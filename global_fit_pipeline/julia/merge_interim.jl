# L1B + L2 in-memory merge → interim NetCDF for SVD/batch retrieval (self-contained).
# Variable lookup follows the same group-aware convention as the global-fit merge utilities.

using NCDatasets

function _find_var_from_dataset(ds, varname::String)
    groups_to_check = Pair{String, Any}[]
    try
        if hasproperty(ds, :group)
            for group_name in keys(ds.group)
                push!(groups_to_check, group_name => ds.group[group_name])
            end
        end
    catch
        nothing
    end
    if isempty(groups_to_check)
        push!(groups_to_check, "" => ds)
    end
    for (_, group) in groups_to_check
        haskey(group, varname) || continue
        var = group[varname]
        try
            dimnames(var)
        catch
            continue
        end
        return (var = var, dims = collect(String.(dimnames(var))))
    end
    if haskey(ds, varname)
        var = ds[varname]
        return (var = var, dims = collect(String.(dimnames(var))))
    end
    error("Variable '$varname' not found in dataset")
end

"""Same search as `_find_var_from_dataset`, but returns `nothing` if the variable is absent (any group or root)."""
function _find_var_from_dataset_optional(ds, varname::String)
    groups_to_check = Pair{String, Any}[]
    try
        if hasproperty(ds, :group)
            for group_name in keys(ds.group)
                push!(groups_to_check, group_name => ds.group[group_name])
            end
        end
    catch
        nothing
    end
    if isempty(groups_to_check)
        push!(groups_to_check, "" => ds)
    end
    for (_, group) in groups_to_check
        haskey(group, varname) || continue
        var = group[varname]
        try
            dimnames(var)
        catch
            continue
        end
        return (var = var, dims = collect(String.(dimnames(var))))
    end
    if haskey(ds, varname)
        var = ds[varname]
        return (var = var, dims = collect(String.(dimnames(var))))
    end
    return nothing
end

const _MERGE_DIM_RENAME = Dict(
    "pixels_per_line" => "pixels",
    "number_of_lines" => "scans",
    "number_of_bands" => "red_wavelength",
    "red_bands" => "red_wavelength",
)

function _to_pixels_scans_bands(arr, dnames)
    out_dims = [get(_MERGE_DIM_RENAME, d, d) for d in dnames]
    i_pix = findfirst(==("pixels"), out_dims)
    i_scan = findfirst(==("scans"), out_dims)
    i_band = findfirst(==("red_wavelength"), out_dims)
    isnothing(i_band) && (i_band = findfirst(x -> occursin("band", lowercase(x)), out_dims))
    perm = (i_pix, i_scan, i_band)
    return permutedims(arr, perm)
end

function _to_pixels_scans_2d(arr, dnames)
    out_dims = [get(_MERGE_DIM_RENAME, d, d) for d in dnames]
    i_pix = findfirst(==("pixels"), out_dims)
    i_scan = findfirst(==("scans"), out_dims)
    perm = (i_pix, i_scan)
    return permutedims(arr, perm)
end

"""
    preprocess_and_merge_in_memory(L1B_path, L2AOP_path, L2BGC_path, output_path)

Read L1B/L2AOP/L2BGC directly, compute TOA in memory, write single interim file.
"""
function preprocess_and_merge_in_memory(
    L1B_path::AbstractString,
    L2AOP_path::AbstractString,
    L2BGC_path::Union{AbstractString, Nothing},
    output_path::AbstractString,
)
    ds_l1b = Dataset(L1B_path)
    ds_aop = Dataset(L2AOP_path)
    try
        rhot_info = _find_var_from_dataset(ds_l1b, "rhot_red")
        rhot_raw = rhot_info.var[:]
        sol_info = _find_var_from_dataset(ds_l1b, "red_solar_irradiance")
        solar_irrad = sol_info.var[:]
        sza_info = _find_var_from_dataset(ds_l1b, "solar_zenith")
        sza = sza_info.var[:]
        earth_sun = ds_l1b.attrib["earth_sun_distance_correction"]

        solar_irrad = reshape(solar_irrad, (1, 1, size(solar_irrad)...))
        sza = reshape(sza, (size(sza)..., 1))
        Rtoa = Float32.(replace(rhot_raw .* solar_irrad .* cosd.(sza) ./ π / earth_sun, missing => NaN))
        Rtoa = _to_pixels_scans_bands(Rtoa, rhot_info.dims)

        wl_info = _find_var_from_dataset(ds_l1b, "red_wavelength")
        wl = Float32.(replace(wl_info.var[:], missing => NaN))

        lat_info = _find_var_from_dataset(ds_l1b, "latitude")
        lat = Float32.(replace(_to_pixels_scans_2d(lat_info.var[:], lat_info.dims), missing => NaN))
        lon_info = _find_var_from_dataset(ds_l1b, "longitude")
        lon = Float32.(replace(_to_pixels_scans_2d(lon_info.var[:], lon_info.dims), missing => NaN))

        wm = nothing
        if haskey(ds_l1b, "watermask")
            wm_info = _find_var_from_dataset(ds_l1b, "watermask")
            wm_raw = wm_info.var[:]
            T = Base.nonmissingtype(eltype(wm_raw))
            wm = T.(replace(_to_pixels_scans_2d(wm_raw, wm_info.dims), missing => typemin(T)))
        end

        # Optional L2 OC_AOP field: used only if [batch_fit].pixel_filter_vars includes "nflh".
        nflh_info = _find_var_from_dataset_optional(ds_aop, "nflh")
        nflh = if nflh_info !== nothing
            Float32.(replace(_to_pixels_scans_2d(nflh_info.var[:], nflh_info.dims), missing => NaN))
        else
            nothing
        end

        chlor_a = nothing
        if L2BGC_path !== nothing && isfile(L2BGC_path)
            ds_bgc = Dataset(L2BGC_path)
            try
                try
                    chl_info = _find_var_from_dataset(ds_bgc, "chlor_a")
                    chlor_a = Float32.(replace(_to_pixels_scans_2d(chl_info.var[:], chl_info.dims), missing => NaN))
                catch
                end
            finally
                close(ds_bgc)
            end
        end

        n_pix, n_scan, n_band = size(Rtoa)
        merged = Dataset(output_path, "c")
        defDim(merged, "pixels", n_pix)
        defDim(merged, "scans", n_scan)
        defDim(merged, "red_wavelength", n_band)

        comp = (shuffle = true, deflatelevel = 4)
        v_spec = defVar(merged, "Rtoa_red", Float32, ("pixels", "scans", "red_wavelength"), fillvalue = Float32(-9999.0); comp...)
        v_spec.attrib["long_name"] = "Top-of-Atmosphere radiance at Red Band"
        v_spec.attrib["units"] = "W/m^2/sr/μm"
        v_spec[:, :, :] = Rtoa

        v_wl = defVar(merged, "red_wavelength", Float32, ("red_wavelength",))
        v_wl[:] = wl
        v_wl.attrib["units"] = "nm"

        v_lat = defVar(merged, "latitude", Float32, ("pixels", "scans"), fillvalue = Float32(NaN); comp...)
        v_lat[:, :] = lat
        v_lon = defVar(merged, "longitude", Float32, ("pixels", "scans"), fillvalue = Float32(NaN); comp...)
        v_lon[:, :] = lon

        if wm !== nothing
            v_wm = defVar(merged, "watermask", eltype(wm), ("pixels", "scans"), fillvalue = typemin(eltype(wm)); comp...)
            v_wm[:, :] = wm
        end

        if nflh !== nothing
            v_nflh = defVar(merged, "nflh", Float32, ("pixels", "scans"), fillvalue = Float32(NaN); comp...)
            v_nflh[:, :] = nflh
        end

        if chlor_a !== nothing
            v_chl = defVar(merged, "chlor_a", Float32, ("pixels", "scans"), fillvalue = Float32(NaN); comp...)
            v_chl[:, :] = chlor_a
        end

        merged.attrib["source_L1B"] = basename(L1B_path)
        merged.attrib["source_L2AOP"] = basename(L2AOP_path)
        merged.attrib["history"] = "Merged by preprocess_and_merge_in_memory (global_fit_pipeline)"

        close(merged)
    finally
        close(ds_l1b)
        close(ds_aop)
    end
    return output_path
end
