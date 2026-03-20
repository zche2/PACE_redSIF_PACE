# =========================================================================================================
# Merge L1B subset and L2AOP subset into a single interim NetCDF for batch retrieval.
# Output has dims (pixels, scans) and vars: Rtoa_red, red_wavelength, latitude, longitude, nflh, watermask.
# =========================================================================================================

using NCDatasets

# Map dim name to output name
const _DIM_RENAME = Dict(
    "pixels_per_line" => "pixels",
    "number_of_lines" => "scans",
    "number_of_bands" => "red_wavelength",
    "red_bands" => "red_wavelength",
)

function _to_pixels_scans_bands(arr, dnames)
    out_dims = [get(_DIM_RENAME, d, d) for d in dnames]
    i_pix = findfirst(==("pixels"), out_dims)
    i_scan = findfirst(==("scans"), out_dims)
    i_band = findfirst(==("red_wavelength"), out_dims)
    isnothing(i_band) && (i_band = findfirst(x -> occursin("band", lowercase(x)), out_dims))
    perm = (i_pix, i_scan, i_band)
    return permutedims(arr, perm)
end

function _to_pixels_scans_2d(arr, dnames)
    out_dims = [get(_DIM_RENAME, d, d) for d in dnames]
    i_pix = findfirst(==("pixels"), out_dims)
    i_scan = findfirst(==("scans"), out_dims)
    perm = (i_pix, i_scan)
    return permutedims(arr, perm)
end

"""
    preprocess_and_merge_in_memory(L1B_path, L2AOP_path, L2BGC_path, output_path)

Read L1B/L2AOP/L2BGC directly, compute TOA in memory, write single interim file.
Skips subset files to save I/O and storage. Uses ~2× granule size in RAM.
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
        # L1B: read and TOA
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

        nflh_info = _find_var_from_dataset(ds_aop, "nflh")
        nflh = Float32.(replace(_to_pixels_scans_2d(nflh_info.var[:], nflh_info.dims), missing => NaN))

        chlor_a = nothing
        if L2BGC_path !== nothing && isfile(L2BGC_path)
            ds_bgc = Dataset(L2BGC_path)
            try
                try
                    chl_info = _find_var_from_dataset(ds_bgc, "chlor_a")
                    chlor_a = Float32.(replace(_to_pixels_scans_2d(chl_info.var[:], chl_info.dims), missing => NaN))
                catch
                    # chlor_a optional for filtering
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

        comp = (shuffle=true, deflatelevel=4)
        v_spec = defVar(merged, "Rtoa_red", Float32, ("pixels", "scans", "red_wavelength"), fillvalue=Float32(-9999.0); comp...)
        v_spec.attrib["long_name"] = "Top-of-Atmosphere radiance at Red Band"
        v_spec.attrib["units"] = "W/m^2/sr/μm"
        v_spec[:, :, :] = Rtoa

        v_wl = defVar(merged, "red_wavelength", Float32, ("red_wavelength",))
        v_wl[:] = wl
        v_wl.attrib["units"] = "nm"

        v_lat = defVar(merged, "latitude", Float32, ("pixels", "scans"), fillvalue=Float32(NaN); comp...)
        v_lat[:, :] = lat
        v_lon = defVar(merged, "longitude", Float32, ("pixels", "scans"), fillvalue=Float32(NaN); comp...)
        v_lon[:, :] = lon

        if wm !== nothing
            v_wm = defVar(merged, "watermask", eltype(wm), ("pixels", "scans"), fillvalue=typemin(eltype(wm)); comp...)
            v_wm[:, :] = wm
        end

        v_nflh = defVar(merged, "nflh", Float32, ("pixels", "scans"), fillvalue=Float32(NaN); comp...)
        v_nflh[:, :] = nflh

        if chlor_a !== nothing
            v_chl = defVar(merged, "chlor_a", Float32, ("pixels", "scans"), fillvalue=Float32(NaN); comp...)
            v_chl[:, :] = chlor_a
        end

        merged.attrib["source_L1B"] = basename(L1B_path)
        merged.attrib["source_L2AOP"] = basename(L2AOP_path)
        merged.attrib["history"] = "Merged by preprocess_and_merge_in_memory"

        close(merged)
    finally
        close(ds_l1b)
        close(ds_aop)
    end
    return output_path
end

function merge_global_fit_inputs(
    L1B_subset_path::AbstractString,
    L2AOP_subset_path::AbstractString,
    output_path::AbstractString;
    L2BGC_subset_path::Union{AbstractString, Nothing}=nothing,
)
    ds_l1b = Dataset(L1B_subset_path)
    ds_aop = Dataset(L2AOP_subset_path)

    try
        # L1B subset has Rtoa_red, red_wavelength, latitude, longitude, watermask (after TOA_radiance post-process)
        # L2AOP has nflh
        # Both must have pixels and scans dims (from rename_dims in subset)
        haskey(ds_l1b, "Rtoa_red") || error("L1B subset missing Rtoa_red")
        haskey(ds_l1b, "red_wavelength") || error("L1B subset missing red_wavelength")
        haskey(ds_aop, "nflh") || error("L2AOP subset missing nflh")

        v_spec = ds_l1b["Rtoa_red"]
        dnames = collect(String.(dimnames(v_spec)))
        n_pix = size(v_spec, findfirst(==("pixels"), dnames))
        n_scan = size(v_spec, findfirst(==("scans"), dnames))
        i_band = findfirst(==("red_wavelength"), dnames)
        isnothing(i_band) && (i_band = findfirst(x -> occursin("band", lowercase(x)), dnames))
        isnothing(i_band) && error("Could not find band dimension in Rtoa_red dims: $dnames")
        n_band = size(v_spec, i_band)

        merged = Dataset(output_path, "c")
        defDim(merged, "pixels", n_pix)
        defDim(merged, "scans", n_scan)
        defDim(merged, "red_wavelength", n_band)

        # Copy spectrum (Rtoa_red) - Float32 halves storage vs Float64, compression reduces file size
        v_spec_out = defVar(merged, "Rtoa_red", Float32, ("pixels", "scans", "red_wavelength"), fillvalue=Float32(-9999.0), shuffle=true, deflatelevel=4)
        v_spec_out.attrib["long_name"] = "Top-of-Atmosphere radiance at Red Band"
        v_spec_out.attrib["units"] = "W/m^2/sr/μm"
        raw = ds_l1b["Rtoa_red"][:]
        # Ensure (pixels, scans, bands) layout for batch_fit
        arr_spec = if dnames == ["pixels", "scans", "red_wavelength"] || dnames == ["pixels", "scans", "band"] ||
           dnames == ["pixels", "scans", "number_of_bands"]
            Float32.(replace(raw, missing => NaN))
        elseif dnames == ["scans", "pixels", "red_wavelength"] || dnames == ["scans", "pixels", "band"] ||
               dnames == ["scans", "pixels", "number_of_bands"]
            Float32.(replace(permutedims(raw, (2, 1, 3)), missing => NaN))
        elseif dnames == ["pixels", "red_wavelength", "scans"] || dnames == ["pixels", "band", "scans"]
            Float32.(replace(permutedims(raw, (1, 3, 2)), missing => NaN))
        elseif dnames == ["scans", "red_wavelength", "pixels"] || dnames == ["scans", "band", "pixels"]
            Float32.(replace(permutedims(raw, (3, 1, 2)), missing => NaN))
        else
            error("Unsupported Rtoa_red dims: $dnames")
        end
        v_spec_out[:, :, :] = arr_spec

        # Copy wavelength - Float32 sufficient for nm
        wl = ds_l1b["red_wavelength"][:]
        v_wl = defVar(merged, "red_wavelength", Float32, ("red_wavelength",))
        v_wl[:] = collect(Float32.(replace(wl, missing => NaN)))
        v_wl.attrib["units"] = get(ds_l1b["red_wavelength"].attrib, "units", "nm")

        # Copy geolocation
        for vn in ["latitude", "longitude"]
            haskey(ds_l1b, vn) || continue
            v_in = ds_l1b[vn]
            raw = v_in[:, :]
            d = collect(String.(dimnames(v_in)))
            arr = if d == ["pixels", "scans"]
                raw
            elseif d == ["scans", "pixels"]
                permutedims(raw, (2, 1))
            else
                error("Unsupported $vn dims: $d")
            end
            v_out = defVar(merged, vn, Float32, ("pixels", "scans"), fillvalue=Float32(NaN), shuffle=true, deflatelevel=4)
            v_out[:, :] = Float32.(replace(arr, missing => NaN))
        end

        # Copy watermask
        if haskey(ds_l1b, "watermask")
            v_in = ds_l1b["watermask"]
            raw = v_in[:, :]
            d = collect(String.(dimnames(v_in)))
            arr = if d == ["pixels", "scans"]
                raw
            elseif d == ["scans", "pixels"]
                permutedims(raw, (2, 1))
            else
                error("Unsupported watermask dims: $d")
            end
            T = Base.nonmissingtype(eltype(arr))
            v_out = defVar(merged, "watermask", T, ("pixels", "scans"), fillvalue=typemin(T), shuffle=true, deflatelevel=4)
            v_out[:, :] = T.(replace(arr, missing => typemin(T)))
        end

        # Copy nflh from L2AOP
        v_in = ds_aop["nflh"]
        raw = v_in[:, :]
        d = collect(String.(dimnames(v_in)))
        arr = if d == ["pixels", "scans"]
            raw
        elseif d == ["scans", "pixels"]
            permutedims(raw, (2, 1))
        else
            error("Unsupported nflh dims: $d")
        end
        v_out = defVar(merged, "nflh", Float32, ("pixels", "scans"), fillvalue=Float32(NaN), shuffle=true, deflatelevel=4)
        v_out[:, :] = Float32.(replace(arr, missing => NaN))

        # Optional: chlor_a from L2BGC for filtering
        if L2BGC_subset_path !== nothing && isfile(L2BGC_subset_path)
            ds_bgc = Dataset(L2BGC_subset_path)
            try
                if haskey(ds_bgc, "chlor_a")
                    v_in = ds_bgc["chlor_a"]
                    raw = v_in[:, :]
                    d = collect(String.(dimnames(v_in)))
                    arr = if d == ["pixels", "scans"]
                        raw
                    elseif d == ["scans", "pixels"]
                        permutedims(raw, (2, 1))
                    else
                        @warn "Skipping chlor_a: unsupported dims $d"
                        arr = nothing
                    end
                    if arr !== nothing
                        v_out = defVar(merged, "chlor_a", Float32, ("pixels", "scans"), fillvalue=Float32(NaN), shuffle=true, deflatelevel=4)
                        v_out[:, :] = Float32.(replace(arr, missing => NaN))
                    end
                end
            finally
                close(ds_bgc)
            end
        end

        # Copy useful global attributes
        for (k, v) in ds_l1b.attrib
            k in ["_NCProperties", "_SuperblockVersion"] && continue
            try
                merged.attrib[k] = v
            catch
                nothing
            end
        end

        merged.attrib["source_L1B"] = basename(L1B_subset_path)
        merged.attrib["source_L2AOP"] = basename(L2AOP_subset_path)
        merged.attrib["history"] = "Merged by merge_global_fit_inputs"

        close(merged)
    finally
        close(ds_l1b)
        close(ds_aop)
    end

    return output_path
end
