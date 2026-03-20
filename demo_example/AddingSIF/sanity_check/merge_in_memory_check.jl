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

function rescale_data(
    x,
    scale_factor,
    offset,
)
    return (x .- offset) .* scale_factor
end

# Find variable and return (var, dims) from dataset (searches groups and root)
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
        return (var=var, dims=collect(String.(dimnames(var))))
    end
    # Also check root
    if haskey(ds, varname)
        var = ds[varname]
        return (var=var, dims=collect(String.(dimnames(var))))
    end
    error("Variable '$varname' not found in dataset")
end

# Read var with optional rescale (scale_factor, add_offset)
function _read_var_rescaled(var)
    data = var[:]
    if haskey(var.attrib, "scale_factor") && haskey(var.attrib, "add_offset")
        data = rescale_data(data, var.attrib["scale_factor"], var.attrib["add_offset"])
    end
    return data
end

L1B_path = "/home/zhe2/data/PACE/L1B_V3/PACE_OCI.20250927T123954.L1B.V3.nc"
L2AOP_path = "/home/zhe2/data/PACE/L2_AOP_V3.1/PACE_OCI.20250927T123954.L2.OC_AOP.V3_1.nc"
L2BGC_path = "/home/zhe2/data/PACE/L2_BGC_V3.1/PACE_OCI.20250927T123954.L2.OC_BGC.V3_1.nc"

ds_l1b = Dataset(L1B_path)
ds_aop = Dataset(L2AOP_path)
ds_bgc = Dataset(L2BGC_path)

rhot_info = _find_var_from_dataset(ds_l1b, "rhot_red")
rhot_raw = _read_var_rescaled(rhot_info.var)
sol_info = _find_var_from_dataset(ds_l1b, "red_solar_irradiance")
solar_irrad = _read_var_rescaled(sol_info.var)
sza_info = _find_var_from_dataset(ds_l1b, "solar_zenith")
sza = _read_var_rescaled(sza_info.var)
earth_sun = ds_l1b.attrib["earth_sun_distance_correction"]
