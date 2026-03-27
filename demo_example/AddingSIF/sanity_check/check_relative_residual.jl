using NCDatasets
using Plots
using Interpolations

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
rhot_raw = rhot_info.var[:]
sol_info = _find_var_from_dataset(ds_l1b, "red_solar_irradiance")
solar_irrad = sol_info.var[:]
sza_info = _find_var_from_dataset(ds_l1b, "solar_zenith")
sza = sza_info.var[:]
earth_sun = ds_l1b.attrib["earth_sun_distance_correction"]

# reshape
solar_irradiance = reshape(solar_irrad, (1, 1, size(solar_irrad)...))
solar_zenith_angle = reshape(sza, (size(sza)..., 1))

# compute TOA radiance
Rtoa = rhot_raw .* solar_irradiance .* cosd.(solar_zenith_angle) ./ π ./ earth_sun
Rtoa = _to_pixels_scans_bands(Rtoa, rhot_info.dims)

# load retrieved SIF (residuals)
addedSIFfile = "/home/zhe2/data/PACE/adding_sif_output/PACE_OCI.20250927T123954.L1B.V3_adding_sif_p5_20260320.nc"
ds_sif = Dataset(addedSIFfile)
rmse = ds_sif["rmse"].var[:]
red_wavelength = ds_sif["red_wavelength"].var[:]

# interpolate Rtoa on red_wavelength, Rtoa is a 3D array (pixels, scans, bands)
wvlen_l1b_info = _find_var_from_dataset(ds_l1b, "red_wavelength")
wvlen_l1b = wvlen_l1b_info.var[:]

# Gridded: interpolate(knots_tuple, array, scheme_tuple)
itp = interpolate(
    (1:size(Rtoa,1), 1:size(Rtoa,2), wvlen_l1b),
    Rtoa,
    (NoInterp(), NoInterp(), Gridded(Linear()))
)

# Query at target wavelengths — integer indexing for dim 1 & 2, nm for dim 3
Rtoa_subset = [itp[i, j, w] for i in axes(Rtoa,1), j in axes(Rtoa,2), w in red_wavelength]

# compute relative residual
relative_residual = rmse ./ Rtoa_subset

# plot
# choose only where sif is retrieved
sif = ds_sif["sif_radiance_678nm"].var[:]
idx_sif = findall(sif .> 0)
# Each ci is CartesianIndex(pixel, scan)
Rtoa_subset_sif = stack([Rtoa_subset[ci[1], ci[2], :] for ci in idx_sif])
Rtoa_full_sif = stack([Rtoa[ci[1], ci[2], :] for ci in idx_sif])
relative_residual_sif = stack([relative_residual[ci[1], ci[2], :] for ci in idx_sif])

# --- figure 1: radiance & subsetted radiance
# shape: (83, n_sif_pixels)  — transpose if needed
p1 = plot(
    red_wavelength, Rtoa_subset_sif[:, 1:10], 
    size=(800, 300), legend=false, dpi=300,
    margin=10Plots.mm,
    xlabel="Wavelength [nm]",
    ylabel="TOA radiance [W/m²/µm/sr]",
    linewidth=2.5,
)
plot!(p1, wvlen_l1b, Rtoa_full_sif[:, 1:10], linestyle=:dash)

# --- figure 2: relative residual
p2 = plot(
    red_wavelength, relative_residual_sif[:, 1:1000:end], 
    size=(800, 300), legend=false, dpi=300,
    margin=10Plots.mm,
    xlabel="Wavelength [nm]",
    ylabel="Relative residual",
    linewidth=1,
    alpha=0.1,
    ylims=(-0.02, 0.02),
    grid=true,
    xticks=640:10:760,
)
