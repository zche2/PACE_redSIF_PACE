using NCDatasets
using Plots
using Interpolations
using Statistics

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

L1B_path = "/home/zhe2/data/PACE/L1B_V3/PACE_OCI.20260125T163441.L1B.V3.nc"
L2AOP_path = "/home/zhe2/data/PACE/L2_AOP_V3.1/PACE_OCI.20260125T163441.L2.OC_AOP.V3_1.nc"
L2BGC_path = "/home/zhe2/data/PACE/L2_BGC_V3.1/PACE_OCI.20260125T163441.L2.OC_BGC.V3_1.nc"

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
vza_info = _find_var_from_dataset(ds_l1b, "sensor_zenith")
vza = vza_info.var[:]

# reshape
solar_irradiance = reshape(solar_irrad, (1, 1, size(solar_irrad)...))
solar_zenith_angle = reshape(sza, (size(sza)..., 1))

# compute TOA radiance
Rtoa = rhot_raw .* solar_irradiance .* cosd.(solar_zenith_angle) ./ π ./ earth_sun
Rtoa = _to_pixels_scans_bands(Rtoa, rhot_info.dims)

# load retrieved SIF (residuals)
addedSIFfile = "/home/zhe2/data/PACE/adding_sif_output/PACE_OCI.20260125T163441.L1B.V3_adding_sif_zero_20260401.nc"
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
    red_wavelength, relative_residual_sif[:, 1:8000:end], 
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

# --- figure 3: ratio of modelled to observed radiance
g = relative_residual_sif .+ 1   # g: vicarious gain
# mean gain
mean_g = mean(g, dims=2)

p3 = plot(
    red_wavelength, g[:, 1:1000:end], 
    size=(800, 300), legend=false, dpi=300,
    margin=10Plots.mm,
    xlabel="Wavelength [nm]",
    ylabel="Vicarious gain (g)",
    ylims=(0.99, 1.01),
    linewidth=1,
    alpha=0.1,
)

plot!(p3, red_wavelength, mean_g, linestyle=:solid, color=:black, label="Mean gain", linewidth=2)

# --- figure 4: mean gain vs sza
sza_stack = stack([sza[ci[1], ci[2]] for ci in idx_sif])

# set sza and vza bins
sza_bins = 20:2.5:45
sza_binned = [findfirst(x -> x >= k, sza_bins) for k in sza_stack]
mean_g_sza = [mean(g[:, findall(sza_binned .== i)], dims=2) for i in 1:length(sza_bins)]
counts = [length(findall(sza_binned .== i)) for i in 1:length(sza_bins)]
# set the color to continuous cmap
n_bins = length(sza_bins)
bin_colors = reshape([cgrad(:viridis)[i] for i in range(0, 1, length=n_bins)], 1, n_bins)
# set labels and convert it to vector
labels = label = reshape([string(v) * "°" for v in sza_bins], 1, n_bins)
p4 = plot(
    red_wavelength, mean_g_sza, 
    color=bin_colors,
    size=(1400, 300),
    dpi=300,
    margin=10Plots.mm,
    linewidth=2,
    alpha=0.8,
    labels=labels,
    legend_position=:outerright,
    legend_title="Solar zenith angle [°]",
)

# --- figure 5: mean gain vs vza
vza_stack = stack([vza[ci[1], ci[2]] for ci in idx_sif])
vza_bins = 20:5:80
vza_binned = [findfirst(x -> x >= k, vza_bins) for k in vza_stack]
mean_g_vza = [mean(g[:, findall(vza_binned .== i)], dims=2) for i in 1:length(vza_bins)]
counts = [length(findall(vza_binned .== i)) for i in 1:length(vza_bins)]
# set the color to continuous cmap
n_bins = length(vza_bins)
bin_colors = reshape([cgrad(:viridis)[i] for i in range(0, 1, length=n_bins)], 1, n_bins)
# set labels and convert it to vector
labels = label = reshape([string(v) * "°" for v in vza_bins], 1, n_bins)
p5 = plot(
    red_wavelength, mean_g_vza, 
    color=bin_colors,
    size=(1400, 300),
    dpi=300,
    margin=10Plots.mm,
    linewidth=2,
    alpha=0.8,
    labels=labels,
    legend_position=:outerright,
    legend_title="View zenith angle [°]",
)

# --- figure 6: mean gain vs. chl
chl_info = _find_var_from_dataset(ds_bgc, "chlor_a")
chl = chl_info.var[:]
chl_stack = [chl[ci[1], ci[2]] for ci in idx_sif]
# remove missing chl (filter both chl_stack and corresponding g columns)
valid_chl_idx = findall(.!ismissing.(chl_stack))
chl_stack = Float64.(chl_stack[valid_chl_idx])
g_chl = g[:, valid_chl_idx]
# bins - log scale
chl_bins = exp10.(range(-1.5, 1, length=10))
chl_binned = [findfirst(x -> x >= k, chl_bins) for k in chl_stack]
mean_g_chl = [mean(g_chl[:, findall(chl_binned .== i)], dims=2) for i in 1:length(chl_bins)]
counts = [length(findall(chl_binned .== i)) for i in 1:length(chl_bins)]
# set the color to continuous cmap
n_bins = length(chl_bins)
bin_colors = reshape([cgrad(:viridis)[i] for i in range(0, 1, length=n_bins)], 1, n_bins)
# set labels and convert it to vector
label = reshape(["$(round(v, digits=2)) [$(cnt)]" for (v, cnt) in zip(chl_bins, counts)], 1, n_bins)
p6 = plot(
    red_wavelength, mean_g_chl, 
    color=bin_colors,
    size=(1600, 300),
    dpi=300,
    margin=10Plots.mm,
    linewidth=2,
    alpha=0.8,
    label=label,
    # legend_columns=4,
    legend_position=:outerright,
    legend_title="Chl concentration [counts]",
)