# =========================================================================================================
# This script pre-processes PACE OCI data, to:
# 1) get TOA radiance
# 2) merge with the corresponding geolocation data (latitude, longitude, solar zenith angle, etc.)
# 3) select valid pixels (the criteria TBD)
# =========================================================================================================

using PACE_SIF
using NCDatasets, Glob
using Plots    # just for fun XD
# include dependencies - to be updated
include("./set_parameters.jl")
include("./pixel_retrieval.jl")

println("=== Starting script with $(Threads.nthreads()) threads ===")

# ===========================================
# preprocess PACE OCI data
# ==========================================
# input dir
L1B_dir = "/home/zhe2/data/PACE/L1B_V3"
L2AOP_dir = "/home/zhe2/data/PACE/L2_AOP"
L2BGC_dir = "/home/zhe2/data/PACE/L2_BGC"
# interim dir 
interim_dir = "/home/zhe2/data/PACE/interim"

# load one scene - to be updated!
L1B_file_pattern = "PACE_OCI.20250130T202059.L1B.V3.nc"
L2AOP_file_pattern = "PACE_OCI.20250130T202059.L2.OC_AOP.V3_0.nc"
L2BGC_file_pattern = "PACE_OCI.20250130T202059.L2.OC_BGC.V3_0.nc"
# L1B_file_pattern = "PACE_OCI.20250928T2224??.L1B.V3.nc"
# L2AOP_file_pattern = "PACE_OCI.20250928T2224??.L2.OC_AOP.V3_1.nc"
# L2BGC_file_pattern = "PACE_OCI.20250928T2224??.L2.OC_BGC.V3_1.nc"
L1B_file = first(glob(L1B_file_pattern, L1B_dir))
L2AOP_file = first(glob(L2AOP_file_pattern, L2AOP_dir))
L2BGC_file = first(glob(L2BGC_file_pattern, L2BGC_dir))
ds_L1B = NCDataset(L1B_file)
ds_L2AOP = NCDataset(L2AOP_file)
ds_L2BGC = NCDataset(L2BGC_file)

# show all groups and variables in the L1B dataset
println("Groups in L1B dataset:")
for group in keys(ds_L1B)
    println("Group: $group")
    println("Variables in this group:")
    for var in keys(ds_L1B[group])
        println("  - $var")
    end
end

# select variables that I want to Use
var_list_L1B = [
    "red_wavelength", "red_solar_irradiance", "watermask",
    "sensor_azimuth", "sensor_zenith", "solar_azimuth", "solar_zenith",
    "quality_flag", "rhot_red"
    ]
var_list_L2AOP = ["nflh", "l2_flags"]
var_list_L2BGC = ["chlor_a"]

# subset the data
function subset_netcdf_dataset(
        filepath::String, 
        selected_vars::Vector{String}, 
        output_path::String; 
        indices::Dict=Dict(), 
        prefix_groups::Bool=false,
        rename_dims=Dict(),
        post_process_func::Union{Function, Nothing}=nothing
    ) 
    # Open input dataset 
    ds = Dataset(filepath)
    
    # Create output filename
    output_name = "subset_" * basename(filepath)
    output_file = joinpath(output_path, output_name)

    # Create output dataset, overwrite if exists
    merged_data = Dataset(output_file, "c")

    try
        # Keep track of defined dimensions
        defined_dims = Set{String}()
        
        for group_name in keys(ds.group)
            group = ds.group[group_name]
            
            for var_name in keys(group)
                # Check if this variable is in the selected list
                if var_name in selected_vars
                    var = group[var_name]
                    
                    # Determine output variable name
                    key = prefix_groups ? "$(group_name)_$(var_name)" : var_name
                    
                    # Get dimensions
                    var_dims = dimnames(var)

                    println("Variable: $var_name")
                    
                    # Build indexing tuple for subsetting
                    idx_tuple = []
                    subset_dims = []
                    subset_sizes = []
                    
                    for (i, dim_name) in enumerate(var_dims)
                        if haskey(indices, dim_name)
                            # Use provided indices
                            push!(idx_tuple, indices[dim_name])
                            push!(subset_dims, dim_name)
                            push!(subset_sizes, length(indices[dim_name]))
                        else
                            # Use all indices
                            push!(idx_tuple, :)
                            push!(subset_dims, dim_name)
                            push!(subset_sizes, size(var, i))
                        end
                    end
                    
                    # Define dimensions in output dataset if not already defined
                    for (dim_name, dim_size) in zip(subset_dims, subset_sizes)
                        # Use renamed dimension if mapping exists
                        output_dim_name = haskey(rename_dims, dim_name) ? rename_dims[dim_name] : dim_name
                        
                        if !(output_dim_name in defined_dims)
                            defDim(merged_data, output_dim_name, dim_size)
                            push!(defined_dims, output_dim_name)
                        end
                    end
                    
                    # Create variable in output dataset
                    var_type = eltype(var)
                    println("  | Original type: ", var_type)

                    nc_type = if Missing <: var_type
                        Base.nonmissingtype(var_type)
                    else
                        var_type
                    end
                    println("  | Output type: ", nc_type)

                    # Create variable
                    fillvalue = haskey(var.attrib, "_FillValue") ? var.attrib["_FillValue"] : nothing
                    renamed_subset_dims = [haskey(rename_dims, d) ? rename_dims[d] : d for d in subset_dims]
                    if fillvalue !== nothing
                        merged_var = defVar(
                            merged_data, 
                            key, 
                            nc_type, 
                            tuple(renamed_subset_dims...),
                            fillvalue = fillvalue
                        )
                    else 
                        merged_var = defVar(
                            merged_data, 
                            key, 
                            nc_type, 
                            tuple(renamed_subset_dims...)
                        )
                    end

                    # Copy attributes, skipping _FillValue and missing_value
                    for (attr_name, attr_value) in var.attrib
                        println("  | Attribute: $attr_name = $attr_value")
                        if attr_name ∉ ["_FillValue", "missing_value"]
                            try
                                merged_var.attrib[attr_name] = attr_value
                            catch e
                                @warn "Could not copy attribute $attr_name for $key: $e"
                            end
                        end
                    end

                    # rescale data by checking whether scale_factor and add_offset exist
                    if haskey(var.attrib, "scale_factor") && haskey(var.attrib, "add_offset")
                        scale_factor = var.attrib["scale_factor"]
                        add_offset = var.attrib["add_offset"]
                        println("  | Rescaling data with scale_factor=$scale_factor and add_offset=$add_offset")
                        # Load original data
                        original_data = var[idx_tuple...]
                        # Rescale data
                        rescaled_data = rescale_data(original_data, scale_factor, add_offset)
                        # Write rescaled data to output variable
                        merged_var[:] = rescaled_data
                        # add a flag attribute to indicate that the data has been rescaled
                        merged_var.attrib["rescaled"] = "true"
                    else
                        # otherwise keep the original data (NCDatasets handles fill values automatically)
                        merged_var[:] = var[idx_tuple...]
                    end
                    println("  ✓ Copied variable: $key ($(subset_sizes))")
                end
            end
        end
        
        # Copy global attributes
        for (attr_name, attr_value) in ds.attrib
            merged_data.attrib[attr_name] = attr_value
        end

        # post process if needed
        if post_process_func !== nothing
            post_process_func(merged_data)
        end
        println("  ✓ Subset created: $output_file")
        
    finally
        close(ds)
        close(merged_data)
    end

    return output_file
end

# functino to rescale data
function rescale_data(
        x,   # input packed data
        scale_factor,  # scale factor
        offset,        # offset
    )
    # Rescale data
    x_rescaled = (x .- offset) ./ scale_factor;
    return x_rescaled
end

# function to calculate Rtoa from radiance and solar irradiance
function TOA_radiance(ds::Dataset)
    println("Creating TOA radiance variable...")

    # Get the NetCDF variable (not the data)
    rhot_red_var = ds["rhot_red"]
    
    # Get dimension names from the variable
    dims = dimnames(rhot_red_var)
    data_type = eltype(rhot_red_var)
    
    # Load the data
    rhot_red = rhot_red_var[:]
    solar_irradiance = ds["red_solar_irradiance"][:]
    solar_zenith_angle = ds["solar_zenith"][:]
    earth_sun_correction = ds.attrib["earth_sun_distance_correction"]   # unitless

    # Ensure dimensions are compatible for broadcasting
    solar_irradiance = reshape(solar_irradiance, (1, 1, size(solar_irradiance)...))
    solar_zenith_angle = reshape(solar_zenith_angle, (size(solar_zenith_angle)..., 1))

    # Shape check
    println("  | Radiance shape: ", size(rhot_red))
    println("  | Solar irradiance shape: ", size(solar_irradiance))
    println("  | Solar zenith angle shape: ", size(solar_zenith_angle))
    println("  | Earth-sun distance correction: ", earth_sun_correction)

    # Calculate TOA radiance
    Rtoa = rhot_red .* solar_irradiance .* cosd.(solar_zenith_angle) ./ π / earth_sun_correction
    
    # Get the non-missing type for NetCDF
    nc_type = Missing <: data_type ? Base.nonmissingtype(data_type) : data_type
    
    # Create new variable with proper dimensions
    Rtoa_var = defVar(ds, "Rtoa_red", nc_type, dims, fillvalue=nc_type(-9999.0))
    
    # Add attributes
    Rtoa_var.attrib["long_name"] = "Top-of-Atmosphere radiance at Red Band"
    Rtoa_var.attrib["units"] = "W/m^2/sr/μm"
    
    # Write data
    Rtoa_var[:] = Rtoa
    
    println("  | TOA radiance variable created with shape: ", size(Rtoa_var[:]))
    
    return ds
end

subset_L1B = subset_netcdf_dataset(
    L1B_file,
    var_list_L1B,
    interim_dir,
    post_process_func = TOA_radiance
)

subset_L2AOP = subset_netcdf_dataset(
    L2AOP_file, 
    var_list_L2AOP, 
    interim_dir,
    rename_dims = Dict(
        "pixels_per_line" => "pixel",
        "number_of_lines" => "scans"
    )
)

subset_L2BGC = subset_netcdf_dataset(
    L2BGC_file,
    var_list_L2BGC,
    interim_dir,
    rename_dims = Dict(
        "pixels_per_line" => "pixel",
        "number_of_lines" => "scans"
    )
)

# show the dimensions of the subsetted dataset
ds_subset_L1B = NCDataset(subset_L1B)
println("Dimensions in subsetted L1B dataset\n $(ds_subset_L1B.dim)")
ds_subset_L2AOP = NCDataset(subset_L2AOP)
println("Dimensions in subsetted L2AOP dataset\n $(ds_subset_L2AOP.dim)")
ds_subset_L2BGC = NCDataset(subset_L2BGC)
println("Dimensions in subsetted L2BGC dataset\n $(ds_subset_L2BGC.dim)")

# ===========================================
# Pixelwise Retrieval
# ==========================================
# set parameters
nflh_threshold = 0.05;           # threshold for valid pixels based on NFLH
DecompositionMethod = :SVD;    # "NMF" or "SVD"
if_log = true;                 # whether to do log-SVD for transmittance
n     = 10;
ranks = 15;
nPC   = ranks;
nSIF  = 2;
nIter = 25;
thr_Converge = 1e-6;

λ_min = 620.0
λ_max = 860.0

λ_remove_min = 750.0
λ_remove_max = 749.0

λ_bl_ref = [607.99, 610.36, 612.73, 615.14, 617.6, 620.06, 622.53, 669.52, 
670.76, 671.99, 673.24, 674.51, 675.73, 676.97, 678.21, 679.45, 
754.3, 779.33, 867.11, 869.61, 872.13]

scale_factor_SIF = 20

configuration = "Configurations: \n" *
          "Wavelength range: $λ_min - $λ_max nm\n" *
          "SNR degradation range: $λ_remove_min - $λ_remove_max nm\n" *
          "SIF scale factor: $scale_factor_SIF\n" *
          "Decomposition method: $DecompositionMethod with if_log=$if_log\n" *
          "NFLH threshold: $nflh_threshold\n" *
          "NMF rank: $rank\n" *
          "Order of polynomials to fit: $n, Number of retrieval PCs: $nPC, SIF PCs: $nSIF\n" *
          "Number of iterations: $nIter, Convergence threshold: $thr_Converge\n"


println(configuration)

params = setup_retrieval_parameters(
    ds_subset_L1B["red_wavelength"][:],
    ds_subset_L1B["red_solar_irradiance"][:],
    λ_min, λ_max, scale_factor_SIF, DecompositionMethod, if_log, n, nPC, nSIF, nIter, thr_Converge,
);

# apply Retrieval_for_Pixel func.

# subset wavelength range
wvlen_index = findall(λ_min .<= ds_subset_L1B["red_wavelength"][:] .<= λ_max);
# collect variables as needed
R_toa = ds_subset_L1B["Rtoa_red"][:, :, wvlen_index];
sza   = ds_subset_L1B["solar_zenith"][:];
vza   = ds_subset_L1B["sensor_zenith"][:];
nflh  = ds_subset_L2AOP["nflh"][:];
chlor_a = ds_subset_L2BGC["chlor_a"][:];
# flag checks whether nflh is above the threshold and whether chlor_a is valid
flag = Matrix{Bool}( .!ismissing.(chlor_a) .& .!ismissing.(nflh) .& (nflh .> nflh_threshold) );

# try a small subset area
test_R_toa = R_toa[1:10, 1:10, :];
test_sza = sza[1:10, 1:10];
test_vza = vza[1:10, 1:10];
test_nflh = nflh[1:10, 1:10];
test_chlor_a = chlor_a[1:10, 1:10];
test_flag = flag[1:10, 1:10];
retrieval_result = process_all_pixels(
    test_R_toa, test_sza, test_vza, test_nflh, test_chlor_a, test_flag, params
);

# test_pixel_idx = (846, 1470)  # (scan, pixel)
# test_R_toa = R_toa[test_pixel_idx..., :]
# test_sza = sza[test_pixel_idx...]
# test_vza = vza[test_pixel_idx...]
# test_nflh = nflh[test_pixel_idx...]
# test_chlor_a = chlor_a[test_pixel_idx...]
# test_flag = flag[test_pixel_idx...]
# retrieval_result = retrieve_pixel(
#     test_R_toa, test_sza, test_vza, test_nflh, test_chlor_a, test_flag, params
# )

# retrieved state vector
# results = process_all_pixels(
#     R_toa, sza, vza, nflh, chlor_a, flag, params
# );

# save metadata: configuration and parameters
# get time
# timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
# save_dir  = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/retrieval_from_realData/"
# save_name = save_dir * "retrieval_results_20250928T222455" * timestamp * ".jld2"

