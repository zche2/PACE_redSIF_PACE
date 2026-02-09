# =========================================================================================================
# This script 
# =========================================================================================================

using PACE_SIF
using NCDatasets
using JLD2

include("./pixel_retrieval.jl")
include("./pre_process.jl")

# ===========================
# Here the function begins
# ===========================

function granule_retrieval(
        L1B_file::String,
        L2AOP_file::String,
        L2BGC_file::String,
        interim_dir::String,
        output_str::String,
        params::RetrievalParams,
        configuration::String;
        var_list_L1B = ["red_wavelength", "red_solar_irradiance", "watermask","latitude", "longitude",
                        "sensor_azimuth", "sensor_zenith", "solar_azimuth", "solar_zenith", "quality_flag", "rhot_red"],
        var_list_L2AOP = ["nflh"], 
        var_list_L2BGC = ["chlor_a"],
        L1B_post_process = TOA_radiance,  # Default function
        L2AOP_post_process = nothing,     # No post-processing
        L2BGC_post_process = nothing      # No post-processing
    )

    # 1) read, preprocess, and subset the L1B, L2AOP, and L2BGC datasets
    println("Starting granule retrieval for:")
    println("  L1B file: $L1B_file")
    println("  L2AOP file: $L2AOP_file")
    println("  L2BGC file: $L2BGC_file\n")


    subset_L1B = subset_netcdf_dataset(
        L1B_file,
        var_list_L1B,
        interim_dir,
        post_process_func = L1B_post_process
    )

    subset_L2AOP = subset_netcdf_dataset(
        L2AOP_file, 
        var_list_L2AOP, 
        interim_dir,
        rename_dims = Dict(
            "pixels_per_line" => "pixel",
            "number_of_lines" => "scans"
        ),
        post_process_func = L2AOP_post_process
    )

    subset_L2BGC = subset_netcdf_dataset(
        L2BGC_file,
        var_list_L2BGC,
        interim_dir,
        rename_dims = Dict(
            "pixels_per_line" => "pixel",
            "number_of_lines" => "scans"
        ),
        post_process_func = L2BGC_post_process
    )

    # show the dimensions of the subsetted dataset
    ds_subset_L1B   = NCDataset(subset_L1B)
    ds_subset_L2AOP = NCDataset(subset_L2AOP)
    ds_subset_L2BGC = NCDataset(subset_L2BGC)

    # 2) apply the pixel retrieval function to each pixel in the granule
    # subset wavelength range
    wvlen_index = findall(λ_min .<= ds_subset_L1B["red_wavelength"][:] .<= λ_max);
    # collect variables as needed
    lat   = ds_subset_L1B["latitude"][:];
    lon   = ds_subset_L1B["longitude"][:];
    R_toa = ds_subset_L1B["Rtoa_red"][:, :, wvlen_index];
    sza   = ds_subset_L1B["solar_zenith"][:];
    vza   = ds_subset_L1B["sensor_zenith"][:];
    nflh  = ds_subset_L2AOP["nflh"][:];
    chlor_a = ds_subset_L2BGC["chlor_a"][:];
    # flag checks whether nflh is above the threshold and whether chlor_a is valid
    flag = Matrix{Bool}( .!ismissing.(chlor_a) .& .!ismissing.(nflh) .& (nflh .> nflh_threshold) );

    # retrieved state vector
    results = process_all_pixels(
        R_toa, sza, vza, nflh, chlor_a, flag, params
    );
    num_valid_pixels = count(!ismissing, results[:, :, 1]);
    println("  Number of total valid pixels: ", length(flag[flag .== true]))
    println("  Number of valid pixels retrieved: $num_valid_pixels")

    # 3) Save the retrieved state vector
    # convert to NCDataset
    retrieved_state_vector_dims = ("pixel", "scan", "state_vector")
    retrieved_state_vector_var  = "retrieved_state_vector"
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    output_name = joinpath(interim_dir, "SIF_$(output_str)_$timestamp.nc")

    ds_output = Dataset(output_name, "c")
    defDim(ds_output, "pixel", size(results, 1))
    defDim(ds_output, "scan", size(results, 2))
    defDim(ds_output, "state_vector", size(results, 3))
    # add var: retrieved state vector, lat, and lon
    retrieved_var = defVar(ds_output, retrieved_state_vector_var, Float64, retrieved_state_vector_dims, fillvalue=-9999.0)
    lat_var = defVar(ds_output, "latitude", Float64, ("pixel", "scan"), fillvalue=-9999.0)
    lon_var = defVar(ds_output, "longitude", Float64, ("pixel", "scan"), fillvalue=-9999.0)
    # add attributes
    retrieved_var.attrib["configuration"] = configuration
    # write data
    retrieved_var[:] = results
    lat_var[:] = lat
    lon_var[:] = lon
    println("  Retrieved state vector saved to: $output_name")
    close(ds_output)


    # # delete interim L2 dataset
    # rm(subset_L2AOP)
    # rm(subset_L2BGC)

end


# just for archive:

# try a small subset area
# test_R_toa = R_toa[1220:1225, 1707:1710, :];
# test_sza = sza[1220:1225, 1707:1710];
# test_vza = vza[1220:1225, 1707:1710];
# test_nflh = nflh[1220:1225, 1707:1710];
# test_chlor_a = chlor_a[1220:1225, 1707:1710];
# test_flag = flag[1220:1225, 1707:1710];
# retrieval_result = process_all_pixels(
#     test_R_toa, test_sza, test_vza, test_nflh, test_chlor_a, test_flag, params
# );
# # how many non missing results do we have?
# num_valid_pixels = count(!ismissing, retrieval_result[:, :, 1])

# try a single pixel
# test_pixel_idx = (1222, 1710)  # (scan, pixel)
# test_R_toa = R_toa[test_pixel_idx..., :]
# test_sza = sza[test_pixel_idx...]
# test_vza = vza[test_pixel_idx...]
# test_nflh = nflh[test_pixel_idx...]
# test_chlor_a = chlor_a[test_pixel_idx...]
# test_flag = flag[test_pixel_idx...]
# retrieval_result = retrieve_pixel(
#     test_R_toa, test_sza, test_vza, test_nflh, test_chlor_a, test_flag, params
# );
