using JLD2
using NCDatasets, Glob
using LinearAlgebra
using Plots
using PACE_SIF
using LegendrePolynomials
using Statistics
using Dates


#=============================================================================
# Configuration
=============================================================================#
dir = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/retrieval_from_realData/20250927/"  # Replace with your actual retrieval directory
output_dir = "/home/zhe2/FraLab/PACE_redSIF_PACE/global_fit/gridded_output"
mkpath(output_dir)

# params file
params_file = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/retrieval_from_realData/20250927/retrieval_params_20250927.jld2"
@load params_file params

# Define global grid (0.1° resolution)
lat_edges = -90.0:0.1:90.0
lon_edges = -180.0:0.1:180.0
lat_centers = (lat_edges[1:end-1] .+ lat_edges[2:end]) ./ 2
lon_centers = (lon_edges[1:end-1] .+ lon_edges[2:end]) ./ 2

#=============================================================================
# Helper Function to Find Corresponding L1B File
=============================================================================#
function find_l1b_file(retrieval_file::String, l1b_dir::String)
    # Extract identifier from retrieval file (adjust pattern as needed)
    base = basename(retrieval_file)
    # find YYYYMMDDTHHMMSS out of the retrieval file name
    m = base[5:19]
    # substitute into "PACE_OCI.YYYYMMDDTHHMMSS.L1B.V3.nc"
    l1b_name = "subset_PACE_OCI." * m * ".L1B.V3.nc"
    l1b_path = joinpath(l1b_dir, l1b_name)
    println("Looking for L1B file: ", l1b_path)
    return isfile(l1b_path) ? l1b_path : nothing
end

#=============================================================================
# Reconstruction Function (from retrieval code)
=============================================================================#
function reconstruct_components_1d(
    x::Vector, params::RetrievalParams; sza=nothing, E=nothing, Rtoa_obs=nothing, if_log=true
)
    nPoly = params.nPoly
    nPC = params.nPC
    nSIF = params.nSIF

    # 1. Reflectance (Legendre polynomial expansion)
    λc = params.λc  # Central wavelength for normalization
    v = [collectPl(λ, lmax=nPoly) for λ in λc]  # Legendre polynomials
    ρ = hcat(v...)' * x[1:(nPoly + 1)]

    # 2. Transmittance
    trans_mat = params.PrinComp[:, 1:nPC]
    β = x[(nPoly + 2):(nPoly + nPC + 1)]

    if if_log
        logT₁ = trans_mat * β
        T₁ = exp.(-logT₁)
    else
        T₁ = trans_mat * β
    end

    # 3. Smoothing parameter and T₂
    smooth_x = 10.0 / (1.0 + exp(-x[nPoly + nPC + 2])) + 1.0

    if if_log
        T₂ = exp.(-smooth_x .* (trans_mat * β))
    else
        T₂ = exp.(smooth_x .* log.(T₁))
    end

    # 4. SIF
    SIF_shape = params.SIFComp[:, 1:nSIF]
    SIF_coeff = x[(nPoly + nPC + 3):(nPoly + nPC + nSIF + 2)]
    SIF = SIF_shape * SIF_coeff

    # 5. TOA Radiance (if angles provided)
    radiance = nothing
    if !isnothing(sza) && !isnothing(E)
        radiance = @. E * cosd(sza) / π * T₂ * ρ + SIF * T₁
    elseif !isnothing(params.E) && !isnothing(sza)
        radiance = @. params.E * cosd(sza) / π * T₂ * ρ + SIF * T₁
    end

    # 6. Residual (if Rtoa_obs provided)
    residual = nothing
    if !isnothing(Rtoa_obs) && !isnothing(radiance)
        residual = Rtoa_obs .- radiance
    end

    return Dict(
        "reflectance" => ρ,
        "T1" => T₁,
        "T2" => T₂,
        "SIF" => SIF,
        "radiance" => radiance,
        "residual" => residual,
        "smooth_param" => smooth_x,
        "beta" => β,
        "gamma" => x[1:(nPoly + 1)],
        "SIF_coeff" => SIF_coeff
    )
end

#=============================================================================
# Main Processing Loop
=============================================================================#

# Get all retrieval files
retrieval_files = glob("SIF_*.nc", dir)
# Initialize global gridded arrays (accumulate across all files)
gridded_SIF_sum = zeros(Float64, length(lat_centers), length(lon_centers))
gridded_residual_sum = zeros(Float64, length(lat_centers), length(lon_centers))
gridded_count = zeros(Int32, length(lat_centers), length(lon_centers))

for file in retrieval_files
    println("Processing: ", basename(file))
    start_time = now()
    
    # Find corresponding L1B file
    l1b_file = find_l1b_file(file, dir)
    if isnothing(l1b_file)
        println("Warning: L1B file not found for ", basename(file))
        continue
    end
    
    # Read retrieval data
    ds = NCDataset(file, "r")
    state_vector = ds["retrieved_state_vector"][:]  # [scan x pixel x n_params]
    latitude = ds["latitude"][:]
    longitude = ds["longitude"][:]
    close(ds)
    
    # Read L1B data
    ds_l1b   = NCDataset(l1b_file, "r")
    red_wavelength = ds_l1b["red_wavelength"][:]
    λ_range = findall(red_wavelength .>= minimum(params.λ) .&& red_wavelength .<= maximum(params.λ))
    Rtoa_obs = ds_l1b["Rtoa_red"][:, :, λ_range]; # [scan x pixel x wavelength]
    sza = ds_l1b["solar_zenith"][:]  # Adjust variable name as needed
    close(ds_l1b)
    
    # Grid the data
    for i in axes(state_vector, 1), j in axes(state_vector, 2)
        lat = latitude[i, j]
        lon = longitude[i, j]
        
        # Find grid indices
        lat_idx = searchsortedfirst(lat_edges, lat) - 1
        lon_idx = searchsortedfirst(lon_edges, lon) - 1
        
        # Check bounds
        if 1 <= lat_idx <= length(lat_centers) && 1 <= lon_idx <= length(lon_centers)
            # Reconstruct components
            result = reconstruct_components_1d(
                state_vector[i, j, :], params,
                E=params.E,
                sza=sza[i, j], Rtoa_obs=Rtoa_obs[i, j, :]
            )
            
            # Extract SIF at 681 nm
            λ_idx = argmin(abs.(params.λ .- 681.0))
            SIF_value = result["SIF"][λ_idx]
            
            # Calculate residual norm
            residual_norm = isnothing(result["residual"]) ? NaN : sqrt(sum(result["residual"].^2))
            
            # Accumulate
            if !ismissing(SIF_value) && !isnan(residual_norm)
                gridded_count[lat_idx, lon_idx] += 1
                gridded_SIF_sum[lat_idx, lon_idx] += SIF_value
                gridded_residual_sum[lat_idx, lon_idx] += residual_norm
            end
        end
    end
    # show time elapse
    end_time = now()
    println("Time taken for ", basename(file), ": ", end_time - start_time)
end

# Calculate mean and create Union{Missing, Float64} arrays
gridded_SIF_mean = Array{Union{Missing, Float64}}(undef, length(lat_centers), length(lon_centers))
gridded_residual_norm = Array{Union{Missing, Float64}}(undef, length(lat_centers), length(lon_centers))

for i in eachindex(lat_centers), j in eachindex(lon_centers)
    if gridded_count[i, j] > 0
        gridded_SIF_mean[i, j] = gridded_SIF_sum[i, j] / gridded_count[i, j]
        gridded_residual_norm[i, j] = gridded_residual_sum[i, j] / gridded_count[i, j]
    else
        gridded_SIF_mean[i, j] = missing
        gridded_residual_norm[i, j] = missing
    end
end

# Save single gridded output file
output_file = joinpath(output_dir, "gridded_all_20250927.nc")
ds_out = NCDataset(output_file, "c")

defDim(ds_out, "lat", length(lat_centers))
defDim(ds_out, "lon", length(lon_centers))

defVar(ds_out, "latitude", lat_centers, ("lat",))
defVar(ds_out, "longitude", lon_centers, ("lon",))
defVar(ds_out, "SIF_mean", gridded_SIF_mean, ("lat", "lon"))
defVar(ds_out, "residual_l2norm", gridded_residual_norm, ("lat", "lon"))
defVar(ds_out, "count", gridded_count, ("lat", "lon"))

close(ds_out)
println("Saved: ", output_file)

println("Gridding complete!")
println("  | Total grid cells with data: ", count(!ismissing, gridded_SIF_mean))
