# Activate project environment
import Pkg
Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")

# Load packages
using JLD2, Interpolations, Revise
using Base.Threads, Dates
using ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics
using Polynomials, Random
using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools
using PACE_SIF

println("Running with $(Threads.nthreads()) threads")

# ===========================================
# CONFIGURATION PARAMETERS
# ===========================================

#====== Wavelength range ======#
λ_min = 620.0
λ_max = 860.0

λ_remove_min = 750.0
λ_remove_max = 749.0

#====== Pseudo observation generation ======#
n_sample = 5000
random_seed = 512

# Polynomial fitting
order = 4

# Baseline wavelengths for fitting
λ_bl_ref = [607.99, 610.36, 612.73, 615.14, 617.6, 620.06, 622.53, 669.52, 
            670.76, 671.99, 673.24, 674.51, 675.73, 676.97, 678.21, 679.45, 
            754.3, 779.33, 867.11, 869.61, 872.13]

# SIF scaling
scale_factor_SIF = 20

# File paths
path_transmittance_summer = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc"
path_transmittance_winter = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc"
path_oci = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/sample_granule_20240830T131442_new_chl.nc"
path_sif_shapes = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/SIF_singular_vector.jld2"
path_snr = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt"

nflh_threshold = 0.05  # Threshold for valid NFLH

#====== retrieval ======#
rank  = 15;
n     = 10;
nPC   = rank;
nSIF  = 2;
nIter = 25;
thr_Converge = 1e-6;

println("\n=== Configuration ===")
println("Wavelength range: $λ_min - $λ_max nm")
println("SNR degradation range: $λ_remove_min - $λ_remove_max nm")
println("Polynomial order: $order")
println("SIF scale factor: $scale_factor_SIF")
println("Number of samples: $n_sample")
println("Random seed: $random_seed")
println("NFLH threshold: $nflh_threshold")
println("NMF rank: $rank")
println("Order of polynomials to fit: $n, Number of retrieval PCs: $nPC, SIF PCs: $nSIF")
println("Number of iterations: $nIter, Convergence threshold: $thr_Converge")
println("====================\n")


# ===========================================
# Load Data
# ===========================================

println("Loading data...")

# Load MERRA2 transmittance data
println("Loading transmittance data...")
summer = Dataset(path_transmittance_summer)
winter = Dataset(path_transmittance_winter)

trans = cat(summer["transmittance"][:, :], winter["transmittance"][:, :], dims=1)
bands = summer["band"][:]

close(summer)
close(winter)
println("Transmittance data loaded: $(size(trans, 1)) spectra")

# Load PACE OCI data
println("Loading PACE OCI data...")
oci = Dataset(path_oci)

red_band = oci["red_wavelength"][:]
nflh = oci["nflh"][:, :]
vza = oci["sensor_zenith"][:, :]
sza = oci["solar_zenith"][:, :]
chlor_a = oci["chlor_a"][:, :]

# Select spectral band
ind = findall(λ_min .< red_band .< λ_max)

E = oci["red_solar_irradiance"][ind]
R_toa = oci["radiance_red"][:, :, ind]
oci_band = red_band[ind]

close(oci)
println("PACE data loaded: Band selected to [$λ_min, $λ_max] nm")

# Interpolate transmittance to OCI bands
println("Interpolating transmittance to OCI bands...")
trans_new = zeros(size(trans, 1), length(oci_band))

@threads for i in 1:size(trans, 1)
    itp_row = LinearInterpolation(bands, trans[i, :], extrapolation_bc=0)
    trans_new[i, :] = itp_row.(oci_band)
end
println("Transmittance interpolated: $(size(trans_new, 1)) spectra")

# Load and interpolate SIF shapes
println("Loading SIF shapes...")
SIF_shape_dict = JLD2.load(path_sif_shapes)

# SVD to SIF
SIF_SVD = Spectral_SVD(
    SIF_shape_dict["SIF_shapes"]*scale_factor_SIF,
    SIF_shape_dict["SIF_wavelen"],
    Float64.(collect(skipmissing(oci_band))),
    if_log = false
)

SIF_new = SIF_SVD.PrinComp * SIF_SVD.Loading
println("SIF shapes SVD completed!")

# Load SNR
println("Loading SNR data...")
lines = readlines(path_snr)
end_header_index = findfirst(x -> x == "/end_header", lines)
data = readdlm(path_snr, String, skipstart=end_header_index)

FPA = data[:, 1]
wvlen = parse.(Float64, data[:, 2])

wv_val = (λ_min .< wvlen .< λ_max)
snr_ind = findall((FPA .== "Red") .& wv_val)

c1 = parse.(Float64, data[snr_ind, 4])
c2 = parse.(Float64, data[snr_ind, 5])

println("Data loading complete!\n")

# Print summary
println("Data Summary:")
println("  Wavelength range: $λ_min - $λ_max nm ($(length(oci_band)) bands)")
println("  TOA radiance: $(size(R_toa)) pixels")
println("  Transmittance: $(size(trans_new, 1)) spectra")
println("  SIF shapes: $(size(SIF_new, 2)) components")
println("  SNR coefficients: $(length(c1)) bands")
println()


# ===========================================
# Generating Pseudo Observations
# ===========================================

# Wavelength setup
λ = oci_band
λc = center_wavelength(oci_band)

# Find transmittance baseline
bl_ind = map(λ_bl_ref -> argmin(abs.(λ .- λ_bl_ref)), λ_bl_ref)

# Extract valid pixels
valid_mask = findall(coalesce.(nflh .< nflh_threshold, false))
R_noSIF = R_toa[valid_mask, :]
sza_noSIF = sza[valid_mask]
vza_noSIF = vza[valid_mask]

# Fit polynomials
n_pixels = size(R_noSIF, 1)
K₀ = hcat(collectPl.(λc[bl_ind], lmax=order)...)'
K₀_recon = hcat(collectPl.(λc, lmax=order)...)'
K₀_inv = (K₀'K₀) \ K₀'

# Preallocate
R_fitted = zeros(n_pixels, length(λc))
coeffs_record = zeros(n_pixels, order + 1)

println("Processing $n_pixels pixels...")

# Parallel fitting
@threads for i in 1:n_pixels
    rhs = R_noSIF[i, bl_ind] .* π ./ (E[bl_ind] .* cosd(sza_noSIF[i]))
    coeffs_record[i, :] = K₀_inv * rhs
    R_fitted[i, :] = (K₀_recon * coeffs_record[i, :]) .* (E .* cosd(sza_noSIF[i]) ./ π)
    
    if i % 1000 == 0
        println("Processed $i / $n_pixels")
    end
end

println("Fitting completed!")

# Generate pseudo observations
println("Generating $n_sample pseudo observations...")

# Random sampling
Random.seed!(random_seed)
nᵨ = n_pixels
nₛ = size(SIF_new, 2)
nₜ = size(trans_new, 1)
indᵨ = rand(1:nᵨ, n_sample)
indₛ = rand(1:nₛ, n_sample)
indₜ₁ = rand(1:nₜ, n_sample)
indₜ₂ = rand(1:nₜ, n_sample)
ind_sza = rand(1:nᵨ, n_sample)
ind_vza = rand(1:nᵨ, n_sample)

# Preallocate
len_λ = length(λ)
ρ_all = zeros(n_sample, len_λ)
μ₁_all = zeros(n_sample)
μ₂_all = zeros(n_sample)
T₁_all = zeros(n_sample, len_λ)
T₂_all = zeros(n_sample, len_λ)
SIF_all = zeros(n_sample, len_λ)
SIF_loadings_all = zeros(n_sample, nSIF)
pseudo_obs_all   = zeros(n_sample, len_λ)

# Generate pseudo data
@threads for i in 1:n_sample
    # ----- rho -----
    ρ_all[i, :] = K₀_recon * coeffs_record[indᵨ[i], :]
    
    # ----- cos(sza) and cos(vza) -----
    μ₁_all[i] = cosd(sza_noSIF[ind_sza[i]])
    μ₂_all[i] = cosd(vza_noSIF[ind_vza[i]])
    
    # ----- Transmittance -----
    σ₁ = @. -log(trans_new[indₜ₁[i], :])
    σ₂ = @. -log(trans_new[indₜ₂[i], :])
    T₁_all[i, :] = @. exp(-σ₁)
    T₂_all[i, :] = @. exp(-σ₁ - σ₂)
    
    # ----- water-leaving SIF -----
    SIF_all[i, :] = SIF_new[:, indₛ[i]]
    SIF_loadings_all[i, :] = SIF_SVD.Loading[1:nSIF, indₛ[i]]
    
    # ----- TOA -----
    pseudo_obs_all[i, :] = @. E / pi * μ₁_all[i] * ρ_all[i, :] * T₂_all[i, :] + SIF_all[i, :] * T₁_all[i, :]

    # ----- noise -----
    stds = sqrt.(c1 .+ c2 .* pseudo_obs_all[i, :])
    pseudo_obs_all[i, :] += randn(len_λ) .* stds
    
    if i % 1000 == 0
        println("Processed $i / $n_sample samples")
    end
end

println("Pseudo observations complete!")



# ===========================================
# Retrieval setup
# ===========================================

#= Principal Components =#

# NMF
HighResNMF = Spectral_NMF(
    trans, 
    bands,
    Float64.(collect(skipmissing(oci_band))); 
    rank=rank
);
# W and H
λ₀ = HighResNMF.band;
W₀ = HighResNMF.Loading;
H₀ = HighResNMF.PrinComp;

# matrics
mean_val  = [round(mean(W₀[:, i]), digits=2) for i in 1:rank];
max_val   = [round(maximum(W₀[:, i]), digits=2) for i in 1:rank];
min_val   = [round(minimum(W₀[:, i]), digits=2) for i in 1:rank];

# s.d. for the loading term
loading_ave_trans = [mean(W₀[:, i]) for i in 1:rank];
loading_var_trans = [var(W₀[:, i]) for i in 1:rank];
println("NMF decomposition complete!")

# SVD 
loading_var_sif   = var(SIF_SVD.Loading[1:nSIF,:], dims=2) .* 2 ;  # 2 as a scale factor?

# MakePriori
diag_values = vcat(
    fill(1e10, n+1),
    loading_var_trans[1:nPC],
    [2.0],
    loading_var_sif[1:nSIF]
);

Sₐ = Diagonal(diag_values);

println("Diagonal terms updated, with diagonal terms: $(diag(Sₐ))")

# remove bands from retrieval evaluation by manually degrading ther SNR
c2_modified    = copy(c2);
wv_degrade_ind = findall((λ .>= λ_remove_min) .& (λ .<= λ_remove_max));
c2_modified[wv_degrade_ind] .= 1e6;  

# Create the retrieval parameters
params = RetrievalParams(
    # Measurement specific
    λ  = oci_band,                   # Wavelength array
    λc = λc, 						 # Centered wavelength
    λ_bl_ind = bl_ind,               # Baseline band indices
    E        = E,                    # Solar irradiance
	c₁       = c1, 					 # PACE SNR 
	c₂       = c2_modified, 			       	
    
    # Forward model settings
    forward_model = forward_model,
    nPoly = n,                       # Degree of Legendre polynomial
    nPC   = nPC,                     # Number of transmittance PCs
    nSIF  = nSIF,                    # Number of SIF PCs
    Sₐ = Sₐ,   					     # Prior covariance
    βₐ = loading_ave_trans,                # Prior state
    PrinComp = HighResNMF.PrinComp',       # Principal components
    SIFComp  = SIF_SVD.PrinComp,           # SIF components
    
    # Iteration settings (optional, have defaults)
    iteration_method = LM_Iteration!,
    nIter = nIter,
    thr_Converge = thr_Converge
)
println("Parameters setup complete!")

# ===========================================
# Retrieval implementation
# ===========================================
Retrieval_all = Vector{Union{Pixel, Missing}}(undef, n_sample)
start_time = now()
println("Starting retrieval at $start_time...")

@threads for i in 1:n_sample
    try
        Retrieval_all[i] = Retrieval_for_Pixel(
            pseudo_obs_all[i, :],
            sza_noSIF[ind_sza[i]],
            vza_noSIF[ind_sza[i]],
            maximum(SIF_all[i, :]),
            1.0,
            1.0,
            params
        )
    catch e
        @warn "Sample $i failed" exception=e
        Retrieval_all[i] = missing
    end
    
    if i % 250 == 0
        println("Retrieved $i / $n_sample")
    end
end

end_time = now()
elapsed_time = end_time - start_time
println("Retrieval complete! Total time elapsed: $elapsed_time")

# save results
version = "v3_1_2"
message = "Configurations: \n" *
          "Wavelength range: $λ_min - $λ_max nm\n" *
          "SNR degradation range: $λ_remove_min - $λ_remove_max nm\n" *
          "Polynomial order: $order\n" *
          "SIF scale factor: $scale_factor_SIF\n" *
          "Number of samples: $n_sample\n" *
          "Random seed: $random_seed\n" *
          "NFLH threshold: $nflh_threshold\n" *
          "NMF rank: $rank\n" *
          "Order of polynomials to fit: $n, Number of retrieval PCs: $nPC, SIF PCs: $nSIF\n" *
          "Number of iterations: $nIter, Convergence threshold: $thr_Converge\n" *
          "Retrieval started at: $start_time, ended at: $end_time, elapsed time: $elapsed_time\n";
ground_truth = (
    psuedo_obs_all = pseudo_obs_all,
    ρ_all = ρ_all,
    T₁_all = T₁_all,
    T₂_all = T₂_all,
    SIF_all = SIF_all,
    SIF_loadings_all = SIF_loadings_all
);
@save "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/retrieval_from_pseudoObs/retrieval_results_$version.jld2" Retrieval_all ground_truth params message