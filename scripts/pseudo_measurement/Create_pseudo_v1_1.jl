#!/usr/bin/env julia

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
# Load Data
# ==========================================

println("Loading data...")

# Load MERRA2 transmittance data
println("Loading transmittance data...")
summer = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc")
winter = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc")

trans = cat(summer["transmittance"][:, :], winter["transmittance"][:, :], dims=1)
bands = summer["band"][:]

close(summer)
close(winter)
println("Transmittance data loaded: $(size(trans, 1)) spectra")

# Load PACE OCI data
println("Loading PACE OCI data...")
oci = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/sample_granule_20240830T131442_new_chl.nc")

red_band = oci["red_wavelength"][:]
nflh = oci["nflh"][:, :]
vza = oci["sensor_zenith"][:, :]
sza = oci["solar_zenith"][:, :]
chlor_a = oci["chlor_a"][:, :]

# Select spectral band
λ_min = 620.0  # Define if not already set
λ_max = 860.0  # Define if not already set
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
SIF_shape_dict = JLD2.load("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/SIF_singular_vector.jld2")

# Create interpolator
itp = interpolate(SIF_shape_dict["SIF_shapes"], (BSpline(Linear()), NoInterp()))
range₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end]
range₂ = 1:size(itp, 2)
sitp = scale(itp, range₁, range₂)
setp0 = extrapolate(sitp, 0)

# Interpolate to OCI bands
SIF_new = reduce(hcat, [setp0.(oci_band, i) for i in range₂])
# Scale SIF
scale_factor_SIF = 20
SIF_new *= scale_factor_SIF
println("SIF shapes scaled by factor $scale_factor_SIF: $(size(SIF_new, 2)) shapes")
# interpolation in the first dimension and no interp. in the second
itp₂    = interpolate(SIF_shape_dict["SIF_U"], (BSpline(Linear()), NoInterp()));
# scale
r₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end];
r₂ = 1:size(itp₂, 2);
sitp₂   = scale(itp₂, r₁, r₂);
# set extrapolation filling value = 0
setp0₂  = extrapolate(sitp₂, 0)
# interpolation
SIF_PC  = reduce(hcat, [setp0₂.(oci_band, i) for i in range₂]); 
println("SIF shape interpolated")

# load SNR
println("Loading SNR data...")
filename = raw"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt";
lines = readlines(filename);
end_header_index = findfirst(x -> x == "/end_header", lines);
data  = readdlm(filename, String, skipstart=end_header_index);

FPA   = data[:, 1];                   # 1st column: band
wvlen = parse.(Float64, data[:, 2]);  # 2nd column: center wavelength

wv_val  = (λ_min .< wvlen .< λ_max);
snr_ind = findall((FPA .== "Red") .& wv_val);

# get c1 and c2 at that range
c1    = parse.(Float64, data[snr_ind, 4]);  # 4th column: c1
c2    = parse.(Float64, data[snr_ind, 5]);  # 5th column: c2

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
# =========================================

# Wavelength setup
λ = oci_band
λc = center_wavelength(oci_band)

# Find transmittance baseline
λ_bl_ref = [607.99, 610.36, 612.73, 615.14, 617.6, 620.06, 622.53, 669.52, 
            670.76, 671.99, 673.24, 674.51, 675.73, 676.97, 678.21, 679.45, 
            754.3, 779.33, 867.11, 869.61, 872.13]
bl_ind = map(λ_bl_ref -> argmin(abs.(λ .- λ_bl_ref)), λ_bl_ref)

# Extract valid pixels
valid_mask = findall(coalesce.(nflh .< 0.005, false))
R_noSIF = R_toa[valid_mask, :]
sza_noSIF = sza[valid_mask]
vza_noSIF = vza[valid_mask]

# Fit polynomials
order = 6
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
    
    if i % 100 == 0
        println("Processed $i / $n_pixels")
    end
end

println("Fitting complete!")

# Generate pseudo observations
n_sample = 5000
println("Generating $n_sample pseudo observations...")

# Random sampling
Random.seed!(42)
nᵨ       = n_pixels;   # also select SZA from
nₛ       = size(SIF_new, 2);
nₜ        = size(trans_new, 1);
indᵨ     = rand(1:nᵨ, n_sample);
indₛ      = rand(1:nₛ, n_sample);
indₜ₁     = rand(1:nₜ, n_sample);
indₜ₂     = rand(1:nₜ, n_sample);
ind_sza  = rand(1:nᵨ, n_sample);   # sza
ind_vza  = rand(1:nᵨ, n_sample);   # vza

# Preallocate
len_λ   = length(λ);
ρ_all   = zeros(n_sample, len_λ);
μ₁_all  = zeros(n_sample);
μ₂_all  = zeros(n_sample);
T₁_all  = zeros(n_sample, len_λ);
T₂_all  = zeros(n_sample, len_λ);
SIF_all = zeros(n_sample, len_λ);
pseudo_obs_all = zeros(n_sample, len_λ);

# Generate pseudo data
@threads for i in 1:n_sample
    # ----- rho -----
    ρ_all[i, :] = K₀_recon * coeffs_record[indᵨ[i], :];
    
    # ----- cos(sza) and cos(vza) -----
    μ₁_all[i] = cosd(sza_noSIF[ind_sza[i]]);
    μ₂_all[i] = cosd(vza_noSIF[ind_vza[i]]);
    
    # ----- Transmittance -----
    σ₁ = @. - 1 / μ₁_all[i] * log( trans_new[indₜ₁[i], :] );
    σ₂ = @. - 1 / μ₂_all[i] * log( trans_new[indₜ₂[i], :] );
    T₁_all[i, :] = @. exp( - σ₁ );
    T₂_all[i, :] = @. exp( - σ₁ - σ₂ );
    
    # ----- water-leaving SIF -----
    SIF_all[i, :] = SIF_new[:, indₛ[i]];
    
    # ----- TOA -----
    pseudo_obs_all[i, :] = @. E / pi * μ₁_all[i] * ρ_all[i, :] * T₂_all[i, :] + SIF_all[i, :] * T₁_all[i, :];

    # ----- noise -----
    stds = sqrt.(c1 .+ c2 .* pseudo_obs_all[i, :]);
    pseudo_obs_all[i, :] += randn(len_λ) .* stds;
    
    if i % 1000 == 0
        println("Processed $i / $n_sample samples")
    end
end

println("Pseudo observations complete!")

# ===========================================
# Retrieval setup
# ===========================================

# NMF Decomposition
rank       = 15;
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
loading_ave = [mean(W₀[:, i]) for i in 1:rank];
@show loading_sd  = [var(W₀[:, i]) for i in 1:rank];
println("NMF decomposition complete!")

# prepare parameters
n     = 10;
nPC   = rank;
nSIF  = 1;
Sₐ   = I(n+nPC+nSIF+2) .+ 0.;
# update diagonal term
for i=1:(n+1)
    Sₐ[i,i] = 1e10;
    # large variance applies no constrain to these polynomial term
end
# \beta
for i=(n+2):(n+nPC+1)
    Sₐ[i,i]  = loading_sd[i - (n+1)];
end
# \gamma
Sₐ[n+nPC+2, n+nPC+2] = 2;
# SIF magnitude
Sₐ[end, end] = 1;
println("Diagonal terms updated, with diagonal terms: $(diag(Sₐ))")

# Create the retrieval parameters
params = RetrievalParams(
    # Measurement specific
    λ  = oci_band,                   # Wavelength array
    λc = λc, 						 # Centered wavelength
    λ_bl_ind = bl_ind,               # Baseline band indices
    E        = E,                    # Solar irradiance
	c₁       = c1, 					 # PACE SNR 
	c₂       = c2, 			       	 # PACE SNR
    
    # Forward model settings
    forward_model = forward_model,
    nPoly = n,                       # Degree of Legendre polynomial
    nPC   = nPC,                     # Number of transmittance PCs
    nSIF  = nSIF,                    # Number of SIF PCs
    Sₐ = Sₐ,   					     # Prior covariance
    βₐ = loading_ave,                # Prior state
    PrinComp = HighResNMF.PrinComp', # Principal components
    SIFComp  = SIF_PC,       # SIF components
    
    # Iteration settings (optional, have defaults)
    iteration_method = LM_Iteration!,
    nIter = 25,
    thr_Converge = 1e-6
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
version = "v1_3"
message = "nPoly=$n, rank=$rank\n Change penalty function for T₁ and T₂ conversion: smooth_x = 10. / (1 + exp( -x[px.nPoly+px.nPC+2]) ) + 1."
@save "/home/zhe2/FraLab/PACE_redSIF_PACE/scripts/pseudo_measurement/retrieval_results_$version.jld2" Retrieval_all pseudo_obs_all ρ_all T₁_all T₂_all SIF_all params message