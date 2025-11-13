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
n_sample = 200;
println("Generating $n_sample pseudo observations...")

# Random sampling
Random.seed!(512)
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

# load LUT
println("Load LUT for cross sections...")
o2_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_O2.jld2";
o2_sitp = read_rescale(o2_jld2);
h2o_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_H2O.jld2";
h2o_sitp = read_rescale(h2o_jld2);
metadata = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01.log";
ν_grid_o2, p_grid_hPa, t_grid = o2_sitp.ranges;
println("LUT loaded.")

# generate SRF
println("Generating spectral response function...")

function SRF_for_pace(
        λ_max,
        λ_min,
        ν_step;
        filename = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_RSRs.nc"
    )

    # read data
    pace   = Dataset(filename, "r");
	wavlen = pace["wavelength"][:];
	RSR    = pace["RSR"];
	band   = pace["bands"];

    ind₁   = findall( λ_min .< wavlen .< λ_max );
	ind₂   = findall( λ_min .< band .< λ_max )
	λ_msr  = wavlen[ind₁];
	MyRSR  = RSR[ind₁, ind₂];

    # output high res wavelength 
    ν_min     = λ_to_ν(λ_max);
    ν_max     = λ_to_ν(λ_min);
    ν_grid    = ν_min:ν_step:ν_max;
    wvlen_out = ν_to_λ.(reverse(collect(ν_grid)));
    println("  SRF ν grid: $ν_min to $ν_max cm⁻¹ with step $ν_step cm⁻¹: $ν_grid ");
    println("  SRF λ grid: $(minimum(wvlen_out)) to $(maximum(wvlen_out)) nm with length $(length(wvlen_out))");

    # construct
    MyKernel = KernelInstrument(
		band[ind₂],
		λ_msr,
		MyRSR,
		wvlen_out,
		ν_grid
	);

    return MyKernel
end

MyKernel = SRF_for_pace(λ_max, λ_min, ν_grid_o2.step.hi);

println("SRF generated.")

# parameters
println("Setting up retrieval parameters...")
nPoly  = 6
nLayer = 1
nSIF   = 1

Sₐ   = I(nPoly+nSIF+13) .* 0.;
# update diagonal term
for i=1:(nPoly+1)
    Sₐ[i,i] = 1e10;
    # large variance applies no constrain to these polynomial term
end
# SIF magnitude
Sₐ[end, end] = 1;
# pressure terms
index_pressure = [1, 4, 7, 10] .+ nPoly .+ 1;
Sₐ[CartesianIndex.(index_pressure, index_pressure)] .= 1e4
# temperature terms
index_temperature = [2, 5, 8, 11] .+ nPoly .+ 1;
Sₐ[CartesianIndex.(index_temperature, index_temperature)] .= 1e3;
# vcd for O2 (log10 transformed)
index_vcd_o2 = [3, 9] .+ nPoly .+ 1;
Sₐ[CartesianIndex.(index_vcd_o2, index_vcd_o2)]   .= 0.01;   # variance for O2 statistically is ~0.0025002
# vcd for H2O (log10 transformed)
index_vcd_h2o = [6, 12] .+ nPoly .+ 1;
Sₐ[CartesianIndex.(index_vcd_h2o, index_vcd_h2o)] .= 0.5;    # variance for H2O statistically is ~0.381829

println("Sₐ =====> Diagonal terms updated, with diagonal terms: $(diag(Sₐ))")

# βₐ and γₐ
βₐ = [800.0; 250.0; 24.6; 800.0; 290.0; 22.6];              # see demo for the stats
γₐ = βₐ;                       
println("βₐ and γₐ ======> Set to climatology (?) values")

# iteration method
thr_Converge = 1e-4
MyIter = (px, MyModel) -> LM_Iteration!(
    px, 
    MyModel;
    thr_Converge = thr_Converge,
    nIter = 20,
    γ_init = 1000.0,
    γ⁺ = 10.0,
    γ⁻ = 2.0,
    max_runtime = 400.0
);

params = RetrievalParams_xSecFit(
    λ  = λ,
    λc = λc,
    E  = E,
    c₁ = c1,
    c₂ = c2,
    o2_sitp  = o2_sitp,
    h2o_sitp = h2o_sitp,
    InstrumentKernel = MyKernel,
    forward_model    = forward_model,
    nPoly  = nPoly,
    nLayer = nLayer,
    nSIF   = nSIF,
    Sₐ = Sₐ,
    βₐ = βₐ,
    γₐ = γₐ,
    SIFComp  = SIF_PC,      
    iteration_method = MyIter,
    thr_Converge = thr_Converge
);
println("Retrieval parameters set.\n")

# ===========================================
# Run Retrieval
# ===========================================

Retrieval_all = Vector{Union{Pixel_xSecFit, Missing}}(undef, n_sample);
start_time = now();
println("Starting retrieval for $n_sample samples at $start_time...")

# single pixel retrieval
# i = 27;
# Retrieval_for_Pixel(
#     pseudo_obs_all[i, :],
#     sza_noSIF[ind_sza[i]],
#     vza_noSIF[ind_sza[i]],
#     maximum(SIF_all[i, :]),
#     1.0,
#     1.0,
#     params
# );


@threads for i in 1:n_sample
    if i % (n_sample÷10) == 0
        println("Now at $i / $n_sample")
    end

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
    
    if i % (n_sample÷10) == 0
        println("Retrieved $i / $n_sample")
    end
end

end_time = now()
elapsed_time = end_time - start_time
println("Retrieval complete! Total time elapsed: $elapsed_time")

# ===========================================
# Save results
# ===========================================
version = "v2_4"
message = "nPoly=$nPoly, cutoff at the edge hasn't been fully solved; take viewing geometry into account."
# exclude heavy fields: interpolators and kernels (they're big)
params_to_save = (
    λ = params.λ,
    λc = params.λc,
    E  = params.E,
    c₁ = params.c₁,
    c₂ = params.c₂,
    nPoly = params.nPoly,
    nLayer = params.nLayer,
    nSIF = params.nSIF,
    Sₐ = params.Sₐ,
    βₐ = params.βₐ,
    γₐ = params.γₐ,
    SIFComp = params.SIFComp,
    thr_Converge = params.thr_Converge
)
@save "/home/zhe2/FraLab/PACE_redSIF_PACE/scripts/pseudo_measurement/retrieval_results_$version.jld2" Retrieval_all pseudo_obs_all ρ_all T₁_all T₂_all SIF_all params_to_save message MyIter