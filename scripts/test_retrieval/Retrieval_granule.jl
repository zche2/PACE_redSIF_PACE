"""
Test the retrieval with granule data
------
Created: 2025-12-03
"""

# Activate project environment
import Pkg
Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")

# Load packages
using JLD2, Interpolations, Revise
using Base.Threads, Dates
using BlockDiagonals, LinearAlgebra
using ForwardDiff, DiffResults, Plots, DelimitedFiles, NCDatasets, Statistics
using Polynomials, Random
using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools
using PACE_SIF

println("Running with $(Threads.nthreads()) threads")

# ===========================================
# path to OCI data
# ===========================================
granule_name = "sample_granule_20250808T204353_new_chl"
path_oci = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/$(granule_name).nc"

nflh_threshold = 0.05  # Threshold for valid NFLH

println("Loading OCI data from $path_oci ...")

# ===========================================
# set parameters
# ==========================================
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

# saving dir
save_dir  = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/retrieval_from_realData/"
timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
save_name = save_dir * "retrieval_for_" * granule_name * "_" * timestamp * ".jld2"
# ===========================================
# load data and set params
# ===========================================

path_transmittance_summer = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc"
path_transmittance_winter = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc"
path_sif_shapes = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/SIF_singular_vector.jld2"
path_snr = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt"


println("Loading data...")

# Load transmittance data
trans, bands = Dataset(path_transmittance_summer) do summer
    Dataset(path_transmittance_winter) do winter
        (cat(summer["transmittance"][:, :], winter["transmittance"][:, :], dims=1),
         summer["band"][:])
    end
end
println("Transmittance: $(size(trans, 1)) spectra")

# Load PACE OCI data
red_band, nflh, vza, sza, chlor_a, E, R_toa = Dataset(path_oci) do oci
    rb = oci["red_wavelength"][:]
    ind = findall(λ_min .< rb .< λ_max)
    (rb, oci["nflh"][:, :], oci["sensor_zenith"][:, :], oci["solar_zenith"][:, :],
     oci["chlor_a"][:, :], oci["red_solar_irradiance"][ind], oci["radiance_red"][:, :, ind])
end
oci_band = red_band[findall(λ_min .< red_band .< λ_max)]
println("PACE data: Band [$λ_min, $λ_max] nm")

# Interpolate transmittance to OCI bands
trans_new = zeros(size(trans, 1), length(oci_band))
@threads for i in 1:size(trans, 1)
    trans_new[i, :] = LinearInterpolation(bands, trans[i, :], extrapolation_bc=0).(oci_band)
end
println("Transmittance interpolated: $(size(trans_new, 1)) spectra")

# Load and process SIF shapes
SIF_shape_dict = JLD2.load(path_sif_shapes)

SIF_SVD = Spectral_SVD(
    SIF_shape_dict["SIF_shapes"]*scale_factor_SIF,
    SIF_shape_dict["SIF_wavelen"],
    Float64.(collect(skipmissing(oci_band))),
    if_log = false
)
println("SIF shapes processed")

# Load SNR coefficients
c1, c2 = let
    data = readdlm(path_snr, String, skipstart=findfirst(==("/end_header"), readlines(path_snr)))
    wvlen = parse.(Float64, data[:, 2])
    idx = findall((data[:, 1] .== "Red") .& (λ_min .< wvlen .< λ_max))
    (parse.(Float64, data[idx, 4]), parse.(Float64, data[idx, 5]))
end

# Print summary
println("Data Summary:")
println("  Wavelength range: $λ_min - $λ_max nm ($(length(oci_band)) bands)")
println("  TOA radiance: $(size(R_toa)) pixels")
println("  Transmittance: $(size(trans_new, 1)) spectra")
println("  SIF shapes: $(size(SIF_SVD.Loading, 2)) components")
println("  SNR coefficients: $(length(c1)) bands")
println()

# ===========================================
# Retrieval set-up
# ===========================================
println("Setting up retrieval parameters...")

λ = oci_band
λc = center_wavelength(oci_band)
bl_ind = map(λ_ref -> argmin(abs.(λ .- λ_ref)), λ_bl_ref)

loading_ave_trans, loading_var_trans, cov_matx, PrinComp = if DecompositionMethod == :NMF
    res = Spectral_NMF(trans, bands, Float64.(collect(skipmissing(oci_band))); 
                       rank=ranks, if_log=if_log)
    ([mean(res.Loading[:, i]) for i in 1:ranks],
     [var(res.Loading[:, i]) for i in 1:ranks],
     cov(res.Loading, dims=1),
     res.PrinComp')
elseif DecompositionMethod == :SVD
    res = Spectral_SVD(Float64.(trans'), bands, 
                       Float64.(collect(skipmissing(oci_band))), if_log=if_log)
    ([mean(res.Loading[i, :]) for i in 1:nPC],
     [var(res.Loading[i, :]) for i in 1:nPC],
     cov(res.Loading[1:nPC, :], dims=2),
     res.PrinComp[:, 1:nPC])
end

println("\t $(DecompositionMethod) decomposition complete!")

# Prior covariance matrix
cov_matx_sif    = cov(SIF_SVD.Loading[1:nSIF, :], dims=2) .* 2;
loading_var_sif = var(SIF_SVD.Loading[1:nSIF, :], dims=2) .* 2
Sₐ = BlockDiagonal([
    diagm(fill(1e10, n+1)),
    cov_matx,
    diagm([2.0]),
    cov_matx_sif
])

# Degrade SNR for specified wavelength range
c2_modified = copy(c2)
c2_modified[findall((λ .>= λ_remove_min) .& (λ .<= λ_remove_max))] .= 1e6

# Create retrieval parameters
params = RetrievalParams(
    λ=oci_band, λc=λc, λ_bl_ind=bl_ind, E=E, c₁=c1, c₂=c2_modified,
    forward_model=(x, px) -> forward_model(x, px; if_log=if_log, return_components=false),
    nPoly=n, nPC=nPC, nSIF=nSIF, Sₐ=Sₐ, βₐ=loading_ave_trans,
    PrinComp=PrinComp, SIFComp=SIF_SVD.PrinComp,
    iteration_method=LM_Iteration!, nIter=nIter, thr_Converge=thr_Converge
)

println("Parameters setup complete!")


# ===========================================
# Run retrieval
# ===========================================
start_time    = now()

valid_indices = findall(coalesce.(nflh .> 0.05, false))

# Preallocate results
results       = Vector{Union{Missing, Pixel}}(missing, length(valid_indices));
# reference value
ref_Rtoa      = Matrix{Union{Missing, AbstractFloat}}(undef, length(valid_indices), length(λ));
ref_nflh      = Vector{Union{Missing, AbstractFloat}}(missing, length(valid_indices));

# start_time
println("Starting retrieval for $(length(valid_indices)) valid pixels...")

# Loop through valid indices only
for (i, idx) in enumerate(valid_indices)
    try
        # Extract data for this pixel
        R_toa_pixel   = R_toa[idx, :]
        sza_pixel     = sza[idx]
        vza_pixel     = vza[idx]
        nflh_pixel    = nflh[idx]
        chlor_a_pixel = chlor_a[idx]

        ref_nflh[i]   = nflh_pixel;
        ref_Rtoa[i,:] = R_toa_pixel;
        
        # Do retrieval
        results[i] = Retrieval_for_Pixel(
            R_toa_pixel, 
            sza_pixel, 
            vza_pixel, 
            nflh_pixel, 
            chlor_a_pixel,  # chlor_a
            1.0,            # flag
            params
        )

        if i % 1000 == 0
            println("Processed $i / $(length(valid_indices)) pixels")
        end
        
    catch e
        @warn "Failed at index $idx: $e"
        results[i] = missing
    end
end

println("Completed $(count(!ismissing, results)) / $(length(valid_indices)) retrievals")
end_time     = now()
elapsed_time = end_time - start_time

println("Retrieval started at: $start_time, ended at: $end_time, elapsed time: $elapsed_time\n")

# ===========================================
# Save results
# ===========================================
println("Saving results to $save_name ...")
@save save_name results ref_nflh ref_Rtoa params configuration
println("Results saved successfully.")
