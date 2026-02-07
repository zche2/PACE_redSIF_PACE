# =========================================================================================================
# This script wraps up parameter set-up for retrieval
# =========================================================================================================


# Load packages
using JLD2, Interpolations, Revise
using Base.Threads, Dates
using BlockDiagonals, LinearAlgebra
using ForwardDiff, DiffResults, Plots, DelimitedFiles, NCDatasets, Statistics
using Polynomials, Random
using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools
using PACE_SIF


function setup_retrieval_parameters(
    red_band,     # red_band from PACE OCI data
    E,            # extraterrestrial solar irradiance at red_band wavelengths
    λ_min, λ_max,           # wavelength range for retrieval
    scale_factor_SIF,       # scaling factor for SIF shapes (to make the variance comparable to transmittance)
    DecompositionMethod,    # :NMF or :SVD
    if_log,       # whether to do log-SVD for transmittance
    n::Int,       # polynomial term
    nPC::Int,     # number of PCs to retain for transmittance
    nSIF::Int,    # number of SIF PCs
    nIter::Int,   # number of iterations for retrieval
    thr_Converge; # convergence threshold for retrieval
    λ_remove_min = 750.0,
    λ_remove_max = 749.0,
    λ_bl_ref = [607.99, 610.36, 612.73, 615.14, 617.6, 620.06, 622.53, 669.52, 
    670.76, 671.99, 673.24, 674.51, 675.73, 676.97, 678.21, 679.45, 
    754.3, 779.33, 867.11, 869.61, 872.13],
    ranks::Int=10,   # rank for NMF decomposition (ignored if SVD is used)
    path_transmittance_summer:: String="/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc", 
    path_transmittance_winter:: String="/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc", 
    path_sif_shapes:: String="/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/SIF_singular_vector.jld2", 
    path_snr:: String="/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt",
)

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
        SIF_shape_dict["SIF_shapes"] * scale_factor_SIF,
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
    println("  Transmittance: $(size(trans_new, 1)) spectra")
    println("  SIF shapes: $(size(SIF_SVD.Loading, 2)) components")
    println("  SNR coefficients: $(length(c1)) bands")
    println()

    # Retrieval set-up
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
    return params

end
