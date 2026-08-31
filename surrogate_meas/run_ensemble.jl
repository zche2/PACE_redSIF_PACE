#!/usr/bin/env julia
# Build 5,000 surrogate measurements + SVD retrieval ensemble.
#
# Usage:
#   julia --project=. -t 8 surrogate_meas/run_ensemble.jl
#   N_SAMPLES=100 julia --project=. -t 4 surrogate_meas/run_ensemble.jl   # smoke test
#   ZERO_SIF=1 REUSE_TRUTH=0 julia --project=. -t 8 surrogate_meas/run_ensemble.jl  # bias floor test
#
# Outputs under surrogate_meas/batch_ensemble/:
#   truth_ensemble.nc, retrieval_ensemble_<tag>.nc, plots_<tag>/*.png
#   With ZERO_SIF=1: truth_ensemble_zeroSIF.nc, *_zeroSIF.nc / plots_*_zeroSIF/

using Base.Threads
using Random
using Statistics
using LinearAlgebra
using TOML
using NCDatasets
using Plots
using JLD2
using ForwardDiff
using SparseArrays
using Dates

const SCRIPT_DIR = @__DIR__
const REPO_ROOT = dirname(SCRIPT_DIR)

# Reuse single-meas helpers (safe: main() only runs when that file is PROGRAM_FILE)
include(joinpath(SCRIPT_DIR, "build_single_meas.jl"))

# ── ensemble settings ─────────────────────────────────────────────────────────
const OUT_DIR = joinpath(SCRIPT_DIR, "batch_ensemble")
const SVD_CONFIG = get(ENV, "SVD_CONFIG", joinpath(SCRIPT_DIR, "configs", "svd_nPC15_npoly5.toml"))
const N_SAMPLES = parse(Int, get(ENV, "N_SAMPLES", "5000"))
const N_ATM_POOL = parse(Int, get(ENV, "N_ATM_POOL", "80"))
const N_CONT_POOL = parse(Int, get(ENV, "N_CONT_POOL", "120"))
const ENSEMBLE_SEED = parse(Int, get(ENV, "ENSEMBLE_SEED", "20260828"))
const REUSE_TRUTH = lowercase(get(ENV, "REUSE_TRUTH", "1")) in ("1", "true", "yes")
const ZERO_SIF = lowercase(get(ENV, "ZERO_SIF", "0")) in ("1", "true", "yes")
const DECAY_RANGE = (1.0, 8.0)          # exponential reflectance-weight decay
const SIF_STRENGTH_RANGE_ENS = ZERO_SIF ? (0.0, 0.0) : (0.0, 0.5)
const RADIANCE_MEAN_RANGE_ENS = (15.0, 30.0)
const N_PLOT_SPECTRA = 40
const GEN_SZA_DEG = 30.0
const GEN_VZA_DEG = 0.0

function _nleg_tag(n_leg::Integer)
    base = "npoly$(n_leg)"
    return ZERO_SIF ? "$(base)_zeroSIF" : base
end

# ── continuum pool from one L1B granule ───────────────────────────────────────

"""Extract up to `n_want` ocean continuum spectra (rescaled) from one L1B file."""
function build_continuum_pool(
    pace_path::AbstractString,
    λ_band::AbstractVector{<:Real};
    n_want::Int=N_CONT_POOL,
    mean_λ_min::Float64,
    mean_λ_max::Float64,
    max_tries::Int=8000,
)
    ds = Dataset(pace_path)
    λ_all = Float64.(ds["red_wavelength"][:])
    E_all = Float64.(ds["red_solar_irradiance"][:])
    n_pix, n_scan, _ = size(ds["radiance_red"])
    has_water = haskey(ds, "watermask")
    mean_ind = findall(mean_λ_min .< λ_all .< mean_λ_max)
    isempty(mean_ind) && error("No bands in mean-radiance window")

    continua = Vector{Vector{Float64}}()
    meta = NamedTuple[]
    tries = 0
    while length(continua) < n_want && tries < max_tries
        tries += 1
        pixel = rand(1:n_pix)
        scan = rand(1:n_scan)
        if has_water
            w = ds["watermask"][pixel, scan]
            (ismissing(w) || Int(w) != 1) && continue
        end
        R_try = Float64.(coalesce.(ds["radiance_red"][pixel, scan, :], NaN))
        sza_try = Float64(coalesce(ds["solar_zenith"][pixel, scan], NaN))
        mean_R = mean(R_try[mean_ind])
        if !(all(isfinite, R_try) && isfinite(sza_try) && sza_try < 85.0 &&
             all(R_try .> 0) && 5.0 <= mean_R <= 80.0)
            continue
        end
        obs_full = (λ=λ_all, E=E_all, R=R_try, sza=sza_try, vza=0.0,
                    pixel=pixel, scan=scan, pace_path=pace_path)
        obs = subset_obs_wavelength(obs_full, REFLECTANCE_λ_MIN, REFLECTANCE_λ_MAX)
        local fit
        try
            fit = fit_legendre_reflectance(
                obs.λ, obs.E, obs.R, obs.sza;
                order=LEGENDRE_ORDER,
                λ_lo=REFLECTANCE_λ_MIN, λ_hi=REFLECTANCE_λ_MAX,
            )
        catch
            continue
        end
        frac_neg = count(<(0), fit.R_fit) / length(fit.R_fit)
        (mean(fit.R_fit) <= 0 || frac_neg >= 0.05) && continue
        # unit-mean continuum in the working window, later rescaled per sample
        ind_m = findall(mean_λ_min .< obs.λ .< mean_λ_max)
        isempty(ind_m) && continue
        R_unit = fit.R_fit ./ mean(fit.R_fit[ind_m])
        R_band = map_spectrum_to_bands(obs.λ, R_unit, λ_band)
        push!(continua, R_band)
        push!(meta, (pixel=pixel, scan=scan, sza=sza_try, mean_R=mean_R))
    end
    close(ds)
    length(continua) >= 10 || error("Continuum pool too small ($(length(continua))); check L1B")
    println("  continuum pool: $(length(continua)) spectra (tries=$tries)")
    return (R=continua, meta=meta)
end

# ── atmospheric τ pool ────────────────────────────────────────────────────────

"""Precompute τ_layer / τ_cum for `n_pool` random MERRA2 profiles."""
function build_atm_pool(
    ctx,
    n_pool::Int=N_ATM_POOL;
)
    ds = Dataset(TRANS_NC)
    n_profiles = ds.dim["profile"]
    idxs = randperm(n_profiles)[1:min(n_pool, n_profiles)]
    pool = Vector{NamedTuple}(undef, length(idxs))
    println("  building atmospheric τ pool ($n_pool profiles)…")
    for (k, ip) in enumerate(idxs)
        prof = load_profile(ds, ip)
        τ_layer, τ_cum = layer_optical_depth(
            ctx.spectral_axis,
            prof.p_full, prof.temp, prof.vcd_dry, prof.vcd_h2o,
            ctx.o2_sitp, ctx.h2o_sitp,
        )
        pool[k] = (
            profile_index=ip,
            τ_layer=τ_layer,
            τ_cum=τ_cum,
            n_layers=size(τ_layer, 1),
            ps_hpa=prof.ps_hpa,
        )
        if k % 10 == 0 || k == length(idxs)
            println("    atm pool $k / $(length(idxs))")
        end
    end
    close(ds)
    return pool
end

"""Build band T1 (atm-only) and T2 (atm×solar) for one pool member + decay."""
function transmittance_from_pool(
    atm,
    decay::Float64,
    T_solar::AbstractVector{<:Real},
    K::AbstractMatrix{<:Real};
    sza_deg::Float64=GEN_SZA_DEG,
    vza_deg::Float64=GEN_VZA_DEG,
)
    amf_down = 1.0 / cosd(sza_deg)
    amf_up = 1.0 / cosd(vza_deg)
    weights = layer_reflectance_weights(atm.n_layers; shape=:exponential, decay=decay)
    τ_2way = weighted_two_way_optical_depth(atm.τ_cum, weights, amf_down, amf_up)
    τ_1way = column_one_way_optical_depth(atm.τ_layer, amf_up)
    T1_atm = exp.(-τ_1way)
    T2_atm = exp.(-τ_2way)
    T2_hres = T2_atm .* T_solar
    T1_band = vec(K * T1_atm)
    T2_band = vec(K * T2_hres)
    return (T1_band=T1_band, T2_band=T2_band, amf_down=amf_down, amf_up=amf_up, weights=weights)
end

# ── SIF library ───────────────────────────────────────────────────────────────

function load_sif_library_matrix(sif_path::AbstractString, λ_dst::AbstractVector{<:Real})
    sif = JLD2.load(MWEF.must_exist(sif_path))
    shapes = Matrix{Float64}(sif["SIF_shapes"])
    λ_ref = Float64.(sif["SIF_wavelen"])
    n_shapes = size(shapes, 2)
    shapes_band = Matrix{Float64}(undef, length(λ_dst), n_shapes)
    for j in 1:n_shapes
        shapes_band[:, j] = map_spectrum_to_bands(λ_ref, shapes[:, j], λ_dst)
    end
    # peak-normalize each column
    for j in 1:n_shapes
        pk = maximum(abs.(shapes_band[:, j]))
        pk > 0 && (shapes_band[:, j] ./= pk)
    end
    return shapes_band
end

# ── generate ensemble ─────────────────────────────────────────────────────────

function generate_ensemble(
    ctx,
    cfg::Dict;
    n_samples::Int=N_SAMPLES,
)
    λ_band = collect(Float64.(ctx.λ))
    n_band = length(λ_band)
    K = ctx.kernel_rsr_out
    λ_min = Float64(get(get(cfg, "spectral", Dict()), "lambda_min_nm", 640.0))
    λ_max = Float64(get(get(cfg, "spectral", Dict()), "lambda_max_nm", 756.0))

    data_cfg = get(cfg, "data", Dict{String, Any}())
    solar_file = get(data_cfg, "solar_file", "solar_merged_20240731_600_33300_100.out")
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)
    T_solar, _ = MWEF.load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
    T_solar = Float64.(T_solar)

    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    pace_file = get(pace_cfg, "pace_file", "sample_granule_20240830T131442_new_chl.nc")
    pace_path = isabspath(pace_file) ? pace_file : joinpath(ctx.paths.base_dir, pace_file)

    println("=== Continuum pool ===")
    cont_pool = build_continuum_pool(
        pace_path, λ_band;
        n_want=N_CONT_POOL, mean_λ_min=λ_min, mean_λ_max=λ_max,
    )

    println("=== Atmospheric pool ===")
    atm_pool = build_atm_pool(ctx, N_ATM_POOL)

    println("=== SIF library ===")
    sif_shapes = load_sif_library_matrix(ctx.paths.sif_path, λ_band)
    n_sif_lib = size(sif_shapes, 2)
    println("  SIF shapes: $n_sif_lib")

    c1 = Float64.(ctx.band_snr_coeffs["c1"])
    c2 = Float64.(ctx.band_snr_coeffs["c2"])
    length(c1) == n_band || error("SNR length mismatch")

    # allocate
    R_toa_clean = Matrix{Float32}(undef, n_band, n_samples)
    R_toa_noisy = Matrix{Float32}(undef, n_band, n_samples)
    R_cont = Matrix{Float32}(undef, n_band, n_samples)
    R_sif_toa = Matrix{Float32}(undef, n_band, n_samples)
    SIF = Matrix{Float32}(undef, n_band, n_samples)
    T1 = Matrix{Float32}(undef, n_band, n_samples)
    T2 = Matrix{Float32}(undef, n_band, n_samples)
    sigma = Matrix{Float32}(undef, n_band, n_samples)

    profile_index = Vector{Int32}(undef, n_samples)
    decay = Vector{Float32}(undef, n_samples)
    sif_library_index = Vector{Int32}(undef, n_samples)
    sif_strength = Vector{Float32}(undef, n_samples)
    cont_pool_index = Vector{Int32}(undef, n_samples)
    mean_radiance_target = Vector{Float32}(undef, n_samples)

    println("=== Generating $n_samples samples ===")
    for i in 1:n_samples
        ia = rand(1:length(atm_pool))
        ic = rand(1:length(cont_pool.R))
        isif = rand(1:n_sif_lib)
        dec = DECAY_RANGE[1] + (DECAY_RANGE[2] - DECAY_RANGE[1]) * rand()
        strength = SIF_STRENGTH_RANGE_ENS[1] +
                   (SIF_STRENGTH_RANGE_ENS[2] - SIF_STRENGTH_RANGE_ENS[1]) * rand()
        if ZERO_SIF
            strength = 0.0
        end
        target = RADIANCE_MEAN_RANGE_ENS[1] +
                 (RADIANCE_MEAN_RANGE_ENS[2] - RADIANCE_MEAN_RANGE_ENS[1]) * rand()

        tr = transmittance_from_pool(atm_pool[ia], dec, T_solar, K)
        R_c = cont_pool.R[ic] .* target
        SIF_i = sif_shapes[:, isif] .* strength
        R_sif = SIF_i .* tr.T1_band
        R_bg = R_c .* tr.T2_band
        R_clean = R_bg .+ R_sif
        σ = noise_std_from_snr(R_clean, c1, c2)
        R_noisy = R_clean .+ randn(n_band) .* σ

        R_toa_clean[:, i] .= Float32.(R_clean)
        R_toa_noisy[:, i] .= Float32.(R_noisy)
        R_cont[:, i] .= Float32.(R_c)
        R_sif_toa[:, i] .= Float32.(R_sif)
        SIF[:, i] .= Float32.(SIF_i)
        T1[:, i] .= Float32.(tr.T1_band)
        T2[:, i] .= Float32.(tr.T2_band)
        sigma[:, i] .= Float32.(σ)

        profile_index[i] = Int32(atm_pool[ia].profile_index)
        decay[i] = Float32(dec)
        sif_library_index[i] = Int32(isif)
        sif_strength[i] = Float32(strength)
        cont_pool_index[i] = Int32(ic)
        mean_radiance_target[i] = Float32(target)

        if i % 500 == 0 || i == n_samples
            println("  generated $i / $n_samples")
        end
    end

    return (
        λ=λ_band,
        R_toa_clean=R_toa_clean,
        R_toa_noisy=R_toa_noisy,
        R_cont=R_cont,
        R_sif_toa=R_sif_toa,
        SIF=SIF,
        T1=T1,
        T2=T2,
        sigma=sigma,
        profile_index=profile_index,
        decay=decay,
        sif_library_index=sif_library_index,
        sif_strength=sif_strength,
        cont_pool_index=cont_pool_index,
        mean_radiance_target=mean_radiance_target,
        pace_path=pace_path,
        sif_path=String(ctx.paths.sif_path),
    )
end

function write_truth_nc(ens, path::AbstractString)
    n_band, n_samp = size(ens.R_toa_noisy)
    isfile(path) && rm(path)
    ds = NCDataset(path, "c")
    defDim(ds, "band", n_band)
    defDim(ds, "sample", n_samp)
    defVar(ds, "wavelength", Float64.(ens.λ), ("band",); attrib=Dict("units"=>"nm"))
    for (name, arr, units) in (
        ("R_toa_clean", ens.R_toa_clean, "W m-2 sr-1 um-1"),
        ("R_toa_noisy", ens.R_toa_noisy, "W m-2 sr-1 um-1"),
        ("R_cont", ens.R_cont, "W m-2 sr-1 um-1"),
        ("R_sif_toa", ens.R_sif_toa, "W m-2 sr-1 um-1"),
        ("SIF", ens.SIF, "W m-2 sr-1 um-1"),
        ("T1", ens.T1, "1"),
        ("T2", ens.T2, "1"),
        ("sigma_noise", ens.sigma, "W m-2 sr-1 um-1"),
    )
        defVar(ds, name, arr, ("band", "sample"); attrib=Dict("units"=>units))
    end
    defVar(ds, "profile_index", ens.profile_index, ("sample",))
    defVar(ds, "decay", ens.decay, ("sample",); attrib=Dict("long_name"=>"exponential reflectance-weight decay"))
    defVar(ds, "sif_library_index", ens.sif_library_index, ("sample",))
    defVar(ds, "sif_strength", ens.sif_strength, ("sample",); attrib=Dict("long_name"=>"peak SIF radiance"))
    defVar(ds, "cont_pool_index", ens.cont_pool_index, ("sample",))
    defVar(ds, "mean_radiance_target", ens.mean_radiance_target, ("sample",))
    ds.attrib["title"] = "Surrogate measurement ensemble"
    ds.attrib["created"] = string(Dates.now())
    ds.attrib["n_samples"] = n_samp
    ds.attrib["gen_sza_deg"] = GEN_SZA_DEG
    ds.attrib["gen_vza_deg"] = GEN_VZA_DEG
    ds.attrib["sif_path"] = ens.sif_path
    ds.attrib["pace_path"] = ens.pace_path
    ds.attrib["ensemble_seed"] = ENSEMBLE_SEED
    ds.attrib["zero_sif"] = Int(ZERO_SIF)
    close(ds)
    println("Wrote truth NetCDF: ", path)
end

function load_truth_nc(path::AbstractString)
    isfile(path) || error("Missing truth NetCDF: $path")
    ds = Dataset(path)
    ens = (
        λ=Float64.(ds["wavelength"][:]),
        R_toa_clean=Float32.(ds["R_toa_clean"][:, :]),
        R_toa_noisy=Float32.(ds["R_toa_noisy"][:, :]),
        R_cont=Float32.(ds["R_cont"][:, :]),
        R_sif_toa=Float32.(ds["R_sif_toa"][:, :]),
        SIF=Float32.(ds["SIF"][:, :]),
        T1=Float32.(ds["T1"][:, :]),
        T2=Float32.(ds["T2"][:, :]),
        sigma=Float32.(ds["sigma_noise"][:, :]),
        profile_index=Int32.(ds["profile_index"][:]),
        decay=Float32.(ds["decay"][:]),
        sif_library_index=Int32.(ds["sif_library_index"][:]),
        sif_strength=Float32.(ds["sif_strength"][:]),
        cont_pool_index=Int32.(ds["cont_pool_index"][:]),
        mean_radiance_target=Float32.(ds["mean_radiance_target"][:]),
        pace_path=String(ds.attrib["pace_path"]),
        sif_path=String(ds.attrib["sif_path"]),
    )
    close(ds)
    println("Loaded truth NetCDF: ", path, "  (n_samples=$(size(ens.R_toa_noisy, 2)))")
    return ens
end

# ── plots ─────────────────────────────────────────────────────────────────────

function plot_ensemble_diagnostics(ens; out_dir::AbstractString)
    mkpath(out_dir)
    λ = ens.λ
    n = size(ens.R_toa_noisy, 2)
    idx = sort(randperm(n)[1:min(N_PLOT_SPECTRA, n)])

    p1 = plot(size=(1100, 450), legend=false, xlabel="Wavelength [nm]", ylabel="Radiance",
              title="Noisy TOA radiance ($N_PLOT_SPECTRA / $n samples)")
    for j in idx
        plot!(p1, λ, ens.R_toa_noisy[:, j]; lw=0.8, alpha=0.45)
    end
    savefig(p1, joinpath(out_dir, "radiance_samples.png"))

    p2 = plot(size=(1100, 450), legend=false, xlabel="Wavelength [nm]", ylabel="SIF radiance",
              title="Surface SIF shapes (scaled; $N_PLOT_SPECTRA samples)")
    for j in idx
        plot!(p2, λ, ens.SIF[:, j]; lw=0.9, alpha=0.5)
    end
    savefig(p2, joinpath(out_dir, "sif_samples.png"))

    p3 = plot(size=(1100, 450), legend=false, xlabel="Wavelength [nm]", ylabel="SIF×T₁",
              title="SIF at TOA ($N_PLOT_SPECTRA samples)")
    for j in idx
        plot!(p3, λ, ens.R_sif_toa[:, j]; lw=0.9, alpha=0.5)
    end
    savefig(p3, joinpath(out_dir, "sif_toa_samples.png"))

    p4 = plot(
        histogram(ens.decay; bins=30, label="decay", xlabel="Decay", ylabel="Count",
                  title="Exponential decay distribution"),
        histogram(ens.sif_strength; bins=30, label="SIF strength", xlabel="Peak SIF",
                  ylabel="Count", title="SIF strength distribution"),
        histogram(ens.mean_radiance_target; bins=30, label="mean R", xlabel="Mean continuum",
                  ylabel="Count", title="Continuum mean target"),
        layout=(1, 3), size=(1200, 350), legend=false,
    )
    savefig(p4, joinpath(out_dir, "parameter_histograms.png"))

    println("Saved diagnostic plots under ", out_dir)
end

# ── retrieval ─────────────────────────────────────────────────────────────────

function run_ensemble_retrieval(ens, svd_cfg::Dict; out_nc::AbstractString)
    λ = collect(Float64.(ens.λ))
    n_band, n_samp = size(ens.R_toa_noisy)
    println("=== SVD retrieval setup ===")
    sh = prepare_svd_retrieval_setup(svd_cfg, λ)
    println("  n_pc=$(sh.n_pc)  n_legendre=$(sh.n_leg)  n_state=$(sh.layout.n_state)  sif_nev=$(size(sh.sif_basis, 2))")
    sh.n_pc == 15 || @warn "Expected n_pc=15, got $(sh.n_pc)"
    println("  using n_legendre=$(sh.n_leg)")

    solar = load_l1b_solar_on_bands(ens.pace_path, λ)
    println("  solar from $(basename(solar.pace_path)), esd=$(round(solar.esd, digits=4)), SZA=$(PSEUDO_RETRIEVAL_SZA)°")

    n_ev = size(sh.sif_basis, 2)
    n_state = sh.layout.n_state
    i678 = argmin(abs.(λ .- 678.2))

    R_fit = Matrix{Float32}(undef, n_band, n_samp)
    sif_toa_ret = Matrix{Float32}(undef, n_band, n_samp)
    resid = Matrix{Float32}(undef, n_band, n_samp)
    sif_coeff_ret = Matrix{Float32}(undef, n_ev, n_samp)
    state = Matrix{Float32}(undef, n_state, n_samp)
    status = fill(Int16(-1), n_samp)
    converged = fill(UInt8(0), n_samp)
    n_steps = fill(Int16(0), n_samp)
    rmse = fill(Float32(NaN), n_samp)
    reduced_chi2 = fill(Float32(NaN), n_samp)
    dof = fill(Float32(NaN), n_samp)
    sif_678_true = Vector{Float32}(undef, n_samp)
    sif_678_ret = fill(Float32(NaN), n_samp)
    sif_toa_mean_true = Vector{Float32}(undef, n_samp)
    sif_toa_mean_ret = fill(Float32(NaN), n_samp)

    println("=== Retrieving $n_samp samples ($(nthreads()) threads) ===")
    t0 = time()
    @threads for i in 1:n_samp
        y = Float64.(ens.R_toa_noisy[:, i])
        try
            ret = run_pseudo_svd_retrieval(y, λ, sh, solar; sza_deg=PSEUDO_RETRIEVAL_SZA)
            R_fit[:, i] .= Float32.(ret.y_fit)
            sif_toa_ret[:, i] .= Float32.(ret.sif_toa_fit)
            resid[:, i] .= Float32.(ret.resid)
            sif_coeff_ret[:, i] .= Float32.(ret.sif_coeff)
            state[:, i] .= Float32.(ret.x)
            status[i] = Int16(ret.stats.status)
            converged[i] = ret.stats.converged ? UInt8(1) : UInt8(0)
            n_steps[i] = Int16(ret.stats.n_steps)
            rmse[i] = Float32(ret.stats.rmse)
            reduced_chi2[i] = Float32(ret.stats.reduced_chi2)
            dof[i] = Float32(ret.stats.dof)
            sif_678_ret[i] = Float32(ret.sif_toa_fit[i678])
            sif_toa_mean_ret[i] = Float32(mean(ret.sif_toa_fit))
        catch e
            status[i] = Int16(4)
            @warn "Retrieval failed" sample=i exception=e
        end
        sif_678_true[i] = ens.R_sif_toa[i678, i]
        sif_toa_mean_true[i] = Float32(mean(ens.R_sif_toa[:, i]))
        if threadid() == 1 && (i % 200 == 0)
            # approximate progress (not exact under threads)
        end
    end
    # sequential progress-friendly recount
    n_ok = count(==(Int16(1)), status)
    println("  done in $(round(time()-t0, digits=1)) s; status=1 (converged): $n_ok / $n_samp")

    isfile(out_nc) && rm(out_nc)
    ds = NCDataset(out_nc, "c")
    defDim(ds, "band", n_band)
    defDim(ds, "sample", n_samp)
    defDim(ds, "sif_ev", n_ev)
    defDim(ds, "state", n_state)
    defVar(ds, "wavelength", Float64.(λ), ("band",); attrib=Dict("units"=>"nm"))
    defVar(ds, "R_fit", R_fit, ("band", "sample"))
    defVar(ds, "sif_toa_ret", sif_toa_ret, ("band", "sample"))
    defVar(ds, "sif_toa_true", Float32.(ens.R_sif_toa), ("band", "sample"))
    defVar(ds, "resid", resid, ("band", "sample"))
    defVar(ds, "sif_coeff_ret", sif_coeff_ret, ("sif_ev", "sample"))
    defVar(ds, "state", state, ("state", "sample"))
    defVar(ds, "status", status, ("sample",))
    defVar(ds, "converged", converged, ("sample",))
    defVar(ds, "n_steps", n_steps, ("sample",))
    defVar(ds, "rmse", rmse, ("sample",))
    defVar(ds, "reduced_chi2", reduced_chi2, ("sample",))
    defVar(ds, "dof", dof, ("sample",))
    defVar(ds, "sif_678_true", sif_678_true, ("sample",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "sif_678_ret", sif_678_ret, ("sample",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "sif_toa_mean_true", sif_toa_mean_true, ("sample",))
    defVar(ds, "sif_toa_mean_ret", sif_toa_mean_ret, ("sample",))
    defVar(ds, "decay", ens.decay, ("sample",))
    defVar(ds, "sif_strength", ens.sif_strength, ("sample",))
    defVar(ds, "profile_index", ens.profile_index, ("sample",))
    ds.attrib["title"] = "Surrogate ensemble SVD retrievals"
    ds.attrib["created"] = string(Dates.now())
    ds.attrib["n_pc"] = sh.n_pc
    ds.attrib["n_legendre"] = sh.n_leg
    ds.attrib["retrieval_sza_deg"] = PSEUDO_RETRIEVAL_SZA
    ds.attrib["svd_config"] = SVD_CONFIG
    ds.attrib["zero_sif"] = Int(ZERO_SIF)
    close(ds)
    println("Wrote retrieval NetCDF: ", out_nc)

    return (
        sif_678_true=sif_678_true,
        sif_678_ret=sif_678_ret,
        sif_toa_mean_true=sif_toa_mean_true,
        sif_toa_mean_ret=sif_toa_mean_ret,
        sif_toa_ret=sif_toa_ret,
        R_fit=R_fit,
        resid=resid,
        state=state,
        status=status,
        rmse=rmse,
        reduced_chi2=reduced_chi2,
        i678=i678,
        n_pc=sh.n_pc,
        n_leg=sh.n_leg,
        PCs=sh.PCs,
        layout=sh.layout,
        log_trans=sh.log_trans,
    )
end

"""α_coeff = 10/(1+e^{-α_raw}) + 1  (same as make_svd_forward_model_λ)."""
alpha_coeff_from_raw(α_raw::Real) = 10.0 / (1.0 + exp(-Float64(α_raw))) + 1.0

"""
Reconstruct T1 (trans_up) and T2 (trans_updown) from SVD state.
  T1 = exp(PCs·c)           [log_transform=true]
  T2 = T1^α_coeff           [α_coeff from α_raw]
"""
function reconstruct_T_from_state(
    state::AbstractMatrix{<:Real},
    PCs::AbstractMatrix{<:Real},
    layout;
    log_transform::Bool=true,
)
    n_state = layout.n_state
    if size(state, 1) == n_state
        X = state
    elseif size(state, 2) == n_state
        X = collect(state')
    else
        error("state size $(size(state)) incompatible with n_state=$n_state")
    end
    n_pc = layout.n_pc
    n_λ = size(PCs, 1)
    n_samp = size(X, 2)
    T1 = Matrix{Float64}(undef, n_λ, n_samp)
    T2 = Matrix{Float64}(undef, n_λ, n_samp)
    α_coeff = Vector{Float64}(undef, n_samp)
    PC = Float64.(PCs[:, 1:n_pc])
    for i in 1:n_samp
        c = @view X[layout.idx_pc, i]
        α_raw = X[first(layout.idx_alpha), i]
        α = alpha_coeff_from_raw(α_raw)
        α_coeff[i] = α
        t1 = log_transform ? exp.(PC * c) : 1.0 .+ PC * c
        t1 = max.(t1, eps(Float64))
        T1[:, i] .= t1
        T2[:, i] .= exp.(α .* log.(t1))
    end
    return (T1=T1, T2=T2, α_coeff=α_coeff)
end

function plot_radiance_residuals(ens, ret, ok, λ, out_dir)
    R = Float64.(ret.resid[:, ok])
    μ = vec(mean(R; dims=2))
    σ = vec(std(R; dims=2))
    p_res = plot(λ, μ; ribbon=σ, label="mean ± 1σ", color=:navy, lw=2,
                 xlabel="Wavelength [nm]", ylabel="obs − fit",
                 title="Radiance residual spectrum", size=(1000, 400), legend=:outerright)
    hline!(p_res, [0.0]; color=:black, ls=:dash, label="")
    savefig(p_res, joinpath(out_dir, "residual_spectrum.png"))

    # residual vs instrument noise for a few examples
    ex = ok[1:min(6, length(ok))]
    plots_r = []
    for i in ex
        σ_i = Float64.(ens.sigma[:, i])
        pk = plot(λ, ret.resid[:, i]; label="obs − fit", color=:navy, lw=1.4,
                  title="sample $i", xlabel="λ [nm]", ylabel="Residual",
                  legend=:outerright, size=(500, 280))
        plot!(pk, λ, σ_i; label="+σ", color=:gray, lw=1.0, ls=:dash)
        plot!(pk, λ, -σ_i; label="-σ", color=:gray, lw=1.0, ls=:dash)
        push!(plots_r, pk)
    end
    savefig(plot(plots_r...; layout=(2, 3), size=(1400, 700)),
            joinpath(out_dir, "residual_examples.png"))

    # normalized residual mean ± σ
    Rn = similar(R)
    for (j, i) in enumerate(ok)
        σ_i = max.(Float64.(ens.sigma[:, i]), eps(Float64))
        Rn[:, j] .= R[:, j] ./ σ_i
    end
    μn = vec(mean(Rn; dims=2))
    σn = vec(std(Rn; dims=2))
    p_n = plot(λ, μn; ribbon=σn, label="mean ± 1σ", color=:darkorange, lw=2,
               xlabel="Wavelength [nm]", ylabel="(obs − fit) / σ",
               title="Noise-normalized residual spectrum", size=(1000, 400), legend=:outerright)
    hline!(p_n, [0.0]; color=:black, ls=:dash, label="")
    savefig(p_n, joinpath(out_dir, "residual_normalized_spectrum.png"))
end

function plot_T_comparison(ens, ret, ok, λ, out_dir)
    (!hasproperty(ret, :state) || !hasproperty(ret, :PCs)) && (@warn "No state/PCs in ret; skip T₁ plots"; return)
    recon = reconstruct_T_from_state(
        Float64.(ret.state), ret.PCs, ret.layout; log_transform=ret.log_trans,
    )
    T1t = Float64.(ens.T1[:, ok])
    T1f = recon.T1[:, ok]
    T2t = Float64.(ens.T2[:, ok])
    T2f = recon.T2[:, ok]
    α = recon.α_coeff[ok]

    m_t = vec(mean(T1t; dims=1))
    m_f = vec(mean(T1f; dims=1))
    A = hcat(ones(length(m_t)), m_t)
    β = A \ m_f
    r2 = 1 - sum((m_f .- A * β) .^ 2) / max(sum((m_f .- mean(m_f)) .^ 2), eps())
    bias = mean(m_f .- m_t)
    rmse = sqrt(mean((m_f .- m_t) .^ 2))
    lims = extrema(vcat(m_t, m_f))
    pad = 0.05 * (lims[2] - lims[1] + eps())
    lims = (lims[1] - pad, lims[2] + pad)

    p1 = scatter(m_t, m_f; ms=2, alpha=0.3, label="n=$(length(ok))",
                 xlabel="True mean T₁", ylabel="Fitted mean T₁ (from PCs)",
                 title="Mean T₁: bias=$(round(bias, digits=4)), RMSE=$(round(rmse, digits=4)), R²=$(round(r2, digits=3))",
                 size=(700, 650), legend=:topleft)
    plot!(p1, [lims[1], lims[2]], [lims[1], lims[2]]; color=:black, ls=:dash, label="1:1")
    xx = range(lims[1], lims[2]; length=50)
    plot!(p1, xx, β[1] .+ β[2] .* xx; color=:crimson, lw=2,
          label="fit: y=$(round(β[1], digits=3))+$(round(β[2], digits=3))x")
    xlims!(p1, lims); ylims!(p1, lims)
    savefig(p1, joinpath(out_dir, "T1_mean_scatter.png"))

    i750 = argmin(abs.(λ .- 750.0))
    i688 = argmin(abs.(λ .- 688.0))
    function _scatter_band(ib, name)
        tt = vec(T1t[ib, :]); ff = vec(T1f[ib, :])
        lim = extrema(vcat(tt, ff))
        p = scatter(tt, ff; ms=1.5, alpha=0.25, label="",
                    xlabel="True T₁", ylabel="Fitted T₁",
                    title="T₁ @ $(round(λ[ib], digits=1)) nm ($name)", size=(550, 520))
        plot!(p, [lim[1], lim[2]], [lim[1], lim[2]]; color=:black, ls=:dash, label="1:1")
        return p
    end
    savefig(plot(_scatter_band(i750, "continuum"), _scatter_band(i688, "O₂-B");
                 layout=(1, 2), size=(1100, 500)),
            joinpath(out_dir, "T1_band_scatter.png"))

    ex = ok[1:min(6, length(ok))]
    plots_t1 = []
    for i in ex
        pk = plot(λ, ens.T1[:, i]; label="True T₁", color=:navy, lw=2,
                  title="sample $i, α=$(round(recon.α_coeff[i], digits=2))",
                  xlabel="λ [nm]", ylabel="T₁", legend=:outerright, size=(500, 280))
        plot!(pk, λ, recon.T1[:, i]; label="Fitted T₁ (PCs)", color=:crimson, lw=1.5, ls=:dash)
        push!(plots_t1, pk)
    end
    savefig(plot(plots_t1...; layout=(2, 3), size=(1400, 700)),
            joinpath(out_dir, "T1_spectra_examples.png"))

    plots_t2 = []
    for i in ex
        pk = plot(λ, ens.T2[:, i]; label="True T₂ (atm×solar)", color=:darkred, lw=2,
                  title="sample $i", xlabel="λ [nm]", ylabel="T₂", legend=:outerright, size=(500, 280))
        plot!(pk, λ, recon.T2[:, i]; label="Fitted T₂ = T₁^α", color=:orange, lw=1.5, ls=:dash)
        push!(plots_t2, pk)
    end
    savefig(plot(plots_t2...; layout=(2, 3), size=(1400, 700)),
            joinpath(out_dir, "T2_spectra_examples.png"))

    savefig(histogram(α; bins=40, label="", xlabel="α_coeff = 10/(1+e^{-α_raw})+1",
                      ylabel="Count", title="Retrieved α coefficient", size=(700, 400)),
            joinpath(out_dir, "alpha_coeff_hist.png"))

    dT1 = T1f .- T1t
    μ = vec(mean(dT1; dims=2))
    σ = vec(std(dT1; dims=2))
    p_tres = plot(λ, μ; ribbon=σ, label="mean ± 1σ", color=:navy, lw=2,
                  xlabel="Wavelength [nm]", ylabel="Fitted − true T₁",
                  title="T₁ residual spectrum", size=(1000, 400), legend=:outerright)
    hline!(p_tres, [0.0]; color=:black, ls=:dash, label="")
    savefig(p_tres, joinpath(out_dir, "T1_residual_spectrum.png"))

    open(joinpath(out_dir, "T1_comparison_summary.txt"), "w") do io
        println(io, "n_converged\t$(length(ok))")
        println(io, "n_pc\t$(ret.n_pc)")
        println(io, "n_legendre\t$(ret.n_leg)")
        println(io, "log_transform\t$(ret.log_trans)")
        println(io, "mean_T1_bias\t$bias")
        println(io, "mean_T1_rmse\t$rmse")
        println(io, "mean_T1_r2\t$r2")
        println(io, "mean_T1_slope\t$(β[2])")
        println(io, "mean_T1_intercept\t$(β[1])")
        println(io, "median_alpha_coeff\t$(median(α))")
        println(io, "mean_alpha_coeff\t$(mean(α))")
    end
    println("T₁ comparison — bias=$(round(bias, digits=4)), RMSE=$(round(rmse, digits=4)), R²=$(round(r2, digits=3)), slope=$(round(β[2], digits=3))")
    println("α_coeff — mean=$(round(mean(α), digits=3)), median=$(round(median(α), digits=3))")
end

function plot_retrieval_comparison(ens, ret; out_dir::AbstractString)
    mkpath(out_dir)
    ok = findall(ret.status .== 1)
    isempty(ok) && (@warn "No converged retrievals to plot"; return)

    t = Float64.(ret.sif_678_true[ok])
    r = Float64.(ret.sif_678_ret[ok])
    # linear fit
    A = hcat(ones(length(t)), t)
    β = A \ r
    r2 = 1 - sum((r .- A * β) .^ 2) / sum((r .- mean(r)) .^ 2)
    bias = mean(r .- t)
    rmse_s = sqrt(mean((r .- t) .^ 2))

    lims = extrema(vcat(t, r))
    pad = 0.05 * (lims[2] - lims[1] + eps())
    lims = (lims[1] - pad, lims[2] + pad)
    p_sc = scatter(t, r; ms=2, alpha=0.35, label="samples (n=$(length(ok)))",
                   xlabel="True SIF×T₁ @ 678 nm", ylabel="Retrieved SIF×T₁ @ 678 nm",
                   title="SIF@678: bias=$(round(bias, digits=3)), RMSE=$(round(rmse_s, digits=3)), R²=$(round(r2, digits=3))",
                   size=(700, 650), legend=:topleft)
    plot!(p_sc, [lims[1], lims[2]], [lims[1], lims[2]]; color=:black, ls=:dash, label="1:1")
    xx = range(lims[1], lims[2]; length=50)
    plot!(p_sc, xx, β[1] .+ β[2] .* xx; color=:crimson, lw=2,
          label="fit: y=$(round(β[1], digits=3))+$(round(β[2], digits=3))x")
    xlims!(p_sc, lims); ylims!(p_sc, lims)
    savefig(p_sc, joinpath(out_dir, "sif678_scatter.png"))

    tm = Float64.(ret.sif_toa_mean_true[ok])
    rm = Float64.(ret.sif_toa_mean_ret[ok])
    p_m = scatter(tm, rm; ms=2, alpha=0.35, label="",
                  xlabel="True mean SIF×T₁", ylabel="Retrieved mean SIF×T₁",
                  title="Band-mean SIF×T₁", size=(650, 600))
    lim2 = extrema(vcat(tm, rm))
    plot!(p_m, [lim2[1], lim2[2]], [lim2[1], lim2[2]]; color=:black, ls=:dash, label="1:1")
    savefig(p_m, joinpath(out_dir, "sif_mean_scatter.png"))

    p_h = plot(
        histogram(ret.rmse[ok]; bins=40, label="RMSE", title="Retrieval RMSE", xlabel="RMSE"),
        histogram(ret.reduced_chi2[ok]; bins=40, label="rχ²", title="Reduced χ²", xlabel="rχ²"),
        histogram(Float64.(ens.decay[ok]), Float64.(ret.sif_678_ret[ok] .- ret.sif_678_true[ok]);
                  bins=30, xlabel="Decay", ylabel="SIF@678 error (ret−true)",
                  title="SIF error vs decay"),
        layout=(1, 3), size=(1200, 380), legend=false,
    )
    savefig(p_h, joinpath(out_dir, "retrieval_diagnostics.png"))

    # example spectra
    ex = ok[1:min(6, length(ok))]
    λ = ens.λ
    plots_ex = []
    for (k, i) in enumerate(ex)
        pk = plot(λ, ens.R_toa_noisy[:, i]; label="noisy obs", color=:gray, lw=1.2,
                  title="sample $i", xlabel="λ [nm]", ylabel="Radiance", legend=:outerright)
        plot!(pk, λ, ens.R_toa_clean[:, i]; label="true clean", color=:black, lw=1.5)
        plot!(pk, λ, ret.R_fit[:, i]; label="fit", color=:crimson, lw=1.2, ls=:dash)
        push!(plots_ex, pk)
    end
    p_ex = plot(plots_ex...; layout=(2, 3), size=(1400, 700))
    savefig(p_ex, joinpath(out_dir, "toa_fit_examples.png"))

    plots_s = []
    for (k, i) in enumerate(ex)
        pk = plot(λ, ens.R_sif_toa[:, i]; label="true SIF×T₁", color=:forestgreen, lw=2,
                  title="sample $i", xlabel="λ [nm]", ylabel="SIF×T₁", legend=:outerright)
        plot!(pk, λ, ret.sif_toa_ret[:, i]; label="retrieved", color=:darkgreen, lw=1.5, ls=:dash)
        push!(plots_s, pk)
    end
    p_sx = plot(plots_s...; layout=(2, 3), size=(1400, 700))
    savefig(p_sx, joinpath(out_dir, "sif_fit_examples.png"))

    # radiance residuals + reconstructed T₁/T₂ (same content as plot_t1_comparison.jl)
    plot_radiance_residuals(ens, ret, ok, λ, out_dir)
    plot_T_comparison(ens, ret, ok, λ, out_dir)

    open(joinpath(out_dir, "retrieval_summary.txt"), "w") do io
        println(io, "n_samples\t$(size(ens.R_toa_noisy, 2))")
        println(io, "n_converged\t$(length(ok))")
        println(io, "n_pc\t$(ret.n_pc)")
        println(io, "n_legendre\t$(ret.n_leg)")
        println(io, "sif678_bias\t$bias")
        println(io, "sif678_rmse\t$rmse_s")
        println(io, "sif678_r2\t$r2")
        println(io, "sif678_slope\t$(β[2])")
        println(io, "sif678_intercept\t$(β[1])")
        println(io, "median_rmse\t$(median(Float64.(ret.rmse[ok])))")
        println(io, "median_rchi2\t$(median(Float64.(ret.reduced_chi2[ok])))")
    end
    println("Saved retrieval comparison plots under ", out_dir)
    println("SIF@678 — bias=$(round(bias, digits=4)), RMSE=$(round(rmse_s, digits=4)), R²=$(round(r2, digits=3)), slope=$(round(β[2], digits=3))")
end

"""Load retrieval NetCDF + rebuild SVD setup so plots can be regenerated without refitting."""
function load_retrieval_for_plots(ret_nc::AbstractString, svd_cfg::Dict, λ::AbstractVector{<:Real})
    ds = Dataset(ret_nc)
    n_state = Int(ds.dim["state"])
    state = Matrix{Float32}(ds["state"][:, :])
    if size(state, 1) != n_state && size(state, 2) == n_state
        state = collect(state')
    end
    ret = (
        sif_678_true=Float32.(ds["sif_678_true"][:]),
        sif_678_ret=Float32.(ds["sif_678_ret"][:]),
        sif_toa_mean_true=Float32.(ds["sif_toa_mean_true"][:]),
        sif_toa_mean_ret=Float32.(ds["sif_toa_mean_ret"][:]),
        sif_toa_ret=Float32.(ds["sif_toa_ret"][:, :]),
        R_fit=Float32.(ds["R_fit"][:, :]),
        resid=Float32.(ds["resid"][:, :]),
        state=state,
        status=Int16.(ds["status"][:]),
        rmse=Float32.(ds["rmse"][:]),
        reduced_chi2=Float32.(ds["reduced_chi2"][:]),
        n_pc=Int(get(ds.attrib, "n_pc", -1)),
        n_leg=Int(get(ds.attrib, "n_legendre", -1)),
    )
    close(ds)
    sh = prepare_svd_retrieval_setup(svd_cfg, collect(Float64.(λ)))
    return merge(ret, (
        PCs=sh.PCs,
        layout=sh.layout,
        log_trans=sh.log_trans,
        n_pc=sh.n_pc,
        n_leg=sh.n_leg,
        i678=argmin(abs.(λ .- 678.2)),
    ))
end

# ── main ──────────────────────────────────────────────────────────────────────

function main_ensemble()
    Random.seed!(ENSEMBLE_SEED)

    println("="^60)
    println("Surrogate ensemble: N_SAMPLES=$N_SAMPLES  N_ATM_POOL=$N_ATM_POOL  N_CONT_POOL=$N_CONT_POOL")
    println("SVD config: $SVD_CONFIG")
    println("REUSE_TRUTH=$REUSE_TRUTH  ZERO_SIF=$ZERO_SIF")
    println("SIF strength range: $SIF_STRENGTH_RANGE_ENS")
    println("="^60)

    isfile(SVD_CONFIG) || error("Missing SVD config: $SVD_CONFIG")
    svd_cfg = TOML.parsefile(SVD_CONFIG)
    n_pc_cfg = Int(get(get(get(svd_cfg, "fit", Dict()), "svd", Dict()), "n_pc", -1))
    n_leg_cfg = Int(get(get(get(svd_cfg, "fit", Dict()), "svd", Dict()), "n_legendre", -1))
    println("Config check: n_pc=$n_pc_cfg  n_legendre=$n_leg_cfg")
    n_pc_cfg > 0 || error("Invalid n_pc in $SVD_CONFIG")
    n_leg_cfg > 0 || error("Invalid n_legendre in $SVD_CONFIG")

    tag = _nleg_tag(n_leg_cfg)
    plot_dir = joinpath(OUT_DIR, "plots_$(tag)")
    ret_nc = joinpath(OUT_DIR, "retrieval_ensemble_$(tag).nc")
    truth_nc = joinpath(OUT_DIR, ZERO_SIF ? "truth_ensemble_zeroSIF.nc" : "truth_ensemble.nc")
    mkpath(OUT_DIR)
    mkpath(plot_dir)
    println("Output tag: $tag")
    println("  truth → $truth_nc")
    println("  plots → $plot_dir")
    println("  retrieval → $ret_nc")

    ENV["PACE_LUT_INTERPOLATION"] = LUT_INTERPOLATION
    plot_only = lowercase(get(ENV, "PLOT_ONLY", "0")) in ("1", "true", "yes")

    # Zero-SIF must not reuse the nominal (non-zero) truth ensemble.
    reuse = REUSE_TRUTH && isfile(truth_nc)
    if ZERO_SIF && REUSE_TRUTH && isfile(truth_nc)
        # sanity: refuse reuse if file still has non-zero SIF strengths
        ds_chk = Dataset(truth_nc)
        smax = maximum(Float64.(ds_chk["sif_strength"][:]))
        close(ds_chk)
        if smax > 0
            @warn "ZERO_SIF set but $truth_nc has max sif_strength=$smax; regenerating"
            reuse = false
        end
    end

    if reuse
        ens = load_truth_nc(truth_nc)
        if ZERO_SIF
            smax = maximum(Float64.(ens.sif_strength))
            smax == 0 || error("Loaded truth is not zero-SIF (max strength=$smax)")
        end
    else
        println("\n=== Preparing MWE context ===")
        ctx = MWEF.prepare_mwe_inputs(CONFIG_PATH)
        mwe_cfg = TOML.parsefile(CONFIG_PATH)
        ens = generate_ensemble(ctx, mwe_cfg; n_samples=N_SAMPLES)
        if ZERO_SIF
            maximum(ens.sif_strength) == 0 || error("ZERO_SIF generation failed: max strength=$(maximum(ens.sif_strength))")
            println("  verified sif_strength ≡ 0 for all $(size(ens.R_toa_noisy, 2)) samples")
        end
        write_truth_nc(ens, truth_nc)
        plot_ensemble_diagnostics(ens; out_dir=plot_dir)
    end

    if plot_only
        isfile(ret_nc) || error("PLOT_ONLY set but missing $ret_nc")
        ret = load_retrieval_for_plots(ret_nc, svd_cfg, ens.λ)
    else
        ret = run_ensemble_retrieval(ens, svd_cfg; out_nc=ret_nc)
    end
    plot_retrieval_comparison(ens, ret; out_dir=plot_dir)

    if ZERO_SIF
        ok = findall(ret.status .== 1)
        if !isempty(ok)
            floor_mean = mean(Float64.(ret.sif_678_ret[ok]))
            floor_med = median(Float64.(ret.sif_678_ret[ok]))
            floor_std = std(Float64.(ret.sif_678_ret[ok]))
            println("Zero-SIF floor @678 — mean=$(round(floor_mean, digits=4)), median=$(round(floor_med, digits=4)), std=$(round(floor_std, digits=4))")
            open(joinpath(plot_dir, "zeroSIF_floor_summary.txt"), "w") do io
                println(io, "n_converged\t$(length(ok))")
                println(io, "sif678_ret_mean\t$floor_mean")
                println(io, "sif678_ret_median\t$floor_med")
                println(io, "sif678_ret_std\t$floor_std")
                println(io, "sif678_true_max\t$(maximum(Float64.(ret.sif_678_true[ok])))")
            end
        end
    end

    println("\nDone. Results in $OUT_DIR ($tag)")
    return (ens=ens, ret=ret)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_ensemble()
end
