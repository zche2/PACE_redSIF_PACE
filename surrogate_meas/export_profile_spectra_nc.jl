#!/usr/bin/env julia
# Rebuild surrogate spectra for profile 5511 from profile_5511_summary.txt and write NetCDF.
#
# Usage:
#   julia --project=. surrogate_meas/export_profile_spectra_nc.jl
#   SUMMARY=surrogate_meas/output/profile_5511_summary.txt julia --project=. surrogate_meas/export_profile_spectra_nc.jl

using Random
using Statistics
using LinearAlgebra
using TOML
using NCDatasets
using Dates
using JLD2

const SCRIPT_DIR = @__DIR__
const REPO_ROOT = dirname(SCRIPT_DIR)
include(joinpath(SCRIPT_DIR, "build_single_meas.jl"))

function _parse_summary(path::AbstractString)
    d = Dict{String,String}()
    for line in eachline(path)
        isempty(strip(line)) && continue
        parts = split(line, '\t')
        length(parts) >= 2 || continue
        d[parts[1]] = parts[2]
    end
    return d
end

function load_l1b_pixel_at(
    pace_path::AbstractString,
    pixel::Int,
    scan::Int;
    fit_λ_min::Float64=REFLECTANCE_λ_MIN,
    fit_λ_max::Float64=REFLECTANCE_λ_MAX,
)
    ds = Dataset(pace_path)
    λ_all = Float64.(ds["red_wavelength"][:])
    E_all = Float64.(ds["red_solar_irradiance"][:])
    R = Float64.(coalesce.(ds["radiance_red"][pixel, scan, :], NaN))
    sza = Float64(coalesce(ds["solar_zenith"][pixel, scan], NaN))
    vza = haskey(ds, "sensor_zenith") ? Float64(coalesce(ds["sensor_zenith"][pixel, scan], NaN)) : 0.0
    close(ds)
    all(isfinite, R) || error("Non-finite radiance at pixel=$pixel scan=$scan")
    obs_full = (
        λ=λ_all, E=E_all, R=R, sza=sza, vza=vza,
        pixel=pixel, scan=scan, pace_path=pace_path,
    )
    return subset_obs_wavelength(obs_full, fit_λ_min, fit_λ_max)
end

function load_sif_shape_at(
    sif_path::AbstractString,
    λ_dst::AbstractVector{<:Real},
    library_index::Int,
    strength::Float64,
)
    sif = JLD2.load(MWEF.must_exist(sif_path))
    shapes = Matrix{Float64}(sif["SIF_shapes"])
    λ_ref = Float64.(sif["SIF_wavelen"])
    1 <= library_index <= size(shapes, 2) || error("sif_library_index out of range")
    shape_band = map_sif_shape_to_bands(λ_ref, shapes[:, library_index], λ_dst)
    peak = maximum(abs.(shape_band))
    peak > 0 || error("SIF shape is all zero")
    SIF = (shape_band ./ peak) .* strength
    return (
        λ=collect(Float64.(λ_dst)),
        SIF=SIF,
        strength=strength,
        library_index=library_index,
        peak_ref=peak,
    )
end

function main_export()
    summary_path = get(ENV, "SUMMARY", joinpath(SCRIPT_DIR, "output", "profile_5511_summary.txt"))
    isfile(summary_path) || error("Missing summary: $summary_path")
    meta = _parse_summary(summary_path)

    profile_index = parse(Int, meta["profile_index"])
    decay = haskey(meta, "reflectance_decay") ? parse(Float64, meta["reflectance_decay"]) : REFLECTANCE_DECAY
    shape = Symbol(get(meta, "reflectance_shape", String(REFLECTANCE_SHAPE)))
    sza_deg = parse(Float64, meta["sza_deg"])
    vza_deg = parse(Float64, meta["vza_deg"])
    amf_down = parse(Float64, meta["amf_down"])
    amf_up = parse(Float64, meta["amf_up"])
    sif_idx = parse(Int, meta["sif_library_index"])
    sif_strength = parse(Float64, meta["sif_strength"])
    l1b_pixel = parse(Int, meta["l1b_pixel"])
    l1b_scan = parse(Int, meta["l1b_scan"])

    println("Exporting spectra for profile $profile_index")
    println("  decay=$decay  shape=$shape  SIF idx=$sif_idx strength=$sif_strength")
    println("  L1B pixel=$l1b_pixel scan=$l1b_scan")

    ENV["PACE_LUT_INTERPOLATION"] = LUT_INTERPOLATION
    Random.seed!(RANDOM_SEED)
    ctx = MWEF.prepare_mwe_inputs(CONFIG_PATH)
    cfg = TOML.parsefile(CONFIG_PATH)
    data_cfg = get(cfg, "data", Dict{String, Any}())
    solar_file = get(data_cfg, "solar_file", "solar_merged_20240731_600_33300_100.out")
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)

    ds = Dataset(TRANS_NC)
    prof = load_profile(ds, profile_index)
    close(ds)

    weights = layer_reflectance_weights(length(prof.temp); shape=shape, decay=decay)
    τ_layer, τ_cum = layer_optical_depth(
        ctx.spectral_axis, prof.p_full, prof.temp, prof.vcd_dry, prof.vcd_h2o,
        ctx.o2_sitp, ctx.h2o_sitp,
    )
    τ_2way = weighted_two_way_optical_depth(τ_cum, weights, amf_down, amf_up)
    T2_atm = exp.(-τ_2way)
    T_solar, _ = MWEF.load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
    T_solar = Float64.(T_solar)
    T_solar ./= maximum(T_solar)
    T2_hres = T2_atm .* T_solar

    K = ctx.kernel_rsr_out
    λ_band = collect(Float64.(ctx.λ))
    T2_band = vec(K * T2_hres)
    T_solar_band = vec(K * T_solar)
    T2_atm_band = vec(K * T2_atm)

    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    pace_file = get(pace_cfg, "pace_file", "sample_granule_20240830T131442_new_chl.nc")
    pace_path = isabspath(pace_file) ? pace_file : joinpath(ctx.paths.base_dir, pace_file)
    λ_min = Float64(get(get(cfg, "spectral", Dict()), "lambda_min_nm", 640.0))
    λ_max = Float64(get(get(cfg, "spectral", Dict()), "lambda_max_nm", 756.0))

    obs = load_l1b_pixel_at(pace_path, l1b_pixel, l1b_scan)
    fit = fit_legendre_reflectance(obs.λ, obs.E, obs.R, obs.sza; order=LEGENDRE_ORDER)
    resc = rescale_mean_radiance(fit.ρ, fit.R_fit; λ=obs.λ, mean_λ_min=λ_min, mean_λ_max=λ_max)
    reflectance = (obs=obs, fit=fit, resc=resc)

    T1 = build_one_way_transmittance(τ_layer, amf_up, K)
    sif_scaled = load_sif_shape_at(ctx.paths.sif_path, λ_band, sif_idx, sif_strength)
    R_sif_toa = sif_scaled.SIF .* T1.T1_band
    sif_nt = (
        λ=sif_scaled.λ,
        SIF=sif_scaled.SIF,
        R_sif_toa=R_sif_toa,
        strength=sif_scaled.strength,
        library_index=sif_scaled.library_index,
        peak_ref=sif_scaled.peak_ref,
    )
    sif = (sif=sif_nt, T1=T1)

    # Deterministic noise for reproducibility (summary used seeded randn originally)
    Random.seed!(RANDOM_SEED + profile_index)
    toa = rebuild_toa_with_noise(
        λ_band, T2_band, reflectance, sif;
        band_snr_coeffs=ctx.band_snr_coeffs,
        out_path=joinpath(OUTPUT_DIR, "toa_profile$(profile_index)_export_tmp.png"),
    )

    svd_cfg = TOML.parsefile(SVD_RETRIEVAL_CONFIG)
    retrieval = pseudo_retrieval_from_surrogate(
        toa, sif, λ_band, ctx.paths.base_dir;
        retrieval_cfg=svd_cfg,
        out_path=joinpath(OUTPUT_DIR, "retrieval_profile$(profile_index)_export_tmp.png"),
    )

    out_nc = joinpath(OUTPUT_DIR, "profile_$(profile_index)_spectra.nc")
    write_surrogate_profile_nc(
        out_nc;
        profile_index=profile_index,
        prof=prof,
        weights=weights,
        reflectance_decay=decay,
        reflectance_shape=String(shape),
        amf_down=amf_down,
        amf_up=amf_up,
        λ_hres=collect(Float64.(ctx.λ_hres)),
        λ_band=λ_band,
        τ_2way=τ_2way,
        T_solar=T_solar,
        T2_atm=T2_atm,
        T2_hres=T2_hres,
        T_solar_band=T_solar_band,
        T2_atm_band=T2_atm_band,
        T2_band=T2_band,
        reflectance=reflectance,
        sif=sif,
        toa=toa,
        retrieval=retrieval,
    )
    println("Done → $out_nc")
    return out_nc
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_export()
end
