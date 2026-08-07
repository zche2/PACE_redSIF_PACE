#!/usr/bin/env julia
# Match-up SIF retrieval: load TROPOMI↔OCI matchup.nc and run the SVD transmittance LM
# from global_svd_fit_pipeline/svd_retrieval (same FM / LM as swath retrieval).
#
# Usage:
#   julia -t 8 TROPOMI_OCI_colocation/sif_retrieval/run_matchup_sif.jl path/to/config.sif_retrieval.toml
#
# Runs twice per match (adjustable via [matchup].spectra):
#   lt_oci_sim_noisy + sza_trop
#   lt_oci_meas      + sza_pace

using Base.Threads
using Dates
using ForwardDiff
using Interpolations
using JLD2
using LinearAlgebra
using NCDatasets
using SparseArrays
using Statistics
using TOML

const _THIS_DIR = @__DIR__
const _COLO_DIR = dirname(_THIS_DIR)
const _REPO_ROOT = dirname(_COLO_DIR)
const _SVD_DIR = joinpath(_REPO_ROOT, "global_svd_fit_pipeline", "svd_retrieval")

include(joinpath(_SVD_DIR, "svd_helpers.jl"))

"""Local copy of `load_sif_basis` (avoids pulling SimplePACEXSecFit / PACE_SIF)."""
function load_sif_basis_matchup(
    sif_path::AbstractString,
    λ::AbstractVector{<:Real};
    nEV::Int = 2,
    normalize::Bool = true,
)
    nEV < 1 && error("nEV must be >= 1")
    isfile(sif_path) || error("SIF file not found: $sif_path")
    sif = JLD2.load(sif_path)
    sif_u = convert.(Float64, sif["SIF_U"])
    λ_ref = collect(Float64.(sif["SIF_wavelen"]))
    size(sif_u, 1) == length(λ_ref) || error("SIF_U first dim must match SIF_wavelen")
    n_use = min(nEV, size(sif_u, 2))
    dλ = diff(λ_ref)
    step = dλ[1]
    tol = max(1e-10, abs(step) * 1e-8)
    use_cubic = all(abs.(dλ .- step) .<= tol)
    λ_knots = use_cubic ? range(λ_ref[1], step = step, length = length(λ_ref)) : nothing
    basis = zeros(Float64, length(λ), n_use)
    for iev in 1:n_use
        if use_cubic
            itp_ev = CubicSplineInterpolation(λ_knots, sif_u[:, iev]; extrapolation_bc = Line())
            vals = itp_ev.(λ)
            vals[(λ .< λ_ref[1]) .| (λ .> λ_ref[end])] .= 0.0
            basis[:, iev] .= vals
        else
            itp_ev = LinearInterpolation(λ_ref, sif_u[:, iev]; extrapolation_bc = 0.0)
            basis[:, iev] .= itp_ev.(λ)
        end
    end
    if normalize
        s = maximum(abs.(basis[:, 1]))
        s > 0 && (basis ./= s)
    end
    return basis
end

function _abspath(p::AbstractString)
    return isabspath(p) ? String(p) : normpath(joinpath(_REPO_ROOT, p))
end

function _resolve_data(cfg::Dict, key::String)
    data = get(cfg, "data", Dict{String, Any}())
    raw = String(data[key])
    p = _abspath(raw)
    isfile(p) && return p
    base = String(get(data, "base_dir", ""))
    isempty(base) && error("Missing file for [data].$key=$raw")
    cand = joinpath(base, raw)
    isfile(cand) || error("Missing file for [data].$key: tried $p and $cand")
    return cand
end

function _l1b_var(ds, group_name::String, var_name::String)
    if haskey(ds.group, group_name)
        g = ds.group[group_name]
        haskey(g, var_name) && return g[var_name]
    end
    haskey(ds, var_name) && return ds[var_name]
    error("L1B variable $group_name/$var_name not found")
end

"""Read OCI red F0 from a reference L1B; interpolate onto matchup band λ."""
function _load_f0_on_bands(l1b_path::AbstractString, λ_bands::Vector{Float64})
    ds = Dataset(l1b_path)
    try
        E0 = vec(Float64.(_l1b_var(ds, "sensor_band_parameters", "red_solar_irradiance")[:]))
        wl = vec(Float64.(_l1b_var(ds, "sensor_band_parameters", "red_wavelength")[:]))
        es = Float64(ds.attrib["earth_sun_distance_correction"])
        order = sortperm(wl)
        itp = LinearInterpolation(wl[order], E0[order]; extrapolation_bc = Flat())
        E0_b = Float64[itp(λ) for λ in λ_bands]
        return E0_b, es
    finally
        close(ds)
    end
end

function _band_mask(λ_all::Vector{Float64}, λ_min::Float64, λ_max::Float64)
    return findall(x -> λ_min < x < λ_max, λ_all)
end

function _setup_bases(cfg::Dict, λ_ctx::Vector{Float64})
    spectral = get(cfg, "spectral", Dict{String, Any}())
    fit = get(cfg, "fit", Dict{String, Any}())
    svd = get(fit, "svd", Dict{String, Any}())
    λ_min = Float64(get(spectral, "lambda_min_nm", 660.0))
    λ_max = Float64(get(spectral, "lambda_max_nm", 720.0))
    n_pc = Int(get(svd, "n_pc", 10))
    n_leg = Int(get(svd, "n_legendre", 5))
    log_trans = Bool(get(svd, "svd_log_transform", true))
    sif_nev = Int(get(spectral, "sif_nev", 1))
    normalize_sif = Bool(get(spectral, "normalize_sif_first_ev", true))

    summer = _resolve_data(cfg, "summer_nc")
    winter = _resolve_data(cfg, "winter_nc")
    sif_path = _resolve_data(cfg, "sif_file")

    svd_basis = load_svd_basis(summer, winter, λ_ctx; λ_min = λ_min, λ_max = λ_max, n_pc = n_pc, log_transform = log_trans)
    PCs = Float64.(svd_basis.PCs[:, 1:n_pc])
    sif_basis = load_sif_basis_matchup(sif_path, λ_ctx; nEV = sif_nev, normalize = normalize_sif)
    layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_leg, n_ev = size(sif_basis, 2))
    i678 = argmin(abs.(λ_ctx .- 678.2))
    sif_678 = vec(sif_basis[i678, :])

    prior_sigma_default = Float64(get(fit, "prior_sigma_default", 1e12))
    prior_min_sigma = Float64(get(fit, "prior_min_sigma", 1e-3))
    sif_sigma = Float64(get(fit, "sif_sigma", 1e12))
    alpha_mean = Float64(get(svd, "alpha_prior_mean", 1.5))
    alpha_sigma = Float64(get(svd, "alpha_prior_sigma", 0.5))
    use_leg01 = Bool(get(svd, "use_legendre01_prior", true))
    leg01_frac = Float64(get(svd, "legendre01_prior_sigma_fraction", 1.0))
    use_leghig = Bool(get(svd, "use_legendre_higher_prior", true))
    leg_higher_sigma = Float64(get(svd, "legendre_higher_sigma", 1.0))
    pc_prior_mode = String(get(svd, "pc_prior_mode", "loading_variance"))
    pc_sigma_scale = Float64(get(svd, "pc_prior_sigma_scale", 1.0))

    x0 = zeros(Float64, layout.n_state)
    x0[first(layout.idx_alpha)] = alpha_mean
    x0[first(layout.idx_legendre)] = 1.0
    prior_sigma = fill(prior_sigma_default, layout.n_state)
    prior_sigma[first(layout.idx_alpha)] = max(alpha_sigma, prior_min_sigma)
    if pc_prior_mode == "loading_variance"
        n_prof = svd_basis.n_profiles
        for k in 1:n_pc
            prior_sigma[k] = max(svd_basis.S[k] / sqrt(Float64(n_prof)) * pc_sigma_scale, prior_min_sigma)
        end
    end
    if use_leghig && length(layout.idx_legendre) >= 3
        for j in 3:length(layout.idx_legendre)
            prior_sigma[layout.idx_legendre[j]] = max(leg_higher_sigma, prior_min_sigma)
        end
    end
    prior_sigma[layout.idx_sif] .= max(sif_sigma, prior_min_sigma)

    z = _normalized_grid(λ_ctx)
    leg_basis = Float64.(_legendre_design_matrix(z, n_leg))
    A01 = hcat(ones(length(z)), z)
    lower = fill(-Inf, layout.n_state)
    upper = fill(Inf, layout.n_state)
    x_scale = ones(Float64, layout.n_state)
    for k in 1:n_pc
        x_scale[k] = prior_sigma[k]
    end
    x_scale[first(layout.idx_alpha)] = prior_sigma[first(layout.idx_alpha)]

    use_band_snr = Bool(get(fit, "use_band_snr", true))
    meas_sigma = Float64(get(fit, "meas_sigma", 0.01))
    snr_path = _resolve_data(cfg, "pace_snr_file")
    band_snr = if use_band_snr
        load_pace_band_snr_coeffs(snr_path, λ_ctx; λ_min = λ_min, λ_max = λ_max)
    else
        nothing
    end

    lm = (
        lambda0 = Float64(get(fit, "lm_lambda0", 1.0)),
        lambda_up = Float64(get(fit, "lm_lambda_up", 5.0)),
        lambda_down = Float64(get(fit, "lm_lambda_down", 0.7)),
        lambda_min = Float64(get(fit, "lm_lambda_min", 1e-8)),
        lambda_max = Float64(get(fit, "lm_lambda_max", 1e8)),
        max_inner = Int(get(fit, "lm_max_inner", 40)),
    )
    conv = (
        dx_rel_tol = Float64(get(fit, "conv_dx_rel_tol", 1e-6)),
        rmse_rel_tol = Float64(get(fit, "conv_rmse_rel_tol", 1e-6)),
        rmse_abs_tol = Float64(get(fit, "conv_rmse_abs_tol", 1e-6)),
        enabled = Bool(get(fit, "conv_stall_enable", true)),
        window = Int(get(fit, "conv_stall_window", 3)),
        redchi2_target = Float64(get(fit, "conv_stall_redchi2_target", 5.0)),
        redchi2_abs_tol = Float64(get(fit, "conv_stall_redchi2_abs_tol", 0.1)),
        redchi2_rel_tol = Float64(get(fit, "conv_stall_redchi2_rel_tol", 0.03)),
        dx_rel_tol_stall = Float64(get(fit, "conv_stall_dx_rel_tol", 5e-3)),
    )
    max_outer = Int(get(get(cfg, "batch_fit", Dict{String, Any}()), "max_outer_steps", 10))
    stall_dx = Float64(get(fit, "conv_stall_dx_rel_tol", 5e-3))

    return (;
        PCs,
        sif_basis,
        sif_678,
        layout,
        x0,
        prior_sigma,
        prior_min_sigma,
        use_leg01,
        leg01_frac,
        A01,
        lower,
        upper,
        x_scale,
        use_band_snr,
        band_snr,
        meas_sigma,
        lm,
        conv,
        max_outer,
        stall_dx,
        log_trans,
        n_pc,
        n_leg,
        λ_min,
        λ_max,
        state_names = svd_state_names(layout),
    )
end

function _retrieve_one_λ!(
    x_out::Vector{Float64},
    λ_ctx::Vector{Float64},
    y_obs::Vector{Float64},
    solar_eff::Vector{Float64},
    sh,
)
    fm, layout = make_svd_forward_model_λ(
        λ_ctx,
        solar_eff,
        sh.PCs,
        sh.sif_basis;
        n_pc = sh.n_pc,
        n_legendre = sh.n_leg,
        log_transform = sh.log_trans,
    )
    jac = (x -> ForwardDiff.jacobian(fm, x))
    return _run_one_svd_retrieval!(
        x_out,
        fm,
        jac,
        layout,
        y_obs,
        sh.x0,
        sh.x0,
        sh.prior_sigma,
        sh.lower,
        sh.upper,
        sh.x_scale,
        sh.use_leg01,
        sh.A01,
        sh.prior_min_sigma,
        sh.leg01_frac,
        sh.use_band_snr,
        sh.band_snr,
        sh.meas_sigma,
        sh.lm,
        sh.conv,
        sh.stall_dx,
        sh.max_outer,
    )
end

function _spectra_specs(cfg::Dict)
    m = get(cfg, "matchup", Dict{String, Any}())
    specs = get(m, "spectra", ["sim_noisy", "meas"])
    out = NamedTuple[]
    for s in specs
        s = String(s)
        if s == "sim_noisy"
            push!(out, (name = "sim_noisy", lt_var = "lt_oci_sim_noisy", sza_var = "sza_trop"))
        elseif s == "meas"
            push!(out, (name = "meas", lt_var = "lt_oci_meas", sza_var = "sza_pace"))
        elseif s == "sim"
            push!(out, (name = "sim", lt_var = "lt_oci_sim", sza_var = "sza_trop"))
        else
            error("Unknown [matchup].spectra entry: $s (use sim_noisy|meas|sim)")
        end
    end
    return out
end

function run_matchup_sif(config_path::AbstractString)
    cfg = TOML.parsefile(config_path)
    mcfg = get(cfg, "matchup", Dict{String, Any}())
    matchup_nc = _abspath(String(mcfg["matchup_nc"]))
    isfile(matchup_nc) || error("matchup_nc not found: $matchup_nc")

    spectral = get(cfg, "spectral", Dict{String, Any}())
    fit = get(cfg, "fit", Dict{String, Any}())
    svd = get(fit, "svd", Dict{String, Any}())
    n_pc = Int(get(svd, "n_pc", 10))
    n_leg = Int(get(svd, "n_legendre", 5))
    # Always name: {output_root}/sif_retrieval_nPC_{n}_npoly_{m}
    out_root = if haskey(mcfg, "output_root")
        _abspath(String(mcfg["output_root"]))
    elseif haskey(mcfg, "output_dir")
        dirname(_abspath(String(mcfg["output_dir"])))
    else
        error("Set [matchup].output_root or output_dir")
    end
    out_dir = joinpath(out_root, "sif_retrieval_nPC_$(n_pc)_npoly_$(n_leg)")
    mkpath(out_dir)
    l1b_ref = _abspath(String(get(mcfg, "l1b_solar_reference", "")))
    isempty(l1b_ref) && error("Set [matchup].l1b_solar_reference to an OCI L1B for F0/esd")
    isfile(l1b_ref) || error("l1b_solar_reference not found: $l1b_ref")
    max_n = Int(get(mcfg, "max_matches", 0))  # 0 = all
    esd_override = get(mcfg, "esd", nothing)

    λ_min = Float64(get(spectral, "lambda_min_nm", 660.0))
    λ_max = Float64(get(spectral, "lambda_max_nm", 720.0))

    println("Threads: ", nthreads())
    println("matchup: ", matchup_nc)
    println("output:  ", out_dir)
    println("λ window: ($λ_min, $λ_max) nm")

    ds = Dataset(matchup_nc)
    λ_all = vec(Float64.(ds["band"][:]))
    ib = _band_mask(λ_all, λ_min, λ_max)
    isempty(ib) && error("No matchup bands in ($λ_min, $λ_max) nm")
    λ_ctx = λ_all[ib]
    n_match_all = Int(ds.dim["match"])
    n_match = max_n > 0 ? min(max_n, n_match_all) : n_match_all
    println("bands kept: ", length(λ_ctx), " / ", length(λ_all), "   matches: ", n_match, " / ", n_match_all)

    # load geometries used for both runs
    lat_pace = Float32.(ds["lat_pace"][1:n_match])
    lon_pace = Float32.(ds["lon_pace"][1:n_match])
    lat_trop = Float32.(ds["lat_trop"][1:n_match])
    lon_trop = Float32.(ds["lon_trop"][1:n_match])
    sza_pace = Float64.(ds["sza_pace"][1:n_match])
    sza_trop = Float64.(ds["sza_trop"][1:n_match])
    vza_pace = Float32.(ds["vza_pace"][1:n_match])
    vza_trop = Float32.(ds["vza_trop"][1:n_match])

    E0_b, esd_l1b = _load_f0_on_bands(l1b_ref, λ_ctx)
    esd = esd_override === nothing ? esd_l1b : Float64(esd_override)
    println("F0 from L1B ref; esd=", esd)

    sh = _setup_bases(cfg, λ_ctx)
    println("state: ", join(sh.state_names, ", "))
    specs = _spectra_specs(cfg)

    for spec in specs
        println("="^72)
        println("Retrieving spectrum=$(spec.name)  lt=$(spec.lt_var)  sza=$(spec.sza_var)")
        # NetCDF dims for lt_* are ("band","match")
        Lt_bm = Float64.(ds[spec.lt_var][ib, 1:n_match])  # (n_band, n_match)
        SZA = spec.sza_var == "sza_trop" ? sza_trop : sza_pace

        n_state = sh.layout.n_state
        X = fill(NaN32, n_match, n_state)
        converged = fill(UInt8(0), n_match)
        status = fill(Int16(-1), n_match)
        n_steps = fill(Int16(0), n_match)
        rmse = fill(NaN32, n_match)
        redchi2 = fill(NaN32, n_match)
        sif_ev1 = fill(NaN32, n_match)
        sif_678 = fill(NaN32, n_match)

        @threads for i in 1:n_match
            y = vec(Lt_bm[:, i])
            if !all(isfinite, y) || !isfinite(SZA[i])
                status[i] = Int16(3)
                continue
            end
            solar_eff = @. E0_b * cosd(SZA[i]) / π / esd
            x_out = zeros(Float64, n_state)
            try
                res = _retrieve_one_λ!(x_out, λ_ctx, y, solar_eff, sh)
                X[i, :] .= Float32.(x_out)
                converged[i] = res.converged ? UInt8(1) : UInt8(0)
                status[i] = res.status
                n_steps[i] = Int16(res.n_steps)
                rmse[i] = Float32(res.rmse)
                redchi2[i] = Float32(res.reduced_chi2)
                sif_c = x_out[sh.layout.idx_sif]
                sif_ev1[i] = Float32(sif_c[1])
                sif_678[i] = Float32(dot(sh.sif_678, sif_c))
            catch
                status[i] = Int16(4)
            end
        end

        out_nc = joinpath(out_dir, "matchup_sif_$(spec.name).nc")
        isfile(out_nc) && rm(out_nc)
        ds_o = Dataset(out_nc, "c")
        defDim(ds_o, "match", n_match)
        defDim(ds_o, "band", length(λ_ctx))
        defDim(ds_o, "state", n_state)
        ds_o.attrib["title"] = "Match-up SVD SIF retrieval ($(spec.name))"
        ds_o.attrib["source_matchup"] = matchup_nc
        ds_o.attrib["config_file"] = abspath(config_path)
        ds_o.attrib["spectrum"] = spec.name
        ds_o.attrib["lt_var"] = spec.lt_var
        ds_o.attrib["sza_var"] = spec.sza_var
        ds_o.attrib["state_names_csv"] = join(sh.state_names, ",")
        ds_o.attrib["lambda_min_nm"] = λ_min
        ds_o.attrib["lambda_max_nm"] = λ_max
        ds_o.attrib["n_pc"] = sh.n_pc
        ds_o.attrib["n_legendre"] = sh.n_leg
        ds_o.attrib["esd"] = esd
        ds_o.attrib["created"] = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS")

        defVar(ds_o, "band", Float32.(λ_ctx), ("band",); attrib = Dict("units" => "nm"))
        defVar(ds_o, "lat_pace", lat_pace, ("match",))
        defVar(ds_o, "lon_pace", lon_pace, ("match",))
        defVar(ds_o, "lat_trop", lat_trop, ("match",))
        defVar(ds_o, "lon_trop", lon_trop, ("match",))
        defVar(ds_o, "sza_pace", Float32.(sza_pace), ("match",); attrib = Dict("units" => "degree"))
        defVar(ds_o, "sza_trop", Float32.(sza_trop), ("match",); attrib = Dict("units" => "degree"))
        defVar(ds_o, "vza_pace", vza_pace, ("match",))
        defVar(ds_o, "vza_trop", vza_trop, ("match",))
        defVar(ds_o, "sza_used", Float32.(SZA), ("match",); attrib = Dict("long_name" => "SZA used in solar_eff for this retrieval"))
        defVar(ds_o, "x_hat", X, ("match", "state"))
        defVar(ds_o, "converged", converged, ("match",))
        defVar(ds_o, "status_code", status, ("match",))
        defVar(ds_o, "n_steps", n_steps, ("match",))
        defVar(ds_o, "rmse", rmse, ("match",))
        defVar(ds_o, "reduced_chi2", redchi2, ("match",))
        defVar(ds_o, "sif_ev1", sif_ev1, ("match",))
        defVar(ds_o, "sif_radiance_678nm", sif_678, ("match",); attrib = Dict("units" => "W m-2 sr-1 um-1"))
        close(ds_o)
        n_ok = count(==(Int16(1)), status)
        println("wrote $out_nc  converged=$n_ok / $n_match")
    end
    close(ds)
    println("Done.")
end

function main(args = ARGS)
    isempty(args) && error("Usage: julia -t N run_matchup_sif.jl path/to/config.sif_retrieval.toml")
    run_matchup_sif(args[1])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
