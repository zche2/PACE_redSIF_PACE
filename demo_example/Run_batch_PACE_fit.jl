#!/usr/bin/env julia

using TOML
using NCDatasets
using LinearAlgebra
using SparseArrays
using Statistics
using Dates

# Reuse forward model + LM/Jacobian utilities without executing the single-pixel main().
include(joinpath(@__DIR__, "Fit_toy_forward_model.jl"))

const MWEF = SimplePACEXSecFitMWEFunctions

@inline function cfg_get(cfg::Dict, section::String, key::String, default)
    haskey(cfg, section) || return default
    sec = cfg[section]
    sec isa Dict || return default
    return get(sec, key, default)
end

function resolve_orbit_files(cfg::Dict)
    base_dir = String(cfg_get(cfg, "data", "base_dir", ""))
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())

    files = String[]
    cfg_files = get(batch_cfg, "pace_files", String[])
    if cfg_files isa AbstractVector && !isempty(cfg_files)
        for f in cfg_files
            fs = String(f)
            push!(files, isabspath(fs) ? fs : joinpath(base_dir, fs))
        end
    else
        one_file = String(get(pace_cfg, "pace_file", ""))
        isempty(one_file) && error("No pace file configured (pace_observation.pace_file)")
        push!(files, isabspath(one_file) ? one_file : joinpath(base_dir, one_file))
    end
    return files
end

function make_output_path(pace_path::AbstractString, cfg::Dict)
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    out_dir_default = joinpath(@__DIR__, "batch_output")
    out_dir_cfg = String(get(batch_cfg, "output_dir", out_dir_default))
    out_dir = isabspath(out_dir_cfg) ? out_dir_cfg : joinpath(@__DIR__, out_dir_cfg)
    mkpath(out_dir)

    suffix = String(get(batch_cfg, "output_suffix", "_retrieval.nc"))
    stem = splitext(basename(pace_path))[1]
    return joinpath(out_dir, stem * suffix)
end

function find_axis_indices(v, wavelength_var::AbstractString)
    dnames = collect(String.(dimnames(v)))
    i_pix = findfirst(==("pixels"), dnames)
    i_scan = findfirst(==("scans"), dnames)
    i_band = findfirst(==(String(wavelength_var)), dnames)
    if isnothing(i_band)
        i_band = findfirst(x -> occursin("band", lowercase(x)), dnames)
    end
    if isnothing(i_pix) || isnothing(i_scan) || isnothing(i_band)
        error("Could not infer (pixels, scans, bands) dims from $(dnames)")
    end
    return (
        dnames = dnames,
        i_pix = i_pix,
        i_scan = i_scan,
        i_band = i_band,
        n_pix = size(v, i_pix),
        n_scan = size(v, i_scan),
        n_band = size(v, i_band),
    )
end

function read_geo_2d(ds::NCDataset, varname::AbstractString, n_pix::Int, n_scan::Int)
    haskey(ds, varname) || return fill(Float32(NaN), n_pix, n_scan)
    v = ds[varname]
    ndims(v) == 2 || error("Expected 2D geolocation variable '$varname', got ndims=$(ndims(v))")
    d = collect(String.(dimnames(v)))
    raw = v[:, :]

    arr = if d == ["pixels", "scans"]
        raw
    elseif d == ["scans", "pixels"]
        permutedims(raw, (2, 1))
    else
        error("Unsupported geolocation dims for '$varname': $d")
    end

    out = Array{Float32}(undef, size(arr)...)
    @inbounds for j in axes(arr, 2), i in axes(arr, 1)
        a = arr[i, j]
        if ismissing(a)
            out[i, j] = Float32(NaN)
        else
            af = Float64(a)
            out[i, j] = isfinite(af) ? Float32(af) : Float32(NaN)
        end
    end
    return out
end

function make_linear_resampler(λ_src_in::AbstractVector{<:Real}, λ_dst::AbstractVector{<:Real})
    λ_src = collect(Float64.(λ_src_in))
    λ_dst_f = collect(Float64.(λ_dst))
    perm = sortperm(λ_src)
    λs = λ_src[perm]

    lo, hi = extrema(λs)
    (minimum(λ_dst_f) >= lo && maximum(λ_dst_f) <= hi) ||
        error("Target λ range [$(minimum(λ_dst_f)), $(maximum(λ_dst_f))] outside source λ range [$lo, $hi]")

    I = Int[]
    J = Int[]
    V = Float64[]
    sizehint!(I, 2 * length(λ_dst_f))
    sizehint!(J, 2 * length(λ_dst_f))
    sizehint!(V, 2 * length(λ_dst_f))

    for (i, λ) in enumerate(λ_dst_f)
        j_hi = searchsortedfirst(λs, λ)
        if j_hi <= 1
            push!(I, i); push!(J, 1); push!(V, 1.0)
        elseif j_hi > length(λs)
            push!(I, i); push!(J, length(λs)); push!(V, 1.0)
        else
            j_lo = j_hi - 1
            λ_lo = λs[j_lo]
            λ_hi = λs[j_hi]
            if λ_hi == λ_lo
                push!(I, i); push!(J, j_lo); push!(V, 1.0)
            else
                w_hi = (λ - λ_lo) / (λ_hi - λ_lo)
                w_lo = 1.0 - w_hi
                push!(I, i); push!(J, j_lo); push!(V, w_lo)
                push!(I, i); push!(J, j_hi); push!(V, w_hi)
            end
        end
    end
    W = sparse(I, J, V, length(λ_dst_f), length(λs))
    return W, perm
end

@inline function copy_sorted_spectrum!(
    y_sorted::AbstractVector{Float64},
    spec_raw::AbstractVector,
    perm::AbstractVector{Int},
)
    @inbounds for i in eachindex(perm)
        v = spec_raw[perm[i]]
        if ismissing(v)
            return false
        end
        vf = Float64(v)
        if !isfinite(vf)
            return false
        end
        y_sorted[i] = vf
    end
    return true
end

function create_output_dataset(
    output_path::AbstractString,
    n_pix::Int,
    n_scan::Int,
    state_names::Vector{String},
    pace_path::AbstractString,
    config_path::AbstractString,
    pixel_range::UnitRange{Int},
    scan_range::UnitRange{Int},
)
    ds = Dataset(output_path, "c")
    defDim(ds, "pixels", n_pix)
    defDim(ds, "scans", n_scan)
    defDim(ds, "state", length(state_names))

    ds.attrib["title"] = "PACE toy retrieval swath output"
    ds.attrib["history"] = "Created " * Dates.format(now(), Dates.DateFormat("yyyy-mm-ddTHH:MM:SS"))
    ds.attrib["input_pace_file"] = String(pace_path)
    ds.attrib["config_file"] = String(config_path)
    ds.attrib["state_names_csv"] = join(state_names, ",")
    ds.attrib["pixel_start"] = first(pixel_range)
    ds.attrib["pixel_end"] = last(pixel_range)
    ds.attrib["scan_start"] = first(scan_range)
    ds.attrib["scan_end"] = last(scan_range)

    v_lat = defVar(ds, "latitude", Float32, ("pixels", "scans"))
    v_lon = defVar(ds, "longitude", Float32, ("pixels", "scans"))
    v_state = defVar(ds, "x_hat", Float32, ("pixels", "scans", "state"))
    v_conv = defVar(ds, "converged", UInt8, ("pixels", "scans"))
    v_status = defVar(ds, "status_code", Int16, ("pixels", "scans"))
    v_steps = defVar(ds, "n_steps", Int16, ("pixels", "scans"))
    v_rmse = defVar(ds, "rmse", Float32, ("pixels", "scans"))
    v_rchi2 = defVar(ds, "reduced_chi2", Float32, ("pixels", "scans"))
    v_obj = defVar(ds, "objective", Float32, ("pixels", "scans"))
    v_sif1 = defVar(ds, "sif_ev1", Float32, ("pixels", "scans"))
    v_dark = defVar(ds, "is_dark", UInt8, ("pixels", "scans"))
    v_ocean = defVar(ds, "is_ocean", UInt8, ("pixels", "scans"))
    v_pixsrc = defVar(ds, "source_pixel_index", Int32, ("pixels",))
    v_scansrc = defVar(ds, "source_scan_index", Int32, ("scans",))

    v_state.attrib["long_name"] = "Retrieved state vector"
    v_conv.attrib["long_name"] = "1 if convergence criterion reached, 0 otherwise"
    v_status.attrib["long_name"] = "0=ok_not_converged, 1=converged, 2=no_accepted_step, 3=invalid_input, 4=model_failure, 5=non_dark_spectrum_skipped, 6=non_ocean_spectrum_skipped"
    v_steps.attrib["long_name"] = "Number of accepted LM outer steps"
    v_rmse.attrib["long_name"] = "Spectral RMSE between observation and model"
    v_rchi2.attrib["long_name"] = "Reduced chi-square: sum((r/sigma)^2)/(n_meas-n_state)"
    v_obj.attrib["long_name"] = "Final MAP objective value"
    v_sif1.attrib["long_name"] = "Retrieved first SIF eigenvector coefficient (sif_ev1)"
    v_dark.attrib["long_name"] = "1 if spectrum passed dark-scene threshold, 0 otherwise"
    v_ocean.attrib["long_name"] = "1 if spectrum passed ocean-mask filter, 0 otherwise"

    v_pixsrc[:] = collect(Int32.(pixel_range))
    v_scansrc[:] = collect(Int32.(scan_range))

    return ds
end

function build_retrieval_core(config_path::AbstractString)
    cfg = TOML.parsefile(config_path)
    fit_cfg = get(cfg, "fit", Dict{String, Any}())
    data_cfg = get(cfg, "data", Dict{String, Any}())

    state_float_type = MWEF.parse_float_type(cfg)
    ctx = MWEF.prepare_mwe_inputs(config_path)

    n_legendre = Int(get(fit_cfg, "n_legendre", 2))
    preallocate_forward = Bool(get(fit_cfg, "preallocate_forward", true))
    preallocate_ad_forward = Bool(get(fit_cfg, "preallocate_ad_forward", false))
    preallocate_jacobian = Bool(get(fit_cfg, "preallocate_jacobian", false))
    use_hybrid_jacobian = Bool(get(fit_cfg, "use_hybrid_jacobian", false))
    conv_dx_rel_tol = Float64(get(fit_cfg, "conv_dx_rel_tol", 1e-6))
    conv_rmse_rel_tol = Float64(get(fit_cfg, "conv_rmse_rel_tol", 1e-6))
    conv_rmse_abs_tol = Float64(get(fit_cfg, "conv_rmse_abs_tol", 1e-6))
    conv_stall_enable = Bool(get(fit_cfg, "conv_stall_enable", true))
    conv_stall_window = Int(get(fit_cfg, "conv_stall_window", 3))
    conv_stall_redchi2_target = Float64(get(fit_cfg, "conv_stall_redchi2_target", 5.0))
    conv_stall_redchi2_abs_tol = Float64(get(fit_cfg, "conv_stall_redchi2_abs_tol", 0.1))
    conv_stall_redchi2_rel_tol = Float64(get(fit_cfg, "conv_stall_redchi2_rel_tol", 0.03))
    conv_stall_dx_rel_tol = Float64(get(fit_cfg, "conv_stall_dx_rel_tol", 5e-3))
    lm_lambda0 = Float64(get(fit_cfg, "lm_lambda0", 1.0))
    lm_lambda_up = Float64(get(fit_cfg, "lm_lambda_up", 5.0))
    lm_lambda_down = Float64(get(fit_cfg, "lm_lambda_down", 0.7))
    lm_lambda_min = Float64(get(fit_cfg, "lm_lambda_min", 1e-8))
    lm_lambda_max = Float64(get(fit_cfg, "lm_lambda_max", 1e8))
    lm_max_inner = Int(get(fit_cfg, "lm_max_inner", 24))
    meas_sigma = Float64(get(fit_cfg, "meas_sigma", 0.01))
    max_outer_steps_default = Int(get(fit_cfg, "n_plot_steps", 12))

    p_prior_hpa = Float64(get(fit_cfg, "p_prior_hpa", 700.0))
    p_sigma_hpa = Float64(get(fit_cfg, "p_sigma_hpa", 200.0))
    t_prior_k = Float64(get(fit_cfg, "t_prior_k", 280.0))
    t_sigma_k = Float64(get(fit_cfg, "t_sigma_k", 20.0))
    vcd_o2_sigma = Float64(get(fit_cfg, "vcd_o2_sigma", 1e23))
    vcd_h2o_sigma = Float64(get(fit_cfg, "vcd_h2o_sigma", 3e22))
    vcd_slope_prior_sigma_factor = Float64(get(fit_cfg, "vcd_slope_prior_sigma_factor", 1.0))
    use_vcd_slope_prior = Bool(get(fit_cfg, "use_vcd_slope_prior", true))
    use_legendre01_prior = Bool(get(fit_cfg, "use_legendre01_prior", true))
    legendre01_prior_sigma_fraction = Float64(get(fit_cfg, "legendre01_prior_sigma_fraction", 1.0))
    sif_sigma = Float64(get(fit_cfg, "sif_sigma", 1e12))
    prior_min_sigma = Float64(get(fit_cfg, "prior_min_sigma", 1e-3))
    prior_sigma_default = Float64(get(fit_cfg, "prior_sigma_default", 1e12))
    legendre_higher_sigma = Float64(get(fit_cfg, "legendre_higher_sigma", 1.0))
    use_legendre_higher_prior = Bool(get(fit_cfg, "use_legendre_higher_prior", true))
    use_pt_constraints = Bool(get(fit_cfg, "use_pt_constraints", true))
    pt_constraint_sigma_mult = Float64(get(fit_cfg, "pt_constraint_sigma_mult", 3.0))

    solar_file = String(get(data_cfg, "solar_file", "solar_merged_20200720_600_33300_100.out"))
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)
    solar_hres, _ = MWEF.load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
    solar_hres = state_float_type.(solar_hres)

    fm = make_forward_model_simple(
        ctx,
        solar_hres;
        n_legendre=n_legendre,
        preallocate_float64=preallocate_forward && state_float_type == Float64,
        preallocate_float32=preallocate_forward && state_float_type == Float32,
        preallocate_other_types=preallocate_ad_forward,
    )

    layout = state_layout_simple(ctx; n_legendre=n_legendre)
    x0 = initial_state_simple(ctx; n_legendre=n_legendre, T=state_float_type)
    jacobian_eval = if use_hybrid_jacobian
        make_hybrid_jacobian_evaluator(
            fm,
            ctx,
            solar_hres,
            layout;
            n_legendre=n_legendre,
        )
    else
        make_jacobian_evaluator(fm, x0; use_preallocated=preallocate_jacobian)
    end

    x_a_base = copy(Float64.(x0))
    prior_sigma_base = fill(prior_sigma_default, length(x0))

    # Priors analogous to the single-spectrum setup, but static for swath runs.
    x_a_base[layout.idx_p_o2_hpa] = p_prior_hpa
    x_a_base[layout.idx_p_h2o_hpa] = p_prior_hpa
    x_a_base[layout.idx_t_o2_k] = t_prior_k
    x_a_base[layout.idx_t_h2o_k] = t_prior_k
    prior_sigma_base[layout.idx_p_o2_hpa] = p_sigma_hpa
    prior_sigma_base[layout.idx_p_h2o_hpa] = p_sigma_hpa
    prior_sigma_base[layout.idx_t_o2_k] = t_sigma_k
    prior_sigma_base[layout.idx_t_h2o_k] = t_sigma_k

    x_a_base[layout.idx_vcd_o2_intercept] = x0[layout.idx_vcd_o2_intercept]
    x_a_base[layout.idx_vcd_h2o_intercept] = x0[layout.idx_vcd_h2o_intercept]
    prior_sigma_base[layout.idx_vcd_o2_intercept] = vcd_o2_sigma
    prior_sigma_base[layout.idx_vcd_h2o_intercept] = vcd_h2o_sigma
    x_a_base[layout.idx_vcd_o2_sif] = x0[layout.idx_vcd_o2_sif]
    x_a_base[layout.idx_vcd_h2o_sif] = x0[layout.idx_vcd_h2o_sif]
    prior_sigma_base[layout.idx_vcd_o2_sif] = vcd_o2_sigma
    prior_sigma_base[layout.idx_vcd_h2o_sif] = vcd_h2o_sigma

    if use_vcd_slope_prior
        x_a_base[layout.idx_vcd_o2_slope] = 0.0
        x_a_base[layout.idx_vcd_h2o_slope] = 0.0
        prior_sigma_base[layout.idx_vcd_o2_slope] = max(vcd_o2_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
        prior_sigma_base[layout.idx_vcd_h2o_slope] = max(vcd_h2o_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
    end

    x_a_base[layout.idx_sif] .= 0.0
    prior_sigma_base[layout.idx_sif] .= max(sif_sigma, prior_min_sigma)
    if use_legendre_higher_prior && length(layout.idx_legendre) >= 3
        for j in 3:length(layout.idx_legendre)
            idx = layout.idx_legendre[j]
            x_a_base[idx] = 0.0
            prior_sigma_base[idx] = max(legendre_higher_sigma, prior_min_sigma)
        end
    end

    S_e_inv = spdiagm(0 => fill(1.0 / (meas_sigma^2), length(ctx.λ)))

    x_scale_base = ones(Float64, length(x0))
    x_scale_base[layout.idx_vcd_o2_intercept] = vcd_o2_sigma
    x_scale_base[layout.idx_vcd_o2_slope] = max(vcd_o2_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
    x_scale_base[layout.idx_vcd_h2o_intercept] = vcd_h2o_sigma
    x_scale_base[layout.idx_vcd_h2o_slope] = max(vcd_h2o_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
    x_scale_base[layout.idx_vcd_o2_sif] = vcd_o2_sigma
    x_scale_base[layout.idx_vcd_h2o_sif] = vcd_h2o_sigma
    x_scale_base[layout.idx_p_o2_hpa] = p_sigma_hpa
    x_scale_base[layout.idx_p_h2o_hpa] = p_sigma_hpa
    x_scale_base[layout.idx_t_o2_k] = t_sigma_k
    x_scale_base[layout.idx_t_h2o_k] = t_sigma_k
    x_scale_base[layout.idx_sif] .= 1.0
    x_scale_base[layout.idx_legendre] .= 1.0

    lower_bounds = fill(-Inf, length(x0))
    upper_bounds = fill(Inf, length(x0))
    if use_pt_constraints
        lower_bounds[layout.idx_p_o2_hpa] = p_prior_hpa - pt_constraint_sigma_mult * p_sigma_hpa
        upper_bounds[layout.idx_p_o2_hpa] = p_prior_hpa + pt_constraint_sigma_mult * p_sigma_hpa
        lower_bounds[layout.idx_p_h2o_hpa] = p_prior_hpa - pt_constraint_sigma_mult * p_sigma_hpa
        upper_bounds[layout.idx_p_h2o_hpa] = p_prior_hpa + pt_constraint_sigma_mult * p_sigma_hpa
        lower_bounds[layout.idx_t_o2_k] = t_prior_k - pt_constraint_sigma_mult * t_sigma_k
        upper_bounds[layout.idx_t_o2_k] = t_prior_k + pt_constraint_sigma_mult * t_sigma_k
        lower_bounds[layout.idx_t_h2o_k] = t_prior_k - pt_constraint_sigma_mult * t_sigma_k
        upper_bounds[layout.idx_t_h2o_k] = t_prior_k + pt_constraint_sigma_mult * t_sigma_k
    end

    z = _normalized_grid(ctx.λ)
    A01 = hcat(ones(length(z)), z)

    return (
        cfg = cfg,
        ctx = ctx,
        fm = fm,
        layout = layout,
        state_names = state_names_simple(ctx; n_legendre=n_legendre),
        x0_base = copy(Float64.(x0)),
        x_a_base = x_a_base,
        prior_sigma_base = prior_sigma_base,
        use_legendre01_prior = use_legendre01_prior,
        legendre01_prior_sigma_fraction = legendre01_prior_sigma_fraction,
        prior_min_sigma = prior_min_sigma,
        A01 = A01,
        jacobian_eval = jacobian_eval,
        S_e_inv = S_e_inv,
        x_scale_base = x_scale_base,
        lower_bounds = lower_bounds,
        upper_bounds = upper_bounds,
        meas_sigma = meas_sigma,
        conv_dx_rel_tol = conv_dx_rel_tol,
        conv_rmse_rel_tol = conv_rmse_rel_tol,
        conv_rmse_abs_tol = conv_rmse_abs_tol,
        conv_stall = (
            enabled = conv_stall_enable,
            window = conv_stall_window,
            redchi2_target = conv_stall_redchi2_target,
            redchi2_abs_tol = conv_stall_redchi2_abs_tol,
            redchi2_rel_tol = conv_stall_redchi2_rel_tol,
            dx_rel_tol = conv_stall_dx_rel_tol,
        ),
        lm = (
            lambda0 = lm_lambda0,
            lambda_up = lm_lambda_up,
            lambda_down = lm_lambda_down,
            lambda_min = lm_lambda_min,
            lambda_max = lm_lambda_max,
            max_inner = lm_max_inner,
            max_outer_default = max_outer_steps_default,
        ),
    )
end

function run_one_retrieval!(
    x_out::Vector{Float64},
    core,
    y_obs::Vector{Float64},
    max_outer_steps::Int,
)
    # Match single-spectrum bootstrap exactly:
    # 1) start from base x0 and static priors,
    # 2) derive Legendre P0/P1 priors from ratio against f(x0).
    x0 = copy(core.x0_base)

    x_a = copy(core.x_a_base)
    prior_sigma = copy(core.prior_sigma_base)
    layout = core.layout

    # Per-spectrum Legendre P0/P1 priors from continuum ratio fit.
    if core.use_legendre01_prior && length(layout.idx_legendre) >= 1
        y_base = core.fm(x0)
        ratio = y_obs ./ max.(abs.(y_base), eps(Float64))
        w = y_obs .- minimum(y_obs)
        w .+= max(maximum(w), 1.0) * 1e-6
        s = sqrt.(w ./ maximum(w))
        c01 = (core.A01 .* s) \ (ratio .* s)

        leg0_idx = first(layout.idx_legendre)
        x_a[leg0_idx] = c01[1]
        prior_sigma[leg0_idx] = max(abs(c01[1]) * core.legendre01_prior_sigma_fraction, core.prior_min_sigma)

        if length(layout.idx_legendre) >= 2
            leg1_idx = layout.idx_legendre[2]
            x_a[leg1_idx] = c01[2]
            prior_sigma[leg1_idx] = max(abs(c01[2]) * core.legendre01_prior_sigma_fraction, core.prior_min_sigma)
        end
    end

    x_curr = copy(x_a)
    S_a_inv = _spdiag_invvar(prior_sigma)
    x_scale = copy(core.x_scale_base)

    y_curr = copy(core.fm(x_curr))
    rmse_prev = sqrt(mean((y_obs .- y_curr) .^ 2))
    dof = max(length(y_obs) - length(x_curr), 1)
    redchi2_hist = Float64[sum(((y_obs .- y_curr) ./ core.meas_sigma) .^ 2) / dof]
    dx_rel_hist = Float64[]
    λ = core.lm.lambda0
    n_acc = 0
    converged = false
    status = Int16(0)

    for _ in 1:max_outer_steps
        step = try
            lm_one_step(
                core.fm,
                x_curr,
                y_obs;
                x_a=x_a,
                S_e_inv=core.S_e_inv,
                S_a_inv=S_a_inv,
                lambda=λ,
                lambda_up=core.lm.lambda_up,
                lambda_down=core.lm.lambda_down,
                lambda_min=core.lm.lambda_min,
                lambda_max=core.lm.lambda_max,
                max_inner=core.lm.max_inner,
                jacobian_eval=core.jacobian_eval,
                x_scale=x_scale,
                lower_bounds=core.lower_bounds,
                upper_bounds=core.upper_bounds,
            )
        catch
            status = Int16(4)
            break
        end

        λ = step.lambda_next
        if !step.accepted
            stalled_conv, _ = _stalled_convergence(
                dx_rel_hist,
                redchi2_hist;
                enabled=core.conv_stall.enabled,
                window=core.conv_stall.window,
                redchi2_target=core.conv_stall.redchi2_target,
                redchi2_abs_tol=core.conv_stall.redchi2_abs_tol,
                redchi2_rel_tol=core.conv_stall.redchi2_rel_tol,
                dx_rel_tol=core.conv_stall.dx_rel_tol,
            )
            if stalled_conv
                converged = true
                status = Int16(1)
            else
                status = Int16(2)
            end
            break
        end

        n_acc += 1
        x_prev = x_curr
        x_curr = step.x_next
        copyto!(y_curr, step.y_next)
        rmse_curr = sqrt(mean((y_obs .- y_curr) .^ 2))

        dx_rel = norm(step.dx) / max(norm(x_prev), eps(Float64))
        push!(dx_rel_hist, dx_rel)
        rmse_abs_change = abs(rmse_curr - rmse_prev)
        rmse_rel_change = rmse_abs_change / max(abs(rmse_prev), eps(Float64))
        rmse_prev = rmse_curr
        push!(redchi2_hist, sum(((y_obs .- y_curr) ./ core.meas_sigma) .^ 2) / dof)

        if dx_rel < core.conv_dx_rel_tol ||
           rmse_rel_change < core.conv_rmse_rel_tol ||
           rmse_abs_change < core.conv_rmse_abs_tol
            converged = true
            status = Int16(1)
            break
        end
    end

    if status == 0 && !converged
        status = Int16(0)  # valid but not converged within max_outer_steps
    end

    x_out .= x_curr
    resid = y_obs .- y_curr
    rmse = sqrt(mean(resid .^ 2))
    rchi2 = sum((resid ./ core.meas_sigma) .^ 2) / dof
    obj = _cost_with_prior(y_obs, y_curr, x_curr, x_a, core.S_e_inv, S_a_inv)

    return (converged = converged, status = status, n_steps = n_acc, rmse = rmse, reduced_chi2 = rchi2, objective = obj)
end

function run_orbit(core, pace_path::AbstractString, config_path::AbstractString)
    cfg = core.cfg
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
    spectrum_var = String(get(pace_cfg, "spectrum_var", "radiance_red"))

    ds = Dataset(pace_path)
    haskey(ds, wavelength_var) || error("Missing wavelength variable '$wavelength_var' in $pace_path")
    haskey(ds, spectrum_var) || error("Missing spectrum variable '$spectrum_var' in $pace_path")
    v_spec = ds[spectrum_var]
    axes = find_axis_indices(v_spec, wavelength_var)

    λ_src = collect(Float64.(ds[wavelength_var][:]))
    W_interp, perm = make_linear_resampler(λ_src, core.ctx.λ)

    n_pix = axes.n_pix
    n_scan = axes.n_scan
    p_start = Int(get(batch_cfg, "pixel_start", 1))
    s_start = Int(get(batch_cfg, "scan_start", 1))
    p_end_cfg = Int(get(batch_cfg, "pixel_end", 0))
    s_end_cfg = Int(get(batch_cfg, "scan_end", 0))
    p_end = p_end_cfg > 0 ? min(p_end_cfg, n_pix) : n_pix
    s_end = s_end_cfg > 0 ? min(s_end_cfg, n_scan) : n_scan
    p_start = clamp(p_start, 1, p_end)
    s_start = clamp(s_start, 1, s_end)
    pixel_range = p_start:p_end
    scan_range = s_start:s_end
    max_outer_steps = Int(get(batch_cfg, "max_outer_steps", core.lm.max_outer_default))
    dark_filter_enabled = Bool(get(batch_cfg, "dark_filter_enabled", false))
    dark_max_radiance = Float64(get(batch_cfg, "dark_max_radiance", 20.0))
    ocean_filter_enabled = Bool(get(batch_cfg, "ocean_filter_enabled", false))
    watermask_var = String(get(batch_cfg, "watermask_var", "watermask"))
    ocean_mask_values_raw = get(batch_cfg, "ocean_mask_values", Any[1])
    ocean_mask_values = Set{Int}(Int(v) for v in ocean_mask_values_raw)
    isempty(ocean_mask_values) && error("batch_fit.ocean_mask_values must contain at least one value")

    lat = read_geo_2d(ds, "latitude", n_pix, n_scan)
    lon = read_geo_2d(ds, "longitude", n_pix, n_scan)
    watermask = if haskey(ds, watermask_var)
        wm = ds[watermask_var]
        ndims(wm) == 2 || error("Expected 2D watermask variable '$watermask_var', got ndims=$(ndims(wm))")
        d = collect(String.(dimnames(wm)))
        raw = wm[:, :]
        if d == ["pixels", "scans"]
            raw
        elseif d == ["scans", "pixels"]
            permutedims(raw, (2, 1))
        else
            error("Unsupported watermask dims for '$watermask_var': $d")
        end
    else
        if ocean_filter_enabled
            error("Ocean filter enabled but variable '$watermask_var' is missing in $pace_path")
        end
        fill(missing, n_pix, n_scan)
    end

    output_path = make_output_path(pace_path, cfg)
    ds_out = create_output_dataset(
        output_path,
        length(pixel_range),
        length(scan_range),
        core.state_names,
        pace_path,
        config_path,
        pixel_range,
        scan_range,
    )

    # Fill static geolocation and defaults.
    ds_out["latitude"][:, :] = lat[pixel_range, scan_range]
    ds_out["longitude"][:, :] = lon[pixel_range, scan_range]

    y_sorted = zeros(Float64, length(perm))
    y_obs = zeros(Float64, length(core.ctx.λ))
    x_tmp = zeros(Float64, core.layout.n_state)

    println("Running swath retrieval for: ", pace_path)
    println("  pixels: ", first(pixel_range), ":", last(pixel_range), " (", length(pixel_range), ")")
    println("  scans:  ", first(scan_range), ":", last(scan_range), " (", length(scan_range), ")")
    println("  n_state: ", core.layout.n_state, "  n_meas: ", length(y_obs))
    println("  dark filter: ", dark_filter_enabled, " (max radiance <= ", dark_max_radiance, ")")
    println("  ocean filter: ", ocean_filter_enabled, " (", watermask_var, " in ", collect(ocean_mask_values), ")")

    for (j_scan_out, j_scan_src) in enumerate(scan_range)
        inds = Any[Colon() for _ in 1:3]
        inds[axes.i_scan] = j_scan_src
        slab = v_spec[inds...]

        slab_is_pix_band = size(slab) == (axes.n_pix, axes.n_band)
        slab_is_band_pix = size(slab) == (axes.n_band, axes.n_pix)
        (slab_is_pix_band || slab_is_band_pix) ||
            error("Unexpected slab size $(size(slab)) for scan=$j_scan_src")

        state_scan = fill(Float32(NaN), length(pixel_range), core.layout.n_state)
        conv_scan = fill(UInt8(0), length(pixel_range))
        status_scan = fill(Int16(3), length(pixel_range))
        steps_scan = fill(Int16(0), length(pixel_range))
        rmse_scan = fill(Float32(NaN), length(pixel_range))
        rchi2_scan = fill(Float32(NaN), length(pixel_range))
        obj_scan = fill(Float32(NaN), length(pixel_range))
        sif1_scan = fill(Float32(NaN), length(pixel_range))
        dark_scan = fill(UInt8(0), length(pixel_range))
        ocean_scan = fill(UInt8(0), length(pixel_range))

        for (i_pix_out, i_pix_src) in enumerate(pixel_range)
            spec_raw = slab_is_pix_band ? view(slab, i_pix_src, :) : view(slab, :, i_pix_src)
            ok = copy_sorted_spectrum!(y_sorted, spec_raw, perm)
            if !ok
                status_scan[i_pix_out] = Int16(3)
                continue
            end

            wm_val = watermask[i_pix_src, j_scan_src]
            is_ocean = !ismissing(wm_val) && (Int(wm_val) in ocean_mask_values)
            ocean_scan[i_pix_out] = is_ocean ? UInt8(1) : UInt8(0)
            if ocean_filter_enabled && !is_ocean
                status_scan[i_pix_out] = Int16(6)
                continue
            end

            mul!(y_obs, W_interp, y_sorted)
            is_dark = maximum(y_obs) <= dark_max_radiance
            dark_scan[i_pix_out] = is_dark ? UInt8(1) : UInt8(0)
            if dark_filter_enabled && !is_dark
                status_scan[i_pix_out] = Int16(5)
                continue
            end

            stats = run_one_retrieval!(x_tmp, core, y_obs, max_outer_steps)

            state_scan[i_pix_out, :] .= Float32.(x_tmp)
            conv_scan[i_pix_out] = stats.converged ? UInt8(1) : UInt8(0)
            status_scan[i_pix_out] = stats.status
            steps_scan[i_pix_out] = Int16(stats.n_steps)
            rmse_scan[i_pix_out] = Float32(stats.rmse)
            rchi2_scan[i_pix_out] = Float32(stats.reduced_chi2)
            obj_scan[i_pix_out] = Float32(stats.objective)
            if length(core.layout.idx_sif) >= 1
                sif1_scan[i_pix_out] = Float32(x_tmp[first(core.layout.idx_sif)])
            end
        end

        ds_out["x_hat"][:, j_scan_out, :] = state_scan
        ds_out["converged"][:, j_scan_out] = conv_scan
        ds_out["status_code"][:, j_scan_out] = status_scan
        ds_out["n_steps"][:, j_scan_out] = steps_scan
        ds_out["rmse"][:, j_scan_out] = rmse_scan
        ds_out["reduced_chi2"][:, j_scan_out] = rchi2_scan
        ds_out["objective"][:, j_scan_out] = obj_scan
        ds_out["sif_ev1"][:, j_scan_out] = sif1_scan
        ds_out["is_dark"][:, j_scan_out] = dark_scan
        ds_out["is_ocean"][:, j_scan_out] = ocean_scan

        n_conv = count(==(UInt8(1)), conv_scan)
        n_dark = count(==(UInt8(1)), dark_scan)
        n_ocean = count(==(UInt8(1)), ocean_scan)
        println(
            "  scan ", j_scan_src,
            " -> ocean ", n_ocean, "/", length(pixel_range),
            " | dark ", n_dark, "/", length(pixel_range),
            " | converged ", n_conv, "/", length(pixel_range),
        )
    end

    close(ds_out)
    close(ds)
    println("Saved swath retrieval to: ", output_path)
    return output_path
end

function main_batch()
    config_path = get(
        ENV,
        "PACE_MWE_CONFIG",
        joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE.toml"),
    )
    core = build_retrieval_core(config_path)
    files = resolve_orbit_files(core.cfg)

    println("Batch retrieval pipeline")
    println("  config: ", config_path)
    println("  n_orbits: ", length(files))
    for f in files
        run_orbit(core, f, config_path)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_batch()
end
