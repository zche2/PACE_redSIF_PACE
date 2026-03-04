#!/usr/bin/env julia
# Compare fit with SIF vs without SIF (SIF prior pinned to zero via tight prior_sigma).
# Similar structure to Benchmark_interpolation_modes.jl; produces comparison plots of final fit and residuals.

using TOML
using LinearAlgebra
using Statistics
using SparseArrays
using Plots

include(joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE_Functions.jl"))
using .SimplePACEXSecFitMWEFunctions

include(joinpath(@__DIR__, "toy_forward_model.jl"))

# Load Fit_toy_forward_model.jl without executing main()
let
    fit_path = joinpath(@__DIR__, "Fit_toy_forward_model.jl")
    src = read(fit_path, String)
    src = replace(src, r"\nmain\(\)\s*$" => "\n")
    # Also skip "if abspath(PROGRAM_FILE) == @__FILE__" block so we don't run main when including
    src = replace(src, r"\nif abspath\(PROGRAM_FILE\)\s*==\s*@__FILE__\s*\n\s*main\(\)\s*\nend\s*$" => "\n")
    Base.include_string(Main, src, fit_path)
end

"""
    _build_setup(config_path, with_sif::Bool)

Build shared context and fit setup. When with_sif=true use normal SIF prior (weak);
when with_sif=false use very tight SIF prior (prior_sigma[layout.idx_sif] = 1e-12) so SIF stays at 0.
Returns NamedTuple with ctx, fm, layout, x_a, S_a_inv, jacobian_eval, y_obs, lm, x_scale, lower_bounds, upper_bounds, use_band_snr.
"""
function _build_setup(config_path::AbstractString, with_sif::Bool)
    cfg = TOML.parsefile(config_path)
    fit_cfg = get(cfg, "fit", Dict{String, Any}())
    data_cfg = get(cfg, "data", Dict{String, Any}())
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())

    ctx = SimplePACEXSecFitMWEFunctions.prepare_mwe_inputs(config_path)
    state_float_type = SimplePACEXSecFitMWEFunctions.parse_float_type(cfg)
    n_legendre = Int(get(fit_cfg, "n_legendre", 2))
    model_variant = Symbol(get(fit_cfg, "model_variant", "standard"))
    preallocate_forward = Bool(get(fit_cfg, "preallocate_forward", true))
    preallocate_ad_forward = Bool(get(fit_cfg, "preallocate_ad_forward", false))
    use_hybrid_jacobian = Bool(get(fit_cfg, "use_hybrid_jacobian", false))
    use_band_snr = Bool(get(fit_cfg, "use_band_snr", true))
    lm_lambda0 = Float64(get(fit_cfg, "lm_lambda0", 1.0))
    lm_lambda_up = Float64(get(fit_cfg, "lm_lambda_up", 2.0))
    lm_lambda_down = Float64(get(fit_cfg, "lm_lambda_down", 0.7))
    lm_lambda_min = Float64(get(fit_cfg, "lm_lambda_min", 1e-8))
    lm_lambda_max = Float64(get(fit_cfg, "lm_lambda_max", 1e8))
    lm_max_inner = Int(get(fit_cfg, "lm_max_inner", 8))
    prior_min_sigma = Float64(get(fit_cfg, "prior_min_sigma", 1e-3))
    prior_sigma_default = Float64(get(fit_cfg, "prior_sigma_default", 1e30))
    use_legendre01_prior = Bool(get(fit_cfg, "use_legendre01_prior", true))
    legendre01_prior_sigma_fraction = Float64(get(fit_cfg, "legendre01_prior_sigma_fraction", 0.2))
    use_legendre_higher_prior = Bool(get(fit_cfg, "use_legendre_higher_prior", true))
    legendre_higher_sigma = Float64(get(fit_cfg, "legendre_higher_sigma", 1.0))
    use_vcd_slope_prior = Bool(get(fit_cfg, "use_vcd_slope_prior", true))
    vcd_slope_prior_sigma_factor = Float64(get(fit_cfg, "vcd_slope_prior_sigma_factor", 1.0))
    p_prior_hpa = Float64(get(fit_cfg, "p_prior_hpa", 700.0))
    p_sigma_hpa = Float64(get(fit_cfg, "p_sigma_hpa", 200.0))
    t_prior_k = Float64(get(fit_cfg, "t_prior_k", 280.0))
    t_sigma_k = Float64(get(fit_cfg, "t_sigma_k", 20.0))
    vcd_o2_sigma = Float64(get(fit_cfg, "vcd_o2_sigma", 1e23))
    vcd_h2o_sigma = Float64(get(fit_cfg, "vcd_h2o_sigma", 3e22))
    use_pt_constraints = Bool(get(fit_cfg, "use_pt_constraints", true))
    pt_constraint_sigma_mult = Float64(get(fit_cfg, "pt_constraint_sigma_mult", 3.0))
    meas_sigma = Float64(get(fit_cfg, "meas_sigma", 0.01))
    conv_dx_rel_tol = Float64(get(fit_cfg, "conv_dx_rel_tol", 1e-6))
    conv_rmse_rel_tol = Float64(get(fit_cfg, "conv_rmse_rel_tol", 1e-6))
    conv_rmse_abs_tol = Float64(get(fit_cfg, "conv_rmse_abs_tol", 1e-6))
    conv_stall_enable = Bool(get(fit_cfg, "conv_stall_enable", true))
    conv_stall_window = Int(get(fit_cfg, "conv_stall_window", 3))
    conv_stall_redchi2_target = Float64(get(fit_cfg, "conv_stall_redchi2_target", 5.0))
    conv_stall_redchi2_abs_tol = Float64(get(fit_cfg, "conv_stall_redchi2_abs_tol", 0.1))
    conv_stall_redchi2_rel_tol = Float64(get(fit_cfg, "conv_stall_redchi2_rel_tol", 0.03))
    conv_stall_dx_rel_tol = Float64(get(fit_cfg, "conv_stall_dx_rel_tol", 5e-3))

    solar_file = get(data_cfg, "solar_file", "solar_merged_20240731_600_33300_100.out")
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)
    solar_hres, _ = SimplePACEXSecFitMWEFunctions.load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
    solar_hres = state_float_type.(solar_hres)

    pace_file = get(pace_cfg, "pace_file", "sample_granule_20240830T131442_new_chl.nc")
    pace_path = isabspath(pace_file) ? pace_file : joinpath(ctx.paths.base_dir, pace_file)
    wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
    spectrum_var = String(get(pace_cfg, "spectrum_var", "radiance_red"))
    pixel_idx = Int(get(pace_cfg, "pixel_index", 600))
    scan_idx = Int(get(pace_cfg, "scan_index", 800))
    y_obs, _ = SimplePACEXSecFitMWEFunctions.load_pace_spectrum_on_grid(
        pace_path,
        ctx.λ;
        pixel_idx=pixel_idx,
        scan_idx=scan_idx,
        wavelength_var=wavelength_var,
        spectrum_var=spectrum_var,
    )
    y_obs = Float64.(y_obs)

    fm = make_forward_model_simple(
        ctx,
        solar_hres;
        n_legendre=n_legendre,
        model_variant=model_variant,
        preallocate_float64=preallocate_forward && state_float_type == Float64,
        preallocate_float32=preallocate_forward && state_float_type == Float32,
        preallocate_other_types=preallocate_ad_forward,
    )
    layout = state_layout_simple(ctx; n_legendre=n_legendre)
    x0 = initial_state_simple(ctx; n_legendre=n_legendre, T=Float64, model_variant=model_variant)
    jacobian_eval = make_jacobian_evaluator(fm, Float64.(x0); use_preallocated=false)

    # Prior setup per element (matching Fit_toy_forward_model.jl)
    x_a = copy(Float64.(x0))
    prior_sigma = fill(prior_sigma_default, length(x0))

    if use_legendre01_prior && length(layout.idx_legendre) >= 1
        leg0_idx = first(layout.idx_legendre)
        y_base = fm(x0)
        ratio = y_obs ./ max.(abs.(y_base), eps(Float64))
        z = _normalized_grid(ctx.λ)
        A01 = hcat(ones(length(z)), z)
        w = y_obs .- minimum(y_obs)
        w .+= max(maximum(w), 1.0) * 1e-6
        s = sqrt.(w ./ maximum(w))
        c01 = (A01 .* s) \ (ratio .* s)
        x_a[leg0_idx] = c01[1]
        prior_sigma[leg0_idx] = max(abs(c01[1]) * legendre01_prior_sigma_fraction, prior_min_sigma)
        if length(layout.idx_legendre) >= 2
            leg1_idx = layout.idx_legendre[2]
            x_a[leg1_idx] = c01[2]
            prior_sigma[leg1_idx] = max(abs(c01[2]) * legendre01_prior_sigma_fraction, prior_min_sigma)
        end
    end

    if use_legendre_higher_prior && length(layout.idx_legendre) >= 3
        for j in 3:length(layout.idx_legendre)
            idx = layout.idx_legendre[j]
            x_a[idx] = 0.0
            prior_sigma[idx] = max(legendre_higher_sigma, prior_min_sigma)
        end
    end

    x_a[layout.idx_p_o2_hpa] = p_prior_hpa
    x_a[layout.idx_p_h2o_hpa] = p_prior_hpa
    x_a[layout.idx_t_o2_k] = t_prior_k
    x_a[layout.idx_t_h2o_k] = t_prior_k
    x_a[layout.idx_vcd_o2_intercept] = x0[layout.idx_vcd_o2_intercept]
    x_a[layout.idx_vcd_h2o_intercept] = x0[layout.idx_vcd_h2o_intercept]
    x_a[layout.idx_vcd_o2_sif] = x0[layout.idx_vcd_o2_sif]
    x_a[layout.idx_vcd_h2o_sif] = x0[layout.idx_vcd_h2o_sif]
    prior_sigma[layout.idx_p_o2_hpa] = p_sigma_hpa
    prior_sigma[layout.idx_p_h2o_hpa] = p_sigma_hpa
    prior_sigma[layout.idx_t_o2_k] = t_sigma_k
    prior_sigma[layout.idx_t_h2o_k] = t_sigma_k
    prior_sigma[layout.idx_vcd_o2_intercept] = vcd_o2_sigma
    prior_sigma[layout.idx_vcd_h2o_intercept] = vcd_h2o_sigma
    prior_sigma[layout.idx_vcd_o2_sif] = vcd_o2_sigma
    prior_sigma[layout.idx_vcd_h2o_sif] = vcd_h2o_sigma

    if use_vcd_slope_prior
        x_a[layout.idx_vcd_o2_slope] = 0.0
        x_a[layout.idx_vcd_h2o_slope] = 0.0
        prior_sigma[layout.idx_vcd_o2_slope] = max(vcd_o2_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
        prior_sigma[layout.idx_vcd_h2o_slope] = max(vcd_h2o_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
    end

    sif_sigma = with_sif ? Float64(get(fit_cfg, "sif_sigma", 1e12)) : 1e-12
    x_a[layout.idx_sif] .= 0.0
    prior_sigma[layout.idx_sif] .= max(sif_sigma, prior_min_sigma)

    if hasproperty(ctx, :sif_prior_cov) && !isnothing(ctx.sif_prior_cov) && with_sif && length(layout.idx_sif) == size(ctx.sif_prior_cov, 1)
        σ = collect(Float64.(prior_sigma))
        @. σ = clamp(abs(σ), 1e-12, 1e100)
        S_a_inv_dense = Matrix(Diagonal(@. 1.0 / (σ^2)))
        S_a_inv_dense[layout.idx_sif, layout.idx_sif] .= inv(ctx.sif_prior_cov)
        S_a_inv = S_a_inv_dense
    else
        S_a_inv = _spdiag_invvar(prior_sigma)
    end

    x_scale = ones(Float64, length(x0))
    x_scale[layout.idx_vcd_o2_intercept] = vcd_o2_sigma
    x_scale[layout.idx_vcd_o2_slope] = max(vcd_o2_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
    x_scale[layout.idx_vcd_h2o_intercept] = vcd_h2o_sigma
    x_scale[layout.idx_vcd_h2o_slope] = max(vcd_h2o_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
    x_scale[layout.idx_vcd_o2_sif] = vcd_o2_sigma
    x_scale[layout.idx_vcd_h2o_sif] = vcd_h2o_sigma
    x_scale[layout.idx_p_o2_hpa] = p_sigma_hpa
    x_scale[layout.idx_p_h2o_hpa] = p_sigma_hpa
    x_scale[layout.idx_t_o2_k] = t_sigma_k
    x_scale[layout.idx_t_h2o_k] = t_sigma_k
    x_scale[layout.idx_sif] .= 1.0
    x_scale[layout.idx_legendre] .= 1.0

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

    lm_params = (
        lambda0=lm_lambda0,
        lambda_up=lm_lambda_up,
        lambda_down=lm_lambda_down,
        lambda_min=lm_lambda_min,
        lambda_max=lm_lambda_max,
        max_inner=lm_max_inner,
    )

    # S_e_inv and dof for convergence checks (redchi2, stall)
    y_init = fm(x_a)
    S_e_inv = if use_band_snr && !isnothing(ctx.band_snr_coeffs)
        make_Se_inv_from_snr(y_init, ctx.band_snr_coeffs)
    else
        spdiagm(0 => fill(1.0 / (meas_sigma^2), length(y_obs)))
    end
    dof = max(length(y_obs) - layout.n_state, 1)
    conv = (
        dx_rel_tol=conv_dx_rel_tol,
        rmse_rel_tol=conv_rmse_rel_tol,
        rmse_abs_tol=conv_rmse_abs_tol,
        stall_enable=conv_stall_enable,
        stall_window=conv_stall_window,
        stall_redchi2_target=conv_stall_redchi2_target,
        stall_redchi2_abs_tol=conv_stall_redchi2_abs_tol,
        stall_redchi2_rel_tol=conv_stall_redchi2_rel_tol,
        stall_dx_rel_tol=conv_stall_dx_rel_tol,
    )

    return (
        ctx=ctx,
        fm=fm,
        layout=layout,
        x_a=x_a,
        S_a_inv=S_a_inv,
        S_e_inv=S_e_inv,
        jacobian_eval=jacobian_eval,
        y_obs=y_obs,
        lm=lm_params,
        x_scale=x_scale,
        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
        use_band_snr=use_band_snr,
        dof=dof,
        conv=conv,
    )
end

"""
    _run_lm(setup, max_steps::Int)

Run LM outer iterations, each with inner iterations (via lm_one_step max_inner).
Tests convergence (dx_rel, rmse change) and stalled-convergence when no step accepted;
stops early when converged. Returns (x_final, y_final, rmse_final, converged, steps_taken).
"""
function _run_lm(setup, max_steps::Int)
    c = setup.conv
    x_curr = copy(setup.x_a)
    y_curr = setup.fm(x_curr)
    λ = setup.lm.lambda0
    dof = setup.dof
    chi2_curr = dot(setup.y_obs .- y_curr, setup.S_e_inv * (setup.y_obs .- y_curr))
    rmse_series = [sqrt(mean((setup.y_obs .- y_curr) .^ 2))]
    redchi2_series = [chi2_curr / dof]
    dx_rel_series = Float64[]
    converged = false
    steps_taken = 0

    for istep in 1:max_steps
        x_prev = copy(x_curr)
        rmse_prev = rmse_series[end]
        step = try
            lm_one_step(
                setup.fm,
                x_curr,
                setup.y_obs;
                x_a=setup.x_a,
                S_a_inv=setup.S_a_inv,
                lambda=λ,
                lambda_up=setup.lm.lambda_up,
                lambda_down=setup.lm.lambda_down,
                lambda_min=setup.lm.lambda_min,
                lambda_max=setup.lm.lambda_max,
                max_inner=setup.lm.max_inner,
                jacobian_eval=setup.jacobian_eval,
                x_scale=setup.x_scale,
                lower_bounds=setup.lower_bounds,
                upper_bounds=setup.upper_bounds,
                use_band_snr=setup.use_band_snr,
                band_snr_coeffs=setup.ctx.band_snr_coeffs,
            )
        catch err
            break
        end
        λ = step.lambda_next
        steps_taken = istep

        if !step.accepted
            # LM did not accept any of the max_inner inner tries (no step improved RMSE).
            stalled_conv, stalled_msg = _stalled_convergence(
                dx_rel_series,
                redchi2_series;
                enabled=c.stall_enable,
                window=c.stall_window,
                redchi2_target=c.stall_redchi2_target,
                redchi2_abs_tol=c.stall_redchi2_abs_tol,
                redchi2_rel_tol=c.stall_redchi2_rel_tol,
                dx_rel_tol=c.stall_dx_rel_tol,
            )
            if stalled_conv
                converged = true
            end
            break
        end

        x_curr = step.x_next
        y_curr = copy(step.y_next)
        push!(rmse_series, sqrt(mean((setup.y_obs .- y_curr) .^ 2)))
        push!(redchi2_series, step.chi2_next / dof)
        dx_rel = norm(step.dx) / max(norm(x_prev), eps(Float64))
        push!(dx_rel_series, dx_rel)

        rmse_curr = rmse_series[end]
        rmse_abs_change = abs(rmse_curr - rmse_prev)
        rmse_rel_change = rmse_abs_change / max(abs(rmse_prev), eps(Float64))
        if dx_rel < c.dx_rel_tol ||
           rmse_rel_change < c.rmse_rel_tol ||
           rmse_abs_change < c.rmse_abs_tol
            converged = true
            break
        end
    end

    y_final = setup.fm(x_curr)
    rmse_final = sqrt(mean((setup.y_obs .- y_final) .^ 2))
    return (
        x_final=x_curr,
        y_final=y_final,
        rmse_final=rmse_final,
        converged=converged,
        steps_taken=steps_taken,
    )
end

function main()
    config_path = get(
        ENV,
        "PACE_MWE_CONFIG",
        joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE_zcheVer.toml"),
    )
    n_steps = parse(Int, get(ENV, "PACE_COMPARE_SIF_NSTEPS", "20"))

    println("Compare SIF fit: with SIF vs without SIF (SIF prior pinned to 0)")
    println("  config: ", config_path)
    println("  LM steps: ", n_steps)

    setup_with = _build_setup(config_path, true)
    setup_without = _build_setup(config_path, false)

    println("  Running fit WITH SIF (weak prior on SIF coeffs)...")
    res_with = _run_lm(setup_with, n_steps)
    println("    RMSE = ", res_with.rmse_final, "  steps = ", res_with.steps_taken, "/", n_steps, "  converged = ", res_with.converged)

    println("  Running fit WITHOUT SIF (SIF prior sigma = 1e-12)...")
    res_without = _run_lm(setup_without, n_steps)
    println("    RMSE = ", res_without.rmse_final, "  steps = ", res_without.steps_taken, "/", n_steps, "  converged = ", res_without.converged)

    λ = collect(Float64, setup_with.ctx.λ)
    y_obs = setup_with.y_obs

    # ---- Plot: final fit (obs + both models) ----
    p_fit = plot(
        λ,
        y_obs;
        label="Observation",
        lw=2,
        color=:black,
        xlabel="Wavelength [nm]",
        ylabel="Radiance",
        title="Fit comparison: with vs without SIF",
    )
    plot!(p_fit, λ, res_with.y_final; label="Model (with SIF), RMSE=$(round(res_with.rmse_final, digits=6))", lw=1.8, color=:blue)
    plot!(p_fit, λ, res_without.y_final; label="Model (no SIF), RMSE=$(round(res_without.rmse_final, digits=6))", lw=1.8, color=:red, ls=:dash)

    # ---- Plot: residuals ----
    res_with_resid = y_obs .- res_with.y_final
    res_without_resid = y_obs .- res_without.y_final
    p_resid = plot(
        λ,
        res_with_resid;
        label="Residual (with SIF)",
        lw=1.8,
        color=:blue,
        xlabel="Wavelength [nm]",
        ylabel="Obs − Model",
        title="Residuals",
    )
    plot!(p_resid, λ, res_without_resid; label="Residual (no SIF)", lw=1.8, color=:red, ls=:dash)
    hline!(p_resid, [0.0]; color=:black, ls=:dot, lw=1, label="")

    # ---- Plot: SIF spectral shape (from "with SIF" retrieval) ----
    sif_coeff = res_with.x_final[setup_with.layout.idx_sif]
    sif_hres = setup_with.ctx.sif_basis_hres * sif_coeff
    λ_hres = collect(Float64, setup_with.ctx.λ_hres)
    sif_max = maximum(sif_hres)
    sif_mean = mean(sif_hres)
    sif_norm = sqrt(sum(sif_hres .^ 2))
    x_lo, x_hi = extrema(λ_hres)
    y_hi = maximum(sif_hres)
    p_sif = plot(
        λ_hres,
        sif_hres;
        lw=2,
        color=:green,
        xlabel="Wavelength [nm]",
        ylabel="SIF [mW/m²/sr/nm]",
        title="Retrieved SIF spectral shape (with-SIF fit)",
        label="SIF",
    )
    # Label magnitude (top-left corner, one block)
    mag_str = "max = $(round(sif_max, digits=4))\nmean = $(round(sif_mean, digits=4))\n‖SIF‖₂ = $(round(sif_norm, digits=4))"
    annotate!(p_sif, x_lo + 0.02 * (x_hi - x_lo), 0.95 * y_hi, text(mag_str, 8, :left, :top))

    # Top two panels share x (λ); SIF panel has its own high-res λ grid
    p = plot(p_fit, p_resid, p_sif; layout=(3, 1), size=(900, 950))
    out_path = get(
        ENV,
        "PACE_COMPARE_SIF_PLOT",
        joinpath(@__DIR__, "compare_SIF_fit.png"),
    )
    savefig(p, out_path)
    println("  Saved plot: ", out_path)

    println("\nSummary:")
    println("  With SIF:    RMSE = ", res_with.rmse_final, "  (converged = ", res_with.converged, ", steps = ", res_with.steps_taken, ")")
    println("  Without SIF: RMSE = ", res_without.rmse_final, "  (converged = ", res_without.converged, ", steps = ", res_without.steps_taken, ")")
    println("  SIF improves RMSE by factor ", round(res_without.rmse_final / res_with.rmse_final, digits=4))
end

main()
