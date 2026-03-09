#!/usr/bin/env julia
# Add known SIF radiance to real PACE data, run retrieval, and check if the
# pseudo SIF is recovered. Prior sigma, x_a, and LM iteration match Fit_toy_forward_model.jl.

using TOML
using LinearAlgebra
using SparseArrays
using Plots

const PSEUDO_DIR = @__DIR__
const DEMO_DIR = joinpath(PSEUDO_DIR, "..")
include(joinpath(DEMO_DIR, "Fit_toy_forward_model.jl"))

const MWEF = SimplePACEXSecFitMWEFunctions

"""
    add_sif_radiance_to_spectrum(y_obs, ctx, sif_coeff_add)

Add convolved SIF radiance to an observed spectrum: y_mod = y_obs + K * (sif_basis_hres * coeff).
Returns y_obs + SIF_at_OCI (both on ctx.λ).
"""
function add_sif_radiance_to_spectrum(
    y_obs::AbstractVector{<:Real},
    ctx,
    sif_coeff_add::AbstractVector{<:Real},
)
    length(y_obs) == length(ctx.λ) || error("y_obs length must match ctx.λ")
    length(sif_coeff_add) == size(ctx.sif_basis_hres, 2) || error("sif_coeff_add length must match SIF basis columns")
    sif_hres = ctx.sif_basis_hres * sif_coeff_add
    sif_at_oci = ctx.kernel_rsr_out * sif_hres
    return y_obs .+ sif_at_oci
end

function main_sif_addition()
    pseudo_cfg_path = joinpath(PSEUDO_DIR, "pseudo_measurement_config.toml")
    pseudo_cfg = TOML.parsefile(pseudo_cfg_path)
    sif_add_cfg = get(pseudo_cfg, "sif_addition", Dict{String, Any}())
    main_config = normpath(joinpath(PSEUDO_DIR, get(get(pseudo_cfg, "pseudo_measurement", Dict()), "main_config_path", "../Simple_PACE_xSecFit_MWE_zcheVer.toml")))
    pace_file = get(sif_add_cfg, "pace_file", "sample_granule_20240830T131442_new_chl.nc")
    pixel_index = Int(get(sif_add_cfg, "pixel_index", 1260))
    scan_index = Int(get(sif_add_cfg, "scan_index", 1700))
    sif_coeff_add = Float64.(get(sif_add_cfg, "sif_coeff_add", [0.5]))
    max_iter = Int(get(sif_add_cfg, "max_iter", 50))

    println("SIF_addition: add known SIF to PACE spectrum and retrieve")
    println("  main_config: ", main_config)
    println("  pace_file: ", pace_file, "  pixel: ", pixel_index, "  scan: ", scan_index)
    println("  SIF coeff to add (truth): ", sif_coeff_add)

    cfg = TOML.parsefile(main_config)
    ctx = MWEF.prepare_mwe_inputs(main_config)
    data_cfg = get(cfg, "data", Dict{String, Any}())
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    base_dir = get(data_cfg, "base_dir", joinpath(DEMO_DIR, "..", "Files_in_use"))
    pace_path = isabspath(pace_file) ? pace_file : joinpath(base_dir, pace_file)
    wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
    spectrum_var = String(get(pace_cfg, "spectrum_var", "radiance_red"))

    y_obs, _ = MWEF.load_pace_spectrum_on_grid(
        pace_path,
        ctx.λ;
        pixel_idx=pixel_index,
        scan_idx=scan_index,
        wavelength_var=wavelength_var,
        spectrum_var=spectrum_var,
    )
    y_obs = Float64.(y_obs)
    println("  Loaded spectrum: length ", length(y_obs), "  mean radiance ", mean(y_obs))

    length(sif_coeff_add) == size(ctx.sif_basis_hres, 2) || error(
        "sif_coeff_add length ($(length(sif_coeff_add))) must match SIF basis columns ($(size(ctx.sif_basis_hres, 2)))"
    )

    y_modified = add_sif_radiance_to_spectrum(y_obs, ctx, sif_coeff_add)
    println("  Modified spectrum: mean radiance ", mean(y_modified))

    # Run retrieval on y_modified using same prior and iteration as Fit_toy_forward_model.jl
    fit_cfg = get(cfg, "fit", Dict{String, Any}())
    data_cfg = get(cfg, "data", Dict{String, Any}())
    solar_file = get(data_cfg, "solar_file", "solar_merged_20240731_600_33300_100.out")
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)
    solar_hres, _ = MWEF.load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
    solar_hres = Float64.(solar_hres)
    n_legendre = Int(get(fit_cfg, "n_legendre", 3))
    model_variant = Symbol(get(fit_cfg, "model_variant", "standard"))
    use_hybrid_jacobian = Bool(get(fit_cfg, "use_hybrid_jacobian", false))
    preallocate_forward = Bool(get(fit_cfg, "preallocate_forward", true))

    fm = make_forward_model_simple(
        ctx,
        solar_hres;
        n_legendre=n_legendre,
        preallocate_float64=preallocate_forward,
        model_variant=model_variant,
    )
    layout = state_layout_simple(ctx; n_legendre=n_legendre)
    x0 = initial_state_simple(ctx; n_legendre=n_legendre, T=Float64, model_variant=model_variant)
    jacobian_eval = if use_hybrid_jacobian
        make_hybrid_jacobian_evaluator(fm, ctx, solar_hres, layout; n_legendre=n_legendre)
    else
        make_jacobian_evaluator(fm, x0; use_preallocated=false)
    end

    # Prior and x_a: same as Fit_toy_forward_model.jl
    prior_sigma_default = Float64(get(fit_cfg, "prior_sigma_default", 1e30))
    prior_min_sigma = Float64(get(fit_cfg, "prior_min_sigma", 1e-3))
    use_legendre01_prior = Bool(get(fit_cfg, "use_legendre01_prior", true))
    legendre01_prior_sigma_fraction = Float64(get(fit_cfg, "legendre01_prior_sigma_fraction", 0.2))
    use_legendre_higher_prior = Bool(get(fit_cfg, "use_legendre_higher_prior", true))
    legendre_higher_sigma = Float64(get(fit_cfg, "legendre_higher_sigma", 1.0))
    p_prior_hpa = Float64(get(fit_cfg, "p_prior_hpa", 700.0))
    p_sigma_hpa = Float64(get(fit_cfg, "p_sigma_hpa", 200.0))
    t_prior_k = Float64(get(fit_cfg, "t_prior_k", 280.0))
    t_sigma_k = Float64(get(fit_cfg, "t_sigma_k", 20.0))
    vcd_o2_sigma = Float64(get(fit_cfg, "vcd_o2_sigma", 1e23))
    vcd_h2o_sigma = Float64(get(fit_cfg, "vcd_h2o_sigma", 3e22))
    use_vcd_slope_prior = Bool(get(fit_cfg, "use_vcd_slope_prior", true))
    vcd_slope_prior_sigma_factor = Float64(get(fit_cfg, "vcd_slope_prior_sigma_factor", 1.0))
    sif_sigma = Float64(get(fit_cfg, "sif_sigma", 1e12))
    use_pt_constraints = Bool(get(fit_cfg, "use_pt_constraints", true))
    pt_constraint_sigma_mult = Float64(get(fit_cfg, "pt_constraint_sigma_mult", 3.0))
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
    lm_lambda_up = Float64(get(fit_cfg, "lm_lambda_up", 2.0))
    lm_lambda_down = Float64(get(fit_cfg, "lm_lambda_down", 0.7))
    lm_lambda_min = Float64(get(fit_cfg, "lm_lambda_min", 1e-8))
    lm_lambda_max = Float64(get(fit_cfg, "lm_lambda_max", 1e8))
    lm_max_inner = Int(get(fit_cfg, "lm_max_inner", 8))
    use_band_snr = Bool(get(fit_cfg, "use_band_snr", true))
    meas_sigma = Float64(get(fit_cfg, "meas_sigma", 0.01))
    n_plot_steps = Int(get(fit_cfg, "n_plot_steps", 20))

    x_a = copy(x0)
    prior_sigma = fill(prior_sigma_default, length(x0))
    if use_legendre01_prior && length(layout.idx_legendre) >= 1
        leg0_idx = first(layout.idx_legendre)
        y_base = fm(x0)
        ratio = y_modified ./ max.(abs.(y_base), eps(Float64))
        z = _normalized_grid(ctx.λ)
        A01 = hcat(ones(length(z)), z)
        w = y_modified .- minimum(y_modified)
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
    x_a[layout.idx_sif] .= 0.0
    prior_sigma[layout.idx_sif] .= max(sif_sigma, prior_min_sigma)

    if hasproperty(ctx, :sif_prior_cov) && !isnothing(ctx.sif_prior_cov) &&
       length(layout.idx_sif) == size(ctx.sif_prior_cov, 1)
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

    # Iteration: same as Fit_toy_forward_model.jl (LM loop with convergence and stalled check)
    x_curr = copy(x_a)
    y_curr = copy(fm(x_curr))
    S_e_inv = if !use_band_snr || isnothing(ctx.band_snr_coeffs)
        spdiagm(0 => fill(1.0 / (meas_sigma^2), length(y_modified)))
    else
        make_Se_inv_from_snr(y_curr, ctx.band_snr_coeffs)
    end
    dof = max(length(y_modified) - layout.n_state, 1)
    dx_rel_series = Float64[]
    redchi2_series = Float64[]
    λ = lm_lambda0
    n_steps = max(n_plot_steps, 0)
    converged = false
    for istep in 1:n_steps
        println("LM step ", istep)
        x_prev = copy(x_curr)
        rmse_prev = sqrt(mean((y_modified .- y_curr) .^ 2))
        step = try
            lm_one_step(
                fm,
                x_curr,
                y_modified;
                x_a=x_a,
                S_a_inv=S_a_inv,
                lambda=λ,
                lambda_up=lm_lambda_up,
                lambda_down=lm_lambda_down,
                lambda_min=lm_lambda_min,
                lambda_max=lm_lambda_max,
                max_inner=lm_max_inner,
                jacobian_eval=jacobian_eval,
                x_scale=x_scale,
                lower_bounds=lower_bounds,
                upper_bounds=upper_bounds,
                use_band_snr=use_band_snr,
                band_snr_coeffs=ctx.band_snr_coeffs,
                meas_sigma=meas_sigma,
            )
        catch err
            println("  LM step failed: ", err)
            break
        end
        λ = step.lambda_next
        if !step.accepted
            stalled_conv, _ = _stalled_convergence(
                dx_rel_series,
                redchi2_series;
                enabled=conv_stall_enable,
                window=conv_stall_window,
                redchi2_target=conv_stall_redchi2_target,
                redchi2_abs_tol=conv_stall_redchi2_abs_tol,
                redchi2_rel_tol=conv_stall_redchi2_rel_tol,
                dx_rel_tol=conv_stall_dx_rel_tol,
            )
            if stalled_conv
                converged = true
            end
            break
        end
        x_curr = step.x_next
        y_curr = copy(step.y_next)
        chi2_curr = step.chi2_next
        push!(redchi2_series, chi2_curr / dof)
        dx_rel = norm(step.dx) / max(norm(x_prev), eps(Float64))
        push!(dx_rel_series, dx_rel)
        rmse_curr = sqrt(mean((y_modified .- y_curr) .^ 2))
        rmse_abs_change = abs(rmse_curr - rmse_prev)
        rmse_rel_change = rmse_abs_change / max(abs(rmse_prev), eps(Float64))
        if dx_rel < conv_dx_rel_tol ||
           rmse_rel_change < conv_rmse_rel_tol ||
           rmse_abs_change < conv_rmse_abs_tol
            converged = true
            break
        end
    end

    sif_retrieved = x_curr[layout.idx_sif]
    sif_hres      = ctx.sif_basis_hres * sif_retrieved   # reconstructed SIF at high resolution
    sif_added_hres = ctx.sif_basis_hres * sif_coeff_add
    y_reconstructed = fm(x_curr)
    rmse = sqrt(mean((y_obs - y_reconstructed) .^ 2))

    println("\n--- SIF addition experiment results ---")
    println("  Added SIF coeff (truth): ", sif_coeff_add)
    println("  Retrieved SIF coeff:     ", sif_retrieved)
    err = sif_retrieved .- sif_coeff_add
    println("  Difference:             ", err)
    println("  Relative error (ev1):    ", length(sif_coeff_add) >= 1 && abs(sif_coeff_add[1]) > 0 ? err[1] / sif_coeff_add[1] : "N/A")

    # plot the fit
    p_fit = plot(
        ctx.λ,
        y_obs;
        label="Observation",
        lw=2,
        color=:black,
    )
    plot!(p_fit, ctx.λ, y_modified; label="Modified spectrum", lw=1.8, color=:blue)
    plot!(p_fit, ctx.λ, y_reconstructed; label="Reconstructed spectrum", lw=1.8, color=:green)
    savefig(p_fit, "sif_addition_fit.png")
    println("  Saved plot to sif_addition_fit.png")
    
    # added SIF (truth) vs. retrieved SIF
    p_sif = plot(
        ctx.λ_hres,
        sif_added_hres;
        label="Added SIF (truth)",
        lw=2,
        color=:black,
    )
    plot!(p_sif, ctx.λ_hres, sif_hres; label="Retrieved SIF", lw=1.8, color=:red)
    savefig(p_sif, "sif_addition_sif_signal.png")
    println("  Saved plot to sif_addition_sif_signal.png")

    # summarize retrieved x
    x_final = x_curr
    state_names = state_names_simple(ctx; n_legendre=n_legendre)
    println("  Final state vector:", x_final)
    for (i, name) in zip(eachindex(x_final), state_names)
        println("    [$i] $name: $(x_final[i])")
    end

    return (sif_added=sif_coeff_add, sif_retrieved=sif_retrieved, y_obs=y_obs, y_modified=y_modified, x_final=x_curr)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_sif_addition()
end
