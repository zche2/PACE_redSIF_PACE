#!/usr/bin/env julia

using TOML
using Statistics
using SparseArrays

include(joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE_Functions.jl"))
using .SimplePACEXSecFitMWEFunctions

include(joinpath(@__DIR__, "toy_forward_model.jl"))

# Load helper functions from `Fit_toy_forward_model.jl` without executing
# `main()`. Done at top level to avoid world-age issues.
let
    fit_path = joinpath(@__DIR__, "Fit_toy_forward_model.jl")
    src = read(fit_path, String)
    src = replace(src, r"\nmain\(\)\s*$" => "\n")
    Base.include_string(Main, src, fit_path)
end

function _time_avg_ms(f::Function, n::Int)
    t0 = time_ns()
    for _ in 1:n
        f()
    end
    return (time_ns() - t0) / 1e6 / n
end

function _build_setup(config_path::AbstractString, mode::Symbol)
    old_mode = get(ENV, "PACE_LUT_INTERPOLATION", nothing)
    ENV["PACE_LUT_INTERPOLATION"] = String(mode)
    t_setup = @elapsed ctx = SimplePACEXSecFitMWEFunctions.prepare_mwe_inputs(config_path)
    if isnothing(old_mode)
        delete!(ENV, "PACE_LUT_INTERPOLATION")
    else
        ENV["PACE_LUT_INTERPOLATION"] = old_mode
    end

    cfg = TOML.parsefile(config_path)
    fit_cfg = get(cfg, "fit", Dict{String, Any}())
    data_cfg = get(cfg, "data", Dict{String, Any}())
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())

    n_legendre = Int(get(fit_cfg, "n_legendre", 2))
    state_float_type = SimplePACEXSecFitMWEFunctions.parse_float_type(cfg)
    preallocate_forward = Bool(get(fit_cfg, "preallocate_forward", true))
    preallocate_ad_forward = Bool(get(fit_cfg, "preallocate_ad_forward", false))
    use_hybrid_jacobian = Bool(get(fit_cfg, "use_hybrid_jacobian", false))
    lm_lambda0 = Float64(get(fit_cfg, "lm_lambda0", 1.0))
    lm_lambda_up = Float64(get(fit_cfg, "lm_lambda_up", 5.0))
    lm_lambda_down = Float64(get(fit_cfg, "lm_lambda_down", 0.7))
    lm_lambda_min = Float64(get(fit_cfg, "lm_lambda_min", 1e-8))
    lm_lambda_max = Float64(get(fit_cfg, "lm_lambda_max", 1e8))
    lm_max_inner = Int(get(fit_cfg, "lm_max_inner", 24))
    meas_sigma = Float64(get(fit_cfg, "meas_sigma", 0.01))

    solar_file = get(data_cfg, "solar_file", "solar_merged_20200720_600_33300_100.out")
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)
    solar_hres, _ = SimplePACEXSecFitMWEFunctions.load_solar_spectrum_on_grid(
        solar_path,
        ctx.λ_hres;
        header_lines=3,
    )
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
    y_obs = state_float_type.(y_obs)

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
        make_jacobian_evaluator(fm, x0; use_preallocated=false)
    end

    # Keep prior model simple for speed-focused benchmark.
    # Use finite diagonal prior widths to keep LM linear systems nonsingular.
    x_a = copy(Float64.(x0))
    prior_sigma = fill(1e2, length(x0))
    S_a_inv = _spdiag_invvar(prior_sigma)
    S_e_inv = spdiagm(0 => fill(1.0 / (meas_sigma^2), length(y_obs)))

    return (
        t_setup_s = t_setup,
        fm = fm,
        jacobian_eval = jacobian_eval,
        x0 = Float64.(x0),
        x_a = x_a,
        y_obs = Float64.(y_obs),
        S_a_inv = S_a_inv,
        S_e_inv = S_e_inv,
        lm = (
            lambda0 = lm_lambda0,
            lambda_up = lm_lambda_up,
            lambda_down = lm_lambda_down,
            lambda_min = lm_lambda_min,
            lambda_max = lm_lambda_max,
            max_inner = lm_max_inner,
        ),
    )
end

function _run_lm_core_10_steps(setup)
    x = copy(setup.x_a)
    λ = setup.lm.lambda0
    accepted = 0
    for _ in 1:10
        step = lm_one_step(
            setup.fm,
            x,
            setup.y_obs;
            x_a=setup.x_a,
            S_e_inv=setup.S_e_inv,
            S_a_inv=setup.S_a_inv,
            lambda=λ,
            lambda_up=setup.lm.lambda_up,
            lambda_down=setup.lm.lambda_down,
            lambda_min=setup.lm.lambda_min,
            lambda_max=setup.lm.lambda_max,
            max_inner=setup.lm.max_inner,
            jacobian_eval=setup.jacobian_eval,
        )
        λ = step.lambda_next
        step.accepted || break
        x = step.x_next
        accepted += 1
    end
    return accepted
end

function benchmark_mode(config_path::AbstractString, mode::Symbol)
    setup = _build_setup(config_path, mode)

    # Warmup
    setup.fm(setup.x0)
    setup.jacobian_eval(setup.x0)
    _run_lm_core_10_steps(setup)

    forward_ms = _time_avg_ms(() -> setup.fm(setup.x0), 500)
    jacobian_ms = _time_avg_ms(() -> setup.jacobian_eval(setup.x0), 60)
    lm10_ms = _time_avg_ms(() -> _run_lm_core_10_steps(setup), 5)
    accepted = _run_lm_core_10_steps(setup)

    return (
        mode = mode,
        setup_s = setup.t_setup_s,
        forward_ms = forward_ms,
        jacobian_ms = jacobian_ms,
        lm10_ms = lm10_ms,
        lm10_accepted_steps = accepted,
    )
end

function main()
    config_path = get(
        ENV,
        "PACE_MWE_CONFIG",
        joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE.toml"),
    )

    results = map(mode -> benchmark_mode(config_path, mode), (:cubic, :linear))

    println("Interpolation benchmark summary")
    for r in results
        println("  mode = ", r.mode)
        println("    setup time [s]: ", r.setup_s)
        println("    forward eval [ms/call]: ", r.forward_ms)
        println("    Jacobian eval [ms/call]: ", r.jacobian_ms)
        println("    LM core (10-step budget) [ms/run]: ", r.lm10_ms)
        println("    LM accepted steps: ", r.lm10_accepted_steps, "/10")
    end

    cubic = only(filter(r -> r.mode == :cubic, results))
    linear = only(filter(r -> r.mode == :linear, results))
    println("Speedup (linear / cubic)")
    println("  forward speedup: ", cubic.forward_ms / linear.forward_ms, "x")
    println("  Jacobian speedup: ", cubic.jacobian_ms / linear.jacobian_ms, "x")
    println("  LM core speedup: ", cubic.lm10_ms / linear.lm10_ms, "x")
end

main()
