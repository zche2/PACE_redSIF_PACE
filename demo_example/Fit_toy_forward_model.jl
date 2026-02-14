#!/usr/bin/env julia

using TOML
using LinearAlgebra
using Statistics
using ForwardDiff
using SparseArrays
using Plots

include(joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE_Functions.jl"))
using .SimplePACEXSecFitMWEFunctions

include(joinpath(@__DIR__, "toy_forward_model.jl"))

"""
    rodgers_eq59_fit(fm, x0, y_obs; kwargs...)

Simple unconstrained nonlinear least-squares:
min_x 0.5 * ||y_obs - fm(x)||^2

Rodgers Eq. 5.9 iteration (explicit prior-state form):
    x_{i+1} = x_a + (S_a^{-1} + K_i' S_ϵ^{-1} K_i)^{-1}
                    K_i' S_ϵ^{-1} [y - f(x_i) + K_i (x_i - x_a)]

For no-prior retrieval set `use_prior=false` (default), which implies `S_a^{-1}=0`.
"""
function rodgers_eq59_fit(
    fm,
    x0::AbstractVector{<:Real},
    y_obs::AbstractVector{<:Real};
    ctx=nothing,
    layout=nothing,
    x_a::AbstractVector{<:Real}=collect(Float64.(x0)),
    use_prior::Bool=false,
    prior_sigma::Union{Nothing, AbstractVector{<:Real}}=nothing,
    meas_sigma::Union{Nothing, AbstractVector{<:Real}}=nothing,
    jacobian_eval::Union{Nothing, Function}=nothing,
    max_iter::Int=12,
    rel_obj_tol::Float64=1e-8,
    rel_step_tol::Float64=1e-8,
    verbose::Bool=true,
)
    x = collect(Float64.(x0))
    x_a = collect(Float64.(x_a))
    y = fm(x)                    # f(x₀)
    r = y_obs .- y               # residual at x₀
    obj = 0.5 * dot(r, r)        # objective at x₀

    if verbose
        println("Initial objective: ", obj, "   RMSE: ", sqrt(mean(r .^ 2)))
    end

    history = NamedTuple[]
    converged = false
    message = "max_iter reached"

    n_state = length(x0)
    if use_prior
        if isnothing(prior_sigma)
            error("use_prior=true requires prior_sigma")
        end
        length(prior_sigma) == n_state || error("prior_sigma length must match x length")
        S_a_inv = spdiagm(0 => @. 1.0 / (prior_sigma^2))
    else
        S_a_inv = spzeros(Float64, n_state, n_state)
    end

    if isnothing(meas_sigma)
        S_e_inv = spdiagm(0 => ones(Float64, length(y_obs)))
    else
        length(meas_sigma) == length(y_obs) || error("meas_sigma length must match y length")
        S_e_inv = spdiagm(0 => @. 1.0 / (meas_sigma^2))
    end

    for iter in 1:max_iter
        # Evaluate Jacobian at current state xᵢ.
        J = isnothing(jacobian_eval) ? ForwardDiff.jacobian(fm, x) : jacobian_eval(x)
        # Rodgers Eq. 5.9 in explicit x_a form:
        # x_{i+1} = x_a + (S_a^{-1} + K_i' S_e^{-1} K_i)^{-1} K_i' S_e^{-1}
        #           [y - f(x_i) + K_i(x_i - x_a)]
        A = S_a_inv + J' * S_e_inv * J
        innovation = y_obs .- y .+ J * (x .- x_a)
        rhs = J' * S_e_inv * innovation
        x_trial = x_a .+ (A \ rhs)

        step_norm = norm(x_trial .- x)
        x_norm = max(norm(x), eps())

        y_trial = try
            fm(x_trial)                # f(xᵢ₊₁)
        catch err
            message = "model evaluation failed after Rodgers step: $(typeof(err))"
            break
        end
        r_trial = y_obs .- y_trial
        obj_trial = 0.5 * dot(r_trial, r_trial)

        push!(history, (
            iter = iter,
            objective = obj,
            rmse = sqrt(mean(r .^ 2)),
            step_norm = step_norm,
            objective_trial = obj_trial,
        ))

        if verbose
            println(
                "iter=", iter,
                " obj=", obj,
                " obj_trial=", obj_trial,
                " rmse=", sqrt(mean(r .^ 2)),
                " |Δx|=", step_norm,
            )
        end

        obj_prev = obj
        x = x_trial
        y = y_trial
        r = r_trial
        obj = obj_trial

        rel_obj = abs(obj_prev - obj) / max(obj_prev, eps())
        rel_step = step_norm / x_norm

        if rel_obj < rel_obj_tol || rel_step < rel_step_tol
            converged = true
            message = "converged"
            break
        end
    end

    return (
        x = x,
        y = y,
        residual = r,
        objective = obj,
        rmse = sqrt(mean(r .^ 2)),
        converged = converged,
        message = message,
        history = history,
    )
end

"""
    rodgers_eq59_one_step(fm, x_prior, y_obs; kwargs...)

Run exactly one Rodgers Eq. 5.9 update from `x_prior` and return
`(x_next, y_prior, y_next, dx, J, objective_prior, objective_next, cond_A)`.
"""
function rodgers_eq59_one_step(
    fm,
    x_prior::AbstractVector{<:Real},
    y_obs::AbstractVector{<:Real};
    ctx=nothing,
    layout=nothing,
    x_a::AbstractVector{<:Real}=collect(Float64.(x_prior)),
    use_prior::Bool=false,
    prior_sigma::Union{Nothing, AbstractVector{<:Real}}=nothing,
    meas_sigma::Union{Nothing, AbstractVector{<:Real}}=nothing,
    jacobian_eval::Union{Nothing, Function}=nothing,
)
    x = collect(Float64.(x_prior))
    x_a = collect(Float64.(x_a))
    y_prior = fm(x)

    n_state = length(x)
    if use_prior
        isnothing(prior_sigma) && error("use_prior=true requires prior_sigma")
        length(prior_sigma) == n_state || error("prior_sigma length must match x length")
        S_a_inv = spdiagm(0 => @. 1.0 / (prior_sigma^2))
    else
        S_a_inv = spzeros(Float64, n_state, n_state)
    end

    if isnothing(meas_sigma)
        S_e_inv = spdiagm(0 => ones(Float64, length(y_obs)))
    else
        length(meas_sigma) == length(y_obs) || error("meas_sigma length must match y length")
        S_e_inv = spdiagm(0 => @. 1.0 / (meas_sigma^2))
    end

    J = isnothing(jacobian_eval) ? ForwardDiff.jacobian(fm, x) : jacobian_eval(x)
    A = S_a_inv + J' * S_e_inv * J
    cond_A = cond(Matrix(A))
    innovation = y_obs .- y_prior .+ J * (x .- x_a)
    rhs = J' * S_e_inv * innovation
    x_next = x_a .+ (A \ rhs)

    dx = x_next .- x
    y_next = fm(x_next)
    r_prior = y_obs .- y_prior
    r_next = y_obs .- y_next
    obj_prior = 0.5 * dot(r_prior, r_prior)
    obj_next = 0.5 * dot(r_next, r_next)

    return (
        x_next = x_next,
        y_prior = y_prior,
        y_next = y_next,
        dx = dx,
        J = J,
        objective_prior = obj_prior,
        objective_next = obj_next,
        cond_A = cond_A,
    )
end

function _spdiag_invvar(sigma::AbstractVector{<:Real})
    σ = collect(Float64.(sigma))
    # Guard against only numerically pathological scales; keep physically
    # meaningful large sigmas (e.g., 1e23 for weak VCD priors) unchanged.
    @. σ = clamp(abs(σ), 1e-12, 1e100)
    return spdiagm(0 => @. 1.0 / (σ^2))
end

function _cost_with_prior(y_obs, y_mod, x, x_a, S_e_inv, S_a_inv)
    r = y_obs .- y_mod
    dx = x .- x_a
    j_obs = 0.5 * dot(r, S_e_inv * r)
    j_pri = 0.5 * dot(dx, S_a_inv * dx)
    return j_obs + j_pri
end

function _apply_box_constraints!(x, lower, upper)
    @inbounds for i in eachindex(x)
        x[i] = clamp(x[i], lower[i], upper[i])
    end
    return x
end

"""
    _stalled_convergence(dx_rel_hist, redchi2_hist; kwargs...)

Heuristic convergence check for cases where LM cannot find a new accepted step,
but the solution is already near-stationary.

Returns `(is_converged, message)`.
"""
function _stalled_convergence(
    dx_rel_hist::AbstractVector{<:Real},
    redchi2_hist::AbstractVector{<:Real};
    enabled::Bool=true,
    window::Int=3,
    redchi2_target::Float64=5.0,
    redchi2_abs_tol::Float64=0.1,
    redchi2_rel_tol::Float64=0.03,
    dx_rel_tol::Float64=5e-3,
)
    enabled || return (false, "stall criterion disabled")
    n = length(redchi2_hist)
    n >= 2 || return (false, "insufficient iteration history")

    k = min(window, n - 1)
    r_old = Float64(redchi2_hist[n - k])
    r_new = Float64(redchi2_hist[n])
    abs_change = abs(r_new - r_old)
    rel_change = abs_change / max(abs(r_old), eps(Float64))
    dx_last = isempty(dx_rel_hist) ? Inf : Float64(dx_rel_hist[end])

    is_good_chi2 = r_new <= redchi2_target
    is_stable_chi2 = (abs_change <= redchi2_abs_tol) || (rel_change <= redchi2_rel_tol)
    is_small_step = dx_last <= dx_rel_tol

    if is_good_chi2 && (is_stable_chi2 || is_small_step)
        msg = "stalled accepted-step search with near-stationary fit " *
              "(red_chi2=$(r_new), Δred_chi2_abs=$(abs_change), Δred_chi2_rel=$(rel_change), dx_rel_last=$(dx_last))"
        return (true, msg)
    end
    return (false, "stall criterion not met")
end

"""
    make_jacobian_evaluator(fm, x_template; use_preallocated=false)

Return a callable `jac_eval(x)` used to evaluate Jacobians of `fm`.
When `use_preallocated=true`, this reuses both a Jacobian output matrix and
`ForwardDiff.JacobianConfig` for fixed-size `Float64` state vectors.
"""
function make_jacobian_evaluator(
    fm,
    x_template::AbstractVector{<:Real};
    use_preallocated::Bool=false,
)
    if !use_preallocated
        return x -> ForwardDiff.jacobian(fm, x)
    end

    x_ref = collect(Float64.(x_template))
    y_ref = fm(x_ref)
    J = zeros(Float64, length(y_ref), length(x_ref))
    cfg = ForwardDiff.JacobianConfig(fm, x_ref)

    function jac_eval(x::AbstractVector{<:Real})
        if length(x) != size(J, 2) || eltype(x) != Float64
            # Conservative fallback for unexpected runtime type/shape changes.
            return ForwardDiff.jacobian(fm, x)
        end
        ForwardDiff.jacobian!(J, fm, x, cfg)
        return J
    end
    return jac_eval
end

"""
    make_hybrid_jacobian_evaluator(fm, ctx, solar_hres, layout; n_legendre)

Hybrid Jacobian:
- Analytic columns: VCD intercept/slope, SIF-path VCDs, continuum scale, SIF coeffs, Legendre coeffs
- Interpolator-derivative columns: p_o2_hpa, t_o2_k, p_h2o_hpa, t_h2o_k

Pressure/temperature derivatives are obtained directly from dualized LUT calls at
each spectral sample, avoiding expensive whole-model AD for these columns.
"""
function make_hybrid_jacobian_evaluator(
    fm,
    ctx,
    solar_hres::AbstractVector{<:Real},
    layout;
    n_legendre::Int,
)
    n_state = layout.n_state
    spectral_axis = collect(Float64.(ctx.spectral_axis))
    z_hres = _normalized_grid(ctx.λ_hres)
    z_lres = _normalized_grid(ctx.λ)
    K = Matrix{Float64}(hasproperty(ctx, :kernel_rsr_out) ? ctx.kernel_rsr_out : ctx.kernel.RSR_out)
    sif_basis = Matrix{Float64}(ctx.sif_basis_hres)
    leg_basis = Matrix{Float64}(_legendre_design_matrix(z_lres, n_legendre))
    solar = collect(Float64.(solar_hres))

    n_hres = length(z_hres)
    n_lres = size(K, 1)
    xs_o2 = zeros(Float64, n_hres)
    xs_h2o = zeros(Float64, n_hres)
    vcd_o2 = zeros(Float64, n_hres)
    vcd_h2o = zeros(Float64, n_hres)
    trans = zeros(Float64, n_hres)
    trans_sif = zeros(Float64, n_hres)
    sif_hres = zeros(Float64, n_hres)
    y_hres = zeros(Float64, n_hres)
    y_lres = zeros(Float64, n_lres)
    poly = zeros(Float64, n_lres)
    dxsdp_o2 = zeros(Float64, n_hres)
    dxsdt_o2 = zeros(Float64, n_hres)
    dxsdp_h2o = zeros(Float64, n_hres)
    dxsdt_h2o = zeros(Float64, n_hres)
    d_hres = zeros(Float64, n_hres)
    d_lres = zeros(Float64, n_lres)
    J = zeros(Float64, n_lres, n_state)

    function apply_hres_column!(col::Int)
        mul!(d_lres, K, d_hres)
        @inbounds @simd for i in 1:n_lres
            J[i, col] = d_lres[i] * poly[i]
        end
        return nothing
    end

    function jac_eval(x::AbstractVector{<:Real})
        if length(x) != n_state || eltype(x) != Float64
            return ForwardDiff.jacobian(fm, x)
        end

        vcd_o2_i = x[layout.idx_vcd_o2_intercept]
        vcd_o2_s = x[layout.idx_vcd_o2_slope]
        vcd_h2o_i = x[layout.idx_vcd_h2o_intercept]
        vcd_h2o_s = x[layout.idx_vcd_h2o_slope]
        p_o2 = x[layout.idx_p_o2_hpa]
        t_o2 = x[layout.idx_t_o2_k]
        p_h2o = x[layout.idx_p_h2o_hpa]
        t_h2o = x[layout.idx_t_h2o_k]
        vcd_o2_sif = x[layout.idx_vcd_o2_sif]
        vcd_h2o_sif = x[layout.idx_vcd_h2o_sif]
        cont = x[layout.idx_continuum_scale]
        sif_coeff = @view x[layout.idx_sif]
        leg_coeff = @view x[layout.idx_legendre]

        o2_pdual = ForwardDiff.Dual{Nothing}(p_o2, ForwardDiff.Partials((1.0, 0.0)))
        o2_tdual = ForwardDiff.Dual{Nothing}(t_o2, ForwardDiff.Partials((0.0, 1.0)))
        h2o_pdual = ForwardDiff.Dual{Nothing}(p_h2o, ForwardDiff.Partials((1.0, 0.0)))
        h2o_tdual = ForwardDiff.Dual{Nothing}(t_h2o, ForwardDiff.Partials((0.0, 1.0)))
        @inbounds for i in eachindex(spectral_axis)
            vo2 = ctx.o2_sitp(spectral_axis[i], o2_pdual, o2_tdual)
            xs_o2[i] = ForwardDiff.value(vo2)
            po2 = ForwardDiff.partials(vo2)
            dxsdp_o2[i] = po2[1]
            dxsdt_o2[i] = po2[2]

            vh2o = ctx.h2o_sitp(spectral_axis[i], h2o_pdual, h2o_tdual)
            xs_h2o[i] = ForwardDiff.value(vh2o)
            ph2o = ForwardDiff.partials(vh2o)
            dxsdp_h2o[i] = ph2o[1]
            dxsdt_h2o[i] = ph2o[2]
        end

        @. vcd_o2 = vcd_o2_i + vcd_o2_s * z_hres
        @. vcd_h2o = vcd_h2o_i + vcd_h2o_s * z_hres
        @. trans = exp(-(vcd_h2o * xs_h2o + vcd_o2 * xs_o2))
        @. trans_sif = exp(-(vcd_h2o_sif * xs_h2o + vcd_o2_sif * xs_o2))

        mul!(sif_hres, sif_basis, sif_coeff)
        @. y_hres = cont * solar * trans + trans_sif * sif_hres
        mul!(y_lres, K, y_hres)
        mul!(poly, leg_basis, leg_coeff)

        # Analytic VCD derivatives.
        @. d_hres = -cont * solar * trans * xs_o2
        apply_hres_column!(layout.idx_vcd_o2_intercept)
        @. d_hres = -cont * solar * trans * (z_hres * xs_o2)
        apply_hres_column!(layout.idx_vcd_o2_slope)
        @. d_hres = -cont * solar * trans * xs_h2o
        apply_hres_column!(layout.idx_vcd_h2o_intercept)
        @. d_hres = -cont * solar * trans * (z_hres * xs_h2o)
        apply_hres_column!(layout.idx_vcd_h2o_slope)

        # Analytic SIF-path VCD derivatives.
        @. d_hres = -trans_sif * sif_hres * xs_o2
        apply_hres_column!(layout.idx_vcd_o2_sif)
        @. d_hres = -trans_sif * sif_hres * xs_h2o
        apply_hres_column!(layout.idx_vcd_h2o_sif)

        # Analytic p/T derivatives using LUT derivatives.
        @. d_hres = -cont * solar * trans * (vcd_o2 * dxsdp_o2) - trans_sif * sif_hres * (vcd_o2_sif * dxsdp_o2)
        apply_hres_column!(layout.idx_p_o2_hpa)
        @. d_hres = -cont * solar * trans * (vcd_o2 * dxsdt_o2) - trans_sif * sif_hres * (vcd_o2_sif * dxsdt_o2)
        apply_hres_column!(layout.idx_t_o2_k)
        @. d_hres = -cont * solar * trans * (vcd_h2o * dxsdp_h2o) - trans_sif * sif_hres * (vcd_h2o_sif * dxsdp_h2o)
        apply_hres_column!(layout.idx_p_h2o_hpa)
        @. d_hres = -cont * solar * trans * (vcd_h2o * dxsdt_h2o) - trans_sif * sif_hres * (vcd_h2o_sif * dxsdt_h2o)
        apply_hres_column!(layout.idx_t_h2o_k)

        # Analytic continuum derivative.
        @. d_hres = solar * trans
        apply_hres_column!(layout.idx_continuum_scale)

        # Analytic SIF coefficients.
        for iev in 1:layout.n_ev
            @. d_hres = trans_sif * sif_basis[:, iev]
            apply_hres_column!(layout.idx_sif[iev])
        end

        # Analytic Legendre coefficients.
        for j in 1:layout.n_leg_coeff
            col = layout.idx_legendre[j]
            @inbounds @simd for i in 1:n_lres
                J[i, col] = y_lres[i] * leg_basis[i, j]
            end
        end

        return J
    end

    return jac_eval
end

"""
One Rodgers-style LM step for the MAP objective:
J(x) = 0.5*(y-f)^T S_e^-1 (y-f) + 0.5*(x-x_a)^T S_a^-1 (x-x_a)

The LM damping is applied as a pre-factor on the prior precision term:
    A(γ) = K' S_e^-1 K + γ S_a^-1
with γ >= 0.

Returns accepted/rejected step, updated damping, and diagnostics.
"""
function lm_one_step(
    fm,
    x_curr::AbstractVector{<:Real},
    y_obs::AbstractVector{<:Real},
    ;
    x_a::AbstractVector{<:Real},
    S_e_inv,
    S_a_inv,
    lambda::Float64,
    lambda_up::Float64=10.0,
    lambda_down::Float64=0.3,
    lambda_min::Float64=1e-8,
    lambda_max::Float64=1e8,
    max_inner::Int=8,
    jacobian_eval::Union{Nothing, Function}=nothing,
    x_scale::Union{Nothing, AbstractVector{<:Real}}=nothing,
    lower_bounds::Union{Nothing, AbstractVector{<:Real}}=nothing,
    upper_bounds::Union{Nothing, AbstractVector{<:Real}}=nothing,
)
    x = collect(Float64.(x_curr))
    y = fm(x)
    r0 = y_obs .- y
    ssr0 = 0.5 * dot(r0, r0)
    rmse0 = sqrt(mean(r0 .^ 2))

    J = isnothing(jacobian_eval) ? ForwardDiff.jacobian(fm, x) : jacobian_eval(x)
    H_obs = J' * S_e_inv * J
    g_obs = J' * S_e_inv * (y_obs .- y)
    g_pri = S_a_inv * (x_a .- x)

    s = isnothing(x_scale) ? ones(Float64, length(x)) : collect(Float64.(x_scale))
    length(s) == length(x) || error("x_scale length must match state length")
    s .= max.(abs.(s), 1e-12)

    # Solve in scaled variables: dx = S * du, S = diag(s).
    S = Diagonal(s)
    Hs_obs = S * H_obs * S
    S_as = S * S_a_inv * S

    λ = clamp(lambda, lambda_min, lambda_max)
    accepted = false
    x_best = x
    y_best = y
    cost_best = ssr0
    dx_best = zeros(Float64, length(x))
    cond_A_best = NaN
    rmse_lin_best = NaN
    rmse_try_best = rmse0
    pred_red_best = NaN
    act_red_best = NaN
    rho_best = NaN
    n_try = 0

    for _ in 1:max_inner
        n_try += 1
        # Rodgers-style LM damping on prior precision only.
        A = Hs_obs + λ * S_as
        cond_A_try = cond(Matrix(A))
        rhs = S * (g_obs .+ λ .* g_pri)
        du = A \ rhs
        dx = S * du
        # Linearized residual prediction at x + dx: r_lin ≈ r0 - K*dx
        r_lin = r0 .- J * dx
        ssr_lin = 0.5 * dot(r_lin, r_lin)
        rmse_lin = sqrt(mean(r_lin .^ 2))
        x_try = x .+ dx
        if !isnothing(lower_bounds) && !isnothing(upper_bounds)
            _apply_box_constraints!(x_try, lower_bounds, upper_bounds)
        end

        y_try = try
            fm(x_try)
        catch
            λ = clamp(λ * lambda_up, lambda_min, lambda_max)
            continue
        end
        r_try = y_obs .- y_try
        ssr_try = 0.5 * dot(r_try, r_try)
        rmse_try = sqrt(mean(r_try .^ 2))
        pred_red = ssr0 - ssr_lin
        act_red = ssr0 - ssr_try
        rho = pred_red > 0 ? act_red / pred_red : -Inf

        # Accept/reject solely by spectral residual improvement.
        if isfinite(rmse_try) && rmse_try < rmse0
            accepted = true
            x_best = x_try
            y_best = y_try
            cost_best = ssr_try
            dx_best = x_best .- x
            cond_A_best = cond_A_try
            rmse_lin_best = rmse_lin
            rmse_try_best = rmse_try
            pred_red_best = pred_red
            act_red_best = act_red
            rho_best = rho
            λ = clamp(λ * lambda_down, lambda_min, lambda_max)
            break
        else
            λ = clamp(λ * lambda_up, lambda_min, lambda_max)
        end
    end

    return (
        x_next = x_best,
        y_prior = y,
        y_next = y_best,
        dx = dx_best,
        cost_prior = ssr0,
        cost_next = cost_best,
        accepted = accepted,
        lambda_next = λ,
        inner_tries = n_try,
        cond_A = cond_A_best,
        rmse_prior = rmse0,
        rmse_linear = rmse_lin_best,
        rmse_next = rmse_try_best,
        pred_reduction = pred_red_best,
        act_reduction = act_red_best,
        rho = rho_best,
        J = J,
    )
end

function main()
    config_path = get(
        ENV,
        "PACE_MWE_CONFIG",
        joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE.toml"),
    )
    cfg = TOML.parsefile(config_path)

    fit_cfg = get(cfg, "fit", Dict{String, Any}())
    data_cfg = get(cfg, "data", Dict{String, Any}())
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    state_float_type = parse_float_type(cfg)
    ctx = prepare_mwe_inputs(config_path)

    n_legendre = Int(get(fit_cfg, "n_legendre", 2))
    preallocate_forward = Bool(get(fit_cfg, "preallocate_forward", true))
    preallocate_ad_forward = Bool(get(fit_cfg, "preallocate_ad_forward", false))
    preallocate_jacobian = Bool(get(fit_cfg, "preallocate_jacobian", false))
    use_hybrid_jacobian = Bool(get(fit_cfg, "use_hybrid_jacobian", false))
    max_iter = Int(get(fit_cfg, "max_iter", 12))
    rel_obj_tol = Float64(get(fit_cfg, "rel_obj_tol", 1e-8))
    rel_step_tol = Float64(get(fit_cfg, "rel_step_tol", 1e-8))
    conv_dx_rel_tol = Float64(get(fit_cfg, "conv_dx_rel_tol", 1e-6))
    conv_rmse_rel_tol = Float64(get(fit_cfg, "conv_rmse_rel_tol", 1e-6))
    conv_rmse_abs_tol = Float64(get(fit_cfg, "conv_rmse_abs_tol", 1e-6))
    conv_stall_enable = Bool(get(fit_cfg, "conv_stall_enable", true))
    conv_stall_window = Int(get(fit_cfg, "conv_stall_window", 3))
    conv_stall_redchi2_target = Float64(get(fit_cfg, "conv_stall_redchi2_target", 5.0))
    conv_stall_redchi2_abs_tol = Float64(get(fit_cfg, "conv_stall_redchi2_abs_tol", 0.1))
    conv_stall_redchi2_rel_tol = Float64(get(fit_cfg, "conv_stall_redchi2_rel_tol", 0.03))
    conv_stall_dx_rel_tol = Float64(get(fit_cfg, "conv_stall_dx_rel_tol", 5e-3))
    use_legendre01_prior = Bool(get(fit_cfg, "use_legendre01_prior", true))
    legendre01_prior_sigma_fraction = Float64(get(fit_cfg, "legendre01_prior_sigma_fraction", 0.2))
    use_legendre_higher_prior = Bool(get(fit_cfg, "use_legendre_higher_prior", true))
    legendre_higher_sigma = Float64(get(fit_cfg, "legendre_higher_sigma", 1.0))
    use_vcd_slope_prior = Bool(get(fit_cfg, "use_vcd_slope_prior", true))
    vcd_slope_prior_sigma_factor = Float64(get(fit_cfg, "vcd_slope_prior_sigma_factor", 1.0))
    sif_sigma = Float64(get(fit_cfg, "sif_sigma", 1e12))
    prior_min_sigma = Float64(get(fit_cfg, "prior_min_sigma", 1e-3))
    prior_sigma_default = Float64(get(fit_cfg, "prior_sigma_default", 1e30))
    n_plot_steps = Int(get(fit_cfg, "n_plot_steps", 3))
    iter_plot_file = String(get(fit_cfg, "iter_plot_file", "toy_fit_iter_spectral.png"))
    iter_log_file = String(get(fit_cfg, "iter_log_file", "toy_fit_iter_log.csv"))
    # LM controls
    lm_lambda0 = Float64(get(fit_cfg, "lm_lambda0", 1.0))
    lm_lambda_up = Float64(get(fit_cfg, "lm_lambda_up", 2.0))
    lm_lambda_down = Float64(get(fit_cfg, "lm_lambda_down", 0.7))
    lm_lambda_min = Float64(get(fit_cfg, "lm_lambda_min", 1e-8))
    lm_lambda_max = Float64(get(fit_cfg, "lm_lambda_max", 1e8))
    lm_max_inner = Int(get(fit_cfg, "lm_max_inner", 8))
    meas_sigma = Float64(get(fit_cfg, "meas_sigma", 0.01))
    # Prior settings requested by user
    p_prior_hpa = Float64(get(fit_cfg, "p_prior_hpa", 700.0))
    p_sigma_hpa = Float64(get(fit_cfg, "p_sigma_hpa", 200.0))
    t_prior_k = Float64(get(fit_cfg, "t_prior_k", 280.0))
    t_sigma_k = Float64(get(fit_cfg, "t_sigma_k", 20.0))
    vcd_o2_sigma = Float64(get(fit_cfg, "vcd_o2_sigma", 1e23))
    vcd_h2o_sigma = Float64(get(fit_cfg, "vcd_h2o_sigma", 3e22))
    # Box constraints around p/T priors
    use_pt_constraints = Bool(get(fit_cfg, "use_pt_constraints", true))
    pt_constraint_sigma_mult = Float64(get(fit_cfg, "pt_constraint_sigma_mult", 3.0))
    meas_sigma > 0.0 || error("fit.meas_sigma must be > 0")

    solar_file = get(data_cfg, "solar_file", "solar_merged_20200720_600_33300_100.out")
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)
    solar_hres, _ = load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
    solar_hres = state_float_type.(solar_hres)

    pace_file = get(pace_cfg, "pace_file", "sample_granule_20240830T131442_new_chl.nc")
    pace_path = isabspath(pace_file) ? pace_file : joinpath(ctx.paths.base_dir, pace_file)
    wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
    spectrum_var = String(get(pace_cfg, "spectrum_var", "radiance_red"))
    pixel_idx = Int(get(pace_cfg, "pixel_index", 600))
    scan_idx = Int(get(pace_cfg, "scan_index", 800))
    y_obs, obs_info = load_pace_spectrum_on_grid(
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
        make_jacobian_evaluator(
            fm,
            x0;
            use_preallocated=preallocate_jacobian,
        )
    end

    # Good first guess for scale from linear least-squares.
    y0 = fm(x0)
    x0[layout.idx_continuum_scale] = dot(y_obs, y0) / dot(y0, y0)

    # Prior setup (Rodgers Eq. 5.9):
    # Set priors for first two Legendre terms (P0,P1) from an envelope-like fit.
    x_a = copy(x0)
    prior_sigma = fill(prior_sigma_default, length(x0))
    use_prior = false
    if use_legendre01_prior && length(layout.idx_legendre) >= 1
        leg0_idx = first(layout.idx_legendre)
        y_base = fm(x0)
        ratio = y_obs ./ max.(abs.(y_base), eps(Float64))
        z = _normalized_grid(ctx.λ)
        A01 = hcat(ones(length(z)), z)

        # Weight by radiance level to emphasize upper-envelope shape.
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
        use_prior = true
    end

    # Moderate priors for higher Legendre terms (P2 and above) to reduce
    # degeneracy with continuum/VCD-slope terms when using many polynomial DOFs.
    if use_legendre_higher_prior && length(layout.idx_legendre) >= 3
        for j in 3:length(layout.idx_legendre)
            idx = layout.idx_legendre[j]
            x_a[idx] = 0.0
            prior_sigma[idx] = max(legendre_higher_sigma, prior_min_sigma)
        end
        use_prior = true
    end

    # Requested priors for p/T and VCD.
    x_a[layout.idx_p_o2_hpa] = p_prior_hpa
    x_a[layout.idx_p_h2o_hpa] = p_prior_hpa
    x_a[layout.idx_t_o2_k] = t_prior_k
    x_a[layout.idx_t_h2o_k] = t_prior_k
    prior_sigma[layout.idx_p_o2_hpa] = p_sigma_hpa
    prior_sigma[layout.idx_p_h2o_hpa] = p_sigma_hpa
    prior_sigma[layout.idx_t_o2_k] = t_sigma_k
    prior_sigma[layout.idx_t_h2o_k] = t_sigma_k

    x_a[layout.idx_vcd_o2_intercept] = x0[layout.idx_vcd_o2_intercept]
    x_a[layout.idx_vcd_h2o_intercept] = x0[layout.idx_vcd_h2o_intercept]
    prior_sigma[layout.idx_vcd_o2_intercept] = vcd_o2_sigma
    prior_sigma[layout.idx_vcd_h2o_intercept] = vcd_h2o_sigma
    # SIF-path VCD priors: shorter path than direct solar beam (initially 50%).
    x_a[layout.idx_vcd_o2_sif] = x0[layout.idx_vcd_o2_sif]
    x_a[layout.idx_vcd_h2o_sif] = x0[layout.idx_vcd_h2o_sif]
    prior_sigma[layout.idx_vcd_o2_sif] = vcd_o2_sigma
    prior_sigma[layout.idx_vcd_h2o_sif] = vcd_h2o_sigma

    # VCD slope priors: mean slope = 0.
    if use_vcd_slope_prior
        x_a[layout.idx_vcd_o2_slope] = 0.0
        x_a[layout.idx_vcd_h2o_slope] = 0.0
        prior_sigma[layout.idx_vcd_o2_slope] = max(vcd_o2_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
        prior_sigma[layout.idx_vcd_h2o_slope] = max(vcd_h2o_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
    end

    # Explicit SIF priors: mean = 0 for all SIF coefficients.
    x_a[layout.idx_sif] .= 0.0
    prior_sigma[layout.idx_sif] .= max(sif_sigma, prior_min_sigma)
    use_prior = true

    S_a_inv = _spdiag_invvar(prior_sigma)
    S_e_inv = spdiagm(0 => fill(1.0 / (meas_sigma^2), length(y_obs)))

    # Parameter scaling for LM updates (conditioning improvement).
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
    x_scale[layout.idx_continuum_scale] = max(abs(x0[layout.idx_continuum_scale]), 10.0)
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

    println("Fitting one PACE spectrum")
    println("  config: ", config_path)
    println("  selected pixel/scan: ", (obs_info.pixel_idx, obs_info.scan_idx))
    println("  n_state: ", layout.n_state, " (nEV=", layout.n_ev, ", n_legendre=", layout.n_legendre, ")")
    println("  state float type: ", state_float_type)
    println("  LUT interpolation mode: ", getproperty(ctx, :lut_interpolation))
    println("  forward preallocation: ", preallocate_forward)
    println("  AD forward preallocation (other numeric types): ", preallocate_ad_forward)
    println("  Jacobian preallocation (ForwardDiff.jacobian!): ", preallocate_jacobian)
    println("  Hybrid Jacobian (analytic + AD for p/T): ", use_hybrid_jacobian)
    println("  initial scale: ", x0[layout.idx_continuum_scale])
    println("  measurement sigma (1σ): ", meas_sigma)
    println(
        "  stalled-convergence check: ",
        conv_stall_enable,
        " (window=", conv_stall_window,
        ", red_chi2_target=", conv_stall_redchi2_target,
        ", red_chi2_abs_tol=", conv_stall_redchi2_abs_tol,
        ", red_chi2_rel_tol=", conv_stall_redchi2_rel_tol,
        ", dx_rel_tol=", conv_stall_dx_rel_tol,
        ")",
    )
    if use_prior
        leg0_idx = first(layout.idx_legendre)
        println("  priors on Legendre coeffs:")
        println("    x_a[leg0] = ", x_a[leg0_idx], "  sigma[leg0] = ", prior_sigma[leg0_idx])
        if length(layout.idx_legendre) >= 2
            leg1_idx = layout.idx_legendre[2]
            println("    x_a[leg1] = ", x_a[leg1_idx], "  sigma[leg1] = ", prior_sigma[leg1_idx])
        end
        println("  priors on VCD:")
        println(
            "    x_a[vcd_o2_intercept] = ", x_a[layout.idx_vcd_o2_intercept],
            "  sigma = ", prior_sigma[layout.idx_vcd_o2_intercept],
        )
        println(
            "    x_a[vcd_h2o_intercept] = ", x_a[layout.idx_vcd_h2o_intercept],
            "  sigma = ", prior_sigma[layout.idx_vcd_h2o_intercept],
        )
        println(
            "    x_a[vcd_o2_sif] = ", x_a[layout.idx_vcd_o2_sif],
            "  sigma = ", prior_sigma[layout.idx_vcd_o2_sif],
        )
        println(
            "    x_a[vcd_h2o_sif] = ", x_a[layout.idx_vcd_h2o_sif],
            "  sigma = ", prior_sigma[layout.idx_vcd_h2o_sif],
        )
        println("  priors on VCD slopes:")
        println(
            "    x_a[vcd_o2_slope] = ", x_a[layout.idx_vcd_o2_slope],
            "  sigma = ", prior_sigma[layout.idx_vcd_o2_slope],
        )
        println(
            "    x_a[vcd_h2o_slope] = ", x_a[layout.idx_vcd_h2o_slope],
            "  sigma = ", prior_sigma[layout.idx_vcd_h2o_slope],
        )
        println("  priors on p/T:")
        println("    p prior = ", p_prior_hpa, " sigma = ", p_sigma_hpa)
        println("    T prior = ", t_prior_k, " sigma = ", t_sigma_k)
        println("  priors on SIF coeffs:")
        println("    x_a[sif_ev*] = 0.0  sigma = ", prior_sigma[first(layout.idx_sif)])
        println("  LM scales:")
        println("    vcd_o2 scale = ", x_scale[layout.idx_vcd_o2_intercept], "  slope scale = ", x_scale[layout.idx_vcd_o2_slope])
        println("    vcd_h2o scale = ", x_scale[layout.idx_vcd_h2o_intercept], "  slope scale = ", x_scale[layout.idx_vcd_h2o_slope])
        println("    vcd_o2_sif scale = ", x_scale[layout.idx_vcd_o2_sif], "  vcd_h2o_sif scale = ", x_scale[layout.idx_vcd_h2o_sif])
        println("    p scale = ", p_sigma_hpa, "  T scale = ", t_sigma_k, "  continuum scale = ", x_scale[layout.idx_continuum_scale])
        if use_pt_constraints
            println(
                "  p/T box constraints: ±", pt_constraint_sigma_mult, "σ ",
                "(p in [", lower_bounds[layout.idx_p_o2_hpa], ", ", upper_bounds[layout.idx_p_o2_hpa], "], ",
                "T in [", lower_bounds[layout.idx_t_o2_k], ", ", upper_bounds[layout.idx_t_o2_k], "])",
            )
        end
    else
        println("  prior: disabled")
    end

    # Multi-step LM from prior state (spectral-space comparison).
    x_curr = copy(x_a)
    y_curr = copy(fm(x_curr))
    dof = max(length(y_obs) - layout.n_state, 1)
    chi2_curr = sum(((y_obs .- y_curr) ./ meas_sigma) .^ 2)
    x_series = [copy(x_curr)]
    y_series = [copy(y_curr)]
    obj_series = [_cost_with_prior(y_obs, y_curr, x_curr, x_a, S_e_inv, S_a_inv)]
    rmse_series = [sqrt(mean((y_obs .- y_curr) .^ 2))]
    chi2_series = [chi2_curr]
    redchi2_series = [chi2_curr / dof]
    dx_norm_series = Float64[]
    dx_rel_series = Float64[]
    cond_series = Float64[]
    rmse_linear_series = Float64[]
    rho_series = Float64[]
    lambda_series = Float64[lm_lambda0]
    accepted_series = Bool[]
    failed_step = 0
    failed_error = ""
    converged = false
    convergence_reason = ""
    λ = lm_lambda0

    n_steps = max(n_plot_steps, 0)
    for istep in 1:n_steps
        x_prev = copy(x_curr)
        rmse_prev = rmse_series[end]
        step = try
            lm_one_step(
                fm,
                x_curr,
                y_obs;
                x_a=x_a,
                S_e_inv=S_e_inv,
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
            )
        catch err
            failed_step = istep
            failed_error = string(typeof(err))
            break
        end
        λ = step.lambda_next
        if !step.accepted
            stalled_conv, stalled_msg = _stalled_convergence(
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
                convergence_reason = stalled_msg
            else
                failed_step = istep
                failed_error = "no accepted LM update after $(lm_max_inner) inner tries"
            end
            break
        end
        x_curr = step.x_next
        push!(x_series, copy(x_curr))
        y_curr = copy(step.y_next)
        push!(y_series, copy(y_curr))
        push!(obj_series, step.cost_next)
        push!(rmse_series, sqrt(mean((y_obs .- y_curr) .^ 2)))
        chi2_curr = sum(((y_obs .- y_curr) ./ meas_sigma) .^ 2)
        push!(chi2_series, chi2_curr)
        push!(redchi2_series, chi2_curr / dof)
        push!(dx_norm_series, norm(step.dx))
        push!(cond_series, step.cond_A)
        push!(rmse_linear_series, step.rmse_linear)
        push!(rho_series, step.rho)
        push!(lambda_series, λ)
        push!(accepted_series, step.accepted)

        dx_rel = norm(step.dx) / max(norm(x_prev), eps(Float64))
        push!(dx_rel_series, dx_rel)
        rmse_curr = rmse_series[end]
        rmse_abs_change = abs(rmse_curr - rmse_prev)
        rmse_rel_change = rmse_abs_change / max(abs(rmse_prev), eps(Float64))
        if dx_rel < conv_dx_rel_tol ||
           rmse_rel_change < conv_rmse_rel_tol ||
           rmse_abs_change < conv_rmse_abs_tol
            converged = true
            convergence_reason = "dx_rel=$(dx_rel), rmse_abs_change=$(rmse_abs_change), rmse_rel_change=$(rmse_rel_change)"
            break
        end
    end

    println()
    println("LM multi-step summary")
    println(
        "  objective prior: ", obj_series[1],
        "   RMSE prior: ", rmse_series[1],
        "   chi2 prior: ", chi2_series[1],
        "   red_chi2 prior: ", redchi2_series[1],
    )
    for istep in 1:length(dx_norm_series)
        println(
            "  step ", istep,
            " objective: ", obj_series[istep + 1],
            "   RMSE: ", rmse_series[istep + 1],
            "   chi2: ", chi2_series[istep + 1],
            "   red_chi2: ", redchi2_series[istep + 1],
            "   |dx|: ", dx_norm_series[istep],
            "   |dx|/|x|: ", dx_rel_series[istep],
            "   cond(A): ", cond_series[istep],
            "   RMSE_lin: ", rmse_linear_series[istep],
            "   rho: ", rho_series[istep],
            "   λ: ", lambda_series[istep + 1],
            "   accepted: ", accepted_series[istep],
        )
    end
    if failed_step > 0
        println("  stopped early at step ", failed_step, " due to model failure: ", failed_error)
    elseif converged
        println("  converged: ", convergence_reason)
    end

    state_names = state_names_simple(ctx; n_legendre=n_legendre)
    println("State vector by step")
    for istep in 0:(length(x_series) - 1)
        dxn = istep == 0 ? 0.0 : dx_norm_series[istep]
        dxr = istep == 0 ? NaN : dx_rel_series[istep]
        cnd = istep == 0 ? NaN : cond_series[istep]
        println(
            "  step ", istep,
            " | obj=", obj_series[istep + 1],
            " rmse=", rmse_series[istep + 1],
            " chi2=", chi2_series[istep + 1],
            " red_chi2=", redchi2_series[istep + 1],
            " |dx|=", dxn,
            " |dx|/|x|=", dxr,
            " cond(A)=", cnd,
            " λ=", lambda_series[istep + 1],
        )
        xk = x_series[istep + 1]
        for j in eachindex(state_names)
            println("    ", state_names[j], " = ", xk[j])
        end
    end

    # Save compact iteration diagnostics and all state elements to CSV.
    log_path = isabspath(iter_log_file) ? iter_log_file : joinpath(@__DIR__, iter_log_file)
    open(log_path, "w") do io
        println(
            io,
            join(
                vcat(
                    [
                        "step",
                        "objective",
                        "rmse",
                        "chi2",
                        "reduced_chi2",
                        "rmse_linear",
                        "rho",
                        "dx_norm",
                        "dx_rel",
                        "cond_A",
                        "lambda",
                        "accepted",
                    ],
                    state_names,
                ),
                ",",
            ),
        )
        for istep in 0:(length(x_series) - 1)
            dxn = istep == 0 ? 0.0 : dx_norm_series[istep]
            dxr = istep == 0 ? NaN : dx_rel_series[istep]
            cnd = istep == 0 ? NaN : cond_series[istep]
            rlin = istep == 0 ? NaN : rmse_linear_series[istep]
            rho = istep == 0 ? NaN : rho_series[istep]
            acc = istep == 0 ? true : accepted_series[istep]
            row = String[
                string(istep),
                string(obj_series[istep + 1]),
                string(rmse_series[istep + 1]),
                string(chi2_series[istep + 1]),
                string(redchi2_series[istep + 1]),
                string(rlin),
                string(rho),
                string(dxn),
                string(dxr),
                string(cnd),
                string(lambda_series[istep + 1]),
                string(acc),
            ]
            append!(row, string.(x_series[istep + 1]))
            println(io, join(row, ","))
        end
    end
    println("  saved iteration log: ", log_path)

    n_steps_done = length(y_series) - 1

    p1 = plot(
        ctx.λ,
        y_obs;
        label="PACE measurement",
        lw=2.8,
        color=:black,
        xlabel="Wavelength [nm]",
        ylabel="Radiance",
        title="Spectral Agreement: Prior + $(n_steps_done) Rodgers Steps",
        size=(1000, 800),
    )

    plot!(p1, ctx.λ, y_series[1]; label="Model @ prior", lw=2.0, color=:steelblue, ls=:dash)
    step_colors = [:firebrick, :darkorange, :forestgreen, :purple, :brown]
    for istep in 1:n_steps_done
        c = step_colors[mod1(istep, length(step_colors))]
        plot!(p1, ctx.λ, y_series[istep + 1]; label="Model @ step $istep", lw=2.0, color=c)
    end

    p2 = plot(
        ctx.λ,
        y_obs .- y_series[1];
        label="Residual @ prior",
        lw=2.0,
        color=:steelblue,
        ls=:dash,
        xlabel="Wavelength [nm]",
        ylabel="Obs - Model",
        title="Spectral Residuals",
    )
    for istep in 1:n_steps_done
        c = step_colors[mod1(istep, length(step_colors))]
        plot!(p2, ctx.λ, y_obs .- y_series[istep + 1]; label="Residual @ step $istep", lw=1.8, color=c)
    end
    hline!(p2, [0.0]; color=:black, ls=:dot, lw=1.0, label="")

    p = plot(p1, p2; layout=(2, 1), link=:x, size=(1000, 800))

    save_path = isabspath(iter_plot_file) ? iter_plot_file : joinpath(@__DIR__, iter_plot_file)
    savefig(p, save_path)
    println("  saved plot: ", save_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
