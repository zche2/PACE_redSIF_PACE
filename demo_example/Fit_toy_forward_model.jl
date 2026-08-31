#!/usr/bin/env julia

using TOML
using LinearAlgebra
using Statistics
using ForwardDiff
using SparseArrays
using Plots
using JLD2

include(joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE_Functions.jl"))
using .SimplePACEXSecFitMWEFunctions

include(joinpath(@__DIR__, "toy_forward_model.jl"))

"""
    make_Se_inv_from_snr(y, band_snr_coeffs)

Compute diagonal measurement error covariance inverse S_e^{-1} from
per-pixel SNR model encoded in `band_snr_coeffs`.

The noise model is:
    σ²(λ) = c₁(λ) + c₂(λ) * y(λ)

`band_snr_coeffs` must be a Dict with keys `"c1"` and `"c2"`, each a
vector of length `n` (one entry per spectral pixel).

Returns a sparse diagonal matrix `S_e_inv` of size `(n, n)`.
"""
function make_Se_inv_from_snr(
    y::AbstractVector{<:Real},
    band_snr_coeffs::Dict,
)
    n = length(y)
    c1 = band_snr_coeffs["c1"]
    c2 = band_snr_coeffs["c2"]
    length(c1) == n || error(
        "band_snr_coeffs c1 length ($(length(c1))) must match spectrum length ($n)"
    )
    length(c2) == n || error(
        "band_snr_coeffs c2 length ($(length(c2))) must match spectrum length ($n)"
    )
    sigma2 = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        sigma2[i] = Float64(c1[i]) + Float64(c2[i]) * y[i]
    end
    return spdiagm(0 => @. 1.0 / sigma2)
end


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

function _get_sigma2(radiance, c1, c2)
    # Simple noise model: σ² = c₁ + (c₂ * radiance)
    return c1 .+ (c2 .* radiance)
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
- Analytic columns: VCD intercept/slope, SIF-path VCDs, SIF coeffs, Legendre coeffs
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
    # z_lres = _normalized_grid(ctx.λ)
    K = Matrix{Float64}(hasproperty(ctx, :kernel_rsr_out) ? ctx.kernel_rsr_out : ctx.kernel.RSR_out)
    sif_basis = Matrix{Float64}(ctx.sif_basis_hres)
    leg_basis = Matrix{Float64}(_legendre_design_matrix(z_hres, n_legendre))
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
    rho = zeros(Float64, n_hres)
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
            J[i, col] = d_lres[i]
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
        mul!(rho, leg_basis, leg_coeff)
        @. y_hres = solar * trans * rho / π + trans_sif * sif_hres

        # Analytic VCD derivatives (solar term includes rho/π).
        @. d_hres = -solar * trans * rho / π * xs_o2
        apply_hres_column!(layout.idx_vcd_o2_intercept)
        @. d_hres = -solar * trans * rho / π * (z_hres * xs_o2)
        apply_hres_column!(layout.idx_vcd_o2_slope)
        @. d_hres = -solar * trans * rho / π * xs_h2o
        apply_hres_column!(layout.idx_vcd_h2o_intercept)
        @. d_hres = -solar * trans * rho / π * (z_hres * xs_h2o)
        apply_hres_column!(layout.idx_vcd_h2o_slope)

        # Analytic SIF-path VCD derivatives (sif_hres is raw, before trans_sif multiply).
        @. d_hres = -trans_sif * sif_hres * xs_o2
        apply_hres_column!(layout.idx_vcd_o2_sif)
        @. d_hres = -trans_sif * sif_hres * xs_h2o
        apply_hres_column!(layout.idx_vcd_h2o_sif)

        # Analytic p/T derivatives using LUT derivatives (solar term includes rho/π).
        @. d_hres = -solar * trans * rho / π * (vcd_o2 * dxsdp_o2) - trans_sif * sif_hres * (vcd_o2_sif * dxsdp_o2)
        apply_hres_column!(layout.idx_p_o2_hpa)
        @. d_hres = -solar * trans * rho / π * (vcd_o2 * dxsdt_o2) - trans_sif * sif_hres * (vcd_o2_sif * dxsdt_o2)
        apply_hres_column!(layout.idx_t_o2_k)
        @. d_hres = -solar * trans * rho / π * (vcd_h2o * dxsdp_h2o) - trans_sif * sif_hres * (vcd_h2o_sif * dxsdp_h2o)
        apply_hres_column!(layout.idx_p_h2o_hpa)
        @. d_hres = -solar * trans * rho / π * (vcd_h2o * dxsdt_h2o) - trans_sif * sif_hres * (vcd_h2o_sif * dxsdt_h2o)
        apply_hres_column!(layout.idx_t_h2o_k)

        # Analytic SIF coefficients.
        for iev in 1:layout.n_ev
            @. d_hres = trans_sif * sif_basis[:, iev]
            apply_hres_column!(layout.idx_sif[iev])
        end

        # Analytic Legendre coefficients: ∂y_hres/∂a_j = solar * trans / π * leg_basis[:, j].
        for j in 1:layout.n_leg_coeff
            @. d_hres = solar * trans / π * (@view leg_basis[:, j])
            apply_hres_column!(layout.idx_legendre[j])
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

`S_e_inv` may be supplied directly, or constructed on-the-fly from SNR
coefficients (`band_snr_coeffs`) or a scalar fallback (`meas_sigma`).

Returns accepted/rejected step, updated damping, and diagnostics.
"""
function lm_one_step(
    fm,
    x_curr::AbstractVector{<:Real},
    y_obs::AbstractVector{<:Real},
    ;
    x_a::AbstractVector{<:Real},
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
    # Measurement noise options (used only when S_e_inv is nothing)
    use_band_snr::Bool=false,
    band_snr_coeffs=nothing,
    meas_sigma::Float64=0.01,
)
    x = collect(Float64.(x_curr))
    y = fm(x)
    r0 = y_obs .- y
    ssr0 = 0.5 * dot(r0, r0)
    rmse0 = sqrt(mean(r0 .^ 2))

    # Build S_e_inv lazily (only when not supplied by caller).
    Se_inv = if use_band_snr && !isnothing(band_snr_coeffs)
        make_Se_inv_from_snr(y, band_snr_coeffs)
    else
        spdiagm(0 => fill(1.0 / meas_sigma^2, length(y_obs)))
    end
    # println("  | Using S_e_inv with diagonal entries in range [", minimum(diag(Se_inv)), ", ", maximum(diag(Se_inv)), "]")

    J = isnothing(jacobian_eval) ? ForwardDiff.jacobian(fm, x) : jacobian_eval(x)
    H_obs = J' * Se_inv * J
    g_obs = J' * Se_inv * (y_obs .- y)
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
        # Rodgers-style LM: prior Hessian baseline (1+λ)S_as; prior gradient unscaled in rhs.
        A = Hs_obs + (1 + λ) * S_as
        cond_A_try = cond(Matrix(A))
        rhs = S * (g_obs .+ g_pri)
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

    chi2_curr = dot(y_obs .- y_best, Se_inv * (y_obs .- y_best))

    return (
        x_next = x_best,
        y_prior = y,
        y_next = y_best,
        dx = dx_best,
        cost_prior = ssr0,
        cost_next = cost_best,
        chi2_next = chi2_curr,
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

"""
One damped Gauss-Newton (MAP) step for the same objective as [`lm_one_step`](@ref):
    A = K' S_e^{-1} K + S_a^{-1}   (scaled: Hs_obs + S_as)
    full step `Δx_full` solves the normal equations; the applied update is `gn_damping * Δx_full`
    (constant learning-rate style damping; use `< 1` to shrink steps for stability).

Unlike LM, there is no inner retry loop; the step is taken whenever `fm(x + gn_damping*Δx_full)` succeeds.
Returns the same NamedTuple shape as `lm_one_step` (`lambda_next` is `NaN`; `inner_tries` is `1`).
"""
function gn_one_step(
    fm,
    x_curr::AbstractVector{<:Real},
    y_obs::AbstractVector{<:Real},
    ;
    x_a::AbstractVector{<:Real},
    S_a_inv,
    jacobian_eval::Union{Nothing, Function}=nothing,
    x_scale::Union{Nothing, AbstractVector{<:Real}}=nothing,
    lower_bounds::Union{Nothing, AbstractVector{<:Real}}=nothing,
    upper_bounds::Union{Nothing, AbstractVector{<:Real}}=nothing,
    use_band_snr::Bool=false,
    band_snr_coeffs=nothing,
    meas_sigma::Float64=0.01,
    gn_damping::Float64=1.0,
)
    gn_damping > 0.0 || error("gn_damping must be positive, got $gn_damping")
    x = collect(Float64.(x_curr))
    y = fm(x)
    r0 = y_obs .- y
    ssr0 = 0.5 * dot(r0, r0)
    rmse0 = sqrt(mean(r0 .^ 2))

    Se_inv = if use_band_snr && !isnothing(band_snr_coeffs)
        make_Se_inv_from_snr(y, band_snr_coeffs)
    else
        spdiagm(0 => fill(1.0 / meas_sigma^2, length(y_obs)))
    end

    J = isnothing(jacobian_eval) ? ForwardDiff.jacobian(fm, x) : jacobian_eval(x)
    H_obs = J' * Se_inv * J
    g_obs = J' * Se_inv * (y_obs .- y)
    g_pri = S_a_inv * (x_a .- x)

    s = isnothing(x_scale) ? ones(Float64, length(x)) : collect(Float64.(x_scale))
    length(s) == length(x) || error("x_scale length must match state length")
    s .= max.(abs.(s), 1e-12)

    S = Diagonal(s)
    Hs_obs = S * H_obs * S
    S_as = S * S_a_inv * S

    # Full GN: A = Hs_obs + S_as (equivalent to LM with γ=1 on prior precision only)
    A = Hs_obs + S_as
    cond_A_try = cond(Matrix(A))
    rhs = S * (g_obs .+ g_pri)
    du = A \ rhs
    dx_full = S * du
    dx = gn_damping .* dx_full
    r_lin = r0 .- J * dx
    ssr_lin = 0.5 * dot(r_lin, r_lin)
    rmse_lin = sqrt(mean(r_lin .^ 2))
    x_try = x .+ dx
    if !isnothing(lower_bounds) && !isnothing(upper_bounds)
        _apply_box_constraints!(x_try, lower_bounds, upper_bounds)
    end

    y_try = fm(x_try)

    r_try = y_obs .- y_try
    ssr_try = 0.5 * dot(r_try, r_try)
    rmse_try = sqrt(mean(r_try .^ 2))
    pred_red = ssr0 - ssr_lin
    act_red = ssr0 - ssr_try
    rho = pred_red > 0 ? act_red / pred_red : -Inf

    chi2_curr = dot(y_obs .- y_try, Se_inv * (y_obs .- y_try))
    dx_best = x_try .- x

    return (
        x_next = x_try,
        y_prior = y,
        y_next = y_try,
        dx = dx_best,
        cost_prior = ssr0,
        cost_next = ssr_try,
        chi2_next = chi2_curr,
        accepted = true,
        lambda_next = NaN,
        inner_tries = 1,
        cond_A = cond_A_try,
        rmse_prior = rmse0,
        rmse_linear = rmse_lin,
        rmse_next = rmse_try,
        pred_reduction = pred_red,
        act_reduction = act_red,
        rho = rho,
        J = J,
    )
end

mutable struct LbfgsbState
    m::Int
    s_hist::Vector{Vector{Float64}}
    y_hist::Vector{Vector{Float64}}
    rho::Vector{Float64}
    # Step cache: quantities computed at x_best of an accepted step, ready to be
    # reused as the starting values of the very next lbfgsb_one_step call (which
    # begins at that same x_best).  Avoids one fm eval + one Jacobian eval per
    # iteration after the first.
    cached_y    ::Union{Nothing, Vector{Float64}}
    cached_Se_inv                                   # sparse or diagonal matrix
    cached_J    ::Union{Nothing, Matrix{Float64}}
    cached_grad ::Union{Nothing, Vector{Float64}}
end

function LbfgsbState(m::Int)
    m > 0 || error("lbfgs_m must be > 0, got $m")
    return LbfgsbState(
        m,
        Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), Float64[],
        nothing, nothing, nothing, nothing,
    )
end

function _lbfgsb_clear_cache!(state::LbfgsbState)
    state.cached_y     = nothing
    state.cached_Se_inv = nothing
    state.cached_J     = nothing
    state.cached_grad  = nothing
    return nothing
end

function _lbfgs_direction(
    grad::Vector{Float64},
    state::LbfgsbState,
)
    q = copy(grad)
    if isempty(state.s_hist)
        return -q
    end
    α = zeros(Float64, length(state.s_hist))
    for i in length(state.s_hist):-1:1
        α[i] = state.rho[i] * dot(state.s_hist[i], q)
        q .-= α[i] .* state.y_hist[i]
    end
    γ = dot(state.s_hist[end], state.y_hist[end]) / max(dot(state.y_hist[end], state.y_hist[end]), 1e-18)
    γ = max(γ, 1e-10)
    r = γ .* q
    for i in eachindex(state.s_hist)
        β = state.rho[i] * dot(state.y_hist[i], r)
        r .+= state.s_hist[i] .* (α[i] - β)
    end
    return -r
end

function _max_step_to_bounds(
    x::Vector{Float64},
    p::Vector{Float64},
    lower::Vector{Float64},
    upper::Vector{Float64},
)
    α_max = Inf
    @inbounds for i in eachindex(x)
        pi = p[i]
        if pi > 1e-30
            α_max = min(α_max, (upper[i] - x[i]) / pi)
        elseif pi < -1e-30
            α_max = min(α_max, (lower[i] - x[i]) / pi)
        end
    end
    return isfinite(α_max) ? max(α_max, 0.0) : 0.0
end

function _lbfgsb_push_memory!(
    state::LbfgsbState,
    s::Vector{Float64},
    y::Vector{Float64},
)
    ys = dot(y, s)
    ys > 1e-12 || return false
    while length(state.s_hist) >= state.m
        popfirst!(state.s_hist)
        popfirst!(state.y_hist)
        popfirst!(state.rho)
    end
    push!(state.s_hist, s)
    push!(state.y_hist, y)
    push!(state.rho, 1.0 / ys)
    return true
end

"""
One bound-constrained L-BFGS-B-style MAP step with Armijo backtracking.

Memory update policy on rejected steps is explicit:
if no acceptable line-search step is found, `(s, y, rho)` history is unchanged.
"""
function lbfgsb_one_step(
        fm,
        x_curr::AbstractVector{<:Real},
        y_obs::AbstractVector{<:Real},
        ;
        x_a::AbstractVector{<:Real},
        S_a_inv,
        jacobian_eval::Union{Nothing, Function}=nothing,
        lower_bounds::Union{Nothing, AbstractVector{<:Real}}=nothing,
        upper_bounds::Union{Nothing, AbstractVector{<:Real}}=nothing,
        use_band_snr::Bool=false,
        band_snr_coeffs=nothing,
        meas_sigma::Float64=0.01,
        lbfgs_state::LbfgsbState,
        lbfgs_c1::Float64=1e-4,
        lbfgs_max_backtrack::Int=20,
        x_scale::Union{Nothing, AbstractVector{<:Real}}=nothing,
    )
    lbfgs_c1 > 0 || error("lbfgs_c1 must be > 0")
    lbfgs_max_backtrack > 0 || error("lbfgs_max_backtrack must be > 0")

    x = collect(Float64.(x_curr))

    # Reuse fm/Se_inv/J/grad computed at the end of the previous accepted step
    # (which ended at this same x) to avoid redundant evaluations.
    local y, Se_inv, J, grad0
    if lbfgs_state.cached_y !== nothing
        y      = lbfgs_state.cached_y
        Se_inv = lbfgs_state.cached_Se_inv
        J      = lbfgs_state.cached_J
        grad0  = lbfgs_state.cached_grad
        _lbfgsb_clear_cache!(lbfgs_state)
    else
        y = fm(x)
        Se_inv = if use_band_snr && !isnothing(band_snr_coeffs)
            make_Se_inv_from_snr(y, band_snr_coeffs)
        else
            spdiagm(0 => fill(1.0 / meas_sigma^2, length(y_obs)))
        end
        J     = isnothing(jacobian_eval) ? ForwardDiff.jacobian(fm, x) : jacobian_eval(x)
        g_obs = J' * Se_inv * (y_obs .- y)
        g_pri = S_a_inv * (x_a .- x)
        grad0 = -(g_obs .+ g_pri)
    end

    r0    = y_obs .- y
    ssr0  = 0.5 * dot(r0, r0)
    rmse0 = sqrt(mean(r0 .^ 2))
    cost0 = _cost_with_prior(y_obs, y, x, x_a, Se_inv, S_a_inv)

    s = isnothing(x_scale) ? ones(Float64, length(x)) : collect(Float64.(x_scale))
    s .= max.(abs.(s), 1e-12)

    lb = isnothing(lower_bounds) ? fill(-Inf, length(x)) : collect(Float64.(lower_bounds))
    ub = isnothing(upper_bounds) ? fill(Inf, length(x)) : collect(Float64.(upper_bounds))

    # Work in scaled u-space (u = x/s) so all parameters have comparable magnitudes.
    # grad_u = ∂C/∂u = s .* grad_x  (chain rule: ∂C/∂u_i = s_i * ∂C/∂x_i)
    grad0_u = grad0 .* s
    p_u = _lbfgs_direction(grad0_u, lbfgs_state)
    if dot(grad0_u, p_u) > 0
        p_u .= -grad0_u
    end
    p = p_u .* s   # convert u-space direction back to x-space step

    α_cap = _max_step_to_bounds(x, p, lb, ub)
    if α_cap <= 0
        chi2_curr = dot(y_obs .- y, Se_inv * (y_obs .- y))
        return (
            x_next = x,
            y_prior = y,
            y_next = y,
            dx = zeros(Float64, length(x)),
            cost_prior = cost0,
            cost_next = cost0,
            chi2_next = chi2_curr,
            accepted = false,
            lambda_next = NaN,
            inner_tries = 1,
            cond_A = NaN,
            rmse_prior = rmse0,
            rmse_linear = rmse0,
            rmse_next = rmse0,
            pred_reduction = 0.0,
            act_reduction = 0.0,
            rho = NaN,
            J = J,
            reject_reason = "zero feasible step to bounds",
        )
    end

    α = min(1.0, 0.99 * α_cap)
    accepted = false
    x_best = x
    y_best = y
    cost_best = cost0
    rmse_best = rmse0
    n_try = 0

    for _ in 1:lbfgs_max_backtrack
        n_try += 1
        x_try = clamp.(x .+ α .* p, lb, ub)
        y_try = try
            fm(x_try)
        catch
            α *= 0.5
            α < 1e-16 * max(α_cap, 1.0) && break
            continue
        end
        cost_try = _cost_with_prior(y_obs, y_try, x_try, x_a, Se_inv, S_a_inv)
        armijo_rhs = cost0 + lbfgs_c1 * dot(grad0, x_try .- x)
        if cost_try <= armijo_rhs + 1e-12 * max(1.0, abs(cost0))
            accepted = true
            x_best = x_try
            y_best = y_try
            cost_best = cost_try
            rmse_best = sqrt(mean((y_obs .- y_try) .^ 2))
            break
        end
        α *= 0.5
        α < 1e-16 * max(α_cap, 1.0) && break
    end

    if !accepted
        chi2_curr = dot(y_obs .- y, Se_inv * (y_obs .- y))
        return (
            x_next = x,
            y_prior = y,
            y_next = y,
            dx = zeros(Float64, length(x)),
            cost_prior = cost0,
            cost_next = cost0,
            chi2_next = chi2_curr,
            accepted = false,
            lambda_next = NaN,
            inner_tries = max(n_try, 1),
            cond_A = NaN,
            rmse_prior = rmse0,
            rmse_linear = rmse0,
            rmse_next = rmse0,
            pred_reduction = 0.0,
            act_reduction = 0.0,
            rho = NaN,
            J = J,
            reject_reason = "no acceptable Armijo step in $(max(n_try, 1)) tries",
        )
    end

    dx = x_best .- x
    r_lin = r0 .- J * dx
    ssr_lin = 0.5 * dot(r_lin, r_lin)
    rmse_lin = sqrt(mean(r_lin .^ 2))
    ssr_try = 0.5 * dot(y_obs .- y_best, y_obs .- y_best)
    pred_red = ssr0 - ssr_lin
    act_red = ssr0 - ssr_try
    rho = pred_red > 0 ? act_red / pred_red : NaN

    J_new = isnothing(jacobian_eval) ? ForwardDiff.jacobian(fm, x_best) : jacobian_eval(x_best)
    # Build Se_inv from y_best so that grad_new matches what the next call will
    # compute when it reuses this cached gradient (next call's Se_inv is also
    # built from y_best = its starting y).
    Se_inv_new = if use_band_snr && !isnothing(band_snr_coeffs)
        make_Se_inv_from_snr(y_best, band_snr_coeffs)
    else
        spdiagm(0 => fill(1.0 / meas_sigma^2, length(y_obs)))
    end
    g_obs_new = J_new' * Se_inv_new * (y_obs .- y_best)
    g_pri_new = S_a_inv * (x_a .- x_best)
    grad_new = -(g_obs_new .+ g_pri_new)
    # Store curvature pairs in u-space so _lbfgs_direction sees consistent scaling.
    du   = dx ./ s
    dy_u = (grad_new .- grad0) .* s
    _lbfgsb_push_memory!(lbfgs_state, du, dy_u)

    # Cache y_best, Se_inv_new, J_new, grad_new for the next call's starting point.
    # copy(J_new) is needed when jacobian_eval writes into a preallocated buffer.
    lbfgs_state.cached_y      = copy(y_best)
    lbfgs_state.cached_Se_inv = Se_inv_new
    lbfgs_state.cached_J      = copy(J_new)
    lbfgs_state.cached_grad   = copy(grad_new)

    chi2_curr = dot(y_obs .- y_best, Se_inv_new * (y_obs .- y_best))
    return (
        x_next = x_best,
        y_prior = y,
        y_next = y_best,
        dx = dx,
        cost_prior = cost0,
        cost_next = cost_best,
        chi2_next = chi2_curr,
        accepted = true,
        lambda_next = NaN,
        inner_tries = n_try,
        cond_A = NaN,
        rmse_prior = rmse0,
        rmse_linear = rmse_lin,
        rmse_next = rmse_best,
        pred_reduction = pred_red,
        act_reduction = act_red,
        rho = rho,
        J = J,
        reject_reason = "",
    )
end

function main(;
    silent::Bool=false,
    return_rmse_series_only::Bool=false,
    return_benchmark::Bool=false,
    fit_method_override::Union{Nothing,Symbol}=nothing,
    x_init::Union{Nothing,AbstractVector{<:Real}}=nothing,
)
    config_path = get(
        ENV,
        "PACE_MWE_CONFIG",
        joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE_zcheVer.toml"),
    )
    cfg = TOML.parsefile(config_path)

    fit_cfg = get(cfg, "fit", Dict{String, Any}())
    data_cfg = get(cfg, "data", Dict{String, Any}())
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    state_float_type = parse_float_type(cfg)
    ctx = prepare_mwe_inputs(config_path)

    model_variant = Symbol(get(fit_cfg, "model_variant", "standard"))
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
    gn_damping = Float64(get(fit_cfg, "gn_damping", 1.0))
    gn_damping > 0.0 || error("fit.gn_damping must be > 0")
    lbfgs_m = Int(get(fit_cfg, "lbfgs_m", 5))
    lbfgs_c1 = Float64(get(fit_cfg, "lbfgs_c1", 1e-4))
    lbfgs_max_backtrack = Int(get(fit_cfg, "lbfgs_max_backtrack", 20))
    lbfgs_m > 0 || error("fit.lbfgs_m must be > 0")
    lbfgs_c1 > 0.0 || error("fit.lbfgs_c1 must be > 0")
    lbfgs_max_backtrack > 0 || error("fit.lbfgs_max_backtrack must be > 0")
    fit_method = Symbol(get(fit_cfg, "fit_method", "lm"))
    fit_method in (:lm, :gn, :lbfgsb) || error("fit.fit_method must be \"lm\", \"gn\", or \"lbfgsb\", got $(repr(string(fit_method)))")
    if fit_method_override !== nothing
        fit_method = fit_method_override
        fit_method in (:lm, :gn, :lbfgsb) || error("fit_method_override must be :lm, :gn, or :lbfgsb")
    end
    # Measurement error settings requested by user
    use_band_snr = Bool(get(fit_cfg, "use_band_snr", true))
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

    solar_file = get(data_cfg, "solar_file", "solar_merged_20240731_600_33300_100.out")
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
        model_variant=model_variant,
    )
    layout = state_layout_simple(ctx; n_legendre=n_legendre)
    x0 = initial_state_simple(ctx; n_legendre=n_legendre, T=state_float_type, model_variant=model_variant)
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
    x_a[layout.idx_vcd_o2_intercept]  = x0[layout.idx_vcd_o2_intercept]
    x_a[layout.idx_vcd_h2o_intercept] = x0[layout.idx_vcd_h2o_intercept]
    x_a[layout.idx_vcd_o2_sif]  = x0[layout.idx_vcd_o2_sif]
    x_a[layout.idx_vcd_h2o_sif] = x0[layout.idx_vcd_h2o_sif]
    prior_sigma[layout.idx_p_o2_hpa] = p_sigma_hpa
    prior_sigma[layout.idx_p_h2o_hpa] = p_sigma_hpa
    prior_sigma[layout.idx_t_o2_k] = t_sigma_k
    prior_sigma[layout.idx_t_h2o_k] = t_sigma_k

    prior_sigma[layout.idx_vcd_o2_intercept] = vcd_o2_sigma
    prior_sigma[layout.idx_vcd_h2o_intercept] = vcd_h2o_sigma
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

    # S_a_inv: use SIF basis covariance from Spectral_SVD when available; else diagonal prior.
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

    if !silent
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
            if hasproperty(ctx, :sif_prior_cov) && !isnothing(ctx.sif_prior_cov) &&
                length(layout.idx_sif) == size(ctx.sif_prior_cov, 1)
                println("    x_a[sif_ev*] = 0.0  prior = covariance block (inv(cov) from SIF shapes)")
                println("    SIF prior covariance matrix:")
                println("    ", ctx.sif_prior_cov)
            else
                println("    x_a[sif_ev*] = 0.0  sigma = ", prior_sigma[first(layout.idx_sif)])
            end
            println("  LM scales:")
            println("    vcd_o2 scale = ", x_scale[layout.idx_vcd_o2_intercept], "  slope scale = ", x_scale[layout.idx_vcd_o2_slope])
            println("    vcd_h2o scale = ", x_scale[layout.idx_vcd_h2o_intercept], "  slope scale = ", x_scale[layout.idx_vcd_h2o_slope])
            println("    vcd_o2_sif scale = ", x_scale[layout.idx_vcd_o2_sif], "  vcd_h2o_sif scale = ", x_scale[layout.idx_vcd_h2o_sif])
            println("    p scale = ", p_sigma_hpa, "  T scale = ", t_sigma_k)
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

        println("\n" * "="^70)
        method_label = fit_method == :gn ? "Gauss-Newton" : (fit_method == :lbfgsb ? "L-BFGS-B" : "Levenberg-Marquardt")
        println(
            "Starting ",
            method_label,
            " optimization using model variant: ",
            model_variant,
        )
        println("="^70)
    end

    # Multi-step LM or GN from prior state (spectral-space comparison).
    n_forward = Ref(0)
    n_jacobian = Ref(0)
    fm_eval = x -> begin
        n_forward[] += 1
        return fm(x)
    end
    jac_eval = x -> begin
        n_jacobian[] += 1
        return jacobian_eval(x)
    end
    x_curr = isnothing(x_init) ? copy(x_a) : copy(x_init)
    y_curr = copy(fm_eval(x_curr))
    S_e_inv = if !use_band_snr || isnothing(ctx.band_snr_coeffs)
        spdiagm(0 => fill(1.0 / (meas_sigma^2), length(y_obs)))
    else
        make_Se_inv_from_snr(y_curr, ctx.band_snr_coeffs)
    end
    dof = max(length(y_obs) - layout.n_state, 1)
    chi2_curr = dot(y_obs .- y_curr, S_e_inv * (y_obs .- y_curr))
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
    lbfgs_state = LbfgsbState(lbfgs_m)
    t_start_ns = time_ns()

    n_steps = max(n_plot_steps, 0)
    for istep in 1:n_steps
        x_prev = copy(x_curr)
        rmse_prev = rmse_series[end]
        step = try
            if fit_method == :gn
                gn_one_step(
                    fm_eval,
                    x_curr,
                    y_obs;
                    x_a=x_a,
                    S_a_inv=S_a_inv,
                    jacobian_eval=jac_eval,
                    x_scale=x_scale,
                    lower_bounds=lower_bounds,
                    upper_bounds=upper_bounds,
                    use_band_snr=use_band_snr,
                    band_snr_coeffs=ctx.band_snr_coeffs,
                    meas_sigma=meas_sigma,
                    gn_damping=gn_damping,
                )
            elseif fit_method == :lbfgsb
                lbfgsb_one_step(
                    fm_eval,
                    x_curr,
                    y_obs;
                    x_a=x_a,
                    S_a_inv=S_a_inv,
                    jacobian_eval=jac_eval,
                    lower_bounds=lower_bounds,
                    upper_bounds=upper_bounds,
                    use_band_snr=use_band_snr,
                    band_snr_coeffs=ctx.band_snr_coeffs,
                    meas_sigma=meas_sigma,
                    lbfgs_state=lbfgs_state,
                    lbfgs_c1=lbfgs_c1,
                    lbfgs_max_backtrack=lbfgs_max_backtrack,
                    x_scale=x_scale,
                )
            else
                lm_one_step(
                    fm_eval,
                    x_curr,
                    y_obs;
                    x_a=x_a,
                    S_a_inv=S_a_inv,
                    lambda=λ,
                    lambda_up=lm_lambda_up,
                    lambda_down=lm_lambda_down,
                    lambda_min=lm_lambda_min,
                    lambda_max=lm_lambda_max,
                    max_inner=lm_max_inner,
                    jacobian_eval=jac_eval,
                    x_scale=x_scale,
                    lower_bounds=lower_bounds,
                    upper_bounds=upper_bounds,
                    use_band_snr=use_band_snr,
                    band_snr_coeffs=ctx.band_snr_coeffs,
                    meas_sigma=meas_sigma,
                )
            end
        catch err
            failed_step = istep
            failed_error = string(typeof(err))
            break
        end
        if fit_method == :lm
            λ = step.lambda_next
        end
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
                failed_error = if fit_method == :gn
                    "Gauss-Newton step rejected (forward model failed at trial state)"
                elseif fit_method == :lbfgsb
                    step.reject_reason
                else
                    "no accepted LM update after $(lm_max_inner) inner tries"
                end
            end
            break
        end
        x_curr = step.x_next
        push!(x_series, copy(x_curr))
        y_curr = copy(step.y_next)
        push!(y_series, copy(y_curr))
        push!(obj_series, step.cost_next)
        push!(rmse_series, sqrt(mean((y_obs .- y_curr) .^ 2)))
        chi2_curr = step.chi2_next
        push!(chi2_series, chi2_curr)
        push!(redchi2_series, chi2_curr / dof)
        push!(dx_norm_series, norm(step.dx))
        push!(cond_series, step.cond_A)
        push!(rmse_linear_series, step.rmse_linear)
        push!(rho_series, step.rho)
        push!(lambda_series, fit_method == :lm ? λ : NaN)
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

    if return_rmse_series_only
        return rmse_series
    end
    elapsed_wall_s = (time_ns() - t_start_ns) / 1e9

    if return_benchmark
        return (
            method = fit_method,
            rmse_series = rmse_series,
            x_final = copy(x_curr),
            x_prior = copy(x_a),
            state_names = state_names_simple(ctx; n_legendre=n_legendre),
            elapsed_wall_s = elapsed_wall_s,
            n_forward = n_forward[],
            n_jacobian = n_jacobian[],
            converged = converged,
            failed_step = failed_step,
            failed_error = failed_error,
            convergence_reason = convergence_reason,
            n_steps_done = length(x_series) - 1,
            wavelength = collect(Float64.(ctx.λ)),
            y_obs = collect(Float64.(y_obs)),
            y_prior = copy(y_series[1]),
            y_final = copy(y_curr),
        )
    end

    println()
    println((fit_method == :gn ? "GN" : fit_method == :lbfgsb ? "L-BFGS-B" : "LM"), " multi-step summary")
    println(
        "  objective prior: ", obj_series[1],
        "   RMSE prior: ", rmse_series[1],
        "   chi2 prior: ", chi2_series[1],
        "   red_chi2 prior: ", redchi2_series[1],
    )
    for istep in eachindex(dx_norm_series)
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
        println("  stopped early at step ", failed_step, ": ", failed_error)
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

        # sanity check
        if model_variant == :standard
            # print ratio of vcd_h2o to vcd_h2o_sif and vcd_o2 to vcd_o2_sif to check if SIF VCDs are reasonable relative to direct-beam VCDs
            vcd_o2_mean  = xk[layout.idx_vcd_o2_intercept] + xk[layout.idx_vcd_o2_slope] * mean(_normalized_grid(ctx.λ))
            vcd_h2o_mean = xk[layout.idx_vcd_h2o_intercept] + xk[layout.idx_vcd_h2o_slope] * mean(_normalized_grid(ctx.λ))
            vcd_o2_ratio = xk[layout.idx_vcd_o2_sif] / vcd_o2_mean
            vcd_h2o_ratio = xk[layout.idx_vcd_h2o_sif] / vcd_h2o_mean
            println("    vcd_o2_sif / vcd_o2_solar = ", vcd_o2_ratio)
            println("    vcd_h2o_sif / vcd_h2o_solar = ", vcd_h2o_ratio)
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
        lw=3.5,
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
        ylims=(-0.2, 0.2),
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

    # =========================================================================
    # Reconstruct high-resolution SIF from final state
    # =========================================================================
    
    println("\n" * "="^70)
    println("Reconstructing high-resolution SIF from final state")
    println("="^70)
    
    x_final = x_series[end]
    
    # Extract SIF coefficients from final state
    sif_coeff_final = x_final[layout.idx_sif]
    println("Final SIF coefficients: ", sif_coeff_final)
    
    # Reconstruct high-resolution SIF spectrum
    sif_hres = ctx.sif_basis_hres * sif_coeff_final
    
    println("High-resolution SIF reconstructed:")
    println("  Grid: $(length(ctx.λ_hres)) points")
    println("  Range: [$(minimum(ctx.λ_hres)), $(maximum(ctx.λ_hres))] nm")
    println("  SIF range: [$(minimum(sif_hres)), $(maximum(sif_hres))]")
    println("  SIF mean: $(mean(sif_hres))")
    println("  SIF std: $(std(sif_hres))")
    
    # Reconstruct other high-resolution components for context
    vcd_o2_i = x_final[layout.idx_vcd_o2_intercept]
    vcd_o2_s = x_final[layout.idx_vcd_o2_slope]
    vcd_h2o_i = x_final[layout.idx_vcd_h2o_intercept]
    vcd_h2o_s = x_final[layout.idx_vcd_h2o_slope]
    p_o2 = x_final[layout.idx_p_o2_hpa]
    t_o2 = x_final[layout.idx_t_o2_k]
    p_h2o = x_final[layout.idx_p_h2o_hpa]
    t_h2o = x_final[layout.idx_t_h2o_k]
    vcd_o2_sif = x_final[layout.idx_vcd_o2_sif]
    vcd_h2o_sif = x_final[layout.idx_vcd_h2o_sif]
    leg_coeff = x_final[layout.idx_legendre]
    
    # Get cross-sections at high resolution
    spectral_axis = collect(Float64.(ctx.spectral_axis))
    z_hres = _normalized_grid(ctx.λ_hres)
    K      = hasproperty(ctx, :kernel_rsr_out) ? ctx.kernel_rsr_out : ctx.kernel.RSR_out
    
    xs_o2_hres = [ctx.o2_sitp(spectral_axis[i], p_o2, t_o2) for i in eachindex(spectral_axis)]
    xs_h2o_hres = [ctx.h2o_sitp(spectral_axis[i], p_h2o, t_h2o) for i in eachindex(spectral_axis)]
    
    # Calculate transmittances at high and low resolution
    vcd_o2_hres = @. vcd_o2_i + vcd_o2_s * z_hres
    vcd_h2o_hres = @. vcd_h2o_i + vcd_h2o_s * z_hres
    trans_hres = @. exp(-(vcd_h2o_hres * xs_h2o_hres + vcd_o2_hres * xs_o2_hres))
    trans_sif_hres = @. exp(-(vcd_h2o_sif * xs_h2o_hres + vcd_o2_sif * xs_o2_hres))
    
    trans_lres = K * trans_hres
    trans_sif_lres = K * trans_sif_hres
    
    # Calculate reflected solar continuum
    leg_basis = Matrix{Float64}(_legendre_design_matrix(z_hres, n_legendre))
    rho_hres = leg_basis * leg_coeff
    solar_rho_hres = @. solar_hres * rho_hres / π
    solar_rho_lres = K * solar_rho_hres
    solar_continuum_hres = @. solar_hres * trans_hres * rho_hres / π
    solar_continuum_lres = K * solar_continuum_hres
    
    # Calculate SIF contribution to TOA radiance at high and low resolution
    sif_contribution_hres = @. trans_sif_hres * sif_hres
    sif_contribution_lres = K * sif_contribution_hres
    
    # Total high-resolution radiance at high and low resolution
    radiance_total_hres = solar_continuum_hres .+ sif_contribution_hres
    radiance_total_lres_postconv = K * radiance_total_hres
    radiance_total_lres_preconv  = solar_continuum_lres .+ sif_contribution_lres
    
    println("\nHigh-resolution components:")
    println("  Solar continuum contribution: $(mean(solar_continuum_hres))")
    println("  SIF contribution: $(mean(sif_contribution_hres))")
    println("  Total radiance: $(mean(radiance_total_hres))")
    println("  SIF fraction: $(mean(sif_contribution_hres) / mean(radiance_total_hres) * 100)%")
    
    # =========================================================================
    # Plot high-resolution SIF and components
    # =========================================================================
    
    p_sif = plot(
        ctx.λ_hres,
        sif_hres;
        xlabel="Wavelength [nm]",
        ylabel="SIF [mW/m²/sr/nm]",
        title="Retrieved High-Resolution SIF Spectrum",
        lw=2,
        color=:red,
        label="SIF",
        legend=:topright,
        size=(1000, 600)
    )
    
    p_trans = plot(
        ctx.λ_hres,
        trans_hres;
        xlabel="Wavelength [nm]",
        ylabel="Transmittance",
        title="Atmospheric Transmittance",
        lw=2,
        color=:blue,
        label="Direct path (T↓↑)",
        legend=:bottomleft,
        alpha=0.5,
        ylims=(0.0, 1.2)
    )
    plot!(p_trans, ctx.λ, trans_lres;
            lw=2, color=:blue, ls=:dash, label="Direct path (T↓↑) (low-res)")
    plot!(p_trans, ctx.λ_hres, trans_sif_hres;
            lw=2, color=:orange, label="SIF path (T↑)", alpha=0.5)
    plot!(p_trans, ctx.λ, trans_sif_lres;
            lw=2, color=:orange, ls=:dash, label="SIF path (T↑) (low-res)")
    
    p_contrib = plot(
        ctx.λ_hres,
        solar_continuum_hres;
        xlabel="Wavelength [nm]",
        ylabel="Radiance",
        title="High-Res Radiance Components",
        lw=2,
        color=:lightblue,
        label="Solar continuum",
        legend=:topright,
        alpha=0.5
    )
    plot!(p_contrib, ctx.λ_hres, sif_contribution_hres;
            lw=2, color=:red, label="SIF contribution")
    plot!(p_contrib, ctx.λ_hres, radiance_total_hres;
            lw=2, color=:silver, ls=:dash, label="Total", alpha=0.5)
    plot!(p_contrib, ctx.λ, radiance_total_lres_postconv;
            lw=2, color=:black, label="Total (low-res) (post-conv)")
    plot!(p_contrib, ctx.λ, solar_rho_lres;
            lw=2, color=:green, label="Solar radiance baseline (low-res)")

    p_sif_all = plot(p_sif, p_trans, p_contrib;
                        layout=(3, 1),
                        size=(1000, 1200))
    
    sif_plot_file = replace(iter_plot_file, ".png" => "_sif_hires.png")
    sif_save_path = isabspath(sif_plot_file) ? sif_plot_file : joinpath(@__DIR__, sif_plot_file)
    savefig(p_sif_all, sif_save_path)
    println("  saved SIF plot: ", sif_save_path)
    
    # =========================================================================
    # Save high-resolution SIF to file
    # =========================================================================

    sif_output_file = replace(iter_log_file, ".csv" => "_sif_hires.jld2")
    sif_output_path = isabspath(sif_output_file) ? sif_output_file : joinpath(@__DIR__, sif_output_file)
    
    @save sif_output_path λ_hres=ctx.λ_hres sif_hres sif_coeff_final trans_hres trans_sif_hres solar_continuum_hres sif_contribution_hres radiance_total_hres
    
    println("  saved SIF data: ", sif_output_path)
    
    # Print SIF statistics at key wavelengths
    println("\nSIF at key wavelengths:")
    key_wavelengths = [685.0, 740.0, 760.0]
    for λ_key in key_wavelengths
        idx = argmin(abs.(ctx.λ_hres .- λ_key))
        λ_actual = ctx.λ_hres[idx]
        sif_val = sif_hres[idx]
        println("  λ = $(round(λ_actual, digits=2)) nm: SIF = $(round(sif_val, digits=4)) mW/m²/sr/nm")
    end
    
    println("="^70)

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
