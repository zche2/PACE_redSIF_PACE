# Self-contained LM + SVD transmittance helpers for global_fit_pipeline SVD batch retrieval.
# Does NOT include Fit_toy_forward_model.jl (lm_one_step and SNR helpers copied).

using DelimitedFiles
using ForwardDiff
using Interpolations
using LinearAlgebra
using NCDatasets
using SparseArrays
using Statistics

# ── PACE baseline SNR file → per-target λ (pseudo_measurement / band-only pipeline) ─

"""Interpolate Red-band SNR coefficients from `PACE_OCI_L1BLUT_baseline_SNR_*.txt` onto `λ_target` nm."""
function load_pace_band_snr_coeffs(
    pace_snr_path::AbstractString,
    λ_target::AbstractVector{Float64};
    λ_min::Float64,
    λ_max::Float64,
)
    isfile(pace_snr_path) || error("SNR file not found: $pace_snr_path")
    snr_lines = readlines(pace_snr_path)
    header_end_idx = findfirst(line -> occursin("/end_header", line), snr_lines)
    skipstart = isnothing(header_end_idx) ? 0 : Int(header_end_idx)
    snr_data = readdlm(pace_snr_path, String; skipstart = skipstart)
    snr_band_full = parse.(Float64, snr_data[:, 2])
    idx = findall((snr_data[:, 1] .== "Red") .& (λ_min .<= snr_band_full .<= λ_max))
    isempty(idx) && error("No Red SNR rows with λ in ($λ_min, $λ_max) nm in $pace_snr_path")
    wl_snr = snr_band_full[idx]
    c1_raw = parse.(Float64, snr_data[idx, 4])
    c2_raw = parse.(Float64, snr_data[idx, 5])
    p = sortperm(wl_snr)
    wl_s = wl_snr[p]
    c1_s = c1_raw[p]
    c2_s = c2_raw[p]
    itp1 = LinearInterpolation(wl_s, c1_s; extrapolation_bc = Flat())
    itp2 = LinearInterpolation(wl_s, c2_s; extrapolation_bc = Flat())
    c1_out = Float64[itp1(λ) for λ in λ_target]
    c2_out = Float64[itp2(λ) for λ in λ_target]
    return Dict("c1" => c1_out, "c2" => c2_out)
end

# ── SNR / LM utilities (from Fit_toy_forward_model.jl) ────────────────────────

function make_Se_inv_from_snr(
    y::AbstractVector{<:Real},
    band_snr_coeffs::Dict,
)
    n = length(y)
    c1 = band_snr_coeffs["c1"]
    c2 = band_snr_coeffs["c2"]
    length(c1) == n || error("band_snr_coeffs c1 length ($(length(c1))) must match spectrum length ($n)")
    length(c2) == n || error("band_snr_coeffs c2 length ($(length(c2))) must match spectrum length ($n)")
    sigma2 = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        sigma2[i] = Float64(c1[i]) + Float64(c2[i]) * y[i]
    end
    return spdiagm(0 => @. 1.0 / sigma2)
end

"""Build measurement-noise inverse covariance from current model radiance `y` (or fixed σ)."""
function _make_Se_inv(
    y::AbstractVector{<:Real},
    use_band_snr::Bool,
    band_snr_coeffs,
    meas_sigma::Float64,
)
    if use_band_snr && !isnothing(band_snr_coeffs)
        return make_Se_inv_from_snr(y, band_snr_coeffs)
    end
    return spdiagm(0 => fill(1.0 / (meas_sigma^2), length(y)))
end

function _spdiag_invvar(sigma::AbstractVector{<:Real})
    σ = collect(Float64.(sigma))
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

function _stalled_convergence(
    dx_rel_hist::AbstractVector{<:Real},
    redchi2_hist::AbstractVector{<:Real};
    enabled::Bool = true,
    window::Int = 3,
    redchi2_target::Float64 = 5.0,
    redchi2_abs_tol::Float64 = 0.1,
    redchi2_rel_tol::Float64 = 0.03,
    dx_rel_tol::Float64 = 5e-3,
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
        msg = "stalled LM search near-stationary (red_chi2=$(r_new))"
        return (true, msg)
    end
    return (false, "stall criterion not met")
end

function lm_one_step(
    fm,
    x_curr::AbstractVector{<:Real},
    y_obs::AbstractVector{<:Real};
    x_a::AbstractVector{<:Real},
    S_a_inv,
    lambda::Float64,
    lambda_up::Float64 = 10.0,
    lambda_down::Float64 = 0.3,
    lambda_min::Float64 = 1e-8,
    lambda_max::Float64 = 1e8,
    max_inner::Int = 8,
    jacobian_eval = nothing,
    x_scale = nothing,
    lower_bounds = nothing,
    upper_bounds = nothing,
    use_band_snr::Bool = false,
    band_snr_coeffs = nothing,
    meas_sigma::Float64 = 0.01,
)
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
        A = Hs_obs + (1 + λ) * S_as
        cond_A_try = cond(Matrix(A))
        rhs = S * (g_obs .+ g_pri)
        du = A \ rhs
        dx = S * du
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
    lm_one_step_from_J(fm, x_curr, y_obs, y_at_x, J_at_x; ...)

Same as [`lm_one_step`](@ref) but uses precomputed `y_at_x = fm(x_curr)` and `J_at_x` Jacobian at `x_curr`
(so the caller can supply a batched Jacobian). The inner loop still calls `fm(x_try)` for trial states.
"""
function lm_one_step_from_J(
    fm,
    x_curr::AbstractVector{<:Real},
    y_obs::AbstractVector{<:Real},
    y_at_x::AbstractVector{<:Real},
    J_at_x::AbstractMatrix{<:Real};
    x_a::AbstractVector{<:Real},
    S_a_inv,
    lambda::Float64,
    lambda_up::Float64 = 10.0,
    lambda_down::Float64 = 0.3,
    lambda_min::Float64 = 1e-8,
    lambda_max::Float64 = 1e8,
    max_inner::Int = 8,
    x_scale = nothing,
    lower_bounds = nothing,
    upper_bounds = nothing,
    use_band_snr::Bool = false,
    band_snr_coeffs = nothing,
    meas_sigma::Float64 = 0.01,
)
    x = collect(Float64.(x_curr))
    y = collect(Float64.(y_at_x))
    r0 = y_obs .- y
    ssr0 = 0.5 * dot(r0, r0)
    rmse0 = sqrt(mean(r0 .^ 2))
    Se_inv = if use_band_snr && !isnothing(band_snr_coeffs)
        make_Se_inv_from_snr(y, band_snr_coeffs)
    else
        spdiagm(0 => fill(1.0 / meas_sigma^2, length(y_obs)))
    end
    J = collect(Float64.(J_at_x))
    H_obs = J' * Se_inv * J
    g_obs = J' * Se_inv * (y_obs .- y)
    g_pri = S_a_inv * (x_a .- x)
    s = isnothing(x_scale) ? ones(Float64, length(x)) : collect(Float64.(x_scale))
    length(s) == length(x) || error("x_scale length must match state length")
    s .= max.(abs.(s), 1e-12)
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
        A = Hs_obs + (1 + λ) * S_as
        cond_A_try = cond(Matrix(A))
        rhs = S * (g_obs .+ g_pri)
        du = A \ rhs
        dx = S * du
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
    # compute current chi2
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

# ── Spectral design (from toy_forward_model.jl) ───────────────────────────────

function _normalized_grid(λ::AbstractVector{<:Real})
    T = eltype(λ) <: AbstractFloat ? eltype(λ) : Float64
    λf = collect(T.(λ))
    λmin, λmax = extrema(λf)
    λspan = λmax - λmin
    if λspan == 0
        return zeros(T, length(λf))
    end
    λcenter = T(0.5) * (λmin + λmax)
    halfspan = T(0.5) * λspan
    return @. (λf - λcenter) / halfspan
end

function _legendre_design_matrix(z::AbstractVector{<:Real}, n_legendre::Int)
    T = eltype(z) <: AbstractFloat ? eltype(z) : Float64
    zt = collect(T.(z))
    n = length(zt)
    A = zeros(T, n, n_legendre + 1)
    A[:, 1] .= one(T)
    n_legendre == 0 && return A
    A[:, 2] .= zt
    for l in 1:(n_legendre - 1)
        @. A[:, l + 2] = (T(2l + 1) * zt * A[:, l + 1] - T(l) * A[:, l]) / T(l + 1)
    end
    return A
end

# ── SVD state / basis ─────────────────────────────────────────────────────────

function svd_state_layout(; n_pc::Int, n_legendre::Int, n_ev::Int)
    n_leg_coeff = n_legendre + 1
    n_alpha = 1
    n_state = n_pc + n_alpha + n_leg_coeff + n_ev
    i0 = n_pc + n_alpha
    return (
        n_pc = n_pc,
        n_alpha = n_alpha,
        n_legendre = n_legendre,
        n_leg_coeff = n_leg_coeff,
        n_ev = n_ev,
        n_state = n_state,
        idx_pc = 1:n_pc,
        idx_alpha = (n_pc + 1):(n_pc + n_alpha),
        idx_legendre = (i0 + 1):(i0 + n_leg_coeff),
        idx_sif = (i0 + n_leg_coeff + 1):n_state,
    )
end

function load_svd_basis(
    summer_nc::String,
    winter_nc::String,
    λ_pc_target::AbstractVector{<:Float64};
    λ_min::Float64,
    λ_max::Float64,
    n_pc::Int,
    log_transform::Bool = false,
)
    isfile(summer_nc) || error("summer_nc not found: $summer_nc")
    isfile(winter_nc) || error("winter_nc not found: $winter_nc")
    trans_s, bands = NCDataset(summer_nc, "r") do ds
        collect(Float64.(ds["transmittance"][:, :])),
        collect(Float64.(ds["band"][:]))
    end
    trans_w = NCDataset(winter_nc, "r") do ds
        collect(Float64.(ds["transmittance"][:, :]))
    end
    trans = vcat(trans_s, trans_w)
    n_profiles = size(trans, 1)
    ind = findall(λ_min .< bands .< λ_max)
    isempty(ind) && error("No transmittance bands found in [$(λ_min), $(λ_max)] nm")
    bands_sel = bands[ind]
    trans_sel = trans[:, ind]
    a = log_transform ? log.(max.(trans_sel, 1e-10)) : (trans_sel .- 1.0)
    F = svd(a')
    U = F.U
    S = F.S
    S_norm = S ./ sum(S) .* 100
    println("SVD basis: $(n_profiles) profiles, $(length(bands_sel)) bands in [$(λ_min), $(λ_max)] nm")
    nλ = length(λ_pc_target)
    PCs = zeros(Float64, nλ, size(U, 2))
    for k in axes(U, 2)
        itp = LinearInterpolation(bands_sel, U[:, k]; extrapolation_bc = Flat())
        PCs[:, k] .= itp.(λ_pc_target)
    end
    return (; PCs = PCs, S = S, S_norm = S_norm, n_profiles = n_profiles, bands_sel = bands_sel, U = U)
end

function svd_state_names(layout)
    names = String[]
    for k in 1:layout.n_pc
        push!(names, "pc_$k")
    end
    push!(names, "alpha")
    for j in 0:layout.n_legendre
        push!(names, "legendre_p$j")
    end
    for j in 1:layout.n_ev
        push!(names, "sif_ev$j")
    end
    return names
end

function make_svd_forward_model_λ(
    λ_obs::AbstractVector{Float64},
    solar_eff::AbstractVector{Float64},
    PCs_obs::Matrix{Float64},
    sif_basis_obs::Matrix{Float64};
    n_pc::Int,
    n_legendre::Int,
    log_transform::Bool,
)
    length(λ_obs) == length(solar_eff) || error("λ and solar_eff length mismatch")
    PCs = Float64.(PCs_obs[:, 1:n_pc])
    SIF = Float64.(sif_basis_obs)
    z_obs = _normalized_grid(λ_obs)
    leg_basis = _legendre_design_matrix(z_obs, n_legendre)
    layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_legendre, n_ev = size(SIF, 2))

    function fm_svd(x::AbstractVector)
        length(x) == layout.n_state || error("SVD state length $(length(x)) != $(layout.n_state)")
        c_vec = @view x[layout.idx_pc]
        alpha_raw = x[first(layout.idx_alpha)]
        alpha_coeff = 10.0 / (1.0 + exp(-alpha_raw)) + 1.0
        leg_coeff = @view x[layout.idx_legendre]
        sif_coeff = @view x[layout.idx_sif]
        trans_up = log_transform ? exp.(PCs * c_vec) : 1.0 .+ PCs * c_vec
        trans_updown = exp.(alpha_coeff .* log.(max.(trans_up, eps(Float64))))
        rho_obs = leg_basis * leg_coeff
        sif_toa = trans_up .* (SIF * sif_coeff)
        return @.(solar_eff * trans_updown * rho_obs / π + sif_toa)
    end
    return fm_svd, layout
end

function _run_one_svd_retrieval!(
    x_out::Vector{Float64},
    fm,
    jac_eval,
    layout,
    y_obs::Vector{Float64},
    x0::Vector{Float64},
    x_a::Vector{Float64},
    prior_sigma::Vector{Float64},
    lower::Vector{Float64},
    upper::Vector{Float64},
    x_scale::Vector{Float64},
    use_leg01_prior::Bool,
    A01::Matrix{Float64},
    prior_min_sigma::Float64,
    legendre01_frac::Float64,
    use_band_snr::Bool,
    band_snr_coeffs,
    meas_sigma::Float64,
    lm::NamedTuple,
    conv::NamedTuple,
    stall_dx_rel_tol::Float64,
    max_outer_steps::Int,
)
    x_a_loc = copy(x_a)
    σ_loc = copy(prior_sigma)
    if use_leg01_prior && length(layout.idx_legendre) >= 1
        y0 = fm(x0)
        ratio = y_obs ./ max.(abs.(y0), eps(Float64))
        z = A01[:, 2]
        w = y_obs .- minimum(y_obs)
        w .+= max(maximum(w), 1.0) * 1e-6
        sv = sqrt.(w ./ maximum(w))
        c01 = (A01 .* sv) \ (ratio .* sv)
        leg0 = first(layout.idx_legendre)
        x_a_loc[leg0] = c01[1]
        σ_loc[leg0] = max(abs(c01[1]) * legendre01_frac, prior_min_sigma)
        if length(layout.idx_legendre) >= 2
            leg1 = layout.idx_legendre[2]
            x_a_loc[leg1] = c01[2]
            σ_loc[leg1] = max(abs(c01[2]) * legendre01_frac, prior_min_sigma)
        end
    end
    S_a_inv = _spdiag_invvar(σ_loc)
    x_curr = copy(x_a_loc)
    y_curr = fm(x_curr)
    # Rebuild S_e from current model radiance whenever reduced χ² is evaluated
    # (matches lm_one_step weighting and final y_curr ≈ y_obs under good fits).
    S_e_inv = _make_Se_inv(y_curr, use_band_snr, band_snr_coeffs, meas_sigma)
    rmse_prev = sqrt(mean((y_obs .- y_curr) .^ 2))
    chi2_curr = dot(y_obs .- y_curr, S_e_inv * (y_obs .- y_curr))
    # Classical m−n DOF for in-loop reduced-χ² history (stall detection).
    dof_loop = max(length(y_obs) - length(x_curr), 1)
    redchi2_hist = Float64[chi2_curr / dof_loop]
    dx_rel_hist = Float64[]
    λ = lm.lambda0
    n_acc = 0
    converged = false
    status = Int16(0)
    for _ in 1:max_outer_steps
        step = try
            lm_one_step(
                fm,
                x_curr,
                y_obs;
                x_a = x_a_loc,
                S_a_inv = S_a_inv,
                lambda = λ,
                lambda_up = lm.lambda_up,
                lambda_down = lm.lambda_down,
                lambda_min = lm.lambda_min,
                lambda_max = lm.lambda_max,
                max_inner = lm.max_inner,
                jacobian_eval = jac_eval,
                x_scale = x_scale,
                lower_bounds = lower,
                upper_bounds = upper,
                use_band_snr = use_band_snr,
                band_snr_coeffs = band_snr_coeffs,
                meas_sigma = meas_sigma,
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
                enabled = conv.enabled,
                window = conv.window,
                redchi2_target = conv.redchi2_target,
                redchi2_abs_tol = conv.redchi2_abs_tol,
                redchi2_rel_tol = conv.redchi2_rel_tol,
                dx_rel_tol = stall_dx_rel_tol,
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
        S_e_inv = _make_Se_inv(y_curr, use_band_snr, band_snr_coeffs, meas_sigma)
        chi2_curr = dot(y_obs .- y_curr, S_e_inv * (y_obs .- y_curr))
        push!(redchi2_hist, chi2_curr / dof_loop)
        if dx_rel < conv.dx_rel_tol ||
           rmse_rel_change < conv.rmse_rel_tol ||
           rmse_abs_change < conv.rmse_abs_tol
            converged = true
            status = Int16(1)
            break
        end
    end
    if status == 0 && !converged
        status = Int16(0)
    end
    x_out .= x_curr
    resid = y_obs .- y_curr
    rmse = sqrt(mean(resid .^ 2))
    # Final reduced χ² / objective always use S_e(y_curr) at the returned state
    S_e_inv = _make_Se_inv(y_curr, use_band_snr, band_snr_coeffs, meas_sigma)
    chi2_curr = dot(resid, S_e_inv * resid)
    obj = _cost_with_prior(y_obs, y_curr, x_curr, x_a_loc, S_e_inv, S_a_inv)
    # obtain the final Jacobian
    J_final = jac_eval(x_curr)
    # compute the final Hessian
    H_obs_final = J_final' * S_e_inv * J_final
    # compute current posterior sigma
    S_post = inv(Matrix(H_obs_final + S_a_inv))
    # Averaging-kernel trace: tr(A), A = S_post * H_obs.
    trace_A = Float64(tr(S_post * H_obs_final))
    # correction to final dof
    dof     = length(y_curr) - trace_A
    rchi2   = chi2_curr / max(dof, eps(Float64))
    return (
        converged = converged, 
        status = status, 
        n_steps = n_acc, 
        rmse = rmse, 
        reduced_chi2 = rchi2, 
        objective = obj,
        S_posterior = S_post,
        dof = dof,
        )
end
