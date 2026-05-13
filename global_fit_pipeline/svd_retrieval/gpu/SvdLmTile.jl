# Tile LM: batched forward + batched FD Jacobian on active columns each outer iteration;
# parallel `lm_one_step_from_J` per active column. Optional CUDA for batched predict/J.

using LinearAlgebra
using SparseArrays
using Statistics

using CUDA

include(joinpath(@__DIR__, "SvdFmBatched.jl"))

function _to_gpu_or_cpu(A::AbstractMatrix{T}, use_cuda::Bool) where {T}
    (use_cuda && CUDA.functional()) || return A
    return CuArray(T.(A))
end

function predict_svd_batched_maybe_gpu(
    x::Matrix{Float64},
    solar_eff::Matrix{Float64},
    PCs_g,
    SIF_g,
    leg_g,
    log_transform::Bool,
    layout,
    use_cuda::Bool,
)
    if use_cuda && CUDA.functional()
        xc = CuArray(x)
        sc = CuArray(solar_eff)
        y = predict_svd_batched(xc, sc, PCs_g, SIF_g, leg_g, log_transform, layout)
        return Array(y)
    end
    return predict_svd_batched(x, solar_eff, PCs_g, SIF_g, leg_g, log_transform, layout)
end

function jacobian_svd_batched_fd_maybe_gpu(
    x::Matrix{Float64},
    solar_eff::Matrix{Float64},
    PCs_g,
    SIF_g,
    leg_g,
    log_transform::Bool,
    layout,
    use_cuda::Bool,
    ε::Float64,
)
    if use_cuda && CUDA.functional()
        n_state, n_tile = size(x)
        nλ = size(PCs_g, 1)
        xc = CuArray(x)
        sc = CuArray(solar_eff)
        y0 = predict_svd_batched(xc, sc, PCs_g, SIF_g, leg_g, log_transform, layout)
        J = CUDA.zeros(Float64, nλ, n_state, n_tile)
        xp = copy(xc)
        for k in 1:n_state
            copyto!(xp, xc)
            @inbounds xp[k, :] .+= ε
            yp = predict_svd_batched(xp, sc, PCs_g, SIF_g, leg_g, log_transform, layout)
            J[:, k, :] .= (yp .- y0) ./ ε
        end
        return Array(J), Array(y0)
    end
    return jacobian_svd_batched_fd(
        x,
        solar_eff,
        PCs_g,
        SIF_g,
        leg_g,
        log_transform,
        layout;
        ε = ε,
    )
end

function _make_fm_column(
    t::Int,
    solar_eff::Matrix{Float64},
    PCs::Matrix{Float64},
    SIF::Matrix{Float64},
    leg_basis::Matrix{Float64},
    log_transform::Bool,
    layout,
)
    return function fm(xv::AbstractVector)
        X = reshape(copy(xv), :, 1)
        se = reshape(view(solar_eff, :, t), :, 1)
        return vec(predict_svd_batched(X, se, PCs, SIF, leg_basis, log_transform, layout))
    end
end

"""
    run_tile_svd_retrieval!(x_out, y_obs, solar_eff, PCs, SIF, leg_basis, log_transform, layout, ...)

`y_obs`, `solar_eff`: `(nλ, K)`. Matches per-column semantics of `_run_one_svd_retrieval!` using
batched `Y`/`J` for active columns each outer iteration.
"""
function run_tile_svd_retrieval!(
    x_out::Matrix{Float64},
    y_obs::Matrix{Float64},
    solar_eff::Matrix{Float64},
    PCs::Matrix{Float64},
    SIF::Matrix{Float64},
    leg_basis::Matrix{Float64},
    log_transform::Bool,
    layout,
    x0::Vector{Float64},
    prior_sigma_template::Vector{Float64},
    lower::Vector{Float64},
    upper::Vector{Float64},
    x_scale_template::Vector{Float64},
    use_leg01::Bool,
    A01::Matrix{Float64},
    prior_min_sigma::Float64,
    leg01_frac::Float64,
    use_band_snr::Bool,
    band_snr_coeffs,
    meas_sigma::Float64,
    lm::NamedTuple,
    conv::NamedTuple,
    stall_dx_rel_tol::Float64,
    max_outer_steps::Int;
    use_cuda::Bool = false,
)
    K = size(y_obs, 2)
    n_state = layout.n_state
    nλ = size(y_obs, 1)
    size(x_out) == (n_state, K) || error("x_out size mismatch")
    x_a_mat = repeat(reshape(x0, :, 1), 1, K)
    σ_mat = repeat(reshape(prior_sigma_template, :, 1), 1, K)
    if use_leg01 && length(layout.idx_legendre) >= 1
        for t in 1:K
            fm0 = _make_fm_column(t, solar_eff, PCs, SIF, leg_basis, log_transform, layout)
            y0 = fm0(x0)
            ratio = y_obs[:, t] ./ max.(abs.(y0), eps(Float64))
            w = y_obs[:, t] .- minimum(y_obs[:, t])
            w .+= max(maximum(w), 1.0) * 1e-6
            sv = sqrt.(w ./ maximum(w))
            c01 = (A01 .* sv) \ (ratio .* sv)
            leg0 = first(layout.idx_legendre)
            x_a_mat[leg0, t] = c01[1]
            σ_mat[leg0, t] = max(abs(c01[1]) * leg01_frac, prior_min_sigma)
            if length(layout.idx_legendre) >= 2
                leg1 = layout.idx_legendre[2]
                x_a_mat[leg1, t] = c01[2]
                σ_mat[leg1, t] = max(abs(c01[2]) * leg01_frac, prior_min_sigma)
            end
        end
    end
    PCs_g = _to_gpu_or_cpu(PCs, use_cuda)
    SIF_g = _to_gpu_or_cpu(SIF, use_cuda)
    leg_g = _to_gpu_or_cpu(leg_basis, use_cuda)

    X = copy(x_a_mat)
    lam = fill(Float64(lm.lambda0), K)
    dof = max(nλ - n_state, 1)
    redchi2_hist = [Float64[] for _ in 1:K]
    dx_rel_hist = [Float64[] for _ in 1:K]
    Y_init = predict_svd_batched_maybe_gpu(X, solar_eff, PCs_g, SIF_g, leg_g, log_transform, layout, use_cuda)
    for t in 1:K
        Se0 = if use_band_snr && band_snr_coeffs !== nothing
            make_Se_inv_from_snr(Y_init[:, t], band_snr_coeffs)
        else
            spdiagm(0 => fill(1.0 / meas_sigma^2, nλ))
        end
        r0 = y_obs[:, t] .- Y_init[:, t]
        chi0 = dot(r0, Se0 * r0)
        push!(redchi2_hist[t], chi0 / dof)
    end
    active = trues(K)
    status = fill(Int16(0), K)
    n_acc = zeros(Int, K)
    converged = falses(K)
    rmse_prev = zeros(K)
    for t in 1:K
        rmse_prev[t] = sqrt(mean((y_obs[:, t] .- Y_init[:, t]) .^ 2))
    end
    ε = sqrt(eps(Float64))
    for _ in 1:max_outer_steps
        any(active) || break
        idx = findall(active)
        Xsub = X[:, idx]
        ssub = solar_eff[:, idx]
        Ysub, Jsub = jacobian_svd_batched_fd_maybe_gpu(
            Xsub,
            ssub,
            PCs_g,
            SIF_g,
            leg_g,
            log_transform,
            layout,
            use_cuda,
            ε,
        )
        X_next = copy(X)
        Base.Threads.@threads for k in eachindex(idx)
            t = idx[k]
            fm_t = _make_fm_column(t, solar_eff, PCs, SIF, leg_basis, log_transform, layout)
            S_a_inv_t = _spdiag_invvar(σ_mat[:, t])
            step = try
                lm_one_step_from_J(
                    fm_t,
                    X[:, t],
                    y_obs[:, t],
                    Ysub[:, k],
                    Jsub[:, :, k];
                    x_a = x_a_mat[:, t],
                    S_a_inv = S_a_inv_t,
                    lambda = lam[t],
                    lambda_up = lm.lambda_up,
                    lambda_down = lm.lambda_down,
                    lambda_min = lm.lambda_min,
                    lambda_max = lm.lambda_max,
                    max_inner = lm.max_inner,
                    x_scale = x_scale_template,
                    lower_bounds = lower,
                    upper_bounds = upper,
                    use_band_snr = use_band_snr,
                    band_snr_coeffs = band_snr_coeffs,
                    meas_sigma = meas_sigma,
                )
            catch
                status[t] = Int16(4)
                active[t] = false
                continue
            end
            lam[t] = step.lambda_next
            if !step.accepted
                stalled_conv, _ = _stalled_convergence(
                    dx_rel_hist[t],
                    redchi2_hist[t];
                    enabled = conv.enabled,
                    window = conv.window,
                    redchi2_target = conv.redchi2_target,
                    redchi2_abs_tol = conv.redchi2_abs_tol,
                    redchi2_rel_tol = conv.redchi2_rel_tol,
                    dx_rel_tol = stall_dx_rel_tol,
                )
                if stalled_conv
                    converged[t] = true
                    status[t] = Int16(1)
                else
                    status[t] = Int16(2)
                end
                active[t] = false
                continue
            end
            n_acc[t] += 1
            x_prev = X[:, t]
            X_next[:, t] .= step.x_next
            y_curr_t = step.y_next
            rmse_curr = sqrt(mean((y_obs[:, t] .- y_curr_t) .^ 2))
            dx_rel = norm(step.dx) / max(norm(x_prev), eps(Float64))
            push!(dx_rel_hist[t], dx_rel)
            rmse_abs_change = abs(rmse_curr - rmse_prev[t])
            rmse_rel_change = rmse_abs_change / max(abs(rmse_prev[t]), eps(Float64))
            rmse_prev[t] = rmse_curr
            Se_t = if use_band_snr && band_snr_coeffs !== nothing
                make_Se_inv_from_snr(y_curr_t, band_snr_coeffs)
            else
                spdiagm(0 => fill(1.0 / meas_sigma^2, nλ))
            end
            chi2_curr = dot(y_obs[:, t] .- y_curr_t, Se_t * (y_obs[:, t] .- y_curr_t))
            push!(redchi2_hist[t], chi2_curr / dof)
            if dx_rel < conv.dx_rel_tol ||
               rmse_rel_change < conv.rmse_rel_tol ||
               rmse_abs_change < conv.rmse_abs_tol
                converged[t] = true
                status[t] = Int16(1)
                active[t] = false
            end
        end
        X .= X_next
    end
    for t in 1:K
        if status[t] == Int16(0) && !converged[t]
            status[t] = Int16(0)
        end
    end
    x_out .= X
    return (
        converged = converged,
        status = status,
        n_steps = n_acc,
        x_a_mat = copy(x_a_mat),
        sigma_mat = copy(σ_mat),
    )
end

"""Optional CPU batched small solves (placeholder for future CUSOLVER batching)."""
function batched_cholesky_solve!(A_stack::AbstractArray{Float64,3}, rhs_mat::AbstractMatrix{Float64})
    nb = size(A_stack, 3)
    out = similar(rhs_mat)
    for b in 1:nb
        out[:, b] .= A_stack[:, :, b] \ rhs_mat[:, b]
    end
    return out
end
