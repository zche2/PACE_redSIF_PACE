# Batched SVD forward model: y[:,t] depends only on x[:,t] (no cross-talk between columns).
# Used for GPU-accelerated tiles and for tests vs scalar `make_svd_forward_model_λ`.

"""
    predict_svd_batched(x, solar_eff, PCs, SIF, leg_basis, log_transform, layout)

- `x`: `(n_state, n_tile)` state matrix
- `solar_eff`: `(nλ, n_tile)` effective solar factor per band and pixel
- `PCs`: `(nλ, n_pc)`, `SIF`: `(nλ, n_ev)`, `leg_basis`: `(nλ, n_leg_coeff)`

Returns `y` with size `(nλ, n_tile)` matching `make_svd_forward_model_λ` / `fm_svd` for each column.
"""
function predict_svd_batched(
    x::AbstractMatrix{T},
    solar_eff::AbstractMatrix{T},
    PCs::AbstractMatrix{T},
    SIF::AbstractMatrix{T},
    leg_basis::AbstractMatrix{T},
    log_transform::Bool,
    layout,
) where {T<:AbstractFloat}
    n_state, n_tile = size(x)
    n_state == layout.n_state || error("x rows $n_state != layout.n_state $(layout.n_state)")
    size(solar_eff, 2) == n_tile || error("solar_eff columns $(size(solar_eff,2)) != n_tile $n_tile")
    nλ = size(PCs, 1)
    size(solar_eff, 1) == nλ || error("solar_eff rows mismatch PCs rows")
    C = x[layout.idx_pc, :]                 # (n_pc, n_tile)
    ar = reshape(x[layout.idx_alpha, :], 1, n_tile)
    α = @. T(10.0) / (T(1.0) + exp(-ar)) + T(1.0)
    L = x[layout.idx_legendre, :]           # (n_leg_coeff, n_tile)
    S = x[layout.idx_sif, :]                # (n_ev, n_tile)
    trans_up = log_transform ? exp.(PCs * C) : T(1.0) .+ PCs * C
    tu = max.(trans_up, eps(T))
    tud = exp.(α .* log.(tu))
    rho = leg_basis * L
    sif_toa = trans_up .* (SIF * S)
    return @. solar_eff * tud * rho / T(π) + sif_toa
end

"""Central-difference Jacobian `(nλ, n_state, n_tile)` using batched forward."""
function jacobian_svd_batched_fd(
    x::AbstractMatrix{T},
    solar_eff::AbstractMatrix{T},
    PCs::AbstractMatrix{T},
    SIF::AbstractMatrix{T},
    leg_basis::AbstractMatrix{T},
    log_transform::Bool,
    layout;
    ε::T = T(sqrt(eps(T))),
) where {T<:AbstractFloat}
    n_state, n_tile = size(x)
    nλ = size(PCs, 1)
    y0 = predict_svd_batched(x, solar_eff, PCs, SIF, leg_basis, log_transform, layout)
    J = zeros(T, nλ, n_state, n_tile)
    xp = copy(x)
    for k in 1:n_state
        @views xp[k, :] .= x[k, :]
        xp[k, :] .+= ε
        yp = predict_svd_batched(xp, solar_eff, PCs, SIF, leg_basis, log_transform, layout)
        J[:, k, :] .= (yp .- y0) ./ ε
        @views xp[k, :] .= x[k, :]
    end
    return J, y0
end
