"""
    make_forward_model_simple(ctx, solar_hres; n_legendre=2)

Build a closure `f(x)` that captures fixed context once:
- `ctx.λ_hres`, `ctx.spectral_axis`
- `ctx.o2_sitp`, `ctx.h2o_sitp`
- `ctx.kernel.RSR_out`
- `ctx.sif_basis_hres`
- `solar_hres`

This is the preferred pattern for performance and automatic differentiation.

State vector layout for `x`:
1. `vcd_o2_intercept`
2. `vcd_o2_slope`
3. `vcd_h2o_intercept`
4. `vcd_h2o_slope`
5. `p_o2_hpa`
6. `t_o2_k`
7. `p_h2o_hpa`
8. `t_h2o_k`
9. `continuum_scale`
10..`9+nEV`: SIF EV coefficients (on high-res grid)
`(10+nEV)`..end: Legendre coefficients `a0..aN` with `N=n_legendre`

Low-res multiplicative polynomial:
`poly(λ) = Σ a_n P_n(z(λ))`, where `z` is λ mapped to [-1, 1].
"""
function make_forward_model_simple(
    ctx,
    solar_hres::AbstractVector{<:Real};
    n_legendre::Int=2,
)
    n_legendre >= 0 || error("n_legendre must be >= 0")

    λ_hres = collect(ctx.λ_hres)
    λ_lres = collect(ctx.λ)
    spectral_axis = collect(ctx.spectral_axis)
    K = ctx.kernel.RSR_out
    sif_basis_hres = ctx.sif_basis_hres

    length(solar_hres) == length(λ_hres) ||
        error("solar_hres length must match ctx.λ_hres length")
    size(K, 2) == length(λ_hres) ||
        error("kernel columns must match ctx.λ_hres length")
    size(sif_basis_hres, 1) == length(λ_hres) ||
        error("sif_basis_hres first dimension must match ctx.λ_hres length")

    layout = state_layout_simple(ctx; n_legendre=n_legendre)
    n_ev = layout.n_ev
    z_hres = _normalized_grid(λ_hres)
    z_lres = _normalized_grid(λ_lres)
    leg_basis = _legendre_design_matrix(z_lres, n_legendre)

    function forward_model_inner(x::AbstractVector)
        length(x) == layout.n_state ||
            error("x length must be exactly $(layout.n_state) (nEV=$(layout.n_ev), n_legendre=$(layout.n_legendre))")

        vcd_o2_intercept = x[layout.idx_vcd_o2_intercept]
        vcd_o2_slope = x[layout.idx_vcd_o2_slope]
        vcd_h2o_intercept = x[layout.idx_vcd_h2o_intercept]
        vcd_h2o_slope = x[layout.idx_vcd_h2o_slope]
        p_o2_hpa = x[layout.idx_p_o2_hpa]
        t_o2_k = x[layout.idx_t_o2_k]
        p_h2o_hpa = x[layout.idx_p_h2o_hpa]
        t_h2o_k = x[layout.idx_t_h2o_k]
        continuum_scale = x[layout.idx_continuum_scale]
        sif_coeff = @view x[layout.idx_sif]
        leg_coeff = @view x[layout.idx_legendre]

        xs_o2 = vec(ctx.o2_sitp(spectral_axis, p_o2_hpa, t_o2_k))
        xs_h2o = vec(ctx.h2o_sitp(spectral_axis, p_h2o_hpa, t_h2o_k))

        vcd_o2_λ = @. vcd_o2_intercept + vcd_o2_slope * z_hres
        vcd_h2o_λ = @. vcd_h2o_intercept + vcd_h2o_slope * z_hres
        trans = @. exp(-(vcd_h2o_λ * xs_h2o + vcd_o2_λ * xs_o2))
        sif_hres = sif_basis_hres * sif_coeff
        y_hres = @. continuum_scale * solar_hres * trans + sif_hres

        y_lres = K * y_hres
        poly = leg_basis * leg_coeff
        return y_lres .* poly
    end

    return forward_model_inner
end

"""
    forward_model_simple(x; ctx, solar_hres, kwargs...)

Convenience wrapper that builds and evaluates the model in one call.
For repeated calls (and AD/Jacobians), prefer:

```julia
fm = make_forward_model_simple(ctx, solar_hres)
y = fm(x)
```
"""
function forward_model_simple(
    x::AbstractVector;
    ctx,
    solar_hres::AbstractVector{<:Real},
    kwargs...,
)
    fm = make_forward_model_simple(ctx, solar_hres; kwargs...)
    return fm(x)
end

"""
    state_layout_simple(ctx; n_legendre=2)

Return a named tuple describing the full state-vector layout.
The SIF part is always linked to `size(ctx.sif_basis_hres, 2)`.
"""
function state_layout_simple(ctx; n_legendre::Int=2)
    n_legendre >= 0 || error("n_legendre must be >= 0")
    n_ev = size(ctx.sif_basis_hres, 2)
    n_leg_coeff = n_legendre + 1
    n_state = 9 + n_ev + n_leg_coeff
    return (
        n_ev = n_ev,
        n_legendre = n_legendre,
        n_leg_coeff = n_leg_coeff,
        n_state = n_state,
        idx_vcd_o2_intercept = 1,
        idx_vcd_o2_slope = 2,
        idx_vcd_h2o_intercept = 3,
        idx_vcd_h2o_slope = 4,
        idx_p_o2_hpa = 5,
        idx_t_o2_k = 6,
        idx_p_h2o_hpa = 7,
        idx_t_h2o_k = 8,
        idx_continuum_scale = 9,
        idx_sif = 10:(9 + n_ev),
        idx_legendre = (10 + n_ev):(9 + n_ev + n_leg_coeff),
    )
end

"""
    state_length_simple(ctx; n_legendre=2)

Return required state vector length for the toy model.
"""
function state_length_simple(ctx; n_legendre::Int=2)
    return state_layout_simple(ctx; n_legendre=n_legendre).n_state
end

"""
    state_names_simple(ctx; n_legendre=2)

Return ordered state names matching `state_layout_simple`.
"""
function state_names_simple(ctx; n_legendre::Int=2)
    layout = state_layout_simple(ctx; n_legendre=n_legendre)
    names = String[
        "vcd_o2_intercept",
        "vcd_o2_slope",
        "vcd_h2o_intercept",
        "vcd_h2o_slope",
        "p_o2_hpa",
        "t_o2_k",
        "p_h2o_hpa",
        "t_h2o_k",
        "continuum_scale",
    ]
    append!(names, ["sif_ev$(i)" for i in 1:layout.n_ev])
    append!(names, ["legendre_p$(i - 1)" for i in 1:layout.n_leg_coeff])
    return names
end

"""
    initial_state_simple(ctx; n_legendre=2)

Build a default initial state matching current model layout.
"""
function initial_state_simple(ctx; n_legendre::Int=2)
    layout = state_layout_simple(ctx; n_legendre=n_legendre)
    x0 = zeros(Float64, layout.n_state)
    x0[layout.idx_vcd_o2_intercept] = 3.981071705534985e24   # VCD_O2 intercept
    x0[layout.idx_vcd_o2_slope] = 0.0                        # VCD_O2 slope
    x0[layout.idx_vcd_h2o_intercept] = 3.981071705534969e22  # VCD_H2O intercept
    x0[layout.idx_vcd_h2o_slope] = 0.0                       # VCD_H2O slope
    x0[layout.idx_p_o2_hpa] = 800.0     # p_o2_hpa
    x0[layout.idx_t_o2_k] = 250.0       # t_o2_k
    x0[layout.idx_p_h2o_hpa] = 780.0    # p_h2o_hpa
    x0[layout.idx_t_h2o_k] = 285.0      # t_h2o_k
    x0[layout.idx_continuum_scale] = 1.0 # continuum_scale
    # SIF coeffs default to zero.
    # For multiplicative Legendre polynomial without +1 term, set P0 coeff = 1
    # so the default polynomial factor is unity.
    x0[first(layout.idx_legendre)] = 1.0
    return x0
end

function _normalized_grid(λ::AbstractVector{<:Real})
    λf = collect(Float64.(λ))
    λmin, λmax = extrema(λf)
    λspan = λmax - λmin
    if λspan == 0
        return zeros(Float64, length(λf))
    end
    λcenter = 0.5 * (λmin + λmax)
    halfspan = 0.5 * λspan
    # Explicit center-normalization: z = 0 at λcenter.
    return @. (λf - λcenter) / halfspan
end

function _legendre_design_matrix(z::AbstractVector{<:Real}, n_legendre::Int)
    n = length(z)
    A = zeros(Float64, n, n_legendre + 1)
    A[:, 1] .= 1.0
    n_legendre == 0 && return A
    A[:, 2] .= z
    for l in 1:(n_legendre - 1)
        # (l+1)P_{l+1} = (2l+1)zP_l - lP_{l-1}
        A[:, l + 2] .= ((2l + 1) .* z .* A[:, l + 1] .- l .* A[:, l]) ./ (l + 1)
    end
    return A
end
