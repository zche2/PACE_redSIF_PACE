using LinearAlgebra: mul!

struct ForwardScratch{T<:Real}
    spectral_axis::Vector{T}
    z_hres::Vector{T}
    solar_hres::Vector{T}
    sif_basis::Matrix{T}
    K::Matrix{T}
    leg_basis::Matrix{T}
    xs_o2::Vector{T}
    xs_h2o::Vector{T}
    vcd_o2::Vector{T}
    vcd_h2o::Vector{T}
    trans::Vector{T}
    trans_sif::Vector{T}
    sif_hres::Vector{T}
    y_hres::Vector{T}
    y_lres::Vector{T}
    poly::Vector{T}
    out::Vector{T}
end

function _build_scratch(
    ::Type{T},
    spectral_axis,
    z_hres,
    solar_hres,
    sif_basis_hres,
    K,
    leg_basis,
) where {T<:Real}
    n_hres = length(z_hres)
    n_lres = size(K, 1)
    return ForwardScratch{T}(
        T.(spectral_axis),
        T.(z_hres),
        T.(solar_hres),
        T.(sif_basis_hres),
        T.(K),
        T.(leg_basis),
        zeros(T, n_hres),
        zeros(T, n_hres),
        zeros(T, n_hres),
        zeros(T, n_hres),
        zeros(T, n_hres),
        zeros(T, n_hres),
        zeros(T, n_hres),
        zeros(T, n_hres),
        zeros(T, n_lres),
        zeros(T, n_lres),
        zeros(T, n_lres),
    )
end

@inline function _eval_xs!(
    dest::AbstractVector,
    sitp,
    spectral_axis::AbstractVector,
    p,
    t,
)
    # In-place broadcast keeps this path allocation-free and easy to read.
    dest .= sitp.(spectral_axis, p, t)
    return dest
end

function _forward_generic(
    ctx,
    spectral_axis,
    z_hres,
    solar_hres,
    sif_basis_hres,
    K,
    leg_basis,
    vcd_o2_intercept,
    vcd_o2_slope,
    vcd_h2o_intercept,
    vcd_h2o_slope,
    p_o2_hpa,
    t_o2_k,
    p_h2o_hpa,
    t_h2o_k,
    vcd_o2_sif,
    vcd_h2o_sif,
    sif_coeff,
    leg_coeff,
)
    xs_o2 = vec(ctx.o2_sitp(spectral_axis, p_o2_hpa, t_o2_k))
    xs_h2o = vec(ctx.h2o_sitp(spectral_axis, p_h2o_hpa, t_h2o_k))

    vcd_o2_λ = @. vcd_o2_intercept + vcd_o2_slope * z_hres
    vcd_h2o_λ = @. vcd_h2o_intercept + vcd_h2o_slope * z_hres
    trans = @. exp(-(vcd_h2o_λ * xs_h2o + vcd_o2_λ * xs_o2))
    trans_sif = @. exp(-(vcd_h2o_sif * xs_h2o + vcd_o2_sif * xs_o2))
    sif_hres = if length(sif_coeff) > 0
        sif_basis_hres * sif_coeff
    else
        zeros(eltype(sif_basis_hres), size(sif_basis_hres, 1))
    end
    y_hres = @. solar_hres * trans + trans_sif * sif_hres

    y_lres = K * y_hres
    poly = leg_basis * leg_coeff
    return y_lres .* poly
end

function _forward_prealloc!(
    sc::ForwardScratch{T},
    ctx,
    vcd_o2_intercept,
    vcd_o2_slope,
    vcd_h2o_intercept,
    vcd_h2o_slope,
    p_o2_hpa,
    t_o2_k,
    p_h2o_hpa,
    t_h2o_k,
    vcd_o2_sif,
    vcd_h2o_sif,
    sif_coeff,
    leg_coeff,
) where {T<:Real}
    p_o2 = T(p_o2_hpa)
    t_o2 = T(t_o2_k)
    p_h2o = T(p_h2o_hpa)
    t_h2o = T(t_h2o_k)
    vcd_o2_sif_t = T(vcd_o2_sif)
    vcd_h2o_sif_t = T(vcd_h2o_sif)
    vcd_o2_i = T(vcd_o2_intercept)
    vcd_o2_s = T(vcd_o2_slope)
    vcd_h2o_i = T(vcd_h2o_intercept)
    vcd_h2o_s = T(vcd_h2o_slope)
    _eval_xs!(sc.xs_o2, ctx.o2_sitp, sc.spectral_axis, p_o2, t_o2)
    _eval_xs!(sc.xs_h2o, ctx.h2o_sitp, sc.spectral_axis, p_h2o, t_h2o)

    @. sc.vcd_o2 = vcd_o2_i + vcd_o2_s * sc.z_hres
    @. sc.vcd_h2o = vcd_h2o_i + vcd_h2o_s * sc.z_hres
    @. sc.trans = exp(-(sc.vcd_h2o * sc.xs_h2o + sc.vcd_o2 * sc.xs_o2))
    @. sc.trans_sif = exp(-(vcd_h2o_sif_t * sc.xs_h2o + vcd_o2_sif_t * sc.xs_o2))

    if length(sif_coeff) > 0
        mul!(sc.sif_hres, sc.sif_basis, sif_coeff)
    else
        fill!(sc.sif_hres, 0.0)  # Zero out SIF when not retrieving it
    end

    @. sc.y_hres = sc.solar_hres * sc.trans + sc.trans_sif * sc.sif_hres

    mul!(sc.y_lres, sc.K, sc.y_hres)
    mul!(sc.poly, sc.leg_basis, leg_coeff)
    @. sc.out = sc.y_lres * sc.poly
    return sc.out
end

"""
    make_forward_model_simple(
        ctx,
        solar_hres;
        n_legendre=2,
        preallocate_float64=false,
        preallocate_float32=false,
        preallocate_other_types=false,
    )

Build a closure `f(x)` that captures fixed context once:
- `ctx.λ_hres`, `ctx.spectral_axis`
- `ctx.o2_sitp`, `ctx.h2o_sitp`
- `ctx.kernel.RSR_out`
- `ctx.sif_basis_hres`
- `solar_hres`

This is the preferred pattern for performance and automatic differentiation.
If `preallocate_float64=true` and/or `preallocate_float32=true`, corresponding
floating-point calls reuse internal work buffers to reduce allocations;
if `preallocate_other_types=true`, additional numeric element types
(such as `ForwardDiff.Dual`) also get cached scratch buffers.

Supports both SIF retrieval (sif_nev > 0) and no-SIF mode (sif_nev = 0).

State vector layout for `x` when sif_nev > 0:
1. `vcd_o2_intercept`
2. `vcd_o2_slope`
3. `vcd_h2o_intercept`
4. `vcd_h2o_slope`
5. `p_o2_hpa`
6. `t_o2_k`
7. `p_h2o_hpa`
8. `t_h2o_k`
9. `vcd_o2_sif` (SIF-path O2 VCD)
10. `vcd_h2o_sif` (SIF-path H2O VCD)
11..`10+sif_nev`: SIF EV coefficients (on high-res grid)
`(11+sif_nev)`..end: Legendre coefficients `a0..aN` with `N=n_legendre`

State vector layout for `x` when sif_nev = 0:
1-8: Same as above (VCDs, p/T)
9..end: Legendre coefficients `a0..aN` with `N=n_legendre`

Low-res multiplicative polynomial:
`poly(λ) = Σ a_n P_n(z(λ))`, where `z` is λ mapped to [-1, 1].
"""
function make_forward_model_simple(
    ctx,
    solar_hres::AbstractVector{<:Real};
    n_legendre::Int=2,
    preallocate_float64::Bool=false,
    preallocate_float32::Bool=false,
    preallocate_other_types::Bool=false,
)
    n_legendre >= 0 || error("n_legendre must be >= 0")

    λ_hres = collect(ctx.λ_hres)
    λ_lres = collect(ctx.λ)
    spectral_axis = collect(ctx.spectral_axis)
    # check spectral_axis
    println("extrema of spectral_axis: $(spectral_axis[1]), $(spectral_axis[end])")
    K = hasproperty(ctx, :kernel_rsr_out) ? ctx.kernel_rsr_out : ctx.kernel.RSR_out
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

    # Type-specific scratch buffers for allocation-light forward evaluations.
    # This closure is stateful and not thread-safe when preallocation is enabled.
    scratch64 = preallocate_float64 ?
        _build_scratch(Float64, spectral_axis, z_hres, solar_hres, sif_basis_hres, K, leg_basis) :
        nothing
    scratch32 = preallocate_float32 ?
        _build_scratch(Float32, spectral_axis, z_hres, solar_hres, sif_basis_hres, K, leg_basis) :
        nothing
    # Optional cache for non-Float32/64 state types (e.g. ForwardDiff.Dual).
    # This keeps AD paths cleaner without forcing the default behavior.
    scratch_other = preallocate_other_types ? Dict{DataType, Any}() : nothing

    function _scratch_for_type(::Type{T}) where {T}
        T === Float64 && return scratch64
        T === Float32 && return scratch32
        preallocate_other_types || return nothing
        @assert !isnothing(scratch_other)
        haskey(scratch_other, T) && return scratch_other[T]
        sc = try
            _build_scratch(T, spectral_axis, z_hres, solar_hres, sif_basis_hres, K, leg_basis)
        catch
            nothing
        end
        scratch_other[T] = sc
        return sc
    end

    function forward_model_inner(x::AbstractVector)
        length(x) == layout.n_state ||
            error("x length must be exactly $(layout.n_state) (sif_nev=$(layout.sif_nev), n_legendre=$(layout.n_legendre))")

        # Always extract core atmospheric parameters
        vcd_o2_intercept = x[layout.idx_vcd_o2_intercept]
        vcd_o2_slope = x[layout.idx_vcd_o2_slope]
        vcd_h2o_intercept = x[layout.idx_vcd_h2o_intercept]
        vcd_h2o_slope = x[layout.idx_vcd_h2o_slope]
        p_o2_hpa = x[layout.idx_p_o2_hpa]
        t_o2_k = x[layout.idx_t_o2_k]
        p_h2o_hpa = x[layout.idx_p_h2o_hpa]
        t_h2o_k = x[layout.idx_t_h2o_k]
        
        # Extract SIF-related parameters only if sif_nev > 0
        if n_ev > 0
            vcd_o2_sif = x[layout.idx_vcd_o2_sif]
            vcd_h2o_sif = x[layout.idx_vcd_h2o_sif]
            sif_coeff = @view x[layout.idx_sif]
        else
            # Set to zero/empty for no-SIF mode
            vcd_o2_sif = zero(eltype(x))
            vcd_h2o_sif = zero(eltype(x))
            sif_coeff = view(x, 1:0)  # Empty view
        end
        
        leg_coeff = @view x[layout.idx_legendre]

        sc = _scratch_for_type(eltype(x))
        if isnothing(sc)
            return _forward_generic(
                ctx,
                spectral_axis,
                z_hres,
                solar_hres,
                sif_basis_hres,
                K,
                leg_basis,
                vcd_o2_intercept,
                vcd_o2_slope,
                vcd_h2o_intercept,
                vcd_h2o_slope,
                p_o2_hpa,
                t_o2_k,
                p_h2o_hpa,
                t_h2o_k,
                vcd_o2_sif,
                vcd_h2o_sif,
                sif_coeff,
                leg_coeff,
            )
        end
        return _forward_prealloc!(
            sc,
            ctx,
            vcd_o2_intercept,
            vcd_o2_slope,
            vcd_h2o_intercept,
            vcd_h2o_slope,
            p_o2_hpa,
            t_o2_k,
            p_h2o_hpa,
            t_h2o_k,
            vcd_o2_sif,
            vcd_h2o_sif,
            sif_coeff,
            leg_coeff,
        )
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
The SIF part is always linked to `ctx.n_ev` (number of SIF eigenvectors).
If `n_ev = 0`, SIF-path VCDs and SIF coefficients are excluded from the state.
"""
function state_layout_simple(ctx; n_legendre::Int=2)
    n_ev = ctx.sif_nev
    
    # Core state: VCD intercepts, slopes, p/T for both gases, SIF-path VCDs
    n_core = 8  # vcd_o2_i, vcd_o2_s, vcd_h2o_i, vcd_h2o_s, p_o2, t_o2, p_h2o, t_h2o
    n_sif_path = 2  # vcd_o2_sif, vcd_h2o_sif
    n_leg_coeff = n_legendre + 1  # Legendre coefficients a0..aN
    # State size depends on whether SIF is included
    if n_ev > 0
        n_state = n_core + n_sif_path + n_ev + n_leg_coeff
    else
        n_state = n_core + n_leg_coeff
    end
    
    # Build index map
    idx = 1
    idx_vcd_o2_intercept = idx; idx += 1
    idx_vcd_o2_slope = idx; idx += 1
    idx_vcd_h2o_intercept = idx; idx += 1
    idx_vcd_h2o_slope = idx; idx += 1
    idx_p_o2_hpa = idx; idx += 1
    idx_t_o2_k = idx; idx += 1
    idx_p_h2o_hpa = idx; idx += 1
    idx_t_h2o_k = idx; idx += 1
    
    # SIF-path VCDs only if retrieving SIF
    if n_ev > 0
        idx_vcd_o2_sif = idx; idx += 1
        idx_vcd_h2o_sif = idx; idx += 1
        idx_sif = idx:(idx + n_ev - 1); idx += n_ev
    else
        idx_vcd_o2_sif = 0  # Placeholder
        idx_vcd_h2o_sif = 0
        idx_sif = 1:0  # Empty range
    end
    
    idx_legendre = idx:(idx + n_legendre)
    
    return (
        n_state = n_state,
        n_ev = n_ev,
        n_legendre = n_legendre,
        n_leg_coeff = n_leg_coeff,
        idx_vcd_o2_intercept = idx_vcd_o2_intercept,
        idx_vcd_o2_slope = idx_vcd_o2_slope,
        idx_vcd_h2o_intercept = idx_vcd_h2o_intercept,
        idx_vcd_h2o_slope = idx_vcd_h2o_slope,
        idx_p_o2_hpa = idx_p_o2_hpa,
        idx_t_o2_k = idx_t_o2_k,
        idx_p_h2o_hpa = idx_p_h2o_hpa,
        idx_t_h2o_k = idx_t_h2o_k,
        idx_vcd_o2_sif = idx_vcd_o2_sif,
        idx_vcd_h2o_sif = idx_vcd_h2o_sif,
        idx_sif = idx_sif,
        idx_legendre = idx_legendre,
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
Supports both SIF retrieval (n_ev > 0) and no-SIF mode (n_ev = 0).
"""
function state_names_simple(ctx; n_legendre::Int=2)
    layout = state_layout_simple(ctx; n_legendre=n_legendre)
    names = Vector{String}(undef, layout.n_state)
    
    names[layout.idx_vcd_o2_intercept] = "vcd_o2_intercept"
    names[layout.idx_vcd_o2_slope] = "vcd_o2_slope"
    names[layout.idx_vcd_h2o_intercept] = "vcd_h2o_intercept"
    names[layout.idx_vcd_h2o_slope] = "vcd_h2o_slope"
    names[layout.idx_p_o2_hpa] = "p_o2_hpa"
    names[layout.idx_t_o2_k] = "t_o2_k"
    names[layout.idx_p_h2o_hpa] = "p_h2o_hpa"
    names[layout.idx_t_h2o_k] = "t_h2o_k"
    
    if layout.n_ev > 0
        names[layout.idx_vcd_o2_sif] = "vcd_o2_sif"
        names[layout.idx_vcd_h2o_sif] = "vcd_h2o_sif"
        
        for i in 1:layout.n_ev
            names[layout.idx_sif[i]] = "sif_ev$(i)"
        end
    end
    
    for i in 1:layout.n_leg_coeff
        names[layout.idx_legendre[i]] = "legendre_p$(i-1)"
    end
    
    return names
end

"""
    initial_state_simple(ctx; n_legendre=2, T::Type{<:AbstractFloat}=Float64)

Build a default initial state matching current model layout.
Supports both SIF retrieval (n_ev > 0) and no-SIF mode (n_ev = 0).
"""
function initial_state_simple(
    ctx;
    n_legendre::Int=2,
    T::Type{<:AbstractFloat}=Float64,
)
    layout = state_layout_simple(ctx; n_legendre=n_legendre)
    x = zeros(T, layout.n_state)
    
    # VCD estimates
    x[layout.idx_vcd_o2_intercept] = T(3.981071705534985e24)   # VCD_O2 intercept
    x[layout.idx_vcd_o2_slope] = T(0.0)                        # VCD_O2 slope
    x[layout.idx_vcd_h2o_intercept] = T(3.981071705534969e22)  # VCD_H2O intercept
    x[layout.idx_vcd_h2o_slope] = T(0.0)                       # VCD_H2O slope
    
    # Pressure and temperature
    x[layout.idx_p_o2_hpa] = T(800.0)     # p_o2_hpa
    x[layout.idx_t_o2_k] = T(250.0)       # t_o2_k
    x[layout.idx_p_h2o_hpa] = T(780.0)    # p_h2o_hpa
    x[layout.idx_t_h2o_k] = T(285.0)      # t_h2o_k
    
    # SIF-path VCDs and coefficients (only if retrieving SIF)
    if layout.n_ev > 0
        # Shorter effective light path for emitted SIF than for incoming solar
        x[layout.idx_vcd_o2_sif] = T(0.5) * x[layout.idx_vcd_o2_intercept]
        x[layout.idx_vcd_h2o_sif] = T(0.5) * x[layout.idx_vcd_h2o_intercept]
        
        # SIF coefficients start at zero
        x[layout.idx_sif] .= T(0.0)
    end
    
    # Legendre polynomial coefficients
    x[layout.idx_legendre] .= T(0.0)
    x[first(layout.idx_legendre)] = T(1.0)  # P0 = 1 (baseline multiplicative factor)
    
    return x
end

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
    # Explicit center-normalization: z = 0 at λcenter.
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
        # (l+1)P_{l+1} = (2l+1)zP_l - lP_{l-1}
        @. A[:, l + 2] = (
            T(2l + 1) * zt * A[:, l + 1] - T(l) * A[:, l]
        ) / T(l + 1)
    end
    return A
end
