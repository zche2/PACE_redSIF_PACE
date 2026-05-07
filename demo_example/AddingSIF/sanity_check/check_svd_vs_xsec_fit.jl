#!/usr/bin/env julia
# Compare LUT cross-section fit vs SVD transmittance fit on a single PACE pixel.
#
# SVD state: [PC coeffs | α | Legendre | SIF] — see `svd_state_layout` and `[fit.svd]` in the TOML.
#
# Config is read from svd_vs_xsec_config.toml (or set SVD_XSEC_CONFIG env var).
# Both retrievals use the same LM optimizer, PACE observation, and solar spectrum.
#
# Outputs in svd_vs_xsec_check/:
#   svd_basis.png              — first n_pc PCs vs wavelength with variance explained
#   rmse_vs_iteration.png      — RMSE trace for both methods on the same axes
#   spectra_comparison.png     — observed vs prior vs final for both methods + residuals
#   sif_trans_components.png   — convolved reflectance-path vs SIF-path: LUT xSec vs SVD
#   transmittance_one_two.png  — LUT: K×T vs K×T² on OCI bands; SVD: T_up / T_up² / T_up^α on OCI bands
#   jacobian_lut_hybrid_vs_forwarddiff.png — LUT: normalized (∂y/∂x_j)·σ_post,j vs λ (hybrid vs ForwardDiff)
#   jacobian_svd_forwarddiff.png          — SVD: normalized (∂y/∂x_j)·σ_post,j vs λ (ForwardDiff)
#   averaging_kernel_lut.png              — LUT averaging kernel A = GK (state×state, final state)
#   averaging_kernel_svd.png              — SVD averaging kernel A = GK (state×state, final state)
#   comparison_summary.csv     — n_iter, n_fwd, n_jac, RMSE, wall time
#
# Run:
#   julia --project=. demo_example/AddingSIF/sanity_check/check_svd_vs_xsec_fit.jl
# or:
#   SVD_XSEC_CONFIG=/path/to/config.toml julia --project=. ...

using TOML
using NCDatasets
using JLD2
using Interpolations
using LinearAlgebra
using SparseArrays
using Statistics
using Plots
using Printf

const _SANITY_DIR   = @__DIR__
const _DEMO_EXAMPLE = joinpath(_SANITY_DIR, "..", "..")
const _OUT_DIR      = joinpath(_SANITY_DIR, "svd_vs_xsec_check")

include(joinpath(_DEMO_EXAMPLE, "Fit_toy_forward_model.jl"))

# ── SVD forward model ─────────────────────────────────────────────────────────

"""
State-vector layout for the SVD transmittance forward model.

  [c_1 .. c_{n_pc} | α | leg_0 .. leg_{n_legendre} | sif_ev_0 .. sif_ev_{n_ev-1}]

`α` scales the upwelling transmittance to a down–up path: `trans_updown = trans_up^α`
(equivalently `exp(α * log(trans_up))`). `α = 1` recovers the PC transmittance only.
"""
function svd_state_layout(; n_pc::Int, n_legendre::Int, n_ev::Int)
    n_leg_coeff = n_legendre + 1
    n_alpha     = 1
    n_state     = n_pc + n_alpha + n_leg_coeff + n_ev
    i0 = n_pc + n_alpha
    return (
        n_pc         = n_pc,
        n_alpha      = n_alpha,
        n_legendre   = n_legendre,
        n_leg_coeff  = n_leg_coeff,
        n_ev         = n_ev,
        n_state      = n_state,
        idx_pc       = 1:n_pc,
        idx_alpha    = (n_pc + 1):(n_pc + n_alpha),
        idx_legendre = (i0 + 1):(i0 + n_leg_coeff),
        idx_sif      = (i0 + n_leg_coeff + 1):n_state,
    )
end

"""
Build the SVD transmittance forward model on the **observation band grid** (same as pseudo
measurements: PCs from the NetCDF live on OCI band centers, interpolated to `ctx.λ`).

State vector: `[c_1..c_{n_pc}, α, leg_0..leg_{n_leg}, sif_ev_0..sif_ev_{n_ev}]`

Per band (length = `length(ctx.λ)`):
  linear:  trans_up = 1 + PCs * c_vec
  log:     trans_up = exp.(PCs * c_vec)
  trans_updown = trans_up^α
  ρ = leg_basis_obs * leg_coeff
  y = (K*solar) * trans_updown * ρ / π + trans_up * (K * sif_basis_hres * c_sif)

This matches compiling the physics in low resolution like the pseudo-measurement pipeline.
It approximates `K * (y_hres(state))` by applying `K` only to `solar` and to the SIF shape
`(sif_basis_hres * c)` while using band-local `trans_up`; narrow bands keep the error small.
"""
function make_svd_forward_model(
    ctx,
    solar_hres::AbstractVector{<:Real},
    PCs_obs::Matrix{<:Real};
    n_pc::Int,
    n_legendre::Int,
    log_transform::Bool = false,
)
    n_pc > 0       || error("n_pc must be > 0")
    n_legendre >= 0 || error("n_legendre must be >= 0")
    size(PCs_obs, 2) >= n_pc || error("PCs_obs has fewer columns than n_pc")

    λ_obs         = collect(Float64.(ctx.λ))
    size(PCs_obs, 1) == length(λ_obs) ||
        error("PCs_obs rows ($(size(PCs_obs, 1))) must match ctx.λ length ($(length(λ_obs)))")
    K             = hasproperty(ctx, :kernel_rsr_out) ? ctx.kernel_rsr_out : ctx.kernel.RSR_out
    sif_basis_hres = Float64.(ctx.sif_basis_hres)
    n_ev          = size(sif_basis_hres, 2)
    solar         = Float64.(collect(solar_hres))
    PCs           = Float64.(PCs_obs[:, 1:n_pc])

    solar_lres = Vector(K * solar)
    sif_lres   = K * sif_basis_hres   # (n_bands × n_ev)

    z_obs     = _normalized_grid(λ_obs)
    leg_basis = _legendre_design_matrix(z_obs, n_legendre)
    layout    = svd_state_layout(; n_pc, n_legendre, n_ev)

    function fm_svd(x::AbstractVector)
        length(x) == layout.n_state ||
            error("SVD state vector length $(length(x)) ≠ expected $(layout.n_state)")
        c_vec        = @view x[layout.idx_pc]
        alpha_coeff  = 10.0 ./ (1.0 .+ exp(-x[first(layout.idx_alpha)])) + 1.0
        leg_coeff    = @view x[layout.idx_legendre]
        sif_coeff    = @view x[layout.idx_sif]

        trans_up     = log_transform ? exp.(PCs * c_vec) : 1.0 .+ PCs * c_vec
        trans_updown = exp.(alpha_coeff .* log.(max.(trans_up, eps(Float64))))
        rho_obs      = leg_basis * leg_coeff
        sif_toa_lres = trans_up .* (sif_lres * sif_coeff)

        return @.(solar_lres * trans_updown * rho_obs / π + sif_toa_lres)
    end

    return fm_svd, layout
end

"""
Split LUT xSec TOA radiance into reflectance-path and SIF-path contributions on the observation grid:

  y_refl_lres = K * (solar * trans * ρ / π)
  y_sif_lres  = K * (trans_sif * (sif_basis * c))

Matches `_forward_generic` in `toy_forward_model.jl`.
"""
function xsec_reflectance_and_sif_lres(
    ctx,
    solar_hres::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    layout,
    n_legendre::Int,
    model_variant::Symbol,
)
    λ_hres        = collect(Float64.(ctx.λ_hres))
    spectral_axis = collect(Float64.(ctx.spectral_axis))
    K             = hasproperty(ctx, :kernel_rsr_out) ? ctx.kernel_rsr_out : ctx.kernel.RSR_out
    sif_basis_hres = Float64.(ctx.sif_basis_hres)
    solar         = Float64.(collect(solar_hres))
    z_hres        = _normalized_grid(λ_hres)
    leg_basis     = _legendre_design_matrix(z_hres, n_legendre)

    vcd_o2_intercept = x[layout.idx_vcd_o2_intercept]
    vcd_o2_slope = x[layout.idx_vcd_o2_slope]
    vcd_h2o_intercept = x[layout.idx_vcd_h2o_intercept]
    vcd_h2o_slope = x[layout.idx_vcd_h2o_slope]
    p_o2_hpa = x[layout.idx_p_o2_hpa]
    t_o2_k = x[layout.idx_t_o2_k]
    p_h2o_hpa = x[layout.idx_p_h2o_hpa]
    t_h2o_k = x[layout.idx_t_h2o_k]
    vcd_o2_sif = x[layout.idx_vcd_o2_sif]
    vcd_h2o_sif = x[layout.idx_vcd_h2o_sif]
    sif_coeff = @view x[layout.idx_sif]
    leg_coeff = @view x[layout.idx_legendre]

    xs_o2  = vec(ctx.o2_sitp(spectral_axis, p_o2_hpa, t_o2_k))
    xs_h2o = vec(ctx.h2o_sitp(spectral_axis, p_h2o_hpa, t_h2o_k))

    vcd_o2_λ  = @. vcd_o2_intercept + vcd_o2_slope * z_hres
    vcd_h2o_λ = @. vcd_h2o_intercept + vcd_h2o_slope * z_hres
    trans     = @. exp(-(vcd_h2o_λ * xs_h2o + vcd_o2_λ * xs_o2))

    trans_sif = if model_variant == :standard
        @. exp(-(vcd_h2o_sif * xs_h2o + vcd_o2_sif * xs_o2))
    elseif model_variant == :vcd_ratio
        vcd_o2_sif_λ  = @. vcd_o2_sif * vcd_o2_λ
        vcd_h2o_sif_λ = @. vcd_h2o_sif * vcd_h2o_λ
        @. exp(-(vcd_h2o_sif_λ * xs_h2o + vcd_o2_sif_λ * xs_o2))
    else
        error("Unsupported model_variant: $model_variant")
    end

    rho_hres        = leg_basis * leg_coeff
    sif_hres_rad    = sif_basis_hres * sif_coeff
    y_refl_hres     = @. solar * trans * rho_hres / π
    y_sif_hres      = @. trans_sif * sif_hres_rad

    return (K * y_refl_hres, K * y_sif_hres)
end

"""LUT gas transmittance `trans` on high-res grid (Beer–Lambert); same as reflectance path factor."""
function xsec_gas_transmittance_hres(ctx, x::AbstractVector{<:Real}, layout)
    λ_hres        = collect(Float64.(ctx.λ_hres))
    spectral_axis = collect(Float64.(ctx.spectral_axis))
    z_hres        = _normalized_grid(λ_hres)

    vcd_o2_intercept = x[layout.idx_vcd_o2_intercept]
    vcd_o2_slope = x[layout.idx_vcd_o2_slope]
    vcd_h2o_intercept = x[layout.idx_vcd_h2o_intercept]
    vcd_h2o_slope = x[layout.idx_vcd_h2o_slope]
    vcd_o2_sif = x[layout.idx_vcd_o2_sif]
    vcd_h2o_sif = x[layout.idx_vcd_h2o_sif]
    p_o2_hpa = x[layout.idx_p_o2_hpa]
    t_o2_k = x[layout.idx_t_o2_k]
    p_h2o_hpa = x[layout.idx_p_h2o_hpa]
    t_h2o_k = x[layout.idx_t_h2o_k]

    xs_o2  = vec(ctx.o2_sitp(spectral_axis, p_o2_hpa, t_o2_k))
    xs_h2o = vec(ctx.h2o_sitp(spectral_axis, p_h2o_hpa, t_h2o_k))

    vcd_o2_λ  = @. vcd_o2_intercept + vcd_o2_slope * z_hres
    vcd_h2o_λ = @. vcd_h2o_intercept + vcd_h2o_slope * z_hres
    trans_updown = @. exp(-(vcd_h2o_λ * xs_h2o + vcd_o2_λ * xs_o2))
    trans_up     = @. exp(-(vcd_h2o_sif * xs_h2o + vcd_o2_sif * xs_o2))

    return (λ_hres, trans_updown, trans_up)
end

"""SVD `trans_up`, solar-path `trans_updown`, on `ctx.λ` (same PC basis as `make_svd_forward_model`)."""
function svd_transmittance_obs(
    ctx,
    PCs_obs::Matrix{<:Real},
    x::AbstractVector{<:Real},
    layout_svd;
    n_pc::Int,
    log_transform::Bool,
)
    λ_obs = collect(Float64.(ctx.λ))
    PCs   = Float64.(PCs_obs[:, 1:n_pc])
    c_vec = @view x[layout_svd.idx_pc]
    alpha_coeff = 10.0 ./ (1.0 .+ exp(-x[first(layout_svd.idx_alpha)])) + 1.0
    trans_up = log_transform ? exp.(PCs * c_vec) : 1.0 .+ PCs * c_vec
    trans_updown = exp.(alpha_coeff .* log.(max.(trans_up, eps(Float64))))
    return (λ_obs, trans_up, trans_updown)
end

"""Split SVD TOA into reflectance-path and SIF-path on the observation grid (matches `make_svd_forward_model`)."""
function svd_reflectance_and_sif_lres(
    ctx,
    solar_hres::AbstractVector{<:Real},
    PCs_obs::Matrix{<:Real},
    x::AbstractVector{<:Real},
    layout_svd;
    n_pc::Int,
    n_legendre::Int,
    log_transform::Bool,
)
    λ_obs          = collect(Float64.(ctx.λ))
    K              = hasproperty(ctx, :kernel_rsr_out) ? ctx.kernel_rsr_out : ctx.kernel.RSR_out
    sif_basis_hres = Float64.(ctx.sif_basis_hres)
    solar          = Float64.(collect(solar_hres))
    PCs            = Float64.(PCs_obs[:, 1:n_pc])

    solar_lres = Vector(K * solar)
    sif_lres = K * sif_basis_hres

    z_obs     = _normalized_grid(λ_obs)
    leg_basis = _legendre_design_matrix(z_obs, n_legendre)

    c_vec       = @view x[layout_svd.idx_pc]
    alpha_coeff = 10.0 ./ (1.0 .+ exp(-x[first(layout_svd.idx_alpha)])) + 1.0
    leg_coeff   = @view x[layout_svd.idx_legendre]
    sif_coeff   = @view x[layout_svd.idx_sif]

    trans_up     = log_transform ? exp.(PCs * c_vec) : 1.0 .+ PCs * c_vec
    trans_updown = exp.(alpha_coeff .* log.(max.(trans_up, eps(Float64))))
    rho_obs      = leg_basis * leg_coeff
    sif_lres_toa = trans_up .* (sif_lres * sif_coeff)

    y_refl_lres = @. solar_lres * trans_updown * rho_obs / π
    return (y_refl_lres, sif_lres_toa)
end

# ── Generic LM retrieval loop ─────────────────────────────────────────────────

"""
Run LM retrieval from `x0`, return a named tuple compatible with the benchmark format
used by `main()` in Fit_toy_forward_model.jl.
"""
function run_lm_retrieval(
    fm,
    jacobian_eval,
    x0::Vector{Float64},
    x_a::Vector{Float64},
    S_a_inv,
    x_scale::Vector{Float64},
    y_obs::Vector{Float64},
    lower_bounds::Vector{Float64},
    upper_bounds::Vector{Float64},
    wavelength::Vector{Float64};
    # LM params
    lm_lambda0::Float64    = 1.0,
    lm_lambda_up::Float64  = 5.0,
    lm_lambda_down::Float64 = 0.7,
    lm_lambda_min::Float64 = 1e-8,
    lm_lambda_max::Float64 = 1e8,
    lm_max_inner::Int      = 24,
    n_steps::Int           = 50,
    # Convergence
    conv_dx_rel_tol::Float64   = 1e-6,
    conv_rmse_rel_tol::Float64 = 1e-6,
    conv_rmse_abs_tol::Float64 = 1e-6,
    conv_stall_enable::Bool    = true,
    conv_stall_window::Int     = 3,
    conv_stall_redchi2_target::Float64  = 5.0,
    conv_stall_redchi2_abs_tol::Float64 = 0.1,
    conv_stall_redchi2_rel_tol::Float64 = 0.03,
    conv_stall_dx_rel_tol::Float64      = 5e-3,
    # SNR
    use_band_snr::Bool = true,
    band_snr_coeffs    = nothing,
    meas_sigma::Float64 = 0.01,
)
    n_forward  = Ref(0)
    n_jacobian = Ref(0)
    fm_eval    = x -> (n_forward[] += 1; fm(x))
    jac_eval   = x -> (n_jacobian[] += 1; jacobian_eval(x))

    x_curr    = copy(x0)
    y_curr    = fm_eval(x_curr)

    S_e_inv = if use_band_snr && !isnothing(band_snr_coeffs)
        make_Se_inv_from_snr(y_curr, band_snr_coeffs)
    else
        spdiagm(0 => fill(1.0 / meas_sigma^2, length(y_obs)))
    end

    dof = max(length(y_obs) - length(x0), 1)
    rmse0 = sqrt(mean((y_obs .- y_curr) .^ 2))

    rmse_series    = Float64[rmse0]
    dx_rel_series  = Float64[]
    redchi2_series = Float64[dot(y_obs .- y_curr, S_e_inv * (y_obs .- y_curr)) / dof]

    y_prior    = copy(y_curr)
    λ          = lm_lambda0
    converged  = false
    conv_reason = ""
    failed_step = 0
    t0         = time_ns()

    for istep in 1:n_steps
        x_prev    = copy(x_curr)
        rmse_prev = rmse_series[end]
        step = try
            lm_one_step(
                fm_eval, x_curr, y_obs;
                x_a            = x_a,
                S_a_inv        = S_a_inv,
                lambda         = λ,
                lambda_up      = lm_lambda_up,
                lambda_down    = lm_lambda_down,
                lambda_min     = lm_lambda_min,
                lambda_max     = lm_lambda_max,
                max_inner      = lm_max_inner,
                jacobian_eval  = jac_eval,
                x_scale        = x_scale,
                lower_bounds   = lower_bounds,
                upper_bounds   = upper_bounds,
                use_band_snr   = use_band_snr,
                band_snr_coeffs = band_snr_coeffs,
                meas_sigma     = meas_sigma,
            )
        catch err
            failed_step = istep
            conv_reason = "lm_one_step threw: " * sprint(showerror, err)
            break
        end

        λ = step.lambda_next

        if !step.accepted
            stalled, msg = _stalled_convergence(
                dx_rel_series, redchi2_series;
                enabled              = conv_stall_enable,
                window               = conv_stall_window,
                redchi2_target       = conv_stall_redchi2_target,
                redchi2_abs_tol      = conv_stall_redchi2_abs_tol,
                redchi2_rel_tol      = conv_stall_redchi2_rel_tol,
                dx_rel_tol           = conv_stall_dx_rel_tol,
            )
            if stalled
                converged   = true
                conv_reason = msg
            else
                failed_step = istep
                conv_reason = "LM outer step rejected (no acceptable inner iterate); " * msg
            end
            break
        end

        x_curr  = step.x_next
        y_curr  = copy(step.y_next)
        rmse_new = sqrt(mean((y_obs .- y_curr) .^ 2))
        push!(rmse_series, rmse_new)

        chi2_new = dot(y_obs .- y_curr, S_e_inv * (y_obs .- y_curr))
        push!(redchi2_series, chi2_new / dof)

        dx_rel = norm(step.dx) / max(norm(x_prev), eps(Float64))
        push!(dx_rel_series, dx_rel)

        rmse_abs = abs(rmse_new - rmse_prev)
        rmse_rel = rmse_abs / max(abs(rmse_prev), eps(Float64))
        if dx_rel < conv_dx_rel_tol || rmse_rel < conv_rmse_rel_tol || rmse_abs < conv_rmse_abs_tol
            converged   = true
            conv_reason = "dx_rel=$(dx_rel), rmse_abs=$(rmse_abs), rmse_rel=$(rmse_rel)"
            break
        end
    end

    elapsed_s = (time_ns() - t0) / 1e9
    return (
        rmse_series       = rmse_series,
        x_final           = copy(x_curr),
        x_prior           = copy(x0),
        elapsed_wall_s    = elapsed_s,
        n_forward         = n_forward[],
        n_jacobian        = n_jacobian[],
        converged         = converged,
        convergence_reason = conv_reason,
        failed_step       = failed_step,
        n_steps_done      = length(rmse_series) - 1,
        wavelength        = wavelength,
        y_obs             = copy(y_obs),
        y_prior           = y_prior,
        y_final           = copy(y_curr),
    )
end

"""
Posterior covariance approximation at the final state:
  S_post = (K' S_e^{-1} K + S_a^{-1})^{-1}
Returns posterior standard deviations σ_post = sqrt(diag(S_post)).
"""
function posterior_sigma_from_jacobian(K::AbstractMatrix, S_e_inv, S_a_inv)
    H = K' * S_e_inv * K + S_a_inv
    S_post = Matrix(inv(H))
    return sqrt.(max.(diag(S_post), 0.0))
end

"""Column-wise Jacobian normalization: J_norm[:,j] = J[:,j] * σ_post[j]."""
function normalize_jacobian_by_sigma(J::AbstractMatrix, σ_post::AbstractVector)
    size(J, 2) == length(σ_post) ||
        error("σ_post length $(length(σ_post)) must match Jacobian columns $(size(J,2))")
    return J .* permutedims(collect(Float64.(σ_post)))
end

"""
Gain matrix and averaging kernel at final state:
  G = (S_a^{-1} + K' S_e^{-1} K)^{-1} K' S_e^{-1}
  A = G K
"""
function gain_and_averaging_kernel(K::AbstractMatrix, S_e_inv, S_a_inv)
    H = S_a_inv + K' * S_e_inv * K
    G = Matrix(H \ (K' * S_e_inv))
    A = Matrix(G * K)
    return G, A
end

# ── Load transmittance NetCDF and build SVD basis ──────────────────────────────

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
        collect(Float64.(ds["transmittance"][:, :])),  # (n_profiles, n_bands)
        collect(Float64.(ds["band"][:]))
    end
    trans_w = NCDataset(winter_nc, "r") do ds
        collect(Float64.(ds["transmittance"][:, :]))
    end
    trans = vcat(trans_s, trans_w)   # (n_profiles_total, n_bands)
    n_profiles = size(trans, 1)

    # Select fitting window
    ind = findall(λ_min .< bands .< λ_max)
    isempty(ind) && error("No transmittance bands found in [$(λ_min), $(λ_max)] nm")
    bands_sel = bands[ind]
    trans_sel = trans[:, ind]         # (n_profiles, n_bands_sel)

    # Center: subtract 1 (or log transform)
    a = log_transform ? log.(max.(trans_sel, 1e-10)) : (trans_sel .- 1.0)

    # SVD: a' is (n_bands_sel, n_profiles)
    F  = svd(a')
    U  = F.U                # (n_bands_sel, n_svs)
    S  = F.S                # (n_svs,) singular values
    S_norm = S ./ sum(S) .* 100  # percent variance

    println("SVD basis: $(n_profiles) profiles, $(length(bands_sel)) bands in [$(λ_min), $(λ_max)] nm")
    println("  Variance explained by first $n_pc PCs: ",
            @sprintf("%.2f", sum(S_norm[1:min(n_pc, length(S_norm))])), "%")
    for k in 1:min(n_pc, length(S_norm))
        @printf("  PC%-2d: %.3f%%  (S_raw=%.4e, σ_prior=%.4e)\n",
                k, S_norm[k], S[k], S[k] / sqrt(n_profiles))
    end

    # Interpolate each PC onto the target grid (OCI band centers = ctx.λ for the forward model)
    nλ   = length(λ_pc_target)
    PCs  = zeros(Float64, nλ, size(U, 2))
    for k in eachindex(axes(U, 2))
        itp = LinearInterpolation(bands_sel, U[:, k]; extrapolation_bc=Flat())
        PCs[:, k] .= itp.(λ_pc_target)
    end

    return (
        PCs        = PCs,          # (nλ, n_svs)
        S          = S,            # singular values (unnormalised)
        S_norm     = S_norm,       # % variance
        n_profiles = n_profiles,
        bands_sel  = bands_sel,
        U          = U,            # (n_bands_sel, n_svs)
    )
end

# ── Plotting helpers ──────────────────────────────────────────────────────────

function _plot_svd_basis(path::String, svd_basis, n_pc::Int)
    λ  = svd_basis.bands_sel
    U  = svd_basis.U
    SN = svd_basis.S_norm
    colors = [:royalblue, :firebrick, :darkgreen, :darkorange, :purple,
              :teal, :brown, :magenta, :gray, :olive]

    p = plot(
        xlabel = "Wavelength [nm]",
        ylabel = "PC amplitude",
        title  = "SVD transmittance basis (first $n_pc PCs)",
        legend = :outertopright,
        size   = (900, 450),
        legendfontsize = 7,
    )
    for k in 1:min(n_pc, size(U, 2))
        plot!(p, λ, U[:, k];
              label  = @sprintf("PC%d (%.2f%%)", k, SN[k]),
              color  = colors[mod1(k, length(colors))],
              lw     = 1.8)
    end
    savefig(p, path)
end

function _plot_rmse(path::String, res_xsec, res_svd)
    p = plot(
        xlabel = "Iteration (0 = prior)",
        ylabel = "RMSE (spectral)",
        title  = "RMSE vs iteration: LUT xSec vs SVD transmittance fit",
        legend = :topright,
        size   = (800, 450),
    )
    plot!(p, 0:(length(res_xsec.rmse_series) - 1), res_xsec.rmse_series;
          label = "LUT xSec (LM)",  color = :royalblue,  lw = 2.0,
          marker = :circle,  markersize = 3)
    plot!(p, 0:(length(res_svd.rmse_series) - 1), res_svd.rmse_series;
          label = "SVD trans. (LM)", color = :firebrick, lw = 2.0,
          marker = :square, markersize = 3)
    savefig(p, path)
end

function _plot_spectra(path::String, res_xsec, res_svd)
    λ     = res_xsec.wavelength
    y_obs = res_xsec.y_obs
    colors = [:royalblue, :firebrick]

    p_spec = plot(
        λ, y_obs;
        label = "Observed", color = :black, lw = 2.0,
        xlabel = "Wavelength [nm]", ylabel = "R_toa",
        title  = "Reconstructed vs observed R_toa",
        legend = :outertopright, legendfontsize = 7,
    )
    plot!(p_spec, λ, res_xsec.y_prior;  label = "Prior (xSec)", color = :gray, lw = 1.2, ls = :dash)
    plot!(p_spec, λ, res_xsec.y_final;  label = "LUT xSec final",   color = colors[1], lw = 1.8)
    plot!(p_spec, λ, res_svd.y_final;   label = "SVD trans. final",  color = colors[2], lw = 1.8)

    p_resid = plot(
        λ, zeros(length(λ));
        label = nothing, color = :black, lw = 1.0, ls = :dot,
        xlabel = "Wavelength [nm]", ylabel = "Residual (obs − fit)",
        title  = "Spectral residuals at final state",
        legend = :outertopright, legendfontsize = 7,
    )
    for (label, res, col) in [
            ("LUT xSec",    res_xsec, colors[1]),
            ("SVD trans.",  res_svd,  colors[2]),
        ]
        plot!(p_resid, λ, y_obs .- res.y_final; label = label, color = col, lw = 1.5)
    end

    p = plot(p_spec, p_resid; layout = (2, 1), size = (950, 850), left_margin = 5Plots.mm)
    savefig(p, path)
end

function _plot_sif_trans_components(
    path::String,
    λ_obs::AbstractVector{<:Real},
    refl_xsec::AbstractVector,
    sif_xsec::AbstractVector,
    refl_svd::AbstractVector,
    sif_svd::AbstractVector,
)
    colors = [:royalblue, :firebrick]

    p_refl = plot(
        xlabel = "Wavelength [nm]",
        ylabel = "Radiance contribution",
        title  = "Convolved solar×T×ρ/π path (reflectance term, excludes SIF)",
        legend = :outertopright,
        legendfontsize = 7,
    )
    plot!(p_refl, λ_obs, refl_xsec; label = "LUT xSec", color = colors[1], lw = 2.0)
    plot!(p_refl, λ_obs, refl_svd;  label = "SVD trans.", color = colors[2], lw = 2.0)

    p_sif = plot(
        xlabel = "Wavelength [nm]",
        ylabel = "Radiance contribution",
        title  = "Convolved SIF-path TOA (K × SIF term)",
        legend = :outertopright,
        legendfontsize = 7,
    )
    plot!(p_sif, λ_obs, sif_xsec; label = "LUT xSec", color = colors[1], lw = 2.0)
    plot!(p_sif, λ_obs, sif_svd;  label = "SVD trans.", color = colors[2], lw = 2.0)

    p = plot(p_refl, p_sif; layout = (2, 1), size = (950, 850), left_margin = 5Plots.mm)
    savefig(p, path)
end

"""LUT: one subplot per state — ∂y/∂x_j vs λ (hybrid vs ForwardDiff)."""
function _plot_jacobian_lut_hybrid_vs_fd(
    path::String,
    λ_obs::AbstractVector{<:Real},
    J_hybrid::AbstractMatrix{<:Real},
    J_fd::AbstractMatrix{<:Real},
    state_names::AbstractVector{String};
    n_cols::Int = 1,
)
    size(J_hybrid) == size(J_fd) ||
        error("Jacobian shapes differ: hybrid $(size(J_hybrid)) vs FD $(size(J_fd))")
    nb, ns = size(J_hybrid)
    length(state_names) == ns ||
        error("state_names length $(length(state_names)) ≠ Jacobian columns $ns")
    length(λ_obs) == nb ||
        error("λ_obs length $(length(λ_obs)) ≠ Jacobian rows (bands) $nb")

    _plot_jacobian_columns_spectra(
        path, λ_obs, state_names, J_hybrid;
        J_compare = J_fd,
        supertitle = "LUT xSec final state — normalized Jacobian (∂y/∂x)·σ_post vs λ",
        label_a = "hybrid norm.",
        label_b = "ForwardDiff norm.",
        n_cols = n_cols,
    )
end

"""SVD: one subplot per state — ∂y/∂x_j vs λ (ForwardDiff)."""
function _plot_jacobian_svd_spectra(
    path::String,
    λ_obs::AbstractVector{<:Real},
    J_svd::AbstractMatrix{<:Real},
    state_names::AbstractVector{String};
    n_cols::Int = 1,
)
    nb, ns = size(J_svd)
    length(state_names) == ns ||
        error("state_names length $(length(state_names)) ≠ Jacobian columns $ns")
    length(λ_obs) == nb ||
        error("λ_obs length $(length(λ_obs)) ≠ Jacobian rows (bands) $nb")

    _plot_jacobian_columns_spectra(
        path, λ_obs, state_names, J_svd;
        J_compare = nothing,
        supertitle = "SVD transmittance final state — normalized Jacobian (∂y/∂x)·σ_post vs λ",
        label_a = "ForwardDiff norm.",
        label_b = "",
        n_cols = n_cols,
    )
end

function _plot_jacobian_columns_spectra(
    path::String,
    λ_obs::AbstractVector{<:Real},
    state_names::AbstractVector{String},
    J_primary::AbstractMatrix{<:Real};
    J_compare::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    supertitle::String,
    label_a::String,
    label_b::String,
    n_cols::Int,
)
    _, ns = size(J_primary)
    n_cols = max(1, min(n_cols, ns))
    n_rows = cld(ns, n_cols)

    plots = []

    dual = !isnothing(J_compare)
    if dual
        size(J_compare) == size(J_primary) || error("Compare Jacobian shape must match primary")
    end

    for j in 1:ns
        col_a = @view J_primary[:, j]

        ylab = "(∂y/∂$(state_names[j]))·σ_post"
        pj = plot(;
            title = state_names[j],
            titlefontsize = 8,
            ylabel = ylab,
            yguidefontsize = 6,
            grid = true,
            legend = j == 1 ? :topright : nothing,
            legendfontsize = 5,
            xlabel = "Wavelength [nm]",
            xlabelfontsize = 7,
            left_margin = 10Plots.mm,
            bottom_margin = 3Plots.mm,
        )

        if dual
            m = maximum(abs.(vcat(Vector(col_a), Vector(J_compare[:, j]))))
            plot!(pj, λ_obs, col_a; label = label_a, color = :royalblue, lw = 1.2)
            plot!(pj, λ_obs, J_compare[:, j]; label = label_b, color = :firebrick, lw = 1.2, ls = :dash)
            plot!(pj; ylims = (-max(m, 1e-30), max(m, 1e-30)))
        else
            ca = Vector(col_a)
            m = maximum(abs.(ca))
            plot!(pj, λ_obs, ca; label = label_a, color = :royalblue, lw = 1.2)
            plot!(pj; ylims = (-max(m, 1e-30), max(m, 1e-30)))
        end
        push!(plots, pj)
    end

    pw = 360 * n_cols + 80
    ph = 175 * n_rows + 75
    p = plot(
        plots...;
        layout       = (n_rows, n_cols),
        size         = (pw, ph),
        plot_title   = supertitle,
        plot_titlefontsize = 11,
        top_margin   = 6Plots.mm,
    )
    savefig(p, path)
end

function _plot_averaging_kernel_heatmap(
    path::String,
    A::AbstractMatrix{<:Real},
    state_names::AbstractVector{String};
    title::String,
)
    n1, n2 = size(A)
    n1 == n2 || error("Averaging kernel must be square, got size $(size(A))")
    length(state_names) == n1 ||
        error("state_names length $(length(state_names)) must match A size $n1")

    h = heatmap(
        1:n2, 1:n1, A;
        xlabel = "True state index (x_true)",
        ylabel = "Retrieved state index (x_ret)",
        xticks = (1:n2, state_names),
        yticks = (1:n1, state_names),
        xrotation = 60,
        tickfontsize = 7,
        color = :balance,
        clims = (-1.0, 1.0),
        colorbar = true,
        left_margin = 12Plots.mm,
        bottom_margin = 10Plots.mm,
        title = title,
    )
    savefig(h, path)
end

function _plot_transmittance_one_two(
    path::String,
    λ_obs::AbstractVector,
    T_up_xsec_lres::AbstractVector,
    T_updown_xsec_lres::AbstractVector,
    T_up_svd::AbstractVector,
    T_updown_svd::AbstractVector,
)
    colors = [:royalblue, :dodgerblue, :firebrick, :orangered, :gray]

    p_lut = plot(
        xlabel = "Wavelength [nm]",
        ylabel = "Transmittance (convolved)",
        title  = "LUT xSec: K×T vs K×T2 on OCI bands (gas, final state)",
        legend = :outertopright,
        legendfontsize = 7,
    )
    plot!(p_lut, λ_obs, T_up_xsec_lres; label = "One-way lres T (SIF)", color = colors[1], lw = 2.0)
    plot!(p_lut, λ_obs, T_updown_xsec_lres; label = "Two-way lres T (solar)", color = colors[2], lw = 2.0, ls = :dash)

    p_svd = plot(
        xlabel = "Wavelength [nm]",
        ylabel = "Transmittance (PC basis on OCI λ)",
        title  = "SVD: T_up, T_up^α on observation bands (final state)",
        legend = :outertopright,
        legendfontsize = 7,
    )
    plot!(p_svd, λ_obs, T_up_svd; label = "One-way T_up (PC basis)", color = colors[3], lw = 2.0)
    plot!(p_svd, λ_obs, T_updown_svd; label = "Solar-path T_up^α", color = colors[5], lw = 1.8, ls = :dot)

    p = plot(p_lut, p_svd; layout = (2, 1), size = (950, 850), left_margin = 5Plots.mm)
    savefig(p, path)
end

function _write_summary(path::String, res_xsec, res_svd)
    open(path, "w") do io
        println(io, "method,n_iter,n_forward,n_jacobian,rmse_prior,rmse_final," *
                    "rmse_reduction_pct,elapsed_s,converged,convergence_reason")
        for (label, res) in [("LUT_xSec_LM", res_xsec), ("SVD_trans_LM", res_svd)]
            rmse0 = res.rmse_series[1]
            rmsef = res.rmse_series[end]
            red   = 100 * (rmse0 - rmsef) / max(abs(rmse0), eps(Float64))
            println(io, join([
                label, res.n_steps_done,
                res.n_forward, res.n_jacobian,
                rmse0, rmsef, red,
                res.elapsed_wall_s,
                res.converged,
                repr(res.convergence_reason),
            ], ","))
        end
    end
end

"""Ordered parameter names matching `svd_state_layout` (index `i` ↔ name `i`)."""
function _svd_state_names(layout)
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
    length(names) == layout.n_state ||
        error("internal: SVD name count $(length(names)) ≠ n_state $(layout.n_state)")
    return names
end

function _print_retrieval_state(
    io::IO,
    title::AbstractString,
    res,
    names::AbstractVector{String},
)
    xp = res.x_prior
    xf = res.x_final
    length(names) == length(xf) ||
        error("state name count $(length(names)) ≠ length(x_final) $(length(xf))")
    dx = xf .- xp
    nrm_dx = norm(dx)
    nrm_xp = norm(xp)
    rel = nrm_dx / max(nrm_xp, eps(Float64))
    println(io)
    println(io, repeat("-", 70))
    println(io, title)
    println(io, "  converged:           ", res.converged)
    println(io, "  convergence_reason:  ", isempty(res.convergence_reason) ? "(none)" : res.convergence_reason)
    println(io, "  failed_step:         ", res.failed_step)
    println(io, "  n_outer_steps:       ", res.n_steps_done)
    @printf(io, "  ||Δx||:              %.6e\n", nrm_dx)
    @printf(io, "  ||Δx|| / ||x_prior||: %.6e\n", rel)
    println(io, "  x_prior → x_final  (Δ = final − prior):")
    w = max(8, maximum(length.(names)))
    for i in eachindex(xf)
        @printf(
            io,
            "    %-*s  % .8e  →  % .8e   (Δ % .4e)\n",
            w, names[i], xp[i], xf[i], dx[i],
        )
    end
end

# ── Main ──────────────────────────────────────────────────────────────────────

function main_compare()
    config_path = get(
        ENV, "SVD_XSEC_CONFIG",
        joinpath(_SANITY_DIR, "svd_vs_xsec_config.toml"),
    )
    isfile(config_path) || error("Config not found: $config_path")
    cfg = TOML.parsefile(config_path)
    mkpath(_OUT_DIR)

    println("="^70)
    println("SVD Transmittance vs LUT Cross-Section Fit Comparison")
    println("Config: $config_path")
    println("="^70)

    # ── Shared setup ────────────────────────────────────────────────────────
    ENV["PACE_MWE_CONFIG"] = config_path
    ctx = prepare_mwe_inputs(config_path)

    data_cfg  = get(cfg, "data",            Dict{String,Any}())
    pace_cfg  = get(cfg, "pace_observation", Dict{String,Any}())
    spec_cfg  = get(cfg, "spectral",        Dict{String,Any}())
    fit_cfg   = get(cfg, "fit",             Dict{String,Any}())
    xsec_cfg  = get(fit_cfg, "xsec",        Dict{String,Any}())
    svd_cfg   = get(fit_cfg, "svd",         Dict{String,Any}())

    float_type     = parse_float_type(cfg)
    meas_sigma     = Float64(get(fit_cfg, "meas_sigma",         0.01))
    sif_nev        = Int(get(spec_cfg, "sif_nev",               1))
    n_steps        = Int(get(fit_cfg, "n_plot_steps",           50))
    kernel_cfg     = get(cfg, "kernel", Dict{String,Any}())
    kernel_wants_snr = Bool(get(kernel_cfg, "use_band_snr", true))
    # Retrieval noise: [fit] use_band_snr overrides; if absent, follow [kernel] (same as prepare_mwe_inputs).
    use_band_snr     = haskey(fit_cfg, "use_band_snr") ? Bool(fit_cfg["use_band_snr"]) : kernel_wants_snr
    if use_band_snr && isnothing(ctx.band_snr_coeffs)
        error(
            "Band SNR requested (use_band_snr=true) but SNR coefficients were not loaded. " *
            "Set [kernel] use_band_snr = true and a valid [data] pace_snr_file (see Simple_PACE_xSecFit_MWE_zcheVer.toml), " *
            "then re-run so prepare_mwe_inputs can populate ctx.band_snr_coeffs.",
        )
    end
    if use_band_snr && !kernel_wants_snr
        error(
            "[fit] use_band_snr=true requires [kernel] use_band_snr=true so the SNR file is read during prepare_mwe_inputs; " *
            "or set [fit] use_band_snr=false to use meas_sigma only.",
        )
    end
    if use_band_snr
        println(
            "Retrieval measurement noise: band-specific SNR (σ² = c1 + c2·y) from ",
            basename(string(ctx.paths.pace_snr_path)),
        )
    else
        println("Retrieval measurement noise: uniform meas_sigma = ", meas_sigma)
    end
    prior_sigma_default = Float64(get(fit_cfg, "prior_sigma_default", 1e12))
    prior_min_sigma     = Float64(get(fit_cfg, "prior_min_sigma",     1e-3))
    sif_sigma           = Float64(get(fit_cfg, "sif_sigma",           1e12))

    # Shared LM controls
    lm_lambda0    = Float64(get(fit_cfg, "lm_lambda0",    1.0))
    lm_lambda_up  = Float64(get(fit_cfg, "lm_lambda_up",  5.0))
    lm_lambda_down = Float64(get(fit_cfg, "lm_lambda_down", 0.7))
    lm_lambda_min = Float64(get(fit_cfg, "lm_lambda_min", 1e-8))
    lm_lambda_max = Float64(get(fit_cfg, "lm_lambda_max", 1e8))
    lm_max_inner  = Int(get(fit_cfg, "lm_max_inner",     24))

    conv_dx_rel_tol   = Float64(get(fit_cfg, "conv_dx_rel_tol",   1e-6))
    conv_rmse_rel_tol = Float64(get(fit_cfg, "conv_rmse_rel_tol", 1e-6))
    conv_rmse_abs_tol = Float64(get(fit_cfg, "conv_rmse_abs_tol", 1e-6))
    conv_stall_enable          = Bool(get(fit_cfg,   "conv_stall_enable",          true))
    conv_stall_window          = Int(get(fit_cfg,    "conv_stall_window",           3))
    conv_stall_redchi2_target  = Float64(get(fit_cfg,"conv_stall_redchi2_target",  5.0))
    conv_stall_redchi2_abs_tol = Float64(get(fit_cfg,"conv_stall_redchi2_abs_tol", 0.1))
    conv_stall_redchi2_rel_tol = Float64(get(fit_cfg,"conv_stall_redchi2_rel_tol", 0.03))
    conv_stall_dx_rel_tol      = Float64(get(fit_cfg,"conv_stall_dx_rel_tol",      5e-3))

    lm_kwargs = (
        lm_lambda0     = lm_lambda0,   lm_lambda_up = lm_lambda_up,
        lm_lambda_down = lm_lambda_down, lm_lambda_min = lm_lambda_min,
        lm_lambda_max  = lm_lambda_max,  lm_max_inner = lm_max_inner,
        n_steps        = n_steps,
        conv_dx_rel_tol   = conv_dx_rel_tol, conv_rmse_rel_tol = conv_rmse_rel_tol,
        conv_rmse_abs_tol = conv_rmse_abs_tol,
        conv_stall_enable = conv_stall_enable,
        conv_stall_window = conv_stall_window,
        conv_stall_redchi2_target  = conv_stall_redchi2_target,
        conv_stall_redchi2_abs_tol = conv_stall_redchi2_abs_tol,
        conv_stall_redchi2_rel_tol = conv_stall_redchi2_rel_tol,
        conv_stall_dx_rel_tol      = conv_stall_dx_rel_tol,
        use_band_snr    = use_band_snr,
        band_snr_coeffs = ctx.band_snr_coeffs,
        meas_sigma      = meas_sigma,
    )

    # Solar spectrum
    solar_file = get(data_cfg, "solar_file", "solar_merged_20240731_600_33300_100.out")
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)
    solar_hres, _ = load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
    solar_hres = float_type.(solar_hres)

    # PACE observation
    pace_file = get(pace_cfg, "pace_file", "sample_granule_20240830T131442_new_chl.nc")
    pace_path = isabspath(pace_file) ? pace_file : joinpath(ctx.paths.base_dir, pace_file)
    pixel_idx  = Int(get(pace_cfg, "pixel_index", 600))
    scan_idx   = Int(get(pace_cfg, "scan_index",  800))
    y_obs, _ = load_pace_spectrum_on_grid(
        pace_path, ctx.λ;
        pixel_idx = pixel_idx, scan_idx = scan_idx,
        wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength")),
        spectrum_var   = String(get(pace_cfg, "spectrum_var",   "radiance_red")),
    )
    y_obs    = float_type.(y_obs)
    λ_hres   = collect(Float64.(ctx.λ_hres))
    λ_obs    = collect(Float64.(ctx.λ))

    # ── LUT xSec fit ────────────────────────────────────────────────────────
    println("\n[1/2] LUT cross-section fit ...")

    n_leg_xsec = Int(get(xsec_cfg, "n_legendre",    3))
    model_variant = Symbol(get(xsec_cfg, "model_variant", "standard"))
    use_hybrid_jac = Bool(get(xsec_cfg, "use_hybrid_jacobian", true))

    fm_xsec = make_forward_model_simple(
        ctx, solar_hres;
        n_legendre       = n_leg_xsec,
        preallocate_float64 = true,
        model_variant    = model_variant,
    )
    layout_xsec = state_layout_simple(ctx; n_legendre = n_leg_xsec)
    x0_xsec     = initial_state_simple(ctx; n_legendre = n_leg_xsec, T = Float64,
                                       model_variant = model_variant)
    jac_xsec = use_hybrid_jac ?
        make_hybrid_jacobian_evaluator(fm_xsec, ctx, solar_hres, layout_xsec;
                                       n_legendre = n_leg_xsec) :
        make_jacobian_evaluator(fm_xsec, x0_xsec)

    # --- xSec priors ---
    vcd_o2_sigma  = Float64(get(xsec_cfg, "vcd_o2_sigma",  1e25))
    vcd_h2o_sigma = Float64(get(xsec_cfg, "vcd_h2o_sigma", 3e22))
    p_prior_hpa   = Float64(get(xsec_cfg, "p_prior_hpa",   700.0))
    p_sigma_hpa   = Float64(get(xsec_cfg, "p_sigma_hpa",   200.0))
    t_prior_k     = Float64(get(xsec_cfg, "t_prior_k",     273.0))
    t_sigma_k     = Float64(get(xsec_cfg, "t_sigma_k",     20.0))
    use_pt_constraints      = Bool(get(xsec_cfg, "use_pt_constraints",     true))
    pt_constraint_sigma_mult = Float64(get(xsec_cfg, "pt_constraint_sigma_mult", 3.0))
    vcd_slope_prior_sigma_factor = Float64(get(xsec_cfg, "vcd_slope_prior_sigma_factor", 1.0))
    leg01_frac_xsec = Float64(get(xsec_cfg, "legendre01_prior_sigma_fraction", 1.0))
    leg_higher_xsec = Float64(get(xsec_cfg, "legendre_higher_sigma", 1.0))
    use_leg01_xsec  = Bool(get(xsec_cfg, "use_legendre01_prior", true))
    use_leghig_xsec = Bool(get(xsec_cfg, "use_legendre_higher_prior", true))

    x_a_xsec     = copy(x0_xsec)
    prior_sigma_xsec = fill(prior_sigma_default, length(x0_xsec))

    # Legendre P0/P1 prior: fit ratio y_obs / y0 with a linear polynomial
    if use_leg01_xsec && length(layout_xsec.idx_legendre) >= 1
        y0 = fm_xsec(x0_xsec)
        ratio = y_obs ./ max.(abs.(y0), eps(Float64))
        z = _normalized_grid(λ_obs)
        A01 = hcat(ones(length(z)), z)
        w = y_obs .- minimum(y_obs); w .+= max(maximum(w), 1.0) * 1e-6
        s = sqrt.(w ./ maximum(w))
        c01 = (A01 .* s) \ (ratio .* s)
        leg0 = first(layout_xsec.idx_legendre)
        x_a_xsec[leg0] = c01[1]
        prior_sigma_xsec[leg0] = max(abs(c01[1]) * leg01_frac_xsec, prior_min_sigma)
        if length(layout_xsec.idx_legendre) >= 2
            leg1 = layout_xsec.idx_legendre[2]
            x_a_xsec[leg1] = c01[2]
            prior_sigma_xsec[leg1] = max(abs(c01[2]) * leg01_frac_xsec, prior_min_sigma)
        end
    end
    if use_leghig_xsec && length(layout_xsec.idx_legendre) >= 3
        for j in 3:length(layout_xsec.idx_legendre)
            prior_sigma_xsec[layout_xsec.idx_legendre[j]] = max(leg_higher_xsec, prior_min_sigma)
        end
    end

    x_a_xsec[layout_xsec.idx_p_o2_hpa]  = p_prior_hpa
    x_a_xsec[layout_xsec.idx_p_h2o_hpa] = p_prior_hpa
    x_a_xsec[layout_xsec.idx_t_o2_k]    = t_prior_k
    x_a_xsec[layout_xsec.idx_t_h2o_k]   = t_prior_k
    prior_sigma_xsec[layout_xsec.idx_p_o2_hpa]  = p_sigma_hpa
    prior_sigma_xsec[layout_xsec.idx_p_h2o_hpa] = p_sigma_hpa
    prior_sigma_xsec[layout_xsec.idx_t_o2_k]    = t_sigma_k
    prior_sigma_xsec[layout_xsec.idx_t_h2o_k]   = t_sigma_k
    prior_sigma_xsec[layout_xsec.idx_vcd_o2_intercept]  = vcd_o2_sigma
    prior_sigma_xsec[layout_xsec.idx_vcd_h2o_intercept] = vcd_h2o_sigma
    prior_sigma_xsec[layout_xsec.idx_vcd_o2_sif]  = vcd_o2_sigma
    prior_sigma_xsec[layout_xsec.idx_vcd_h2o_sif] = vcd_h2o_sigma
    prior_sigma_xsec[layout_xsec.idx_vcd_o2_slope]  = max(vcd_o2_sigma  * vcd_slope_prior_sigma_factor, prior_min_sigma)
    prior_sigma_xsec[layout_xsec.idx_vcd_h2o_slope] = max(vcd_h2o_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
    prior_sigma_xsec[layout_xsec.idx_sif] .= max(sif_sigma, prior_min_sigma)
    S_a_inv_xsec = _spdiag_invvar(prior_sigma_xsec)

    x_scale_xsec = ones(Float64, length(x0_xsec))
    x_scale_xsec[layout_xsec.idx_vcd_o2_intercept] = vcd_o2_sigma
    x_scale_xsec[layout_xsec.idx_vcd_o2_slope]     = max(vcd_o2_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
    x_scale_xsec[layout_xsec.idx_vcd_h2o_intercept] = vcd_h2o_sigma
    x_scale_xsec[layout_xsec.idx_vcd_h2o_slope]     = max(vcd_h2o_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
    x_scale_xsec[layout_xsec.idx_vcd_o2_sif]  = vcd_o2_sigma
    x_scale_xsec[layout_xsec.idx_vcd_h2o_sif] = vcd_h2o_sigma
    x_scale_xsec[layout_xsec.idx_p_o2_hpa]    = p_sigma_hpa
    x_scale_xsec[layout_xsec.idx_p_h2o_hpa]   = p_sigma_hpa
    x_scale_xsec[layout_xsec.idx_t_o2_k]      = t_sigma_k
    x_scale_xsec[layout_xsec.idx_t_h2o_k]     = t_sigma_k

    lower_xsec = fill(-Inf, length(x0_xsec))
    upper_xsec = fill( Inf, length(x0_xsec))
    if use_pt_constraints
        for idx in [layout_xsec.idx_p_o2_hpa, layout_xsec.idx_p_h2o_hpa]
            lower_xsec[idx] = p_prior_hpa - pt_constraint_sigma_mult * p_sigma_hpa
            upper_xsec[idx] = p_prior_hpa + pt_constraint_sigma_mult * p_sigma_hpa
        end
        for idx in [layout_xsec.idx_t_o2_k, layout_xsec.idx_t_h2o_k]
            lower_xsec[idx] = t_prior_k - pt_constraint_sigma_mult * t_sigma_k
            upper_xsec[idx] = t_prior_k + pt_constraint_sigma_mult * t_sigma_k
        end
    end

    res_xsec = run_lm_retrieval(
        fm_xsec, jac_xsec, Float64.(x_a_xsec), Float64.(x_a_xsec),
        S_a_inv_xsec, x_scale_xsec,
        Float64.(y_obs), lower_xsec, upper_xsec, λ_obs;
        lm_kwargs...,
    )
    @printf("  xSec: %d iter, %d fwd, %d jac, RMSE %.4e → %.4e (%.1f%% reduction)\n",
            res_xsec.n_steps_done, res_xsec.n_forward, res_xsec.n_jacobian,
            res_xsec.rmse_series[1], res_xsec.rmse_series[end],
            100*(res_xsec.rmse_series[1]-res_xsec.rmse_series[end]) /
                max(abs(res_xsec.rmse_series[1]), eps(Float64)))
    _print_retrieval_state(
        stdout, "LUT cross-section (LM) — state vector",
        res_xsec, state_names_simple(ctx; n_legendre = n_leg_xsec),
    )

    # ── SVD transmittance fit ────────────────────────────────────────────────
    println("\n[2/2] SVD transmittance fit ...")

    n_leg_svd   = Int(get(svd_cfg, "n_legendre",        3))
    n_pc        = Int(get(svd_cfg, "n_pc",              5))
    log_trans   = Bool(get(svd_cfg, "svd_log_transform", false))
    pc_prior_mode  = String(get(svd_cfg, "pc_prior_mode",  "loading_variance"))
    pc_sigma_scale = Float64(get(svd_cfg, "pc_prior_sigma_scale", 1.0))
    λ_min = Float64(get(cfg["spectral"], "lambda_min_nm", 640.0))
    λ_max = Float64(get(cfg["spectral"], "lambda_max_nm", 756.0))
    summer_nc = svd_cfg["summer_nc"]
    winter_nc = svd_cfg["winter_nc"]

    svd_basis = load_svd_basis(
        summer_nc, winter_nc, λ_obs;
        λ_min = λ_min, λ_max = λ_max, n_pc = n_pc, log_transform = log_trans,
    )

    fm_svd, layout_svd = make_svd_forward_model(
        ctx, solar_hres, svd_basis.PCs;
        n_pc = n_pc, n_legendre = n_leg_svd, log_transform = log_trans,
    )
    n_ev_svd = layout_svd.n_ev
    jac_svd  = make_jacobian_evaluator(fm_svd, zeros(layout_svd.n_state))

    # --- SVD priors ---
    leg01_frac_svd = Float64(get(svd_cfg, "legendre01_prior_sigma_fraction", 1.0))
    leg_higher_svd = Float64(get(svd_cfg, "legendre_higher_sigma", 1.0))
    use_leg01_svd  = Bool(get(svd_cfg, "use_legendre01_prior", true))
    use_leghig_svd = Bool(get(svd_cfg, "use_legendre_higher_prior", true))

    alpha_mean   = Float64(get(svd_cfg, "alpha_prior_mean",  1.0))
    alpha_sigma  = Float64(get(svd_cfg, "alpha_prior_sigma", 0.3))
    alpha_sigma  = max(alpha_sigma, prior_min_sigma)

    x0_svd       = zeros(Float64, layout_svd.n_state)
    x0_svd[first(layout_svd.idx_alpha)]    = alpha_mean   # trans_updown = trans_up^α; α=1 matches basis T
    x0_svd[first(layout_svd.idx_legendre)] = 1.0          # P0 = 1 → unit continuum
    x_a_svd      = copy(x0_svd)
    prior_sigma_svd = fill(prior_sigma_default, layout_svd.n_state)
    prior_sigma_svd[first(layout_svd.idx_alpha)] = alpha_sigma

    # PC priors: from loading variance or uninformative
    if pc_prior_mode == "loading_variance"
        n_prof = svd_basis.n_profiles
        for k in 1:n_pc
            sigma_k = svd_basis.S[k] / sqrt(Float64(n_prof)) * pc_sigma_scale
            prior_sigma_svd[k] = max(sigma_k, prior_min_sigma)
        end
    end  # else: leave as prior_sigma_default (uninformative)

    # Legendre P0/P1 prior for SVD: ratio fit at x0
    if use_leg01_svd && length(layout_svd.idx_legendre) >= 1
        y0 = fm_svd(x0_svd)
        ratio = Float64.(y_obs) ./ max.(abs.(y0), eps(Float64))
        z = _normalized_grid(λ_obs)
        A01 = hcat(ones(length(z)), z)
        w = Float64.(y_obs) .- minimum(Float64.(y_obs))
        w .+= max(maximum(w), 1.0) * 1e-6
        sv = sqrt.(w ./ maximum(w))
        c01 = (A01 .* sv) \ (ratio .* sv)
        leg0 = first(layout_svd.idx_legendre)
        x_a_svd[leg0] = c01[1]; x0_svd[leg0] = c01[1]
        prior_sigma_svd[leg0] = max(abs(c01[1]) * leg01_frac_svd, prior_min_sigma)
        if length(layout_svd.idx_legendre) >= 2
            leg1 = layout_svd.idx_legendre[2]
            x_a_svd[leg1] = c01[2]; x0_svd[leg1] = c01[2]
            prior_sigma_svd[leg1] = max(abs(c01[2]) * leg01_frac_svd, prior_min_sigma)
        end
    end
    if use_leghig_svd && length(layout_svd.idx_legendre) >= 3
        for j in 3:length(layout_svd.idx_legendre)
            prior_sigma_svd[layout_svd.idx_legendre[j]] = max(leg_higher_svd, prior_min_sigma)
        end
    end
    prior_sigma_svd[layout_svd.idx_sif] .= max(sif_sigma, prior_min_sigma)
    S_a_inv_svd = _spdiag_invvar(prior_sigma_svd)

    x_scale_svd = ones(Float64, layout_svd.n_state)
    for k in 1:n_pc
        x_scale_svd[k] = prior_sigma_svd[k]
    end
    x_scale_svd[first(layout_svd.idx_alpha)] = alpha_sigma

    lower_svd = fill(-Inf, layout_svd.n_state)
    upper_svd = fill( Inf, layout_svd.n_state)

    res_svd = run_lm_retrieval(
        fm_svd, jac_svd, Float64.(x0_svd), Float64.(x_a_svd),
        S_a_inv_svd, x_scale_svd,
        Float64.(y_obs), lower_svd, upper_svd, λ_obs;
        lm_kwargs...,
    )
    @printf("  SVD:  %d iter, %d fwd, %d jac, RMSE %.4e → %.4e (%.1f%% reduction)\n",
            res_svd.n_steps_done, res_svd.n_forward, res_svd.n_jacobian,
            res_svd.rmse_series[1], res_svd.rmse_series[end],
            100*(res_svd.rmse_series[1]-res_svd.rmse_series[end]) /
                max(abs(res_svd.rmse_series[1]), eps(Float64)))
    _print_retrieval_state(
        stdout, "SVD transmittance (LM) — state vector",
        res_svd, _svd_state_names(layout_svd),
    )

    # ── Print summary ────────────────────────────────────────────────────────
    println("\n" * "="^70)
    @printf("%-18s  %6s  %6s  %6s  %10s  %7s\n",
            "Method", "n_iter", "n_fwd", "n_jac", "RMSE_f", "red%")
    println("-"^60)
    for (label, res) in [("LUT xSec (LM)", res_xsec), ("SVD trans. (LM)", res_svd)]
        rmse0 = res.rmse_series[1]; rmsef = res.rmse_series[end]
        @printf("%-18s  %6d  %6d  %6d  %10.4e  %6.1f%%\n",
                label, res.n_steps_done, res.n_forward, res.n_jacobian,
                rmsef, 100*(rmse0-rmsef)/max(abs(rmse0),eps(Float64)))
    end

    # ── Write outputs ────────────────────────────────────────────────────────
    basis_png       = joinpath(_OUT_DIR, "svd_basis.png")
    rmse_png        = joinpath(_OUT_DIR, "rmse_vs_iteration.png")
    spectra_png     = joinpath(_OUT_DIR, "spectra_comparison.png")
    sif_trans_png   = joinpath(_OUT_DIR, "sif_trans_components.png")
    trans_12_png    = joinpath(_OUT_DIR, "transmittance_one_two.png")
    jac_png         = joinpath(_OUT_DIR, "jacobian_lut_hybrid_vs_forwarddiff.png")
    jac_svd_png     = joinpath(_OUT_DIR, "jacobian_svd_forwarddiff.png")
    ak_lut_png      = joinpath(_OUT_DIR, "averaging_kernel_lut.png")
    ak_svd_png      = joinpath(_OUT_DIR, "averaging_kernel_svd.png")
    summary_csv     = joinpath(_OUT_DIR, "comparison_summary.csv")
    n_cols          = 3

    x_lut_final = Float64.(res_xsec.x_final)
    S_e_inv_lut = if use_band_snr && !isnothing(ctx.band_snr_coeffs)
        make_Se_inv_from_snr(res_xsec.y_final, ctx.band_snr_coeffs)
    else
        spdiagm(0 => fill(1.0 / meas_sigma^2, length(y_obs)))
    end
    jac_hybrid_eval = make_hybrid_jacobian_evaluator(
        fm_xsec, ctx, solar_hres, layout_xsec; n_legendre = n_leg_xsec,
    )
    jac_fd_lut = make_jacobian_evaluator(fm_xsec, x_lut_final)
    J_hyb = Matrix(jac_hybrid_eval(x_lut_final))
    J_fd_lut = Matrix(jac_fd_lut(x_lut_final))
    σ_post_lut = posterior_sigma_from_jacobian(J_fd_lut, S_e_inv_lut, S_a_inv_xsec)
    J_hyb_norm = normalize_jacobian_by_sigma(J_hyb, σ_post_lut)
    J_fd_lut_norm = normalize_jacobian_by_sigma(J_fd_lut, σ_post_lut)
    jac_names_lut = state_names_simple(ctx; n_legendre = n_leg_xsec)
    _plot_jacobian_lut_hybrid_vs_fd(jac_png, λ_obs, J_hyb_norm, J_fd_lut_norm, jac_names_lut; n_cols = n_cols)
    jac_rmse_lut = sqrt(mean((J_hyb .- J_fd_lut) .^ 2))
    @printf("\nLUT Jacobian at final state — hybrid vs ForwardDiff element RMSE: %.6e\n", jac_rmse_lut)
    @printf("LUT posterior sigma: min %.3e, median %.3e, max %.3e\n",
            minimum(σ_post_lut), median(σ_post_lut), maximum(σ_post_lut))
    G_lut, A_lut = gain_and_averaging_kernel(J_fd_lut, S_e_inv_lut, S_a_inv_xsec)
    _plot_averaging_kernel_heatmap(
        ak_lut_png, A_lut, jac_names_lut;
        title = "LUT averaging kernel A = G K (final state)",
    )
    Abar_lut = mean(diag(A_lut))
    @printf("LUT averaging kernel diag: min %.3e, mean %.3e, max %.3e\n",
            minimum(diag(A_lut)), mean(diag(A_lut)), maximum(diag(A_lut)))
    @printf("LUT info metric (Ā = mean(diag(A))): %.6f\n", Abar_lut)

    x_svd_final = Float64.(res_svd.x_final)
    S_e_inv_svd = if use_band_snr && !isnothing(ctx.band_snr_coeffs)
        make_Se_inv_from_snr(res_svd.y_final, ctx.band_snr_coeffs)
    else
        spdiagm(0 => fill(1.0 / meas_sigma^2, length(y_obs)))
    end
    jac_fd_svd = make_jacobian_evaluator(fm_svd, x_svd_final)
    J_svd = Matrix(jac_fd_svd(x_svd_final))
    σ_post_svd = posterior_sigma_from_jacobian(J_svd, S_e_inv_svd, S_a_inv_svd)
    J_svd_norm = normalize_jacobian_by_sigma(J_svd, σ_post_svd)
    jac_names_svd = _svd_state_names(layout_svd)
    _plot_jacobian_svd_spectra(jac_svd_png, λ_obs, J_svd_norm, jac_names_svd; n_cols=n_cols)
    @printf("SVD Jacobian at final state — │∂y/∂x│ mean %.6e max %.6e\n",
            mean(abs.(J_svd)), maximum(abs.(J_svd)))
    @printf("SVD posterior sigma: min %.3e, median %.3e, max %.3e\n",
            minimum(σ_post_svd), median(σ_post_svd), maximum(σ_post_svd))
    G_svd, A_svd = gain_and_averaging_kernel(J_svd, S_e_inv_svd, S_a_inv_svd)
    _plot_averaging_kernel_heatmap(
        ak_svd_png, A_svd, jac_names_svd;
        title = "SVD averaging kernel A = G K (final state)",
    )
    Abar_svd = mean(diag(A_svd))
    @printf("SVD averaging kernel diag: min %.3e, mean %.3e, max %.3e\n",
            minimum(diag(A_svd)), mean(diag(A_svd)), maximum(diag(A_svd)))
    @printf("SVD info metric (Ā = mean(diag(A))): %.6f\n", Abar_svd)

    refl_xsec, sif_xsec_band = xsec_reflectance_and_sif_lres(
        ctx, solar_hres, res_xsec.x_final, layout_xsec, n_leg_xsec, model_variant,
    )
    refl_svd, sif_svd_band = svd_reflectance_and_sif_lres(
        ctx, solar_hres, svd_basis.PCs, res_svd.x_final, layout_svd;
        n_pc = n_pc, n_legendre = n_leg_svd, log_transform = log_trans,
    )

    λ_hres_T, T_up_xsec_final, T_updown_xsec_final = xsec_gas_transmittance_hres(ctx, res_xsec.x_final, layout_xsec)
    K_T = hasproperty(ctx, :kernel_rsr_out) ? ctx.kernel_rsr_out : ctx.kernel.RSR_out
    T_up_xsec_lres  = Vector(K_T * T_up_xsec_final)
    T_updown_xsec_lres = Vector(K_T * T_updown_xsec_final)
    _, T_up_svd, T_updown_svd = svd_transmittance_obs(
        ctx, svd_basis.PCs, res_svd.x_final, layout_svd;
        n_pc = n_pc, log_transform = log_trans,
    )

    _plot_svd_basis(basis_png, svd_basis, n_pc)
    _plot_rmse(rmse_png, res_xsec, res_svd)
    _plot_spectra(spectra_png, res_xsec, res_svd)
    _plot_sif_trans_components(sif_trans_png, λ_obs, refl_xsec, sif_xsec_band, refl_svd, sif_svd_band)
    _plot_transmittance_one_two(trans_12_png, λ_obs, T_up_xsec_lres, T_updown_xsec_lres, T_up_svd, T_updown_svd)
    _write_summary(summary_csv, res_xsec, res_svd)

    println("\nOutputs written to: $_OUT_DIR/")
    for f in [basis_png, rmse_png, spectra_png, sif_trans_png, trans_12_png,
              jac_png, jac_svd_png, ak_lut_png, ak_svd_png, summary_csv]
        println("  ", basename(f))
    end
end

main_compare()
