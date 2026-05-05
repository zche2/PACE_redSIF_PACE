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
Build the SVD transmittance forward model closure.

State vector: `[c_1..c_{n_pc}, α, leg_0..leg_{n_leg}, sif_ev_0..sif_ev_{n_ev}]`

Forward model (high-res):
  linear:  trans_up = 1 + PCs_hres * c_vec          (SVD on transmittance)
  log:     trans_up = exp.(PCs_hres * c_vec)         (SVD on log-transmittance)
  trans_updown = trans_up^α = exp.(α * log(trans_up))   (path / air-mass scaling)
  ρ(λ)   = leg_basis * leg_coeff          (multiplicative continuum)
  SIF(λ)  = trans_up * sif_basis * sif_coeff  (SIF uses upwelling transmittance)
  y_hres  = solar * trans_updown * ρ / π + SIF
  y_lres  = K * y_hres
"""
function make_svd_forward_model(
    ctx,
    solar_hres::AbstractVector{<:Real},
    PCs_hres::Matrix{<:Real};
    n_pc::Int,
    n_legendre::Int,
    log_transform::Bool = false,
)
    n_pc > 0       || error("n_pc must be > 0")
    n_legendre >= 0 || error("n_legendre must be >= 0")
    size(PCs_hres, 2) >= n_pc || error("PCs_hres has fewer columns than n_pc")

    λ_hres        = collect(Float64.(ctx.λ_hres))
    K             = hasproperty(ctx, :kernel_rsr_out) ? ctx.kernel_rsr_out : ctx.kernel.RSR_out
    sif_basis_hres = Float64.(ctx.sif_basis_hres)
    n_ev          = size(sif_basis_hres, 2)
    solar         = Float64.(solar_hres)
    PCs           = Float64.(PCs_hres[:, 1:n_pc])    # (n_hres, n_pc)

    z_hres    = _normalized_grid(λ_hres)
    leg_basis = _legendre_design_matrix(z_hres, n_legendre)  # (n_hres, n_leg_coeff)
    layout    = svd_state_layout(; n_pc, n_legendre, n_ev)

    function fm_svd(x::AbstractVector)
        length(x) == layout.n_state ||
            error("SVD state vector length $(length(x)) ≠ expected $(layout.n_state)")
        c_vec        = @view x[layout.idx_pc]
        alpha_coeff  = x[first(layout.idx_alpha)]   # scalar: trans_updown = trans_up^alpha_coeff
        leg_coeff    = @view x[layout.idx_legendre]
        sif_coeff    = @view x[layout.idx_sif]

        # Reconstruct transmittance from PCs; path scaling trans_updown = trans_up^alpha_coeff
        trans_up     = log_transform ? exp.(PCs * c_vec) : 1.0 .+ PCs * c_vec  # (n_hres,)
        trans_updown = exp.(alpha_coeff .* log.(max.(trans_up, eps(Float64))))

        # Continuum: multiplicative polynomial
        rho_hres = leg_basis * leg_coeff                   # (n_hres,)

        # SIF: attenuated by same transmittance
        sif_hres = trans_up .* (sif_basis_hres * sif_coeff)  # (n_hres,)

        # TOA radiance
        y_hres = @. solar * trans_updown * rho_hres / π + sif_hres
        return K * y_hres
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

"""Split SVD TOA radiance into reflectance-path and SIF-path (same splitting as `make_svd_forward_model`)."""
function svd_reflectance_and_sif_lres(
    ctx,
    solar_hres::AbstractVector{<:Real},
    PCs_hres::Matrix{<:Real},
    x::AbstractVector{<:Real},
    layout_svd;
    n_pc::Int,
    n_legendre::Int,
    log_transform::Bool,
)
    λ_hres         = collect(Float64.(ctx.λ_hres))
    K              = hasproperty(ctx, :kernel_rsr_out) ? ctx.kernel_rsr_out : ctx.kernel.RSR_out
    sif_basis_hres = Float64.(ctx.sif_basis_hres)
    solar          = Float64.(collect(solar_hres))
    PCs            = Float64.(PCs_hres[:, 1:n_pc])

    z_hres    = _normalized_grid(λ_hres)
    leg_basis = _legendre_design_matrix(z_hres, n_legendre)

    c_vec       = @view x[layout_svd.idx_pc]
    alpha_coeff = x[first(layout_svd.idx_alpha)]
    leg_coeff   = @view x[layout_svd.idx_legendre]
    sif_coeff   = @view x[layout_svd.idx_sif]

    trans_up     = log_transform ? exp.(PCs * c_vec) : 1.0 .+ PCs * c_vec
    trans_updown = exp.(alpha_coeff .* log.(max.(trans_up, eps(Float64))))
    rho_hres     = leg_basis * leg_coeff
    sif_hres     = trans_up .* (sif_basis_hres * sif_coeff)

    y_refl_hres = @. solar * trans_updown * rho_hres / π

    return (K * y_refl_hres, K * sif_hres)
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

# ── Load transmittance NetCDF and build SVD basis ──────────────────────────────

function load_svd_basis(
    summer_nc::String,
    winter_nc::String,
    λ_hres::AbstractVector{<:Float64};
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

    # Interpolate each PC to the high-res wavelength grid
    n_hres   = length(λ_hres)
    PCs_hres = zeros(Float64, n_hres, size(U, 2))
    for k in eachindex(axes(U, 2))
        itp = LinearInterpolation(bands_sel, U[:, k]; extrapolation_bc=Flat())
        PCs_hres[:, k] .= itp.(λ_hres)
    end

    return (
        PCs_hres   = PCs_hres,     # (n_hres, n_svs)
        S          = S,             # singular values (unnormalised)
        S_norm     = S_norm,        # % variance
        n_profiles = n_profiles,
        bands_sel  = bands_sel,
        U          = U,             # (n_bands_sel, n_svs)
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
    use_band_snr   = haskey(fit_cfg, "use_band_snr") ? Bool(fit_cfg["use_band_snr"]) : kernel_wants_snr
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
        summer_nc, winter_nc, λ_hres;
        λ_min = λ_min, λ_max = λ_max, n_pc = n_pc, log_transform = log_trans,
    )

    fm_svd, layout_svd = make_svd_forward_model(
        ctx, solar_hres, svd_basis.PCs_hres;
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
    summary_csv     = joinpath(_OUT_DIR, "comparison_summary.csv")

    refl_xsec, sif_xsec_band = xsec_reflectance_and_sif_lres(
        ctx, solar_hres, res_xsec.x_final, layout_xsec, n_leg_xsec, model_variant,
    )
    refl_svd, sif_svd_band = svd_reflectance_and_sif_lres(
        ctx, solar_hres, svd_basis.PCs_hres, res_svd.x_final, layout_svd;
        n_pc = n_pc, n_legendre = n_leg_svd, log_transform = log_trans,
    )

    _plot_svd_basis(basis_png, svd_basis, n_pc)
    _plot_rmse(rmse_png, res_xsec, res_svd)
    _plot_spectra(spectra_png, res_xsec, res_svd)
    _plot_sif_trans_components(sif_trans_png, λ_obs, refl_xsec, sif_xsec_band, refl_svd, sif_svd_band)
    _write_summary(summary_csv, res_xsec, res_svd)

    println("\nOutputs written to: $_OUT_DIR/")
    for f in [basis_png, rmse_png, spectra_png, sif_trans_png, summary_csv]
        println("  ", basename(f))
    end
end

main_compare()
