#!/usr/bin/env julia
# Plot true T1 vs fitted T1 from ensemble NetCDFs.
# Fitted T1 = trans_up from PCs; α maps T1 → T2 via T2 = T1^α_coeff.
#
#   julia --project=. surrogate_meas/plot_t1_comparison.jl

using Statistics
using LinearAlgebra
using TOML
using NCDatasets
using Plots

const SCRIPT_DIR = @__DIR__
const REPO_ROOT = dirname(SCRIPT_DIR)
const OUT_DIR = joinpath(SCRIPT_DIR, "batch_ensemble")
const PLOT_DIR = joinpath(OUT_DIR, "plots")
const SVD_CONFIG = joinpath(SCRIPT_DIR, "configs", "svd_nPC15_npoly5.toml")
const TRUTH_NC = joinpath(OUT_DIR, "truth_ensemble.nc")
const RET_NC = joinpath(OUT_DIR, "retrieval_ensemble.nc")

include(joinpath(REPO_ROOT, "global_svd_fit_pipeline", "svd_retrieval", "svd_helpers.jl"))

"""α_coeff = 10/(1+e^{-α_raw}) + 1  (same as make_svd_forward_model_λ)."""
alpha_coeff_from_raw(α_raw::Real) = 10.0 / (1.0 + exp(-Float64(α_raw))) + 1.0

"""
Reconstruct T1 (trans_up) and T2 (trans_updown) from SVD state.
  T1 = exp(PCs·c)           [log_transform=true]
  T2 = T1^α_coeff           [α_coeff from α_raw]
"""
function reconstruct_T_from_state(
    state::AbstractMatrix{<:Real},  # (n_state, n_sample) or (n_sample, n_state)
    PCs::AbstractMatrix{<:Real},
    layout;
    log_transform::Bool=true,
)
    # accept either orientation
    n_state = layout.n_state
    if size(state, 1) == n_state
        X = state
    elseif size(state, 2) == n_state
        X = collect(state')
    else
        error("state size $(size(state)) incompatible with n_state=$n_state")
    end
    n_pc = layout.n_pc
    n_λ = size(PCs, 1)
    n_samp = size(X, 2)
    T1 = Matrix{Float64}(undef, n_λ, n_samp)
    T2 = Matrix{Float64}(undef, n_λ, n_samp)
    α_coeff = Vector{Float64}(undef, n_samp)
    PC = Float64.(PCs[:, 1:n_pc])
    for i in 1:n_samp
        c = @view X[layout.idx_pc, i]
        α_raw = X[first(layout.idx_alpha), i]
        α = alpha_coeff_from_raw(α_raw)
        α_coeff[i] = α
        t1 = log_transform ? exp.(PC * c) : 1.0 .+ PC * c
        t1 = max.(t1, eps(Float64))
        T1[:, i] .= t1
        T2[:, i] .= exp.(α .* log.(t1))
    end
    return (T1=T1, T2=T2, α_coeff=α_coeff)
end

function _as_band_sample(A, n_band, n_samp)
    if size(A) == (n_band, n_samp)
        return Matrix{Float64}(A)
    elseif size(A) == (n_samp, n_band)
        return Matrix{Float64}(A')
    else
        error("Unexpected array size $(size(A)); expected ($n_band,$n_samp) or ($n_samp,$n_band)")
    end
end

function main()
    mkpath(PLOT_DIR)
    isfile(TRUTH_NC) || error("Missing $TRUTH_NC")
    isfile(RET_NC) || error("Missing $RET_NC")
    isfile(SVD_CONFIG) || error("Missing $SVD_CONFIG")

    cfg = TOML.parsefile(SVD_CONFIG)
    svd_cfg = get(get(cfg, "fit", Dict()), "svd", Dict())
    data_cfg = get(cfg, "data", Dict{String, Any}())
    spectral = get(cfg, "spectral", Dict{String, Any}())
    n_pc = Int(get(svd_cfg, "n_pc", 15))
    n_leg = Int(get(svd_cfg, "n_legendre", 5))
    log_trans = Bool(get(svd_cfg, "svd_log_transform", true))
    λ_min = Float64(get(spectral, "lambda_min_nm", 640.0))
    λ_max = Float64(get(spectral, "lambda_max_nm", 756.0))
    base_dir = String(get(data_cfg, "base_dir", ""))
    summer = String(get(data_cfg, "summer_nc", ""))
    winter = String(get(data_cfg, "winter_nc", ""))
    summer_nc = isabspath(summer) ? summer : joinpath(base_dir, summer)
    winter_nc = isabspath(winter) ? winter : joinpath(base_dir, winter)

    ds_t = Dataset(TRUTH_NC)
    ds_r = Dataset(RET_NC)
    λ = Float64.(ds_t["wavelength"][:])
    n_band = length(λ)
    n_samp = Int(ds_t.dim["sample"])
    T1_true = _as_band_sample(ds_t["T1"][:], n_band, n_samp)
    T2_true = _as_band_sample(ds_t["T2"][:], n_band, n_samp)
    status = Int16.(ds_r["status"][:])
    state = ds_r["state"][:]
    # state may be (sample, state) from NetCDF
    if ndims(state) == 2 && size(state, 2) == Int(ds_r.dim["state"])
        state_mat = Matrix{Float64}(state')  # -> (n_state, n_sample)
    else
        state_mat = Matrix{Float64}(state)
    end
    close(ds_t)
    close(ds_r)

    println("Loading SVD PCs (n_pc=$n_pc, log_transform=$log_trans)…")
    basis = load_svd_basis(summer_nc, winter_nc, λ; λ_min=λ_min, λ_max=λ_max, n_pc=n_pc, log_transform=log_trans)
    layout = svd_state_layout(; n_pc=n_pc, n_legendre=n_leg, n_ev=1)
    recon = reconstruct_T_from_state(state_mat, basis.PCs, layout; log_transform=log_trans)

    ok = findall(status .== 1)
    println("Converged samples: $(length(ok)) / $n_samp")

    T1t = T1_true[:, ok]
    T1f = recon.T1[:, ok]
    T2t = T2_true[:, ok]
    T2f = recon.T2[:, ok]
    α = recon.α_coeff[ok]

    # ── mean-T1 scatter ───────────────────────────────────────────────────────
    m_t = vec(mean(T1t; dims=1))
    m_f = vec(mean(T1f; dims=1))
    A = hcat(ones(length(m_t)), m_t)
    β = A \ m_f
    r2 = 1 - sum((m_f .- A * β) .^ 2) / max(sum((m_f .- mean(m_f)) .^ 2), eps())
    bias = mean(m_f .- m_t)
    rmse = sqrt(mean((m_f .- m_t) .^ 2))
    lims = extrema(vcat(m_t, m_f))
    pad = 0.05 * (lims[2] - lims[1] + eps())
    lims = (lims[1] - pad, lims[2] + pad)

    p1 = scatter(m_t, m_f; ms=2, alpha=0.3, label="n=$(length(ok))",
                 xlabel="True mean T₁", ylabel="Fitted mean T₁ (from PCs)",
                 title="Mean T₁: bias=$(round(bias, digits=4)), RMSE=$(round(rmse, digits=4)), R²=$(round(r2, digits=3))",
                 size=(700, 650), legend=:topleft)
    plot!(p1, [lims[1], lims[2]], [lims[1], lims[2]]; color=:black, ls=:dash, label="1:1")
    xx = range(lims[1], lims[2]; length=50)
    plot!(p1, xx, β[1] .+ β[2] .* xx; color=:crimson, lw=2,
          label="fit: y=$(round(β[1], digits=3))+$(round(β[2], digits=3))x")
    xlims!(p1, lims); ylims!(p1, lims)
    savefig(p1, joinpath(PLOT_DIR, "T1_mean_scatter.png"))

    # ── T1 at continuum (~750 nm) and O2-B (~688 nm) ─────────────────────────
    i750 = argmin(abs.(λ .- 750.0))
    i688 = argmin(abs.(λ .- 688.0))
    function _scatter_band(ib, name)
        tt = vec(T1t[ib, :]); ff = vec(T1f[ib, :])
        lim = extrema(vcat(tt, ff))
        p = scatter(tt, ff; ms=1.5, alpha=0.25, label="",
                    xlabel="True T₁", ylabel="Fitted T₁",
                    title="T₁ @ $(round(λ[ib], digits=1)) nm ($name)", size=(550, 520))
        plot!(p, [lim[1], lim[2]], [lim[1], lim[2]]; color=:black, ls=:dash, label="1:1")
        return p
    end
    p_b = plot(_scatter_band(i750, "continuum"), _scatter_band(i688, "O₂-B");
               layout=(1, 2), size=(1100, 500))
    savefig(p_b, joinpath(PLOT_DIR, "T1_band_scatter.png"))

    # ── spectral overlays ─────────────────────────────────────────────────────
    ex = ok[1:min(6, length(ok))]
    plots_t1 = []
    for i in ex
        pk = plot(λ, T1_true[:, i]; label="True T₁", color=:navy, lw=2,
                  title="sample $i, α=$(round(recon.α_coeff[i], digits=2))",
                  xlabel="λ [nm]", ylabel="T₁", legend=:outerright, size=(500, 280))
        plot!(pk, λ, recon.T1[:, i]; label="Fitted T₁ (PCs)", color=:crimson, lw=1.5, ls=:dash)
        push!(plots_t1, pk)
    end
    savefig(plot(plots_t1...; layout=(2, 3), size=(1400, 700)), joinpath(PLOT_DIR, "T1_spectra_examples.png"))

    # ── T2 for reference: true includes solar; fitted is atm-only^α ───────────
    # Compare atm-proxy: true T2 / mean(T2 in continuum) shape vs fitted — or
    # show T2_fit vs true T2 only as spectral examples with a note.
    plots_t2 = []
    for i in ex
        pk = plot(λ, T2_true[:, i]; label="True T₂ (atm×solar)", color=:darkred, lw=2,
                  title="sample $i", xlabel="λ [nm]", ylabel="T₂", legend=:outerright, size=(500, 280))
        plot!(pk, λ, recon.T2[:, i]; label="Fitted T₂ = T₁^α", color=:orange, lw=1.5, ls=:dash)
        push!(plots_t2, pk)
    end
    savefig(plot(plots_t2...; layout=(2, 3), size=(1400, 700)), joinpath(PLOT_DIR, "T2_spectra_examples.png"))

    pα = histogram(α; bins=40, label="", xlabel="α_coeff = 10/(1+e^{-α_raw})+1",
                   ylabel="Count", title="Retrieved α coefficient", size=(700, 400))
    savefig(pα, joinpath(PLOT_DIR, "alpha_coeff_hist.png"))

    # spectral mean ± std of residual
    dT1 = T1f .- T1t
    μ = vec(mean(dT1; dims=2))
    σ = vec(std(dT1; dims=2))
    p_res = plot(λ, μ; ribbon=σ, label="mean ± 1σ", color=:navy, lw=2,
                 xlabel="Wavelength [nm]", ylabel="Fitted − true T₁",
                 title="T₁ residual spectrum", size=(1000, 400), legend=:outerright)
    hline!(p_res, [0.0]; color=:black, ls=:dash, label="")
    savefig(p_res, joinpath(PLOT_DIR, "T1_residual_spectrum.png"))

    open(joinpath(PLOT_DIR, "T1_comparison_summary.txt"), "w") do io
        println(io, "n_converged\t$(length(ok))")
        println(io, "n_pc\t$n_pc")
        println(io, "log_transform\t$log_trans")
        println(io, "mean_T1_bias\t$bias")
        println(io, "mean_T1_rmse\t$rmse")
        println(io, "mean_T1_r2\t$r2")
        println(io, "mean_T1_slope\t$(β[2])")
        println(io, "mean_T1_intercept\t$(β[1])")
        println(io, "median_alpha_coeff\t$(median(α))")
        println(io, "mean_alpha_coeff\t$(mean(α))")
    end

    println("Saved T₁ comparison plots under $PLOT_DIR")
    println("Mean T₁ — bias=$(round(bias, digits=4)), RMSE=$(round(rmse, digits=4)), R²=$(round(r2, digits=3)), slope=$(round(β[2], digits=3))")
    println("α_coeff — mean=$(round(mean(α), digits=3)), median=$(round(median(α), digits=3))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
