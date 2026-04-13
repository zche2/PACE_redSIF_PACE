#!/usr/bin/env julia
# Hybrid retrieval: GN until it breaks, then L-BFGS-B continues from that state.
#
# Produces in hybrid_check/:
#   rmse_vs_iteration.png  — combined RMSE trace with phase boundary marked
#   spectra_comparison.png — observed vs prior vs each method final vs hybrid final
#   hybrid_summary.csv     — per-phase and overall metrics
#
# Run:
#   julia --project=. demo_example/AddingSIF/sanity_check/check_hybrid_gn_lbfgsb.jl

using Plots
using Printf
using Statistics

const _SANITY_DIR   = @__DIR__
const _DEMO_EXAMPLE = joinpath(_SANITY_DIR, "..", "..")
const _FIT_SCRIPT   = joinpath(_DEMO_EXAMPLE, "Fit_toy_forward_model.jl")
const _OUT_DIR      = joinpath(_SANITY_DIR, "hybrid_check")

include(_FIT_SCRIPT)

# ── run helpers ───────────────────────────────────────────────────────────────

function _run(method_sym::Symbol; x_init=nothing)
    main(;
        silent            = true,
        return_benchmark  = true,
        fit_method_override = method_sym,
        x_init            = x_init,
    )
end

# ── plotting ──────────────────────────────────────────────────────────────────

"""
Single RMSE-vs-iteration plot showing four traces:
  • Pure GN          (gray dashed)
  • Pure L-BFGS-B    (blue dashed)
  • Pure LM          (green dashed, reference)
  • Hybrid GN→LBFGSB (orange solid, bold)
A vertical dashed line marks the GN/L-BFGS-B handover iteration.
"""
function _plot_rmse(
    path::String,
    res_gn,
    res_lbfgsb,
    res_lm,
    res_hybrid_gn,
    res_hybrid_lbfgsb,
)
    n_gn    = length(res_gn.rmse_series)          # includes iteration 0
    # combined hybrid trace: GN iter 0..k, then L-BFGS-B iter k+1..end
    # (drop the duplicated iteration-0 of the second phase)
    hybrid_rmse = vcat(res_hybrid_gn.rmse_series, res_hybrid_lbfgsb.rmse_series[2:end])
    handover    = n_gn - 1                         # iteration index (0-based) of handover

    p = plot(
        xlabel    = "Iteration (0 = prior)",
        ylabel    = "RMSE (spectral)",
        title     = "RMSE vs iteration: GN → L-BFGS-B hybrid",
        legend    = :topright,
        size      = (900, 500),
    )

    # reference pure methods
    plot!(p, 0:(length(res_lm.rmse_series)-1),     res_lm.rmse_series;
          label="Pure LM",      color=:darkgreen, lw=1.5, ls=:dash, marker=:circle,    markersize=3)
    plot!(p, 0:(length(res_gn.rmse_series)-1),     res_gn.rmse_series;
          label="Pure GN",      color=:gray,      lw=1.5, ls=:dash, marker=:square,    markersize=3)
    plot!(p, 0:(length(res_lbfgsb.rmse_series)-1), res_lbfgsb.rmse_series;
          label="Pure L-BFGS-B",color=:royalblue, lw=1.5, ls=:dash, marker=:diamond,  markersize=3)

    # hybrid trace
    plot!(p, 0:(length(hybrid_rmse)-1), hybrid_rmse;
          label="Hybrid GN→L-BFGS-B", color=:darkorange, lw=2.5, ls=:solid,
          marker=:utriangle, markersize=4)

    # handover boundary
    vline!(p, [handover]; color=:darkorange, ls=:dot, lw=1.5, label="GN→LBFGSB handover")
    annotate!(p, handover + 0.3, maximum(hybrid_rmse) * 0.97,
              text("handover\niter $(handover)", 8, :left, :darkorange))

    savefig(p, path)
    return hybrid_rmse, handover
end

"""
Two-panel spectral comparison: observed vs prior vs pure methods vs hybrid.
"""
function _plot_spectra(
    path::String,
    res_gn, res_lbfgsb, res_lm, res_hybrid_lbfgsb,
)
    λ     = res_lm.wavelength
    y_obs = res_lm.y_obs

    p_spec = plot(
        λ, y_obs;
        label="Observed", color=:black, lw=2.0,
        xlabel="Wavelength [nm]", ylabel="R_toa",
        title="Spectra: pure methods vs hybrid (final state)",
        legend=:outertopright, legendfontsize=7,
    )
    plot!(p_spec, λ, res_lm.y_prior; label="Prior", color=:gray, lw=1.5, ls=:dash)
    plot!(p_spec, λ, res_lm.y_final;      label="LM final",         color=:darkgreen,  lw=1.2, ls=:dot)
    plot!(p_spec, λ, res_gn.y_final;      label="GN final",         color=:gray,       lw=1.2, ls=:dot)
    plot!(p_spec, λ, res_lbfgsb.y_final;  label="L-BFGS-B final",   color=:royalblue,  lw=1.2, ls=:dot)
    plot!(p_spec, λ, res_hybrid_lbfgsb.y_final; label="Hybrid final", color=:darkorange, lw=2.0, ls=:solid)

    p_resid = plot(
        λ, zeros(length(λ));
        label=nothing, color=:black, lw=1.0, ls=:dot,
        xlabel="Wavelength [nm]", ylabel="Residual (obs − fit)",
        title="Spectral residuals", legend=:outertopright, legendfontsize=7,
    )
    for (label, res, col) in [
            ("LM",          res_lm,              :darkgreen),
            ("GN",          res_gn,              :gray),
            ("L-BFGS-B",    res_lbfgsb,          :royalblue),
            ("Hybrid",      res_hybrid_lbfgsb,   :darkorange),
        ]
        plot!(p_resid, λ, y_obs .- res.y_final; label=label, color=col, lw=1.5)
    end

    p = plot(p_spec, p_resid; layout=(2,1), size=(1000, 850), left_margin=5Plots.mm)
    savefig(p, path)
end

# ── CSV ───────────────────────────────────────────────────────────────────────

function _write_summary(
    path::String,
    res_gn, res_lbfgsb, res_lm,
    res_hybrid_gn, res_hybrid_lbfgsb,
    hybrid_rmse,
)
    open(path, "w") do io
        println(io, "label,n_iter,n_forward,n_jacobian,rmse_prior,rmse_final,rmse_reduction_pct,elapsed_s,notes")
        function row(label, res; extra_fwd=0, extra_jac=0, extra_s=0.0, notes="")
            rmse0 = res.rmse_series[1]
            rmsef = res.rmse_series[end]
            red   = 100 * (rmse0 - rmsef) / max(abs(rmse0), eps(Float64))
            println(io, join([label, res.n_steps_done,
                              res.n_forward + extra_fwd,
                              res.n_jacobian + extra_jac,
                              rmse0, rmsef, red,
                              res.elapsed_wall_s + extra_s, notes], ","))
        end
        row("Pure LM",       res_lm)
        row("Pure GN",       res_gn)
        row("Pure L-BFGS-B", res_lbfgsb)
        # hybrid: combine counts from both phases
        n_iter_hybrid = res_hybrid_gn.n_steps_done + res_hybrid_lbfgsb.n_steps_done
        rmse_prior    = res_hybrid_gn.rmse_series[1]
        rmse_final    = hybrid_rmse[end]
        red           = 100 * (rmse_prior - rmse_final) / max(abs(rmse_prior), eps(Float64))
        elapsed       = res_hybrid_gn.elapsed_wall_s + res_hybrid_lbfgsb.elapsed_wall_s
        println(io, join(["Hybrid GN→LBFGSB", n_iter_hybrid,
                          res_hybrid_gn.n_forward  + res_hybrid_lbfgsb.n_forward,
                          res_hybrid_gn.n_jacobian + res_hybrid_lbfgsb.n_jacobian,
                          rmse_prior, rmse_final, red, elapsed,
                          "GN phase: $(res_hybrid_gn.n_steps_done) iter; " *
                          "LBFGSB phase: $(res_hybrid_lbfgsb.n_steps_done) iter"], ","))
    end
end

# ── main ──────────────────────────────────────────────────────────────────────

function main_hybrid()
    config_path = get(
        ENV, "PACE_MWE_CONFIG",
        joinpath(_DEMO_EXAMPLE, "Simple_PACE_xSecFit_MWE_zcheVer.toml"),
    )
    isfile(config_path) || error("Config not found: $config_path")
    ENV["PACE_MWE_CONFIG"] = config_path
    mkpath(_OUT_DIR)

    println("="^70)
    println("Hybrid GN → L-BFGS-B fine-tuning comparison")
    println("Config: ", config_path)
    println("="^70)

    # ── pure baselines ────────────────────────────────────────────────────────
    println("\n[1/5] Pure LM ...")
    res_lm      = _run(:lm)
    println("[2/5] Pure GN ...")
    res_gn      = _run(:gn)
    println("[3/5] Pure L-BFGS-B ...")
    res_lbfgsb  = _run(:lbfgsb)

    # ── hybrid: GN phase (same as pure GN) ───────────────────────────────────
    println("[4/5] Hybrid — GN phase ...")
    res_hybrid_gn = res_gn      # reuse; GN runs identically from prior

    # ── hybrid: L-BFGS-B warm-started from GN's final x ─────────────────────
    println("[5/5] Hybrid — L-BFGS-B phase (warm-started from GN final state) ...")
    res_hybrid_lbfgsb = _run(:lbfgsb; x_init = res_hybrid_gn.x_final)

    # ── print summary ─────────────────────────────────────────────────────────
    println("\n" * "="^70)
    @printf("%-22s  %6s  %6s  %6s  %10s  %8s\n",
            "Method", "n_iter", "n_fwd", "n_jac", "rmse_f", "reduction")
    println("-"^70)
    for (label, res) in [
            ("Pure LM",       res_lm),
            ("Pure GN",       res_gn),
            ("Pure L-BFGS-B", res_lbfgsb),
        ]
        rmse0 = res.rmse_series[1]; rmsef = res.rmse_series[end]
        red   = 100*(rmse0-rmsef)/max(abs(rmse0),eps(Float64))
        @printf("%-22s  %6d  %6d  %6d  %10.4e  %7.1f%%\n",
                label, res.n_steps_done, res.n_forward, res.n_jacobian, rmsef, red)
    end
    # hybrid combined
    n_gn_iters    = res_hybrid_gn.n_steps_done
    hybrid_rmse_f = res_hybrid_lbfgsb.rmse_series[end]
    rmse0         = res_hybrid_gn.rmse_series[1]
    red           = 100*(rmse0 - hybrid_rmse_f)/max(abs(rmse0), eps(Float64))
    n_fwd_total   = res_hybrid_gn.n_forward  + res_hybrid_lbfgsb.n_forward
    n_jac_total   = res_hybrid_gn.n_jacobian + res_hybrid_lbfgsb.n_jacobian
    n_iter_total  = res_hybrid_gn.n_steps_done + res_hybrid_lbfgsb.n_steps_done
    @printf("%-22s  %6d  %6d  %6d  %10.4e  %7.1f%%   (GN: %d iter, LBFGSB: %d iter)\n",
            "Hybrid GN→L-BFGS-B",
            n_iter_total, n_fwd_total, n_jac_total, hybrid_rmse_f, red,
            n_gn_iters, res_hybrid_lbfgsb.n_steps_done)

    # ── outputs ───────────────────────────────────────────────────────────────
    rmse_png    = joinpath(_OUT_DIR, "rmse_vs_iteration.png")
    spectra_png = joinpath(_OUT_DIR, "spectra_comparison.png")
    summary_csv = joinpath(_OUT_DIR, "hybrid_summary.csv")

    hybrid_rmse, _ = _plot_rmse(rmse_png,    res_gn, res_lbfgsb, res_lm,
                                              res_hybrid_gn, res_hybrid_lbfgsb)
    _plot_spectra(spectra_png, res_gn, res_lbfgsb, res_lm, res_hybrid_lbfgsb)
    _write_summary(summary_csv, res_gn, res_lbfgsb, res_lm,
                   res_hybrid_gn, res_hybrid_lbfgsb, hybrid_rmse)

    println("\nOutputs written to: ", _OUT_DIR, "/")
    for f in [rmse_png, spectra_png, summary_csv]
        println("  ", basename(f))
    end
end

main_hybrid()
