#!/usr/bin/env julia
# Compare LM vs GN vs L-BFGS-B on one single-pixel toy retrieval (no AddingSIF pipeline).
# Uses demo_example/Fit_toy_forward_model.jl and [fit] from the MWE TOML (PACE_MWE_CONFIG).

using Plots
using Printf

const _SANITY_DIR = @__DIR__
const _DEMO_EXAMPLE = joinpath(_SANITY_DIR, "..", "..")
const _FIT_SCRIPT = joinpath(_DEMO_EXAMPLE, "Fit_toy_forward_model.jl")
const _OUT_DIR = joinpath(_SANITY_DIR, "iteration_method_check")

include(_FIT_SCRIPT)

function _run_method(method_sym::Symbol)
    return main(; silent=true, return_benchmark=true, fit_method_override=method_sym)
end

function _write_metrics_csv(path::String, labels::Vector{String}, results::Vector)
    open(path, "w") do io
        println(io, "method,rmse_prior,rmse_final,rmse_reduction_frac,elapsed_wall_s,n_forward,n_jacobian,converged,failed_step,failed_error,convergence_reason,n_steps_done")
        for (label, res) in zip(labels, results)
            rmse0 = res.rmse_series[1]
            rmsef = res.rmse_series[end]
            redfrac = (rmse0 - rmsef) / max(abs(rmse0), eps(Float64))
            println(
                io,
                join(
                    [
                        label,
                        string(rmse0),
                        string(rmsef),
                        string(redfrac),
                        string(res.elapsed_wall_s),
                        string(res.n_forward),
                        string(res.n_jacobian),
                        string(res.converged),
                        string(res.failed_step),
                        repr(res.failed_error),
                        repr(res.convergence_reason),
                        string(res.n_steps_done),
                    ],
                    ",",
                ),
            )
        end
    end
end

function _write_final_state_csv(path::String, labels::Vector{String}, results::Vector)
    state_names = results[1].state_names
    open(path, "w") do io
        println(io, join(vcat(["method"], state_names), ","))
        for (label, res) in zip(labels, results)
            row = vcat([label], string.(res.x_final))
            println(io, join(row, ","))
        end
    end
end

function _plot_rmse(path::String, labels::Vector{String}, results::Vector)
    markers = [:circle, :square, :diamond]
    p = plot(
        xlabel="Iteration (0 = prior state)",
        ylabel="RMSE (spectral)",
        title="Toy single-pixel retrieval: RMSE vs iteration",
        legend=:topright,
    )
    for (i, (label, res)) in enumerate(zip(labels, results))
        plot!(
            p,
            0:(length(res.rmse_series) - 1),
            res.rmse_series;
            label=label,
            marker=markers[mod1(i, length(markers))],
            markersize=3,
            lw=1.8,
        )
    end
    savefig(p, path)
end

"""
Two-panel spectral comparison across all three retrieval methods.

Top panel   — observed R_toa + prior-state spectrum + final reconstructed
              spectrum for each method, all on the same axes.
Bottom panel — residuals (observed − reconstructed) at the final state for
               each method, plus a zero reference line.
"""
function _plot_spectra(path::String, labels::Vector{String}, results::Vector)
    λ     = results[1].wavelength
    y_obs = results[1].y_obs

    method_colors  = [:royalblue, :darkgreen, :firebrick]
    method_markers = [:circle, :square, :diamond]

    # ── top panel: observed + prior + all finals ──────────────────────────
    p_spec = plot(
        λ, y_obs;
        label      = "Observed",
        color      = :black,
        lw         = 2.0,
        ls         = :solid,
        xlabel     = "Wavelength [nm]",
        ylabel     = "R_toa",
        title      = "Reconstructed vs observed R_toa",
        legend     = :outertopright,
        legendfontsize = 7,
    )
    plot!(
        p_spec, λ, results[1].y_prior;
        label  = "Prior",
        color  = :gray,
        lw     = 1.5,
        ls     = :dash,
    )
    for (i, (label, res)) in enumerate(zip(labels, results))
        plot!(
            p_spec, λ, res.y_final;
            label  = label,
            color  = method_colors[mod1(i, length(method_colors))],
            lw     = 1.5,
            ls     = :solid,
        )
    end

    # ── bottom panel: residuals ───────────────────────────────────────────
    p_resid = plot(
        λ, zeros(length(λ));
        label      = nothing,
        color      = :black,
        lw         = 1.0,
        ls         = :dot,
        xlabel     = "Wavelength [nm]",
        ylabel     = "Residual (obs − fit)",
        title      = "Spectral residuals at final state",
        legend     = :outertopright,
        legendfontsize = 7,
    )
    for (i, (label, res)) in enumerate(zip(labels, results))
        resid = y_obs .- res.y_final
        plot!(
            p_resid, λ, resid;
            label  = label,
            color  = method_colors[mod1(i, length(method_colors))],
            lw     = 1.5,
            marker = method_markers[mod1(i, length(method_markers))],
            markersize = 2,
        )
    end

    p = plot(p_spec, p_resid; layout=(2, 1), size=(900, 800), left_margin=5Plots.mm)
    savefig(p, path)
end

function _plot_cost(path::String, labels::Vector{String}, results::Vector)
    x = 1:length(labels)
    wall = [r.elapsed_wall_s for r in results]
    n_fwd = [r.n_forward for r in results]
    n_jac = [r.n_jacobian for r in results]

    p1 = bar(
        x,
        wall;
        xticks=(x, labels),
        ylabel="Wall time [s]",
        title="Computational cost (wall time)",
        legend=false,
    )
    p2 = bar(
        x,
        [n_fwd n_jac];
        xticks=(x, labels),
        label=["n_forward" "n_jacobian"],
        ylabel="Evaluation count",
        title="Computational cost (evaluations)",
        legend=:topright,
    )
    p = plot(p1, p2; layout=(2, 1), size=(900, 800))
    savefig(p, path)
end

function main_sanity()
    config_path = get(ENV, "PACE_MWE_CONFIG", joinpath(_DEMO_EXAMPLE, "Simple_PACE_xSecFit_MWE_zcheVer.toml"))
    isfile(config_path) || error("Config not found: $config_path")
    ENV["PACE_MWE_CONFIG"] = config_path
    mkpath(_OUT_DIR)

    methods = [:lm, :gn, :lbfgsb]
    labels = ["LM", "GN", "L-BFGS-B"]
    println("Single-pixel toy retrieval (no SIF injection): ", config_path)

    results = Any[]
    for (sym, label) in zip(methods, labels)
        println("Running ", label, " ...")
        push!(results, _run_method(sym))
    end

    for (label, res) in zip(labels, results)
        rmse0 = res.rmse_series[1]
        rmsef = res.rmse_series[end]
        redfrac = (rmse0 - rmsef) / max(abs(rmse0), eps(Float64))
        @printf(
            "%-8s rmse0=%.6e  rmsef=%.6e  reduction=%.2f%%  time=%.3fs  n_fwd=%d  n_jac=%d\n",
            label,
            rmse0,
            rmsef,
            100 * redfrac,
            res.elapsed_wall_s,
            res.n_forward,
            res.n_jacobian,
        )
    end

    rmse_png    = joinpath(_OUT_DIR, "rmse_vs_iteration.png")
    cost_png    = joinpath(_OUT_DIR, "computational_cost.png")
    spectra_png = joinpath(_OUT_DIR, "spectra_comparison.png")
    metrics_csv = joinpath(_OUT_DIR, "method_metrics.csv")
    state_csv   = joinpath(_OUT_DIR, "final_state_by_method.csv")

    _plot_rmse(rmse_png, labels, results)
    _plot_cost(cost_png, labels, results)
    _plot_spectra(spectra_png, labels, results)
    _write_metrics_csv(metrics_csv, labels, results)
    _write_final_state_csv(state_csv, labels, results)

    println("Wrote ", rmse_png)
    println("Wrote ", cost_png)
    println("Wrote ", spectra_png)
    println("Wrote ", metrics_csv)
    println("Wrote ", state_csv)
end

main_sanity()
