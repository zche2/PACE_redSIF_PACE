#!/usr/bin/env julia
# Compare L-BFGS-B performance across different memory sizes (lbfgs_m).
# Metrics: RMSE trace, objective trace, wall time, n_forward, n_jacobian.
#
# Run from repo root:
#   julia --project=. demo_example/AddingSIF/sanity_check/check_lbfgsb_m.jl

using Plots
using Printf
using TOML

const _SANITY_DIR  = @__DIR__
const _DEMO_EXAMPLE = joinpath(_SANITY_DIR, "..", "..")
const _FIT_SCRIPT   = joinpath(_DEMO_EXAMPLE, "Fit_toy_forward_model.jl")
const _OUT_DIR      = joinpath(_SANITY_DIR, "lbfgsb_m_check")

include(_FIT_SCRIPT)

# ── helpers ──────────────────────────────────────────────────────────────────

"""
Write a patched copy of the base TOML with `fit.lbfgs_m` overridden,
call `main(; silent, return_benchmark, fit_method_override=:lbfgsb)`,
then restore the previous ENV value.
"""
function _run_lbfgsb_m(base_config_path::String, m::Int)
    cfg = TOML.parsefile(base_config_path)
    cfg["fit"]["lbfgs_m"] = m
    cfg["fit"]["fit_method"] = "lbfgsb"

    tmp = tempname() * ".toml"
    open(tmp, "w") do io
        TOML.print(io, cfg)
    end

    prev = get(ENV, "PACE_MWE_CONFIG", nothing)
    ENV["PACE_MWE_CONFIG"] = tmp
    result = try
        main(; silent=true, return_benchmark=true, fit_method_override=:lbfgsb)
    finally
        if isnothing(prev)
            delete!(ENV, "PACE_MWE_CONFIG")
        else
            ENV["PACE_MWE_CONFIG"] = prev
        end
        rm(tmp; force=true)
    end
    return result
end

# ── CSV writers ──────────────────────────────────────────────────────────────

function _write_metrics_csv(path::String, m_vals::Vector{Int}, results::Vector)
    open(path, "w") do io
        println(io, "lbfgs_m,rmse_prior,rmse_final,rmse_reduction_frac," *
                    "elapsed_wall_s,n_forward,n_jacobian,converged," *
                    "failed_step,failed_error,convergence_reason,n_steps_done")
        for (m, res) in zip(m_vals, results)
            rmse0   = res.rmse_series[1]
            rmsef   = res.rmse_series[end]
            redfrac = (rmse0 - rmsef) / max(abs(rmse0), eps(Float64))
            println(io, join([
                string(m),
                string(rmse0), string(rmsef), string(redfrac),
                string(res.elapsed_wall_s),
                string(res.n_forward), string(res.n_jacobian),
                string(res.converged),
                string(res.failed_step), repr(res.failed_error),
                repr(res.convergence_reason),
                string(res.n_steps_done),
            ], ","))
        end
    end
end

# ── plots ────────────────────────────────────────────────────────────────────

function _plot_rmse_traces(path::String, m_vals::Vector{Int}, results::Vector)
    markers = [:circle, :square, :diamond, :utriangle, :star5]
    colors  = [:royalblue, :darkgreen, :firebrick, :darkorange, :purple]
    p = plot(
        xlabel = "Iteration (0 = prior state)",
        ylabel = "RMSE (spectral)",
        title  = "L-BFGS-B RMSE vs iteration (varying m)",
        legend = :topright,
    )
    for (i, (m, res)) in enumerate(zip(m_vals, results))
        plot!(
            p,
            0:(length(res.rmse_series) - 1),
            res.rmse_series;
            label    = "m = $m",
            marker   = markers[mod1(i, length(markers))],
            markersize = 3,
            lw       = 1.8,
            color    = colors[mod1(i, length(colors))],
        )
    end
    savefig(p, path)
end

function _plot_cost(path::String, m_vals::Vector{Int}, results::Vector)
    labels  = ["m=$m" for m in m_vals]
    x       = collect(1:length(m_vals))
    wall    = [r.elapsed_wall_s for r in results]
    n_fwd   = Float64[r.n_forward  for r in results]
    n_jac   = Float64[r.n_jacobian for r in results]
    n_steps = Float64[r.n_steps_done for r in results]

    p_wall = bar(
        x, wall;
        xticks = (x, labels),
        ylabel = "Wall time [s]",
        title  = "Wall time by m",
        legend = false,
        color  = :steelblue,
    )

    # Manual grouped bar: offset two series by ±0.2 with bar width 0.35
    bw = 0.35
    p_evals = bar(
        x .- 0.2, n_fwd;
        bar_width = bw,
        xticks    = (x, labels),
        label     = "n_forward",
        ylabel    = "Evaluation count",
        title     = "Evaluations by m",
        legend    = :topright,
        color     = :royalblue,
    )
    bar!(
        p_evals,
        x .+ 0.2, n_jac;
        bar_width = bw,
        label     = "n_jacobian",
        color     = :firebrick,
    )

    p_steps = bar(
        x, n_steps;
        xticks = (x, labels),
        ylabel = "Iterations accepted",
        title  = "Iterations done by m",
        legend = false,
        color  = :darkorange,
    )
    p = plot(p_wall, p_evals, p_steps; layout=(3, 1), size=(900, 1000))
    savefig(p, path)
end

function _plot_evals_per_iter(path::String, m_vals::Vector{Int}, results::Vector)
    labels  = ["m=$m" for m in m_vals]
    x       = collect(1:length(m_vals))
    fwd_per = [r.n_forward  / max(r.n_steps_done, 1) for r in results]
    jac_per = [r.n_jacobian / max(r.n_steps_done, 1) for r in results]

    bw = 0.35
    p = bar(
        x .- 0.2, fwd_per;
        bar_width = bw,
        xticks    = (x, labels),
        label     = "fwd/iter",
        ylabel    = "Evaluations per accepted iteration",
        title     = "L-BFGS-B cost per iteration (varying m)",
        legend    = :topright,
        color     = :royalblue,
        size      = (800, 500),
    )
    bar!(
        p,
        x .+ 0.2, jac_per;
        bar_width = bw,
        label     = "jac/iter",
        color     = :firebrick,
    )
    savefig(p, path)
end

function _plot_rmse_vs_cost(path::String, m_vals::Vector{Int}, results::Vector)
    rmsef  = [r.rmse_series[end] for r in results]
    n_jac  = [r.n_jacobian for r in results]
    wall   = [r.elapsed_wall_s for r in results]

    p1 = scatter(
        n_jac, rmsef;
        xlabel     = "n_jacobian (total)",
        ylabel     = "Final RMSE",
        title      = "Final RMSE vs Jacobian cost",
        legend     = false,
        series_annotations = ["  m=$m" for m in m_vals],
        markersize = 6,
        color      = :royalblue,
    )
    p2 = scatter(
        wall, rmsef;
        xlabel     = "Wall time [s]",
        ylabel     = "Final RMSE",
        title      = "Final RMSE vs Wall time",
        legend     = false,
        series_annotations = ["  m=$m" for m in m_vals],
        markersize = 6,
        color      = :darkgreen,
    )
    p = plot(p1, p2; layout=(1, 2), size=(1000, 450))
    savefig(p, path)
end

# ── main ─────────────────────────────────────────────────────────────────────

function main_m_check()
    config_path = get(
        ENV,
        "PACE_MWE_CONFIG",
        joinpath(_DEMO_EXAMPLE, "Simple_PACE_xSecFit_MWE_zcheVer.toml"),
    )
    isfile(config_path) || error("Config not found: $config_path")
    mkpath(_OUT_DIR)

    m_vals = [3, 5, 10]

    println("L-BFGS-B memory-size comparison (m ∈ $m_vals)")
    println("Config: ", config_path)
    println()

    results = Any[]
    for m in m_vals
        println("Running L-BFGS-B with m = $m ...")
        push!(results, _run_lbfgsb_m(config_path, m))
    end

    println()
    println("─"^75)
    @printf("%-8s  %-12s  %-12s  %-10s  %-6s  %-7s  %-7s  %-6s\n",
            "m", "rmse_prior", "rmse_final", "reduction", "time_s",
            "n_fwd", "n_jac", "n_iter")
    println("─"^75)
    for (m, res) in zip(m_vals, results)
        rmse0   = res.rmse_series[1]
        rmsef   = res.rmse_series[end]
        redfrac = (rmse0 - rmsef) / max(abs(rmse0), eps(Float64))
        @printf(
            "m=%-6d  %-12.6e  %-12.6e  %8.2f%%  %6.2fs  %-7d  %-7d  %-6d\n",
            m, rmse0, rmsef, 100 * redfrac,
            res.elapsed_wall_s, res.n_forward, res.n_jacobian, res.n_steps_done,
        )
    end
    println("─"^75)

    rmse_png    = joinpath(_OUT_DIR, "rmse_vs_iteration.png")
    cost_png    = joinpath(_OUT_DIR, "computational_cost.png")
    piter_png   = joinpath(_OUT_DIR, "evals_per_iter.png")
    scatter_png = joinpath(_OUT_DIR, "rmse_vs_cost.png")
    metrics_csv = joinpath(_OUT_DIR, "lbfgsb_m_metrics.csv")

    _plot_rmse_traces(rmse_png,   m_vals, results)
    _plot_cost(cost_png,          m_vals, results)
    _plot_evals_per_iter(piter_png, m_vals, results)
    _plot_rmse_vs_cost(scatter_png, m_vals, results)
    _write_metrics_csv(metrics_csv, m_vals, results)

    println()
    println("Outputs written to ", _OUT_DIR, "/")
    for f in [rmse_png, cost_png, piter_png, scatter_png, metrics_csv]
        println("  ", basename(f))
    end
end

main_m_check()
