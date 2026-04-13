#!/usr/bin/env julia
# Benchmark comparison: LM vs GN vs L-BFGS-B on a single-pixel toy retrieval.
#
# Runs each method N_TRIALS times (default 5) to build a timing distribution,
# then produces:
#   benchmark_check/
#     timing_distribution.png   — violin/box-style strip plot of wall times per method
#     timing_log.csv             — raw trial wall times for all methods
#     rmse_vs_iteration.png      — RMSE trace per iteration (first trial only)
#     spectra_comparison.png     — reconstructed vs observed + residuals (first trial)
#     benchmark_summary.csv      — mean/median/std/min/max wall time + eval counts
#
# Run from repo root:
#   julia --project=. demo_example/AddingSIF/sanity_check/check_benchmark_methods.jl
#
# Override number of trials:
#   N_TRIALS=10 julia --project=. ...

using Plots
using Printf
using Statistics

const _SANITY_DIR   = @__DIR__
const _DEMO_EXAMPLE = joinpath(_SANITY_DIR, "..", "..")
const _FIT_SCRIPT   = joinpath(_DEMO_EXAMPLE, "Fit_toy_forward_model.jl")
const _OUT_DIR      = joinpath(_SANITY_DIR, "benchmark_check")

include(_FIT_SCRIPT)

# ── helpers ──────────────────────────────────────────────────────────────────

"""
Run `method_sym` once and return (result_namedtuple, wall_time_s).
Wall time is measured externally with `time_ns()` so it is independent of
the internal timer inside `main` (which measures only the solver loop).
"""
function _timed_run(method_sym::Symbol)
    t0 = time_ns()
    res = main(; silent=true, return_benchmark=true, fit_method_override=method_sym)
    wall_s = (time_ns() - t0) / 1e9
    return res, wall_s
end

# ── CSV writers ───────────────────────────────────────────────────────────────

function _write_timing_log(
    path::String,
    labels::Vector{String},
    all_times::Vector{Vector{Float64}},
)
    n_trials = maximum(length.(all_times))
    open(path, "w") do io
        println(io, join(labels, ","))
        for t in 1:n_trials
            row = [
                t <= length(all_times[i]) ? string(all_times[i][t]) : ""
                for i in eachindex(labels)
            ]
            println(io, join(row, ","))
        end
    end
end

function _write_benchmark_summary(
    path::String,
    labels::Vector{String},
    all_times::Vector{Vector{Float64}},
    results_first::Vector,          # benchmark result from trial 1 of each method
)
    open(path, "w") do io
        println(io, "method,n_trials,mean_s,median_s,std_s,min_s,max_s," *
                    "n_forward,n_jacobian,n_steps,rmse_final,rmse_reduction_pct," *
                    "converged,convergence_reason")
        for (i, (label, times, res)) in enumerate(zip(labels, all_times, results_first))
            rmse0   = res.rmse_series[1]
            rmsef   = res.rmse_series[end]
            redfrac = (rmse0 - rmsef) / max(abs(rmse0), eps(Float64))
            println(io, join([
                label,
                string(length(times)),
                string(mean(times)),
                string(median(times)),
                string(length(times) > 1 ? std(times) : 0.0),
                string(minimum(times)),
                string(maximum(times)),
                string(res.n_forward),
                string(res.n_jacobian),
                string(res.n_steps_done),
                string(rmsef),
                string(100 * redfrac),
                string(res.converged),
                repr(res.convergence_reason),
            ], ","))
        end
    end
end

# ── plots ─────────────────────────────────────────────────────────────────────

"""
Strip + box timing distribution plot (pure Plots.jl, no StatsPlots).
Each method gets a column; individual trial times are shown as dots,
the mean as a horizontal bar.
"""
function _plot_timing_distribution(
    path::String,
    labels::Vector{String},
    all_times::Vector{Vector{Float64}},
)
    method_colors = [:royalblue, :darkgreen, :firebrick]
    x_positions   = collect(1:length(labels))

    p = plot(
        xlabel    = "Method",
        ylabel    = "Wall time [s]",
        title     = "Timing distribution (all trials)",
        legend    = false,
        xticks    = (x_positions, labels),
        size      = (700, 500),
    )

    for (i, (times, col)) in enumerate(zip(all_times, method_colors))
        xi = x_positions[i]
        # Jitter x positions slightly so overlapping dots are visible
        jitter = (rand(length(times)) .- 0.5) .* 0.15
        scatter!(
            p,
            fill(xi, length(times)) .+ jitter,
            times;
            color      = col,
            markersize = 6,
            markerstrokewidth = 0,
            alpha      = 0.7,
        )
        # Mean bar
        m = mean(times)
        plot!(
            p,
            [xi - 0.25, xi + 0.25],
            [m, m];
            color = col,
            lw    = 3,
        )
        # Annotate mean
        annotate!(
            p,
            xi, m * 1.04,
            text(@sprintf("μ=%.2fs", m), 8, :center, col),
        )
    end

    savefig(p, path)
end

function _plot_rmse(path::String, labels::Vector{String}, results::Vector)
    method_colors  = [:royalblue, :darkgreen, :firebrick]
    method_markers = [:circle, :square, :diamond]
    p = plot(
        xlabel    = "Iteration (0 = prior state)",
        ylabel    = "RMSE (spectral)",
        title     = "RMSE vs iteration (trial 1)",
        legend    = :topright,
    )
    for (i, (label, res)) in enumerate(zip(labels, results))
        plot!(
            p,
            0:(length(res.rmse_series) - 1),
            res.rmse_series;
            label      = label,
            color      = method_colors[mod1(i, length(method_colors))],
            marker     = method_markers[mod1(i, length(method_markers))],
            markersize = 3,
            lw         = 1.8,
        )
    end
    savefig(p, path)
end

function _plot_spectra(path::String, labels::Vector{String}, results::Vector)
    λ     = results[1].wavelength
    y_obs = results[1].y_obs
    method_colors  = [:royalblue, :darkgreen, :firebrick]
    method_markers = [:circle, :square, :diamond]

    p_spec = plot(
        λ, y_obs;
        label          = "Observed",
        color          = :black,
        lw             = 2.0,
        xlabel         = "Wavelength [nm]",
        ylabel         = "R_toa",
        title          = "Reconstructed vs observed R_toa (trial 1)",
        legend         = :outertopright,
        legendfontsize = 7,
    )
    plot!(p_spec, λ, results[1].y_prior; label="Prior", color=:gray, lw=1.5, ls=:dash)
    for (i, (label, res)) in enumerate(zip(labels, results))
        plot!(
            p_spec, λ, res.y_final;
            label = label,
            color = method_colors[mod1(i, length(method_colors))],
            lw    = 1.5,
        )
    end

    p_resid = plot(
        λ, zeros(length(λ));
        label          = nothing,
        color          = :black,
        lw             = 1.0,
        ls             = :dot,
        xlabel         = "Wavelength [nm]",
        ylabel         = "Residual (obs − fit)",
        title          = "Spectral residuals at final state (trial 1)",
        legend         = :outertopright,
        legendfontsize = 7,
    )
    for (i, (label, res)) in enumerate(zip(labels, results))
        resid = y_obs .- res.y_final
        plot!(
            p_resid, λ, resid;
            label      = label,
            color      = method_colors[mod1(i, length(method_colors))],
            lw         = 1.5,
            marker     = method_markers[mod1(i, length(method_markers))],
            markersize = 2,
        )
    end

    p = plot(p_spec, p_resid; layout=(2, 1), size=(900, 800), left_margin=5Plots.mm)
    savefig(p, path)
end

# ── main ──────────────────────────────────────────────────────────────────────

function main_benchmark()
    config_path = get(
        ENV,
        "PACE_MWE_CONFIG",
        joinpath(_DEMO_EXAMPLE, "Simple_PACE_xSecFit_MWE_zcheVer.toml"),
    )
    isfile(config_path) || error("Config not found: $config_path")
    ENV["PACE_MWE_CONFIG"] = config_path
    mkpath(_OUT_DIR)

    n_trials = parse(Int, get(ENV, "N_TRIALS", "5"))

    methods = [:lm, :gn, :lbfgsb]
    labels  = ["LM", "GN", "L-BFGS-B"]

    println("="^70)
    println("Benchmark: LM vs GN vs L-BFGS-B  ($n_trials trials each)")
    println("Config: ", config_path)
    println("="^70)

    # all_times[method_idx] = Vector of wall times across trials
    all_times      = [Float64[] for _ in methods]
    # results_first[method_idx] = benchmark result from trial 1 (for plots)
    results_first  = Any[nothing for _ in methods]

    for (mi, (sym, label)) in enumerate(zip(methods, labels))
        println("\n── $label ($n_trials trials) ──")
        for t in 1:n_trials
            res, wall_s = _timed_run(sym)
            push!(all_times[mi], wall_s)
            if t == 1
                results_first[mi] = res
            end
            rmsef   = res.rmse_series[end]
            rmse0   = res.rmse_series[1]
            redfrac = (rmse0 - rmsef) / max(abs(rmse0), eps(Float64))
            @printf(
                "  trial %2d: wall=%.3fs  solver=%.3fs  rmse_f=%.4e  red=%.1f%%  n_jac=%d\n",
                t,
                wall_s,
                res.elapsed_wall_s,
                rmsef,
                100 * redfrac,
                res.n_jacobian,
            )
        end
        times = all_times[mi]
        @printf(
            "  → mean=%.3fs  median=%.3fs  std=%.3fs  [%.3f, %.3f]\n",
            mean(times), median(times),
            length(times) > 1 ? std(times) : 0.0,
            minimum(times), maximum(times),
        )
    end

    println("\n" * "="^70)
    println("Summary (trial 1 eval counts):")
    @printf("%-10s  %6s  %6s  %6s  %8s  %8s\n",
            "Method", "n_fwd", "n_jac", "n_iter", "rmse_f", "wall_s(μ)")
    println("-"^60)
    for (i, (label, res)) in enumerate(zip(labels, results_first))
        @printf("%-10s  %6d  %6d  %6d  %8.4e  %8.3f\n",
                label,
                res.n_forward, res.n_jacobian, res.n_steps_done,
                res.rmse_series[end],
                mean(all_times[i]))
    end

    # ── write outputs ─────────────────────────────────────────────────────────
    timing_png    = joinpath(_OUT_DIR, "timing_distribution.png")
    rmse_png      = joinpath(_OUT_DIR, "rmse_vs_iteration.png")
    spectra_png   = joinpath(_OUT_DIR, "spectra_comparison.png")
    timing_csv    = joinpath(_OUT_DIR, "timing_log.csv")
    summary_csv   = joinpath(_OUT_DIR, "benchmark_summary.csv")

    _plot_timing_distribution(timing_png,  labels, all_times)
    _plot_rmse(rmse_png,                   labels, results_first)
    _plot_spectra(spectra_png,             labels, results_first)
    _write_timing_log(timing_csv,          labels, all_times)
    _write_benchmark_summary(summary_csv,  labels, all_times, results_first)

    println("\nOutputs written to: ", _OUT_DIR, "/")
    for f in [timing_png, rmse_png, spectra_png, timing_csv, summary_csv]
        println("  ", basename(f))
    end
end

main_benchmark()
