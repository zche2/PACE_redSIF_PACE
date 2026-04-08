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

    rmse_png = joinpath(_OUT_DIR, "rmse_vs_iteration.png")
    cost_png = joinpath(_OUT_DIR, "computational_cost.png")
    metrics_csv = joinpath(_OUT_DIR, "method_metrics.csv")
    state_csv = joinpath(_OUT_DIR, "final_state_by_method.csv")

    _plot_rmse(rmse_png, labels, results)
    _plot_cost(cost_png, labels, results)
    _write_metrics_csv(metrics_csv, labels, results)
    _write_final_state_csv(state_csv, labels, results)

    println("Wrote ", rmse_png)
    println("Wrote ", cost_png)
    println("Wrote ", metrics_csv)
    println("Wrote ", state_csv)
end

main_sanity()
#!/usr/bin/env julia
# Compare LM vs Gauss-Newton on the same single-pixel toy retrieval (no AddingSIF pipeline).
# Uses demo_example/Fit_toy_forward_model.jl and [fit] from the MWE TOML (PACE_MWE_CONFIG).

using DelimitedFiles
using Plots

const _SANITY_DIR = @__DIR__
const _DEMO_EXAMPLE = joinpath(_SANITY_DIR, "..", "..")
const _FIT_SCRIPT = joinpath(_DEMO_EXAMPLE, "Fit_toy_forward_model.jl")

include(_FIT_SCRIPT)

config_path = get(ENV, "PACE_MWE_CONFIG", joinpath(_DEMO_EXAMPLE, "Simple_PACE_xSecFit_MWE_zcheVer.toml"))
isfile(config_path) || error("Config not found: $config_path")
ENV["PACE_MWE_CONFIG"] = config_path

println("Single-pixel toy retrieval (no SIF injection): ", config_path)
println("Running GN ...")
rmse_gn = main(; silent=true, return_rmse_series_only=true, fit_method_override=:gn)
println("Running LM ...")
rmse_lm = main(; silent=true, return_rmse_series_only=true, fit_method_override=:lm)

println("rmse_lm = ", rmse_lm)
println("rmse_gn = ", rmse_gn)

# Plot RMSE vs iteration (0-based step index on x-axis)
p = plot(
    0:(length(rmse_lm) - 1),
    rmse_lm;
    label="LM",
    marker=:circle,
    markersize=3,
    lw=1.5,
    xlabel="Iteration (0 = prior state)",
    ylabel="RMSE (spectral)",
    title="Toy single-pixel retrieval: RMSE vs iteration",
    legend=:topright,
)
plot!(
    p,
    0:(length(rmse_gn) - 1),
    rmse_gn;
    label="GN",
    marker=:square,
    markersize=3,
    lw=1.5,
)
png_name = "rmse_vs_iteration.png"
savefig(p, png_name)
println("Wrote ", png_name)

println("Done. LM length=$(length(rmse_lm)), GN length=$(length(rmse_gn)) (includes prior at index 0).")
