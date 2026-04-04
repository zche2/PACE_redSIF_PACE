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
