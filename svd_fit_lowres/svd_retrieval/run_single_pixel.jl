#!/usr/bin/env julia
# Single-pixel SVD retrieval + spectra diagnostic plot.
#
# Usage:
#   julia --project=. svd_fit_lowres/svd_retrieval/run_single_pixel.jl path/to/global_fit_pipeline.toml
#
# Requires [svd_retrieval] granule_id plus pixel/scan (or pixel_start=pixel_end, scan_start=scan_end).
# Writes a PNG under output_dir: observed/model on left axis; SIF (TOA) and residual on twinx.

using JLD2
using Plots
using TOML

const _SVD_SCRIPT_DIR = @__DIR__
const _PIPE_DIR = dirname(_SVD_SCRIPT_DIR)
const _REPO_ROOT = dirname(_PIPE_DIR)

include(joinpath(_SVD_SCRIPT_DIR, "run_svd_fit.jl"))

function _resolve_pixel_scan(svd::Dict)
    pixel = if haskey(svd, "pixel")
        Int(svd["pixel"])
    elseif haskey(svd, "pixel_start") || haskey(svd, "pixel_end")
        p0 = Int(get(svd, "pixel_start", get(svd, "pixel_end")))
        p1 = Int(get(svd, "pixel_end", p0))
        p0 == p1 || error("single-pixel mode requires pixel_start == pixel_end (got $p0:$p1); or set pixel = N")
        p0
    else
        error("Set [svd_retrieval].pixel or pixel_start/pixel_end for run_single_pixel.jl")
    end
    scan = if haskey(svd, "scan")
        Int(svd["scan"])
    elseif haskey(svd, "scan_start") || haskey(svd, "scan_end")
        s0 = Int(get(svd, "scan_start", get(svd, "scan_end")))
        s1 = Int(get(svd, "scan_end", s0))
        s0 == s1 || error("single-pixel mode requires scan_start == scan_end (got $s0:$s1); or set scan = N")
        s0
    else
        error("Set [svd_retrieval].scan or scan_start/scan_end for run_single_pixel.jl")
    end
    return pixel, scan
end

function _plot_single_pixel_spectra(r; title_str::AbstractString, out_png::AbstractString)
    λ = r.λ
    plt = plot(
        λ, r.y_obs;
        label = "observed",
        xlabel = "Wavelength (nm)",
        ylabel = "Radiance",
        legend = :topright,
        title = title_str,
        linewidth = 2,
        size = (900, 500),
    )
    plot!(plt, λ, r.y_mod; label = "model", linewidth = 2)
    # twin y (Plots.twinx): SIF TOA contribution + residual share the same secondary axis
    plt2 = twinx(plt)
    plot!(plt2, λ, r.sif_toa; label = "SIF (TOA)", linewidth = 1.5, linestyle = :dash, color = :green)
    plot!(plt2, λ, r.residual; label = "residual", linewidth = 1.5, linestyle = :dot, color = :gray)
    ylabel!(plt2, "SIF / residual")
    mkpath(dirname(out_png))
    savefig(plt, out_png)
    return out_png
end

function run_single_pixel(pipeline_path::AbstractString)
    pipe = TOML.parsefile(pipeline_path)
    svd = get(pipe, "svd_retrieval", Dict{String, Any}())
    gid = String(get(svd, "granule_id", ""))
    isempty(gid) && error("[svd_retrieval].granule_id required for run_single_pixel.jl")
    pixel, scan = _resolve_pixel_scan(svd)

    env = _prepare_run_env!(pipe, pipeline_path)
    L1B_path, L2AOP_path, L2BGC_path = _resolve_granule_paths(env.L1B_dir, env.L2AOP_dir, env.L2BGC_dir, gid)
    isfile(L1B_path) || error("L1B not found: $L1B_path")
    isfile(L2AOP_path) || error("L2AOP not found: $L2AOP_path")

    mkpath(env.interim_dir)
    mkpath(env.output_dir)
    interim_path = joinpath(env.interim_dir, "interim_$(gid)_single.nc")
    println("Single-pixel SVD: granule=$gid pixel=$pixel scan=$scan")
    println("  pipeline: ", env.pipeline_path)
    preprocess_and_merge_in_memory(
        L1B_path,
        L2AOP_path,
        isfile(L2BGC_path) ? L2BGC_path : nothing,
        interim_path,
    )

    r = retrieve_svd_single_pixel(interim_path, L1B_path, env.retrieval_cfg, pixel, scan)
    rm(interim_path; force = true)

    stem = "svd_single_$(gid)_p$(pixel)_s$(scan)"
    jld_path = joinpath(env.output_dir, "$(stem).jld2")
    png_path = joinpath(env.output_dir, "$(stem)_spectra.png")
    JLD2.jldsave(
        jld_path;
        λ = r.λ,
        y_obs = r.y_obs,
        y_mod = r.y_mod,
        residual = r.residual,
        sif_toa = r.sif_toa,
        sif_shape = r.sif_shape,
        x_hat = r.x_hat,
        status_code = r.status_code,
        converged = r.converged,
        reduced_chi2 = r.reduced_chi2,
        sif_radiance_678nm = r.sif_radiance_678nm,
        latitude = r.latitude,
        longitude = r.longitude,
        pixel = r.pixel,
        scan = r.scan,
        granule_id = gid,
    )
    title_str = "$(gid)  p=$(pixel) s=$(scan)  SIF@678=$(round(r.sif_radiance_678nm; digits=4))  χ²ᵣ=$(round(r.reduced_chi2; digits=3))  status=$(r.status_code)"
    _plot_single_pixel_spectra(r; title_str = title_str, out_png = png_path)
    println("  status=$(r.status_code) converged=$(r.converged) sif678=$(r.sif_radiance_678nm) redχ²=$(r.reduced_chi2)")
    println("  wrote ", jld_path)
    println("  wrote ", png_path)
    return (; jld_path, png_path, result = r)
end

function main()
    pipeline_path = get(ENV, "PACE_SVD_PIPELINE_CONFIG", normpath(joinpath(_PIPE_DIR, "global_fit_pipeline.toml")))
    if !isempty(ARGS)
        pipeline_path = ARGS[1]
    end
    isfile(pipeline_path) || error("Pipeline config not found: $pipeline_path")
    run_single_pixel(pipeline_path)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
