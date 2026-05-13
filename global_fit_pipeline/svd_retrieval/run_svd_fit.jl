#!/usr/bin/env julia
# Multi-granule SVD retrieval: read global_fit_pipeline.toml (single file: [general], [svd_retrieval],
# plus [data], [spectral], [fit], …). No demo_example paths.
#
# Usage:
#   julia -t 8 global_fit_pipeline/svd_retrieval/run_svd_fit.jl [path/to/global_fit_pipeline.toml]

using Base.Threads
using Dates
using Glob
using TOML

const _SVD_SCRIPT_DIR = @__DIR__
const _PIPE_DIR = dirname(_SVD_SCRIPT_DIR)
const _REPO_ROOT = dirname(_PIPE_DIR)

include(joinpath(_SVD_SCRIPT_DIR, "svd_granule.jl"))

"""Top-level TOML sections consumed only for orchestration (stripped before the retrieval dict)."""
const _PIPELINE_ORCHESTRATION_ROOTS = ("download", "general", "svd_retrieval")

function _merge_config(cfg_base::Dict, cfg_overlay::Dict)
    out = copy(cfg_base)
    for (k, v) in cfg_overlay
        if haskey(out, k) && out[k] isa Dict && v isa Dict
            out[k] = _merge_config(out[k], v)
        else
            out[k] = v
        end
    end
    return out
end

function _resolve_repo_path(p::AbstractString)
    isabspath(p) ? String(p) : normpath(joinpath(_REPO_ROOT, p))
end

"""Everything except [download], [general], [svd_retrieval] → in-memory retrieval dict for the MWE stack."""
function _retrieval_cfg_from_pipeline(pipe::Dict)
    out = Dict{String, Any}()
    for (k, v) in pipe
        k in _PIPELINE_ORCHESTRATION_ROOTS && continue
        out[k] = v
    end
    out
end

"""In-memory retrieval config: same overrides as before (no `merged_svd_retrieval.toml` on disk)."""
function _effective_retrieval_cfg_dict(
    retrieval_cfg::Dict{String, Any},
    output_dir::AbstractString,
    parallel_pixels::Bool,
)
    # Per-granule pixel LM threading (`julia -t N`, batch_fit.use_threads). Granule-level `@threads` is not used
    # (concurrent NetCDF/HDF5 in one process corrupts the C library).
    pixel_threads = parallel_pixels
    return _merge_config(
        retrieval_cfg,
        Dict{String, Any}(
            "batch_fit" => merge(
                get(retrieval_cfg, "batch_fit", Dict{String, Any}()),
                Dict{String, Any}(
                    "output_dir" => output_dir,
                    "use_threads" => pixel_threads,
                ),
            ),
            "pace_observation" => merge(
                get(retrieval_cfg, "pace_observation", Dict{String, Any}()),
                Dict{String, Any}("spectrum_var" => "Rtoa_red", "wavelength_var" => "red_wavelength"),
            ),
        ),
    )
end

function _pixel_scan_kw(svd::Dict)
    has_ps = haskey(svd, "pixel_start") || haskey(svd, "pixel_end")
    has_ss = haskey(svd, "scan_start") || haskey(svd, "scan_end")
    pixel_range = if has_ps
        p0 = Int(get(svd, "pixel_start", 1))
        p1 = Int(get(svd, "pixel_end", p0))
        p0:p1
    else
        nothing
    end
    scan_range = if has_ss
        s0 = Int(get(svd, "scan_start", 1))
        s1 = Int(get(svd, "scan_end", s0))
        s0:s1
    else
        nothing
    end
    return (; pixel_range, scan_range)
end

function _extract_granule_id(L1B_basename::AbstractString)
    stem = splitext(L1B_basename)[1]
    m = match(r"^PACE_OCI\.(.+)\.L1B\.V3$", stem)
    return m !== nothing ? String(m.captures[1]) : stem
end

function _resolve_granule_paths(
    L1B_dir::AbstractString,
    L2AOP_dir::AbstractString,
    L2BGC_dir::AbstractString,
    granule_id::AbstractString,
)
    L1B_path = joinpath(L1B_dir, "PACE_OCI.$(granule_id).L1B.V3.nc")
    L2AOP_path = joinpath(L2AOP_dir, "PACE_OCI.$(granule_id).L2.OC_AOP.V3_1.nc")
    L2BGC_path = joinpath(L2BGC_dir, "PACE_OCI.$(granule_id).L2.OC_BGC.V3_1.nc")
    return L1B_path, L2AOP_path, L2BGC_path
end

function _parse_yyyymmdd(s::AbstractString)
    str = String(s)
    m = match(r"^(\d{4})(\d{2})(\d{2})$", str)
    m === nothing && error("Invalid calendar day (expect YYYYMMDD): $(repr(str))")
    y, mo, d = parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3])
    return Date(y, mo, d)
end

"""Inclusive calendar days from start to end, stepping by `step_days` (>= 1)."""
function _expand_dates_range(start_s::AbstractString, end_s::AbstractString, step_days::Int)
    step_days >= 1 || error("[svd_retrieval].date_interval_days must be >= 1, got $step_days")
    d0 = _parse_yyyymmdd(start_s)
    d1 = _parse_yyyymmdd(end_s)
    d1 >= d0 || error("[svd_retrieval].date_end must be on or after date_start ($start_s .. $end_s)")
    out = String[]
    d = d0
    inc = Day(step_days)
    while d <= d1
        push!(out, Dates.format(d, dateformat"yyyymmdd"))
        d += inc
    end
    return out
end

"""
Resolve `mode=dates` day list: prefer non-empty `dates = [...]`; else `date_start` + `date_end`
with optional `date_interval_days` (default 1).
"""
function _resolve_dates_mode_list(svd::Dict)
    dates_raw = get(svd, "dates", nothing)
    if dates_raw !== nothing && isa(dates_raw, AbstractVector) && !isempty(dates_raw)
        return [String(d) for d in dates_raw]
    end
    ds = String(get(svd, "date_start", ""))
    de = String(get(svd, "date_end", ""))
    if !isempty(ds) && !isempty(de)
        step = Int(get(svd, "date_interval_days", 1))
        return _expand_dates_range(ds, de, step)
    end
    error(
        "mode=dates requires [svd_retrieval].dates = [...] or date_start and date_end " *
            "(optional date_interval_days, default 1 calendar day)",
    )
end

function _prepare_run_env!(pipe::Dict, pipeline_path::AbstractString)
    gen = get(pipe, "general", Dict{String, Any}())
    gen isa Dict || error("[general] must be a table")
    svd = get(pipe, "svd_retrieval", Dict{String, Any}())
    isempty(svd) && error("[svd_retrieval] table missing")

    L1B_dir = _resolve_repo_path(String(get(gen, "L1B_dir", "")))
    L2AOP_dir = _resolve_repo_path(String(get(gen, "L2AOP_dir", "")))
    L2BGC_dir = _resolve_repo_path(String(get(gen, "L2BGC_dir", "")))
    isempty(L1B_dir) && error("[general].L1B_dir is required")
    isempty(L2AOP_dir) && error("[general].L2AOP_dir is required")
    isempty(L2BGC_dir) && error("[general].L2BGC_dir is required")

    output_dir = _resolve_repo_path(String(get(svd, "output_dir", joinpath(_REPO_ROOT, "svd_retrieval_output"))))
    interim_dir = _resolve_repo_path(String(get(svd, "interim_dir", joinpath(_REPO_ROOT, "svd_retrieval_interim"))))
    parallel_granules = Bool(get(svd, "parallel_granules", false))
    parallel_pixels = Bool(get(svd, "parallel_pixels", true))
    use_in_memory_merge = Bool(get(svd, "use_in_memory_merge", true))

    retrieval_cfg = _retrieval_cfg_from_pipeline(pipe)
    isempty(retrieval_cfg) && error(
        "Pipeline TOML must define retrieval sections (e.g. [data], [spectral], [kernel], [fit], [batch_fit]) " *
            "besides [general], [download], [svd_retrieval]",
    )

    retrieval_eff = _effective_retrieval_cfg_dict(retrieval_cfg, output_dir, parallel_pixels)
    bf = get(retrieval_eff, "batch_fit", Dict{String, Any}())
    use_gpu_merged = haskey(svd, "use_gpu") ? Bool(svd["use_gpu"]) : Bool(get(bf, "use_gpu", false))
    gpu_tile_merged =
        haskey(svd, "gpu_tile_pixels") ? Int(svd["gpu_tile_pixels"]) : Int(get(bf, "gpu_tile_pixels", 256))
    retrieval_eff["batch_fit"] = merge(
        bf,
        Dict{String, Any}(
            "use_gpu" => use_gpu_merged,
            "gpu_tile_pixels" => max(1, gpu_tile_merged),
        ),
    )
    mkpath(output_dir)
    mkpath(interim_dir)
    if parallel_granules && Threads.nthreads() > 1
        @warn "parallel_granules=true does not enable multi-threaded granules (unsafe with NetCDF/HDF5). Granules run one at a time; use parallel_pixels=true and `julia -t N`, or launch separate Julia processes per granule for throughput."
    end
    pipeline_abspath = abspath(String(pipeline_path))
    return (;
        L1B_dir,
        L2AOP_dir,
        L2BGC_dir,
        output_dir,
        interim_dir,
        retrieval_cfg = retrieval_eff,
        pipeline_path = pipeline_abspath,
        use_in_memory_merge,
        parallel_granules,
    )
end

function run_svd_global_fit_date(pipeline_path::AbstractString; date_override::Union{Nothing, AbstractString} = nothing)
    pipe = TOML.parsefile(pipeline_path)
    svd = get(pipe, "svd_retrieval", Dict{String, Any}())
    if date_override !== nothing
        svd = merge(svd, Dict{String, Any}("date" => String(date_override)))
    end
    date = String(get(svd, "date", ""))
    isempty(date) && error("[svd_retrieval].date required when mode=date (or pass date_override=)")
    pipe_run = date_override === nothing ? pipe : merge(pipe, Dict{String, Any}("svd_retrieval" => svd))
    env = _prepare_run_env!(pipe_run, pipeline_path)
    pr = _pixel_scan_kw(svd)
    pattern = "PACE_OCI.$(date)T*.L1B.V3.nc"
    L1B_files = Glob.glob(pattern, env.L1B_dir)
    sort!(L1B_files)
    if isempty(L1B_files)
        @warn "No L1B files for date $date in $(env.L1B_dir) (pattern $pattern)"
        return String[]
    end
    granules = [(f, _extract_granule_id(basename(f))) for f in L1B_files]
    println(
        "SVD global fit (date=$date): $(length(granules)) granule(s); granule concurrency=sequential (NetCDF/HDF5-safe)",
    )
    println("  pipeline config: ", env.pipeline_path)
    results_per_k = Vector{Union{Nothing,String}}(nothing, length(granules))
    for k in eachindex(granules)
        L1B_path, gid = granules[k]
        L2AOP_path, L2BGC_path = _resolve_granule_paths(env.L1B_dir, env.L2AOP_dir, env.L2BGC_dir, gid)
        isfile(L1B_path) || (@warn "L1B missing: $L1B_path"; continue)
        isfile(L2AOP_path) || (@warn "L2AOP missing: $L2AOP_path"; continue)
        run_svd_granule(
            L1B_path,
            L2AOP_path,
            isfile(L2BGC_path) ? L2BGC_path : nothing,
            env.output_dir,
            env.interim_dir,
            env.retrieval_cfg,
            env.pipeline_path;
            use_in_memory_merge = env.use_in_memory_merge,
            pixel_range = pr.pixel_range,
            scan_range = pr.scan_range,
        )
        results_per_k[k] = L1B_path
    end
    return String[r for r in results_per_k if r !== nothing]
end

function run_svd_global_fit_dates(pipeline_path::AbstractString)
    cfg = TOML.parsefile(pipeline_path)
    svd = get(cfg, "svd_retrieval", Dict{String, Any}())
    dates = _resolve_dates_mode_list(svd)
    println("SVD global fit (dates mode): $(length(dates)) calendar day(s): ", join(dates, ", "))
    all_results = String[]
    for d in dates
        append!(all_results, run_svd_global_fit_date(pipeline_path; date_override = d))
    end
    return all_results
end

function run_svd_global_fit_granule(pipeline_path::AbstractString, granule_id::AbstractString)
    pipe = TOML.parsefile(pipeline_path)
    env = _prepare_run_env!(pipe, pipeline_path)
    svd = get(pipe, "svd_retrieval", Dict{String, Any}())
    pr = _pixel_scan_kw(svd)
    L1B_path, L2AOP_path, L2BGC_path = _resolve_granule_paths(env.L1B_dir, env.L2AOP_dir, env.L2BGC_dir, granule_id)
    isfile(L1B_path) || error("L1B not found: $L1B_path")
    isfile(L2AOP_path) || error("L2AOP not found: $L2AOP_path")
    println("SVD global fit (granule=$granule_id)")
    println("  pipeline config: ", env.pipeline_path)
    run_svd_granule(
        L1B_path,
        L2AOP_path,
        isfile(L2BGC_path) ? L2BGC_path : nothing,
        env.output_dir,
        env.interim_dir,
        env.retrieval_cfg,
        env.pipeline_path;
        use_in_memory_merge = env.use_in_memory_merge,
        pixel_range = pr.pixel_range,
        scan_range = pr.scan_range,
    )
    return [L1B_path]
end

function main()
    pipeline_path = get(ENV, "PACE_SVD_PIPELINE_CONFIG", normpath(joinpath(_PIPE_DIR, "global_fit_pipeline.toml")))
    if !isempty(ARGS)
        pipeline_path = ARGS[1]
    end
    isfile(pipeline_path) || error("Pipeline config not found: $pipeline_path")
    cfg = TOML.parsefile(pipeline_path)
    svd = get(cfg, "svd_retrieval", Dict{String, Any}())
    mode = String(get(svd, "mode", "date"))
    if mode == "date"
        run_svd_global_fit_date(pipeline_path)
    elseif mode == "dates"
        run_svd_global_fit_dates(pipeline_path)
    elseif mode == "granule"
        gid = String(get(svd, "granule_id", ""))
        isempty(gid) && error("[svd_retrieval].granule_id required when mode=granule")
        run_svd_global_fit_granule(pipeline_path, gid)
    else
        error("Unknown [svd_retrieval].mode=$(repr(mode)); use date, dates, or granule")
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
