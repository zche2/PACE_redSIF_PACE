#!/usr/bin/env julia
# Global fit pipeline: raw PACE L1B+L2 → TOA radiance → batch retrieval.
# Modes: date (all granules for a date) or granule (single granule by ID).
# Parallel: granule-level (Threads.@threads) and pixel-level (via batch_fit).
#
# Usage: julia -t N demo_example/global_fit/Run_global_fit.jl [config_path]
#        (N = number of threads)

using Base.Threads
using TOML
using Glob

const _GLOBAL_FIT_DIR = @__DIR__
const _DEMO_DIR = joinpath(_GLOBAL_FIT_DIR, "..")

include(joinpath(_GLOBAL_FIT_DIR, "granule_retrieval.jl"))

# ----- Granule ID and path resolution -----
function _extract_granule_id(L1B_basename::AbstractString)
    stem = splitext(L1B_basename)[1]
    m = match(r"^PACE_OCI\.(.+)\.L1B\.V3$", stem)
    return m !== nothing ? m.captures[1] : stem
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

# ----- Mode 1: Date-based (all granules for a date) -----
function run_global_fit_date(config_path::AbstractString)
    base_config = normpath(joinpath(_DEMO_DIR, "Simple_PACE_xSecFit_MWE_zcheVer.toml"))
    cfg_base = isfile(base_config) ? TOML.parsefile(base_config) : Dict{String, Any}()
    cfg_overlay = TOML.parsefile(config_path)
    cfg = _merge_config(cfg_base, cfg_overlay)
    gf = get(cfg, "global_fit", Dict{String, Any}())
    date = String(get(gf, "date", ""))
    isempty(date) && error("[global_fit].date required when mode=date")
    L1B_dir = String(get(gf, "L1B_dir", ""))
    L2AOP_dir = String(get(gf, "L2AOP_dir", ""))
    L2BGC_dir = String(get(gf, "L2BGC_dir", ""))
    output_dir = String(get(gf, "output_dir", joinpath(_DEMO_DIR, "global_fit_output")))
    interim_dir = String(get(gf, "interim_dir", joinpath(_DEMO_DIR, "global_fit_interim")))
    parallel_granules = Bool(get(gf, "parallel_granules", true))
    use_in_memory_merge = Bool(get(gf, "use_in_memory_merge", !parallel_granules))

    L1B_dir = isabspath(L1B_dir) ? L1B_dir : joinpath(_DEMO_DIR, L1B_dir)
    L2AOP_dir = isabspath(L2AOP_dir) ? L2AOP_dir : joinpath(_DEMO_DIR, L2AOP_dir)
    L2BGC_dir = isabspath(L2BGC_dir) ? L2BGC_dir : joinpath(_DEMO_DIR, L2BGC_dir)
    output_dir = isabspath(output_dir) ? output_dir : joinpath(_DEMO_DIR, output_dir)
    interim_dir = isabspath(interim_dir) ? interim_dir : joinpath(_DEMO_DIR, interim_dir)

    mkpath(interim_dir)
    cfg["batch_fit"] = merge(get(cfg, "batch_fit", Dict()), Dict("output_dir" => output_dir, "use_threads" => batch_use_threads))
    cfg["pace_observation"] = merge(get(cfg, "pace_observation", Dict()), Dict("spectrum_var" => "Rtoa_red", "wavelength_var" => "red_wavelength"))
    merged_config = joinpath(interim_dir, "merged_config.toml")
    open(merged_config, "w") do io
        TOML.print(io, cfg)
    end

    pattern = "PACE_OCI.$(date)T*.L1B.V3.nc"
    L1B_files = Glob.glob(pattern, L1B_dir)
    sort!(L1B_files)

    if isempty(L1B_files)
        @warn "No L1B files found for date $date in $L1B_dir (pattern: $pattern)"
        return String[]
    end

    granules = [(f, _extract_granule_id(basename(f))) for f in L1B_files]
    n_threads = parallel_granules ? nthreads() : 1
    # When parallel_granules: use threads for granules, disable pixel-level threads to avoid oversubscription
    batch_use_threads = !parallel_granules && Bool(get(get(cfg, "batch_fit", Dict()), "use_threads", true))

    println("Global fit (date mode): $date")
    println("  granules: ", length(granules))
    println("  parallel_granules: ", parallel_granules, " (threads: ", n_threads, ")")
    println("  batch use_threads (pixel-level): ", batch_use_threads)
    println("  output_dir: ", output_dir)

    results = String[]
    if parallel_granules && n_threads > 1
        Threads.@threads for i in eachindex(granules)
            L1B_path, granule_id = granules[i]
            L2AOP_path, L2BGC_path = _resolve_granule_paths(L1B_dir, L2AOP_dir, L2BGC_dir, granule_id)
            isfile(L1B_path) || (@warn "L1B missing: $L1B_path"; continue)
            isfile(L2AOP_path) || (@warn "L2AOP missing: $L2AOP_path"; continue)
            run_granule_retrieval(
                L1B_path,
                L2AOP_path,
                isfile(L2BGC_path) ? L2BGC_path : nothing,
                output_dir,
                interim_dir,
                merged_config;
                use_threads=false,
                keep_interim=false,
                shared_config=true,
                use_in_memory_merge=use_in_memory_merge,
            )
            push!(results, L1B_path)
        end
    else
        for (L1B_path, granule_id) in granules
            L2AOP_path, L2BGC_path = _resolve_granule_paths(L1B_dir, L2AOP_dir, L2BGC_dir, granule_id)
            isfile(L1B_path) || (@warn "L1B missing: $L1B_path"; continue)
            isfile(L2AOP_path) || (@warn "L2AOP missing: $L2AOP_path"; continue)
            run_granule_retrieval(
                L1B_path,
                L2AOP_path,
                isfile(L2BGC_path) ? L2BGC_path : nothing,
                output_dir,
                interim_dir,
                merged_config;
                use_threads=batch_use_threads,
                keep_interim=false,
                shared_config=true,
                use_in_memory_merge=use_in_memory_merge,
            )
            push!(results, L1B_path)
        end
    end

    return results
end

# ----- Mode 2: Single granule by ID -----
function run_global_fit_granule(config_path::AbstractString, granule_id::AbstractString)
    base_config = normpath(joinpath(_DEMO_DIR, "Simple_PACE_xSecFit_MWE_zcheVer.toml"))
    cfg_base = isfile(base_config) ? TOML.parsefile(base_config) : Dict{String, Any}()
    cfg_overlay = TOML.parsefile(config_path)
    cfg = _merge_config(cfg_base, cfg_overlay)
    gf = get(cfg, "global_fit", Dict{String, Any}())
    L1B_dir = String(get(gf, "L1B_dir", ""))
    L2AOP_dir = String(get(gf, "L2AOP_dir", ""))
    L2BGC_dir = String(get(gf, "L2BGC_dir", ""))
    output_dir = String(get(gf, "output_dir", joinpath(_DEMO_DIR, "global_fit_output")))
    interim_dir = String(get(gf, "interim_dir", joinpath(_DEMO_DIR, "global_fit_interim")))

    L1B_dir = isabspath(L1B_dir) ? L1B_dir : joinpath(_DEMO_DIR, L1B_dir)
    L2AOP_dir = isabspath(L2AOP_dir) ? L2AOP_dir : joinpath(_DEMO_DIR, L2AOP_dir)
    L2BGC_dir = isabspath(L2BGC_dir) ? L2BGC_dir : joinpath(_DEMO_DIR, L2BGC_dir)
    output_dir = isabspath(output_dir) ? output_dir : joinpath(_DEMO_DIR, output_dir)
    interim_dir = isabspath(interim_dir) ? interim_dir : joinpath(_DEMO_DIR, interim_dir)

    use_in_memory_merge = Bool(get(gf, "use_in_memory_merge", true))

    mkpath(interim_dir)
    cfg["batch_fit"] = merge(get(cfg, "batch_fit", Dict()), Dict("output_dir" => output_dir, "use_threads" => true))
    cfg["pace_observation"] = merge(get(cfg, "pace_observation", Dict()), Dict("spectrum_var" => "Rtoa_red", "wavelength_var" => "red_wavelength"))
    merged_config = joinpath(interim_dir, "merged_config_granule.toml")
    open(merged_config, "w") do io
        TOML.print(io, cfg)
    end

    L1B_path, L2AOP_path, L2BGC_path = _resolve_granule_paths(L1B_dir, L2AOP_dir, L2BGC_dir, granule_id)
    isfile(L1B_path) || error("L1B not found: $L1B_path")
    isfile(L2AOP_path) || error("L2AOP not found: $L2AOP_path")

    println("Global fit (granule mode): $granule_id")
    run_granule_retrieval(
        L1B_path,
        L2AOP_path,
        isfile(L2BGC_path) ? L2BGC_path : nothing,
        output_dir,
        interim_dir,
        merged_config;
        use_threads=true,
        keep_interim=false,
        shared_config=true,
        use_in_memory_merge=use_in_memory_merge,
    )
    return [L1B_path]
end

# Deep merge: cfg_overlay overrides cfg_base
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

# ----- Main entry -----
function main()
    config_path = get(
        ENV,
        "PACE_GLOBAL_FIT_CONFIG",
        normpath(joinpath(_DEMO_DIR, "global_fit", "global_fit_config.toml")),
    )
    if !isempty(ARGS)
        config_path = ARGS[1]
    end
    isfile(config_path) || error("Config not found: $config_path")

    cfg = TOML.parsefile(config_path)
    gf = get(cfg, "global_fit", Dict{String, Any}())
    mode = String(get(gf, "mode", "date"))
    granule_id = String(get(gf, "granule_id", ""))

    if mode == "date"
        run_global_fit_date(config_path)
    elseif mode == "granule"
        isempty(granule_id) && error("[global_fit].granule_id required when mode=granule")
        run_global_fit_granule(config_path, granule_id)
    else
        error("Unknown mode: $mode (use 'date' or 'granule')")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
