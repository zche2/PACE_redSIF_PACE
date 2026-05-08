#!/usr/bin/env julia
# CLI for PACE–TROPOMI L1B coincidence matching.
#   julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl [config.toml]
#   MATCH_CONFIG=/path/to/config.toml julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl
#   julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl --fetch-remote [config.toml]
#   julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl --self-test
#   julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl --help
using TOML
using Dates
using Logging
using Printf
using Statistics

include(joinpath(@__DIR__, "MatchCoincidence.jl"))
using .MatchCoincidence
include(joinpath(@__DIR__, "RemoteFetch.jl"))
using .RemoteFetch

function print_usage(io::IO = stdout)
    prog = "toolbox/pace_tropomi_coincidence/run_match.jl"
    print(io, """
Usage:
  julia --project=. $prog [--fetch-remote] [CONFIG.toml]
  julia --project=. $prog --self-test
  julia --project=. $prog --help

Arguments:
  CONFIG.toml       TOML config (paths, match, output, optional remote.*).
                    If omitted: use env MATCH_CONFIG if set, else toolbox/pace_tropomi_coincidence/config.toml (beside run_match.jl).

Flags:
  --fetch-remote    Download enabled [remote.tropomi] / [remote.pace] granules, then match.
  --self-test       Run built-in tests only (no config file).
  --help, -h        Show this text.

Environment:
  MATCH_CONFIG      When non-empty, config path (overrides positional CONFIG.toml).

More detail: toolbox/pace_tropomi_coincidence/README.md and config.example.toml
""")
end

function _sym_band(s::AbstractString)::Symbol
    u = uppercase(strip(s))
    u == "BD5" && return :BD5
    u == "BD6" && return :BD6
    error("tropomi_band must be BD5 or BD6, got $(repr(s))")
end

function _missing_toml_window_key(x)::Bool
    x === nothing && return true
    isa(x, AbstractString) && isempty(strip(x)) && return true
    return false
end

function _sym_multi(s::AbstractString)::Symbol
    x = lowercase(strip(s))
    x == "skip" && return :skip
    x == "first" && return :first
    x == "error" && return :error
    error("multi_pace_timestamp must be skip|first|error, got $(repr(s))")
end

function load_config(path::AbstractString)::Tuple{
    MatchConfig,
    String,
    Union{Nothing,String},
    Union{Nothing,RemoteTropomiGESDISCConfig},
    Union{Nothing,RemotePaceEarthAccessConfig},
}
    cfg_toml = TOML.parsefile(path)
    paths = get(cfg_toml, "paths", Dict{String,Any}())
    match_c = get(cfg_toml, "match", Dict{String,Any}())
    out_c = get(cfg_toml, "output", Dict{String,Any}())

    pace_dirs = String[get(paths, "pace_dirs", String[])...]
    tropomi_dirs = String[get(paths, "tropomi_dirs", String[])...]
    pace_manifest = get(paths, "pace_manifest", nothing)
    tropomi_manifest = get(paths, "tropomi_manifest", nothing)
    pace_manifest = pace_manifest === nothing || isempty(string(pace_manifest)) ? nothing : String(pace_manifest)
    tropomi_manifest = tropomi_manifest === nothing || isempty(string(tropomi_manifest)) ? nothing : String(tropomi_manifest)

    ws_raw = get(match_c, "utc_window_start", nothing)
    we_raw = get(match_c, "utc_window_end", nothing)
    ws_miss = _missing_toml_window_key(ws_raw)
    we_miss = _missing_toml_window_key(we_raw)
    ws_miss != we_miss && error("[match] utc_window_start and utc_window_end must both be set or both omitted")
    utc_ws = nothing
    utc_we = nothing
    if !ws_miss
        utc_ws = parse_match_window_endpoint(ws_raw, false)
        utc_we = parse_match_window_endpoint(we_raw, true)
    end
    stride = Int(haskey(match_c, "utc_day_stride") ? match_c["utc_day_stride"] : get(match_c, "utc_interval_days", 1))

    cfg = MatchConfig(;
        pace_dirs,
        tropomi_dirs,
        pace_manifest,
        tropomi_manifest,
        pace_glob = String(get(match_c, "pace_glob", "**/PACE_OCI.*.L1B.V3.nc")),
        tropomi_band = _sym_band(String(get(match_c, "tropomi_band", "BD6"))),
        strict_pace_parse = Bool(get(match_c, "strict_pace_parse", false)),
        multi_pace_timestamp = _sym_multi(String(get(match_c, "multi_pace_timestamp", "skip"))),
        inclusive_endpoints = Bool(get(match_c, "inclusive_endpoints", true)),
        neighbor_days = Int(get(match_c, "neighbor_days", 0)),
        use_netcdf_time_for_pace = Bool(get(match_c, "use_netcdf_time_for_pace", false)),
        pace_basename_strip = String[get(match_c, "pace_basename_strip", String[])...],
        utc_window_start = utc_ws,
        utc_window_end = utc_we,
        utc_day_stride = stride,
    )

    csv_path = String(get(out_c, "csv_path", "coincidences.csv"))
    jsonl_raw = get(out_c, "jsonl_path", nothing)
    jsonl_path = jsonl_raw === nothing || isempty(string(jsonl_raw)) ? nothing : String(jsonl_raw)

    rt, rp = load_remote_configs(cfg_toml, cfg.tropomi_band)
    return cfg, csv_path, jsonl_path, rt, rp
end

function parse_fetch_remote(argv)::Tuple{Bool,Vector{String}}
    use_fetch = false
    rest = String[]
    for a in argv
        if a == "--fetch-remote"
            use_fetch = true
        else
            push!(rest, a)
        end
    end
    return use_fetch, rest
end

function main()
    argv = ARGS
    if any(x -> x == "--help" || x == "-h", argv)
        print_usage()
        exit(0)
    end
    if any(==("--self-test"), argv)
        ok = run_self_tests()
        ok &= RemoteFetch.run_remote_self_tests()
        exit(ok ? 0 : 1)
    end

    fetch_remote, argv_rest = parse_fetch_remote(argv)

    config_path = get(ENV, "MATCH_CONFIG", "")
    if isempty(config_path)
        config_path = length(argv_rest) >= 1 ? argv_rest[1] : joinpath(@__DIR__, "config.toml")
    end
    isfile(config_path) || error("Config not found: $config_path")

    cfg, csv_path, jsonl_path, rt, rp = load_config(config_path)
    println("PACE–TROPOMI coincidence matching")
    println("  config: ", config_path)
    println("  tropomi_band: ", cfg.tropomi_band)
    println("  neighbor_days: ", cfg.neighbor_days)
    if cfg.utc_window_start !== nothing
        println("  utc_window: ", cfg.utc_window_start, " … ", cfg.utc_window_end,
                " (utc_day_stride=", cfg.utc_day_stride, ")")
    end

    if fetch_remote
        active_t = rt !== nothing && rt.enabled
        active_p = rp !== nothing && rp.enabled
        if !active_t && !active_p
            @warn "--fetch-remote was set but neither [remote.tropomi] nor [remote.pace] has enabled=true"
        end
        fetch_remote_if_requested!(cfg, rt, rp; do_fetch = true)
    end

    rows = find_coincidences(cfg)
    isdir(dirname(abspath(csv_path))) || mkpath(dirname(abspath(csv_path)))
    write_coincidences_csv(csv_path, rows)
    println("  wrote CSV: ", csv_path, " (", length(rows), " pairs)")

    if jsonl_path !== nothing
        isdir(dirname(abspath(jsonl_path))) || mkpath(dirname(abspath(jsonl_path)))
        write_coincidences_jsonl(jsonl_path, rows)
        println("  wrote JSONL: ", jsonl_path)
    end

    if !isempty(rows)
        margins = [r.containment_margin_s for r in rows]
        @printf("  containment_margin_s: min %.4f  median %.4f  max %.4f\n",
                minimum(margins), median(Float64.(margins)), maximum(margins))
    end
end

main()
