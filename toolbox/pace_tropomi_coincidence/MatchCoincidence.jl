"""
PACE–TROPOMI L1B coincidence matching: parse granule times from basenames,
bucket by UTC date, pair when PACE granule instant lies in TROPOMI [T_start, T_end].
"""
module MatchCoincidence

using Dates
using Glob
using JSON
using NCDatasets
using Printf

export MatchConfig,
    parse_match_window_endpoint,
    parse_tropomi_granule,
    parse_pace_granule,
    collect_local_paths,
    read_manifest_paths,
    find_coincidences,
    write_coincidences_csv,
    write_coincidences_jsonl,
    containment_margin_seconds,
    run_self_tests

# --- Config --------------------------------------------------------------------

Base.@kwdef struct MatchConfig
    pace_dirs::Vector{String}       = String[]
    tropomi_dirs::Vector{String}    = String[]
    pace_manifest::Union{Nothing,String}    = nothing
    tropomi_manifest::Union{Nothing,String} = nothing
    pace_glob::String               = "**/PACE_OCI.*.L1B.V3.nc"
    tropomi_band::Symbol            = :BD6   # :BD5 or :BD6
    strict_pace_parse::Bool         = false
    """When basename has multiple `YYYYMMDDTHHMMSS`: `:skip`, `:first`, `:error`"""
    multi_pace_timestamp::Symbol    = :skip
    inclusive_endpoints::Bool        = true
    """Also compare PACE against TROPOMI whose T_start falls on adjacent UTC days (±k)."""
    neighbor_days::Int             = 0
    use_netcdf_time_for_pace::Bool  = false
    """Strip these substrings from basename before PACE parse (optional)."""
    pace_basename_strip::Vector{String} = String[]
    """Inclusive UTC filter on granule times; both must be set or both `nothing` (no filter)."""
    utc_window_start::Union{Nothing,DateTime} = nothing
    utc_window_end::Union{Nothing,DateTime} = nothing
    """Keep only granules whose anchor UTC date (`Date(t_pace)` / TROPOMI `Date(t_start)`) satisfies
    `mod(day_index_from_window_start, utc_day_stride)==0` with `utc_day_stride ≥ 1` (default 1 = every day)."""
    utc_day_stride::Int = 1
end

# --- Regex ---------------------------------------------------------------------

const RE_TROPOMI_L1B_RA = r"L1B_RA_BD([56])_(\d{8}T\d{6})_(\d{8}T\d{6})_(\d+)_"
"""Canonical PACE stem fragment; allow extra characters after `.V3` (renamed files)."""
const RE_PACE_CANONICAL = r"PACE_OCI\.(\d{8}T\d{6})\.L1B\.V3"
const RE_GRANULE_TIME = r"\d{8}T\d{6}"

function _parse_yyyymmddTHHMMSS(s::AbstractString)::DateTime
    DateTime(
        parse(Int, s[1:4]),
        parse(Int, s[5:6]),
        parse(Int, s[7:8]),
        parse(Int, s[10:11]),
        parse(Int, s[12:13]),
        parse(Int, s[14:15]),
    )
end

"""Parse TROPOMI L1B RA basename; returns (band_char, t_start, t_end, orbit_str) or nothing."""
function parse_tropomi_granule(path::AbstractString)
    base = basename(path)
    m = match(RE_TROPOMI_L1B_RA, base)
    m === nothing && return nothing
    band = m.captures[1]
    t0 = _parse_yyyymmddTHHMMSS(m.captures[2])
    t1 = _parse_yyyymmddTHHMMSS(m.captures[3])
    orbit = m.captures[4]
    if t1 < t0
        t0, t1 = t1, t0
    end
    return (; band, t_start = t0, t_end = t1, orbit)
end

function _strip_basename_tags(base::AbstractString, strips::Vector{String})
    s = base
    for tag in strips
        isempty(tag) && continue
        s = replace(s, tag => "")
    end
    return s
end

"""
Parse PACE granule UTC instant from path basename.
Returns `(t_pace, parse_mode)` where `parse_mode` is `:canonical`, `:fallback_substring`, or `:netcdf`,
or `(nothing, :failed)` / `(nothing, :skipped_multi)`.
"""
function parse_pace_granule(path::AbstractString, cfg::MatchConfig)
    raw_base = basename(path)
    stem, _ext = splitext(raw_base)
    stem = _strip_basename_tags(stem, cfg.pace_basename_strip)

    mc = match(RE_PACE_CANONICAL, stem)
    if mc !== nothing
        return (_parse_yyyymmddTHHMMSS(mc.captures[1]), :canonical)
    end
    cfg.strict_pace_parse && return (nothing, :failed)

    tokens = [m.match for m in eachmatch(Regex(RE_GRANULE_TIME.pattern), stem)]
    if isempty(tokens)
        if cfg.use_netcdf_time_for_pace && isfile(path)
            tnc = _try_pace_time_netcdf(path)
            tnc !== nothing && return (tnc, :netcdf)
        end
        return (nothing, :failed)
    elseif length(tokens) == 1
        return (_parse_yyyymmddTHHMMSS(tokens[1]), :fallback_substring)
    else
        if cfg.multi_pace_timestamp === :first
            return (_parse_yyyymmddTHHMMSS(tokens[1]), :fallback_substring)
        elseif cfg.multi_pace_timestamp === :error
            error("PACE basename has multiple granule-like timestamps: $raw_base")
        else
            return (nothing, :skipped_multi)
        end
    end
end

function _try_pace_time_netcdf(path::AbstractString)::Union{Nothing,DateTime}
    try
        ds = Dataset(path)
        try
            for attr in ("time_coverage_start", "time_coverage_begin")
                if haskey(ds.attrib, attr)
                    return parse_iso8601_safe(string(ds.attrib[attr]))
                end
            end
        finally
            close(ds)
        end
    catch
        return nothing
    end
    return nothing
end

"""Parse TOML `[match] utc_window_*`: datetime string, or date-only `yyyy-mm-dd` (start → 00:00:00, end → 23:59:59 UTC)."""
function parse_match_window_endpoint(x, is_end::Bool)::DateTime
    x isa DateTime && return x
    x isa Date && return is_end ? DateTime(x) + Day(1) - Second(1) : DateTime(x)
    s = strip(string(x))
    isempty(s) && error("utc_window: empty value")
    if occursin('T', s) || occursin('t', s)
        return parse_iso8601_safe(replace(s, 't' => 'T'))
    end
    d = Date(s)
    return is_end ? DateTime(d) + Day(1) - Second(1) : DateTime(d)
end

"""Lightweight ISO8601-ish parse for typical NASA attributes."""
function parse_iso8601_safe(s::AbstractString)::DateTime
    s = strip(s)
    if endswith(s, 'Z')
        s = s[1:prevind(s, lastindex(s))]
    end
    if occursin('T', s)
        date_part, time_part = split(s, 'T'; limit = 2)
        y, mo, d = parse.(Int, split(date_part, '-'))
        rest = split(time_part, '.'; limit = 2)[1]
        hms = split(rest, ':'; limit = 3)
        h = parse(Int, hms[1])
        mi = length(hms) >= 2 ? parse(Int, hms[2]) : 0
        sec = length(hms) >= 3 ? parse(Int, hms[3]) : 0
        return DateTime(y, mo, d, h, mi, sec)
    end
    error("Unrecognized time string: $s")
end

"""Minimum distance from `t` to nearest endpoint of [lo, hi] in seconds (≥0 if inside)."""
function containment_margin_seconds(
    t::DateTime,
    lo::DateTime,
    hi::DateTime;
    inclusive::Bool = true,
)::Float64
    a = min(lo, hi)
    b = max(lo, hi)
    dt_lo = float_seconds(t - a)
    dt_hi = float_seconds(b - t)
    inside = if inclusive
        (t >= a && t <= b)
    else
        (t > a && t < b)
    end
    if !inside
        return -min(abs(dt_lo), abs(dt_hi))
    end
    return min(dt_lo, dt_hi)
end

float_seconds(x::Millisecond) = x.value / 1000.0

function _intervals_overlap_inclusive(a0::DateTime, a1::DateTime, b0::DateTime, b1::DateTime)::Bool
    a_lo, a_hi = minmax(a0, a1)
    b_lo, b_hi = minmax(b0, b1)
    return !(a_hi < b_lo || b_hi < a_lo)
end

function _pace_in_tropomi_window(
    t_pace::DateTime,
    t_start::DateTime,
    t_end::DateTime;
    inclusive::Bool,
)::Bool
    if inclusive
        return t_start <= t_pace <= t_end
    else
        return t_start < t_pace < t_end
    end
end

# --- Path collection -----------------------------------------------------------

"""Glob under `prefix`: matches files in the root of `prefix` and in subdirectories (Glob.jl `**` does not match depth-0 alone)."""
function collect_local_paths(dirs::Vector{String}, pattern::AbstractString)::Vector{String}
    out = String[]
    pat = startswith(pattern, "**/") ? pattern[4:end] : pattern
    for d in dirs
        isdir(d) || continue
        prefix = abspath(d)
        for p in glob(pat, prefix)
            isfile(p) && push!(out, abspath(p))
        end
        for p in glob(joinpath("**", pat), prefix)
            isfile(p) && push!(out, abspath(p))
        end
    end
    unique!(sort!(out))
    return out
end

function read_manifest_paths(manifest_path::AbstractString)::Vector{String}
    lines = String[]
    open(manifest_path) do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            startswith(s, '#') && continue
            push!(lines, abspath(expanduser(s)))
        end
    end
    return lines
end

function _filter_tropomi_band(paths::Vector{String}, band::Symbol)::Vector{String}
    tag = band === :BD5 ? "BD5" : band === :BD6 ? "BD6" : error("tropomi_band must be :BD5 or :BD6")
    filter(p -> occursin("L1B_RA_$(tag)_", basename(p)), paths)
end

# --- Records & matching --------------------------------------------------------

struct TropomiRecord
    path::String
    t_start::DateTime
    t_end::DateTime
    orbit::String
    bucket_date::Date   # UTC date of t_start
end

struct PaceRecord
    path::String
    t_pace::DateTime
    bucket_date::Date
    parse_mode::Symbol
end

struct CoincidenceRow
    pace_path::String
    tropomi_path::String
    pace_t::DateTime
    tropomi_t_start::DateTime
    tropomi_t_end::DateTime
    containment_margin_s::Float64
    tropomi_orbit::String
    pace_parse_mode::Symbol
end

function _build_tropomi_records(paths::Vector{String}, cfg::MatchConfig)::Vector{TropomiRecord}
    recs = TropomiRecord[]
    for p in paths
        info = parse_tropomi_granule(p)
        info === nothing && continue
        wanted = cfg.tropomi_band === :BD5 ? "5" : "6"
        info.band != wanted && continue
        push!(recs, TropomiRecord(p, info.t_start, info.t_end, info.orbit, Date(info.t_start)))
    end
    return recs
end

function _build_pace_records(paths::Vector{String}, cfg::MatchConfig)::Vector{PaceRecord}
    recs = PaceRecord[]
    for p in paths
        t, mode = parse_pace_granule(p, cfg)
        (t === nothing || mode === :failed || mode === :skipped_multi || mode === :ambiguous) && continue
        push!(recs, PaceRecord(p, t, Date(t), mode))
    end
    return recs
end

function _tropomi_index_by_date(recs::Vector{TropomiRecord})::Dict{Date,Vector{TropomiRecord}}
    d = Dict{Date,Vector{TropomiRecord}}()
    for r in recs
        push!(get!(Vector{TropomiRecord}, d, r.bucket_date), r)
    end
    return d
end

function _validate_utc_window!(cfg::MatchConfig)
    ws = cfg.utc_window_start
    we = cfg.utc_window_end
    if (ws === nothing) != (we === nothing)
        error("match.utc_window_start and utc_window_end must both be set or both omitted")
    end
    if ws !== nothing && we !== nothing && ws > we
        error("match.utc_window_start must be ≤ utc_window_end")
    end
    cfg.utc_day_stride >= 1 || error("match.utc_day_stride must be ≥ 1")
end

function _filter_pace_by_window(recs::Vector{PaceRecord}, cfg::MatchConfig)::Vector{PaceRecord}
    ws = cfg.utc_window_start
    we = cfg.utc_window_end
    (ws === nothing || we === nothing) && return recs
    anchor = Date(ws)
    stride = max(1, cfg.utc_day_stride)
    filter(recs) do pr
        pr.t_pace < ws && return false
        pr.t_pace > we && return false
        stride == 1 && return true
        return mod((Date(pr.t_pace) - anchor).value, stride) == 0
    end
end

function _filter_tropomi_by_window(recs::Vector{TropomiRecord}, cfg::MatchConfig)::Vector{TropomiRecord}
    ws = cfg.utc_window_start
    we = cfg.utc_window_end
    (ws === nothing || we === nothing) && return recs
    anchor = Date(ws)
    stride = max(1, cfg.utc_day_stride)
    filter(recs) do tr
        !_intervals_overlap_inclusive(tr.t_start, tr.t_end, ws, we) && return false
        stride == 1 && return true
        return mod((tr.bucket_date - anchor).value, stride) == 0
    end
end

function find_coincidences(cfg::MatchConfig)::Vector{CoincidenceRow}
    pace_paths = String[]
    trop_paths = String[]

    if cfg.pace_manifest !== nothing
        append!(pace_paths, read_manifest_paths(cfg.pace_manifest))
    end
    append!(pace_paths, collect_local_paths(cfg.pace_dirs, cfg.pace_glob))

    if cfg.tropomi_manifest !== nothing
        append!(trop_paths, read_manifest_paths(cfg.tropomi_manifest))
    end
    for pat in ("*L1B_RA_BD6*", "*L1B_RA_BD5*")
        append!(trop_paths, collect_local_paths(cfg.tropomi_dirs, pat))
    end
    unique!(sort!(pace_paths))
    unique!(sort!(trop_paths))
    trop_paths = _filter_tropomi_band(trop_paths, cfg.tropomi_band)

    _validate_utc_window!(cfg)
    trop_recs = _filter_tropomi_by_window(_build_tropomi_records(trop_paths, cfg), cfg)
    pace_recs = _filter_pace_by_window(_build_pace_records(pace_paths, cfg), cfg)
    idx = _tropomi_index_by_date(trop_recs)

    rows = CoincidenceRow[]
    k = max(0, cfg.neighbor_days)
    for pr in pace_recs
        dates = [pr.bucket_date + Day(δ) for δ in (-k):k]
        for d in dates
            for tr in get(idx, d, TropomiRecord[])
                if _pace_in_tropomi_window(pr.t_pace, tr.t_start, tr.t_end; inclusive = cfg.inclusive_endpoints)
                    margin = containment_margin_seconds(
                        pr.t_pace, tr.t_start, tr.t_end; inclusive = cfg.inclusive_endpoints,
                    )
                    push!(rows, CoincidenceRow(
                        pr.path, tr.path, pr.t_pace, tr.t_start, tr.t_end,
                        margin, tr.orbit, pr.parse_mode,
                    ))
                end
            end
        end
    end
    return rows
end

function write_coincidences_csv(path::AbstractString, rows::Vector{CoincidenceRow})
    open(path, "w") do io
        println(io, join([
            "pace_path", "tropomi_path", "pace_t", "tropomi_t_start", "tropomi_t_end",
            "containment_margin_s", "tropomi_orbit", "pace_parse_mode",
        ], ','))
        for r in rows
            println(io, join([
                csv_escape(r.pace_path),
                csv_escape(r.tropomi_path),
                csv_escape(string(r.pace_t)),
                csv_escape(string(r.tropomi_t_start)),
                csv_escape(string(r.tropomi_t_end)),
                string(r.containment_margin_s),
                csv_escape(r.tropomi_orbit),
                csv_escape(string(r.pace_parse_mode)),
            ], ','))
        end
    end
end

function write_coincidences_jsonl(path::AbstractString, rows::Vector{CoincidenceRow})
    open(path, "w") do io
        for r in rows
            println(io, JSON.json(Dict(
                "pace_path" => r.pace_path,
                "tropomi_path" => r.tropomi_path,
                "pace_t" => string(r.pace_t),
                "tropomi_t_start" => string(r.tropomi_t_start),
                "tropomi_t_end" => string(r.tropomi_t_end),
                "containment_margin_s" => r.containment_margin_s,
                "tropomi_orbit" => r.tropomi_orbit,
                "pace_parse_mode" => string(r.pace_parse_mode),
            )))
        end
    end
end

function csv_escape(s::AbstractString)::String
    if occursin(r"[\",\n]", s)
        return '"' * replace(s, "\"" => "\"\"") * '"'
    end
    return s
end

# --- Self-tests ------------------------------------------------------------------

function run_self_tests()::Bool
    ok = true

    tpath = "S5P_OFFL_L1B_RA_BD6_20250130T160818_20250130T174948_37828_03_020100_20250130T193517.nc"
    info = parse_tropomi_granule(tpath)
    info === nothing && (println("FAIL tropomi parse"); return false)
    info.t_start == DateTime(2025, 1, 30, 16, 8, 18) || (ok = false)
    info.t_end == DateTime(2025, 1, 30, 17, 49, 48) || (ok = false)

    cfg = MatchConfig(; strict_pace_parse = false, multi_pace_timestamp = :skip)
    tc, pm = parse_pace_granule("PACE_OCI.20250130T170000.L1B.V3.nc", cfg)
    (tc == DateTime(2025, 1, 30, 17, 0, 0) && pm === :canonical) || (ok = false)

    tc_suf, pm_suf = parse_pace_granule("my_prefix_PACE_OCI.20250130T170000.L1B.V3_suffix.nc", cfg)
    (tc_suf == DateTime(2025, 1, 30, 17, 0, 0) && pm_suf === :canonical) || (ok = false)

    tr = DateTime(2025, 1, 30, 16, 0, 0)
    t0 = DateTime(2025, 1, 30, 15, 0, 0)
    t1 = DateTime(2025, 1, 30, 18, 0, 0)
    _pace_in_tropomi_window(tr, t0, t1; inclusive = true) || (ok = false)
    !_pace_in_tropomi_window(tr - Hour(2), t0, t1; inclusive = true) || (ok = false)

    cfg_m = MatchConfig(; strict_pace_parse = false, multi_pace_timestamp = :skip)
    ta, ma = parse_pace_granule("backup_20250130T120000_copy_20250130T170000.nc", cfg_m)
    (ta === nothing && ma === :skipped_multi) || (ok = false)

    # Endpoints: t_pace == T_end inclusive
    te = DateTime(2025, 1, 30, 17, 49, 48)
    _pace_in_tropomi_window(te, info.t_start, info.t_end; inclusive = true) || (ok = false)
    !_pace_in_tropomi_window(te, info.t_start, info.t_end; inclusive = false) || (ok = false)

    if ok
        # End-to-end: two named files, one expected pair (PACE instant inside TROPOMI window)
        d = mktempdir()
        pdir = joinpath(d, "pace")
        tdir = joinpath(d, "tropo")
        mkpath(pdir)
        mkpath(tdir)
        write(joinpath(pdir, "PACE_OCI.20250130T170000.L1B.V3.nc"), "")
        write(joinpath(tdir, "S5P_OFFL_L1B_RA_BD6_20250130T160818_20250130T174948_37828_03_020100_20250130T193517.nc"), "")
        icfg = MatchConfig(; pace_dirs = [pdir], tropomi_dirs = [tdir], tropomi_band = :BD6, neighbor_days = 0)
        irows = find_coincidences(icfg)
        if length(irows) != 1
            ok = false
        else
            r = irows[1]
            r.pace_t == DateTime(2025, 1, 30, 17, 0, 0) || (ok = false)
            r.tropomi_orbit == "37828" || (ok = false)
        end
    end

    if ok
        ws_d = parse_match_window_endpoint("2025-01-07", false)
        we_d = parse_match_window_endpoint("2025-01-07", true)
        (ws_d == DateTime(2025, 1, 7, 0, 0, 0) && we_d == DateTime(2025, 1, 7, 23, 59, 59)) || (ok = false)
    end

    if ok
        # UTC window drops PACE granules outside [start,end]
        d = mktempdir()
        pdir = joinpath(d, "pace")
        tdir = joinpath(d, "tropo")
        mkpath(pdir)
        mkpath(tdir)
        write(joinpath(pdir, "PACE_OCI.20250130T100000.L1B.V3.nc"), "")
        write(joinpath(pdir, "PACE_OCI.20250130T170000.L1B.V3.nc"), "")
        write(joinpath(tdir, "S5P_OFFL_L1B_RA_BD6_20250130T160818_20250130T174948_37828_03_020100_20250130T193517.nc"), "")
        wcfg = MatchConfig(;
            pace_dirs = [pdir],
            tropomi_dirs = [tdir],
            tropomi_band = :BD6,
            utc_window_start = DateTime(2025, 1, 30, 12, 0, 0),
            utc_window_end = DateTime(2025, 1, 30, 18, 0, 0),
        )
        wrows = find_coincidences(wcfg)
        (length(wrows) == 1 && wrows[1].pace_t == DateTime(2025, 1, 30, 17, 0, 0)) || (ok = false)
    end

    if ok
        # utc_day_stride: anchor is Date(window_start); Jan 31 excluded when stride==2
        d = mktempdir()
        pdir = joinpath(d, "pace")
        tdir = joinpath(d, "tropo")
        mkpath(pdir)
        mkpath(tdir)
        write(joinpath(pdir, "PACE_OCI.20250130T170000.L1B.V3.nc"), "")
        write(joinpath(pdir, "PACE_OCI.20250131T170000.L1B.V3.nc"), "")
        write(joinpath(tdir, "S5P_OFFL_L1B_RA_BD6_20250130T160818_20250130T174948_37828_03_020100_20250130T193517.nc"), "")
        write(joinpath(tdir, "S5P_OFFL_L1B_RA_BD6_20250131T160818_20250131T174948_37828_03_020100_20250131T193517.nc"), "")
        scfg = MatchConfig(;
            pace_dirs = [pdir],
            tropomi_dirs = [tdir],
            tropomi_band = :BD6,
            utc_window_start = DateTime(2025, 1, 30, 0, 0, 0),
            utc_window_end = DateTime(2025, 1, 31, 23, 59, 59),
            utc_day_stride = 2,
        )
        srows = find_coincidences(scfg)
        (length(srows) == 1 && srows[1].pace_t == DateTime(2025, 1, 30, 17, 0, 0)) || (ok = false)
    end

    if ok
        println("MatchCoincidence self-tests: OK")
    else
        println("MatchCoincidence self-tests: FAILED")
    end
    return ok
end

end # module
