"""
Remote listing/download helpers for GES DISC (TROPOMI L1B RA) and Earthdata PACE via earthaccess (Python).
Included after MatchCoincidence in run_match.jl.
"""
module RemoteFetch

using Base64
using Dates
using Downloads
using JSON
using Printf

import Main.MatchCoincidence as MC

export RemoteTropomiGESDISCConfig,
    RemotePaceEarthAccessConfig,
    gesdisc_default_base_url,
    load_remote_configs,
    fetch_remote_if_requested!,
    run_remote_self_tests

# --- Config structs -------------------------------------------------------------

Base.@kwdef struct RemoteTropomiGESDISCConfig
    enabled::Bool = false
    """Empty uses gesdisc_default_base_url(band)."""
    base_url::String = ""
    date_start::Date
    date_end::Date
    cache_dir::String
    """Maximum number of .nc files to download (nothing = no limit)."""
    max_downloads::Union{Nothing,Int} = nothing
    """Optional curl/wget-style cookie jar path for authenticated sessions."""
    cookie_file::Union{Nothing,String} = nothing
    timeout_s::Int = 120
end

Base.@kwdef struct RemotePaceEarthAccessConfig
    enabled::Bool = false
    short_name::String = "PACE_OCI_L1B_SCI"
    version::String = "3"
    temporal_start::DateTime
    temporal_end::DateTime
    """If false, `earthaccess.search_data` is called without `bounding_box` (global / temporal-only; can return many granules)."""
    use_bounding_box::Bool = true
    bbox_w::Float64 = -180.0
    bbox_s::Float64 = -90.0
    bbox_e::Float64 = 180.0
    bbox_n::Float64 = 90.0
    cache_dir::String
    python_exe::String = "python3"
end

function _parse_cfg_datetime(x)::DateTime
    x isa DateTime && return x
    return DateTime(String(x), dateformat"yyyy-mm-ddTHH:MM:SS")
end

"""Default GES DISC HTTPS root for L1B RA; `hir_collection` is the DAAC folder token after `S5P_L1B_RA_BD{5,6}_` (e.g. `HiR.3`, `HiR2`, `HiR.2` — use the exact name shown on the server)."""
function gesdisc_default_base_url(band::Symbol; hir_collection::AbstractString = "HiR.3")::String
    tag = band === :BD5 ? "BD5" : band === :BD6 ? "BD6" : error("band must be :BD5 or :BD6")
    hc = String(strip(hir_collection))
    isempty(hc) && error("hir_collection must be non-empty")
    occursin('/', hc) && error("hir_collection must not contain '/'")
    return "https://tropomi.gesdisc.eosdis.nasa.gov/data/S5P_TROPOMI_Level1B/S5P_L1B_RA_$(tag)_$(hc)/"
end

"""Doc placeholders like `/path/to/...` make `mkpath` try to create `/path` under the real root and usually fail."""
function _assert_reasonable_cache_dir(label::AbstractString, raw::AbstractString, resolved::AbstractString)
    s = lowercase(strip(raw))
    r = lowercase(resolved)
    bad = occursin("/path/to", s) || r == "/path" || startswith(r, "/path/")
    bad && error(
        "$(label) cache_dir must be a writable folder you choose (documentation placeholders such as " *
        "`/path/to/...` resolve under the real filesystem root `/path` and normally fail permission). " *
        "Use e.g. `~/cache/pace_tropomi_coincidence/tropomi` — tilde expands via expanduser — or `./cache/…` relative to cwd. Got: $(repr(raw))",
    )
    return nothing
end

function _present_nonblank(d::AbstractDict, k::AbstractString)::Bool
    haskey(d, k) || return false
    v = d[k]
    v === nothing && return false
    isa(v, AbstractString) && isempty(strip(v)) && return false
    return true
end

"""When `[remote.*]` omits dates, use `[match]` `utc_window_start`/`utc_window_end` (same interval as coincidence filter)."""
function _shared_utc_window_from_match(cfg_toml::AbstractDict)::Union{Nothing,Tuple{DateTime,DateTime}}
    match_c = get(cfg_toml, "match", Dict{String,Any}())
    ws_raw = get(match_c, "utc_window_start", nothing)
    we_raw = get(match_c, "utc_window_end", nothing)
    ws_miss = ws_raw === nothing || (isa(ws_raw, AbstractString) && isempty(strip(ws_raw)))
    we_miss = we_raw === nothing || (isa(we_raw, AbstractString) && isempty(strip(we_raw)))
    ws_miss && we_miss && return nothing
    ws_miss != we_miss && error("[match] utc_window_start and utc_window_end must both be set")
    ws = MC.parse_match_window_endpoint(ws_raw, false)
    we = MC.parse_match_window_endpoint(we_raw, true)
    return (ws, we)
end

function load_remote_configs(
    cfg_toml::AbstractDict,
    tropomi_band::Symbol,
)::Tuple{Union{Nothing,RemoteTropomiGESDISCConfig},Union{Nothing,RemotePaceEarthAccessConfig}}
    remote = get(cfg_toml, "remote", Dict{String,Any}())
    rt = get(remote, "tropomi", nothing)
    rp = get(remote, "pace", nothing)

    tropomi_cfg = if rt !== nothing && Bool(get(rt, "enabled", false))
        ds, de = if _present_nonblank(rt, "date_start") && _present_nonblank(rt, "date_end")
            (Date(String(rt["date_start"])), Date(String(rt["date_end"])))
        elseif _present_nonblank(rt, "date_start") || _present_nonblank(rt, "date_end")
            error("[remote.tropomi] set both date_start and date_end, or omit both and set [match] utc_window_*")
        else
            mw = _shared_utc_window_from_match(cfg_toml)
            mw === nothing &&
                error("[remote.tropomi] enabled: set date_start/date_end or [match] utc_window_start/utc_window_end")
            (Date(mw[1]), Date(mw[2]))
        end
        base = String(get(rt, "base_url", ""))
        hir = String(get(rt, "gesdisc_hir", get(rt, "hir_collection", "HiR.3")))
        isempty(base) && (base = gesdisc_default_base_url(tropomi_band; hir_collection = hir))
        endswith(base, '/') || (base *= '/')
        mx = get(rt, "max_downloads", nothing)
        mx = mx === nothing || mx === "" ? nothing : Int(mx)
        ck = get(rt, "cookie_file", nothing)
        ck = ck === nothing || isempty(string(ck)) ? nothing : String(ck)
        cdt = get(rt, "cache_dir", nothing)
        (cdt === nothing || (isa(cdt, AbstractString) && isempty(strip(cdt)))) &&
            error("[remote.tropomi] cache_dir is required")
        cache_dir_t = abspath(expanduser(String(cdt)))
        _assert_reasonable_cache_dir("[remote.tropomi]", String(cdt), cache_dir_t)
        RemoteTropomiGESDISCConfig(;
            enabled = true,
            base_url = base,
            date_start = ds,
            date_end = de,
            cache_dir = cache_dir_t,
            max_downloads = mx,
            cookie_file = ck,
            timeout_s = Int(get(rt, "timeout_s", 120)),
        )
    else
        nothing
    end

    pace_cfg = if rp !== nothing && Bool(get(rp, "enabled", false))
        ts, te = if _present_nonblank(rp, "temporal_start") && _present_nonblank(rp, "temporal_end")
            (_parse_cfg_datetime(rp["temporal_start"]), _parse_cfg_datetime(rp["temporal_end"]))
        elseif _present_nonblank(rp, "temporal_start") || _present_nonblank(rp, "temporal_end")
            error("[remote.pace] set both temporal_start and temporal_end, or omit both and set [match] utc_window_*")
        else
            mw = _shared_utc_window_from_match(cfg_toml)
            mw === nothing &&
                error("[remote.pace] enabled: set temporal_start/temporal_end or [match] utc_window_start/utc_window_end")
            (mw[1], mw[2])
        end
        use_bbox = Bool(get(rp, "use_bounding_box", true))
        bb = Float64[-180.0, -90.0, 180.0, 90.0]
        if use_bbox
            bb = Float64[get(rp, "bounding_box", bb)...]
            length(bb) == 4 || error("remote.pace bounding_box must have 4 numbers [west,south,east,north]")
        end
        cdp = get(rp, "cache_dir", nothing)
        (cdp === nothing || (isa(cdp, AbstractString) && isempty(strip(cdp)))) &&
            error("[remote.pace] cache_dir is required")
        cache_dir_p = abspath(expanduser(String(cdp)))
        _assert_reasonable_cache_dir("[remote.pace]", String(cdp), cache_dir_p)
        RemotePaceEarthAccessConfig(;
            enabled = true,
            short_name = String(get(rp, "short_name", "PACE_OCI_L1B_SCI")),
            version = String(get(rp, "version", "3")),
            temporal_start = ts,
            temporal_end = te,
            use_bounding_box = use_bbox,
            bbox_w = bb[1],
            bbox_s = bb[2],
            bbox_e = bb[3],
            bbox_n = bb[4],
            cache_dir = cache_dir_p,
            python_exe = String(get(rp, "python", "python3")),
        )
    else
        nothing
    end

    return tropomi_cfg, pace_cfg
end

# --- HTTP helpers --------------------------------------------------------------

function _basic_auth_headers()::Vector{Pair{String,String}}
    u = get(ENV, "EARTHDATA_USERNAME", "")
    p = get(ENV, "EARTHDATA_PASSWORD", "")
    if !isempty(u) && !isempty(p)
        return ["Authorization" => "Basic $(base64encode(u * ":" * p))"]
    end
    return Pair{String,String}[]
end

function fetch_body(url::AbstractString; headers = _basic_auth_headers(), timeout::Real = 120.0)::String
    io = IOBuffer()
    r = Downloads.request(url; headers, timeout = timeout, output = io)
    (200 <= r.status < 300) || error("HTTP $(r.status) GET $url")
    data = take!(io)
    !isempty(data) && return String(data)
    # Older stdlib `Response` types had `.body`; newer Julia (e.g. 1.12) streams to `output` only.
    if hasfield(typeof(r), :body)
        b = getfield(r, :body)
        b !== nothing && !isempty(b) && return String(b)
    end
    error("Downloads.request returned empty response for GET $url")
end

function download_to_path(url::AbstractString, dest::AbstractString; cookie_file = nothing, timeout::Real = 600.0)
    mkpath(dirname(abspath(dest)))
    headers = _basic_auth_headers()
    try
        Downloads.download(url, dest; headers, timeout = timeout)
    catch e
        cookie_file === nothing && rethrow(e)
        isfile(cookie_file) || rethrow(e)
        run(`curl -f -L --connect-timeout 60 --max-time $timeout -b $(cookie_file) -c $(cookie_file) -o $(dest) $(url)`)
    end
    return dest
end

const RE_HREF = r"href\s*=\s*\"([^\"]+)\""

"""Apache-style directory listing href values."""
function extract_hrefs(html::AbstractString)::Vector{String}
    out = String[]
    for m in eachmatch(RE_HREF, html)
        push!(out, replace(m.captures[1], "&amp;" => "&"))
    end
    return unique!(out)
end

function resolve_href(base_url::AbstractString, href::AbstractString)::String
    if startswith(href, "http://") || startswith(href, "https://")
        return href
    end
    startswith(href, "javascript:") && return href
    m = match(r"^(https?://[^/]+)", base_url)
    m === nothing && error("resolve_href: invalid base_url $base_url")
    origin = m.match
    if startswith(href, "/")
        return origin * href
    end
    base_no_frag = first(split(base_url, '#'))
    prefix = if endswith(base_no_frag, '/')
        base_no_frag
    else
        i = findlast('/', base_no_frag)
        i === nothing ? base_no_frag * "/" : base_no_frag[1:i]
    end
    return prefix * href
end

"""Local mirror path under cache_dir preserving path relative to base_url origin path."""
function mirror_local_path(cache_dir::AbstractString, base_url::AbstractString, file_url::AbstractString)::String
    bu = first(split(base_url, '#'))
    fu = first(split(file_url, '#'))
    m1 = match(r"^https?://[^/]+(.*)$", bu)
    m2 = match(r"^https?://[^/]+(.*)$", fu)
    m1 === nothing && error("mirror_local_path: bad base_url")
    m2 === nothing && error("mirror_local_path: bad file_url")
    p1 = String(m1.captures[1])
    p2 = String(m2.captures[1])
    startswith(p2, p1) || return joinpath(cache_dir, basename(p2))
    rel = p2[length(p1)+1:end]
    rel = lstrip(rel, '/')
    isempty(rel) && return joinpath(cache_dir, "download.nc")
    return joinpath(cache_dir, rel)
end

function _should_follow_dir(href::AbstractString)::Bool
    occursin("Parent Directory", href) && return false
    href == "../" && return false
    startswith(href, "./") && return false
    endswith(href, '/') || return false
    return true
end

function _band_tag(band::Symbol)::String
    band === :BD5 ? "BD5" : band === :BD6 ? "BD6" : error("bad band")
end

function _nc_in_band_and_dates(url::AbstractString, band::Symbol, ds::Date, de::Date)::Bool
    low = lowercase(url)
    endswith(low, ".nc") || return false
    btag = _band_tag(band)
    occursin("L1B_RA_$(btag)_", basename(url)) || return false
    info = MC.parse_tropomi_granule(basename(split(url, '?')[1]))
    info === nothing && return false
    d = Date(info.t_start)
    return ds <= d <= de
end

"""
Breadth-first crawl of GES DISC-style directory pages; collect HTTPS URLs to .nc granules
whose TROPOMI filename parses and whose granule start date lies in [date_start, date_end].
"""
function collect_gesdisc_nc_urls(
    base_url::AbstractString,
    band::Symbol,
    ds::Date,
    de::Date;
    max_urls::Union{Nothing,Int} = nothing,
    timeout_s::Real = 120.0,
)::Vector{String}
    seen_page = Set{String}()
    q = String[base_url]
    out = String[]
    while !isempty(q)
        url = popfirst!(q)
        url = string(first(split(url, '#')))
        url in seen_page && continue
        push!(seen_page, url)
        html = fetch_body(url; timeout = timeout_s)
        for href in extract_hrefs(html)
            (startswith(href, "mailto:") || startswith(href, "javascript:")) && continue
            abs_url = resolve_href(url, href)
            if endswith(lowercase(abs_url), ".nc")
                if _nc_in_band_and_dates(abs_url, band, ds, de)
                    push!(out, abs_url)
                    unique!(out)
                    if max_urls !== nothing && length(out) >= max_urls
                        return out
                    end
                end
            elseif _should_follow_dir(href)
                # prune year directories when possible: .../2025/
                ym = match(r"/(\d{4})/?$", abs_url)
                if ym !== nothing
                    y = parse(Int, ym.captures[1])
                    if y < year(ds) || y > year(de)
                        continue
                    end
                end
                push!(q, abs_url)
            end
        end
    end
    return unique!(out)
end

function fetch_tropomi_gesdisc!(
    cfg::MC.MatchConfig,
    opt::RemoteTropomiGESDISCConfig,
)::Int
    isdir(opt.cache_dir) || mkpath(opt.cache_dir)
    urls = collect_gesdisc_nc_urls(
        opt.base_url,
        cfg.tropomi_band,
        opt.date_start,
        opt.date_end;
        max_urls = opt.max_downloads,
        timeout_s = opt.timeout_s,
    )
    println("  [remote.tropomi] discovered ", length(urls), " .nc URLs under date range")
    n_new = 0
    for url in urls
        dest = mirror_local_path(opt.cache_dir, opt.base_url, url)
        if isfile(dest) && filesize(dest) > 0
            continue
        end
        try
            download_to_path(url, dest; cookie_file = opt.cookie_file, timeout = opt.timeout_s * 5)
            n_new += 1
            @printf("    downloaded: %s\n", basename(dest))
        catch e
            @warn "[remote.tropomi] failed download" url exception = (e, catch_backtrace())
        end
    end
    push!(cfg.tropomi_dirs, opt.cache_dir)
    unique!(cfg.tropomi_dirs)
    return n_new
end

function fetch_pace_earthaccess!(cfg::MC.MatchConfig, opt::RemotePaceEarthAccessConfig)::Int
    isdir(opt.cache_dir) || mkpath(opt.cache_dir)
    script = joinpath(@__DIR__, "fetch_pace_earthaccess.py")
    isfile(script) || error("Missing script: $script")
    payload = Dict{String,Any}(
        "short_name" => opt.short_name,
        "version" => opt.version,
        "temporal_start" => Dates.format(opt.temporal_start, dateformat"yyyy-mm-ddTHH:MM:SS"),
        "temporal_end" => Dates.format(opt.temporal_end, dateformat"yyyy-mm-ddTHH:MM:SS"),
        "use_bounding_box" => opt.use_bounding_box,
        "cache_dir" => opt.cache_dir,
    )
    if opt.use_bounding_box
        payload["bounding_box"] = [opt.bbox_w, opt.bbox_s, opt.bbox_e, opt.bbox_n]
    end
    tmp = joinpath(mktempdir(), "pace_earthaccess.json")
    write(tmp, JSON.json(payload))
    out = try
        read(`$(opt.python_exe) $script $tmp`, String)
    catch e
        error("earthaccess fetch failed (is earthaccess installed for $(opt.python_exe)?): $e")
    end
    data = JSON.parse(out)
    get(data, "ok", false) || error("fetch_pace_earthaccess: $(get(data, "error", out))")
    paths = get(data, "paths", String[])
    println("  [remote.pace] earthaccess reported ", length(paths), " local paths")
    push!(cfg.pace_dirs, opt.cache_dir)
    unique!(cfg.pace_dirs)
    return length(paths)
end

"""
If `do_fetch` is true, run enabled remote downloads and append cache dirs to `cfg`.
"""
function fetch_remote_if_requested!(
    cfg::MC.MatchConfig,
    rt::Union{Nothing,RemoteTropomiGESDISCConfig},
    rp::Union{Nothing,RemotePaceEarthAccessConfig};
    do_fetch::Bool,
)
    do_fetch || return
    println("Remote fetch (--fetch-remote)")
    if rt !== nothing && rt.enabled
        fetch_tropomi_gesdisc!(cfg, rt)
    end
    if rp !== nothing && rp.enabled
        fetch_pace_earthaccess!(cfg, rp)
    end
end

# --- Self-tests ----------------------------------------------------------------

function run_remote_self_tests()::Bool
    ok = true
    html = """<a href="2025/">2025/</a><a href="/data/root/file.nc">file.nc</a>"""
    hs = extract_hrefs(html)
    length(hs) >= 2 || (ok = false)
    base = "https://tropomi.gesdisc.eosdis.nasa.gov/data/S5P_TROPOMI_Level1B/S5P_L1B_RA_BD6_HiR.3/"
    u1 = resolve_href(base, "326/")
    occursin("326/", u1) || (ok = false)
    u2 = resolve_href(base, "/data/other/granule.nc")
    startswith(u2, "https://tropomi.gesdisc.eosdis.nasa.gov/data/other/granule.nc") || (ok = false)

    mpath = mirror_local_path("/cache", base, base * "2025/326/S5P_x_L1B_RA_BD6_20251126T120000_20251126T130000_1.nc")
    occursin("326", mpath) && occursin(".nc", mpath) || (ok = false)

    if ok
        println("RemoteFetch self-tests: OK")
    else
        println("RemoteFetch self-tests: FAILED")
    end
    return ok
end

end # module RemoteFetch
