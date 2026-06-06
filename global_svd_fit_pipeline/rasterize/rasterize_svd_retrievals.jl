"""
rasterize_svd_retrievals.jl

Bin PACE SVD retrieval granules into a regular lat/lon grid using centered
rolling time windows.  One NetCDF output file is written per window containing:
  • sif_radiance_678nm(lon, lat)  — unweighted mean of valid soundings
  • counts(lon, lat)              — number of valid soundings per cell

Usage:
  julia --project=. global_fit_pipeline/rasterize/rasterize_svd_retrievals.jl \\
        global_fit_pipeline/rasterize/rasterize.example.toml

Input files must match:  interim_<YYYYMMDDTHHmmss>_svd_retrieval_*.nc
Output files are named:  sif678_raster_<YYYYMMDD>_<YYYYMMDD>.nc
Filters: [rasterize.filters].status_codes = [1, ...] (or legacy valid_status_only).
         [rasterize.filters].l2_flags_reject = "ALL" | ["FLAG1","FLAG2"] | omit
output_dir is created with mkpath if missing; each window prints granule file count.
"""

using Dates
using NCDatasets
using ProgressMeter
using TOML


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

function _die(msg::String)
    println(stderr, "ERROR: ", msg)
    exit(1)
end

function _get(d::Dict, key::String, default)
    v = get(d, key, nothing)
    return (v === nothing) ? default : v
end

struct RasterConfig
    input_dir::String
    output_dir::String
    l2aop_dir::String        # required when exclude_missing_nflh or l2_flags_reject_mask != 0
    start_date::Date
    end_date::Date
    chunk_frequency_days::Int
    half_chunk_days::Int
    resolution::Int          # 180/resolution = degrees per cell
    status_codes::Union{Nothing, Set{Int16}}  # nothing → no status filter
    exclude_dark::Bool
    exclude_ocean::Bool
    ocean_only::Bool       # keep only is_ocean == 1 (mutually exclusive with exclude_ocean)
    exclude_missing_nflh::Bool
    nflh_var::String
    l2_flags_reject_mask::UInt32  # 0x00000000 = disabled; 0xFFFFFFFF = require flag==0
    max_sif::Float64
    max_chi2::Float64        # Inf → disabled
end

function parse_config(path::String)::RasterConfig
    cfg = TOML.parsefile(path)
    r = get(cfg, "rasterize", nothing)
    r isa Dict || _die("Missing [rasterize] table in $path")

    input_dir  = _get(r, "input_dir",  "")
    output_dir = _get(r, "output_dir", "")
    l2aop_dir  = String(_get(r, "l2aop_dir", ""))
    isempty(input_dir)  && _die("[rasterize] input_dir is required")
    isempty(output_dir) && _die("[rasterize] output_dir is required")

    t0_s = _get(r, "start_date", "")
    t1_s = _get(r, "end_date",   "")
    (isempty(t0_s) || isempty(t1_s)) && _die("[rasterize] start_date and end_date are required")
    t0 = Date(t0_s)
    t1 = Date(t1_s)
    t0 <= t1 || _die("start_date must be ≤ end_date")

    freq = Int(_get(r, "chunk_frequency_days", 16))
    half = Int(_get(r, "half_chunk_days",      8))

    grid = get(r, "grid", Dict{String,Any}())
    res  = Int(_get(grid, "resolution", 180))
    res > 0 || _die("[rasterize.grid] resolution must be > 0")

    filt = get(r, "filters", Dict{String,Any}())
    status_codes = _parse_status_codes(filt)
    exclude_dark  = Bool(_get(filt, "exclude_dark",  true))
    ocean_only    = Bool(_get(filt, "ocean_only",    false))
    exclude_ocean = Bool(_get(filt, "exclude_ocean", !ocean_only))
    ocean_only && exclude_ocean &&
        _die("[rasterize.filters] ocean_only and exclude_ocean cannot both be true")
    exclude_missing_nflh = Bool(_get(filt, "exclude_missing_nflh", false))
    nflh_var = String(_get(filt, "nflh_var", "nflh"))
    exclude_missing_nflh && isempty(l2aop_dir) &&
        _die("[rasterize] l2aop_dir is required when [rasterize.filters].exclude_missing_nflh = true")
    l2_flags_reject_mask = _parse_l2_flags_mask(filt)
    l2_flags_reject_mask != 0 && isempty(l2aop_dir) &&
        _die("[rasterize] l2aop_dir is required when [rasterize.filters].l2_flags_reject is set")
    max_sif  = Float64(_get(filt, "max_sif",  Inf))
    max_chi2 = Float64(_get(filt, "max_chi2", Inf))

    return RasterConfig(
        input_dir, output_dir, l2aop_dir,
        t0, t1, freq, half, res,
        status_codes, exclude_dark, exclude_ocean, ocean_only,
        exclude_missing_nflh, nflh_var,
        l2_flags_reject_mask,
        max_sif, max_chi2,
    )
end

"""Parse `[rasterize.filters].status_codes` or legacy `valid_status_only`."""
function _parse_status_codes(filt::Dict)
    if haskey(filt, "status_codes")
        raw = filt["status_codes"]
        raw isa AbstractVector || _die("[rasterize.filters] status_codes must be a list of integers")
        isempty(raw) && _die("[rasterize.filters] status_codes must not be empty")
        return Set{Int16}(Int16.(raw))
    end
    # Legacy bool: true → converged only (status 1); false → no status filter
    if Bool(_get(filt, "valid_status_only", true))
        return Set{Int16}([Int16(1)])
    end
    return nothing
end

# Bit positions of named L2 AOP flags (flag_meanings order from PACE OCI L2 files).
const _L2_FLAG_BITS = Dict{String, Int}(
    "ATMFAIL"    => 0,   "LAND"       => 1,   "PRODWARN"   => 2,   "HIGLINT"    => 3,
    "HILT"       => 4,   "HISATZEN"   => 5,   "COASTZ"     => 6,
    "STRAYLIGHT" => 8,   "CLDICE"     => 9,   "COCCOLITH"  => 10,  "TURBIDW"    => 11,
    "HISOLZEN"   => 12,  "LOWLW"      => 14,  "CHLFAIL"    => 15,  "NAVWARN"    => 16,
    "ABSAER"     => 17,  "MAXAERITER" => 19,  "MODGLINT"   => 20,  "CHLWARN"    => 21,
    "ATMWARN"    => 22,  "OPSHAL"     => 23,  "SEAICE"     => 24,  "NAVFAIL"    => 25,
    "FILTER"     => 26,  "BOWTIEDEL"  => 28,  "HIPOL"      => 29,  "PRODFAIL"   => 30,
)

"""Parse `[rasterize.filters].l2_flags_reject` into a UInt32 bitmask.

  - `"ALL"` or `true`       → 0xFFFFFFFF  (reject any pixel with any flag set; only flag==0 passes)
  - `["CLDICE", "LAND"]`    → bitmask of those named flags
  - `[9, 1]`                → bitmask from raw bit positions
  - missing / `false` / `""` → 0x00000000 (disabled)
"""
function _parse_l2_flags_mask(filt::Dict)::UInt32
    raw = get(filt, "l2_flags_reject", nothing)
    raw === nothing && return UInt32(0)
    if raw isa Bool
        return raw ? typemax(UInt32) : UInt32(0)
    end
    if raw isa String
        s = uppercase(strip(raw))
        isempty(s) && return UInt32(0)
        s == "ALL" && return typemax(UInt32)
        _die("[rasterize.filters] l2_flags_reject string must be \"ALL\" (got \"$raw\")")
    end
    if raw isa AbstractVector
        mask = UInt32(0)
        for entry in raw
            if entry isa Integer
                (0 <= entry <= 31) || _die("[rasterize.filters] l2_flags_reject bit $entry out of range 0–31")
                mask |= UInt32(1) << entry
            elseif entry isa String
                bit = get(_L2_FLAG_BITS, uppercase(strip(entry)), nothing)
                bit === nothing && _die("[rasterize.filters] unknown l2_flags_reject flag name: \"$entry\"\n" *
                    "  Known names: $(join(sort(collect(keys(_L2_FLAG_BITS))), ", "))")
                mask |= UInt32(1) << bit
            else
                _die("[rasterize.filters] l2_flags_reject entries must be strings or integers, got $(typeof(entry))")
            end
        end
        return mask
    end
    _die("[rasterize.filters] l2_flags_reject must be \"ALL\", a list of flag names/bit integers, or omitted")
end


# ---------------------------------------------------------------------------
# Granule discovery
# ---------------------------------------------------------------------------

# Matches: interim_20250702T123456_svd_retrieval_*.nc
const _GRANULE_RE = r"^interim_(\d{8})T\d{6}_svd_retrieval_.*\.nc$"
const _GRANULE_ID_RE = r"^interim_(\d{8}T\d{6})_svd_retrieval_"
const _L2AOP_FNAME_RE = r"^PACE_OCI\.(\d{8}T\d{6})\.L2\.OC_AOP.*\.nc$"

function _granule_id_from_retrieval_path(fpath::String)
    m = match(_GRANULE_ID_RE, basename(fpath))
    return m === nothing ? nothing : String(m.captures[1])
end

"""Match L2 granule ids on all but the last two time digits (seconds), e.g. `…211157` ↔ `…211148`."""
function _granule_id_coarse_prefix(gid::String)
    length(gid) >= 3 ? gid[1:(end - 2)] : gid
end

function _l2aop_granule_id_from_fname(fname::String)
    m = match(_L2AOP_FNAME_RE, fname)
    return m === nothing ? nothing : String(m.captures[1])
end

function _resolve_l2aop_path(l2aop_dir::String, granule_id::String)
    exact = joinpath(l2aop_dir, "PACE_OCI.$(granule_id).L2.OC_AOP.V3_1.nc")
    isfile(exact) && return exact

    prefix = _granule_id_coarse_prefix(granule_id)
    candidates = Tuple{String, String}[]  # (path, l2_granule_id)
    for fname in readdir(l2aop_dir)
        endswith(fname, ".nc") || continue
        l2_id = _l2aop_granule_id_from_fname(fname)
        l2_id === nothing && continue
        if _granule_id_coarse_prefix(l2_id) == prefix
            push!(candidates, (joinpath(l2aop_dir, fname), l2_id))
        end
    end
    isempty(candidates) && return nothing
    if length(candidates) == 1
        path, l2_id = candidates[1]
        l2_id != granule_id &&
            @info "L2 AOP matched on coarse granule id (ignore last 2 time digits)" retrieval_id=granule_id l2_id=l2_id file=basename(path)
        return path
    end
    # Several granules share the same minute: pick closest seconds (last two digits of time)
    best = argmin(candidates) do (_, l2_id)
        abs(parse(Int, granule_id[(end - 1):end]) - parse(Int, l2_id[(end - 1):end]))
    end
    path, l2_id = candidates[best]
    @info "L2 AOP coarse granule match (multiple candidates)" retrieval_id=granule_id l2_id=l2_id file=basename(path)
    return path
end

function _find_ncvar_optional(ds, varname::String)
    groups_to_check = Pair{String, Any}[]
    try
        if hasproperty(ds, :group)
            for group_name in keys(ds.group)
                push!(groups_to_check, group_name => ds.group[group_name])
            end
        end
    catch
        nothing
    end
    if isempty(groups_to_check)
        push!(groups_to_check, "" => ds)
    end
    for (_, group) in groups_to_check
        haskey(group, varname) || continue
        var = group[varname]
        try
            dimnames(var)
        catch
            continue
        end
        return var
    end
    return haskey(ds, varname) ? ds[varname] : nothing
end

"""Read 2D field as (pixels, scans) from root or child groups."""
function _read_pixels_scans_2d(ds, varname::String)
    v = _find_ncvar_optional(ds, varname)
    v === nothing && return nothing
    ndims(v) == 2 || error("Expected 2D '$varname', got ndims=$(ndims(v))")
    d = collect(String.(dimnames(v)))
    raw = Array(v)
    arr = if d == ["pixels", "scans"] || d == ["pixels_per_line", "number_of_lines"]
        raw
    elseif d == ["scans", "pixels"] || d == ["number_of_lines", "pixels_per_line"]
        permutedims(raw, (2, 1))
    else
        i_pix = findfirst(x -> occursin("pixel", lowercase(x)), d)
        i_scan = findfirst(x -> occursin("line", lowercase(x)) || x == "scans", d)
        if i_pix !== nothing && i_scan !== nothing
            perm = (i_pix, i_scan)
            permutedims(raw, perm)
        else
            error("Unsupported dims for '$varname': $d")
        end
    end
    return arr
end

function _nc_fill_value(v)
    fv = get(v.attrib, "_FillValue", nothing)
    fv === nothing && return nothing
    return fv isa AbstractArray ? first(fv) : fv
end

function _nflh_value_present(v, fillv)
    ismissing(v) && return false
    x = Float64(v)
    !isfinite(x) && return false
    fillv !== nothing && x == Float64(fillv) && return false
    return true
end

"""`BitMatrix` (pixels × scans) aligned with retrieval swath variables."""
function _read_l2aop_masks(
    retrieval_ds,
    l2aop_path::String,
    nflh_var::String,
    need_nflh::Bool,
    l2_flags_reject_mask::UInt32,
)
    need_flags = l2_flags_reject_mask != 0

    nflh_raw, flags_raw = NCDatasets.Dataset(l2aop_path, "r") do l2ds
        # --- nFLH ---
        nflh_mat = nothing
        if need_nflh
            v = _find_ncvar_optional(l2ds, nflh_var)
            if v !== nothing
                fillv = _nc_fill_value(v)
                raw = _read_pixels_scans_2d(l2ds, nflh_var)
                if raw !== nothing
                    present = falses(size(raw)...)
                    @inbounds for j in axes(raw, 2), i in axes(raw, 1)
                        present[i, j] = _nflh_value_present(raw[i, j], fillv)
                    end
                    nflh_mat = present
                end
            end
        end

        # --- l2_flags ---
        flags_mat = nothing
        if need_flags
            v = _find_ncvar_optional(l2ds, "l2_flags")
            if v !== nothing
                raw = _read_pixels_scans_2d(l2ds, "l2_flags")
                if raw !== nothing
                    ok = trues(size(raw)...)
                    @inbounds for j in axes(raw, 2), i in axes(raw, 1)
                        ok[i, j] = (UInt32(raw[i, j]) & l2_flags_reject_mask) == 0
                    end
                    flags_mat = ok
                end
            end
        end

        (nflh_mat, flags_mat)
    end

    # Nothing to align
    nflh_raw === nothing && flags_raw === nothing && return nothing, nothing

    # Build source index arrays (same for both masks — compute once)
    ref_size = nflh_raw !== nothing ? size(nflh_raw) : size(flags_raw)
    n_pix_l2, n_scan_l2 = ref_size
    src_pix = haskey(retrieval_ds, "source_pixel_index") ?
        Vector{Int}(Array(retrieval_ds["source_pixel_index"])) : collect(1:n_pix_l2)
    src_scan = haskey(retrieval_ds, "source_scan_index") ?
        Vector{Int}(Array(retrieval_ds["source_scan_index"])) : collect(1:n_scan_l2)
    n_pix_r  = length(src_pix)
    n_scan_r = length(src_scan)

    function _align(mat)
        mat === nothing && return nothing
        out = falses(n_pix_r, n_scan_r)
        @inbounds for j in 1:n_scan_r, i in 1:n_pix_r
            ip = src_pix[i];  js = src_scan[j]
            if 1 <= ip <= size(mat, 1) && 1 <= js <= size(mat, 2)
                out[i, j] = mat[ip, js]
            end
        end
        return out
    end

    return _align(nflh_raw), _align(flags_raw)
end

"""Return Dict{Date, Vector{String}}: sensing date → list of matching file paths."""
function discover_granules(input_dir::String)::Dict{Date, Vector{String}}
    isdir(input_dir) || _die("input_dir not found: $input_dir")
    date_map = Dict{Date, Vector{String}}()
    for fname in readdir(input_dir)
        m = match(_GRANULE_RE, fname)
        m === nothing && continue
        d = Date(m.captures[1], "yyyymmdd")
        fpath = joinpath(input_dir, fname)
        push!(get!(date_map, d, String[]), fpath)
    end
    return date_map
end


# ---------------------------------------------------------------------------
# Grid accumulation
# ---------------------------------------------------------------------------

struct Grid
    sif::Matrix{Float64}    # n_lon × n_lat, running mean
    counts::Matrix{Int32}   # n_lon × n_lat
    n_lon::Int
    n_lat::Int
    res_step::Float64       # degrees per cell
end

function Grid(resolution::Int)
    n_lon = 2 * resolution
    n_lat = resolution
    res_step = 180.0 / resolution
    Grid(zeros(Float64, n_lon, n_lat), zeros(Int32, n_lon, n_lat), n_lon, n_lat, res_step)
end

@inline function _cell(g::Grid, lon::Float64, lat::Float64)
    i = clamp(floor(Int, (lon + 180.0) / g.res_step) + 1, 1, g.n_lon)
    j = clamp(floor(Int, (lat +  90.0) / g.res_step) + 1, 1, g.n_lat)
    return i, j
end

function accumulate_file!(g::Grid, fpath::String, cfg::RasterConfig)
    NCDatasets.Dataset(fpath, "r") do ds
        lat_raw = _read_pixels_scans_2d(ds, "latitude")
        lat_raw === nothing && error("Missing latitude in $(basename(fpath))")
        lon = Float64.(_read_pixels_scans_2d(ds, "longitude"))
        sif = Float64.(_read_pixels_scans_2d(ds, "sif_radiance_678nm"))
        stat = _read_pixels_scans_2d(ds, "status_code")
        dark = _read_pixels_scans_2d(ds, "is_dark")
        ocean = _read_pixels_scans_2d(ds, "is_ocean")
        chi2 = Float64.(_read_pixels_scans_2d(ds, "reduced_chi2"))
        lat = Float64.(lat_raw)

        flags_ok = nothing
        if cfg.l2_flags_reject_mask != 0
            gid = _granule_id_from_retrieval_path(fpath)
            gid === nothing && error("Cannot parse granule id from $(basename(fpath))")
            l2_path = _resolve_l2aop_path(cfg.l2aop_dir, gid)
            l2_path === nothing && error("L2 AOP not found for granule $gid in $(cfg.l2aop_dir)")
            _, flags_ok = _read_l2aop_masks(
                ds, l2_path, cfg.nflh_var,
                false, cfg.l2_flags_reject_mask,
            )
            cfg.l2_flags_reject_mask != 0 && flags_ok === nothing &&
                @warn "l2_flags not found in L2 AOP file — flag filter skipped" file=basename(l2_path)
        end

        for idx in eachindex(lat)
            (isnan(lat[idx]) || isnan(lon[idx]) || isnan(sif[idx])) && continue
            flags_ok !== nothing && !flags_ok[idx] && continue
            if cfg.status_codes !== nothing && !(Int16(stat[idx]) in cfg.status_codes)
                continue
            end
            cfg.exclude_dark && dark[idx] != 0 && continue
            if cfg.ocean_only
                ocean[idx] == 0 && continue
            elseif cfg.exclude_ocean
                ocean[idx] != 0 && continue
            end
            abs(sif[idx]) > cfg.max_sif  && continue
            !isnan(chi2[idx]) && chi2[idx] > cfg.max_chi2 && continue

            i, j = _cell(g, lon[idx], lat[idx])
            g.counts[i, j] += Int32(1)
            n = g.counts[i, j]
            g.sif[i, j] += (sif[idx] - g.sif[i, j]) / n
        end
    end
end

function accumulate_grid(
    files::Vector{String},
    cfg::RasterConfig;
    progress_desc::String = "granule files",
)
    g = Grid(cfg.resolution)
    n_read = 0
    prog = Progress(length(files); desc = progress_desc)
    for fpath in files
        try
            accumulate_file!(g, fpath, cfg)
            n_read += 1
        catch e
            @warn "Skipping file (read error)" file=basename(fpath) exception=(e, catch_backtrace())
        end
        next!(prog; showvalues = [(:file, basename(fpath))])
    end
    finish!(prog)
    return g, n_read
end


# ---------------------------------------------------------------------------
# NetCDF output
# ---------------------------------------------------------------------------

function save_netcdf(g::Grid, win_start::Date, win_end::Date,
                     output_dir::String, cfg::RasterConfig)
    mkpath(output_dir)
    fname = "sif678_raster_$(Dates.format(win_start, "yyyymmdd"))_$(Dates.format(win_end, "yyyymmdd")).nc"
    fpath = joinpath(output_dir, fname)

    res_step = g.res_step
    lon_axis = collect((-180.0 + res_step/2) : res_step : 180.0)
    lat_axis = collect(( -90.0 + res_step/2) : res_step :  90.0)

    # Build SIF array: cells with no observations → NaN
    sif_out = Float32.(g.sif)
    sif_out[g.counts .== 0] .= NaN32

    NCDatasets.Dataset(fpath, "c") do ds
        ds.attrib["Conventions"]  = "CF-1.8"
        ds.attrib["title"]        = "PACE SVD SIF 678nm rasterized"
        ds.attrib["window_start"] = string(win_start)
        ds.attrib["window_end"]   = string(win_end)
        ds.attrib["input_dir"]    = cfg.input_dir
        ds.attrib["resolution"]   = cfg.resolution
        ds.attrib["res_step_deg"] = res_step
        if cfg.status_codes === nothing
            ds.attrib["status_codes"] = "all"
        else
            ds.attrib["status_codes"] = join(string.(sort(collect(cfg.status_codes))), ",")
        end
        ds.attrib["ocean_only"]    = cfg.ocean_only ? "true" : "false"
        ds.attrib["exclude_ocean"] = cfg.exclude_ocean ? "true" : "false"
        ds.attrib["exclude_missing_nflh"] = cfg.exclude_missing_nflh ? "true" : "false"
        if cfg.exclude_missing_nflh
            ds.attrib["l2aop_dir"] = cfg.l2aop_dir
            ds.attrib["nflh_var"]  = cfg.nflh_var
        end
        if cfg.l2_flags_reject_mask != 0
            ds.attrib["l2_flags_reject_mask"] = string(cfg.l2_flags_reject_mask, base=16, pad=8)
            if cfg.l2_flags_reject_mask == typemax(UInt32)
                ds.attrib["l2_flags_reject"] = "ALL (l2_flags == 0 required)"
            else
                names = [k for (k, v) in _L2_FLAG_BITS if (cfg.l2_flags_reject_mask >> v) & 1 == 1]
                ds.attrib["l2_flags_reject"] = join(sort(names), ",")
            end
            isempty(cfg.l2aop_dir) || (ds.attrib["l2aop_dir"] = cfg.l2aop_dir)
        end
        ds.attrib["history"]      = "Created " * string(now())

        defDim(ds, "lon", g.n_lon)
        defDim(ds, "lat", g.n_lat)

        v_lon = defVar(ds, "lon", Float32, ("lon",),
            attrib = ["units" => "degrees_east",
                      "long_name" => "longitude",
                      "axis" => "X"])
        v_lon[:] = Float32.(lon_axis)

        v_lat = defVar(ds, "lat", Float32, ("lat",),
            attrib = ["units" => "degrees_north",
                      "long_name" => "latitude",
                      "axis" => "Y"])
        v_lat[:] = Float32.(lat_axis)

        v_sif = defVar(ds, "sif_radiance_678nm", Float32, ("lon", "lat"),
            attrib = ["long_name"    => "Mean SIF radiance at 678 nm",
                      "units"        => "W m-2 sr-1 um-1",
                      "_FillValue"   => NaN32,
                      "coordinates"  => "lon lat"])
        v_sif[:, :] = sif_out

        v_cnt = defVar(ds, "counts", Int32, ("lon", "lat"),
            attrib = ["long_name"   => "Number of valid soundings",
                      "units"       => "1",
                      "_FillValue"  => Int32(0),
                      "coordinates" => "lon lat"])
        v_cnt[:, :] = g.counts
    end

    return fpath
end


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    length(ARGS) == 1 || _die("usage: julia rasterize_svd_retrievals.jl <config.toml>")
    cfg_path = ARGS[1]
    isfile(cfg_path) || _die("config not found: $cfg_path")

    cfg = parse_config(cfg_path)
    mkpath(cfg.output_dir)

    status_desc = cfg.status_codes === nothing ? "all status codes" :
        "status_code ∈ {$(join(string.(sort(collect(cfg.status_codes))), ", "))}"
    @info "Rasterize SVD retrievals" input_dir=cfg.input_dir output_dir=cfg.output_dir start=cfg.start_date end_=cfg.end_date resolution=cfg.resolution status_filter=status_desc

    date_map = discover_granules(cfg.input_dir)
    n_granules = sum(length(v) for v in values(date_map))
    @info "Discovered granules" n_dates=length(date_map) n_files=n_granules
    println("Input: $(cfg.input_dir)  ($(n_granules) granule file(s) over $(length(date_map)) date(s))")
    println("Output: $(cfg.output_dir)  (created if missing)")
    println("Status filter: $(status_desc)")
    if cfg.exclude_missing_nflh
        println("nFLH filter: exclude pixels with missing $(cfg.nflh_var) (L2 AOP under $(cfg.l2aop_dir))")
    end
    if cfg.l2_flags_reject_mask != 0
        if cfg.l2_flags_reject_mask == typemax(UInt32)
            println("l2_flags filter: ALL — only l2_flags == 0 pixels pass")
        else
            names = sort([k for (k, v) in _L2_FLAG_BITS if (cfg.l2_flags_reject_mask >> v) & 1 == 1])
            println("l2_flags filter: reject if any of {$(join(names, ", "))} bits set  (mask=0x$(string(cfg.l2_flags_reject_mask, base=16, pad=8)))")
        end
    end

    n_windows  = 0
    n_written  = 0

    for center in cfg.start_date : Day(cfg.chunk_frequency_days) : cfg.end_date
        win_start = center - Day(cfg.half_chunk_days)
        win_end   = center + Day(cfg.half_chunk_days)

        files = String[]
        for (d, fs) in date_map
            win_start <= d <= win_end && append!(files, fs)
        end
        sort!(files)

        n_windows += 1
        win_label = "$(Dates.format(win_start, "yyyy-mm-dd")) – $(Dates.format(win_end, "yyyy-mm-dd"))"
        if isempty(files)
            println("[$n_windows] Window $win_label: 0 granule files — skipping")
            @info "Window has no granules — skipping" window_start=win_start window_end=win_end
            continue
        end

        println("[$n_windows] Window $win_label: $(length(files)) granule file(s)")
        @info "Processing window" window_start=win_start window_end=win_end n_files=length(files)
        prog_label = "window $n_windows ($(length(files)) files)"
        g, n_read = accumulate_grid(files, cfg; progress_desc = prog_label)
        if n_read < length(files)
            println("  read $(n_read)/$(length(files)) file(s) ($(length(files) - n_read) skipped due to errors)")
        end

        n_valid = sum(g.counts .> 0)
        n_obs   = sum(g.counts)
        println("  $(n_obs) valid soundings → $(n_valid) grid cell(s) with data")
        @info "Window complete" valid_cells=n_valid total_obs=n_obs n_files_read=n_read
        if n_obs == 0
            hint = if cfg.ocean_only
                "ocean_only=true but no is_ocean==1 pixels passed other filters; check retrieval watermask / batch_fit.ocean_mask_values"
            else
                "no pixels passed filters (status_codes, exclude_dark, l2_flags_reject, max_sif, max_chi2, …)"
            end
            @warn "Window produced zero observations" window_start=win_start window_end=win_end hint=hint
        end

        fpath = save_netcdf(g, win_start, win_end, cfg.output_dir, cfg)
        println("  wrote $(fpath)")
        @info "Wrote" file=basename(fpath) path=fpath
        n_written += 1
    end

    println("Done: $(n_written) raster file(s) written to $(cfg.output_dir)")
    @info "Done" n_windows=n_windows n_written=n_written
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
