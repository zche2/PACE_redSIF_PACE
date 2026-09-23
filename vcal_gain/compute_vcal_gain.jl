#!/usr/bin/env julia
# Vicarious gain for a given date.
#
#   julia --project=. vcal_gain/compute_vcal_gain.jl vcal_gain/vcal_gain.example.toml
#   DATE=20260125 julia --project=. vcal_gain/compute_vcal_gain.jl ...
#
# Formula (matches demo_example/AddingSIF/sanity_check/check_relative_residual.jl):
#   vcal_gain = 1 + residual / Rtoa
#
# residual sources:
#   - AddingSIF: spectral 3D `rmse` on the retrieval NetCDF
#   - SVD (2D scalar `rmse`): rebuild ŷ = FM(x_hat) and residual = ŷ − Rtoa from L1B
#
# Pixels with missing/non-finite L2 AOP `nflh` are skipped (gain left NaN).

using Dates
using Glob
using Interpolations
using JLD2
using LinearAlgebra
using NCDatasets
using Statistics
using TOML

const _VCAL_DIR = @__DIR__
const _REPO_ROOT = dirname(_VCAL_DIR)
const _SVD_HELPERS = joinpath(_REPO_ROOT, "global_svd_fit_pipeline", "svd_retrieval", "svd_helpers.jl")
include(_SVD_HELPERS)

const _DIM_RENAME = Dict(
    "pixels_per_line" => "pixels",
    "number_of_lines" => "scans",
    "number_of_bands" => "wavelength",
    "red_bands" => "wavelength",
    "red_wavelength" => "wavelength",
)

const _DEFAULT_SVD_CONFIG = joinpath(
    _REPO_ROOT,
    "svd_configs",
    "global_fit_new_Sep26",
    "multitask.new.Jan25global_fit_pipeline.npoly3nPC15.toml",
)

# ---------------------------------------------------------------------------
# NC helpers
# ---------------------------------------------------------------------------

function _find_var(ds, varname::String)
    groups = Pair{String, Any}[]
    try
        if hasproperty(ds, :group)
            for gname in keys(ds.group)
                push!(groups, gname => ds.group[gname])
            end
        end
    catch
    end
    isempty(groups) && push!(groups, "" => ds)
    for (_, g) in groups
        haskey(g, varname) || continue
        v = g[varname]
        try
            dimnames(v)
        catch
            continue
        end
        return (var = v, dims = collect(String.(dimnames(v))))
    end
    if haskey(ds, varname)
        v = ds[varname]
        return (var = v, dims = collect(String.(dimnames(v))))
    end
    return nothing
end

function _to_pixels_scans_2d(arr, dnames)
    out = [get(_DIM_RENAME, d, d) for d in dnames]
    i_pix = findfirst(==("pixels"), out)
    i_scan = findfirst(==("scans"), out)
    (i_pix !== nothing && i_scan !== nothing) ||
        error("Need pixels/scans dims; got $dnames")
    return permutedims(arr, (i_pix, i_scan))
end

function _to_pixels_scans_bands(arr, dnames)
    out = [get(_DIM_RENAME, d, d) for d in dnames]
    i_pix = findfirst(==("pixels"), out)
    i_scan = findfirst(==("scans"), out)
    i_band = findfirst(==("wavelength"), out)
    isnothing(i_band) && (i_band = findfirst(x -> occursin("band", lowercase(x)), out))
    (i_pix !== nothing && i_scan !== nothing && i_band !== nothing) ||
        error("Need pixels/scans/wavelength dims; got $dnames")
    return permutedims(arr, (i_pix, i_scan, i_band))
end

function _as_float64_nan(A)
    out = similar(A, Float64)
    @inbounds for i in eachindex(A)
        v = A[i]
        out[i] = (v === missing || (v isa Number && !isfinite(v))) ? NaN : Float64(v)
    end
    return out
end

# ---------------------------------------------------------------------------
# Path resolution (flat or YYYY/MM/DD; V3_1 / V3_2)
# ---------------------------------------------------------------------------

function _day_subdir(root::AbstractString, granule_id::AbstractString)
    length(granule_id) >= 8 || return nothing
    d = joinpath(root, granule_id[1:4], granule_id[5:6], granule_id[7:8])
    return isdir(d) ? d : nothing
end

function _search_dirs(root::AbstractString, granule_id::AbstractString)
    dirs = String[root]
    day = _day_subdir(root, granule_id)
    day !== nothing && pushfirst!(dirs, day)
    return dirs
end

function _resolve_product(root::AbstractString, granule_id::AbstractString, kind::Symbol)
    # kind ∈ (:l1b, :l2aop, :l2bgc)
    prefixes = if kind === :l1b
        ["PACE_OCI.$(granule_id).L1B."]
    elseif kind === :l2aop
        ["PACE_OCI.$(granule_id).L2.OC_AOP."]
    elseif kind === :l2bgc
        ["PACE_OCI.$(granule_id).L2.OC_BGC."]
    else
        error("unknown product kind $kind")
    end
    for dir in _search_dirs(root, granule_id)
        isdir(dir) || continue
        for fname in readdir(dir)
            endswith(fname, ".nc") || continue
            any(startswith(fname, p) for p in prefixes) && return joinpath(dir, fname)
        end
    end
    # Coarse match: ignore last 2 seconds digits
    length(granule_id) >= 3 || return nothing
    prefix = granule_id[1:(end - 2)]
    candidates = String[]
    for dir in _search_dirs(root, granule_id)
        isdir(dir) || continue
        for fname in readdir(dir)
            endswith(fname, ".nc") || continue
            any(startswith(fname, p) for p in
                (kind === :l1b ? ["PACE_OCI.$(prefix)"] :
                 kind === :l2aop ? ["PACE_OCI.$(prefix)"] :
                 ["PACE_OCI.$(prefix)"])) || continue
            occursin(kind === :l1b ? ".L1B." :
                     kind === :l2aop ? ".L2.OC_AOP." : ".L2.OC_BGC.", fname) || continue
            push!(candidates, joinpath(dir, fname))
        end
    end
    isempty(candidates) && return nothing
    return sort(candidates)[1]
end

function _granule_id_from_retrieval(path::AbstractString)
    b = basename(path)
    m = match(r"PACE_OCI\.(\d{8}T\d{6})", b)
    m !== nothing && return String(m.captures[1])
    m = match(r"interim_(\d{8}T\d{6})", b)
    m !== nothing && return String(m.captures[1])
    return nothing
end

function _discover_retrievals(retrieval_dir::AbstractString, date::AbstractString, pattern::AbstractString)
    pat = replace(pattern, "{date}" => date)
    files = String[]
    # Prefer date subfolder retrieval_dir/YYYY/MM/DD (SVD / PACE layout).
    day_dir = nothing
    if length(date) == 8
        day_dir = joinpath(retrieval_dir, date[1:4], date[5:6], date[7:8])
        if isdir(day_dir)
            append!(files, Glob.glob(pat, day_dir))
        end
    end
    # Also allow flat retrieval_dir (AddingSIF-style).
    append!(files, Glob.glob(pat, retrieval_dir))
    sort!(unique!(files))
    return files, day_dir
end

# ---------------------------------------------------------------------------
# Read L2 fields aligned to retrieval subset
# ---------------------------------------------------------------------------

function _read_l2_field(path::AbstractString, varname::AbstractString)
    NCDataset(path) do ds
        info = _find_var(ds, varname)
        info === nothing && error("$varname not found in $path")
        raw = _as_float64_nan(Array(info.var))
        return _to_pixels_scans_2d(raw, info.dims)
    end
end

"""Subset full-swath L2 map to the retrieval window using pixel_start/end attrs (1-based inclusive)."""
function _subset_to_retrieval(full::AbstractMatrix, ret_ds; n_pix::Int, n_scan::Int)
    ps = Int(get(ret_ds.attrib, "pixel_start", 1))
    pe = Int(get(ret_ds.attrib, "pixel_end", size(full, 1)))
    ss = Int(get(ret_ds.attrib, "scan_start", 1))
    se = Int(get(ret_ds.attrib, "scan_end", size(full, 2)))
    pe = min(pe, size(full, 1))
    se = min(se, size(full, 2))
    sub = full[ps:pe, ss:se]
    size(sub, 1) == n_pix && size(sub, 2) == n_scan && return sub
    # Fallback: if sizes already match full product, assume retrieval is full swath
    if size(full, 1) == n_pix && size(full, 2) == n_scan
        return full
    end
    error("Cannot align L2 $(size(full)) to retrieval ($n_pix, $n_scan) with " *
          "pixel=[$ps,$pe] scan=[$ss,$se] → $(size(sub))")
end

# ---------------------------------------------------------------------------
# Rtoa from L1B
# ---------------------------------------------------------------------------

function _rtoa_from_l1b(l1b_path::AbstractString, λ_ret::AbstractVector{<:Real};
                        pixel_start::Int, pixel_end::Int, scan_start::Int, scan_end::Int)
    NCDataset(l1b_path) do ds
        rhot_info = _find_var(ds, "rhot_red")
        sol_info = _find_var(ds, "red_solar_irradiance")
        sza_info = _find_var(ds, "solar_zenith")
        wl_info = _find_var(ds, "red_wavelength")
        (rhot_info === nothing || sol_info === nothing ||
         sza_info === nothing || wl_info === nothing) &&
            error("L1B missing rhot_red / solar / sza / wavelength: $l1b_path")

        rhot = _as_float64_nan(Array(rhot_info.var))
        rhot = _to_pixels_scans_bands(rhot, rhot_info.dims)
        F0 = vec(_as_float64_nan(Array(sol_info.var)))
        sza = _as_float64_nan(Array(sza_info.var))
        sza = _to_pixels_scans_2d(sza, sza_info.dims)
        wl = vec(_as_float64_nan(Array(wl_info.var)))
        es = Float64(get(ds.attrib, "earth_sun_distance_correction", 1.0))

        pe = min(pixel_end, size(rhot, 1))
        se = min(scan_end, size(rhot, 2))
        rhot = rhot[pixel_start:pe, scan_start:se, :]
        sza = sza[pixel_start:pe, scan_start:se]

        # Lt = rhot * F0 * cos(SZA) / (π * es)
        Rtoa = similar(rhot)
        for b in axes(rhot, 3)
            @. Rtoa[:, :, b] = rhot[:, :, b] * F0[b] * cosd(sza) / (π * es)
        end

        # Interpolate onto retrieval wavelengths
        itp = interpolate(
            (1:size(Rtoa, 1), 1:size(Rtoa, 2), wl),
            Rtoa,
            (NoInterp(), NoInterp(), Gridded(Linear())),
        )
        out = Array{Float64}(undef, size(Rtoa, 1), size(Rtoa, 2), length(λ_ret))
        for j in axes(out, 2), i in axes(out, 1), k in eachindex(λ_ret)
            out[i, j, k] = itp[i, j, Float64(λ_ret[k])]
        end
        return out
    end
end

"""Band-native Rtoa + solar in (λ_min, λ_max), matching SVD interim wavelength selection."""
function _l1b_rtoa_window(l1b_path::AbstractString;
                          λ_min::Float64, λ_max::Float64,
                          pixel_start::Int, pixel_end::Int,
                          scan_start::Int, scan_end::Int)
    NCDataset(l1b_path) do ds
        rhot_info = _find_var(ds, "rhot_red")
        sol_info = _find_var(ds, "red_solar_irradiance")
        sza_info = _find_var(ds, "solar_zenith")
        wl_info = _find_var(ds, "red_wavelength")
        (rhot_info === nothing || sol_info === nothing ||
         sza_info === nothing || wl_info === nothing) &&
            error("L1B missing rhot_red / solar / sza / wavelength: $l1b_path")

        wl = vec(_as_float64_nan(Array(wl_info.var)))
        F0 = vec(_as_float64_nan(Array(sol_info.var)))
        length(wl) == length(F0) || error("L1B wavelength/F0 length mismatch")
        es = Float64(get(ds.attrib, "earth_sun_distance_correction", 1.0))

        perm = sortperm(wl)
        wl_sorted = wl[perm]
        slice = findall(x -> λ_min < x < λ_max, wl_sorted)
        isempty(slice) && error("No L1B red bands with λ_min < λ < λ_max ($λ_min, $λ_max)")
        band_idx = perm[slice]   # original band indices, ascending λ
        λ = wl[band_idx]
        F0_win = F0[band_idx]

        rhot = _as_float64_nan(Array(rhot_info.var))
        rhot = _to_pixels_scans_bands(rhot, rhot_info.dims)
        sza = _to_pixels_scans_2d(_as_float64_nan(Array(sza_info.var)), sza_info.dims)

        pe = min(pixel_end, size(rhot, 1))
        se = min(scan_end, size(rhot, 2))
        rhot = rhot[pixel_start:pe, scan_start:se, band_idx]
        sza = sza[pixel_start:pe, scan_start:se]

        Rtoa = similar(rhot)
        for b in axes(rhot, 3)
            @. Rtoa[:, :, b] = rhot[:, :, b] * F0_win[b] * cosd(sza) / (π * es)
        end
        return λ, Rtoa, F0_win, es, sza
    end
end

# ---------------------------------------------------------------------------
# SVD residual rebuild
# ---------------------------------------------------------------------------

function _resolve_svd_config(ret_ds, svd_config_override::AbstractString)
    candidates = String[]
    isempty(svd_config_override) || push!(candidates, svd_config_override)
    attr = String(get(ret_ds.attrib, "config_file", ""))
    isempty(attr) || push!(candidates, attr)
    # Legacy attrib path → Sep26 replacement
    if !isempty(attr) && occursin("/global_fit_new/", attr) && !occursin("global_fit_new_Sep26", attr)
        push!(candidates, replace(attr, "/global_fit_new/" => "/global_fit_new_Sep26/"))
    end
    push!(candidates, _DEFAULT_SVD_CONFIG)
    for c in candidates
        isfile(c) && return c
    end
    error("No usable SVD config. Tried:\n  " * join(candidates, "\n  "))
end

"""Minimal SIF-basis loader (same math as SimplePACEXSecFitMWEFunctions.load_sif_basis)."""
function _load_sif_basis(sif_path::AbstractString, λ::AbstractVector{<:Real};
                         nEV::Int = 1, normalize::Bool = true)
    sif = JLD2.load(sif_path)
    sif_u = convert.(Float64, sif["SIF_U"])
    λ_ref = collect(Float64.(sif["SIF_wavelen"]))
    size(sif_u, 1) == length(λ_ref) || error("SIF_U first dim must match SIF_wavelen")
    n_use = min(nEV, size(sif_u, 2))
    dλ = diff(λ_ref)
    step = dλ[1]
    tol = max(1e-10, abs(step) * 1e-8)
    use_cubic = all(abs.(dλ .- step) .<= tol)
    λ_knots = use_cubic ? range(λ_ref[1], step = step, length = length(λ_ref)) : nothing
    basis = zeros(Float64, length(λ), n_use)
    for iev in 1:n_use
        if use_cubic
            itp_ev = CubicSplineInterpolation(λ_knots, sif_u[:, iev]; extrapolation_bc = Line())
            vals = itp_ev.(λ)
            vals[(λ .< λ_ref[1]) .| (λ .> λ_ref[end])] .= 0.0
            basis[:, iev] .= vals
        else
            itp_ev = LinearInterpolation(λ_ref, sif_u[:, iev], extrapolation_bc = 0.0)
            basis[:, iev] .= itp_ev.(λ)
        end
    end
    if normalize
        s = maximum(abs.(basis[:, 1]))
        s > 0 && (basis ./= s)
    end
    return basis
end

function _abs_data_path(base_dir::AbstractString, rel::AbstractString)
    isabspath(rel) ? String(rel) : joinpath(base_dir, rel)
end

"""
solar_eff scale relative to F0*μ0/es.

- `1.0`: current code (π only inside FM) — `solar_eff = F0*μ0/es`
- `1/π`: pre-fix double-π retrievals — `solar_eff = F0*μ0/(π*es)` while FM still `/π`
"""
function _detect_svd_solar_scale(
    x_hat::Array{Float64, 3},
    Rtoa::Array{Float64, 3},
    sza::Matrix{Float64},
    F0::Vector{Float64},
    es::Float64,
    PCs_n::Matrix{Float64},
    sif_basis::Matrix{Float64},
    leg_basis::Matrix{Float64},
    layout;
    log_transform::Bool,
    converged::Union{Nothing, AbstractMatrix},
    n_sample::Int = 64,
)
    n_pix, n_scan, n_λ = size(Rtoa)
    candidates = Float64[1.0, 1 / π]
    # Collect sample pixel indices (prefer converged)
    idxs = Tuple{Int, Int}[]
    for j in 1:n_scan, i in 1:n_pix
        if converged !== nothing
            cv = converged[i, j]
            (cv == true || cv == 1 || cv == 1.0) || continue
        end
        xok = true
        @inbounds for k in 1:layout.n_state
            if !isfinite(x_hat[i, j, k])
                xok = false
                break
            end
        end
        xok || continue
        yok = true
        @inbounds for k in 1:n_λ
            if !isfinite(Rtoa[i, j, k]) || abs(Rtoa[i, j, k]) < 1e-12
                yok = false
                break
            end
        end
        yok || continue
        isfinite(sza[i, j]) || continue
        push!(idxs, (i, j))
        length(idxs) >= n_sample && break
    end
    isempty(idxs) && return 1.0, "once_in_fm (no sample; default)"

    x_buf = Vector{Float64}(undef, layout.n_state)
    solar_eff = Vector{Float64}(undef, n_λ)
    yhat = Vector{Float64}(undef, n_λ)
    best_scale = 1.0
    best_med = Inf
    best_label = ""
    for scale in candidates
        abs_rel = Float64[]
        for (i, j) in idxs
            @inbounds for k in 1:layout.n_state
                x_buf[k] = x_hat[i, j, k]
            end
            μ0 = cosd(sza[i, j])
            @. solar_eff = F0 * μ0 / es * scale
            c_vec = @view x_buf[layout.idx_pc]
            alpha_raw = x_buf[first(layout.idx_alpha)]
            alpha_coeff = 10.0 / (1.0 + exp(-alpha_raw)) + 1.0
            leg_coeff = @view x_buf[layout.idx_legendre]
            sif_coeff = @view x_buf[layout.idx_sif]
            trans_up = log_transform ? exp.(PCs_n * c_vec) : 1.0 .+ PCs_n * c_vec
            trans_updown = exp.(alpha_coeff .* log.(max.(trans_up, eps(Float64))))
            rho_obs = leg_basis * leg_coeff
            sif_toa = trans_up .* (sif_basis * sif_coeff)
            @. yhat = solar_eff * trans_updown * rho_obs / π + sif_toa
            # median |yhat/yobs - 1|
            rels = abs.(yhat ./ (@view Rtoa[i, j, :]) .- 1)
            push!(abs_rel, median(rels))
        end
        med = median(abs_rel)
        if med < best_med
            best_med = med
            best_scale = scale
            best_label = isapprox(scale, 1 / π) ?
                "also_in_solar_eff (double-π era; med|ŷ/y-1|=$(round(med; digits=4)))" :
                "once_in_fm (med|ŷ/y-1|=$(round(med; digits=4)))"
        end
    end
    return best_scale, best_label
end

function _rebuild_svd_residual!(
    residual::Array{Float64, 3},
    x_hat::Array{Float64, 3},
    Rtoa::Array{Float64, 3},
    sza::Matrix{Float64},
    F0::Vector{Float64},
    es::Float64,
    λ::Vector{Float64},
    PCs::Matrix{Float64},
    sif_basis::Matrix{Float64};
    n_pc::Int,
    n_legendre::Int,
    log_transform::Bool,
    valid_mask::BitMatrix,
    converged::Union{Nothing, AbstractMatrix} = nothing,
    solar_eff_scale::Float64 = 1.0,
)
    n_pix, n_scan, n_λ = size(Rtoa)
    size(residual) == (n_pix, n_scan, n_λ) || error("residual size mismatch")
    size(x_hat, 1) == n_pix && size(x_hat, 2) == n_scan || error("x_hat spatial size mismatch")
    fill!(residual, NaN)

    # Hoist FM pieces (solar_eff varies per pixel; π applied once here as in make_svd_forward_model_λ)
    PCs_n = Float64.(PCs[:, 1:n_pc])
    SIF = Float64.(sif_basis)
    layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_legendre, n_ev = size(SIF, 2))
    size(x_hat, 3) == layout.n_state ||
        error("x_hat state length $(size(x_hat, 3)) != layout.n_state $(layout.n_state)")
    leg_basis = _legendre_design_matrix(_normalized_grid(λ), n_legendre)
    solar_eff = Vector{Float64}(undef, n_λ)
    yhat = Vector{Float64}(undef, n_λ)
    x_buf = Vector{Float64}(undef, layout.n_state)

    n_done = 0
    for j in 1:n_scan, i in 1:n_pix
        valid_mask[i, j] || continue
        if converged !== nothing
            cv = converged[i, j]
            (cv == true || cv == 1 || cv == 1.0) || continue
        end
        @inbounds for k in 1:layout.n_state
            x_buf[k] = x_hat[i, j, k]
        end
        all(isfinite, x_buf) || continue
        yobs_ok = true
        @inbounds for k in 1:n_λ
            if !isfinite(Rtoa[i, j, k])
                yobs_ok = false
                break
            end
        end
        yobs_ok || continue
        μ0 = cosd(sza[i, j])
        isfinite(μ0) || continue
        @. solar_eff = F0 * μ0 / es * solar_eff_scale

        c_vec = @view x_buf[layout.idx_pc]
        alpha_raw = x_buf[first(layout.idx_alpha)]
        alpha_coeff = 10.0 / (1.0 + exp(-alpha_raw)) + 1.0
        leg_coeff = @view x_buf[layout.idx_legendre]
        sif_coeff = @view x_buf[layout.idx_sif]
        trans_up = log_transform ? exp.(PCs_n * c_vec) : 1.0 .+ PCs_n * c_vec
        trans_updown = exp.(alpha_coeff .* log.(max.(trans_up, eps(Float64))))
        rho_obs = leg_basis * leg_coeff
        sif_toa = trans_up .* (SIF * sif_coeff)
        @. yhat = solar_eff * trans_updown * rho_obs / π + sif_toa

        @inbounds for k in 1:n_λ
            residual[i, j, k] = yhat[k] - Rtoa[i, j, k]
        end
        n_done += 1
    end
    return n_done
end

# ---------------------------------------------------------------------------
# Core: one granule
# ---------------------------------------------------------------------------

function compute_vcal_gain_granule(
    retrieval_path::AbstractString,
    l2aop_path::AbstractString,
    l2bgc_path::Union{AbstractString, Nothing},
    l1b_path::Union{AbstractString, Nothing},
    output_dir::AbstractString;
    svd_config::AbstractString = "",
)
    mkpath(output_dir)
    granule_id = _granule_id_from_retrieval(retrieval_path)
    granule_id === nothing && error("Cannot parse granule id from $(basename(retrieval_path))")

    NCDataset(retrieval_path) do ds
        rmse_info = _find_var(ds, "rmse")
        rmse_info === nothing && error("No `rmse` in $retrieval_path")
        rmse_ndims = ndims(rmse_info.var)

        residual = nothing
        λ = nothing
        Rtoa = nothing
        n_pix = n_scan = n_λ = 0
        residual_source = ""
        svd_config_used = ""

        if rmse_ndims == 3
            residual = _to_pixels_scans_bands(_as_float64_nan(Array(rmse_info.var)), rmse_info.dims)
            n_pix, n_scan, n_λ = size(residual)
            wl_info = _find_var(ds, "red_wavelength")
            wl_info === nothing && (wl_info = _find_var(ds, "wavelength"))
            wl_info === nothing && error("No wavelength variable in $retrieval_path")
            λ = vec(_as_float64_nan(Array(wl_info.var)))
            length(λ) == n_λ || error("wavelength length $(length(λ)) ≠ residual bands $n_λ")

            rtoa_info = _find_var(ds, "Rtoa_red")
            Rtoa = if rtoa_info !== nothing && ndims(rtoa_info.var) == 3
                _to_pixels_scans_bands(_as_float64_nan(Array(rtoa_info.var)), rtoa_info.dims)
            else
                l1b_path === nothing && error("Rtoa_red missing and no L1B path for $granule_id")
                ps = Int(get(ds.attrib, "pixel_start", 1))
                pe = Int(get(ds.attrib, "pixel_end", n_pix))
                ss = Int(get(ds.attrib, "scan_start", 1))
                se = Int(get(ds.attrib, "scan_end", n_scan))
                _rtoa_from_l1b(l1b_path, λ; pixel_start = ps, pixel_end = pe, scan_start = ss, scan_end = se)
            end
            size(Rtoa) == size(residual) ||
                error("Rtoa size $(size(Rtoa)) ≠ residual $(size(residual))")
            residual_source = "retrieval_rmse_3d"

        elseif rmse_ndims == 2
            # SVD scalar rmse → rebuild spectral residual from x_hat + L1B
            x_info = _find_var(ds, "x_hat")
            x_info === nothing &&
                error("SVD retrieval lacks `x_hat`; cannot rebuild spectral residual: $retrieval_path")
            l1b_path === nothing &&
                error("SVD residual rebuild requires L1B for $granule_id")

            x_hat = _as_float64_nan(Array(x_info.var))
            # Expect (pixels, scans, state) or permute if needed
            xd = [get(_DIM_RENAME, d, d) for d in x_info.dims]
            if xd == ["pixels", "scans", "state"] || (length(xd) == 3 && xd[1] == "pixels" && xd[2] == "scans")
                # already (pix, scan, state)
            elseif findfirst(==("pixels"), xd) !== nothing
                ip = findfirst(==("pixels"), xd)
                isc = findfirst(==("scans"), xd)
                ist = findfirst(d -> d == "state" || occursin("state", d), xd)
                ist === nothing && (ist = findfirst(i -> i != ip && i != isc, 1:3))
                x_hat = permutedims(x_hat, (ip, isc, ist))
            else
                error("Unexpected x_hat dims=$(x_info.dims)")
            end
            n_pix, n_scan = size(x_hat, 1), size(x_hat, 2)

            cfg_path = _resolve_svd_config(ds, svd_config)
            svd_config_used = cfg_path
            pipe = TOML.parsefile(cfg_path)
            spectral_cfg = get(pipe, "spectral", Dict{String, Any}())
            fit_cfg = get(pipe, "fit", Dict{String, Any}())
            svd_cfg = get(fit_cfg, "svd", Dict{String, Any}())
            data_cfg = get(pipe, "data", Dict{String, Any}())
            λ_min = Float64(get(spectral_cfg, "lambda_min_nm", 640.0))
            λ_max = Float64(get(spectral_cfg, "lambda_max_nm", 756.0))
            n_pc = Int(get(svd_cfg, "n_pc", 15))
            n_leg = Int(get(svd_cfg, "n_legendre", 3))
            log_trans = Bool(get(svd_cfg, "svd_log_transform", false))
            sif_nev = Int(get(spectral_cfg, "sif_nev", 1))
            normalize_sif = Bool(get(spectral_cfg, "normalize_sif_first_ev", true))
            base_dir = String(get(data_cfg, "base_dir", ""))
            summer_abs = _abs_data_path(base_dir, String(data_cfg["summer_nc"]))
            winter_abs = _abs_data_path(base_dir, String(data_cfg["winter_nc"]))
            sif_path = _abs_data_path(base_dir, String(data_cfg["sif_file"]))

            ps = Int(get(ds.attrib, "pixel_start", 1))
            pe = Int(get(ds.attrib, "pixel_end", n_pix))
            ss = Int(get(ds.attrib, "scan_start", 1))
            se = Int(get(ds.attrib, "scan_end", n_scan))
            λ, Rtoa, F0, es, sza = _l1b_rtoa_window(
                l1b_path; λ_min = λ_min, λ_max = λ_max,
                pixel_start = ps, pixel_end = pe, scan_start = ss, scan_end = se,
            )
            size(Rtoa, 1) == n_pix && size(Rtoa, 2) == n_scan ||
                error("L1B Rtoa spatial $(size(Rtoa)[1:2]) ≠ retrieval ($n_pix, $n_scan)")
            n_λ = length(λ)
            size(x_hat, 3) == n_pc + 1 + (n_leg + 1) + sif_nev ||
                @warn "x_hat state length $(size(x_hat,3)) vs expected n_pc+α+leg+sif" n_pc n_leg sif_nev

            println("    SVD rebuild: config=$(basename(cfg_path))  nλ=$n_λ  " *
                    "λ∈($λ_min,$λ_max)  n_pc=$n_pc")
            svd_basis = load_svd_basis(
                summer_abs, winter_abs, λ;
                λ_min = λ_min, λ_max = λ_max, n_pc = n_pc, log_transform = log_trans,
            )
            PCs = Float64.(svd_basis.PCs[:, 1:n_pc])
            sif_basis = _load_sif_basis(sif_path, λ; nEV = sif_nev, normalize = normalize_sif)

            # L2 nflh mask needed before rebuild loop — read early
            nflh_early = _subset_to_retrieval(
                _read_l2_field(l2aop_path, "nflh"), ds; n_pix = n_pix, n_scan = n_scan,
            )
            valid_nflh = isfinite.(nflh_early)
            conv = nothing
            if (info = _find_var(ds, "converged")) !== nothing
                conv = Array(info.var)
                try
                    conv = _to_pixels_scans_2d(_as_float64_nan(conv), info.dims)
                catch
                    conv = reshape(Float64.(conv), n_pix, n_scan)
                end
            end

            layout0 = svd_state_layout(; n_pc = n_pc, n_legendre = n_leg, n_ev = size(sif_basis, 2))
            PCs_n0 = Float64.(PCs[:, 1:n_pc])
            leg0 = _legendre_design_matrix(_normalized_grid(λ), n_leg)
            solar_scale, solar_label = _detect_svd_solar_scale(
                x_hat, Rtoa, sza, F0, es, PCs_n0, Float64.(sif_basis), leg0, layout0;
                log_transform = log_trans, converged = conv,
            )
            println("    SVD rebuild: solar_eff convention → $solar_label")

            residual = fill(NaN, n_pix, n_scan, n_λ)
            n_fm = _rebuild_svd_residual!(
                residual, x_hat, Rtoa, sza, F0, es, λ, PCs, sif_basis;
                n_pc = n_pc, n_legendre = n_leg, log_transform = log_trans,
                valid_mask = valid_nflh, converged = conv,
                solar_eff_scale = solar_scale,
            )
            println("    SVD rebuild: FM evaluated on $n_fm pixels (valid nflh ∩ converged)")
            residual_source = "svd_fm_xhat_minus_l1b_rtoa"

            # Stash early nflh for write path below
            lat_info = _find_var(ds, "latitude")
            lon_info = _find_var(ds, "longitude")
            lat = lat_info === nothing ? fill(NaN, n_pix, n_scan) :
                  _to_pixels_scans_2d(_as_float64_nan(Array(lat_info.var)), lat_info.dims)
            lon = lon_info === nothing ? fill(NaN, n_pix, n_scan) :
                  _to_pixels_scans_2d(_as_float64_nan(Array(lon_info.var)), lon_info.dims)

            sif = nothing
            if (info = _find_var(ds, "sif_radiance_678nm")) !== nothing
                sif = _to_pixels_scans_2d(_as_float64_nan(Array(info.var)), info.dims)
            end

            chlor = fill(NaN, n_pix, n_scan)
            if l2bgc_path !== nothing && isfile(l2bgc_path)
                chlor = _subset_to_retrieval(
                    _read_l2_field(l2bgc_path, "chlor_a"), ds; n_pix = n_pix, n_scan = n_scan,
                )
            end

            n_valid = count(valid_nflh)
            n_skip = n_pix * n_scan - n_valid
            gain = fill(NaN, n_pix, n_scan, n_λ)
            @inbounds for j in 1:n_scan, i in 1:n_pix
                valid_nflh[i, j] || continue
                for k in 1:n_λ
                    r = Rtoa[i, j, k]
                    (!isfinite(r) || abs(r) < 1e-12) && continue
                    res = residual[i, j, k]
                    isfinite(res) || continue
                    gain[i, j, k] = 1.0 + res / r
                end
            end

            return _write_vcal_gain(
                retrieval_path, l2aop_path, l2bgc_path, l1b_path, output_dir, granule_id,
                λ, lat, lon, gain, nflh_early, chlor, sif, conv,
                n_valid, n_skip, residual_source, svd_config_used;
                svd_solar_eff_note = solar_label,
            )
        else
            error("Unexpected `rmse` ndims=$rmse_ndims in $retrieval_path")
        end

        # --- AddingSIF / spectral-rmse path ---
        lat_info = _find_var(ds, "latitude")
        lon_info = _find_var(ds, "longitude")
        lat = lat_info === nothing ? fill(NaN, n_pix, n_scan) :
              _to_pixels_scans_2d(_as_float64_nan(Array(lat_info.var)), lat_info.dims)
        lon = lon_info === nothing ? fill(NaN, n_pix, n_scan) :
              _to_pixels_scans_2d(_as_float64_nan(Array(lon_info.var)), lon_info.dims)

        sif = nothing
        if (info = _find_var(ds, "sif_radiance_678nm")) !== nothing
            sif = _to_pixels_scans_2d(_as_float64_nan(Array(info.var)), info.dims)
        end
        conv = nothing
        if (info = _find_var(ds, "converged")) !== nothing
            conv = Array(info.var)
            try
                conv = _to_pixels_scans_2d(_as_float64_nan(conv), info.dims)
            catch
                conv = reshape(Float64.(conv), n_pix, n_scan)
            end
        end

        nflh = _subset_to_retrieval(
            _read_l2_field(l2aop_path, "nflh"), ds; n_pix = n_pix, n_scan = n_scan,
        )
        chlor = fill(NaN, n_pix, n_scan)
        if l2bgc_path !== nothing && isfile(l2bgc_path)
            chlor = _subset_to_retrieval(
                _read_l2_field(l2bgc_path, "chlor_a"), ds; n_pix = n_pix, n_scan = n_scan,
            )
        end

        valid_nflh = isfinite.(nflh)
        n_valid = count(valid_nflh)
        n_skip = n_pix * n_scan - n_valid

        gain = fill(NaN, n_pix, n_scan, n_λ)
        @inbounds for j in 1:n_scan, i in 1:n_pix
            valid_nflh[i, j] || continue
            for k in 1:n_λ
                r = Rtoa[i, j, k]
                (!isfinite(r) || abs(r) < 1e-12) && continue
                res = residual[i, j, k]
                isfinite(res) || continue
                gain[i, j, k] = 1.0 + res / r
            end
        end

        return _write_vcal_gain(
            retrieval_path, l2aop_path, l2bgc_path, l1b_path, output_dir, granule_id,
            λ, lat, lon, gain, nflh, chlor, sif, conv,
            n_valid, n_skip, residual_source, svd_config_used,
        )
    end
end

function _write_vcal_gain(
    retrieval_path, l2aop_path, l2bgc_path, l1b_path, output_dir, granule_id,
    λ, lat, lon, gain, nflh, chlor, sif, conv,
    n_valid, n_skip, residual_source, svd_config_used;
    svd_solar_eff_note::AbstractString = "",
)
    n_pix, n_scan, n_λ = size(gain)
    stem = replace(splitext(basename(retrieval_path))[1], r"^PACE_OCI\." => "")
    out_path = joinpath(output_dir, "vcal_gain_$(stem).nc")
    isfile(out_path) && rm(out_path)
    NCDataset(out_path, "c") do out
        defDim(out, "pixels", n_pix)
        defDim(out, "scans", n_scan)
        defDim(out, "wavelength", n_λ)

        defVar(out, "wavelength", Float64, ("wavelength",); attrib = Dict(
            "units" => "nm", "long_name" => "retrieval band center"))
        defVar(out, "latitude", Float32, ("pixels", "scans"); attrib = Dict(
            "units" => "degrees_north"))
        defVar(out, "longitude", Float32, ("pixels", "scans"); attrib = Dict(
            "units" => "degrees_east"))
        defVar(out, "vcal_gain", Float32, ("pixels", "scans", "wavelength"); attrib = Dict(
            "long_name" => "vicarious gain = 1 + residual/Rtoa",
            "comment" => "NaN where nflh missing/non-finite or Rtoa≈0"),
            fillvalue = Float32(NaN))
        defVar(out, "nflh", Float32, ("pixels", "scans"); attrib = Dict(
            "units" => "W m-2 sr-1 um-1",
            "long_name" => "L2 AOP normalized fluorescence line height"),
            fillvalue = Float32(NaN))
        defVar(out, "chlor_a", Float32, ("pixels", "scans"); attrib = Dict(
            "units" => "mg m-3", "long_name" => "L2 BGC chlorophyll-a"),
            fillvalue = Float32(NaN))

        out["wavelength"][:] = λ
        out["latitude"][:, :] = Float32.(lat)
        out["longitude"][:, :] = Float32.(lon)
        out["vcal_gain"][:, :, :] = Float32.(gain)
        out["nflh"][:, :] = Float32.(nflh)
        out["chlor_a"][:, :] = Float32.(chlor)

        if sif !== nothing
            defVar(out, "sif_radiance_678nm", Float32, ("pixels", "scans");
                fillvalue = Float32(NaN))
            out["sif_radiance_678nm"][:, :] = Float32.(sif)
        end
        if conv !== nothing
            defVar(out, "converged", Float32, ("pixels", "scans");
                fillvalue = Float32(NaN))
            out["converged"][:, :] = Float32.(conv)
        end

        out.attrib["title"] = "PACE vicarious gain (spectral)"
        out.attrib["granule_id"] = granule_id
        out.attrib["gain_formula"] = "vcal_gain = 1 + residual/Rtoa"
        out.attrib["residual_source"] = residual_source
        out.attrib["retrieval_file"] = retrieval_path
        out.attrib["l2aop_file"] = l2aop_path
        l2bgc_path !== nothing && (out.attrib["l2bgc_file"] = l2bgc_path)
        l1b_path !== nothing && (out.attrib["l1b_file"] = l1b_path)
        !isempty(svd_config_used) && (out.attrib["svd_config"] = svd_config_used)
        !isempty(svd_solar_eff_note) && (out.attrib["svd_solar_eff_convention"] = svd_solar_eff_note)
        out.attrib["n_pixels_valid_nflh"] = n_valid
        out.attrib["n_pixels_skipped_missing_nflh"] = n_skip
        out.attrib["history"] = "Created $(Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))"
    end
    println("  wrote $out_path  valid_nflh=$n_valid  skipped_missing_nflh=$n_skip")
    return out_path
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main(args)
    length(args) >= 1 || error("usage: julia compute_vcal_gain.jl <config.toml>")
    cfg_path = args[1]
    isfile(cfg_path) || error("config not found: $cfg_path")
    cfg = TOML.parsefile(cfg_path)
    v = get(cfg, "vcal_gain", cfg)

    date = get(ENV, "DATE", String(get(v, "date", "")))
    isempty(date) && error("Set [vcal_gain].date or DATE=YYYYMMDD")
    occursin(r"^\d{8}$", date) || error("date must be YYYYMMDD, got $date")

    retrieval_dir = String(get(v, "retrieval_dir", ""))
    L1B_dir = String(get(v, "L1B_dir", ""))
    L2AOP_dir = String(get(v, "L2AOP_dir", ""))
    L2BGC_dir = String(get(v, "L2BGC_dir", ""))
    output_dir = String(get(v, "output_dir", joinpath(_VCAL_DIR, "output")))
    retrieval_glob = String(get(v, "retrieval_glob", "interim_{date}T*_svd_retrieval*.nc"))
    first_only = Bool(get(v, "first_only", false))
    svd_config = String(get(v, "svd_config", ""))

    isempty(retrieval_dir) && error("[vcal_gain].retrieval_dir required")
    isempty(L2AOP_dir) && error("[vcal_gain].L2AOP_dir required")
    isdir(retrieval_dir) || error("retrieval_dir not found: $retrieval_dir")
    isdir(L2AOP_dir) || error("L2AOP_dir not found: $L2AOP_dir")

    files, day_dir = _discover_retrievals(retrieval_dir, date, retrieval_glob)
    filter!(f -> occursin("adding_sif", basename(f)) || occursin("retrieval", basename(f)) ||
                  occursin("PACE_OCI.", basename(f)) || occursin("interim_", basename(f)), files)
    if isempty(files)
        searched = day_dir !== nothing && isdir(day_dir) ? day_dir : retrieval_dir
        error("No retrieval files for date $date under $searched (glob=$retrieval_glob). " *
              "Expected layout like \$retrieval_dir/YYYY/MM/DD/interim_\$(date)T*_svd_retrieval*.nc")
    end

    # Prefer spectral-rmse products; keep scalar-rmse SVD files (rebuild path).
    usable = String[]
    for f in files
        ok = false
        try
            NCDataset(f) do ds
                info = _find_var(ds, "rmse")
                info === nothing && return
                nd = ndims(info.var)
                if nd == 3
                    ok = true
                elseif nd == 2
                    ok = _find_var(ds, "x_hat") !== nothing
                end
            end
        catch
            ok = false
        end
        ok && push!(usable, f)
    end
    isempty(usable) && error(
        "Found $(length(files)) file(s) for $date but none have 3D spectral `rmse` " *
        "or SVD `x_hat` for residual rebuild.")
    first_only && (usable = usable[1:1])

    println("Vicarious gain for date $date")
    println("  retrievals: $(length(usable))")
    println("  output_dir: $output_dir")
    !isempty(svd_config) && println("  svd_config: $svd_config")

    written = String[]
    for f in usable
        gid = _granule_id_from_retrieval(f)
        gid === nothing && (@warn "skip (no granule id)" file = basename(f); continue)
        l2aop = _resolve_product(L2AOP_dir, gid, :l2aop)
        l2aop === nothing && (@warn "L2 AOP missing — skip" granule = gid; continue)
        l2bgc = isempty(L2BGC_dir) ? nothing : _resolve_product(L2BGC_dir, gid, :l2bgc)
        l1b = isempty(L1B_dir) ? nothing : _resolve_product(L1B_dir, gid, :l1b)
        println("  granule $gid")
        println("    retrieval: $f")
        println("    L2AOP: $l2aop")
        l2bgc !== nothing && println("    L2BGC: $l2bgc")
        l1b !== nothing && println("    L1B: $l1b")
        out = compute_vcal_gain_granule(
            f, l2aop, l2bgc, l1b, output_dir; svd_config = svd_config,
        )
        push!(written, out)
    end
    println("Done: $(length(written)) file(s)")
    return written
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
