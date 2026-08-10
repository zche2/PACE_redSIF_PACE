#!/usr/bin/env julia
# SVD transmittance swath retrieval on interim global_fit NetCDF (Rtoa_red + red_wavelength).
# Band-native grid (pseudo_measurement style): no high-res λ, no O₂/H₂O LUT or kernel regeneration.
# Parallel over pixels (Threads); thread-safe ForwardDiff Jacobians (fresh per call).
#
# Usage (internal): run_svd_orbit_full_nc_parallel(interim_nc, l1b_path, retrieval_cfg, pipeline_config_path) — see run_svd_fit.jl

using Base.Threads
using Dates
using ForwardDiff
using LinearAlgebra
using Glob
using NCDatasets
using ProgressMeter
using SparseArrays
using Statistics
using TOML

const _SVD_DIR = @__DIR__
const _PIPE_DIR = dirname(_SVD_DIR)
const _REPO_ROOT = dirname(_PIPE_DIR)
const _PIPELINE_JULIA_DIR = joinpath(_PIPE_DIR, "julia")

include(joinpath(_PIPELINE_JULIA_DIR, "sif_basis.jl"))
include(joinpath(_PIPELINE_JULIA_DIR, "pace_io.jl"))
include(joinpath(_SVD_DIR, "svd_helpers.jl"))

const STATUS_PIXEL_FILTER_SKIPPED = Int16(7)

"""Default: show progress when stdout is a TTY (off under nohup / redirected logs)."""
function _default_show_progress()
    return stdout isa Base.TTY
end

function _granule_id_from_interim_path(interim_nc::AbstractString)
    stem = splitext(basename(interim_nc))[1]
    startswith(stem, "interim_") && return stem[8:end]
    return stem
end

"""Find a variable in NetCDF4 child groups or at dataset root (PACE L1B often uses e.g. `sensor_band_parameters/red_solar_irradiance`)."""
function _find_var_in_dataset(ds, varname::String)
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
    haskey(ds, varname) || error("Variable '$varname' not found in L1B (searched child groups and root)")
    return ds[varname]
end

"""Group-aware variable lookup; returns `nothing` if absent (PACE L1B: e.g. `geolocation_data/watermask`)."""
function _find_var_in_dataset_optional(ds, varname::String)
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

function _read_watermask_2d(ds::NCDataset, varname::AbstractString, n_pix::Int, n_scan::Int)
    v = _find_var_in_dataset_optional(ds, String(varname))
    v === nothing && return fill(missing, n_pix, n_scan)
    ndims(v) == 2 || error("Expected 2D watermask variable '$varname', got ndims=$(ndims(v))")
    d = collect(String.(dimnames(v)))
    raw = v[:, :]
    arr = if d == ["pixels", "scans"]
        raw
    elseif d == ["scans", "pixels"]
        permutedims(raw, (2, 1))
    else
        error("Unsupported watermask dims for '$varname': $d")
    end
    out = Matrix{Union{Missing, Int}}(undef, size(arr)...)
    fillv = get(v.attrib, "_FillValue", nothing)
    if fillv === nothing && haskey(ds.attrib, "watermask_fill_value")
        fillv = parse(Int, String(ds.attrib["watermask_fill_value"]))
    end
    if fillv === nothing
        T = eltype(arr)
        fillv = T <: Unsigned ? typemax(T) : typemin(T)
    else
        fillv = fillv isa AbstractArray ? first(fillv) : fillv
    end
    @inbounds for j in axes(arr, 2), i in axes(arr, 1)
        a = arr[i, j]
        if ismissing(a) || a == fillv
            out[i, j] = missing
        else
            out[i, j] = Int(a)
        end
    end
    return out
end

@inline function _is_ocean_pixel(wm_val, ocean_mask_values::Set{Int})
    return !ismissing(wm_val) && Int(wm_val) in ocean_mask_values
end

# ----- NC / axis helpers (same layout as batch_fit Run_batch_full_nc) -----

function _find_axis_indices(v, wavelength_var::AbstractString)
    dnames = collect(String.(dimnames(v)))
    i_pix = findfirst(==("pixels"), dnames)
    i_scan = findfirst(==("scans"), dnames)
    i_band = findfirst(==(String(wavelength_var)), dnames)
    if isnothing(i_band)
        i_band = findfirst(x -> occursin("band", lowercase(x)), dnames)
    end
    if isnothing(i_pix) || isnothing(i_scan) || isnothing(i_band)
        error("Could not infer (pixels, scans, bands) dims from $(dnames)")
    end
    return (
        dnames = dnames,
        i_pix = i_pix,
        i_scan = i_scan,
        i_band = i_band,
        n_pix = size(v, i_pix),
        n_scan = size(v, i_scan),
        n_band = size(v, i_band),
    )
end

function _read_geo_2d(ds::NCDataset, varname::AbstractString, n_pix::Int, n_scan::Int)
    haskey(ds, varname) || return fill(Float32(NaN), n_pix, n_scan)
    v = ds[varname]
    ndims(v) == 2 || error("Expected 2D geolocation variable '$varname', got ndims=$(ndims(v))")
    d = collect(String.(dimnames(v)))
    raw = v[:, :]
    arr = if d == ["pixels", "scans"]
        raw
    elseif d == ["scans", "pixels"]
        permutedims(raw, (2, 1))
    else
        error("Unsupported geolocation dims for '$varname': $d")
    end
    out = Array{Float32}(undef, size(arr)...)
    @inbounds for j in axes(arr, 2), i in axes(arr, 1)
        a = arr[i, j]
        if ismissing(a)
            out[i, j] = Float32(NaN)
        else
            af = Float64(a)
            out[i, j] = isfinite(af) ? Float32(af) : Float32(NaN)
        end
    end
    return out
end

@inline function _copy_sorted_spectrum!(
    y_sorted::AbstractVector{Float64},
    spec_raw::AbstractVector,
    perm::AbstractVector{Int},
)
    @inbounds for i in eachindex(perm)
        v = spec_raw[perm[i]]
        if ismissing(v)
            return false
        end
        vf = Float64(v)
        # Interim Rtoa_red uses _FillValue ≈ -9999 for missing L1B inputs
        if !isfinite(vf) || vf <= -9000.0
            return false
        end
        y_sorted[i] = vf
    end
    return true
end

function _clamp_inclusive_range(r::UnitRange{Int}, n::Int)
    a = clamp(first(r), 1, n)
    b = clamp(last(r), 1, n)
    a <= b || error("Empty pixel/scan window after clamp to 1:$n (got $(first(r)):$(last(r)))")
    return a:b
end

function build_pixel_eligible_mask(
    ds::NCDataset,
    n_pix::Int,
    n_scan::Int,
    pixel_filter_vars::AbstractVector,
)
    eligible = trues(n_pix, n_scan)
    isempty(pixel_filter_vars) && return eligible
    for varname in pixel_filter_vars
        name = String(varname)
        if !haskey(ds, name)
            @warn "Pixel filter variable '$name' not in interim file; skipping this filter (merge may omit optional L2 fields)."
            continue
        end
        v = ds[name]
        ndims(v) == 2 || error("Pixel filter variable '$name' must be 2D, got ndims=$(ndims(v))")
        d = collect(String.(dimnames(v)))
        raw = v[:, :]
        arr = if d == ["pixels", "scans"]
            raw
        elseif d == ["scans", "pixels"]
            permutedims(raw, (2, 1))
        else
            error("Filter variable '$name' must have dims (pixels, scans) or (scans, pixels), got $d")
        end
        size(arr) == (n_pix, n_scan) || error(
            "Filter variable '$name' size $(size(arr)) does not match (n_pix=$n_pix, n_scan=$n_scan)",
        )
        for j in 1:n_scan, i in 1:n_pix
            eligible[i, j] = eligible[i, j] && !ismissing(arr[i, j])
        end
    end
    return eligible
end

function _create_output_dataset(
    output_path::AbstractString,
    n_pix::Int,
    n_scan::Int,
    n_sif_ev::Int,
    state_names::Vector{String},
    pace_path::AbstractString,
    config_path::AbstractString,
    pixel_range::UnitRange{Int},
    scan_range::UnitRange{Int};
    compression_level::Int = 1,
)
    ds = Dataset(output_path, "c")
    defDim(ds, "pixels", n_pix)
    defDim(ds, "scans", n_scan)
    defDim(ds, "state", length(state_names))
    defDim(ds, "sif_nev", n_sif_ev)
    ds.attrib["title"] = "PACE SVD transmittance retrieval swath output"
    ds.attrib["retrieval_type"] = "svd"
    ds.attrib["history"] = "Created " * Dates.format(now(), Dates.DateFormat("yyyy-mm-ddTHH:MM:SS"))
    ds.attrib["input_pace_file"] = String(pace_path)
    ds.attrib["config_file"] = String(config_path)
    ds.attrib["state_names_csv"] = join(state_names, ",")
    ds.attrib["pixel_start"] = first(pixel_range)
    ds.attrib["pixel_end"] = last(pixel_range)
    ds.attrib["scan_start"] = first(scan_range)
    ds.attrib["scan_end"] = last(scan_range)
    comp = (shuffle = true, deflatelevel = compression_level)
    defVar(ds, "latitude", Float32, ("pixels", "scans"); comp...)
    defVar(ds, "longitude", Float32, ("pixels", "scans"); comp...)
    defVar(ds, "x_hat", Float32, ("pixels", "scans", "state"); comp...)
    defVar(ds, "converged", UInt8, ("pixels", "scans"); comp...)
    defVar(ds, "status_code", Int16, ("pixels", "scans"); comp...)
    defVar(ds, "n_steps", Int16, ("pixels", "scans"); comp...)
    defVar(ds, "rmse", Float32, ("pixels", "scans"); comp...)
    defVar(ds, "reduced_chi2", Float32, ("pixels", "scans"); comp...)
    defVar(ds, "objective", Float32, ("pixels", "scans"); comp...)
    defVar(ds, "sif_ev1", Float32, ("pixels", "scans"); comp...)
    defVar(ds, "sif_coeffs", Float32, ("pixels", "scans", "sif_nev"); comp...)
    defVar(ds, "sif_radiance_678nm", Float32, ("pixels", "scans"); comp...)
    defVar(ds, "is_dark", UInt8, ("pixels", "scans"); comp...)
    defVar(ds, "is_ocean", UInt8, ("pixels", "scans"); comp...)
    defVar(ds, "source_pixel_index", Int32, ("pixels",))
    defVar(ds, "source_scan_index", Int32, ("scans",))
    ds["source_pixel_index"][:] = collect(Int32.(pixel_range))
    ds["source_scan_index"][:] = collect(Int32.(scan_range))
    return ds
end

function _svd_output_dir(cfg::AbstractDict)
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    out_dir_cfg = String(get(batch_cfg, "output_dir", joinpath(_REPO_ROOT, "svd_retrieval_output")))
    return isabspath(out_dir_cfg) ? out_dir_cfg : joinpath(_REPO_ROOT, out_dir_cfg)
end

function _make_svd_output_path(interim_path::AbstractString, cfg::Dict)
    out_dir = _svd_output_dir(cfg)
    mkpath(out_dir)
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    suffix = String(get(batch_cfg, "output_suffix_parallel", "_svd_retrieval_full_parallel.nc"))
    stem = splitext(basename(interim_path))[1]
    return joinpath(out_dir, stem * suffix)
end

"""Expected retrieval NetCDF for a granule (uses `[batch_fit].output_dir` and `output_suffix_parallel`)."""
function svd_expected_output_path(granule_id::AbstractString, interim_dir::AbstractString, cfg::AbstractDict)
    interim_path = joinpath(interim_dir, "interim_$(granule_id).nc")
    return _make_svd_output_path(interim_path, cfg)
end

"""
Find an existing retrieval file for `granule_id` under `[batch_fit].output_dir`.
Matches `interim_<granule_id>_*.nc` (any suffix). Returns `nothing` if none.
"""
function svd_find_existing_retrieval(granule_id::AbstractString, cfg::AbstractDict)
    out_dir = _svd_output_dir(cfg)
    isdir(out_dir) || return nothing
    pattern = "interim_$(granule_id)_*.nc"
    files = Glob.glob(pattern, out_dir)
    isempty(files) && return nothing
    return String(sort(files)[1])
end

function _read_l1b_extras(
    l1b_path::AbstractString,
    perm::Vector{Int},
    n_pix::Int,
    n_scan::Int,
)
    ds = Dataset(l1b_path)
    try
        v_e0 = _find_var_in_dataset(ds, "red_solar_irradiance")
        E0 = vec(Float64.(v_e0[:]))
        length(E0) == length(perm) || error("L1B red_solar_irradiance length $(length(E0)) != perm length $(length(perm))")
        E0_sorted = E0[perm]
        es = Float64(ds.attrib["earth_sun_distance_correction"])
        v_sza = _find_var_in_dataset(ds, "solar_zenith")
        sza_raw = v_sza[:]
        sza = Float64.(replace(sza_raw, missing => NaN))
        sza_arr = if ndims(sza) == 3
            dropdims(sza, dims = 3)
        elseif ndims(sza) == 2
            collect(sza)
        else
            error("Unexpected solar_zenith ndims=$(ndims(sza))")
        end
        if size(sza_arr) == (n_scan, n_pix)
            sza_arr = collect(sza_arr')
        end
        size(sza_arr) == (n_pix, n_scan) ||
            error("solar_zenith size $(size(sza_arr)) does not match interim (pixels=$n_pix, scans=$n_scan)")
        return E0_sorted, es, sza_arr
    finally
        close(ds)
    end
end

mutable struct SvdSwathShared
    λ_ctx::Vector{Float64}
    PCs::Matrix{Float64}
    sif_basis::Matrix{Float64}
    sif_basis_678::Vector{Float64}
    leg_basis::Matrix{Float64}
    n_pc::Int
    n_legendre::Int
    log_transform::Bool
    layout
    x0::Vector{Float64}
    prior_sigma_template::Vector{Float64}
    alpha_mean::Float64
    alpha_sigma::Float64
    use_leg01::Bool
    leg01_frac::Float64
    use_leghig::Bool
    leg_higher_sigma::Float64
    pc_prior_mode::String
    pc_sigma_scale::Float64
    svd_S::Vector{Float64}
    n_profiles::Int
    prior_min_sigma::Float64
    sif_sigma::Float64
    prior_sigma_default::Float64
    A01::Matrix{Float64}
    lower::Vector{Float64}
    upper::Vector{Float64}
    x_scale_template::Vector{Float64}
    use_band_snr::Bool
    band_snr_coeffs
    meas_sigma::Float64
    lm::NamedTuple
    conv::NamedTuple
    solar_band_ctx::Vector{Float64}
    es::Float64
    sza::Array{Float64,2}
    slice_sorted::Vector{Int}
    perm::Vector{Int}
    state_names::Vector{String}
end

function _build_svd_swath_shared(interim_nc::AbstractString, l1b_path::AbstractString, retrieval_cfg::AbstractDict)
    cfg = retrieval_cfg
    ds_i = Dataset(interim_nc)
    wl_raw = Float64.(ds_i["red_wavelength"][:])
    v_spec = ds_i["Rtoa_red"]
    ax = _find_axis_indices(v_spec, "red_wavelength")
    n_pix_i, n_scan_i = ax.n_pix, ax.n_scan
    close(ds_i)

    perm = sortperm(wl_raw)
    λ_sorted = wl_raw[perm]
    spectral_cfg = get(cfg, "spectral", Dict{String, Any}())
    λ_min = Float64(get(spectral_cfg, "lambda_min_nm", 640.0))
    λ_max = Float64(get(spectral_cfg, "lambda_max_nm", 756.0))
    slice_sorted = findall(x -> λ_min < x < λ_max, λ_sorted)
    isempty(slice_sorted) &&
        error("No red bands with λ_min < λ < λ_max ($(λ_min), $(λ_max)) nm after sorting interim red_wavelength")
    λ_ctx = λ_sorted[slice_sorted]

    E0_sorted, es, sza = _read_l1b_extras(l1b_path, perm, n_pix_i, n_scan_i)
    solar_band_ctx = E0_sorted[slice_sorted]

    fit_cfg = get(cfg, "fit", Dict{String, Any}())
    svd_cfg = get(fit_cfg, "svd", Dict{String, Any}())
    isempty(svd_cfg) && error("Retrieval config must contain [fit.svd] (n_pc, n_legendre, …); transmittance NetCDFs in [data].summer_nc / winter_nc (or legacy [fit.svd]).")
    data_cfg = get(cfg, "data", Dict{String, Any}())
    # Legacy: allow transmittance paths under [fit.svd] if [data] omits them
    data_cfg_eff = Dict{String, Any}(data_cfg)
    if !haskey(data_cfg_eff, "summer_nc") && haskey(svd_cfg, "summer_nc")
        data_cfg_eff["summer_nc"] = svd_cfg["summer_nc"]
    end
    if !haskey(data_cfg_eff, "winter_nc") && haskey(svd_cfg, "winter_nc")
        data_cfg_eff["winter_nc"] = svd_cfg["winter_nc"]
    end
    paths = resolve_svd_data_paths(data_cfg_eff)
    (paths.summer_nc !== nothing && paths.winter_nc !== nothing) ||
        error("Set [data].summer_nc and [data].winter_nc (or legacy [fit.svd] summer_nc / winter_nc)")
    paths.sif_path === nothing &&
        error("Set [data].sif_file (JLD2 SIF basis); SVD path does not use LUT/kernel preparation.")
    summer_abs = paths.summer_nc
    winter_abs = paths.winter_nc
    sif_path = paths.sif_path
    isfile(sif_path) || error("SIF file not found: $sif_path")
    n_pc = Int(get(svd_cfg, "n_pc", 5))
    n_leg = Int(get(svd_cfg, "n_legendre", 3))
    log_trans = Bool(get(svd_cfg, "svd_log_transform", false))
    sif_nev = Int(get(spectral_cfg, "sif_nev", 1))
    normalize_sif = Bool(get(spectral_cfg, "normalize_sif_first_ev", true))
    svd_basis = load_svd_basis(summer_abs, winter_abs, λ_ctx; λ_min = λ_min, λ_max = λ_max, n_pc = n_pc, log_transform = log_trans)
    PCs = Float64.(svd_basis.PCs[:, 1:n_pc])
    sif_basis = load_sif_basis(sif_path, λ_ctx; nEV = sif_nev, normalize = normalize_sif)
    size(sif_basis, 1) == length(λ_ctx) ||
        error("sif_basis rows $(size(sif_basis,1)) != length(λ_ctx) $(length(λ_ctx))")
    n_ev = size(sif_basis, 2)
    layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_leg, n_ev = n_ev)
    i678 = argmin(abs.(λ_ctx .- 678.2))
    sif_basis_678 = vec(sif_basis[i678, :])

    prior_sigma_default = Float64(get(fit_cfg, "prior_sigma_default", 1e12))
    prior_min_sigma = Float64(get(fit_cfg, "prior_min_sigma", 1e-3))
    sif_sigma = Float64(get(fit_cfg, "sif_sigma", 1e12))
    alpha_mean = Float64(get(svd_cfg, "alpha_prior_mean", 1.0))
    alpha_sigma = Float64(get(svd_cfg, "alpha_prior_sigma", 0.3))
    use_leg01 = Bool(get(svd_cfg, "use_legendre01_prior", true))
    leg01_frac = Float64(get(svd_cfg, "legendre01_prior_sigma_fraction", 1.0))
    use_leghig = Bool(get(svd_cfg, "use_legendre_higher_prior", true))
    leg_higher_sigma = Float64(get(svd_cfg, "legendre_higher_sigma", 1.0))
    pc_prior_mode = String(get(svd_cfg, "pc_prior_mode", "loading_variance"))
    pc_sigma_scale = Float64(get(svd_cfg, "pc_prior_sigma_scale", 1.0))

    x0 = zeros(Float64, layout.n_state)
    x0[first(layout.idx_alpha)] = alpha_mean
    x0[first(layout.idx_legendre)] = 1.0
    prior_sigma = fill(prior_sigma_default, layout.n_state)
    prior_sigma[first(layout.idx_alpha)] = max(alpha_sigma, prior_min_sigma)
    if pc_prior_mode == "loading_variance"
        n_prof = svd_basis.n_profiles
        for k in 1:n_pc
            sigma_k = svd_basis.S[k] / sqrt(Float64(n_prof)) * pc_sigma_scale
            prior_sigma[k] = max(sigma_k, prior_min_sigma)
        end
    end
    if use_leghig && length(layout.idx_legendre) >= 3
        for j in 3:length(layout.idx_legendre)
            prior_sigma[layout.idx_legendre[j]] = max(leg_higher_sigma, prior_min_sigma)
        end
    end
    prior_sigma[layout.idx_sif] .= max(sif_sigma, prior_min_sigma)

    z = _normalized_grid(λ_ctx)
    leg_basis = Float64.(_legendre_design_matrix(z, n_leg))
    A01 = hcat(ones(length(z)), z)
    lower = fill(-Inf, layout.n_state)
    upper = fill(Inf, layout.n_state)
    x_scale = ones(Float64, layout.n_state)
    for k in 1:n_pc
        x_scale[k] = prior_sigma[k]
    end
    x_scale[first(layout.idx_alpha)] = prior_sigma[first(layout.idx_alpha)]

    kernel_cfg = get(cfg, "kernel", Dict{String, Any}())
    use_band_snr = haskey(fit_cfg, "use_band_snr") ? Bool(fit_cfg["use_band_snr"]) : Bool(get(kernel_cfg, "use_band_snr", true))
    meas_sigma = Float64(get(fit_cfg, "meas_sigma", 0.01))
    pace_snr_path = paths.pace_snr_path
    if pace_snr_path === nothing
        # Default relative name under base_dir when key omitted
        pace_snr_path = resolve_data_path("PACE_OCI_L1BLUT_baseline_SNR_1.1.txt", paths.base_dir)
    end
    band_snr_coeffs = if use_band_snr
        load_pace_band_snr_coeffs(pace_snr_path, λ_ctx; λ_min = λ_min, λ_max = λ_max)
    else
        nothing
    end

    lm = (
        lambda0 = Float64(get(fit_cfg, "lm_lambda0", 1.0)),
        lambda_up = Float64(get(fit_cfg, "lm_lambda_up", 5.0)),
        lambda_down = Float64(get(fit_cfg, "lm_lambda_down", 0.7)),
        lambda_min = Float64(get(fit_cfg, "lm_lambda_min", 1e-8)),
        lambda_max = Float64(get(fit_cfg, "lm_lambda_max", 1e8)),
        max_inner = Int(get(fit_cfg, "lm_max_inner", 24)),
    )
    conv = (
        dx_rel_tol = Float64(get(fit_cfg, "conv_dx_rel_tol", 1e-6)),
        rmse_rel_tol = Float64(get(fit_cfg, "conv_rmse_rel_tol", 1e-6)),
        rmse_abs_tol = Float64(get(fit_cfg, "conv_rmse_abs_tol", 1e-6)),
        enabled = Bool(get(fit_cfg, "conv_stall_enable", true)),
        window = Int(get(fit_cfg, "conv_stall_window", 3)),
        redchi2_target = Float64(get(fit_cfg, "conv_stall_redchi2_target", 5.0)),
        redchi2_abs_tol = Float64(get(fit_cfg, "conv_stall_redchi2_abs_tol", 0.1)),
        redchi2_rel_tol = Float64(get(fit_cfg, "conv_stall_redchi2_rel_tol", 0.03)),
        dx_rel_tol_stall = Float64(get(fit_cfg, "conv_stall_dx_rel_tol", 5e-3)),
    )

    names = svd_state_names(layout)
    println(
        "SVD band-only setup (no LUT): ",
        length(λ_ctx),
        " red bands in (λ_min, λ_max) = (",
        λ_min,
        ", ",
        λ_max,
        ") nm; ",
        length(perm),
        " raw bands sorted.",
    )
    return SvdSwathShared(
        λ_ctx, PCs, sif_basis, sif_basis_678, leg_basis, n_pc, n_leg, log_trans, layout, x0, prior_sigma,
        alpha_mean, alpha_sigma, use_leg01, leg01_frac, use_leghig, leg_higher_sigma,
        pc_prior_mode, pc_sigma_scale, collect(Float64.(svd_basis.S)), svd_basis.n_profiles,
        prior_min_sigma, sif_sigma, prior_sigma_default, A01, lower, upper, x_scale,
        use_band_snr, band_snr_coeffs, meas_sigma, lm, conv, solar_band_ctx, es, sza, slice_sorted, perm,
        names,
    )
end

function _process_one_pixel_svd!(
    i_pix_out::Int,
    i_pix_src::Int,
    j_scan_src::Int,
    sh::SvdSwathShared,
    buf,
    slab,
    slab_is_pix_band::Bool,
    eligible,
    watermask,
    ocean_filter_enabled::Bool,
    ocean_mask_values::Set{Int},
    dark_filter_enabled::Bool,
    dark_max_radiance::Float64,
    max_outer_steps::Int,
    n_sif_ev::Int,
    state_scan,
    conv_scan,
    status_scan,
    steps_scan,
    rmse_scan,
    rchi2_scan,
    obj_scan,
    sif1_scan,
    sif_coeffs_scan,
    sif_678_scan,
    dark_scan,
    ocean_scan,
)
    if !eligible[i_pix_src, j_scan_src]
        status_scan[i_pix_out] = STATUS_PIXEL_FILTER_SKIPPED
        return
    end
    wm_val = watermask[i_pix_src, j_scan_src]
    is_ocean = _is_ocean_pixel(wm_val, ocean_mask_values)
    ocean_scan[i_pix_out] = is_ocean ? UInt8(1) : UInt8(0)
    if ocean_filter_enabled && !is_ocean
        status_scan[i_pix_out] = Int16(6)
        return
    end
    spec_raw = slab_is_pix_band ? view(slab, i_pix_src, :) : view(slab, :, i_pix_src)
    ok = _copy_sorted_spectrum!(buf.y_sorted, spec_raw, sh.perm)
    if !ok
        status_scan[i_pix_out] = Int16(3)
        return
    end
    buf.y_obs .= buf.y_sorted[sh.slice_sorted]
    is_dark = maximum(buf.y_obs) <= dark_max_radiance
    dark_scan[i_pix_out] = is_dark ? UInt8(1) : UInt8(0)
    if dark_filter_enabled && !is_dark
        status_scan[i_pix_out] = Int16(5)
        return
    end
    cosz = sh.sza[i_pix_src, j_scan_src]
    if !isfinite(cosz)
        status_scan[i_pix_out] = Int16(3)
        return
    end
    solar_eff = @. sh.solar_band_ctx * cosd(cosz) / π / sh.es
    fm, layout = make_svd_forward_model_λ(
        sh.λ_ctx,
        solar_eff,
        sh.PCs,
        sh.sif_basis;
        n_pc = sh.n_pc,
        n_legendre = sh.n_legendre,
        log_transform = sh.log_transform,
    )
    jac_eval = (x -> ForwardDiff.jacobian(fm, x))
    x_a = copy(sh.x0)
    σ_prior = copy(sh.prior_sigma_template)
    if sh.use_leg01 && length(layout.idx_legendre) >= 1
        y0 = fm(sh.x0)
        ratio = buf.y_obs ./ max.(abs.(y0), eps(Float64))
        z = sh.A01[:, 2]
        w = buf.y_obs .- minimum(buf.y_obs)
        w .+= max(maximum(w), 1.0) * 1e-6
        sv = sqrt.(w ./ maximum(w))
        c01 = (sh.A01 .* sv) \ (ratio .* sv)
        leg0 = first(layout.idx_legendre)
        x_a[leg0] = c01[1]
        σ_prior[leg0] = max(abs(c01[1]) * sh.leg01_frac, sh.prior_min_sigma)
        if length(layout.idx_legendre) >= 2
            leg1 = layout.idx_legendre[2]
            x_a[leg1] = c01[2]
            σ_prior[leg1] = max(abs(c01[2]) * sh.leg01_frac, sh.prior_min_sigma)
        end
    end
    stats = _run_one_svd_retrieval!(
        buf.x_tmp,
        fm,
        jac_eval,
        layout,
        buf.y_obs,
        copy(sh.x0),
        x_a,
        σ_prior,
        sh.lower,
        sh.upper,
        sh.x_scale_template,
        sh.use_leg01,
        sh.A01,
        sh.prior_min_sigma,
        sh.leg01_frac,
        sh.use_band_snr,
        sh.band_snr_coeffs,
        sh.meas_sigma,
        sh.lm,
        (
            dx_rel_tol = sh.conv.dx_rel_tol,
            rmse_rel_tol = sh.conv.rmse_rel_tol,
            rmse_abs_tol = sh.conv.rmse_abs_tol,
            enabled = sh.conv.enabled,
            window = sh.conv.window,
            redchi2_target = sh.conv.redchi2_target,
            redchi2_abs_tol = sh.conv.redchi2_abs_tol,
            redchi2_rel_tol = sh.conv.redchi2_rel_tol,
        ),
        sh.conv.dx_rel_tol_stall,
        max_outer_steps,
    )
    state_scan[i_pix_out, :] .= Float32.(buf.x_tmp)
    conv_scan[i_pix_out] = stats.converged ? UInt8(1) : UInt8(0)
    status_scan[i_pix_out] = stats.status
    steps_scan[i_pix_out] = Int16(stats.n_steps)
    rmse_scan[i_pix_out] = Float32(stats.rmse)
    rchi2_scan[i_pix_out] = Float32(stats.reduced_chi2)
    obj_scan[i_pix_out] = Float32(stats.objective)
    sif_coeff = buf.x_tmp[layout.idx_sif]
    if length(sif_coeff) >= 1
        sif1_scan[i_pix_out] = Float32(sif_coeff[1])
    end
    sif_coeffs_scan[i_pix_out, :] .= Float32.(sif_coeff)
    sif_678_scan[i_pix_out] = Float32(dot(sh.sif_basis_678, sif_coeff))
    return
end

function run_svd_orbit_full_nc_parallel(
    interim_nc::AbstractString,
    l1b_path::AbstractString,
    retrieval_cfg::AbstractDict,
    pipeline_config_path::AbstractString;
    pixel_range::Union{Nothing,UnitRange{Int}} = nothing,
    scan_range::Union{Nothing,UnitRange{Int}} = nothing,
)
    sh = _build_svd_swath_shared(interim_nc, l1b_path, retrieval_cfg)
    cfg = retrieval_cfg
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    use_threads = Bool(get(batch_cfg, "use_threads", true))
    n_threads = use_threads ? Threads.nthreads() : 1
    # Buffers indexed by threadid(); maxthreadid() can exceed nthreads() (e.g. spare thread slots).
    n_buf = use_threads ? max(1, Base.Threads.maxthreadid()) : 1
    prefetch_input = Bool(get(batch_cfg, "prefetch_input", true))
    compression_level = Int(get(batch_cfg, "output_compression_level", 1))
    deferred_write = Bool(get(batch_cfg, "deferred_write", true))
    show_progress = Bool(get(batch_cfg, "show_progress", _default_show_progress()))

    wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
    spectrum_var = String(get(pace_cfg, "spectrum_var", "Rtoa_red"))
    ds = Dataset(interim_nc)
    haskey(ds, wavelength_var) || error("Missing wavelength variable '$wavelength_var' in $interim_nc")
    haskey(ds, spectrum_var) || error("Missing spectrum variable '$spectrum_var' in $interim_nc")
    v_spec = ds[spectrum_var]
    axes_info = _find_axis_indices(v_spec, wavelength_var)
    n_pix = axes_info.n_pix
    n_scan = axes_info.n_scan

    pixel_range = pixel_range === nothing ? (1:n_pix) : _clamp_inclusive_range(pixel_range, n_pix)
    scan_range = scan_range === nothing ? (1:n_scan) : _clamp_inclusive_range(scan_range, n_scan)

    max_outer_steps = Int(get(batch_cfg, "max_outer_steps", 10))
    dark_filter_enabled = Bool(get(batch_cfg, "dark_filter_enabled", false))
    dark_max_radiance = Float64(get(batch_cfg, "dark_max_radiance", 20.0))
    ocean_filter_enabled = Bool(get(batch_cfg, "ocean_filter_enabled", false))
    watermask_var = String(get(batch_cfg, "watermask_var", "watermask"))
    ocean_mask_values_raw = get(batch_cfg, "ocean_mask_values", Any[1])
    ocean_mask_values = Set{Int}(Int(v) for v in ocean_mask_values_raw)
    isempty(ocean_mask_values) && error("batch_fit.ocean_mask_values must contain at least one value")
    pixel_filter_vars = get(batch_cfg, "pixel_filter_vars", String[])
    pixel_filter_vars = isa(pixel_filter_vars, AbstractVector) ?
        String.(pixel_filter_vars) : String[String(pixel_filter_vars)]
    eligible = build_pixel_eligible_mask(ds, n_pix, n_scan, pixel_filter_vars)
    lat = _read_geo_2d(ds, "latitude", n_pix, n_scan)
    lon = _read_geo_2d(ds, "longitude", n_pix, n_scan)
    watermask = _read_watermask_2d(ds, watermask_var, n_pix, n_scan)
    if ocean_filter_enabled && _find_var_in_dataset_optional(ds, watermask_var) === nothing
        error("Ocean filter enabled but variable '$watermask_var' is missing in $interim_nc")
    end

    n_sif_ev = sh.layout.n_ev
    output_path = _make_svd_output_path(interim_nc, cfg)
    ds_out = _create_output_dataset(
        output_path,
        length(pixel_range),
        length(scan_range),
        n_sif_ev,
        sh.state_names,
        interim_nc,
        pipeline_config_path,
        pixel_range,
        scan_range;
        compression_level = compression_level,
    )
    ds_out["latitude"][:, :] = lat[pixel_range, scan_range]
    ds_out["longitude"][:, :] = lon[pixel_range, scan_range]

    n_band = length(sh.perm)
    buffers = [
        (
            y_sorted = zeros(Float64, n_band),
            y_obs = zeros(Float64, length(sh.λ_ctx)),
            x_tmp = zeros(Float64, sh.layout.n_state),
        ) for _ in 1:n_buf
    ]

    println("SVD full-swath retrieval: ", interim_nc)
    println("  threads: ", n_threads)
    println("  pixels: ", first(pixel_range), ":", last(pixel_range), " (", length(pixel_range), ")")
    println("  scans:  ", first(scan_range), ":", last(scan_range), " (", length(scan_range), ")")
    println("  output: ", output_path)
    println("  prefetch_input: ", prefetch_input, "  deferred_write: ", deferred_write, "  compression_level: ", compression_level)

    # --- I/O acceleration: prefetch full spectrum array ---
    # Replaces ~N_scan serial NCDatasets reads with a single bulk read.
    spectra_prefetched = if prefetch_input
        t_pf = @elapsed begin
            raw = v_spec[:, :, :]
            raw isa Array ? raw : Array(raw)
        end
        @info "Input prefetch done" size_MB=round(sizeof(raw)/1e6; digits=1) t_s=round(t_pf; digits=2)
        raw
    else
        nothing
    end

    # --- I/O acceleration: deferred write buffers ---
    # Accumulate all scan outputs in memory; write once after the loop.
    n_pix_out  = length(pixel_range)
    n_scan_out = length(scan_range)
    n_state    = sh.layout.n_state
    buf_state    = deferred_write ? fill(Float32(NaN), n_pix_out, n_scan_out, n_state) : nothing
    buf_conv     = deferred_write ? fill(UInt8(0),     n_pix_out, n_scan_out)          : nothing
    buf_status   = deferred_write ? fill(Int16(3),     n_pix_out, n_scan_out)          : nothing
    buf_steps    = deferred_write ? fill(Int16(0),     n_pix_out, n_scan_out)          : nothing
    buf_rmse     = deferred_write ? fill(Float32(NaN), n_pix_out, n_scan_out)          : nothing
    buf_rchi2    = deferred_write ? fill(Float32(NaN), n_pix_out, n_scan_out)          : nothing
    buf_obj      = deferred_write ? fill(Float32(NaN), n_pix_out, n_scan_out)          : nothing
    buf_sif1     = deferred_write ? fill(Float32(NaN), n_pix_out, n_scan_out)          : nothing
    buf_sifcoeff = deferred_write ? fill(Float32(NaN), n_pix_out, n_scan_out, n_sif_ev) : nothing
    buf_sif678   = deferred_write ? fill(Float32(NaN), n_pix_out, n_scan_out)          : nothing
    buf_dark     = deferred_write ? fill(UInt8(0),     n_pix_out, n_scan_out)          : nothing
    buf_ocean    = deferred_write ? fill(UInt8(0),     n_pix_out, n_scan_out)          : nothing

    t_read_total    = 0.0
    t_compute_total = 0.0
    t_write_total   = 0.0

    gid = _granule_id_from_interim_path(interim_nc)
    scan_prog = if show_progress
        Progress(length(scan_range); desc = "scans $gid")
    else
        nothing
    end

    for (j_scan_out, j_scan_src) in enumerate(scan_range)
        inds = Any[Colon() for _ in 1:3]
        inds[axes_info.i_scan] = j_scan_src
        t_read_total += if prefetch_input
            @elapsed slab = spectra_prefetched[inds...]
        else
            @elapsed slab = v_spec[inds...]
        end
        slab_is_pix_band = size(slab) == (axes_info.n_pix, axes_info.n_band)
        slab_is_band_pix = size(slab) == (axes_info.n_band, axes_info.n_pix)
        (slab_is_pix_band || slab_is_band_pix) ||
            error("Unexpected slab size $(size(slab)) for scan=$j_scan_src")

        state_scan = fill(Float32(NaN), length(pixel_range), sh.layout.n_state)
        conv_scan = fill(UInt8(0), length(pixel_range))
        status_scan = fill(Int16(3), length(pixel_range))
        steps_scan = fill(Int16(0), length(pixel_range))
        rmse_scan = fill(Float32(NaN), length(pixel_range))
        rchi2_scan = fill(Float32(NaN), length(pixel_range))
        obj_scan = fill(Float32(NaN), length(pixel_range))
        sif1_scan = fill(Float32(NaN), length(pixel_range))
        sif_coeffs_scan = fill(Float32(NaN), length(pixel_range), n_sif_ev)
        sif_678_scan = fill(Float32(NaN), length(pixel_range))
        dark_scan = fill(UInt8(0), length(pixel_range))
        ocean_scan = fill(UInt8(0), length(pixel_range))

        t_compute_total += @elapsed if use_threads && n_threads > 1
            Threads.@threads for i_pix_out in 1:length(pixel_range)
                tid = Threads.threadid()
                i_pix_src = pixel_range[i_pix_out]
                b = buffers[tid]
                _process_one_pixel_svd!(
                    i_pix_out,
                    i_pix_src,
                    j_scan_src,
                    sh,
                    b,
                    slab,
                    slab_is_pix_band,
                    eligible,
                    watermask,
                    ocean_filter_enabled,
                    ocean_mask_values,
                    dark_filter_enabled,
                    dark_max_radiance,
                    max_outer_steps,
                    n_sif_ev,
                    state_scan,
                    conv_scan,
                    status_scan,
                    steps_scan,
                    rmse_scan,
                    rchi2_scan,
                    obj_scan,
                    sif1_scan,
                    sif_coeffs_scan,
                    sif_678_scan,
                    dark_scan,
                    ocean_scan,
                )
            end
        else
            for (i_pix_out, i_pix_src) in enumerate(pixel_range)
                _process_one_pixel_svd!(
                    i_pix_out,
                    i_pix_src,
                    j_scan_src,
                    sh,
                    buffers[1],
                    slab,
                    slab_is_pix_band,
                    eligible,
                    watermask,
                    ocean_filter_enabled,
                    ocean_mask_values,
                    dark_filter_enabled,
                    dark_max_radiance,
                    max_outer_steps,
                    n_sif_ev,
                    state_scan,
                    conv_scan,
                    status_scan,
                    steps_scan,
                    rmse_scan,
                    rchi2_scan,
                    obj_scan,
                    sif1_scan,
                    sif_coeffs_scan,
                    sif_678_scan,
                    dark_scan,
                    ocean_scan,
                )
            end
        end  # @elapsed compute

        t_write_total += @elapsed if deferred_write
            buf_state[:, j_scan_out, :]    = state_scan
            buf_conv[:, j_scan_out]        = conv_scan
            buf_status[:, j_scan_out]      = status_scan
            buf_steps[:, j_scan_out]       = steps_scan
            buf_rmse[:, j_scan_out]        = rmse_scan
            buf_rchi2[:, j_scan_out]       = rchi2_scan
            buf_obj[:, j_scan_out]         = obj_scan
            buf_sif1[:, j_scan_out]        = sif1_scan
            buf_sifcoeff[:, j_scan_out, :] = sif_coeffs_scan
            buf_sif678[:, j_scan_out]      = sif_678_scan
            buf_dark[:, j_scan_out]        = dark_scan
            buf_ocean[:, j_scan_out]       = ocean_scan
        else
            ds_out["x_hat"][:, j_scan_out, :]         = state_scan
            ds_out["converged"][:, j_scan_out]         = conv_scan
            ds_out["status_code"][:, j_scan_out]       = status_scan
            ds_out["n_steps"][:, j_scan_out]           = steps_scan
            ds_out["rmse"][:, j_scan_out]              = rmse_scan
            ds_out["reduced_chi2"][:, j_scan_out]      = rchi2_scan
            ds_out["objective"][:, j_scan_out]         = obj_scan
            ds_out["sif_ev1"][:, j_scan_out]           = sif1_scan
            ds_out["sif_coeffs"][:, j_scan_out, :]     = sif_coeffs_scan
            ds_out["sif_radiance_678nm"][:, j_scan_out] = sif_678_scan
            ds_out["is_dark"][:, j_scan_out]           = dark_scan
            ds_out["is_ocean"][:, j_scan_out]          = ocean_scan
        end
        if scan_prog !== nothing
            next!(scan_prog; showvalues = [(:scan, j_scan_src)])
        end
    end
    if scan_prog !== nothing
        finish!(scan_prog)
    end

    # Bulk write all accumulated results in one pass (only when deferred_write=true).
    if deferred_write
        t_bulk_write = @elapsed begin
            ds_out["x_hat"][:, :, :]            = buf_state
            ds_out["converged"][:, :]            = buf_conv
            ds_out["status_code"][:, :]          = buf_status
            ds_out["n_steps"][:, :]              = buf_steps
            ds_out["rmse"][:, :]                 = buf_rmse
            ds_out["reduced_chi2"][:, :]         = buf_rchi2
            ds_out["objective"][:, :]            = buf_obj
            ds_out["sif_ev1"][:, :]              = buf_sif1
            ds_out["sif_coeffs"][:, :, :]        = buf_sifcoeff
            ds_out["sif_radiance_678nm"][:, :]   = buf_sif678
            ds_out["is_dark"][:, :]              = buf_dark
            ds_out["is_ocean"][:, :]             = buf_ocean
        end
        @info "Bulk write done" t_s=round(t_bulk_write; digits=2)
        t_write_total += t_bulk_write
    end

    t_total = t_read_total + t_compute_total + t_write_total
    safe_total = max(t_total, eps())
    @info "Scan loop timing breakdown" scans=length(scan_range) t_read_s=round(t_read_total; digits=2) t_compute_s=round(t_compute_total; digits=2) t_write_s=round(t_write_total; digits=2) pct_read=round(100*t_read_total/safe_total; digits=1) pct_compute=round(100*t_compute_total/safe_total; digits=1) pct_write=round(100*t_write_total/safe_total; digits=1)
    close(ds_out)
    close(ds)
    println("Saved SVD retrieval to: ", output_path)
    return output_path
end

"""
    retrieve_svd_single_pixel(interim_nc, l1b_path, retrieval_cfg, pixel, scan)

CPU retrieval for one (pixel, scan). Returns a NamedTuple with λ, y_obs, y_mod,
residual, sif_toa (TOA SIF contribution), scalars, and x_hat.
"""
function retrieve_svd_single_pixel(
    interim_nc::AbstractString,
    l1b_path::AbstractString,
    retrieval_cfg::AbstractDict,
    pixel::Int,
    scan::Int,
)
    sh = _build_svd_swath_shared(interim_nc, l1b_path, retrieval_cfg)
    cfg = retrieval_cfg
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
    spectrum_var = String(get(pace_cfg, "spectrum_var", "Rtoa_red"))
    max_outer_steps = Int(get(batch_cfg, "max_outer_steps", 10))
    dark_filter_enabled = Bool(get(batch_cfg, "dark_filter_enabled", false))
    dark_max_radiance = Float64(get(batch_cfg, "dark_max_radiance", 20.0))
    ocean_filter_enabled = Bool(get(batch_cfg, "ocean_filter_enabled", false))
    watermask_var = String(get(batch_cfg, "watermask_var", "watermask"))
    ocean_mask_values = Set{Int}(Int(v) for v in get(batch_cfg, "ocean_mask_values", Any[1]))
    pixel_filter_vars = get(batch_cfg, "pixel_filter_vars", String[])
    pixel_filter_vars = isa(pixel_filter_vars, AbstractVector) ?
        String.(pixel_filter_vars) : String[String(pixel_filter_vars)]

    ds = Dataset(interim_nc)
    v_spec = ds[spectrum_var]
    axes_info = _find_axis_indices(v_spec, wavelength_var)
    n_pix, n_scan = axes_info.n_pix, axes_info.n_scan
    (1 <= pixel <= n_pix) || error("pixel=$pixel outside 1:$n_pix")
    (1 <= scan <= n_scan) || error("scan=$scan outside 1:$n_scan")

    eligible = build_pixel_eligible_mask(ds, n_pix, n_scan, pixel_filter_vars)
    watermask = _read_watermask_2d(ds, watermask_var, n_pix, n_scan)
    lat = _read_geo_2d(ds, "latitude", n_pix, n_scan)
    lon = _read_geo_2d(ds, "longitude", n_pix, n_scan)

    inds = Any[Colon() for _ in 1:3]
    inds[axes_info.i_scan] = scan
    slab = v_spec[inds...]
    slab_is_pix_band = size(slab) == (axes_info.n_pix, axes_info.n_band)
    slab_is_band_pix = size(slab) == (axes_info.n_band, axes_info.n_pix)
    (slab_is_pix_band || slab_is_band_pix) || error("Unexpected slab size $(size(slab))")

    buf = (
        y_sorted = zeros(Float64, length(sh.perm)),
        y_obs = zeros(Float64, length(sh.λ_ctx)),
        x_tmp = zeros(Float64, sh.layout.n_state),
    )
    state_scan = fill(Float32(NaN), 1, sh.layout.n_state)
    conv_scan = fill(UInt8(0), 1)
    status_scan = fill(Int16(3), 1)
    steps_scan = fill(Int16(0), 1)
    rmse_scan = fill(Float32(NaN), 1)
    rchi2_scan = fill(Float32(NaN), 1)
    obj_scan = fill(Float32(NaN), 1)
    sif1_scan = fill(Float32(NaN), 1)
    sif_coeffs_scan = fill(Float32(NaN), 1, sh.layout.n_ev)
    sif_678_scan = fill(Float32(NaN), 1)
    dark_scan = fill(UInt8(0), 1)
    ocean_scan = fill(UInt8(0), 1)

    _process_one_pixel_svd!(
        1, pixel, scan, sh, buf, slab, slab_is_pix_band, eligible, watermask,
        ocean_filter_enabled, ocean_mask_values, dark_filter_enabled, dark_max_radiance,
        max_outer_steps, sh.layout.n_ev,
        state_scan, conv_scan, status_scan, steps_scan, rmse_scan, rchi2_scan, obj_scan,
        sif1_scan, sif_coeffs_scan, sif_678_scan, dark_scan, ocean_scan,
    )
    close(ds)

    x_hat = Float64.(vec(state_scan[1, :]))
    y_obs = copy(buf.y_obs)
    layout = sh.layout
    cosz = sh.sza[pixel, scan]
    solar_eff = @. sh.solar_band_ctx * cosd(cosz) / π / sh.es
    fm, _ = make_svd_forward_model_λ(
        sh.λ_ctx, solar_eff, sh.PCs, sh.sif_basis;
        n_pc = sh.n_pc, n_legendre = sh.n_legendre, log_transform = sh.log_transform,
    )
    y_mod = fill(NaN, length(y_obs))
    sif_toa = fill(NaN, length(y_obs))
    sif_shape = fill(NaN, length(y_obs))
    if status_scan[1] in (Int16(0), Int16(1), Int16(2)) && all(isfinite, x_hat)
        y_mod .= fm(x_hat)
        c_vec = @view x_hat[layout.idx_pc]
        sif_coeff = @view x_hat[layout.idx_sif]
        trans_up = sh.log_transform ? exp.(sh.PCs * c_vec) : 1.0 .+ sh.PCs * c_vec
        sif_shape .= sh.sif_basis * sif_coeff
        sif_toa .= trans_up .* sif_shape
    end
    residual = y_obs .- y_mod

    return (
        λ = copy(sh.λ_ctx),
        y_obs = y_obs,
        y_mod = y_mod,
        residual = residual,
        sif_toa = sif_toa,
        sif_shape = sif_shape,
        x_hat = x_hat,
        status_code = Int(status_scan[1]),
        converged = conv_scan[1] == 0x01,
        n_steps = Int(steps_scan[1]),
        rmse = Float64(rmse_scan[1]),
        reduced_chi2 = Float64(rchi2_scan[1]),
        objective = Float64(obj_scan[1]),
        sif_radiance_678nm = Float64(sif_678_scan[1]),
        latitude = Float64(lat[pixel, scan]),
        longitude = Float64(lon[pixel, scan]),
        pixel = pixel,
        scan = scan,
        interim_nc = String(interim_nc),
        l1b_path = String(l1b_path),
    )
end
