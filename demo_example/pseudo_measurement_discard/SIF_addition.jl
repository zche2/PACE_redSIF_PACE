#!/usr/bin/env julia
# Add known SIF radiance to real PACE data, run retrieval, and check if the
# pseudo SIF is recovered. Supports scene-mode (L1B+L2, parallel) and
# single-pixel mode. SIF shapes from jld2, magnitude 0.5, SNR-consistent noise.

using Base.Threads
using TOML
using LinearAlgebra
using SparseArrays
using Plots
using JLD2
using Interpolations
using NCDatasets
using Glob

const PSEUDO_DIR = @__DIR__
const DEMO_DIR = joinpath(PSEUDO_DIR, "..")
const GLOBAL_FIT_DIR = joinpath(DEMO_DIR, "global_fit")
const BATCH_FIT_DIR = joinpath(DEMO_DIR, "batch_fit")

include(joinpath(DEMO_DIR, "Fit_toy_forward_model.jl"))
include(joinpath(GLOBAL_FIT_DIR, "pre_process.jl"))
include(joinpath(GLOBAL_FIT_DIR, "merge_inputs.jl"))
include(joinpath(BATCH_FIT_DIR, "Run_batch_full_nc.jl"))

const MWEF = SimplePACEXSecFitMWEFunctions
const STATUS_PIXEL_FILTER_SKIPPED = Int16(7)

# Subset/merge constants (from granule_retrieval) for use_in_memory_merge=false path
const L1B_RENAME_DIMS = Dict(
    "pixels_per_line" => "pixels",
    "number_of_lines" => "scans",
    "number_of_bands" => "red_wavelength",
    "red_bands" => "red_wavelength",
)
const L2_RENAME_DIMS = Dict(
    "pixels_per_line" => "pixels",
    "number_of_lines" => "scans",
)
const L1B_VARS = [
    "red_wavelength",
    "red_solar_irradiance",
    "watermask",
    "latitude",
    "longitude",
    "rhot_red",
]
const L2AOP_VARS = ["nflh"]
const L2BGC_VARS = ["chlor_a"]

using Dates

# ----- Output dataset for SIF addition (with and without SIF retrievals) -----
function _create_sif_addition_output_dataset(
    output_path::AbstractString,
    n_pix::Int,
    n_scan::Int,
    n_sif_ev::Int,
    state_names::Vector{String},
    pace_path::AbstractString,
    config_path::AbstractString,
    pixel_range::AbstractRange{Int},
    scan_range::AbstractRange{Int},
)
    ds = Dataset(output_path, "c")
    defDim(ds, "pixels", n_pix)
    defDim(ds, "scans", n_scan)
    defDim(ds, "state", length(state_names))
    defDim(ds, "sif_nev", n_sif_ev)
    ds.attrib["title"] = "PACE SIF addition experiment (with and without SIF retrievals)"
    ds.attrib["history"] = "Created " * Dates.format(now(), Dates.DateFormat("yyyy-mm-ddTHH:MM:SS"))
    ds.attrib["input_pace_file"] = String(pace_path)
    ds.attrib["config_file"] = String(config_path)
    ds.attrib["state_names_csv"] = join(state_names, ",")
    ds.attrib["pixel_start"] = first(pixel_range)
    ds.attrib["pixel_end"] = last(pixel_range)
    ds.attrib["scan_start"] = first(scan_range)
    ds.attrib["scan_end"] = last(scan_range)
    comp = (shuffle=true, deflatelevel=4)
    defVar(ds, "latitude", Float32, ("pixels", "scans"); comp...)
    defVar(ds, "longitude", Float32, ("pixels", "scans"); comp...)
    defVar(ds, "sif_added_678nm", Float32, ("pixels", "scans"); comp...)
    for suffix in ["with_sif", "without_sif"]
        defVar(ds, "x_hat_$(suffix)", Float32, ("pixels", "scans", "state"); comp...)
        defVar(ds, "converged_$(suffix)", UInt8, ("pixels", "scans"); comp...)
        defVar(ds, "status_code_$(suffix)", Int16, ("pixels", "scans"); comp...)
        defVar(ds, "n_steps_$(suffix)", Int16, ("pixels", "scans"); comp...)
        defVar(ds, "rmse_$(suffix)", Float32, ("pixels", "scans"); comp...)
        defVar(ds, "reduced_chi2_$(suffix)", Float32, ("pixels", "scans"); comp...)
        defVar(ds, "objective_$(suffix)", Float32, ("pixels", "scans"); comp...)
        defVar(ds, "sif_ev1_$(suffix)", Float32, ("pixels", "scans"); comp...)
        defVar(ds, "sif_coeffs_$(suffix)", Float32, ("pixels", "scans", "sif_nev"); comp...)
        defVar(ds, "sif_radiance_678nm_$(suffix)", Float32, ("pixels", "scans"); comp...)
    end
    defVar(ds, "is_dark", UInt8, ("pixels", "scans"); comp...)
    defVar(ds, "is_ocean", UInt8, ("pixels", "scans"); comp...)
    defVar(ds, "source_pixel_index", Int32, ("pixels",))
    defVar(ds, "source_scan_index", Int32, ("scans",))
    ds["source_pixel_index"][:] = collect(Int32.(pixel_range))
    ds["source_scan_index"][:] = collect(Int32.(scan_range))
    ds["sif_added_678nm"].attrib["long_name"] = "Added SIF radiance at 678 nm (truth)"
    return ds
end

# ----- Config merge (from Run_global_fit) -----
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

# ----- Granule path resolution (from global_fit) -----
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

# ----- SIF shapes at lres -----
"""
    load_sif_shapes_lres(ctx, sif_path; magnitude=0.5)

Load SIF_shapes from jld2, interpolate each to ctx.λ (lres), scale to given magnitude.
Returns Vector{Vector{Float64}} of length n_shapes.
Magnitude: scale by value at 678 nm (or max if 678 not in range).
"""
function load_sif_shapes_lres(ctx, sif_path::AbstractString; magnitude::Float64=0.5)
    sif = JLD2.load(sif_path)
    haskey(sif, "SIF_shapes") && haskey(sif, "SIF_wavelen") || error("SIF_shapes and SIF_wavelen required in $sif_path")
    shapes_raw = Float64.(sif["SIF_shapes"])
    λ_ref = collect(Float64.(sif["SIF_wavelen"]))
    # shapes_raw is (n_wavelength × n_samples)
    n_bands = size(shapes_raw, 1)
    n_shapes = size(shapes_raw, 2)
    n_bands == length(λ_ref) || error("SIF_shapes first dim must match SIF_wavelen")

    λ_target = collect(Float64.(ctx.λ))
    idx_678 = argmin(abs.(λ_target .- 678.2))
    sif_shapes_lres = Vector{Vector{Float64}}(undef, n_shapes)

    for k in 1:n_shapes
        itp = LinearInterpolation(λ_ref, shapes_raw[:, k], extrapolation_bc=0.0)
        s_lres = Float64.(itp.(λ_target))
        # Scale so magnitude at 678 nm = magnitude (or max if 678 not in range)
        s_val = abs(s_lres[idx_678] > 0 ? s_lres[idx_678] : maximum(abs.(s_lres)))
        s_val > 0 || (s_val = 1.0)
        scale = magnitude / s_val
        sif_shapes_lres[k] = s_lres .* scale
    end
    return sif_shapes_lres
end

"""
    add_sif_radiance_to_spectrum(y_obs, ctx, sif_coeff_add)

Add convolved SIF radiance to an observed spectrum: y_mod = y_obs + K * (sif_basis_hres * coeff).
Returns y_obs + SIF_at_OCI (both on ctx.λ). Used for single-pixel basis mode.
"""
function add_sif_radiance_to_spectrum(
    y_obs::AbstractVector{<:Real},
    ctx,
    sif_coeff_add::AbstractVector{<:Real},
)
    length(y_obs) == length(ctx.λ) || error("y_obs length must match ctx.λ")
    length(sif_coeff_add) == size(ctx.sif_basis_hres, 2) || error("sif_coeff_add length must match SIF basis columns")
    sif_hres = ctx.sif_basis_hres * sif_coeff_add
    sif_at_oci = ctx.kernel_rsr_out * sif_hres
    return y_obs .+ sif_at_oci
end

"""
    add_sif_lres_and_noise!(y_obs, sif_lres, band_snr_coeffs, add_noise; meas_sigma=0.01)

Add sif_lres to y_obs. If add_noise, add Gaussian noise with sigma from band_snr_coeffs.
Sigma computed from sif_lres (the signal we add); y_obs already has its own noise.
Modifies y_obs in place.
"""
function add_sif_lres_and_noise!(
    y_obs::Vector{Float64},
    sif_lres::Vector{Float64},
    band_snr_coeffs::Union{Dict, Nothing},
    add_noise::Bool;
    meas_sigma::Float64=0.01,
)
    length(y_obs) == length(sif_lres) || error("y_obs length must match sif_lres")
    n = length(y_obs)
    sigma = if add_noise && n > 0 && band_snr_coeffs !== nothing
        c1 = band_snr_coeffs["c1"]
        c2 = band_snr_coeffs["c2"]
        length(c1) == n || error("band_snr c1 length mismatch")
        [sqrt(max(Float64(c1[i]) + Float64(c2[i]) * sif_lres[i], 1e-20)) for i in 1:n]
    elseif add_noise && n > 0
        fill(meas_sigma, n)
    else
        Float64[]
    end
    y_obs .+= sif_lres
    if add_noise && !isempty(sigma)
        y_obs .+= randn(n) .* sigma
    end
    return y_obs
end

# ----- Process one pixel (two retrievals: with and without SIF) -----
function _process_one_pixel_sif_addition!(
    i_pix_out::Int,
    i_pix_src::Int,
    j_scan_src::Int,
    core,
    buf,
    slab,
    slab_is_pix_band::Bool,
    perm,
    W_interp,
    eligible,
    watermask,
    ocean_filter_enabled::Bool,
    ocean_mask_values::Set{Int},
    dark_filter_enabled::Bool,
    dark_max_radiance::Float64,
    max_outer_steps::Int,
    n_sif_ev::Int,
    state_scan_with_sif,
    conv_scan_with_sif,
    status_scan_with_sif,
    steps_scan_with_sif,
    rmse_scan_with_sif,
    rchi2_scan_with_sif,
    obj_scan_with_sif,
    sif1_scan_with_sif,
    sif_coeffs_scan_with_sif,
    sif_678_scan_with_sif,
    state_scan_without_sif,
    conv_scan_without_sif,
    status_scan_without_sif,
    steps_scan_without_sif,
    rmse_scan_without_sif,
    rchi2_scan_without_sif,
    obj_scan_without_sif,
    sif1_scan_without_sif,
    sif_coeffs_scan_without_sif,
    sif_678_scan_without_sif,
    sif_added_678_scan,
    dark_scan,
    ocean_scan,
    sif_shapes_lres::Vector{Vector{Float64}},
    add_snr_noise::Bool,
    n_pix::Int,
    meas_sigma::Float64,
)
    if !eligible[i_pix_src, j_scan_src]
        status_scan_with_sif[i_pix_out] = STATUS_PIXEL_FILTER_SKIPPED
        status_scan_without_sif[i_pix_out] = STATUS_PIXEL_FILTER_SKIPPED
        return
    end
    spec_raw = slab_is_pix_band ? view(slab, i_pix_src, :) : view(slab, :, i_pix_src)
    ok = _copy_sorted_spectrum!(buf.y_sorted, spec_raw, perm)
    if !ok
        status_scan_with_sif[i_pix_out] = Int16(3)
        status_scan_without_sif[i_pix_out] = Int16(3)
        return
    end
    wm_val = watermask[i_pix_src, j_scan_src]
    is_ocean = !ismissing(wm_val) && (Int(wm_val) in ocean_mask_values)
    ocean_scan[i_pix_out] = is_ocean ? UInt8(1) : UInt8(0)
    if ocean_filter_enabled && !is_ocean
        status_scan_with_sif[i_pix_out] = Int16(6)
        status_scan_without_sif[i_pix_out] = Int16(6)
        return
    end
    mul!(buf.y_obs, W_interp, buf.y_sorted)
    is_dark = maximum(buf.y_obs) <= dark_max_radiance
    dark_scan[i_pix_out] = is_dark ? UInt8(1) : UInt8(0)
    if dark_filter_enabled && !is_dark
        status_scan_with_sif[i_pix_out] = Int16(5)
        status_scan_without_sif[i_pix_out] = Int16(5)
        return
    end

    # Retrieval 1: without SIF (original spectrum)
    stats = _run_one_retrieval!(buf.x_tmp, core, buf.y_obs, max_outer_steps)
    state_scan_without_sif[i_pix_out, :] .= Float32.(buf.x_tmp)
    conv_scan_without_sif[i_pix_out] = stats.converged ? UInt8(1) : UInt8(0)
    status_scan_without_sif[i_pix_out] = stats.status
    steps_scan_without_sif[i_pix_out] = Int16(stats.n_steps)
    rmse_scan_without_sif[i_pix_out] = Float32(stats.rmse)
    rchi2_scan_without_sif[i_pix_out] = Float32(stats.reduced_chi2)
    obj_scan_without_sif[i_pix_out] = Float32(stats.objective)
    sif_coeff = buf.x_tmp[core.layout.idx_sif]
    if length(sif_coeff) >= 1
        sif1_scan_without_sif[i_pix_out] = Float32(sif_coeff[1])
    end
    sif_coeffs_scan_without_sif[i_pix_out, :] .= Float32.(sif_coeff)
    sif_678_scan_without_sif[i_pix_out] = Float32(dot(core.sif_basis_678, sif_coeff))

    # Add SIF + noise, then retrieval 2: with SIF
    linear_idx = (j_scan_src - 1) * n_pix + i_pix_src
    shape_idx = (linear_idx - 1) % length(sif_shapes_lres) + 1
    sif_lres = sif_shapes_lres[shape_idx]
    idx_678 = argmin(abs.(collect(core.ctx.λ) .- 678.2))
    sif_added_678_scan[i_pix_out] = Float32(sif_lres[idx_678])
    add_sif_lres_and_noise!(
        buf.y_obs,
        sif_lres,
        core.ctx.band_snr_coeffs,
        add_snr_noise;
        meas_sigma=meas_sigma,
    )

    stats = _run_one_retrieval!(buf.x_tmp, core, buf.y_obs, max_outer_steps)
    state_scan_with_sif[i_pix_out, :] .= Float32.(buf.x_tmp)
    conv_scan_with_sif[i_pix_out] = stats.converged ? UInt8(1) : UInt8(0)
    status_scan_with_sif[i_pix_out] = stats.status
    steps_scan_with_sif[i_pix_out] = Int16(stats.n_steps)
    rmse_scan_with_sif[i_pix_out] = Float32(stats.rmse)
    rchi2_scan_with_sif[i_pix_out] = Float32(stats.reduced_chi2)
    obj_scan_with_sif[i_pix_out] = Float32(stats.objective)
    sif_coeff = buf.x_tmp[core.layout.idx_sif]
    if length(sif_coeff) >= 1
        sif1_scan_with_sif[i_pix_out] = Float32(sif_coeff[1])
    end
    sif_coeffs_scan_with_sif[i_pix_out, :] .= Float32.(sif_coeff)
    sif_678_scan_with_sif[i_pix_out] = Float32(dot(core.sif_basis_678, sif_coeff))
    return
end

# ----- Scene-mode: run SIF addition + retrieval on full swath -----
function run_sif_addition_scene(
    cores::Vector,
    pace_path::AbstractString,
    config_path::AbstractString,
    sif_shapes_lres::Vector{Vector{Float64}},
    add_snr_noise::Bool,
    nflh_nonzero::Bool,
)
    core = cores[1]
    cfg = core.cfg
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    fit_cfg = get(cfg, "fit", Dict{String, Any}())
    use_threads = Bool(get(batch_cfg, "use_threads", true))
    n_threads = use_threads ? Threads.nthreads() : 1
    meas_sigma = Float64(get(fit_cfg, "meas_sigma", 0.01))

    wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
    spectrum_var = String(get(pace_cfg, "spectrum_var", "Rtoa_red"))
    ds = Dataset(pace_path)
    haskey(ds, wavelength_var) || error("Missing wavelength variable '$wavelength_var' in $pace_path")
    haskey(ds, spectrum_var) || error("Missing spectrum variable '$spectrum_var' in $pace_path")
    v_spec = ds[spectrum_var]
    axes_info = _find_axis_indices(v_spec, wavelength_var)
    λ_src = collect(Float64.(ds[wavelength_var][:]))
    W_interp, perm = _make_linear_resampler(λ_src, core.ctx.λ)
    n_pix = axes_info.n_pix
    n_scan = axes_info.n_scan
    pixel_range = 1:n_pix
    scan_range = 1:n_scan
    sif_add = get(cfg, "sif_addition", Dict{String, Any}())
    scan_start = get(sif_add, "scan_start", nothing)
    scan_end = get(sif_add, "scan_end", nothing)
    if scan_start !== nothing && scan_end !== nothing
        try
            s1 = Int(scan_start)
            s2 = Int(scan_end)
            scan_range = max(1, s1):min(n_scan, s2)
        catch
            @warn "Invalid scan_start/scan_end, using full swath"
        end
    end
    max_outer_steps = Int(get(batch_cfg, "max_outer_steps", core.lm.max_outer_default))
    dark_filter_enabled = Bool(get(batch_cfg, "dark_filter_enabled", false))
    dark_max_radiance = Float64(get(batch_cfg, "dark_max_radiance", 20.0))
    ocean_filter_enabled = Bool(get(batch_cfg, "ocean_filter_enabled", false))
    watermask_var = String(get(batch_cfg, "watermask_var", "watermask"))
    ocean_mask_values_raw = get(batch_cfg, "ocean_mask_values", Any[1])
    ocean_mask_values = Set{Int}(Int(v) for v in ocean_mask_values_raw)
    isempty(ocean_mask_values) && error("batch_fit.ocean_mask_values must contain at least one value")
    pixel_filter_vars = get(batch_cfg, "pixel_filter_vars", ["nflh"])
    pixel_filter_vars = isa(pixel_filter_vars, AbstractVector) ? String.(pixel_filter_vars) : String[String(pixel_filter_vars)]
    eligible = build_pixel_eligible_mask(ds, n_pix, n_scan, pixel_filter_vars)
    if nflh_nonzero && haskey(ds, "nflh")
        nflh = ds["nflh"][:, :]
        d = collect(String.(dimnames(ds["nflh"])))
        arr = d == ["pixels", "scans"] ? nflh : permutedims(nflh, (2, 1))
        for j in 1:n_scan, i in 1:n_pix
            v = arr[i, j]
            if eligible[i, j] && (ismissing(v) || Float64(v) == 0.0)
                eligible[i, j] = false
            end
        end
    end

    lat = _read_geo_2d(ds, "latitude", n_pix, n_scan)
    lon = _read_geo_2d(ds, "longitude", n_pix, n_scan)
    watermask = if haskey(ds, watermask_var)
        wm = ds[watermask_var]
        ndims(wm) == 2 || error("Expected 2D watermask variable '$watermask_var', got ndims=$(ndims(wm))")
        d = collect(String.(dimnames(wm)))
        raw = wm[:, :]
        d == ["pixels", "scans"] ? raw : permutedims(raw, (2, 1))
    else
        fill(missing, n_pix, n_scan)
    end

    output_path = joinpath(dirname(pace_path), splitext(basename(pace_path))[1] * "_sif_addition_retrieval.nc")
    batch_out = get(batch_cfg, "output_dir", nothing)
    if batch_out !== nothing && !isempty(String(batch_out))
        out_dir = isabspath(batch_out) ? batch_out : joinpath(DEMO_DIR, batch_out)
        mkpath(out_dir)
        output_path = joinpath(out_dir, splitext(basename(pace_path))[1] * "_sif_addition_retrieval.nc")
    end
    n_sif_ev = length(core.layout.idx_sif)
    ds_out = _create_sif_addition_output_dataset(
        output_path,
        length(pixel_range),
        length(scan_range),
        n_sif_ev,
        core.state_names,
        pace_path,
        config_path,
        pixel_range,
        scan_range,
    )
    ds_out["latitude"][:, :] = lat[pixel_range, scan_range]
    ds_out["longitude"][:, :] = lon[pixel_range, scan_range]

    n_eligible = count(eligible)
    println("SIF addition scene: ", pace_path)
    println("  threads: ", n_threads)
    println("  pixels: ", first(pixel_range), ":", last(pixel_range), " (", length(pixel_range), ")")
    println("  scans:  ", first(scan_range), ":", last(scan_range), " (", length(scan_range), ")")
    println("  eligible (nflh != 0): ", n_eligible)
    println("  SIF shapes: ", length(sif_shapes_lres), " | add_snr_noise: ", add_snr_noise)

    buffers = [(
        y_sorted = zeros(Float64, length(perm)),
        y_obs = zeros(Float64, length(core.ctx.λ)),
        x_tmp = zeros(Float64, core.layout.n_state),
    ) for _ in 1:n_threads]

    for (j_scan_out, j_scan_src) in enumerate(scan_range)
        inds = Any[Colon() for _ in 1:3]
        inds[axes_info.i_scan] = j_scan_src
        slab = v_spec[inds...]
        slab_is_pix_band = size(slab) == (axes_info.n_pix, axes_info.n_band)
        slab_is_band_pix = size(slab) == (axes_info.n_band, axes_info.n_pix)
        (slab_is_pix_band || slab_is_band_pix) || error("Unexpected slab size $(size(slab)) for scan=$j_scan_src")

        state_scan_with_sif = fill(Float32(NaN), length(pixel_range), core.layout.n_state)
        conv_scan_with_sif = fill(UInt8(0), length(pixel_range))
        status_scan_with_sif = fill(Int16(3), length(pixel_range))
        steps_scan_with_sif = fill(Int16(0), length(pixel_range))
        rmse_scan_with_sif = fill(Float32(NaN), length(pixel_range))
        rchi2_scan_with_sif = fill(Float32(NaN), length(pixel_range))
        obj_scan_with_sif = fill(Float32(NaN), length(pixel_range))
        sif1_scan_with_sif = fill(Float32(NaN), length(pixel_range))
        sif_coeffs_scan_with_sif = fill(Float32(NaN), length(pixel_range), n_sif_ev)
        sif_678_scan_with_sif = fill(Float32(NaN), length(pixel_range))
        state_scan_without_sif = fill(Float32(NaN), length(pixel_range), core.layout.n_state)
        conv_scan_without_sif = fill(UInt8(0), length(pixel_range))
        status_scan_without_sif = fill(Int16(3), length(pixel_range))
        steps_scan_without_sif = fill(Int16(0), length(pixel_range))
        rmse_scan_without_sif = fill(Float32(NaN), length(pixel_range))
        rchi2_scan_without_sif = fill(Float32(NaN), length(pixel_range))
        obj_scan_without_sif = fill(Float32(NaN), length(pixel_range))
        sif1_scan_without_sif = fill(Float32(NaN), length(pixel_range))
        sif_coeffs_scan_without_sif = fill(Float32(NaN), length(pixel_range), n_sif_ev)
        sif_678_scan_without_sif = fill(Float32(NaN), length(pixel_range))
        sif_added_678_scan = fill(Float32(NaN), length(pixel_range))
        dark_scan = fill(UInt8(0), length(pixel_range))
        ocean_scan = fill(UInt8(0), length(pixel_range))

        if use_threads && n_threads > 1
            Threads.@threads for i_pix_out in 1:length(pixel_range)
                tid = Threads.threadid()
                i_pix_src = pixel_range[i_pix_out]
                c = cores[tid]
                b = buffers[tid]
                _process_one_pixel_sif_addition!(
                    i_pix_out, i_pix_src, j_scan_src, c, b, slab, slab_is_pix_band,
                    perm, W_interp, eligible, watermask,
                    ocean_filter_enabled, ocean_mask_values,
                    dark_filter_enabled, dark_max_radiance, max_outer_steps, n_sif_ev,
                    state_scan_with_sif, conv_scan_with_sif, status_scan_with_sif, steps_scan_with_sif,
                    rmse_scan_with_sif, rchi2_scan_with_sif, obj_scan_with_sif,
                    sif1_scan_with_sif, sif_coeffs_scan_with_sif, sif_678_scan_with_sif,
                    state_scan_without_sif, conv_scan_without_sif, status_scan_without_sif, steps_scan_without_sif,
                    rmse_scan_without_sif, rchi2_scan_without_sif, obj_scan_without_sif,
                    sif1_scan_without_sif, sif_coeffs_scan_without_sif, sif_678_scan_without_sif,
                    sif_added_678_scan, dark_scan, ocean_scan,
                    sif_shapes_lres, add_snr_noise, n_pix, meas_sigma,
                )
            end
        else
            for (i_pix_out, i_pix_src) in enumerate(pixel_range)
                _process_one_pixel_sif_addition!(
                    i_pix_out, i_pix_src, j_scan_src, cores[1], buffers[1], slab, slab_is_pix_band,
                    perm, W_interp, eligible, watermask,
                    ocean_filter_enabled, ocean_mask_values,
                    dark_filter_enabled, dark_max_radiance, max_outer_steps, n_sif_ev,
                    state_scan_with_sif, conv_scan_with_sif, status_scan_with_sif, steps_scan_with_sif,
                    rmse_scan_with_sif, rchi2_scan_with_sif, obj_scan_with_sif,
                    sif1_scan_with_sif, sif_coeffs_scan_with_sif, sif_678_scan_with_sif,
                    state_scan_without_sif, conv_scan_without_sif, status_scan_without_sif, steps_scan_without_sif,
                    rmse_scan_without_sif, rchi2_scan_without_sif, obj_scan_without_sif,
                    sif1_scan_without_sif, sif_coeffs_scan_without_sif, sif_678_scan_without_sif,
                    sif_added_678_scan, dark_scan, ocean_scan,
                    sif_shapes_lres, add_snr_noise, n_pix, meas_sigma,
                )
            end
        end

        ds_out["x_hat_with_sif"][:, j_scan_out, :] = state_scan_with_sif
        ds_out["converged_with_sif"][:, j_scan_out] = conv_scan_with_sif
        ds_out["status_code_with_sif"][:, j_scan_out] = status_scan_with_sif
        ds_out["n_steps_with_sif"][:, j_scan_out] = steps_scan_with_sif
        ds_out["rmse_with_sif"][:, j_scan_out] = rmse_scan_with_sif
        ds_out["reduced_chi2_with_sif"][:, j_scan_out] = rchi2_scan_with_sif
        ds_out["objective_with_sif"][:, j_scan_out] = obj_scan_with_sif
        ds_out["sif_ev1_with_sif"][:, j_scan_out] = sif1_scan_with_sif
        ds_out["sif_coeffs_with_sif"][:, j_scan_out, :] = sif_coeffs_scan_with_sif
        ds_out["sif_radiance_678nm_with_sif"][:, j_scan_out] = sif_678_scan_with_sif
        ds_out["x_hat_without_sif"][:, j_scan_out, :] = state_scan_without_sif
        ds_out["converged_without_sif"][:, j_scan_out] = conv_scan_without_sif
        ds_out["status_code_without_sif"][:, j_scan_out] = status_scan_without_sif
        ds_out["n_steps_without_sif"][:, j_scan_out] = steps_scan_without_sif
        ds_out["rmse_without_sif"][:, j_scan_out] = rmse_scan_without_sif
        ds_out["reduced_chi2_without_sif"][:, j_scan_out] = rchi2_scan_without_sif
        ds_out["objective_without_sif"][:, j_scan_out] = obj_scan_without_sif
        ds_out["sif_ev1_without_sif"][:, j_scan_out] = sif1_scan_without_sif
        ds_out["sif_coeffs_without_sif"][:, j_scan_out, :] = sif_coeffs_scan_without_sif
        ds_out["sif_radiance_678nm_without_sif"][:, j_scan_out] = sif_678_scan_without_sif
        ds_out["sif_added_678nm"][:, j_scan_out] = sif_added_678_scan
        ds_out["is_dark"][:, j_scan_out] = dark_scan
        ds_out["is_ocean"][:, j_scan_out] = ocean_scan
        n_conv_ws = count(==(UInt8(1)), conv_scan_with_sif)
        n_conv_wo = count(==(UInt8(1)), conv_scan_without_sif)
        n_skip = count(==(STATUS_PIXEL_FILTER_SKIPPED), status_scan_with_sif)
        println("  scan ", j_scan_src, " -> filter_skip ", n_skip, " | converged with_sif ", n_conv_ws, " | without_sif ", n_conv_wo, "/", length(pixel_range))
    end
    close(ds_out)
    close(ds)
    println("Saved SIF addition retrieval to: ", output_path)
    return output_path
end

# ----- One-granule SIF addition (mirror run_granule_retrieval) -----
"""
Process one granule: preprocess_and_merge_in_memory → run_sif_addition_scene.
- use_threads: if true, enable pixel-level parallelism
- keep_interim: if false, remove interim file after retrieval
- shared_config: if true, use config_path directly (no per-granule tmp file)
- use_in_memory_merge: if true, read L1B/L2 directly and write single interim
"""
function run_sif_addition_one_granule(
    L1B_path::AbstractString,
    L2AOP_path::AbstractString,
    L2BGC_path::Union{AbstractString, Nothing},
    output_dir::AbstractString,
    interim_dir::AbstractString,
    config_path::AbstractString;
    use_threads::Bool=true,
    keep_interim::Bool=false,
    shared_config::Bool=false,
    use_in_memory_merge::Bool=true,
)
    mkpath(interim_dir)
    mkpath(output_dir)

    stem = splitext(basename(L1B_path))[1]
    m = match(r"^PACE_OCI\.(.+)\.L1B\.V3$", stem)
    granule_id = m !== nothing ? m.captures[1] : stem
    interim_path = joinpath(interim_dir, "interim_$(granule_id).nc")

    subset_L1B = subset_L2AOP = subset_L2BGC = nothing
    if use_in_memory_merge
        preprocess_and_merge_in_memory(L1B_path, L2AOP_path, L2BGC_path, interim_path)
    else
        subset_L1B = subset_netcdf_dataset(
            L1B_path,
            L1B_VARS,
            interim_dir,
            rename_dims=L1B_RENAME_DIMS,
            post_process_func=TOA_radiance,
        )
        subset_L2AOP = subset_netcdf_dataset(
            L2AOP_path,
            L2AOP_VARS,
            interim_dir,
            rename_dims=L2_RENAME_DIMS,
        )
        subset_L2BGC = L2BGC_path !== nothing && isfile(L2BGC_path) ? subset_netcdf_dataset(
            L2BGC_path,
            L2BGC_VARS,
            interim_dir,
            rename_dims=L2_RENAME_DIMS,
        ) : nothing
        merge_global_fit_inputs(
            subset_L1B,
            subset_L2AOP,
            interim_path,
            L2BGC_subset_path=subset_L2BGC,
        )
    end

    cfg = TOML.parsefile(config_path)
    retrieval_config = config_path
    if !shared_config
        batch_cfg = merge(get(cfg, "batch_fit", Dict{String, Any}()), Dict("output_dir" => output_dir, "use_threads" => use_threads))
        cfg["batch_fit"] = batch_cfg
        cfg["pace_observation"] = merge(
            get(cfg, "pace_observation", Dict{String, Any}()),
            Dict("spectrum_var" => "Rtoa_red", "wavelength_var" => "red_wavelength"),
        )
        cfg["batch_fit"]["pixel_filter_vars"] = get(get(cfg, "batch_fit", Dict()), "pixel_filter_vars", ["nflh"])
        tmp_config = joinpath(interim_dir, "tmp_config_$(granule_id).toml")
        open(tmp_config, "w") do io
            TOML.print(io, cfg)
        end
        retrieval_config = tmp_config
    end

    sif_add_cfg = get(cfg, "sif_addition", Dict{String, Any}())
    data_cfg = get(cfg, "data", Dict{String, Any}())
    base_dir = get(data_cfg, "base_dir", joinpath(DEMO_DIR, "..", "Files_in_use"))
    base_dir = isabspath(base_dir) ? base_dir : normpath(joinpath(DEMO_DIR, base_dir))
    sif_magnitude = Float64(get(sif_add_cfg, "sif_magnitude", 0.5))
    sif_jld2 = get(sif_add_cfg, "sif_jld2_path", "")
    if isempty(sif_jld2)
        sif_file = String(get(data_cfg, "sif_file", "SIF_singular_vector.jld2"))
        data_base = String(get(data_cfg, "base_dir", base_dir))
        sif_jld2 = isabspath(data_base) ? joinpath(data_base, sif_file) : joinpath(DEMO_DIR, data_base, sif_file)
    else
        sif_jld2 = isabspath(sif_jld2) ? sif_jld2 : joinpath(base_dir, sif_jld2)
    end
    add_snr_noise = Bool(get(sif_add_cfg, "add_snr_noise", true))
    nflh_nonzero = Bool(get(sif_add_cfg, "nflh_nonzero", true))

    ctx = MWEF.prepare_mwe_inputs(retrieval_config)
    sif_shapes_lres = load_sif_shapes_lres(ctx, sif_jld2; magnitude=sif_magnitude)
    n_threads = use_threads ? Threads.nthreads() : 1
    cores = [_build_retrieval_core(retrieval_config) for _ in 1:n_threads]

    run_sif_addition_scene(cores, interim_path, retrieval_config, sif_shapes_lres, add_snr_noise, nflh_nonzero)

    if !shared_config
        rm(retrieval_config, force=true)
    end
    if !keep_interim
        if !use_in_memory_merge
            subset_L1B !== nothing && rm(subset_L1B, force=true)
            subset_L2AOP !== nothing && rm(subset_L2AOP, force=true)
            subset_L2BGC !== nothing && rm(subset_L2BGC, force=true)
        end
        rm(interim_path, force=true)
    end

    return output_dir
end

# ----- Mode 1: Date-based (all granules for a date) -----
function run_sif_addition_date(config_path::AbstractString)
    pseudo_cfg = TOML.parsefile(config_path)
    main_config_path = normpath(joinpath(dirname(config_path), get(get(pseudo_cfg, "pseudo_measurement", Dict()), "main_config_path", "../Simple_PACE_xSecFit_MWE_zcheVer.toml")))
    cfg_base = isfile(main_config_path) ? TOML.parsefile(main_config_path) : Dict{String, Any}()
    cfg_overlay = pseudo_cfg
    cfg = _merge_config(cfg_base, cfg_overlay)

    sif_add = get(cfg, "sif_addition", Dict{String, Any}())
    date = String(get(sif_add, "date", ""))
    isempty(date) && error("[sif_addition].date required when mode=date")
    L1B_dir = String(get(sif_add, "L1B_dir", ""))
    L2AOP_dir = String(get(sif_add, "L2AOP_dir", ""))
    L2BGC_dir = String(get(sif_add, "L2BGC_dir", ""))
    output_dir = String(get(sif_add, "output_dir", joinpath(PSEUDO_DIR, "sif_addition_output")))
    interim_dir = String(get(sif_add, "interim_dir", joinpath(PSEUDO_DIR, "sif_addition_interim")))
    parallel_granules = Bool(get(sif_add, "parallel_granules", true))
    use_in_memory_merge = Bool(get(sif_add, "use_in_memory_merge", !parallel_granules))

    L1B_dir = isabspath(L1B_dir) ? L1B_dir : joinpath(DEMO_DIR, L1B_dir)
    L2AOP_dir = isabspath(L2AOP_dir) ? L2AOP_dir : joinpath(DEMO_DIR, L2AOP_dir)
    L2BGC_dir = isabspath(L2BGC_dir) ? L2BGC_dir : joinpath(DEMO_DIR, L2BGC_dir)
    output_dir = isabspath(output_dir) ? output_dir : joinpath(DEMO_DIR, output_dir)
    interim_dir = isabspath(interim_dir) ? interim_dir : joinpath(DEMO_DIR, interim_dir)

    mkpath(interim_dir)
    batch_use_threads = !parallel_granules && Bool(get(get(cfg, "batch_fit", Dict()), "use_threads", true))
    cfg["batch_fit"] = merge(get(cfg, "batch_fit", Dict()), Dict(
        "output_dir" => output_dir,
        "use_threads" => batch_use_threads,
        "pixel_filter_vars" => get(sif_add, "pixel_filter_vars", ["nflh"]),
    ))
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

    println("SIF addition (date mode): $date")
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
            run_sif_addition_one_granule(
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
            run_sif_addition_one_granule(
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
function run_sif_addition_granule(config_path::AbstractString, granule_id::AbstractString)
    pseudo_cfg = TOML.parsefile(config_path)
    main_config_path = normpath(joinpath(dirname(config_path), get(get(pseudo_cfg, "pseudo_measurement", Dict()), "main_config_path", "../Simple_PACE_xSecFit_MWE_zcheVer.toml")))
    cfg_base = isfile(main_config_path) ? TOML.parsefile(main_config_path) : Dict{String, Any}()
    cfg_overlay = pseudo_cfg
    cfg = _merge_config(cfg_base, cfg_overlay)

    sif_add = get(cfg, "sif_addition", Dict{String, Any}())
    L1B_dir = String(get(sif_add, "L1B_dir", ""))
    L2AOP_dir = String(get(sif_add, "L2AOP_dir", ""))
    L2BGC_dir = String(get(sif_add, "L2BGC_dir", ""))
    output_dir = String(get(sif_add, "output_dir", joinpath(PSEUDO_DIR, "sif_addition_output")))
    interim_dir = String(get(sif_add, "interim_dir", joinpath(PSEUDO_DIR, "sif_addition_interim")))
    use_in_memory_merge = Bool(get(sif_add, "use_in_memory_merge", true))

    L1B_dir = isabspath(L1B_dir) ? L1B_dir : joinpath(DEMO_DIR, L1B_dir)
    L2AOP_dir = isabspath(L2AOP_dir) ? L2AOP_dir : joinpath(DEMO_DIR, L2AOP_dir)
    L2BGC_dir = isabspath(L2BGC_dir) ? L2BGC_dir : joinpath(DEMO_DIR, L2BGC_dir)
    output_dir = isabspath(output_dir) ? output_dir : joinpath(DEMO_DIR, output_dir)
    interim_dir = isabspath(interim_dir) ? interim_dir : joinpath(DEMO_DIR, interim_dir)

    mkpath(interim_dir)
    cfg["batch_fit"] = merge(get(cfg, "batch_fit", Dict()), Dict(
        "output_dir" => output_dir,
        "use_threads" => true,
        "pixel_filter_vars" => get(sif_add, "pixel_filter_vars", ["nflh"]),
    ))
    cfg["pace_observation"] = merge(get(cfg, "pace_observation", Dict()), Dict("spectrum_var" => "Rtoa_red", "wavelength_var" => "red_wavelength"))
    merged_config = joinpath(interim_dir, "merged_config_granule.toml")
    open(merged_config, "w") do io
        TOML.print(io, cfg)
    end

    L1B_path, L2AOP_path, L2BGC_path = _resolve_granule_paths(L1B_dir, L2AOP_dir, L2BGC_dir, granule_id)
    isfile(L1B_path) || error("L1B not found: $L1B_path")
    isfile(L2AOP_path) || error("L2AOP not found: $L2AOP_path")

    println("SIF addition (granule mode): $granule_id")
    run_sif_addition_one_granule(
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

# ----- Main entry: date, granule, or single-pixel mode -----
function main_sif_addition()
    pseudo_cfg_path = joinpath(PSEUDO_DIR, "pseudo_measurement_config.toml")
    pseudo_cfg = TOML.parsefile(pseudo_cfg_path)
    sif_add_cfg = get(pseudo_cfg, "sif_addition", Dict{String, Any}())
    main_config = normpath(joinpath(PSEUDO_DIR, get(get(pseudo_cfg, "pseudo_measurement", Dict()), "main_config_path", "../Simple_PACE_xSecFit_MWE_zcheVer.toml")))
    cfg = TOML.parsefile(main_config)
    data_cfg = get(cfg, "data", Dict{String, Any}())
    base_dir = get(data_cfg, "base_dir", joinpath(DEMO_DIR, "..", "Files_in_use"))
    base_dir = isabspath(base_dir) ? base_dir : normpath(joinpath(DEMO_DIR, base_dir))

    mode = String(get(sif_add_cfg, "mode", "granule"))
    mode in ("date", "granule", "single_pixel") || error("[sif_addition] mode must be 'date', 'granule', or 'single_pixel', got: $mode")
    granule_id = String(get(sif_add_cfg, "granule_id", ""))
    date = get(sif_add_cfg, "date", "")

    if mode == "date"
        run_sif_addition_date(pseudo_cfg_path)
    elseif mode == "granule"
        isempty(granule_id) && error("[sif_addition].granule_id required when mode=granule")
        run_sif_addition_granule(pseudo_cfg_path, granule_id)
    else  # mode == "single_pixel"
        pace_file = get(sif_add_cfg, "pace_file", "sample_granule_20240830T131442_new_chl.nc")
        pixel_index = Int(get(sif_add_cfg, "pixel_index", 1260))
        scan_index = Int(get(sif_add_cfg, "scan_index", 1700))
        sif_coeff_add = Float64.(get(sif_add_cfg, "sif_coeff_add", [0.5]))
        max_iter = Int(get(sif_add_cfg, "max_iter", 50))

        println("SIF_addition: add known SIF to PACE spectrum and retrieve (single-pixel)")
        println("  main_config: ", main_config)
        println("  pace_file: ", pace_file, "  pixel: ", pixel_index, "  scan: ", scan_index)
        println("  SIF coeff to add (truth): ", sif_coeff_add)

        ctx = MWEF.prepare_mwe_inputs(main_config)
        pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
        pace_path = isabspath(pace_file) ? pace_file : joinpath(base_dir, pace_file)
        wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
        spectrum_var = String(get(pace_cfg, "spectrum_var", "radiance_red"))

        y_obs, _ = MWEF.load_pace_spectrum_on_grid(
            pace_path, ctx.λ;
            pixel_idx=pixel_index, scan_idx=scan_index,
            wavelength_var=wavelength_var, spectrum_var=spectrum_var,
        )
        y_obs = Float64.(y_obs)
        println("  Loaded spectrum: length ", length(y_obs), "  mean radiance ", mean(y_obs))

        length(sif_coeff_add) == size(ctx.sif_basis_hres, 2) || error(
            "sif_coeff_add length ($(length(sif_coeff_add))) must match SIF basis columns ($(size(ctx.sif_basis_hres, 2)))"
        )

        y_modified = add_sif_radiance_to_spectrum(y_obs, ctx, sif_coeff_add)
        println("  Modified spectrum: mean radiance ", mean(y_modified))

        # Run two retrievals: without SIF (y_obs) and with SIF (y_modified)
        fit_cfg = get(cfg, "fit", Dict{String, Any}())
        solar_file = get(data_cfg, "solar_file", "solar_merged_20240731_600_33300_100.out")
        solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)
        solar_hres, _ = MWEF.load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
        solar_hres = Float64.(solar_hres)
        n_legendre = Int(get(fit_cfg, "n_legendre", 3))
        model_variant = Symbol(get(fit_cfg, "model_variant", "standard"))
        use_hybrid_jacobian = Bool(get(fit_cfg, "use_hybrid_jacobian", false))
        preallocate_forward = Bool(get(fit_cfg, "preallocate_forward", true))

        fm = make_forward_model_simple(ctx, solar_hres; n_legendre=n_legendre, preallocate_float64=preallocate_forward, model_variant=model_variant)
        layout = state_layout_simple(ctx; n_legendre=n_legendre)
        x0 = initial_state_simple(ctx; n_legendre=n_legendre, T=Float64, model_variant=model_variant)
        jacobian_eval = use_hybrid_jacobian ? make_hybrid_jacobian_evaluator(fm, ctx, solar_hres, layout; n_legendre=n_legendre) : make_jacobian_evaluator(fm, x0; use_preallocated=false)

        prior_sigma_default = Float64(get(fit_cfg, "prior_sigma_default", 1e30))
        prior_min_sigma = Float64(get(fit_cfg, "prior_min_sigma", 1e-3))
        use_legendre01_prior = Bool(get(fit_cfg, "use_legendre01_prior", true))
        legendre01_prior_sigma_fraction = Float64(get(fit_cfg, "legendre01_prior_sigma_fraction", 0.2))
        use_legendre_higher_prior = Bool(get(fit_cfg, "use_legendre_higher_prior", true))
        legendre_higher_sigma = Float64(get(fit_cfg, "legendre_higher_sigma", 1.0))
        p_prior_hpa = Float64(get(fit_cfg, "p_prior_hpa", 700.0))
        p_sigma_hpa = Float64(get(fit_cfg, "p_sigma_hpa", 200.0))
        t_prior_k = Float64(get(fit_cfg, "t_prior_k", 280.0))
        t_sigma_k = Float64(get(fit_cfg, "t_sigma_k", 20.0))
        vcd_o2_sigma = Float64(get(fit_cfg, "vcd_o2_sigma", 1e23))
        vcd_h2o_sigma = Float64(get(fit_cfg, "vcd_h2o_sigma", 3e22))
        use_vcd_slope_prior = Bool(get(fit_cfg, "use_vcd_slope_prior", true))
        vcd_slope_prior_sigma_factor = Float64(get(fit_cfg, "vcd_slope_prior_sigma_factor", 1.0))
        sif_sigma = Float64(get(fit_cfg, "sif_sigma", 1e12))
        use_pt_constraints = Bool(get(fit_cfg, "use_pt_constraints", true))
        pt_constraint_sigma_mult = Float64(get(fit_cfg, "pt_constraint_sigma_mult", 3.0))
        conv_dx_rel_tol = Float64(get(fit_cfg, "conv_dx_rel_tol", 1e-6))
        conv_rmse_rel_tol = Float64(get(fit_cfg, "conv_rmse_rel_tol", 1e-6))
        conv_rmse_abs_tol = Float64(get(fit_cfg, "conv_rmse_abs_tol", 1e-6))
        conv_stall_enable = Bool(get(fit_cfg, "conv_stall_enable", true))
        conv_stall_window = Int(get(fit_cfg, "conv_stall_window", 3))
        conv_stall_redchi2_target = Float64(get(fit_cfg, "conv_stall_redchi2_target", 5.0))
        conv_stall_redchi2_abs_tol = Float64(get(fit_cfg, "conv_stall_redchi2_abs_tol", 0.1))
        conv_stall_redchi2_rel_tol = Float64(get(fit_cfg, "conv_stall_redchi2_rel_tol", 0.03))
        conv_stall_dx_rel_tol = Float64(get(fit_cfg, "conv_stall_dx_rel_tol", 5e-3))
        lm_lambda0 = Float64(get(fit_cfg, "lm_lambda0", 1.0))
        lm_lambda_up = Float64(get(fit_cfg, "lm_lambda_up", 2.0))
        lm_lambda_down = Float64(get(fit_cfg, "lm_lambda_down", 0.7))
        lm_lambda_min = Float64(get(fit_cfg, "lm_lambda_min", 1e-8))
        lm_lambda_max = Float64(get(fit_cfg, "lm_lambda_max", 1e8))
        lm_max_inner = Int(get(fit_cfg, "lm_max_inner", 8))
        use_band_snr = Bool(get(fit_cfg, "use_band_snr", true))
        meas_sigma = Float64(get(fit_cfg, "meas_sigma", 0.01))
        n_plot_steps = Int(get(fit_cfg, "n_plot_steps", 20))

        x_a = copy(x0)
        prior_sigma = fill(prior_sigma_default, length(x0))
        if use_legendre01_prior && length(layout.idx_legendre) >= 1
            leg0_idx = first(layout.idx_legendre)
            y_base = fm(x0)
            ratio = y_modified ./ max.(abs.(y_base), eps(Float64))
            z = _normalized_grid(ctx.λ)
            A01 = hcat(ones(length(z)), z)
            w = y_modified .- minimum(y_modified)
            w .+= max(maximum(w), 1.0) * 1e-6
            s = sqrt.(w ./ maximum(w))
            c01 = (A01 .* s) \ (ratio .* s)
            x_a[leg0_idx] = c01[1]
            prior_sigma[leg0_idx] = max(abs(c01[1]) * legendre01_prior_sigma_fraction, prior_min_sigma)
            if length(layout.idx_legendre) >= 2
                leg1_idx = layout.idx_legendre[2]
                x_a[leg1_idx] = c01[2]
                prior_sigma[leg1_idx] = max(abs(c01[2]) * legendre01_prior_sigma_fraction, prior_min_sigma)
            end
        end
        if use_legendre_higher_prior && length(layout.idx_legendre) >= 3
            for j in 3:length(layout.idx_legendre)
                idx = layout.idx_legendre[j]
                x_a[idx] = 0.0
                prior_sigma[idx] = max(legendre_higher_sigma, prior_min_sigma)
            end
        end
        x_a[layout.idx_p_o2_hpa] = p_prior_hpa
        x_a[layout.idx_p_h2o_hpa] = p_prior_hpa
        x_a[layout.idx_t_o2_k] = t_prior_k
        x_a[layout.idx_t_h2o_k] = t_prior_k
        x_a[layout.idx_vcd_o2_intercept] = x0[layout.idx_vcd_o2_intercept]
        x_a[layout.idx_vcd_h2o_intercept] = x0[layout.idx_vcd_h2o_intercept]
        x_a[layout.idx_vcd_o2_sif] = x0[layout.idx_vcd_o2_sif]
        x_a[layout.idx_vcd_h2o_sif] = x0[layout.idx_vcd_h2o_sif]
        prior_sigma[layout.idx_p_o2_hpa] = p_sigma_hpa
        prior_sigma[layout.idx_p_h2o_hpa] = p_sigma_hpa
        prior_sigma[layout.idx_t_o2_k] = t_sigma_k
        prior_sigma[layout.idx_t_h2o_k] = t_sigma_k
        prior_sigma[layout.idx_vcd_o2_intercept] = vcd_o2_sigma
        prior_sigma[layout.idx_vcd_h2o_intercept] = vcd_h2o_sigma
        prior_sigma[layout.idx_vcd_o2_sif] = vcd_o2_sigma
        prior_sigma[layout.idx_vcd_h2o_sif] = vcd_h2o_sigma
        if use_vcd_slope_prior
            x_a[layout.idx_vcd_o2_slope] = 0.0
            x_a[layout.idx_vcd_h2o_slope] = 0.0
            prior_sigma[layout.idx_vcd_o2_slope] = max(vcd_o2_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
            prior_sigma[layout.idx_vcd_h2o_slope] = max(vcd_h2o_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
        end
        x_a[layout.idx_sif] .= 0.0
        prior_sigma[layout.idx_sif] .= max(sif_sigma, prior_min_sigma)

        if hasproperty(ctx, :sif_prior_cov) && !isnothing(ctx.sif_prior_cov) && length(layout.idx_sif) == size(ctx.sif_prior_cov, 1)
            σ = collect(Float64.(prior_sigma))
            @. σ = clamp(abs(σ), 1e-12, 1e100)
            S_a_inv_dense = Matrix(Diagonal(@. 1.0 / (σ^2)))
            S_a_inv_dense[layout.idx_sif, layout.idx_sif] .= inv(ctx.sif_prior_cov)
            S_a_inv = S_a_inv_dense
        else
            S_a_inv = _spdiag_invvar(prior_sigma)
        end

        x_scale = ones(Float64, length(x0))
        x_scale[layout.idx_vcd_o2_intercept] = vcd_o2_sigma
        x_scale[layout.idx_vcd_o2_slope] = max(vcd_o2_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
        x_scale[layout.idx_vcd_h2o_intercept] = vcd_h2o_sigma
        x_scale[layout.idx_vcd_h2o_slope] = max(vcd_h2o_sigma * vcd_slope_prior_sigma_factor, prior_min_sigma)
        x_scale[layout.idx_vcd_o2_sif] = vcd_o2_sigma
        x_scale[layout.idx_vcd_h2o_sif] = vcd_h2o_sigma
        x_scale[layout.idx_p_o2_hpa] = p_sigma_hpa
        x_scale[layout.idx_p_h2o_hpa] = p_sigma_hpa
        x_scale[layout.idx_t_o2_k] = t_sigma_k
        x_scale[layout.idx_t_h2o_k] = t_sigma_k
        x_scale[layout.idx_sif] .= 1.0
        x_scale[layout.idx_legendre] .= 1.0

        lower_bounds = fill(-Inf, length(x0))
        upper_bounds = fill(Inf, length(x0))
        if use_pt_constraints
            lower_bounds[layout.idx_p_o2_hpa] = p_prior_hpa - pt_constraint_sigma_mult * p_sigma_hpa
            upper_bounds[layout.idx_p_o2_hpa] = p_prior_hpa + pt_constraint_sigma_mult * p_sigma_hpa
            lower_bounds[layout.idx_p_h2o_hpa] = p_prior_hpa - pt_constraint_sigma_mult * p_sigma_hpa
            upper_bounds[layout.idx_p_h2o_hpa] = p_prior_hpa + pt_constraint_sigma_mult * p_sigma_hpa
            lower_bounds[layout.idx_t_o2_k] = t_prior_k - pt_constraint_sigma_mult * t_sigma_k
            upper_bounds[layout.idx_t_o2_k] = t_prior_k + pt_constraint_sigma_mult * t_sigma_k
            lower_bounds[layout.idx_t_h2o_k] = t_prior_k - pt_constraint_sigma_mult * t_sigma_k
            upper_bounds[layout.idx_t_h2o_k] = t_prior_k + pt_constraint_sigma_mult * t_sigma_k
        end

        _run_single_pixel_lm = (y_target) -> begin
            x_curr = copy(x_a)
            y_curr = copy(fm(x_curr))
            S_e_inv = (!use_band_snr || isnothing(ctx.band_snr_coeffs)) ? spdiagm(0 => fill(1.0 / (meas_sigma^2), length(y_target))) : make_Se_inv_from_snr(y_curr, ctx.band_snr_coeffs)
            dof = max(length(y_target) - layout.n_state, 1)
            dx_rel_series = Float64[]
            redchi2_series = Float64[]
            λ = lm_lambda0
            for istep in 1:max(n_plot_steps, 1)
                x_prev = copy(x_curr)
                rmse_prev = sqrt(mean((y_target .- y_curr) .^ 2))
                step = try
                    lm_one_step(fm, x_curr, y_target; x_a=x_a, S_a_inv=S_a_inv, lambda=λ, lambda_up=lm_lambda_up, lambda_down=lm_lambda_down, lambda_min=lm_lambda_min, lambda_max=lm_lambda_max, max_inner=lm_max_inner, jacobian_eval=jacobian_eval, x_scale=x_scale, lower_bounds=lower_bounds, upper_bounds=upper_bounds, use_band_snr=use_band_snr, band_snr_coeffs=ctx.band_snr_coeffs, meas_sigma=meas_sigma)
                catch err
                    println("  LM step failed: ", err)
                    break
                end
                λ = step.lambda_next
                if !step.accepted
                    stalled_conv, _ = _stalled_convergence(dx_rel_series, redchi2_series; enabled=conv_stall_enable, window=conv_stall_window, redchi2_target=conv_stall_redchi2_target, redchi2_abs_tol=conv_stall_redchi2_abs_tol, redchi2_rel_tol=conv_stall_redchi2_rel_tol, dx_rel_tol=conv_stall_dx_rel_tol)
                    stalled_conv && break
                end
                x_curr = step.x_next
                y_curr = copy(step.y_next)
                push!(redchi2_series, step.chi2_next / dof)
                push!(dx_rel_series, norm(step.dx) / max(norm(x_prev), eps(Float64)))
                rmse_curr = sqrt(mean((y_target .- y_curr) .^ 2))
                if norm(step.dx) / max(norm(x_prev), eps(Float64)) < conv_dx_rel_tol || abs(rmse_curr - rmse_prev) / max(abs(rmse_prev), eps(Float64)) < conv_rmse_rel_tol
                    break
                end
            end
            return x_curr
        end

        println("  Retrieval 1: without SIF (original spectrum)")
        x_without_sif = _run_single_pixel_lm(y_obs)
        println("  Retrieval 2: with SIF (modified spectrum)")
        x_with_sif = _run_single_pixel_lm(y_modified)

        sif_retrieved_with_sif = x_with_sif[layout.idx_sif]
        sif_retrieved_without_sif = x_without_sif[layout.idx_sif]
        sif_hres_with_sif = ctx.sif_basis_hres * sif_retrieved_with_sif
        sif_hres_without_sif = ctx.sif_basis_hres * sif_retrieved_without_sif
        sif_added_hres = ctx.sif_basis_hres * sif_coeff_add
        y_reconstructed_with_sif = fm(x_with_sif)
        y_reconstructed_without_sif = fm(x_without_sif)

        println("\n--- SIF addition experiment results (single-pixel) ---")
        println("  Added SIF coeff (truth): ", sif_coeff_add)
        println("  Retrieved WITH SIF:     ", sif_retrieved_with_sif)
        println("  Retrieved WITHOUT SIF:  ", sif_retrieved_without_sif)
        err_ws = sif_retrieved_with_sif .- sif_coeff_add
        println("  Difference (with SIF):  ", err_ws)
        println("  Relative error ev1:     ", length(sif_coeff_add) >= 1 && abs(sif_coeff_add[1]) > 0 ? err_ws[1] / sif_coeff_add[1] : "N/A")

        p_fit = plot(ctx.λ, y_obs; label="Observation", lw=2, color=:black)
        plot!(p_fit, ctx.λ, y_modified; label="Modified (with SIF)", lw=1.8, color=:blue)
        plot!(p_fit, ctx.λ, y_reconstructed_with_sif; label="Reconstructed (with SIF)", lw=1.8, color=:green)
        plot!(p_fit, ctx.λ, y_reconstructed_without_sif; label="Reconstructed (without SIF)", lw=1.5, color=:orange, ls=:dash)
        savefig(p_fit, joinpath(PSEUDO_DIR, "sif_addition_fit.png"))
        println("  Saved plot to sif_addition_fit.png")

        p_sif = plot(ctx.λ_hres, sif_added_hres; label="Added SIF (truth)", lw=2, color=:black)
        plot!(p_sif, ctx.λ_hres, sif_hres_with_sif; label="Retrieved (with SIF)", lw=1.8, color=:red)
        plot!(p_sif, ctx.λ_hres, sif_hres_without_sif; label="Retrieved (without SIF)", lw=1.5, color=:orange, ls=:dash)
        savefig(p_sif, joinpath(PSEUDO_DIR, "sif_addition_sif_signal.png"))
        println("  Saved plot to sif_addition_sif_signal.png")

        # Store both results to JLD2 for programmatic access
        out_jld2 = joinpath(PSEUDO_DIR, "sif_addition_single_pixel_results.jld2")
        JLD2.@save out_jld2 sif_coeff_add sif_retrieved_with_sif sif_retrieved_without_sif x_with_sif x_without_sif y_obs y_modified
        println("  Saved results to ", out_jld2)

        return (sif_added=sif_coeff_add, sif_retrieved_with_sif=sif_retrieved_with_sif, sif_retrieved_without_sif=sif_retrieved_without_sif, y_obs=y_obs, y_modified=y_modified, x_final_with_sif=x_with_sif, x_final_without_sif=x_without_sif)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_sif_addition()
end
