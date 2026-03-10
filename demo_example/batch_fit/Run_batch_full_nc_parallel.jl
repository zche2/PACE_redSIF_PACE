#!/usr/bin/env julia
# Parallel batch retrieval over the full .nc swath.
# Thread-based pixel parallelization: processes pixels within each scan in parallel.
# Uses per-thread retrieval cores (forward model + Jacobian) for thread safety.
# Same config and output format as Run_batch_full_nc.jl; writes to a different
# output file by default (_retrieval_full_parallel.nc) to avoid overwriting.
#
# Usage: julia -t N demo_example/batch_fit/Run_batch_full_nc_parallel.jl
#        (N = number of threads, e.g. 8)

using Base.Threads
using TOML
using NCDatasets
using LinearAlgebra
using SparseArrays
using Statistics
using Dates

const _BATCH_FIT_DIR = @__DIR__
const _DEMO_DIR = joinpath(_BATCH_FIT_DIR, "..")

# Include shared batch logic (Fit_toy_forward_model, helpers, _build_retrieval_core, etc.)
include(joinpath(_BATCH_FIT_DIR, "Run_batch_full_nc.jl"))

const STATUS_PIXEL_FILTER_SKIPPED = Int16(7)

# ----- Parallel-specific output path -----
function _make_output_path_parallel(pace_path::AbstractString, cfg::Dict)
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    out_dir_default = joinpath(_DEMO_DIR, "batch_output")
    out_dir_cfg = String(get(batch_cfg, "output_dir", out_dir_default))
    out_dir = isabspath(out_dir_cfg) ? out_dir_cfg : joinpath(_DEMO_DIR, out_dir_cfg)
    mkpath(out_dir)
    suffix = String(get(batch_cfg, "output_suffix_parallel", "_retrieval_full_parallel.nc"))
    stem = splitext(basename(pace_path))[1]
    return joinpath(out_dir, stem * suffix)
end

# ----- Process one pixel (thread-safe: uses thread-local core and buffers) -----
function _process_one_pixel!(
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
    spec_raw = slab_is_pix_band ? view(slab, i_pix_src, :) : view(slab, :, i_pix_src)
    ok = _copy_sorted_spectrum!(buf.y_sorted, spec_raw, perm)
    if !ok
        status_scan[i_pix_out] = Int16(3)
        return
    end
    wm_val = watermask[i_pix_src, j_scan_src]
    is_ocean = !ismissing(wm_val) && (Int(wm_val) in ocean_mask_values)
    ocean_scan[i_pix_out] = is_ocean ? UInt8(1) : UInt8(0)
    if ocean_filter_enabled && !is_ocean
        status_scan[i_pix_out] = Int16(6)
        return
    end
    mul!(buf.y_obs, W_interp, buf.y_sorted)
    is_dark = maximum(buf.y_obs) <= dark_max_radiance
    dark_scan[i_pix_out] = is_dark ? UInt8(1) : UInt8(0)
    if dark_filter_enabled && !is_dark
        status_scan[i_pix_out] = Int16(5)
        return
    end
    stats = _run_one_retrieval!(buf.x_tmp, core, buf.y_obs, max_outer_steps)
    state_scan[i_pix_out, :] .= Float32.(buf.x_tmp)
    conv_scan[i_pix_out] = stats.converged ? UInt8(1) : UInt8(0)
    status_scan[i_pix_out] = stats.status
    steps_scan[i_pix_out] = Int16(stats.n_steps)
    rmse_scan[i_pix_out] = Float32(stats.rmse)
    rchi2_scan[i_pix_out] = Float32(stats.reduced_chi2)
    obj_scan[i_pix_out] = Float32(stats.objective)
    sif_coeff = buf.x_tmp[core.layout.idx_sif]
    if length(sif_coeff) >= 1
        sif1_scan[i_pix_out] = Float32(sif_coeff[1])
    end
    sif_coeffs_scan[i_pix_out, :] .= Float32.(sif_coeff)
    sif_678_scan[i_pix_out] = Float32(dot(core.sif_basis_678, sif_coeff))
    return
end

# ----- Full-swath orbit run (parallel) -----
function run_orbit_full_nc_parallel(
    cores::Vector,
    pace_path::AbstractString,
    config_path::AbstractString,
)
    core = cores[1]
    cfg = core.cfg
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    use_threads = Bool(get(batch_cfg, "use_threads", true))
    n_threads = use_threads ? Threads.nthreads() : 1

    wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
    spectrum_var = String(get(pace_cfg, "spectrum_var", "radiance_red"))
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
    max_outer_steps = Int(get(batch_cfg, "max_outer_steps", core.lm.max_outer_default))
    dark_filter_enabled = Bool(get(batch_cfg, "dark_filter_enabled", false))
    dark_max_radiance = Float64(get(batch_cfg, "dark_max_radiance", 20.0))
    ocean_filter_enabled = Bool(get(batch_cfg, "ocean_filter_enabled", false))
    watermask_var = String(get(batch_cfg, "watermask_var", "watermask"))
    ocean_mask_values_raw = get(batch_cfg, "ocean_mask_values", Any[1])
    ocean_mask_values = Set{Int}(Int(v) for v in ocean_mask_values_raw)
    isempty(ocean_mask_values) && error("batch_fit.ocean_mask_values must contain at least one value")
    pixel_filter_vars = get(batch_cfg, "pixel_filter_vars", ["nflh"])
    pixel_filter_vars = isa(pixel_filter_vars, AbstractVector) ?
        String.(pixel_filter_vars) : String[String(pixel_filter_vars)]
    eligible = build_pixel_eligible_mask(ds, n_pix, n_scan, pixel_filter_vars)
    lat = _read_geo_2d(ds, "latitude", n_pix, n_scan)
    lon = _read_geo_2d(ds, "longitude", n_pix, n_scan)
    watermask = if haskey(ds, watermask_var)
        wm = ds[watermask_var]
        ndims(wm) == 2 || error("Expected 2D watermask variable '$watermask_var', got ndims=$(ndims(wm))")
        d = collect(String.(dimnames(wm)))
        raw = wm[:, :]
        if d == ["pixels", "scans"]
            raw
        elseif d == ["scans", "pixels"]
            permutedims(raw, (2, 1))
        else
            error("Unsupported watermask dims for '$watermask_var': $d")
        end
    else
        if ocean_filter_enabled
            error("Ocean filter enabled but variable '$watermask_var' is missing in $pace_path")
        end
        fill(missing, n_pix, n_scan)
    end

    output_path = _make_output_path_parallel(pace_path, cfg)
    n_sif_ev = length(core.layout.idx_sif)
    ds_out = _create_output_dataset(
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
    println("Running full-swath retrieval (parallel) for: ", pace_path)
    println("  threads: ", n_threads)
    println("  pixels: ", first(pixel_range), ":", last(pixel_range), " (", length(pixel_range), ")")
    println("  scans:  ", first(scan_range), ":", last(scan_range), " (", length(scan_range), ")")
    println("  pixel_filter_vars: ", pixel_filter_vars, " -> ", n_eligible, " eligible pixels")
    println("  n_state: ", core.layout.n_state, "  n_meas: ", length(core.ctx.λ))
    println("  dark filter: ", dark_filter_enabled, " (max radiance <= ", dark_max_radiance, ")")
    println("  ocean filter: ", ocean_filter_enabled, " (", watermask_var, " in ", collect(ocean_mask_values), ")")

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
        (slab_is_pix_band || slab_is_band_pix) ||
            error("Unexpected slab size $(size(slab)) for scan=$j_scan_src")

        state_scan = fill(Float32(NaN), length(pixel_range), core.layout.n_state)
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

        if use_threads && n_threads > 1
            Threads.@threads for i_pix_out in 1:length(pixel_range)
                tid = Threads.threadid()
                i_pix_src = pixel_range[i_pix_out]
                c = cores[tid]
                b = buffers[tid]
                _process_one_pixel!(
                    i_pix_out,
                    i_pix_src,
                    j_scan_src,
                    c,
                    b,
                    slab,
                    slab_is_pix_band,
                    perm,
                    W_interp,
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
                _process_one_pixel!(
                    i_pix_out,
                    i_pix_src,
                    j_scan_src,
                    cores[1],
                    buffers[1],
                    slab,
                    slab_is_pix_band,
                    perm,
                    W_interp,
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
        end

        ds_out["x_hat"][:, j_scan_out, :] = state_scan
        ds_out["converged"][:, j_scan_out] = conv_scan
        ds_out["status_code"][:, j_scan_out] = status_scan
        ds_out["n_steps"][:, j_scan_out] = steps_scan
        ds_out["rmse"][:, j_scan_out] = rmse_scan
        ds_out["reduced_chi2"][:, j_scan_out] = rchi2_scan
        ds_out["objective"][:, j_scan_out] = obj_scan
        ds_out["sif_ev1"][:, j_scan_out] = sif1_scan
        ds_out["sif_coeffs"][:, j_scan_out, :] = sif_coeffs_scan
        ds_out["sif_radiance_678nm"][:, j_scan_out] = sif_678_scan
        ds_out["is_dark"][:, j_scan_out] = dark_scan
        ds_out["is_ocean"][:, j_scan_out] = ocean_scan
        n_conv = count(==(UInt8(1)), conv_scan)
        n_dark = count(==(UInt8(1)), dark_scan)
        n_ocean = count(==(UInt8(1)), ocean_scan)
        n_skip = count(==(STATUS_PIXEL_FILTER_SKIPPED), status_scan)
        println(
            "  scan ", j_scan_src,
            " -> filter_skip ", n_skip,
            " | ocean ", n_ocean, "/", length(pixel_range),
            " | dark ", n_dark, "/", length(pixel_range),
            " | converged ", n_conv, "/", length(pixel_range),
        )
    end
    close(ds_out)
    close(ds)
    println("Saved full-swath retrieval to: ", output_path)
    return output_path
end

function main_full_nc_parallel()
    config_path = get(
        ENV,
        "PACE_MWE_CONFIG",
        normpath(joinpath(_DEMO_DIR, "Simple_PACE_xSecFit_MWE_zcheVer.toml")),
    )
    cfg = TOML.parsefile(config_path)
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    use_threads = Bool(get(batch_cfg, "use_threads", true))
    n_threads = use_threads ? Threads.nthreads() : 1

    cores = [_build_retrieval_core(config_path) for _ in 1:n_threads]
    files = _resolve_orbit_files(cores[1].cfg)
    println("Batch full-NC retrieval (parallel, pixel filter: batch_fit.pixel_filter_vars)")
    println("  config: ", config_path)
    println("  threads: ", n_threads)
    println("  n_orbits: ", length(files))
    for f in files
        run_orbit_full_nc_parallel(cores, f, config_path)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_full_nc_parallel()
end
