# =========================================================================================================
# Parallel batch retrieval with SIF addition.
# Thread-based pixel parallelization: add SIF+noise per pixel, then run retrieval.
# Used by granule_retrieval_adding_sif.jl (like Run_batch_full_nc_parallel for granule_retrieval.jl).
# Output includes spectrally resolved rmse, sif_added_678nm.
# =========================================================================================================

using Base.Threads
using TOML
using NCDatasets
using JLD2
using Interpolations
using LinearAlgebra
using SparseArrays
using Statistics
using Dates

const _BATCH_FIT_DIR = @__DIR__
const _DEMO_DIR = joinpath(_BATCH_FIT_DIR, "..")

include(joinpath(_BATCH_FIT_DIR, "Run_batch_full_nc.jl"))

# ----- SIF helpers (from granule_retrieval_adding_sif.jl) -----
"""
    load_sif_shapes_lres(ctx, sif_path; magnitude=0.5)

Load SIF_shapes from jld2, interpolate each to ctx.λ (lres), scale to given magnitude.
Returns Vector{Vector{Float64}} of length n_shapes.
"""
function load_sif_shapes_lres(ctx, sif_path::AbstractString; magnitude::Float64=0.5)
    sif = JLD2.load(sif_path)
    haskey(sif, "SIF_shapes") && haskey(sif, "SIF_wavelen") || error("SIF_shapes and SIF_wavelen required in $sif_path")
    shapes_raw = Float64.(sif["SIF_shapes"])
    λ_ref = collect(Float64.(sif["SIF_wavelen"]))
    n_bands = size(shapes_raw, 1)
    n_shapes = size(shapes_raw, 2)
    n_bands == length(λ_ref) || error("SIF_shapes first dim must match SIF_wavelen")

    λ_target = collect(Float64.(ctx.λ))
    idx_678 = argmin(abs.(λ_target .- 678.2))
    sif_shapes_lres = Vector{Vector{Float64}}(undef, n_shapes)

    for k in 1:n_shapes
        itp = LinearInterpolation(λ_ref, shapes_raw[:, k], extrapolation_bc=0.0)
        s_lres = Float64.(itp.(λ_target))
        s_val = abs(s_lres[idx_678] > 0 ? s_lres[idx_678] : maximum(abs.(s_lres)))
        s_val > 0 || (s_val = 1.0)
        scale = magnitude / s_val
        sif_shapes_lres[k] = s_lres .* scale
    end
    return sif_shapes_lres
end

"""
    add_sif_and_noise!(y_obs, sif_lres, band_snr_coeffs, add_noise)

Add sif_lres to y_obs. If add_noise, add Gaussian noise with sigma = sqrt(c1 + c2 * sif_lres).
Modifies y_obs in place.
"""
function add_sif_and_noise!(
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

# ----- Output dataset (AddingSIF-specific) -----
function _create_adding_sif_output_dataset(
    output_path::AbstractString,
    n_pix::Int,
    n_scan::Int,
    n_sif_ev::Int,
    n_red_wl::Int,
    state_names::Vector{String},
    red_wavelength::Vector{Float64},
    pace_path::AbstractString,
    config_path::AbstractString,
    pixel_range::UnitRange{Int},
    scan_range::UnitRange{Int};
    store_rtoa::Bool=true,
)
    println("Creating output at: ", output_path)
    ds = Dataset(output_path, "c")
    defDim(ds, "pixels", n_pix)
    defDim(ds, "scans", n_scan)
    defDim(ds, "state", length(state_names))
    defDim(ds, "sif_nev", n_sif_ev)
    defDim(ds, "red_wavelength", n_red_wl)
    ds.attrib["title"] = "PACE AddingSIF retrieval output (spectrally resolved rmse)"
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
    v_state = defVar(ds, "x_hat", Float32, ("pixels", "scans", "state"); comp...)
    v_conv = defVar(ds, "converged", UInt8, ("pixels", "scans"); comp...)
    v_sif1 = defVar(ds, "sif_ev1", Float32, ("pixels", "scans"); comp...)
    v_sif_coeffs = defVar(ds, "sif_coeffs", Float32, ("pixels", "scans", "sif_nev"); comp...)
    v_sif_678 = defVar(ds, "sif_radiance_678nm", Float32, ("pixels", "scans"); comp...)
    v_sif_added = defVar(ds, "sif_added_678nm", Float32, ("pixels", "scans"); comp...)
    v_rmse = defVar(ds, "rmse", Float32, ("pixels", "scans", "red_wavelength"); comp...)
    if store_rtoa
        v_rtoa = defVar(ds, "Rtoa_red", Float32, ("pixels", "scans", "red_wavelength"); comp...)
        v_rtoa.attrib["long_name"] = "Observed TOA radiance (after SIF+noise) per band"
    end
    v_red_wl = defVar(ds, "red_wavelength", Float32, ("red_wavelength",); comp...)
    v_state.attrib["long_name"] = "Retrieved state vector"
    v_conv.attrib["long_name"] = "1 if convergence criterion reached, 0 otherwise"
    v_sif1.attrib["long_name"] = "Retrieved first SIF eigenvector coefficient"
    v_sif_coeffs.attrib["long_name"] = "All retrieved SIF eigenvector coefficients"
    v_sif_678.attrib["long_name"] = "Reconstructed SIF radiance at 678.2 nm"
    v_sif_added.attrib["long_name"] = "Truth: added SIF at 678 nm"
    v_rmse.attrib["long_name"] = "Spectrally resolved residual (y_obs - y_fit) per band"
    v_red_wl[:] = Float32.(red_wavelength)
    return ds
end

function _make_adding_sif_output_path(interim_path::AbstractString, cfg::Dict)
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())
    add_cfg = get(cfg, "adding_sif", Dict{String, Any}())
    out_dir = String(get(add_cfg, "output_dir", get(batch_cfg, "output_dir", "AddingSIF_output")))
    out_dir = isabspath(out_dir) ? out_dir : joinpath(_DEMO_DIR, out_dir)
    mkpath(out_dir)
    suffix = String(get(batch_cfg, "output_suffix_parallel", "_adding_sif.nc"))
    stem = splitext(basename(interim_path))[1]
    m = match(r"^interim_(.+)$", stem)
    granule_id = m !== nothing ? m.captures[1] : stem
    return joinpath(out_dir, "PACE_OCI.$(granule_id).L1B.V3" * suffix)
end

# ----- Process one pixel (thread-safe: SIF+noise then retrieval) -----
function _process_one_pixel_adding_sif!(
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
    sif_shapes_lres,
    n_shapes::Int,
    n_pix::Int,
    band_snr_coeffs,
    add_snr_noise::Bool,
    idx_678::Int,
    max_outer_steps::Int,
    store_rtoa::Bool,
    state_scan,
    conv_scan,
    sif1_scan,
    sif_coeffs_scan,
    sif_678_scan,
    sif_added_scan,
    rmse_scan,
    rtoa_scan,
)
    if !eligible[i_pix_src, j_scan_src]
        return
    end
    spec_raw = slab_is_pix_band ? view(slab, i_pix_src, :) : view(slab, :, i_pix_src)
    ok = _copy_sorted_spectrum!(buf.y_sorted, spec_raw, perm)
    if !ok
        return
    end
    mul!(buf.y_obs, W_interp, buf.y_sorted)

    shape_idx = ((j_scan_src - 1) * n_pix + i_pix_src - 1) % n_shapes + 1
    sif_lres = sif_shapes_lres[shape_idx]
    sif_added_scan[i_pix_out] = Float32(sif_lres[idx_678])
    add_sif_and_noise!(buf.y_obs, sif_lres, band_snr_coeffs, add_snr_noise)
    if store_rtoa && rtoa_scan !== nothing
        rtoa_scan[i_pix_out, :] .= Float32.(buf.y_obs)
    end

    stats = _run_one_retrieval!(buf.x_tmp, core, buf.y_obs, max_outer_steps)
    y_fit = core.fm(buf.x_tmp)
    residual = buf.y_obs .- y_fit

    state_scan[i_pix_out, :] .= Float32.(buf.x_tmp)
    conv_scan[i_pix_out] = stats.converged ? UInt8(1) : UInt8(0)
    sif_coeff = buf.x_tmp[core.layout.idx_sif]
    if length(sif_coeff) >= 1
        sif1_scan[i_pix_out] = Float32(sif_coeff[1])
    end
    sif_coeffs_scan[i_pix_out, :] .= Float32.(sif_coeff)
    sif_678_scan[i_pix_out] = Float32(dot(core.sif_basis_678, sif_coeff))
    rmse_scan[i_pix_out, :] .= Float32.(residual)
    return
end

# ----- Full orbit run (parallel, with SIF addition) -----
function run_orbit_full_nc_adding_sif_parallel(
    cores::Vector,
    interim_path::AbstractString,
    config_path::AbstractString,
)
    core = cores[1]
    cfg = core.cfg
    add_cfg = get(cfg, "adding_sif", Dict{String, Any}())
    data_cfg = get(cfg, "data", Dict{String, Any}())
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    batch_cfg = get(cfg, "batch_fit", Dict{String, Any}())

    use_threads = Bool(get(batch_cfg, "use_threads", true))
    n_threads = use_threads ? Threads.nthreads() : 1

    sif_magnitude = Float64(get(add_cfg, "sif_magnitude", 0.5))
    add_snr_noise = Bool(get(add_cfg, "add_snr_noise", true))
    store_rtoa = Bool(get(add_cfg, "store_rtoa", true))

    sif_jld2_path = String(get(add_cfg, "sif_jld2_path", ""))
    if isempty(sif_jld2_path)
        sif_jld2_path = core.ctx.paths.sif_path
    else
        base_dir = String(get(data_cfg, "base_dir", ""))
        sif_jld2_path = isabspath(sif_jld2_path) ? sif_jld2_path : joinpath(base_dir, sif_jld2_path)
    end
    isfile(sif_jld2_path) || error("SIF shapes file not found: $sif_jld2_path")

    sif_shapes_lres = load_sif_shapes_lres(core.ctx, sif_jld2_path; magnitude=sif_magnitude)
    n_shapes = length(sif_shapes_lres)
    band_snr_coeffs = core.ctx.band_snr_coeffs
    λ_target = collect(Float64.(core.ctx.λ))
    idx_678 = argmin(abs.(λ_target .- 678.2))

    wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
    spectrum_var = String(get(pace_cfg, "spectrum_var", "Rtoa_red"))

    ds = Dataset(interim_path)
    haskey(ds, spectrum_var) || error("Missing spectrum variable '$spectrum_var'")
    v_spec = ds[spectrum_var]
    axes_info = _find_axis_indices(v_spec, wavelength_var)
    λ_src = collect(Float64.(ds[wavelength_var][:]))
    W_interp, perm = _make_linear_resampler(λ_src, core.ctx.λ)

    n_pix = axes_info.n_pix
    n_scan = axes_info.n_scan
    pixel_range = 1:n_pix
    scan_range = 1:n_scan

    p_start = Int(get(add_cfg, "pixel_start", 1))
    p_end_cfg = Int(get(add_cfg, "pixel_end", 0))
    s_start = Int(get(add_cfg, "scan_start", 1))
    s_end_cfg = Int(get(add_cfg, "scan_end", 0))
    p_end = p_end_cfg > 0 ? min(p_end_cfg, n_pix) : n_pix
    s_end = s_end_cfg > 0 ? min(s_end_cfg, n_scan) : n_scan
    p_start = clamp(p_start, 1, p_end)
    s_start = clamp(s_start, 1, s_end)
    pixel_range = p_start:p_end
    scan_range = s_start:s_end

    pixel_filter_vars = get(batch_cfg, "pixel_filter_vars", ["nflh"])
    pixel_filter_vars = isa(pixel_filter_vars, AbstractVector) ? String.(pixel_filter_vars) : String[String(pixel_filter_vars)]
    eligible = build_pixel_eligible_mask(ds, n_pix, n_scan, pixel_filter_vars)

    lat = _read_geo_2d(ds, "latitude", n_pix, n_scan)
    lon = _read_geo_2d(ds, "longitude", n_pix, n_scan)
    max_outer_steps = Int(get(batch_cfg, "max_outer_steps", core.lm.max_outer_default))
    n_sif_ev = length(core.layout.idx_sif)
    n_red_wl = length(core.ctx.λ)

    output_path = _make_adding_sif_output_path(interim_path, cfg)
    ds_out = _create_adding_sif_output_dataset(
        output_path,
        length(pixel_range),
        length(scan_range),
        n_sif_ev,
        n_red_wl,
        core.state_names,
        λ_target,
        interim_path,
        config_path,
        pixel_range,
        scan_range;
        store_rtoa=store_rtoa,
    )
    ds_out["latitude"][:, :] = lat[pixel_range, scan_range]
    ds_out["longitude"][:, :] = lon[pixel_range, scan_range]

    buffers = [(
        y_sorted = zeros(Float64, length(perm)),
        y_obs = zeros(Float64, length(core.ctx.λ)),
        x_tmp = zeros(Float64, core.layout.n_state),
    ) for _ in 1:n_threads]

    println("AddingSIF retrieval (parallel): ", interim_path)
    println("  threads: ", n_threads)
    println("  pixels: ", first(pixel_range), ":", last(pixel_range), " (", length(pixel_range), ")  scans: ", first(scan_range), ":", last(scan_range), " (", length(scan_range), ")")
    println("  sif_shapes: ", n_shapes, "  magnitude: ", sif_magnitude, "  add_snr_noise: ", add_snr_noise, "  store_rtoa: ", store_rtoa)

    for (j_scan_out, j_scan_src) in enumerate(scan_range)
        inds = Any[Colon() for _ in 1:3]
        inds[axes_info.i_scan] = j_scan_src
        slab = v_spec[inds...]
        slab_is_pix_band = size(slab) == (axes_info.n_pix, axes_info.n_band)
        slab_is_band_pix = size(slab) == (axes_info.n_band, axes_info.n_pix)
        (slab_is_pix_band || slab_is_band_pix) || error("Unexpected slab size for scan=$j_scan_src")

        state_scan = fill(Float32(NaN), length(pixel_range), core.layout.n_state)
        conv_scan = fill(UInt8(0), length(pixel_range))
        sif1_scan = fill(Float32(NaN), length(pixel_range))
        sif_coeffs_scan = fill(Float32(NaN), length(pixel_range), n_sif_ev)
        sif_678_scan = fill(Float32(NaN), length(pixel_range))
        sif_added_scan = fill(Float32(NaN), length(pixel_range))
        rmse_scan = fill(Float32(NaN), length(pixel_range), n_red_wl)
        rtoa_scan = store_rtoa ? fill(Float32(NaN), length(pixel_range), n_red_wl) : nothing

        if use_threads && n_threads > 1
            Threads.@threads for i_pix_out in 1:length(pixel_range)
                tid = Threads.threadid()
                i_pix_src = pixel_range[i_pix_out]
                c = cores[tid]
                b = buffers[tid]
                _process_one_pixel_adding_sif!(
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
                    sif_shapes_lres,
                    n_shapes,
                    n_pix,
                    band_snr_coeffs,
                    add_snr_noise,
                    idx_678,
                    max_outer_steps,
                    store_rtoa,
                    state_scan,
                    conv_scan,
                    sif1_scan,
                    sif_coeffs_scan,
                    sif_678_scan,
                    sif_added_scan,
                    rmse_scan,
                    rtoa_scan,
                )
            end
        else
            for (i_pix_out, i_pix_src) in enumerate(pixel_range)
                _process_one_pixel_adding_sif!(
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
                    sif_shapes_lres,
                    n_shapes,
                    n_pix,
                    band_snr_coeffs,
                    add_snr_noise,
                    idx_678,
                    max_outer_steps,
                    store_rtoa,
                    state_scan,
                    conv_scan,
                    sif1_scan,
                    sif_coeffs_scan,
                    sif_678_scan,
                    sif_added_scan,
                    rmse_scan,
                    rtoa_scan,
                )
            end
        end

        ds_out["x_hat"][:, j_scan_out, :] = state_scan
        ds_out["converged"][:, j_scan_out] = conv_scan
        ds_out["sif_ev1"][:, j_scan_out] = sif1_scan
        ds_out["sif_coeffs"][:, j_scan_out, :] = sif_coeffs_scan
        ds_out["sif_radiance_678nm"][:, j_scan_out] = sif_678_scan
        ds_out["sif_added_678nm"][:, j_scan_out] = sif_added_scan
        ds_out["rmse"][:, j_scan_out, :] = rmse_scan
        if store_rtoa
            ds_out["Rtoa_red"][:, j_scan_out, :] = rtoa_scan
        end
        n_conv = count(==(UInt8(1)), conv_scan)
        println("  scan ", j_scan_src, " -> converged ", n_conv, "/", length(pixel_range))
    end

    close(ds_out)
    close(ds)
    println("Saved AddingSIF output to: ", output_path)
    return output_path
end
