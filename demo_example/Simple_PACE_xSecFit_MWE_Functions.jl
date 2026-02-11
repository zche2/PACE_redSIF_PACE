module SimplePACEXSecFitMWEFunctions

using TOML
using JLD2
using Interpolations
using NCDatasets
using DelimitedFiles
using PACE_SIF

export read_mwe_config,
       resolve_paths,
       load_sif_basis,
       load_solar_spectrum_on_grid,
       load_pace_spectrum_on_grid,
       infer_lut_spectral_axis,
       build_highres_wavelength_grid,
       make_kernel_spectral_axis,
       subset_kernel,
       build_kernel_from_rsr_nc,
       kernel_quality_metrics,
       compare_kernels,
       prepare_mwe_inputs

"""
    cfg_get(cfg, section, key, default)

Safely fetch a config value from nested TOML dictionary structure.
Returns `default` if section/key is missing or malformed.
"""
function cfg_get(cfg::Dict, section::String, key::String, default)
    haskey(cfg, section) || return default
    sec = cfg[section]
    sec isa Dict || return default
    return get(sec, key, default)
end

"""
    must_exist(path)

Validate that a file exists and return the original path.
Throws an error with a readable message if missing.
"""
function must_exist(path::AbstractString)
    isfile(path) || error("Missing required file: $path")
    return path
end

"""
    read_mwe_config(config_path)

Read and parse the TOML config file for the MWE.
"""
function read_mwe_config(config_path::AbstractString)
    return TOML.parsefile(config_path)
end

"""
    resolve_paths(cfg)

Resolve all path references from config and environment overrides.
`PACE_DATA_DIR` and `PACE_XSEC_RUN` take precedence over TOML values.
"""
function resolve_paths(cfg::Dict)
    base_dir = get(ENV, "PACE_DATA_DIR", cfg_get(cfg, "data", "base_dir", "/Users/cfranken/data"))
    xsec_run = get(ENV, "PACE_XSEC_RUN", cfg_get(cfg, "data", "xsec_run", "wavelength-run"))

    xsec_dir_cfg = cfg_get(cfg, "data", "xsec_dir", "")
    xsec_dir = if !isempty(xsec_dir_cfg)
        xsec_dir_cfg
    else
        d1 = joinpath(base_dir, "interp_xSection", xsec_run)
        d2 = joinpath(base_dir, xsec_run)
        isdir(d1) ? d1 : (isdir(d2) ? d2 : d1)
    end

    kernel_file = cfg_get(cfg, "data", "kernel_file", "KernelInstrument.jld2")
    kernel_path = isabspath(kernel_file) ? kernel_file : joinpath(base_dir, kernel_file)

    sif_file = cfg_get(cfg, "data", "sif_file", "SIF_singular_vector.jld2")
    sif_path = if isabspath(sif_file)
        sif_file
    else
        direct = joinpath(base_dir, sif_file)
        alt = joinpath(base_dir, "reference_spectra", sif_file)
        isfile(direct) ? direct : alt
    end

    o2_file = cfg_get(cfg, "data", "o2_file", "$(xsec_run)_O2.jld2")
    h2o_file = cfg_get(cfg, "data", "h2o_file", "$(xsec_run)_H2O.jld2")
    o2_path = isabspath(o2_file) ? o2_file : joinpath(xsec_dir, o2_file)
    h2o_path = isabspath(h2o_file) ? h2o_file : joinpath(xsec_dir, h2o_file)

    pace_rsr_file = cfg_get(cfg, "data", "pace_rsr_file", "PACE_OCI_RSRs.nc")
    pace_rsr_path = isabspath(pace_rsr_file) ? pace_rsr_file : joinpath(base_dir, pace_rsr_file)

    regen_kernel_file = cfg_get(cfg, "data", "regenerated_kernel_file", "KernelInstrument_regenerated.jld2")
    regen_kernel_path = isabspath(regen_kernel_file) ? regen_kernel_file : joinpath(base_dir, regen_kernel_file)

    return (
        base_dir = base_dir,
        xsec_run = xsec_run,
        xsec_dir = xsec_dir,
        kernel_path = kernel_path,
        sif_path = sif_path,
        o2_path = o2_path,
        h2o_path = h2o_path,
        pace_rsr_path = pace_rsr_path,
        regen_kernel_path = regen_kernel_path,
    )
end

"""
    load_sif_basis(sif_path, λ; nEV=2, normalize=true)

Interpolate `SIF_U[:, 1:nEV]` onto wavelength grid `λ`.
Uses cubic spline interpolation when `SIF_wavelen` is evenly spaced,
and falls back to linear interpolation otherwise.
Returns matrix with size `(length(λ), nEV)`.
"""
function load_sif_basis(
    sif_path::AbstractString,
    λ::AbstractVector{<:Real};
    nEV::Int=2,
    normalize::Bool=true,
)
    nEV < 1 && error("nEV must be >= 1")
    sif = JLD2.load(must_exist(sif_path))

    sif_u = convert.(Float64, sif["SIF_U"])
    λ_ref = collect(Float64.(sif["SIF_wavelen"]))

    size(sif_u, 1) == length(λ_ref) ||
        error("SIF_U first dimension must match SIF_wavelen length")

    n_available = size(sif_u, 2)
    n_use = min(nEV, n_available)
    n_use < nEV && @warn "Requested nEV=$nEV but only $n_available available; using $n_use"

    # CubicSplineInterpolation requires regularly spaced knot ranges.
    dλ = diff(λ_ref)
    step = dλ[1]
    tol = max(1e-10, abs(step) * 1e-8)
    use_cubic = all(abs.(dλ .- step) .<= tol)
    if use_cubic
        λ_knots = range(λ_ref[1], step=step, length=length(λ_ref))
    else
        @warn "SIF_wavelen is not evenly spaced; using linear interpolation."
    end

    # Build basis matrix EV-by-EV for clarity.
    basis = zeros(Float64, length(λ), n_use)
    for iev in 1:n_use
        if use_cubic
            itp_ev = CubicSplineInterpolation(λ_knots, sif_u[:, iev]; extrapolation_bc=Line())
            vals = itp_ev.(λ)
            vals[(λ .< λ_ref[1]) .| (λ .> λ_ref[end])] .= 0.0
            basis[:, iev] .= vals
        else
            itp_ev = LinearInterpolation(λ_ref, sif_u[:, iev], extrapolation_bc=0.0)
            basis[:, iev] .= itp_ev.(λ)
        end
    end

    if normalize
        s = maximum(abs.(basis[:, 1]))
        s > 0 && (basis ./= s)
    end
    return basis
end

"""
    load_solar_spectrum_on_grid(solar_path, λ_target; header_lines=3)

Read the ASCII solar file with columns:
1) wavenumber [cm^-1]
2) solar/transmittance value

Then interpolate onto `λ_target` (nm) and return `(E_target, info)`, where:
- `E_target`: solar values on `λ_target`
- `info`: metadata tuple with source ranges and sample count

The interpolation is strict: it throws if `λ_target` is outside source range.
"""
function load_solar_spectrum_on_grid(
    solar_path::AbstractString,
    λ_target::AbstractVector{<:Real};
    header_lines::Int=3,
)
    tbl = readdlm(must_exist(solar_path), Float64, skipstart=header_lines)
    size(tbl, 2) >= 2 || error("Solar file must contain at least 2 numeric columns.")

    ν_src = tbl[:, 1]
    E_src = tbl[:, 2]
    λ_src = 1e7 ./ ν_src

    p = sortperm(λ_src)
    λ_src_sorted = λ_src[p]
    E_src_sorted = E_src[p]

    λt = collect(Float64.(λ_target))
    λt_min, λt_max = extrema(λt)
    λs_min, λs_max = extrema(λ_src_sorted)
    (λt_min >= λs_min && λt_max <= λs_max) ||
        error("Target λ range [$λt_min, $λt_max] is outside solar λ range [$λs_min, $λs_max].")

    itp = LinearInterpolation(λ_src_sorted, E_src_sorted, extrapolation_bc=Throw())
    E_target = itp.(λt)

    info = (
        n_source = length(λ_src_sorted),
        λ_source_range_nm = (λs_min, λs_max),
        λ_target_range_nm = (λt_min, λt_max),
    )
    return E_target, info
end

"""
    load_pace_spectrum_on_grid(
        pace_path,
        λ_target;
        pixel_idx=1,
        scan_idx=1,
        wavelength_var="red_wavelength",
        spectrum_var="radiance_red",
    )

Read one PACE spectrum from NetCDF and interpolate it onto `λ_target` (nm).

Expected common layout is:
- `wavelength_var`: 1D wavelength vector (e.g. `red_wavelength`)
- `spectrum_var`: 3D spectral cube (e.g. `radiance_red[pixels, scans, red_bands]`)

Returns `(y_target, info)` where `y_target` is the interpolated spectrum and
`info` contains selected pixel/scan, source dimensions, and wavelength ranges.
"""
function load_pace_spectrum_on_grid(
    pace_path::AbstractString,
    λ_target::AbstractVector{<:Real};
    pixel_idx::Int=1,
    scan_idx::Int=1,
    wavelength_var::AbstractString="red_wavelength",
    spectrum_var::AbstractString="radiance_red",
)
    ds = Dataset(must_exist(pace_path))
    haskey(ds, wavelength_var) || error("PACE file missing variable '$wavelength_var'")
    haskey(ds, spectrum_var) || error("PACE file missing variable '$spectrum_var'")

    λ_src = collect(Float64.(ds[wavelength_var][:]))
    v = ds[spectrum_var]
    dnames = collect(String.(dimnames(v)))

    y_src = if ndims(v) == 1
        collect(Float64.(v[:]))
    elseif ndims(v) == 3
        i_pixel = findfirst(==("pixels"), dnames)
        i_scan = findfirst(==("scans"), dnames)
        i_band = findfirst(in(("red_bands", "bands", "wavelengths")), dnames)
        if isnothing(i_pixel) || isnothing(i_scan) || isnothing(i_band)
            # Fall back to the expected [pixel, scan, band] order.
            i_pixel, i_scan, i_band = 1, 2, 3
        end

        n_pixel = size(v, i_pixel)
        n_scan = size(v, i_scan)
        1 <= pixel_idx <= n_pixel || error("pixel_idx=$pixel_idx outside 1:$n_pixel")
        1 <= scan_idx <= n_scan || error("scan_idx=$scan_idx outside 1:$n_scan")

        inds = Any[Colon() for _ in 1:ndims(v)]
        inds[i_pixel] = pixel_idx
        inds[i_scan] = scan_idx
        inds[i_band] = Colon()
        collect(Float64.(v[inds...]))
    else
        close(ds)
        error("Unsupported '$spectrum_var' dimensions: ndims=$(ndims(v)); expected 1D or 3D")
    end
    close(ds)

    length(y_src) == length(λ_src) ||
        error("Length mismatch: $spectrum_var has $(length(y_src)) spectral samples, $wavelength_var has $(length(λ_src))")

    p = sortperm(λ_src)
    λ_src_sorted = λ_src[p]
    y_src_sorted = y_src[p]

    λt = collect(Float64.(λ_target))
    λt_min, λt_max = extrema(λt)
    λs_min, λs_max = extrema(λ_src_sorted)
    (λt_min >= λs_min && λt_max <= λs_max) ||
        error("Target λ range [$λt_min, $λt_max] is outside PACE λ range [$λs_min, $λs_max].")

    itp = LinearInterpolation(λ_src_sorted, y_src_sorted, extrapolation_bc=Throw())
    y_target = itp.(λt)

    info = (
        spectrum_var = String(spectrum_var),
        wavelength_var = String(wavelength_var),
        spectrum_dimnames = Tuple(dnames),
        pixel_idx = pixel_idx,
        scan_idx = scan_idx,
        λ_source_range_nm = (λs_min, λs_max),
        λ_target_range_nm = (λt_min, λt_max),
    )
    return y_target, info
end

"""
    infer_lut_spectral_axis(o2_sitp)

Infer whether LUT first axis is wavelength or wavenumber.
Returns `(spectral_axis_sorted, λ_sorted, axis_unit)`, with `λ_sorted` increasing.
"""
function infer_lut_spectral_axis(o2_sitp)
    spec_raw = collect(o2_sitp.ranges[1])
    is_wavenumber_axis = maximum(spec_raw) > 3000.0
    λ_raw = is_wavenumber_axis ? PACE_SIF.ν_to_λ.(spec_raw) : spec_raw
    axis_unit = is_wavenumber_axis ? "wavenumber_cm^-1" : "wavelength_nm"

    if λ_raw[1] <= λ_raw[end]
        spectral_axis = spec_raw
        λ_hres = λ_raw
    else
        spectral_axis = reverse(spec_raw)
        λ_hres = reverse(λ_raw)
    end

    return spectral_axis, λ_hres, axis_unit
end

"""
    build_highres_wavelength_grid(λ_lut_sorted; λ_min, λ_max, delta_lambda_nm, strict=true)

Construct the high-resolution wavelength grid used by the forward model.
If `strict=true`, requested bounds must lie within LUT wavelength support.
"""
function build_highres_wavelength_grid(
    λ_lut_sorted::AbstractVector{<:Real};
    λ_min::Float64,
    λ_max::Float64,
    delta_lambda_nm::Float64,
    strict::Bool=true,
)
    delta_lambda_nm > 0 || error("delta_lambda_nm must be > 0")
    λ_max > λ_min || error("lambda_max_nm must be > lambda_min_nm")

    lo_lut, hi_lut = extrema(collect(Float64.(λ_lut_sorted)))
    lo_req, hi_req = λ_min, λ_max

    if strict
        (lo_req >= lo_lut && hi_req <= hi_lut) ||
            error("Requested λ range [$lo_req, $hi_req] nm is outside LUT support [$lo_lut, $hi_lut] nm")
    else
        lo_req = clamp(lo_req, lo_lut, hi_lut)
        hi_req = clamp(hi_req, lo_lut, hi_lut)
    end

    n = floor(Int, (hi_req - lo_req) / delta_lambda_nm + 1e-12) + 1
    λ_hres = lo_req .+ delta_lambda_nm .* collect(0:(n - 1))
    if λ_hres[end] < hi_req - 1e-10
        push!(λ_hres, hi_req)
    end
    return λ_hres
end

"""
    make_kernel_spectral_axis(o2_sitp, λ_hres)

Build a spectral axis aligned to `λ_hres` for LUT interpolation.
Returns `(spectral_axis, axis_unit)`.
"""
function make_kernel_spectral_axis(
    o2_sitp,
    λ_hres::AbstractVector{<:Real},
)
    spec_grid_raw = collect(o2_sitp.ranges[1])
    is_wavenumber_axis = maximum(spec_grid_raw) > 3000.0
    axis_unit = is_wavenumber_axis ? "wavenumber_cm^-1" : "wavelength_nm"

    spec_eval = is_wavenumber_axis ? PACE_SIF.λ_to_ν.(λ_hres) : collect(λ_hres)
    lo, hi = extrema(spec_grid_raw)
    # Allow tiny numerical overshoots, but fail for real range mismatch.
    eps = 1e-8
    if minimum(spec_eval) < lo - eps || maximum(spec_eval) > hi + eps
        error("High-res grid maps outside LUT spectral bounds: [$lo, $hi]")
    end
    spec_eval = clamp.(spec_eval, lo, hi)

    if (spec_grid_raw[1] <= spec_grid_raw[end]) != (spec_eval[1] <= spec_eval[end])
        spec_eval = reverse(spec_eval)
    end
    return spec_eval, axis_unit
end

"""
    subset_kernel(kernel, λmin, λmax; spectral_axis=kernel.ν_grid)
"""
function subset_kernel(
    kernel,
    λmin::Float64,
    λmax::Float64;
    spectral_axis=kernel.ν_grid,
)
    idx = findall(λmin .< kernel.band .< λmax)
    isempty(idx) && error("No kernel bands in [$λmin, $λmax] nm")
    length(spectral_axis) == length(kernel.wvlen_out) ||
        error("spectral_axis length must match kernel.wvlen_out")

    band = kernel.band[idx]
    rsr = kernel.RSR[:, idx]
    return PACE_SIF.KernelInstrument(
        band,
        kernel.wvlen,
        rsr,
        kernel.wvlen_out,
        spectral_axis,
    )
end

"""
    build_kernel_from_rsr_nc(pace_rsr_path, λ_hres, spectral_axis; λ_min, λ_max, clip_negative)

Generate `KernelInstrument` directly from PACE RSR NetCDF and chosen high-res grid.
"""
function build_kernel_from_rsr_nc(
    pace_rsr_path::AbstractString,
    λ_hres::AbstractVector{<:Real},
    spectral_axis::AbstractVector{<:Real};
    λ_min::Float64,
    λ_max::Float64,
    clip_negative::Bool=true,
)
    ds = Dataset(must_exist(pace_rsr_path))
    wavlen = collect(Float64.(ds["wavelength"][:]))
    band = collect(Float64.(ds["bands"][:]))
    rsr_all = collect(Float64.(ds["RSR"][:, :]))
    close(ds)

    λ_lo = max(λ_min, minimum(λ_hres))
    λ_hi = min(λ_max, maximum(λ_hres))

    idx_w = findall(minimum(λ_hres) .< wavlen .< maximum(λ_hres))
    idx_b = findall(λ_lo .< band .< λ_hi)
    isempty(idx_w) && error("No RSR wavelength samples overlap λ_hres range")
    isempty(idx_b) && error("No RSR bands overlap requested range [$λ_lo, $λ_hi] nm")

    wavlen_sub = wavlen[idx_w]
    band_sub = band[idx_b]
    rsr_sub = rsr_all[idx_w, idx_b]
    if clip_negative
        rsr_sub = max.(rsr_sub, 0.0)
    end

    kernel = PACE_SIF.KernelInstrument(
        band_sub,
        wavlen_sub,
        rsr_sub,
        collect(Float64.(λ_hres)),
        collect(Float64.(spectral_axis)),
    )

    return kernel, (
        n_wavelength_samples = length(wavlen_sub),
        n_bands = length(band_sub),
        λ_range_nm = (minimum(band_sub), maximum(band_sub)),
    )
end

"""
    kernel_quality_metrics(kernel)

Compute kernel sanity metrics:
- row-sum stats
- negative weight count
- weighted-center offset from nominal band center
- peak-location offset from nominal band center
"""
function kernel_quality_metrics(kernel)
    row_sum = vec(sum(kernel.RSR_out, dims=2))
    n_negative = count(<(0.0), kernel.RSR_out)
    total_weights = length(kernel.RSR_out)

    λ = collect(Float64.(kernel.wvlen_out))
    band = collect(Float64.(kernel.band))
    n_band = length(band)

    center = zeros(Float64, n_band)
    peak = zeros(Float64, n_band)
    for i in 1:n_band
        row = kernel.RSR_out[i, :]
        s = sum(row)
        center[i] = s > 0 ? sum(row .* λ) / s : NaN
        _, j = findmax(row)
        peak[i] = λ[j]
    end

    center_offset = center .- band
    peak_offset = peak .- band

    return (
        n_band = n_band,
        band_nm = band,
        center_nm = center,
        peak_nm = peak,
        row_sum_min = minimum(row_sum),
        row_sum_max = maximum(row_sum),
        row_sum_mean = sum(row_sum) / length(row_sum),
        negative_weight_fraction = n_negative / total_weights,
        center_offset_mean_nm = sum(center_offset) / length(center_offset),
        center_offset_abs_max_nm = maximum(abs.(center_offset)),
        peak_offset_abs_max_nm = maximum(abs.(peak_offset)),
    )
end

"""
    compare_kernels(reference_kernel, candidate_kernel; band_tol_nm=1e-3)

Compare two kernels on shared bands by interpolating reference rows onto
candidate wavelength grid and reporting row-wise differences.
"""
function compare_kernels(
    reference_kernel,
    candidate_kernel;
    band_tol_nm::Float64=1e-3,
)
    ref_band = collect(Float64.(reference_kernel.band))
    cand_band = collect(Float64.(candidate_kernel.band))
    ref_λ = collect(Float64.(reference_kernel.wvlen_out))
    cand_λ = collect(Float64.(candidate_kernel.wvlen_out))

    ref_idx = Int[]
    cand_idx = Int[]
    band_delta = Float64[]
    for (i, b) in enumerate(cand_band)
        j = argmin(abs.(ref_band .- b))
        δ = abs(ref_band[j] - b)
        if δ <= band_tol_nm
            push!(cand_idx, i)
            push!(ref_idx, j)
            push!(band_delta, δ)
        end
    end

    n_common = length(cand_idx)
    n_common > 0 || return (
        n_common = 0,
        band_max_delta_nm = Inf,
        l1_mean = Inf,
        l1_max = Inf,
        l2_mean = Inf,
        l2_max = Inf,
        max_abs_mean = Inf,
        max_abs_max = Inf,
    )

    l1 = zeros(Float64, n_common)
    l2 = zeros(Float64, n_common)
    maxabs = zeros(Float64, n_common)
    for k in 1:n_common
        i = cand_idx[k]
        j = ref_idx[k]

        row_cand = vec(candidate_kernel.RSR_out[i, :])
        interp_ref = LinearInterpolation(ref_λ, vec(reference_kernel.RSR_out[j, :]), extrapolation_bc=0.0)
        row_ref = interp_ref.(cand_λ)

        sc = sum(row_cand)
        sr = sum(row_ref)
        sc > 0 && (row_cand ./= sc)
        sr > 0 && (row_ref ./= sr)

        d = row_cand .- row_ref
        l1[k] = sum(abs.(d))
        l2[k] = sqrt(sum(d .^ 2))
        maxabs[k] = maximum(abs.(d))
    end

    return (
        n_common = n_common,
        band_max_delta_nm = maximum(band_delta),
        l1_mean = sum(l1) / n_common,
        l1_max = maximum(l1),
        l2_mean = sum(l2) / n_common,
        l2_max = maximum(l2),
        max_abs_mean = sum(maxabs) / n_common,
        max_abs_max = maximum(maxabs),
    )
end

"""
    prepare_mwe_inputs(config_path)

Loads config, LUTs, SIF basis, and kernel.

Kernel flow:
- always regenerates from PACE RSR NetCDF on LUT-aligned high-res grid
- optionally compares against stored `kernel_file`
- optionally saves regenerated kernel
"""
function prepare_mwe_inputs(config_path::AbstractString)
    cfg = read_mwe_config(config_path)
    paths = resolve_paths(cfg)

    λ_min = Float64(cfg_get(cfg, "spectral", "lambda_min_nm", 680.0))
    λ_max = Float64(cfg_get(cfg, "spectral", "lambda_max_nm", 770.0))
    delta_λ = Float64(cfg_get(cfg, "spectral", "delta_lambda_nm", 0.005))
    sif_nev = Int(cfg_get(cfg, "spectral", "sif_nev", 2))
    normalize_sif = Bool(cfg_get(cfg, "spectral", "normalize_sif_first_ev", true))

    clip_negative_rsr = Bool(cfg_get(cfg, "kernel", "clip_negative_rsr", true))
    compare_with_stored = Bool(cfg_get(cfg, "kernel", "compare_with_stored", true))
    band_match_tol_nm = Float64(cfg_get(cfg, "kernel", "band_match_tol_nm", 0.01))
    save_regenerated = Bool(cfg_get(cfg, "kernel", "save_regenerated", false))
    use_effective_centers = Bool(cfg_get(cfg, "kernel", "use_effective_band_centers", false))

    o2_sitp = PACE_SIF.read_rescale(must_exist(paths.o2_path))
    h2o_sitp = PACE_SIF.read_rescale(must_exist(paths.h2o_path))
    _, λ_lut, axis_unit = infer_lut_spectral_axis(o2_sitp)
    λ_hres = build_highres_wavelength_grid(
        λ_lut;
        λ_min=λ_min,
        λ_max=λ_max,
        delta_lambda_nm=delta_λ,
        strict=true,
    )
    spectral_axis, _ = make_kernel_spectral_axis(o2_sitp, λ_hres)

    regenerated_kernel, generation_info = build_kernel_from_rsr_nc(
        must_exist(paths.pace_rsr_path),
        λ_hres,
        spectral_axis;
        λ_min=λ_min,
        λ_max=λ_max,
        clip_negative=clip_negative_rsr,
    )
    regenerated_quality = kernel_quality_metrics(regenerated_kernel)

    if save_regenerated
        MyKernel = regenerated_kernel
        @save paths.regen_kernel_path MyKernel
    end

    comparison = nothing
    if compare_with_stored && isfile(paths.kernel_path)
        stored_kernel = let
            @load paths.kernel_path MyKernel
            MyKernel
        end
        stored_kernel = subset_kernel(stored_kernel, λ_min, λ_max)
        comparison = compare_kernels(
            stored_kernel,
            regenerated_kernel;
            band_tol_nm=band_match_tol_nm,
        )
    end

    # --- SIF basis on full high-resolution forward-model grid ---
    sif_basis_hres = load_sif_basis(
        must_exist(paths.sif_path),
        λ_hres;
        nEV=sif_nev,
        normalize=normalize_sif,
    )

    # Convolve high-resolution SIF basis to nominal low-resolution PACE bands.
    band_nominal = collect(Float64.(regenerated_kernel.band))
    sif_basis_lres_nominal = regenerated_kernel.RSR_out * sif_basis_hres

    λ = use_effective_centers ?
        collect(Float64.(regenerated_quality.center_nm)) :
        band_nominal
    λ_source = use_effective_centers ? "kernel_weighted_center" : "pace_band_center"

    # If using effective centers, remap low-res SIF basis from nominal centers.
    if use_effective_centers
        sif_basis_lres = zeros(Float64, length(λ), size(sif_basis_lres_nominal, 2))
        for iev in 1:size(sif_basis_lres_nominal, 2)
            itp_ev = LinearInterpolation(
                band_nominal,
                vec(sif_basis_lres_nominal[:, iev]),
                extrapolation_bc=Line(),
            )
            sif_basis_lres[:, iev] .= itp_ev.(λ)
        end
    else
        sif_basis_lres = sif_basis_lres_nominal
    end

    λc = PACE_SIF.center_wavelength(collect(Union{Missing, Float64}, λ))

    return (
        cfg = cfg,
        paths = paths,
        λ = λ,
        λc = λc,
        kernel = regenerated_kernel,
        regenerated_kernel = regenerated_kernel,
        regenerated_kernel_quality = regenerated_quality,
        generation_info = generation_info,
        kernel_comparison = comparison,
        lambda_source = λ_source,
        delta_lambda_nm = delta_λ,
        o2_sitp = o2_sitp,
        h2o_sitp = h2o_sitp,
        sif_basis_hres = sif_basis_hres,
        sif_basis = sif_basis_lres,
        axis_unit = axis_unit,
        spectral_axis = spectral_axis,
        λ_hres = λ_hres,
    )
end

end # module
