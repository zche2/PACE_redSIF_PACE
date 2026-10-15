# Generate an ensemble of TOA Stokes-I spectra over a LAMBERTIAN ocean surface, with
# aerosol optics from Google-RT/GCHPIO (TOMAS bin-resolved Mie) instead of the YAML species.
#
# Same design as generate_test_spec_lambertian.jl: same seed → same MERRA T/p/q profiles,
# GCHP ocean columns, SIF, geometry and Lambertian albedos, the same OCI convolution and SNR
# noise model, and the same output layout (plus `aod` on OCI bands). Only the aerosol optics differ:
#
#   generate_test_spec_lambertian.jl   one lognormal per YAML species (column-mean wet radius,
#                                      pure-component RI), Mie at the band ends, one phase
#                                      function per species for the whole column
#   this script                        GCHPIO reads the TOMAS state (15 bins × species mass, number,
#                                      aerosol water) of the same GCHP cell, reconstructs the
#                                      within-bin size distribution, mixes RI per bin and layer,
#                                      and runs Mie per (bin, layer) → per-layer τ, ω̃, Greek
#
# Per GCHP column (cached in memory and in OPTICS_CACHE, so a resumed run does not redo Mie):
#   - GCHPIO.Aerosols.compute_aerosol_optics_band on N_QUERY wavelengths spanning the RT grid
#     (Mie at N_NODES Chebyshev nodes, cubic-spline in between; GCHPIO auto-raises the node count
#     for bands wider than 0.15 µm, e.g. 7 nodes for the YAML's 623–850 nm grid)
#   - per GCHP layer: τ(λ) at the query λ, band-mean ω̃ and τ_sca-weighted band-mean Greek coefs
#   Hygroscopic growth uses the GCHP cell's own T/q (as the species path does via RadiWet).
#
# Per sample:
#   - GCHP layers (72, TOA→BOA) are mapped onto the MERRA RT p_half by log-pressure overlap
#     (each GCHP layer's τ conserved; layers outside the RT grid go to the top/bottom RT layer)
#   - RT layer j: τ_j(ν) = Σ_i W_ji τ_i(λ) (linear in λ between query points), ω̃_j and Greek_j
#     τ_sca-weighted over the contributing GCHP layers
#   - one vSmartMOM aerosol mode per active RT layer (diagonal τ_aer), as in
#     Google-RT SURF.build_layer_resolved_model, but injected into the BatchContext model so gas
#     absorption is not rebuilt; phase truncation follows the YAML `truncation:` (override with
#     AER_TRUNCATION=dbge|none). Note ocean_lambertian_0929.yaml sets NoTruncation().
#
#   julia --project=surrogate_meas/full_RT_construction/gchpio_env -t 16 \
#         surrogate_meas/full_RT_construction/generate_test_spec_lambertian_gchpio.jl
#   N_SAMPLES=2 OUT_NC=/path/test.nc julia --project=... -t 16 ...
#   N_SAMPLES=1000 N_RUN=3 OUT_NC=/path/test.nc ...   (first 3 samples of the 1000-sample design,
#                                                       i.e. the same draws as the YAML ensemble)
#   MIE_ARCH=GPU ...           (per-bin Mie on CUDA; default CPU threads, start julia with -t N)
#   OCEAN_YAML=/path/other.yaml OUT_NC=/path/other_out.nc ...   (new OUT_NC, or it resumes the old file)
#
# Environment: gchpio_env/ pins vSmartMOM to GCHPIO's rev (feat/mie-node-api 2cd6934 = main + 21).
# The noise realisation can differ from the YAML ensemble if any package draws from the global RNG.

using vSmartMOM
using vSmartMOM.SolarModel
using vSmartMOM.Scattering: Stokes_I, NAI2, NoTruncation, δBGE, GreekCoefs, AerosolOptics,
    truncate_phase
using vSmartMOM.CoreRT: AerosolState, Optics, SolverConfig, RTModel
using GCHPIO
using GCHPIO.Aerosols: compute_aerosol_optics_band
using NCDatasets
using JLD2
using Interpolations
using Random
using Statistics
using DelimitedFiles

include(joinpath(@__DIR__, "..", "..", "src", "tools", "Instrument.jl"))

# Data files and configs/ are not in git; default to the main checkout's copies.
const REPO_RT = "/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/full_RT_construction"
const OCEAN_YAML = get(ENV, "OCEAN_YAML",
    joinpath(REPO_RT, "..", "configs", "ocean_lambertian_0929.yaml"))
const PACE_RSR_NC = get(ENV, "PACE_RSR_NC",
    "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/PACE_OCI_RSRs.nc")
const SIF_LIB = get(ENV, "SIF_LIB",
    "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/SIF_singular_vector.jld2")
const OCEAN_COLS_NC = get(ENV, "OCEAN_COLS_NC",
    joinpath(REPO_RT, "output_aerosol_profiles", "gchp_ocean_columns_n500.nc"))
const MERRA_SEA_TPQ_NC = get(ENV, "MERRA_SEA_TPQ_NC",
    joinpath(REPO_RT, "merra2_sea_Tpq_columns_n1500.nc"))
const SNR_FILE = get(ENV, "SNR_FILE",
    "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt")
const OUT_NC = get(ENV, "OUT_NC",
    joinpath(REPO_RT, "output_new_Lambertian", "rt_toa_ensemble_lambertian_gchpio.nc"))
# GCHP file the ocean columns were sampled from (its `source` attribute) unless overridden.
const GCHP_NC = get(ENV, "GCHP_NC", "")
const OPTICS_CACHE = get(ENV, "OPTICS_CACHE", replace(OUT_NC, r"\.nc$" => "") * "_aerosol_optics.jld2")

const N_SAMPLES = parse(Int, get(ENV, "N_SAMPLES", "500"))
# Compute at most N_RUN samples in this invocation (0 = all); the design is still drawn for N_SAMPLES.
const N_RUN = parse(Int, get(ENV, "N_RUN", "0"))
const ENABLE_AEROSOLS = parse(Bool, get(ENV, "ENABLE_AEROSOLS", "true"))
const SEED = parse(Int, get(ENV, "ENSEMBLE_SEED", "20260913"))
const SIF_λ = 678.2
const SIF_STRENGTH = (0.0, 0.5)          # W m⁻² sr⁻¹ μm⁻¹ at SIF_λ
const SZA_RANGE = (5.0, 70.0)            # deg
const VZA_RANGE = (0.0, 60.0)            # deg; OCI-like swath
const VAZ_RANGE = (0.0, 180.0)           # deg; relative azimuth (glint vs dark)
const WIND_SPEED_RANGE = (0.0, 10.0)     # m/s
const WHITECAP_ALBEDO_RANGE = (0.1, 0.5)  # 0-1
const INCLUDE_WHITECAPS_RANGE = (0, 1)    # 0 or 1
const LAMBERT_ALBEDO_RANGE = (0.005, 0.025)   # Lambertian ocean reflectance, U(lo, hi)
const ALBEDO_SEED_OFFSET = 1                  # albedo RNG seed = SEED + offset
# Sea Tpq NC was built with PROFILE_STRIDE=3 → 24 layers / 25 half-levels.
const PROFILE_STRIDE = parse(Int, get(ENV, "PROFILE_STRIDE", "3"))
# Force CPU if GPUs are busy: ARCH=CPU julia generate_test_spec.jl
const ARCH = uppercase(get(ENV, "ARCH", "DEFAULT"))

# GCHPIO aerosol optics
const N_NODES = parse(Int, get(ENV, "N_NODES", "3"))      # Mie λ-nodes (GCHPIO may raise it)
const N_QUERY = parse(Int, get(ENV, "N_QUERY", "24"))     # λ where τ is evaluated over the RT grid
const MIE_ARCH = uppercase(get(ENV, "MIE_ARCH", "CPU"))   # CPU (threads) or GPU
const AER_TRUNCATION = lowercase(get(ENV, "AER_TRUNCATION", "yaml"))   # yaml | dbge | none

const h = 6.62607015e-34
const c_light = 299792458.0

# Lambertian surfaces inject SurfaceSIF through vSmartMOM's built-in
# `inject_surface_SIF!` (same factor-2 convention as the Cox-Munk hook in
# generate_test_spec.jl), so no override is needed here.

"New observation geometry + surface sharing the updated atmosphere/optics."
function with_geometry(model, params, sza, vza, vaz, surf)
    FT = params.float_type
    geom = vSmartMOM.CoreRT.ObsGeometry{FT}(
        FT(sza), FT[vza], FT[vaz], params.obs_alt,
    )
    qp = vSmartMOM.CoreRT.rt_set_streams(
        params.quadrature_type, params.l_trunc, geom,
        params.polarization_type, array_type(params.architecture),
    )
    surfaces = [surf for _ in 1:length(model.surfaces)]
    return vSmartMOM.CoreRT.RTModel(
        model.architecture, model.solver, model.numerics,
        geom, qp, model.atmosphere, model.optics, surfaces, model.sources,
    )
end

make_lambertian(FT, albedo) = vSmartMOM.CoreRT.LambertianSurfaceScalar(FT(albedo))

"Per-sample Lambertian albedos from their own RNG (does not touch the design stream)."
function draw_albedo(n)
    rng = MersenneTwister(SEED + ALBEDO_SEED_OFFSET)
    lo, hi = LAMBERT_ALBEDO_RANGE
    return lo .+ (hi - lo) .* rand(rng, n)
end

"""Load sea-surface MERRA T/p/q sample `i` (TOA→BOA, p_half increasing).

`merra2_sea_Tpq_columns_*.nc` is written (layer|half, sample) but NCDatasets
exposes arrays as (sample, layer|half).
"""
function sea_tpq_profile(ds, i)
    T = Float64.(ds["T"][i, :])
    q = Float64.(ds["q"][i, :])
    p_half = Float64.(ds["p_half"][i, :])
    ps = Float64(ds["ps"][i])
    length(p_half) == length(T) + 1 || error(
        "Sea Tpq sample $i: n_half=$(length(p_half)) ≠ n_layer+1=$(length(T)+1)")
    p_half[1] < p_half[end] || error(
        "Sea Tpq sample $i: p_half must be TOA→BOA (increasing); got $(p_half[1]) → $(p_half[end])")
    all(diff(p_half) .> 0) || error("Sea Tpq sample $i: p_half not strictly increasing")
    return (T=T, q=max.(q, 0.0), p_half=p_half, ps=ps)
end

"""Install sample-1 sea T/p/q on `params` so BatchContext allocates the right Nz."""
function apply_sea_atmosphere!(params, sea_ds)
    FT = params.float_type
    src0 = sea_tpq_profile(sea_ds, 1)
    params.p = FT.(src0.p_half)
    params.T = FT.(src0.T)
    params.q = FT.(src0.q)
    params.profile_reduction_n = -1   # already at the target layer count
    n_layer = length(src0.T)
    file_stride = Int(get(sea_ds.attrib, "profile_stride", PROFILE_STRIDE))
    file_stride == PROFILE_STRIDE || @warn "Sea Tpq NC profile_stride=$file_stride " *
        "≠ PROFILE_STRIDE=$PROFILE_STRIDE (using file layers)"
    println("RT atmosphere: $n_layer layers from $MERRA_SEA_TPQ_NC " *
            "(TOA→BOA, p_half $(round(src0.p_half[1]; digits=3)) → " *
            "$(round(src0.p_half[end]; digits=1)) hPa)")
    return src0.p_half
end

function oci_kernel(λ_hres, ν_asc)
    ds = NCDataset(PACE_RSR_NC)
    wavlen = collect(Float64.(ds["wavelength"][:]))
    band = collect(Float64.(ds["bands"][:]))
    rsr = collect(Float64.(ds["RSR"][:, :]))
    close(ds)
    λ_lo, λ_hi = extrema(λ_hres)
    idx_w = findall(λ_lo .< wavlen .< λ_hi)
    idx_b = findall(λ_lo .< band .< λ_hi)
    isempty(idx_b) && error("No OCI bands inside the RT grid")
    return Instrument.KernelInstrument(
        band[idx_b], wavlen[idx_w], max.(rsr[idx_w, idx_b], 0.0),
        collect(Float64.(λ_hres)), collect(Float64.(ν_asc)),
    )
end

function snr_coeffs(λ_oci)
    lines = readlines(SNR_FILE)
    header_end = findfirst(line -> occursin("/end_header", line), lines)
    data = readdlm(SNR_FILE, String; skipstart=isnothing(header_end) ? 0 : header_end)
    red = findall(data[:, 1] .== "Red")
    isempty(red) && error("No Red rows in $SNR_FILE")
    λ = parse.(Float64, data[red, 2])
    c1 = parse.(Float64, data[red, 4])
    c2 = parse.(Float64, data[red, 5])
    p = sortperm(λ)
    λs, c1s, c2s = λ[p], c1[p], c2[p]
    # Drop duplicate LUT wavelengths (keep last) so LinearInterpolation knots are unique.
    keep = trues(length(λs))
    for i in 2:length(λs)
        if λs[i] == λs[i - 1]
            keep[i - 1] = false
        end
    end
    λs, c1s, c2s = λs[keep], c1s[keep], c2s[keep]
    issorted(λs; lt=<) || error("SNR wavelengths not strictly increasing after dedupe")
    itp1 = LinearInterpolation(λs, c1s; extrapolation_bc=Flat())
    itp2 = LinearInterpolation(λs, c2s; extrapolation_bc=Flat())
    return itp1.(λ_oci), itp2.(λ_oci)
end

"Library shapes as water-leaving radiance with I(678 nm) = 1, in wavenumber order."
function unit_sif_shapes(ν, λ_lib, shapes)
    λ_lib = Float64.(λ_lib)
    length(λ_lib) >= 2 || error("SIF library needs ≥2 wavelengths")
    all(diff(λ_lib) .> 0) || error(
        "SIF library wavelengths must be strictly increasing (no duplicate knots)")
    dλ = λ_lib[2] - λ_lib[1]
    # Uniform grid → range knots (guaranteed unique). Non-uniform → use λ_lib as-is.
    knots = all(x -> abs(x - dλ) ≤ 1e-12 * max(abs(dλ), 1.0), diff(λ_lib)) ?
        range(λ_lib[1]; step=dλ, length=length(λ_lib)) : λ_lib
    λ_model = 1e7 ./ ν
    nspec = length(ν)
    nlib = size(shapes, 2)
    out = zeros(nspec, nlib)
    for j in 1:nlib
        itp = CubicSplineInterpolation(knots, shapes[:, j]; extrapolation_bc=Line())
        y = itp.(λ_model)
        y[(λ_model .< λ_lib[1]) .| (λ_model .> λ_lib[end])] .= 0.0
        y678 = itp(SIF_λ)
        y678 != 0 || error("SIF shape $j is zero at $(SIF_λ) nm")
        out[:, j] .= y ./ y678
    end
    return out
end

function draw_design(n, n_profiles, n_sif, n_cols)
    wc_lo, wc_hi = INCLUDE_WHITECAPS_RANGE
    return (
        profile_index = rand(1:n_profiles, n),
        sif_index = rand(1:n_sif, n),
        sif_678 = SIF_STRENGTH[1] .+ (SIF_STRENGTH[2] - SIF_STRENGTH[1]) .* rand(n),
        sza = SZA_RANGE[1] .+ (SZA_RANGE[2] - SZA_RANGE[1]) .* rand(n),
        vza = VZA_RANGE[1] .+ (VZA_RANGE[2] - VZA_RANGE[1]) .* rand(n),
        vaz = VAZ_RANGE[1] .+ (VAZ_RANGE[2] - VAZ_RANGE[1]) .* rand(n),
        wind_speed = WIND_SPEED_RANGE[1] .+ (WIND_SPEED_RANGE[2] - WIND_SPEED_RANGE[1]) .* rand(n),
        whitecap_albedo = WHITECAP_ALBEDO_RANGE[1] .+
            (WHITECAP_ALBEDO_RANGE[2] - WHITECAP_ALBEDO_RANGE[1]) .* rand(n),
        include_whitecaps = rand(wc_lo:wc_hi, n),
        column_index = ENABLE_AEROSOLS ? rand(1:n_cols, n) : zeros(Int, n),
    )
end

function create_output(path, λ_oci, p_full, design, albedo, n_layer)
    mkpath(dirname(path))
    n = length(design.sza)
    n_band = length(λ_oci)
    ds = NCDataset(path, "c")
    defDim(ds, "band", n_band)
    defDim(ds, "sample", n)
    defDim(ds, "layer", n_layer)
    defVar(ds, "wavelength", Float64, ("band",); attrib=Dict(
        "units" => "nm", "long_name" => "OCI band center"))
    defVar(ds, "radiance_clean", Float32, ("band", "sample"); attrib=Dict(
        "units" => "W m-2 sr-1 um-1", "long_name" => "TOA Stokes I, OCI-convolved, no noise"),
        fillvalue=Float32(NaN))
    defVar(ds, "radiance_noisy", Float32, ("band", "sample"); attrib=Dict(
        "units" => "W m-2 sr-1 um-1", "long_name" => "TOA Stokes I plus OCI SNR white noise"),
        fillvalue=Float32(NaN))
    defVar(ds, "sigma_noise", Float32, ("band", "sample"); attrib=Dict(
        "units" => "W m-2 sr-1 um-1", "long_name" => "σ from σ² = c1 + c2·R"))
    defVar(ds, "sif_waterleaving", Float32, ("band", "sample"); attrib=Dict(
        "units" => "W m-2 sr-1 um-1", "long_name" => "water-leaving SIF on OCI bands"))
    defVar(ds, "sif_678", Float32, ("sample",); attrib=Dict(
        "units" => "W m-2 sr-1 um-1", "long_name" => "water-leaving SIF at 678 nm"))
    defVar(ds, "sif_library_index", Int32, ("sample",))
    defVar(ds, "profile_index", Int32, ("sample",); attrib=Dict(
        "long_name" => "1-based index in $MERRA_SEA_TPQ_NC"))
    defVar(ds, "sza", Float32, ("sample",); attrib=Dict("units" => "degree"))
    defVar(ds, "vza", Float32, ("sample",); attrib=Dict("units" => "degree"))
    defVar(ds, "vaz", Float32, ("sample",); attrib=Dict(
        "units" => "degree", "long_name" => "relative azimuth, vSmartMOM convention"))
    defVar(ds, "surface_albedo", Float32, ("sample",); attrib=Dict(
        "units" => "1", "long_name" => "Lambertian ocean surface albedo (used in the RT)"))
    defVar(ds, "wind_speed", Float32, ("sample",); attrib=Dict(
        "units" => "m s-1", "long_name" => "Cox-Munk wind draw (NOT used; kept for alignment)"))
    defVar(ds, "whitecap_albedo", Float32, ("sample",); attrib=Dict(
        "units" => "1", "long_name" => "Cox-Munk whitecap albedo draw (NOT used)"))
    defVar(ds, "include_whitecaps", Int8, ("sample",); attrib=Dict(
        "long_name" => "Cox-Munk whitecap switch draw (NOT used)"))
    defVar(ds, "column_index", Int32, ("sample",); attrib=Dict(
        "long_name" => "1-based GEOS-Chem ocean column index; 0 if aerosols disabled"))
    defVar(ds, "p", Float32, ("layer", "sample"); attrib=Dict(
        "units" => "hPa", "long_name" => "full-level pressure used in the RT"))
    defVar(ds, "T", Float32, ("layer", "sample"); attrib=Dict("units" => "K"))
    defVar(ds, "q", Float32, ("layer", "sample"); attrib=Dict(
        "units" => "kg kg-1", "long_name" => "specific humidity used in the RT"))
    ds["wavelength"][:] = λ_oci
    ds["sif_678"][:] = Float32.(design.sif_678)
    ds["sif_library_index"][:] = Int32.(design.sif_index)
    ds["profile_index"][:] = Int32.(design.profile_index)
    ds["sza"][:] = Float32.(design.sza)
    ds["vza"][:] = Float32.(design.vza)
    ds["vaz"][:] = Float32.(design.vaz)
    ds["surface_albedo"][:] = Float32.(albedo)
    ds["wind_speed"][:] = Float32.(design.wind_speed)
    ds["whitecap_albedo"][:] = Float32.(design.whitecap_albedo)
    ds["include_whitecaps"][:] = Int8.(design.include_whitecaps)
    ds["column_index"][:] = Int32.(design.column_index)
    ds.attrib["n_completed"] = 0
    ds.attrib["n_samples"] = n
    ds.attrib["seed"] = SEED
    ds.attrib["sif_lambda_nm"] = SIF_λ
    ds.attrib["enable_aerosols"] = Int8(ENABLE_AEROSOLS)
    ds.attrib["ocean_columns_nc"] = OCEAN_COLS_NC
    ds.attrib["ocean_yaml"] = OCEAN_YAML
    ds.attrib["merra_sea_tpq_nc"] = MERRA_SEA_TPQ_NC
    ds.attrib["profile_stride"] = PROFILE_STRIDE
    ds.attrib["vertical_order"] = "TOA→BOA (p_half increasing)"
    ds.attrib["radiance"] = "TOA Stokes I from Lambertian surface + SurfaceSIF, OCI-convolved"
    ds.attrib["surface"] = "LambertianSurfaceScalar"
    ds.attrib["lambert_albedo_range"] = collect(LAMBERT_ALBEDO_RANGE)
    ds.attrib["albedo_seed"] = SEED + ALBEDO_SEED_OFFSET
    ds.attrib["pressure_template_hpa"] = join(string.(p_full), ",")
    return ds
end


# ─────────────────────────────────────────────────────────────────────────────
# GCHPIO aerosol optics
# ─────────────────────────────────────────────────────────────────────────────

"Ocean-column table: GCHP cube-sphere indices (0-based in the file) and lat/lon."
function column_table()
    NCDataset(OCEAN_COLS_NC) do d
        gchp = isempty(GCHP_NC) ? String(d.attrib["source"]) : GCHP_NC
        (i_face=Int.(d["i_face"][:]), i_y=Int.(d["i_y"][:]), i_x=Int.(d["i_x"][:]),
         lat=Float64.(d["lat"][:]), lon=Float64.(d["lon"][:]), gchp=gchp)
    end
end

"Truncation applied to the injected aerosol phase functions."
function aerosol_truncation(params)
    FT = params.float_type
    l_cap = params.stream_l_cap > 0 ? params.stream_l_cap : params.l_trunc
    AER_TRUNCATION == "dbge" && return δBGE{FT}(l_cap, FT(params.Δ_angle))
    AER_TRUNCATION == "none" && return NoTruncation()
    AER_TRUNCATION == "yaml" || error("AER_TRUNCATION must be yaml, dbge or none")
    t = params.truncation
    t isa vSmartMOM.Scattering.AutoTruncation && return δBGE{FT}(l_cap, FT(params.Δ_angle))
    return t
end

_greek_matrix(g) = hcat(g.α, g.β, g.γ, g.δ, g.ϵ, g.ζ)

"""GCHP-layer aerosol optics for ocean column `icol`, from the TOMAS state of its GCHP cell.

Returns `(; p_half, λq, τ, ω, G)` on the native GCHP layers (TOA→BOA, hPa):
`τ[iz, k]` extinction optical depth at `λq[k]` (µm), `ω[iz]` band-mean single-scattering albedo,
`G[iz]` (L × 6: α β γ δ ϵ ζ) τ_sca-weighted band-mean Greek coefficients (0 × 6 if no aerosol).
"""
function compute_column_optics(gf, cols, icol, λq)
    iface, iy, ix = cols.i_face[icol], cols.i_y[icol], cols.i_x[icol]
    scene = scene_at(gf, ix + 1, iy + 1, iface + 1)
    (abs(scene.lat - cols.lat[icol]) < 0.05 && abs(mod(scene.lon - cols.lon[icol] + 180, 360) - 180) < 0.05) ||
        error("Column $icol: GCHP cell ($(scene.lat), $(scene.lon)) ≠ file ($(cols.lat[icol]), $(cols.lon[icol]))")
    scene.aerosols === nothing && error("GCHP file has no TOMAS aerosols for column $icol")
    mie_arch = MIE_ARCH == "GPU" ? vSmartMOM.Architectures.GPU() : nothing
    band = compute_aerosol_optics_band(
        scene.aerosols, λq, Stokes_I{Float64}(), NAI2();
        p_half=scene.p_half, T=scene.T, q=scene.q,
        n_nodes=N_NODES, architecture=mie_arch, trunc_type=NoTruncation())
    nlay, nq = length(scene.T), length(λq)
    τ = zeros(nlay, nq)
    ω = zeros(nlay)
    G = [zeros(0, 6) for _ in 1:nlay]
    for iz in 1:nlay
        τ_sum = 0.0
        sca_sum = 0.0
        acc = zeros(0, 6)
        for k in 1:nq
            o = band.per_λ[k].optics[iz]
            o === nothing && continue
            τk = Float64(band.per_λ[k].τ_layer[iz])
            τk > 0 || continue
            τ[iz, k] = τk
            w = τk * Float64(o.ω̃)
            τ_sum += τk
            sca_sum += w
            g = _greek_matrix(o.greek_coefs)
            if size(g, 1) > size(acc, 1)
                acc = vcat(acc, zeros(size(g, 1) - size(acc, 1), 6))
            end
            @views acc[1:size(g, 1), :] .+= w .* g
        end
        if sca_sum > 0
            ω[iz] = sca_sum / τ_sum
            G[iz] = acc ./ sca_sum
        end
    end
    return (; p_half=Float64.(scene.p_half), λq=collect(Float64, λq), τ, ω, G,
              node_λs=collect(Float64, band.node_λs), lat=Float64(scene.lat), lon=Float64(scene.lon))
end

"Column optics from the in-memory cache, then OPTICS_CACHE, else computed (and stored in both)."
function column_optics!(mem, gf, cols, icol, λq)
    haskey(mem, icol) && return mem[icol]
    key = "col$(icol)_n$(N_NODES)_q$(length(λq))_$(round(λq[1]; digits=5))-$(round(λq[end]; digits=5))"
    col = isfile(OPTICS_CACHE) ? jldopen(f -> haskey(f, key) ? f[key] : nothing, OPTICS_CACHE, "r") : nothing
    if col === nothing
        t = @elapsed col = compute_column_optics(gf, cols, icol, λq)
        jldopen(f -> (f[key] = col), OPTICS_CACHE, "a+")
        println("    GCHPIO optics column $icol: $(round(t; digits=1)) s, " *
                "Mie nodes $(round.(col.node_λs; digits=4)) µm, column τ $(round(sum(col.τ[:, 1]); digits=4))" *
                " → $(round(sum(col.τ[:, end]); digits=4)) ($(round(λq[1]; digits=3))→$(round(λq[end]; digits=3)) µm)")
    end
    mem[icol] = col
    return col
end

"""Fraction of each source layer (columns) in each RT layer (rows), by log-pressure overlap.

Both half-level vectors are TOA→BOA in hPa. Every source column sums to 1, so τ is conserved
layer by layer; source layers outside the RT grid go to the nearest (top or bottom) RT layer."""
function overlap_weights(ph_src, ph_rt)
    ls, lr = log.(ph_src), log.(ph_rt)
    ns, nr = length(ls) - 1, length(lr) - 1
    W = zeros(nr, ns)
    for i in 1:ns
        a, b = max(ls[i], lr[1]), min(ls[i + 1], lr[end])
        if b > a
            for j in 1:nr
                o = min(b, lr[j + 1]) - max(a, lr[j])
                o > 0 && (W[j, i] = o)
            end
            W[:, i] ./= sum(W[:, i])
        else
            W[clamp(searchsortedlast(lr, ls[i]), 1, nr), i] = 1.0
        end
    end
    return W
end

"""RTModel = `model` (gases, Rayleigh, geometry from the BatchContext) + GCHPIO aerosols on the
RT layers given by `p_half_rt`: one aerosol mode per active RT layer (diagonal τ_aer)."""
function inject_aerosols(model, params, col, p_half_rt, trunc)
    FT = params.float_type
    nlay = length(p_half_rt) - 1
    W = overlap_weights(col.p_half, p_half_rt)
    τq = W * col.τ                                      # RT layer × query λ
    τ̄_src = vec(mean(col.τ; dims=2))
    sca = W .* (τ̄_src .* col.ω)'                         # band-mean τ_sca, RT layer × GCHP layer
    ext = W .* τ̄_src'
    active = [j for j in 1:nlay if sum(sca[j, :]) > 0]
    isempty(active) && return model

    aos = AerosolOptics{FT}[]
    for j in active
        L = maximum(size(col.G[i], 1) for i in eachindex(col.G) if sca[j, i] > 0)
        g = zeros(L, 6)
        for i in eachindex(col.G)
            sca[j, i] > 0 || continue
            @views g[1:size(col.G[i], 1), :] .+= sca[j, i] .* col.G[i]
        end
        g ./= sum(sca[j, :])
        ao = AerosolOptics{FT}(
            greek_coefs=GreekCoefs((FT.(g[:, c]) for c in 1:6)...),
            ω̃=FT(sum(sca[j, :]) / sum(ext[j, :])), k=FT(sum(ext[j, :])), fᵗ=zero(FT))
        push!(aos, trunc isa δBGE && length(ao.greek_coefs.β) > trunc.l_max ?
            truncate_phase(trunc, ao) : truncate_phase(NoTruncation(), ao))
    end

    n_bands = length(model.atmosphere.spec_bands)
    aerosol_optics = [aos for _ in 1:n_bands]
    τ_aer = map(1:n_bands) do ib
        ν = Float64.(collect(model.atmosphere.spec_bands[ib]))
        λ_um = 1e4 ./ ν
        τb = zeros(FT, length(aos), length(ν), nlay)
        for (m, j) in enumerate(active)
            itp = LinearInterpolation(col.λq, τq[j, :]; extrapolation_bc=Flat())
            τb[m, :, j] .= FT.(max.(itp.(λ_um), 0.0))
        end
        τb
    end

    l_max = [maximum(length(ao.greek_coefs.β) for ao in aerosol_optics[ib]) for ib in 1:n_bands]
    comps = [vSmartMOM.CoreRT._band_components(params, aerosol_optics, model.sources, ib)
             for ib in 1:n_bands]
    m_max = vSmartMOM.CoreRT._derive_m_max_bands_via_traits(
        l_max, params.max_m, comps, model.quad_points.Nstreams)
    sc = model.solver
    solver = SolverConfig{FT, typeof(sc.polarization_type), typeof(sc.quadrature_type)}(
        sc.polarization_type, sc.quadrature_type, m_max, m_max .+ 1, l_max,
        sc.l_trunc, sc.Δ_angle, sc.depol, sc.use_component_traits)
    optics = Optics(model.optics.rayleigh, AerosolState(aerosol_optics, τ_aer),
                    model.optics.τ_abs, model.optics.τ_rayl)
    return RTModel(model.architecture, solver, model.numerics, model.geometry,
                   model.quad_points, model.atmosphere, optics, model.surfaces, model.sources)
end

"Column AOD of `col` interpolated to wavelengths `λ_nm`."
function column_aod(col, λ_nm)
    itp = LinearInterpolation(col.λq, vec(sum(col.τ; dims=1)); extrapolation_bc=Flat())
    return itp.(λ_nm ./ 1e3)
end

function main()
    println("Building RT model from $OCEAN_YAML  (aerosols=$(ENABLE_AEROSOLS), GCHPIO optics)")
    params = parameters_from_yaml(OCEAN_YAML)
    if ARCH == "CPU"
        params.architecture = CPU()
        println("Architecture forced to CPU() via ARCH=CPU")
    elseif ARCH == "GPU"
        params.architecture = GPU()
        println("Architecture forced to GPU() via ARCH=GPU")
    end
    surf0 = only(params.brdf)
    surf0 isa vSmartMOM.CoreRT.LambertianSurfaceScalar ||
        error("Expected LambertianSurfaceScalar in OCEAN_YAML surface:; got $(typeof(surf0))")
    FT_s = typeof(surf0.albedo)
    trunc = aerosol_truncation(params)
    # The base model is gases + Rayleigh only; GCHPIO aerosols are injected per sample.
    params.scattering_params = nothing

    isfile(MERRA_SEA_TPQ_NC) || error("Missing sea T/p/q file: $MERRA_SEA_TPQ_NC")
    sea = NCDataset(MERRA_SEA_TPQ_NC)
    n_profiles = Int(sea.dim["sample"])
    # Replace YAML grid with pre-reduced sea MERRA T/p/q (TOA→BOA).
    p_template = apply_sea_atmosphere!(params, sea)
    n_layer = length(p_template) - 1

    ctx = vSmartMOM.CoreRT.BatchContext(params)
    ν = params.spec_bands[1]
    n_stokes = params.polarization_type.n
    n_to_radiance = @. 100 * h * c_light * ν
    F_sol = SolarModel.default_solar_spectrum_at_earth(ν)[:, 2]
    F₀ = zeros(n_stokes, length(ν))
    F₀[1, :] .= F_sol
    println("polarization: $(typeof(params.polarization_type))  n_stokes=$n_stokes")

    λ_hres = 1e7 ./ reverse(ν)
    kernel = oci_kernel(λ_hres, reverse(ν))
    λ_oci = collect(Float64.(kernel.band))
    c1, c2 = snr_coeffs(λ_oci)
    println("OCI bands: $(length(λ_oci)) in $(extrema(λ_oci)) nm")

    n_cols = 0
    cols = nothing
    gf = nothing
    λq = collect(range(1e4 / maximum(ν), 1e4 / minimum(ν); length=N_QUERY))   # µm
    col_mem = Dict{Int, Any}()
    if ENABLE_AEROSOLS
        isfile(OCEAN_COLS_NC) || error("Missing ocean columns file: $OCEAN_COLS_NC")
        cols = column_table()
        n_cols = length(cols.lat)
        gf = GCHPFile(cols.gchp; aerosol_scheme=:auto)
        println("GCHPIO: $(cols.gchp)  ($n_cols columns)  N_NODES=$N_NODES  N_QUERY=$N_QUERY " *
                "($(round(λq[1]; digits=4))–$(round(λq[end]; digits=4)) µm)  Mie on $MIE_ARCH  " *
                "truncation $(typeof(trunc))  cache $OPTICS_CACHE")
    else
        println("Aerosols disabled; Rayleigh-only optics")
    end

    sif_file = jldopen(SIF_LIB)
    λ_lib = Float64.(sif_file["SIF_wavelen"])
    shapes = Float64.(sif_file["SIF_shapes"])
    close(sif_file)
    sif_unit = unit_sif_shapes(ν, λ_lib, shapes)   # I(678)=1, wavenumber order

    ds, start_i = if isfile(OUT_NC)
        existing = NCDataset(OUT_NC, "a")
        done = Int(get(existing.attrib, "n_completed", 0))
        stored_n = Int(existing.dim["sample"])
        stored_seed = Int(get(existing.attrib, "seed", -1))
        stored_layers = Int(existing.dim["layer"])
        stored_stride = Int(get(existing.attrib, "profile_stride", -1))
        layers_ok = stored_layers == n_layer &&
            (stored_stride == -1 || stored_stride == PROFILE_STRIDE)
        if !layers_ok
            close(existing)
            error("Existing $OUT_NC has layer=$stored_layers (stride=$(stored_stride)), " *
                  "but this run uses layer=$n_layer (PROFILE_STRIDE=$PROFILE_STRIDE). " *
                  "Remove it or set a new OUT_NC.")
        elseif stored_n == N_SAMPLES && stored_seed == SEED && done < N_SAMPLES
            println("Resuming $OUT_NC at sample $(done + 1) / $N_SAMPLES")
            existing, done + 1
        elseif done >= N_SAMPLES && stored_n == N_SAMPLES && stored_seed == SEED
            println("Already complete: $OUT_NC")
            close(existing)
            close(sea)
            return
        else
            close(existing)
            error("Existing $OUT_NC does not match N_SAMPLES=$N_SAMPLES seed=$SEED. Remove it or set OUT_NC.")
        end
    else
        Random.seed!(SEED)
        design = draw_design(N_SAMPLES, n_profiles, size(sif_unit, 2), max(n_cols, 1))
        albedo = draw_albedo(N_SAMPLES)
        out = create_output(OUT_NC, λ_oci, 0.5 .* (p_template[1:end-1] .+ p_template[2:end]),
                            design, albedo, n_layer)
        defVar(out, "aod", Float32, ("band", "sample"); attrib=Dict(
            "units" => "1", "long_name" => "GCHPIO column aerosol optical depth at the OCI band centre"),
            fillvalue=Float32(NaN))
        out.attrib["aerosol_optics"] = "Google-RT GCHPIO TOMAS bin-resolved Mie, layer-resolved"
        out.attrib["gchp_file"] = ENABLE_AEROSOLS ? cols.gchp : ""
        out.attrib["gchpio_n_nodes"] = N_NODES
        out.attrib["gchpio_n_query"] = N_QUERY
        out.attrib["aerosol_truncation"] = string(typeof(trunc))
        out.attrib["aerosol_optics_cache"] = OPTICS_CACHE
        println("Writing $OUT_NC  ($N_SAMPLES samples)")
        out, 1
    end

    t0 = time()
    stop_i = N_RUN > 0 ? min(N_SAMPLES, start_i + N_RUN - 1) : N_SAMPLES
    for i in start_i:stop_i
        ip = Int(ds["profile_index"][i])
        isif = Int(ds["sif_library_index"][i])
        strength = Float64(ds["sif_678"][i])
        sza = Float64(ds["sza"][i])
        vza = Float64(ds["vza"][i])
        vaz = Float64(ds["vaz"][i])
        alb = Float64(ds["surface_albedo"][i])
        icol = Int(ds["column_index"][i])

        src = sea_tpq_profile(sea, ip)
        length(src.T) == n_layer || error(
            "Sea Tpq sample $ip has $(length(src.T)) layers, expected $n_layer")
        p_half, T, q = src.p_half, src.T, src.q
        vSmartMOM.CoreRT.update_model!(ctx; T=T, p_half=p_half, q=q)

        model = ctx.model
        if ENABLE_AEROSOLS
            col = column_optics!(col_mem, gf, cols, icol, λq)
            model = inject_aerosols(ctx.model, params, col, p_half, trunc)
            ds["aod"][:, i] = Float32.(column_aod(col, λ_oci))
        end

        surf = make_lambertian(FT_s, alb)
        scene = with_geometry(model, params, sza, vza, vaz, surf)

        I_wl = sif_unit[:, isif] .* strength
        SIF₀ = zeros(n_stokes, length(ν))
        SIF₀[1, :] .= π .* I_wl ./ n_to_radiance
        sources = SolarBeam(F₀=F₀) + SurfaceSIF(SIF₀=SIF₀)
        R, = rt_run(scene; sources=sources)

        R_I = reverse(Array(R)[1, 1, :] .* n_to_radiance)
        I_oci = vec(kernel.RSR_out * R_I)
        sif_oci = vec(kernel.RSR_out * reverse(I_wl))
        σ = sqrt.(c1 .+ c2 .* max.(I_oci, 0.0))
        noisy = I_oci .+ randn(length(I_oci)) .* σ

        ds["radiance_clean"][:, i] = Float32.(I_oci)
        ds["radiance_noisy"][:, i] = Float32.(noisy)
        ds["sigma_noise"][:, i] = Float32.(σ)
        ds["sif_waterleaving"][:, i] = Float32.(sif_oci)
        ds["p"][:, i] = Float32.(0.5 .* (p_half[1:end-1] .+ p_half[2:end]))
        ds["T"][:, i] = Float32.(T)
        ds["q"][:, i] = Float32.(q)
        ds.attrib["n_completed"] = i
        sync(ds)
        if params.architecture isa GPU
            GC.gc(false)
        end

        dt = time() - t0
        rate = dt / (i - start_i + 1)
        eta = rate * (stop_i - i)
        aer_tag = ENABLE_AEROSOLS ? " col=$icol" : ""
        println("  $i / $N_SAMPLES  albedo=$(round(alb; digits=4))$(aer_tag)  " *
                "$(round(rate, digits=1)) s/sample  ETA $(round(eta / 60, digits=1)) min")
    end
    close(ds)
    close(sea)
    println(stop_i == N_SAMPLES ? "Done: $OUT_NC" : "Stopped after sample $stop_i (N_RUN=$N_RUN): $OUT_NC")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
