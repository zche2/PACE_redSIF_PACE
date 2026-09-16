# Generate an ensemble of TOA Stokes-I spectra with the test_run.jl RT path.
#
# Sampling per sample:
#   - p, T, q from a random MERRA-2 profile, downsampled every PROFILE_STRIDE
#     half-levels (~24 layers for stride=3) before RT
#   - SIF shape from SIF_shapes, magnitude at 678 nm in SIF_STRENGTH
#   - geometry: SZA / VZA / VAZ
#   - Cox-Munk: wind_speed, whitecap_albedo, include_whitecaps
#   - aerosols (if ENABLE_AEROSOLS): random GEOS-Chem ocean column
#   - white noise from the PACE OCI SNR model after OCI convolution
#
# Resume: if OUT_NC already exists with the same seed and N_SAMPLES, unfinished
# samples are filled in. One 1 cm⁻¹ Cox-Munk run is slow; this job checkpoints
# after every sample.
#
#   julia surrogate_meas/full_RT_construction/generate_test_spec.jl
#   N_SAMPLES=2 julia surrogate_meas/full_RT_construction/generate_test_spec.jl
#   ENABLE_AEROSOLS=false N_SAMPLES=2 julia ...

using Pkg
Pkg.activate("/home/zhe2/FraLab/vSmartMOM.jl")
using vSmartMOM
using vSmartMOM.SolarModel
using NCDatasets
using JLD2
using Interpolations
using Random
using Statistics
using DelimitedFiles

include(joinpath(@__DIR__, "..", "..", "src", "tools", "Instrument.jl"))
include(joinpath(@__DIR__, "ocean_column_aerosols.jl"))

const OCEAN_YAML = joinpath(@__DIR__, "..", "configs", "ocean_coxmunk_0912.yaml")
const PACE_RSR_NC = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/PACE_OCI_RSRs.nc"
const SIF_LIB = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/SIF_singular_vector.jld2"
const OCEAN_COLS_NC = joinpath(@__DIR__, "output_aerosol_profiles", "geoschem_ocean_columns_n500.nc")
const TRANS_NC = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc"
const SNR_FILE = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt"
const OUT_NC = get(ENV, "OUT_NC", joinpath(@__DIR__, "output", "rt_toa_ensemble.nc"))

const N_SAMPLES = parse(Int, get(ENV, "N_SAMPLES", "1000"))
const ENABLE_AEROSOLS = parse(Bool, get(ENV, "ENABLE_AEROSOLS", "true"))
const SEED = parse(Int, get(ENV, "ENSEMBLE_SEED", "20260913"))
const SIF_λ = 678.0
const SIF_STRENGTH = (0.0, 0.5)          # W m⁻² sr⁻¹ μm⁻¹ at SIF_λ
const SZA_RANGE = (5.0, 70.0)            # deg
const VZA_RANGE = (0.0, 60.0)            # deg; OCI-like swath
const VAZ_RANGE = (0.0, 180.0)           # deg; relative azimuth (glint vs dark)
const WIND_SPEED_RANGE = (0.0, 10.0)     # m/s
const WHITECAP_ALBEDO_RANGE = (0.1, 0.5)  # 0-1
const INCLUDE_WHITECAPS_RANGE = (0, 1)    # 0 or 1
# MERRA has 72 layers / 73 half-levels; stride 3 → ~24 RT layers.
const PROFILE_STRIDE = parse(Int, get(ENV, "PROFILE_STRIDE", "3"))

const h = 6.62607015e-34
const c_light = 299792458.0

# Isotropic water-leaving SIF on top of the Cox-Munk glint BRDF. Factor 2 = (1/π)·2π.
function vSmartMOM.CoreRT.surface_source_contribute!(
        prep::vSmartMOM.CoreRT.PreparedSurfaceSIF,
        ::vSmartMOM.CoreRT.CoxMunkSurface,
        surface_added_layer, m::Integer, pol_type, architecture)
    m == 0 || return nothing
    iszero(prep.SIF₀) && return nothing
    FT = eltype(surface_added_layer.j₀⁻)
    Nquad = size(surface_added_layer.j₀⁻, 1) ÷ pol_type.n
    surface_added_layer.j₀⁻[:, 1, :] .+= FT(2) .* array_type(architecture)(repeat(FT.(prep.SIF₀), Nquad))
    return nothing
end

"New observation geometry + Cox-Munk surface sharing the updated atmosphere/optics."
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

function make_coxmunk(surf0, wind_speed, whitecap_albedo, include_whitecaps)
    FT = typeof(surf0.wind_speed)
    return vSmartMOM.CoreRT.CoxMunkSurface{FT}(
        wind_speed = FT(wind_speed),
        n_water = surf0.n_water,
        whitecap_albedo = FT(whitecap_albedo),
        include_whitecaps = Bool(include_whitecaps),
        shadowing = surf0.shadowing,
    )
end

"""Push GEOS-Chem ocean column `icol` into `ctx` aerosol optics (fixed species count)."""
function apply_ocean_column!(ctx, params, icol)
    load_ocean_column_aerosols!(params, OCEAN_COLS_NC, icol; yaml_path=OCEAN_YAML, aod_min=0.0)
    rt_list = params.scattering_params.rt_aerosols
    length(rt_list) == ctx.n_aerosols || error(
        "Column $icol has $(length(rt_list)) aerosols but BatchContext expects $(ctx.n_aerosols)")
    for i in 1:ctx.n_aerosols
        rta = rt_list[i]
        vSmartMOM.CoreRT.update_aerosol_loading!(
            ctx, i; τ_ref=rta.τ_ref, profile_dist=rta.profile)
        vSmartMOM.CoreRT.update_aerosol_microphysics!(
            ctx, i, rta.aerosol; τ_ref=rta.τ_ref)
    end
    return nothing
end

function merra_profile(ds, i)
    T = Float64.(ds["temperature"][i, :])
    q = Float64.(ds["q"][i, :])
    ps = Float64(ds["pressure"][i])
    ak = Float64.(ds.attrib["ak"])
    bk = Float64.(ds.attrib["bk"])
    p_half = (ak .+ bk .* (ps * 100.0)) ./ 100.0   # hPa
    return (T=T, q=q, p_half=p_half, ps=ps)
end

"""Keep every `stride`-th MERRA half-level, always including the surface (~24 layers for stride=3)."""
function downsample_half_levels(p_half; stride::Int=PROFILE_STRIDE)
    stride >= 1 || error("PROFILE_STRIDE must be ≥ 1")
    n = length(p_half)
    idx = collect(1:stride:n)
    idx[end] == n || push!(idx, n)
    return p_half[idx]
end

"Map a MERRA-2 profile onto the model's half-level count. Surface pressure is the MERRA value."
function resample_profile(src, p_half_template)
    p_src = src.p_half
    p_src[1] < p_src[end] || error("MERRA p_half must increase downward")
    logp_t = log.(p_half_template)
    f = (logp_t .- logp_t[1]) ./ (logp_t[end] - logp_t[1])
    p_half = exp.(log(p_src[1]) .+ f .* (log(p_src[end]) - log(p_src[1])))
    p_full_src = 0.5 .* (p_src[1:end-1] .+ p_src[2:end])
    p_full = 0.5 .* (p_half[1:end-1] .+ p_half[2:end])
    itpT = LinearInterpolation(log.(p_full_src), src.T; extrapolation_bc=Flat())
    itpQ = LinearInterpolation(log.(p_full_src), src.q; extrapolation_bc=Flat())
    T = itpT.(log.(p_full))
    q = max.(itpQ.(log.(p_full)), 0.0)
    return p_half, T, q
end

"""Replace YAML p/T/q with a MERRA half-level grid downsampled by PROFILE_STRIDE."""
function apply_downsampled_atmosphere!(params, merra_ds)
    FT = params.float_type
    src0 = merra_profile(merra_ds, 1)
    p_half_ds = downsample_half_levels(src0.p_half)
    p_half, T, q = resample_profile(src0, p_half_ds)
    params.p = FT.(p_half)
    params.T = FT.(T)
    params.q = FT.(q)
    params.profile_reduction_n = -1   # already at the target layer count
    n_layer = length(T)
    println("RT atmosphere: $(length(src0.T)) MERRA layers → $n_layer layers " *
            "(every $PROFILE_STRIDE half-levels)")
    return p_half
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
    itp1 = LinearInterpolation(λ[p], c1[p]; extrapolation_bc=Flat())
    itp2 = LinearInterpolation(λ[p], c2[p]; extrapolation_bc=Flat())
    return itp1.(λ_oci), itp2.(λ_oci)
end

"Library shapes as water-leaving radiance with I(678 nm) = 1, in wavenumber order."
function unit_sif_shapes(ν, λ_lib, shapes)
    dλ = λ_lib[2] - λ_lib[1]
    knots = range(λ_lib[1]; step=dλ, length=length(λ_lib))
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

function create_output(path, λ_oci, p_full, design, n_layer)
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
        "long_name" => "1-based index in $TRANS_NC"))
    defVar(ds, "sza", Float32, ("sample",); attrib=Dict("units" => "degree"))
    defVar(ds, "vza", Float32, ("sample",); attrib=Dict("units" => "degree"))
    defVar(ds, "vaz", Float32, ("sample",); attrib=Dict(
        "units" => "degree", "long_name" => "relative azimuth, vSmartMOM convention"))
    defVar(ds, "wind_speed", Float32, ("sample",); attrib=Dict(
        "units" => "m s-1", "long_name" => "Cox-Munk 10-m wind speed"))
    defVar(ds, "whitecap_albedo", Float32, ("sample",); attrib=Dict(
        "units" => "1", "long_name" => "Cox-Munk whitecap Lambertian albedo"))
    defVar(ds, "include_whitecaps", Int8, ("sample",); attrib=Dict(
        "long_name" => "1 = whitecaps on, 0 = off"))
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
    ds.attrib["profile_stride"] = PROFILE_STRIDE
    ds.attrib["radiance"] = "TOA Stokes I from Cox-Munk + SurfaceSIF, OCI-convolved"
    ds.attrib["pressure_template_hpa"] = join(string.(p_full), ",")
    return ds
end

function main()
    println("Building RT model from $OCEAN_YAML  (aerosols=$(ENABLE_AEROSOLS))")
    params = parameters_from_yaml(OCEAN_YAML)
    surf0 = only(params.brdf)
    surf0 isa vSmartMOM.CoreRT.CoxMunkSurface ||
        error("Expected CoxMunkSurface in OCEAN_YAML surface:; got $(typeof(surf0))")

    merra = NCDataset(TRANS_NC)
    n_profiles = merra.dim["profile"]
    # Replace YAML ~33-layer grid with MERRA half-levels downsampled every PROFILE_STRIDE.
    p_template = apply_downsampled_atmosphere!(params, merra)
    n_layer = length(p_template) - 1

    n_cols = 0
    if ENABLE_AEROSOLS
        isfile(OCEAN_COLS_NC) || error("Missing ocean columns file: $OCEAN_COLS_NC")
        ds_cols = NCDataset(OCEAN_COLS_NC)
        n_cols = Int(ds_cols.dim["sample"])
        close(ds_cols)
        # Seed aerosols so BatchContext allocates a fixed species count (aod_min=0).
        load_ocean_column_aerosols!(params, OCEAN_COLS_NC, 1; yaml_path=OCEAN_YAML, aod_min=0.0)
    else
        params.scattering_params = nothing
        println("Aerosols disabled; Rayleigh-only optics")
    end

    ctx = vSmartMOM.CoreRT.BatchContext(params)
    ν = params.spec_bands[1]
    n_to_radiance = @. 100 * h * c_light * ν
    F_sol = SolarModel.default_solar_spectrum_at_earth(ν)[:, 2]
    F₀ = zeros(4, length(ν))
    F₀[1, :] .= F_sol

    λ_hres = 1e7 ./ reverse(ν)
    kernel = oci_kernel(λ_hres, reverse(ν))
    λ_oci = collect(Float64.(kernel.band))
    c1, c2 = snr_coeffs(λ_oci)
    println("OCI bands: $(length(λ_oci)) in $(extrema(λ_oci)) nm")

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
        if stored_n == N_SAMPLES && stored_seed == SEED && done < N_SAMPLES
            println("Resuming $OUT_NC at sample $(done + 1) / $N_SAMPLES")
            existing, done + 1
        elseif done >= N_SAMPLES && stored_n == N_SAMPLES && stored_seed == SEED
            println("Already complete: $OUT_NC")
            close(existing)
            close(merra)
            return
        else
            close(existing)
            error("Existing $OUT_NC does not match N_SAMPLES=$N_SAMPLES seed=$SEED. Remove it or set OUT_NC.")
        end
    else
        Random.seed!(SEED)
        design = draw_design(N_SAMPLES, n_profiles, size(sif_unit, 2), max(n_cols, 1))
        out = create_output(OUT_NC, λ_oci, 0.5 .* (p_template[1:end-1] .+ p_template[2:end]), design, n_layer)
        println("Writing $OUT_NC  ($N_SAMPLES samples)")
        out, 1
    end

    last_col = 0
    t0 = time()
    for i in start_i:N_SAMPLES
        ip = Int(ds["profile_index"][i])
        isif = Int(ds["sif_library_index"][i])
        strength = Float64(ds["sif_678"][i])
        sza = Float64(ds["sza"][i])
        vza = Float64(ds["vza"][i])
        vaz = Float64(ds["vaz"][i])
        wind = Float64(ds["wind_speed"][i])
        wc_alb = Float64(ds["whitecap_albedo"][i])
        wc_on = Bool(Int(ds["include_whitecaps"][i]))
        icol = Int(ds["column_index"][i])

        p_half, T, q = resample_profile(merra_profile(merra, ip), p_template)
        vSmartMOM.CoreRT.update_model!(ctx; T=T, p_half=p_half, q=q)

        if ENABLE_AEROSOLS && icol != last_col
            apply_ocean_column!(ctx, params, icol)
            last_col = icol
        end

        surf = make_coxmunk(surf0, wind, wc_alb, wc_on)
        scene = with_geometry(ctx.model, params, sza, vza, vaz, surf)

        I_wl = sif_unit[:, isif] .* strength
        SIF₀ = zeros(4, length(ν))
        SIF₀[1, :] .= π .* I_wl ./ n_to_radiance
        sources = SolarBeam(F₀=F₀) + SurfaceSIF(SIF₀=SIF₀)
        R, = rt_run(scene; sources=sources)

        R_I = reverse(R[1, 1, :] .* n_to_radiance)
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

        dt = time() - t0
        rate = dt / (i - start_i + 1)
        eta = rate * (N_SAMPLES - i)
        aer_tag = ENABLE_AEROSOLS ? " col=$icol" : ""
        println("  $i / $N_SAMPLES  ws=$(round(wind; digits=1)) wc=$(wc_on)$(aer_tag)  " *
                "$(round(rate, digits=1)) s/sample  ETA $(round(eta / 60, digits=1)) min")
    end
    close(ds)
    close(merra)
    println("Done: $OUT_NC")
end

main()
