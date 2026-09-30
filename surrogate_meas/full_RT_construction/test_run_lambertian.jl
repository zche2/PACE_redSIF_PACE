# Single TOA test spectrum with a Lambertian ocean surface (same atmosphere, aerosol
# column, SIF source, geometry and OCI convolution as `test_run.jl`).
#
# The surface comes straight from the YAML `surface:` entry (default
# `configs/ocean_lambertian_0929.yaml`, LambertianSurfaceScalar). If the YAML holds a
# Cox-Munk surface instead, the same wind/whitecap overrides and SIF hook as
# `test_run.jl` are applied, so both surfaces can be run through one script.
#
# Usage (from the repo root):
#   julia -t 8 surrogate_meas/full_RT_construction/test_run_lambertian.jl
# Environment overrides:
#   OCEAN_YAML=...          surface/atmosphere template
#   ENABLE_AEROSOLS=true    include the GCHP aerosol column (false → Rayleigh only)
#   COLUMN_INDEX=100        aerosol column in OCEAN_COLS_NC
#   OUT_DIR=...             output directory (default: ./output_new_Lambertian)

using Pkg
Pkg.activate("/home/zhe2/FraLab/vSmartMOM.jl")
using vSmartMOM
using vSmartMOM.SolarModel
using NCDatasets
using JLD2
using Interpolations
using Dates

# Same kernel as the SVD pipeline (`build_kernel_from_rsr_nc` + `KernelInstrument`).
include(joinpath(@__DIR__, "..", "..", "src", "tools", "Instrument.jl"))
include(joinpath(@__DIR__, "ocean_column_aerosols.jl"))

const REPO = "/home/zhe2/FraLab/PACE_redSIF_PACE"
const PACE_RSR_NC = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/PACE_OCI_RSRs.nc"
# configs/ is gitignored, so default to the main checkout's copy
const OCEAN_YAML = get(ENV, "OCEAN_YAML",
    joinpath(REPO, "surrogate_meas", "configs", "ocean_lambertian_0929.yaml"))
# 1) SIF
const SIF_LIB = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/SIF_singular_vector.jld2"
const SIF_LIBRARY_INDEX = 1          # column of SIF_shapes
const SIF_PEAK = 0.3                 # water-leaving radiance at SIF_λ, W m⁻² sr⁻¹ μm⁻¹
const SIF_λ = 678.0                  # nm; scale the shape here, not at its maximum
# 2) Aerosol columns
const ENABLE_AEROSOLS = parse(Bool, get(ENV, "ENABLE_AEROSOLS", "true"))
const OCEAN_COLS_NC = joinpath(REPO, "surrogate_meas", "full_RT_construction",
                               "output_aerosol_profiles", "gchp_ocean_columns_n500.nc")
const COLUMN_INDEX = parse(Int, get(ENV, "COLUMN_INDEX", "100"))
# 3) Cox-Munk overrides (used only if the YAML surface is Cox-Munk)
const INCLUDE_WHITECAPS = true
const WHITECAP_ALBEDO = 0.22
const WIND_SPEED = 5.0
const OUT_DIR = get(ENV, "OUT_DIR", joinpath(@__DIR__, "output_new_Lambertian"))

# conversion factor from photon flux to radiance: E = hc/λ = 100⋅h⋅c⋅ν (ν in cm⁻¹)
const h = 6.62607015e-34   # J⋅s
const c = 299792458.0      # m/s

params = parameters_from_yaml(OCEAN_YAML)
surf0 = only(params.brdf)
is_coxmunk = surf0 isa vSmartMOM.CoreRT.CoxMunkSurface
if is_coxmunk
    FTs = typeof(surf0.wind_speed)
    params.brdf = [vSmartMOM.CoreRT.CoxMunkSurface{FTs}(
        wind_speed = FTs(WIND_SPEED), n_water = surf0.n_water,
        whitecap_albedo = FTs(WHITECAP_ALBEDO), include_whitecaps = INCLUDE_WHITECAPS,
        shadowing = surf0.shadowing)]
    println("Surface: Cox-Munk, wind=$(WIND_SPEED) m/s, whitecaps=$(INCLUDE_WHITECAPS), albedo=$(WHITECAP_ALBEDO)")
elseif surf0 isa vSmartMOM.CoreRT.LambertianSurfaceScalar
    println("Surface: Lambertian, albedo=$(surf0.albedo)")
else
    println("Surface: $(typeof(surf0))")
end
if ENABLE_AEROSOLS
    load_ocean_column_aerosols!(params, OCEAN_COLS_NC, COLUMN_INDEX; yaml_path=OCEAN_YAML)
    println("Aerosols: GCHP column $(COLUMN_INDEX) from $(basename(OCEAN_COLS_NC))")
else
    params.scattering_params = nothing   # Rayleigh only
    println("Aerosols disabled")
end
model = model_from_parameters(params)
ν     = params.spec_bands[1]
n_to_radiance = @. 100 * h * c * ν   # J per photon

# solar beam
F_sol = SolarModel.default_solar_spectrum_at_earth(ν)[:, 2]
n_stokes = params.polarization_type.n
F₀ = zeros(n_stokes, length(ν))
F₀[1, :] .= F_sol

# SIF source (same scaling as test_run.jl)
sif_file = jldopen(SIF_LIB)
λ_lib = Float64.(sif_file["SIF_wavelen"])
sif_shape = Float64.(sif_file["SIF_shapes"][:, SIF_LIBRARY_INDEX])
close(sif_file)
λ_model = 1e7 ./ ν
dλ = λ_lib[2] - λ_lib[1]
itp_sif = CubicSplineInterpolation(
    range(λ_lib[1]; step=dλ, length=length(λ_lib)), sif_shape; extrapolation_bc=Line())
I_wl = itp_sif.(λ_model)
I_wl[(λ_model .< λ_lib[1]) .| (λ_model .> λ_lib[end])] .= 0.0
sif_at_ref = itp_sif(SIF_λ)
sif_at_ref != 0 || error("SIF shape is zero at $(SIF_λ) nm")
I_wl .*= SIF_PEAK / sif_at_ref
SIF₀ = zeros(n_stokes, length(ν))
SIF₀[1, :] .= π .* I_wl ./ n_to_radiance

# Cox-Munk only: inject isotropic water-leaving SIF on the m=0 term (as in test_run.jl).
# Lambertian surfaces already do this via vSmartMOM's built-in `inject_surface_SIF!`.
if is_coxmunk
    @eval function vSmartMOM.CoreRT.surface_source_contribute!(
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
end

R, = rt_run(model; sources = SolarBeam(F₀ = F₀) + SurfaceSIF(SIF₀ = SIF₀))

λ_hres = 1e7 ./ reverse(ν)
R_λ = reverse(R, dims=3) .* reshape(reverse(n_to_radiance), 1, 1, :)   # W m⁻² sr⁻¹ μm⁻¹

# OCI RSR convolution (Instrument.jl)
ds = NCDataset(PACE_RSR_NC)
wavlen = collect(Float64.(ds["wavelength"][:]))
band = collect(Float64.(ds["bands"][:]))
rsr_all = collect(Float64.(ds["RSR"][:, :]))
close(ds)
λ_lo, λ_hi = extrema(λ_hres)
idx_w = findall(λ_lo .< wavlen .< λ_hi)
idx_b = findall(λ_lo .< band .< λ_hi)
kernel = Instrument.KernelInstrument(
    band[idx_b], wavlen[idx_w], max.(rsr_all[idx_w, idx_b], 0.0),
    collect(Float64.(λ_hres)), collect(Float64.(reverse(ν))))
λ_oci = kernel.band
nview, npol, _ = size(R_λ)
R_oci = zeros(nview, npol, length(λ_oci))
for iv in 1:nview, ip in 1:npol
    R_oci[iv, ip, :] = kernel.RSR_out * vec(R_λ[iv, ip, :])
end

# output
mkpath(OUT_DIR)
_num_tag(x; digits=3) = replace(string(round(x; digits=digits)), "." => "p")
surf = only(params.brdf)
surf_tag = is_coxmunk ? "coxmunk_ws$(_num_tag(surf.wind_speed; digits=1))" :
           surf isa vSmartMOM.CoreRT.LambertianSurfaceScalar ? "lamb$(_num_tag(surf.albedo; digits=3))" :
           lowercase(string(nameof(typeof(surf))))
aer_tag = ENABLE_AEROSOLS ? "col$(lpad(COLUMN_INDEX, 3, '0'))" : "noaerosol"
out_nc = joinpath(OUT_DIR, "toa_$(surf_tag)_$(aer_tag)_siflib$(SIF_LIBRARY_INDEX)_peak$(_num_tag(SIF_PEAK))_$(_num_tag(SIF_λ; digits=1))nm.nc")
isfile(out_nc) && rm(out_nc)
NCDataset(out_nc, "c") do ds_out
    defDim(ds_out, "view", nview)
    defDim(ds_out, "pol", npol)
    defDim(ds_out, "band_hres", length(λ_hres))
    defDim(ds_out, "band_oci", length(λ_oci))
    defVar(ds_out, "wavelength_hres", Float64.(λ_hres), ("band_hres",); attrib=Dict("units" => "nm"))
    defVar(ds_out, "wavelength_oci", Float64.(λ_oci), ("band_oci",); attrib=Dict("units" => "nm"))
    defVar(ds_out, "radiance_hres", Float32.(R_λ), ("view", "pol", "band_hres");
           attrib=Dict("units" => "W m-2 sr-1 um-1", "long_name" => "TOA upwelling radiance (hi-res)"))
    defVar(ds_out, "radiance_oci", Float32.(R_oci), ("view", "pol", "band_oci");
           attrib=Dict("units" => "W m-2 sr-1 um-1", "long_name" => "TOA upwelling radiance (OCI RSR)"))
    defVar(ds_out, "solar_irradiance_hres", Float32.(reverse(F_sol .* n_to_radiance)), ("band_hres",);
           attrib=Dict("units" => "W m-2 um-1"))
    defVar(ds_out, "sif_waterleaving_hres", Float32.(reverse(I_wl)), ("band_hres",);
           attrib=Dict("units" => "W m-2 sr-1 um-1"))
    defVar(ds_out, "sza", Float64(params.sza), (); attrib=Dict("units" => "degree"))
    defVar(ds_out, "vza", Float64.(params.vza), ("view",); attrib=Dict("units" => "degree"))
    defVar(ds_out, "vaz", Float64.(params.vaz), ("view",); attrib=Dict("units" => "degree"))
    ds_out.attrib["title"] = "vSmartMOM TOA spectra ($(surf_tag))"
    ds_out.attrib["created"] = string(Dates.now())
    ds_out.attrib["ocean_yaml"] = OCEAN_YAML
    ds_out.attrib["surface"] = string(surf)
    ds_out.attrib["enable_aerosols"] = Int8(ENABLE_AEROSOLS)
    ds_out.attrib["ocean_columns_nc"] = OCEAN_COLS_NC
    ds_out.attrib["column_index"] = Int32(COLUMN_INDEX)
    ds_out.attrib["sif_peak"] = Float64(SIF_PEAK)
    ds_out.attrib["sif_lambda_nm"] = Float64(SIF_λ)
end
println("Wrote $out_nc")
