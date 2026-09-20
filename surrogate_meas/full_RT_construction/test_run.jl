using Pkg
Pkg.activate("/home/zhe2/FraLab/vSmartMOM.jl")
using vSmartMOM
using Plots
using vSmartMOM.SolarModel
using NCDatasets
using JLD2
using Interpolations
using Dates

# Same kernel as the SVD pipeline (`build_kernel_from_rsr_nc` + `KernelInstrument`).
include(joinpath(@__DIR__, "..", "..", "src", "tools", "Instrument.jl"))
include(joinpath(@__DIR__, "ocean_column_aerosols.jl"))

const PACE_RSR_NC = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/PACE_OCI_RSRs.nc"
const OCEAN_YAML = joinpath(@__DIR__, "..", "configs", "ocean_coxmunk_0912.yaml")
# 1) SIF
const SIF_LIB = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/SIF_singular_vector.jld2"
const SIF_LIBRARY_INDEX = 1          # column of SIF_shapes
const SIF_PEAK = 0.3                 # water-leaving radiance at SIF_λ, W m⁻² sr⁻¹ μm⁻¹
const SIF_λ = 678.0                  # nm; scale the shape here, not at its maximum
# 2) Aerosol columns
const ENABLE_AEROSOLS = false
const OCEAN_COLS_NC = joinpath(@__DIR__, "output_aerosol_profiles", "gchp_ocean_columns_n500.nc")
const COLUMN_INDEX = 100
# 3) Whitecap and wind speed
const INCLUDE_WHITECAPS = true
const WHITECAP_ALBEDO = 0.22
const WIND_SPEED = 5.0

# conversion factor from photon flux to radiance
h = 6.62607015e-34   # J⋅s
c = 299792458.0      # m/s
    # λ_nm = 1e7/ν, λ_m = λ_nm⋅1e-9 = 0.01/ν
    # E = hc/λ_m = 100⋅h⋅c⋅ν

# config + aerosols from n500 ocean column
params = parameters_from_yaml(OCEAN_YAML)
# Override Cox-Munk surface from script params (YAML is only a template).
surf0 = only(params.brdf)
surf0 isa vSmartMOM.CoreRT.CoxMunkSurface ||
    error("Expected CoxMunkSurface in YAML surface:; got $(typeof(surf0))")
FT = typeof(surf0.wind_speed)
params.brdf = [vSmartMOM.CoreRT.CoxMunkSurface{FT}(
    wind_speed = FT(WIND_SPEED),
    n_water = surf0.n_water,
    whitecap_albedo = FT(WHITECAP_ALBEDO),
    include_whitecaps = INCLUDE_WHITECAPS,
    shadowing = surf0.shadowing,
)]
println("Cox-Munk: wind=$(WIND_SPEED) m/s, whitecaps=$(INCLUDE_WHITECAPS), albedo=$(WHITECAP_ALBEDO)")
if ENABLE_AEROSOLS
    load_ocean_column_aerosols!(params, OCEAN_COLS_NC, COLUMN_INDEX; yaml_path=OCEAN_YAML)
else
    # Drop YAML scattering stub so model_from_parameters builds Rayleigh-only.
    params.scattering_params = nothing
    println("Aerosols disabled (ENABLE_AEROSOLS=false); ignoring YAML scattering: / ocean columns")
end
model = model_from_parameters(params)
ν     = params.spec_bands[1]
n_to_radiance = @. 100 * h * c * ν   # J per photon; ν in cm⁻¹

# solar beam
F_sol = SolarModel.default_solar_spectrum_at_earth(ν)[:, 2] 

n_stokes = params.polarization_type.n
F₀    = zeros(n_stokes, length(ν))   # 1 for Stokes_I(), 4 for Stokes_IQUV()
F₀[1, :] .= F_sol

# SIF source. SIF_shapes is a shape library (640–850 nm), not a physical radiance.
sif_file = jldopen(SIF_LIB)
λ_lib = Float64.(sif_file["SIF_wavelen"])
sif_shape = Float64.(sif_file["SIF_shapes"][:, SIF_LIBRARY_INDEX])
close(sif_file)
1 ≤ SIF_LIBRARY_INDEX || error("SIF_LIBRARY_INDEX must be ≥ 1")

λ_model = 1e7 ./ ν
dλ = λ_lib[2] - λ_lib[1]
itp_sif = CubicSplineInterpolation(
    range(λ_lib[1]; step=dλ, length=length(λ_lib)), sif_shape; extrapolation_bc=Line(),
)
I_wl = itp_sif.(λ_model)                 # W m⁻² sr⁻¹ μm⁻¹ after the scale below
I_wl[(λ_model .< λ_lib[1]) .| (λ_model .> λ_lib[end])] .= 0.0
sif_at_ref = itp_sif(SIF_λ)
sif_at_ref != 0 || error("SIF shape is zero at $(SIF_λ) nm")
I_wl .*= SIF_PEAK / sif_at_ref          # I_wl(SIF_λ) = SIF_PEAK, not the spectral maximum

# SurfaceSIF.SIF₀ is hemispheric irradiance; the solver divides by π to get
# Lambertian radiance. Keep the same per-μm photon unit as F_sol, not mW/cm⁻¹,
# or the two sources cannot be added.
SIF₀ = zeros(n_stokes, length(ν))        # unpolarized: I only (Q=U=V unused if n=1)
SIF₀[1, :] .= π .* I_wl ./ n_to_radiance
println("SIF library index $SIF_LIBRARY_INDEX, I_wl($(SIF_λ) nm) = $(SIF_PEAK) W m⁻² sr⁻¹ μm⁻¹, spectral max $(maximum(I_wl))")

# Cox-Munk is only the glint BRDF. Water-leaving SIF is isotropic and is added
# on the same m=0 j₀⁻ term the Lambertian path uses. Factor 2 = (1/π)·2π.
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

# add to model
model_source = SolarBeam(F₀ = F₀) + SurfaceSIF(SIF₀ = SIF₀)

# RT run
R, T, ieR, ieT, hdr, bhr_uw, bhr_dw = rt_run(model; sources = model_source)

# shape of each variables
println("Shape of R (TOA upwelling stokes field): $(size(R))")
    # dim: n_viewing_geometry x n_pol x n_spec_bands

# Ascending wavelength. ν and R are ascending wavenumber, so reverse both.
λ_hres = 1e7 ./ reverse(ν)
n_to_radiance_λ = reverse(n_to_radiance)
R_λ = reverse(R, dims=3) .* reshape(n_to_radiance_λ, 1, 1, :)   # W m⁻² sr⁻¹ μm⁻¹

# OCI RSR: clip negatives, interpolate onto λ_hres, row-normalize (Instrument.jl).
ds = NCDataset(PACE_RSR_NC)
wavlen = collect(Float64.(ds["wavelength"][:]))
band = collect(Float64.(ds["bands"][:]))
rsr_all = collect(Float64.(ds["RSR"][:, :]))   # NCDatasets: (wavelength, band)
close(ds)
λ_lo, λ_hi = extrema(λ_hres)
idx_w = findall(λ_lo .< wavlen .< λ_hi)
idx_b = findall(λ_lo .< band .< λ_hi)
isempty(idx_b) && error("No OCI bands inside the RT grid [$(λ_lo), $(λ_hi)] nm")
kernel = Instrument.KernelInstrument(
    band[idx_b], wavlen[idx_w], max.(rsr_all[idx_w, idx_b], 0.0),
    collect(Float64.(λ_hres)), collect(Float64.(reverse(ν))),
)
λ_oci = kernel.band
println("OCI kernel: $(size(kernel.RSR_out)) (bands × hi-res samples), bands $(extrema(λ_oci)) nm")

nview, npol, _ = size(R_λ)
R_oci = zeros(nview, npol, length(λ_oci))
for iv in 1:nview, ip in 1:npol
    R_oci[iv, ip, :] = kernel.RSR_out * vec(R_λ[iv, ip, :])
end

# Identify this run in the filename + NetCDF attributes.
const OUT_DIR = joinpath(@__DIR__, "output")
mkpath(OUT_DIR)
surf = only(params.brdf)
surf isa vSmartMOM.CoreRT.CoxMunkSurface ||
    error("Expected CoxMunkSurface in params.brdf; got $(typeof(surf))")
wind_speed = Float64(surf.wind_speed)
include_whitecaps = Bool(surf.include_whitecaps)
whitecap_albedo = Float64(surf.whitecap_albedo)
shadowing = Bool(surf.shadowing)
# Sanity: NC tags follow the script overrides, not the YAML defaults.
wind_speed ≈ WIND_SPEED || @warn "wind_speed mismatch" surf=wind_speed script=WIND_SPEED
include_whitecaps == INCLUDE_WHITECAPS || @warn "include_whitecaps mismatch"
whitecap_albedo ≈ WHITECAP_ALBEDO || @warn "whitecap_albedo mismatch"

_num_tag(x; digits=3) = replace(string(round(x; digits=digits)), "." => "p")
peak_tag = _num_tag(SIF_PEAK)
λ_tag = _num_tag(SIF_λ; digits=1)
ws_tag = _num_tag(wind_speed; digits=1)
wc_on_tag = include_whitecaps ? "true" : "false"
wc_alb_tag = _num_tag(whitecap_albedo; digits=2)
aer_suf = ENABLE_AEROSOLS ? "" : "_noaerosol"
out_nc = joinpath(OUT_DIR,
    "toa_col$(lpad(COLUMN_INDEX, 3, '0'))_siflib$(SIF_LIBRARY_INDEX)_peak$(peak_tag)_$(λ_tag)nm" *
    "_ws$(ws_tag)_wc$(wc_on_tag)_wc$(wc_alb_tag)$(aer_suf).nc")
isfile(out_nc) && rm(out_nc)

ds_out = NCDataset(out_nc, "c")
defDim(ds_out, "view", nview)
defDim(ds_out, "pol", npol)
defDim(ds_out, "band_hres", length(λ_hres))
defDim(ds_out, "band_oci", length(λ_oci))

defVar(ds_out, "wavelength_hres", Float64.(λ_hres), ("band_hres",);
       attrib=Dict("units" => "nm", "long_name" => "high-resolution wavelength (ascending)"))
defVar(ds_out, "wavelength_oci", Float64.(λ_oci), ("band_oci",);
       attrib=Dict("units" => "nm", "long_name" => "OCI band centers"))
defVar(ds_out, "radiance_hres", Float32.(R_λ), ("view", "pol", "band_hres");
       attrib=Dict("units" => "W m-2 sr-1 um-1",
                   "long_name" => "TOA upwelling Stokes radiance (hi-res)",
                   "pol_order" => "I,Q,U,V"))
defVar(ds_out, "radiance_oci", Float32.(R_oci), ("view", "pol", "band_oci");
       attrib=Dict("units" => "W m-2 sr-1 um-1",
                   "long_name" => "TOA upwelling Stokes radiance (OCI RSR)",
                   "pol_order" => "I,Q,U,V"))
defVar(ds_out, "sif_waterleaving_hres", Float32.(reverse(I_wl)), ("band_hres",);
       attrib=Dict("units" => "W m-2 sr-1 um-1",
                   "long_name" => "water-leaving SIF on hi-res grid (ascending λ)"))
defVar(ds_out, "sza", Float64(params.sza), ();
       attrib=Dict("units" => "degree"))
defVar(ds_out, "vza", Float64.(params.vza), ("view",);
       attrib=Dict("units" => "degree"))
defVar(ds_out, "vaz", Float64.(params.vaz), ("view",);
       attrib=Dict("units" => "degree"))

ds_out.attrib["title"] = "vSmartMOM Cox-Munk TOA spectra"
ds_out.attrib["created"] = string(Dates.now())
ds_out.attrib["ocean_yaml"] = OCEAN_YAML
ds_out.attrib["enable_aerosols"] = Int8(ENABLE_AEROSOLS)
ds_out.attrib["ocean_columns_nc"] = OCEAN_COLS_NC
ds_out.attrib["column_index"] = Int32(COLUMN_INDEX)
ds_out.attrib["sif_library"] = SIF_LIB
ds_out.attrib["sif_library_index"] = Int32(SIF_LIBRARY_INDEX)
ds_out.attrib["sif_peak"] = Float64(SIF_PEAK)
ds_out.attrib["sif_peak_units"] = "W m-2 sr-1 um-1"
ds_out.attrib["sif_lambda_nm"] = Float64(SIF_λ)
ds_out.attrib["pace_rsr_nc"] = PACE_RSR_NC
ds_out.attrib["wind_speed"] = wind_speed
ds_out.attrib["wind_speed_units"] = "m s-1"
ds_out.attrib["include_whitecaps"] = Int8(include_whitecaps)
ds_out.attrib["whitecap_albedo"] = whitecap_albedo
ds_out.attrib["shadowing"] = Int8(shadowing)
ds_out.attrib["run_id"] = basename(out_nc)
close(ds_out)
println("Wrote $out_nc")

# Plot
stokes = ("I", "Q", "U", "V")
p_iquv = plot(xlabel="Wavelength (nm)", ylabel="Radiance (W m⁻² sr⁻¹ μm⁻¹)",
              title="TOA upwelling Stokes field")
for ip in 1:4
    plot!(p_iquv, λ_hres, R_λ[1, ip, :], label="$(stokes[ip]) (hi-res)", color=ip, lw=1)
    plot!(p_iquv, λ_oci, R_oci[1, ip, :], label="$(stokes[ip]) (OCI)", color=ip,
          lw=1.5, ls=:dash, marker=:circle, ms=2)
end
display(p_iquv)

# sun-sensor geometry
p_geo = plot(xlabel="Wavelength (nm)", ylabel="Radiance (W m⁻² sr⁻¹ μm⁻¹)",
             title="TOA upwelling Stokes I")
for iv in 1:3
    plot!(p_geo, λ_hres, R_λ[iv, 1, :], label="geo $iv (hi-res)", color=iv, lw=1)
    plot!(p_geo, λ_oci, R_oci[iv, 1, :], label="geo $iv (OCI)", color=iv,
          lw=1.5, ls=:dash, marker=:circle, ms=2)
end
display(p_geo)

# does the ratio change with viewing geometry?
p_ratio = plot(xlabel="Wavelength (nm)", ylabel="Ratio",
               title="Ratio of TOA upwelling Stokes I")
plot!(p_ratio, λ_hres, R_λ[1, 1, :] ./ R_λ[2, 1, :], label="geo1/geo2 (hi-res)", color=1, lw=1)
plot!(p_ratio, λ_oci, R_oci[1, 1, :] ./ R_oci[2, 1, :], label="geo1/geo2 (OCI)", color=1,
      lw=1.5, ls=:dash, marker=:circle, ms=2)
plot!(p_ratio, λ_hres, R_λ[2, 1, :] ./ R_λ[3, 1, :], label="geo2/geo3 (hi-res)", color=2, lw=1)
plot!(p_ratio, λ_oci, R_oci[2, 1, :] ./ R_oci[3, 1, :], label="geo2/geo3 (OCI)", color=2,
      lw=1.5, ls=:dash, marker=:circle, ms=2)
display(p_ratio)

