# TOA spectra for every combination of Rayleigh scale × aerosol scale × SIF intensity × (vza, vaz).
#
# The scene is OCEAN_YAML as written: sza, the paired (vza[i], vaz[i]) views, T/p/q profile, Lambertian
# albedo, gases, spectral grid, streams and truncation. Scan axes, from YAML keys (or env overrides):
#   rayleigh_scaling: [...]   RAYLEIGH_SCALES   factor on the Rayleigh optical depth of every layer
#   aerosol_scaling:  [...]   AEROSOL_SCALES    factor on the GCHPIO column-AER_COLUMN optical depth
#                                               (ω̃ and phase function unchanged); AER_COLUMN=0 → no aerosol
#   SIF_intensity:    [...]   SIF_INTENSITIES   water-leaving SIF at 678 nm (W m⁻² sr⁻¹ µm⁻¹),
#                                               library shape SIF_INDEX
#
# RT is linear in the SurfaceSIF source, so each (Rayleigh, aerosol) pair needs two rt_run calls:
# no SIF (R₀) and SIF with I(678) = 1 (R₁); every intensity is R₀ + SIF₆₇₈·(R₁ − R₀). One explicit
# run at the largest intensity checks this. One rt_run returns all views.
# SUNLIGHT_ON=false drops the solar beam (SIF-only TOA radiance).
#
#   julia --project=surrogate_meas/full_RT_construction/gchpio_env \
#         surrogate_meas/full_RT_construction/generate_rayleigh_scaling_spectra.jl
#   ARCH=CPU AER_COLUMN=77 SCAN_OUT_NC=/path/out.nc julia --project=... ...

include(joinpath(@__DIR__, "generate_test_spec_lambertian_gchpio.jl"))

const SCAN_YAML = get(ENV, "SCAN_YAML",
    joinpath(REPO_RT, "..", "configs", "ocean_lambertian_1010_wavnum0p05_NoRayl.yaml"))
const SCAN_SIF_INDEX = parse(Int, get(ENV, "SIF_INDEX", "1"))
const AER_COLUMN = parse(Int, get(ENV, "AER_COLUMN", "100"))         # 0 = no aerosol
const SUNLIGHT_ON = parse(Bool, get(ENV, "SUNLIGHT_ON", "true"))     # false = solar beam off
const SCAN_OUT_NC = get(ENV, "SCAN_OUT_NC", joinpath(REPO_RT, "output_rayleigh_scaling",
    "rt_toa_rayleigh_scaling_" * splitext(basename(SCAN_YAML))[1] *
    "_" * (SUNLIGHT_ON ? "sunlight_on" : "sunlight_off") *
    "_AerColumn_" * string(AER_COLUMN) * "_grid.nc"))

"List from ENV[env] (comma-separated), else the first YAML line `<key>: [...]` for any of `keys`."
function yaml_list(path, env, keys...)
    haskey(ENV, env) && return parse.(Float64, strip.(split(ENV[env], ",")))
    pat = Regex("^\\s*(?:" * join(keys, "|") * ")\\s*:\\s*\\[(.*)\\]")
    for line in eachline(path)
        m = match(pat, line)
        m === nothing || return parse.(Float64, strip.(split(m[1], ",")))
    end
    error("No `$(first(keys)): [...]` in $path; set $env=...")
end

"Column optics with the aerosol optical depth scaled by `a` (ω̃ and Greek coefficients unchanged)."
scale_column(col, a) = merge(col, (τ = a .* col.τ,))

function scan_main()
    ray_scales = yaml_list(SCAN_YAML, "RAYLEIGH_SCALES", "rayleigh_scaling", "raleigh_scaling")
    sif_vals = yaml_list(SCAN_YAML, "SIF_INTENSITIES", "SIF_intensity", "sif_intensity")
    aer_scales = yaml_list(SCAN_YAML, "AEROSOL_SCALES", "aerosol_scaling")
    if AER_COLUMN == 0
        @warn "AER_COLUMN=0: no aerosol, aerosol_scaling $(aer_scales) collapsed to [0]"
        aer_scales = [0.0]
    end
    params = parameters_from_yaml(SCAN_YAML)
    if ARCH == "CPU"
        params.architecture = CPU()
    elseif ARCH == "GPU"
        params.architecture = GPU()
    end
    vza, vaz = Float64.(params.vza), Float64.(params.vaz)
    length(vza) == length(vaz) || error("vza has $(length(vza)) angles but vaz has $(length(vaz)); they are paired")
    println("Rayleigh scales $ray_scales, aerosol scales $aer_scales (column $AER_COLUMN), " *
            "SIF $sif_vals, views $(collect(zip(vza, vaz))), SUNLIGHT_ON=$SUNLIGHT_ON")
    surf0 = only(params.brdf)
    surf0 isa vSmartMOM.CoreRT.LambertianSurfaceScalar ||
        error("Expected LambertianSurfaceScalar in $SCAN_YAML surface:; got $(typeof(surf0))")
    trunc = aerosol_truncation(params)
    params.scattering_params = nothing        # YAML aerosol is a τ_ref = 0 placeholder
    model = model_from_parameters(params)

    ν = params.spec_bands[1]
    n_stokes = params.polarization_type.n
    n_to_radiance = @. 100 * h * c_light * ν
    F₀ = zeros(n_stokes, length(ν))
    F₀[1, :] .= SolarModel.default_solar_spectrum_at_earth(ν)[:, 2]
    λ_hres = 1e7 ./ reverse(ν)
    kernel = oci_kernel(λ_hres, reverse(ν))
    λ_oci = collect(Float64.(kernel.band))
    c1, c2 = snr_coeffs(λ_oci)

    sif_file = jldopen(SIF_LIB)
    sif_unit = unit_sif_shapes(ν, Float64.(sif_file["SIF_wavelen"]), Float64.(sif_file["SIF_shapes"]))
    close(sif_file)
    I_unit = sif_unit[:, SCAN_SIF_INDEX]                             # I(678) = 1, wavenumber order

    p_half = Float64.(model.atmosphere.profile.p_half)
    n_layer = size(model.optics.τ_rayl[1], 2)
    aod_oci = fill(NaN, length(λ_oci))
    col = nothing
    if AER_COLUMN > 0
        length(p_half) == n_layer + 1 || error("p_half has $(length(p_half)) levels for $n_layer layers")
        cols = column_table()
        gf = GCHPFile(cols.gchp; aerosol_scheme=:auto)
        λq = collect(range(1e4 / maximum(ν), 1e4 / minimum(ν); length=N_QUERY))
        col = compute_column_optics(gf, cols, AER_COLUMN, λq)
        aod_oci = column_aod(col, λ_oci)
        println("GCHPIO column $AER_COLUMN: AOD $(round(aod_oci[1]; digits=4)) → $(round(aod_oci[end]; digits=4))")
    end

    # TOA radiance (ascending λ × view) for model `m` with water-leaving SIF amplitude `amp`
    function toa(m, amp)
        SIF₀ = zeros(n_stokes, length(ν))
        SIF₀[1, :] .= amp .* π .* I_unit ./ n_to_radiance
        sources = SUNLIGHT_ON ? SolarBeam(F₀=F₀) + SurfaceSIF(SIF₀=SIF₀) : SurfaceSIF(SIF₀=SIF₀)
        R, = rt_run(m; sources=sources)
        return reverse(permutedims(Array(R)[:, 1, :] .* n_to_radiance'); dims=1)
    end

    τ0 = deepcopy(model.optics.τ_rayl)
    τ_col_hres = reverse(vec(sum(τ0[1]; dims=2)))                 # ascending λ
    n_hres, n_band, n_v = length(ν), length(λ_oci), length(vza)
    n_r, n_a, n_f = length(ray_scales), length(aer_scales), length(sif_vals)
    R0 = zeros(n_hres, n_v, n_r, n_a)        # no SIF
    R1 = zeros(n_hres, n_v, n_r, n_a)        # SIF with I(678) = 1
    for (ia, a) in enumerate(aer_scales), (ir, s) in enumerate(ray_scales)
        for ib in eachindex(τ0)
            model.optics.τ_rayl[ib] .= s .* τ0[ib]
        end
        m = col === nothing ? model : inject_aerosols(model, params, scale_column(col, a), p_half, trunc)
        t = @elapsed begin
            R0[:, :, ir, ia] .= toa(m, 0.0)
            R1[:, :, ir, ia] .= toa(m, 1.0)
        end
        println("  aerosol ×$a  Rayleigh ×$s  $(round(t; digits=1)) s")
        if ia == 1 && ir == 1
            amp = maximum(sif_vals)
            lin = R0[:, :, 1, 1] .+ amp .* (R1[:, :, 1, 1] .- R0[:, :, 1, 1])
            err = maximum(abs.(toa(m, amp) .- lin)) / maximum(abs.(R1[:, :, 1, 1] .- R0[:, :, 1, 1]))
            println("  SIF linearity check at SIF₆₇₈=$amp: max |explicit − linear| / max SIF signal = $(err)")
            err < 1e-6 || @warn "SIF response is not linear to 1e-6 (got $err)"
        end
    end
    for ib in eachindex(τ0)
        model.optics.τ_rayl[ib] .= τ0[ib]
    end

    conv(x) = reshape(kernel.RSR_out * reshape(x, n_hres, :), n_band, size(x)[2:end]...)
    O0, O1 = conv(R0), conv(R1)
    dO = O1 .- O0
    I_oci = Array{Float64}(undef, n_band, n_f, n_v, n_r, n_a)
    R_hres = Array{Float32}(undef, n_hres, n_f, n_v, n_r, n_a)
    for (k, amp) in enumerate(sif_vals)
        I_oci[:, k, :, :, :] .= O0 .+ amp .* dO
        R_hres[:, k, :, :, :] .= R0 .+ amp .* (R1 .- R0)
    end
    sif_toa = I_oci .- reshape(O0, n_band, 1, n_v, n_r, n_a)
    σ = sqrt.(c1 .+ c2 .* max.(I_oci, 0.0))
    sif_wl = (kernel.RSR_out * reverse(I_unit)) .* sif_vals'
    aod = aod_oci .* aer_scales'

    mkpath(dirname(SCAN_OUT_NC))
    isfile(SCAN_OUT_NC) && rm(SCAN_OUT_NC)
    NCDataset(SCAN_OUT_NC, "c") do ds
        defDim(ds, "band", n_band); defDim(ds, "spectral", n_hres); defDim(ds, "sif", n_f)
        defDim(ds, "view", n_v); defDim(ds, "rayleigh", n_r); defDim(ds, "aerosol", n_a)
        unit = "W m-2 sr-1 um-1"
        put(name, T, dims, data, attrib) = defVar(ds, name, T.(data), dims; attrib=attrib)
        put("rayleigh_scale", Float64, ("rayleigh",), ray_scales,
            Dict("long_name" => "factor applied to the Rayleigh optical depth of every layer"))
        put("aerosol_scale", Float64, ("aerosol",), aer_scales,
            Dict("long_name" => "factor applied to the GCHPIO aerosol optical depth (0: no aerosol)"))
        put("sif_678", Float64, ("sif",), sif_vals,
            Dict("units" => unit, "long_name" => "water-leaving SIF at $(SIF_λ) nm (high-res)"))
        put("vza", Float64, ("view",), vza, Dict("units" => "degree"))
        put("vaz", Float64, ("view",), vaz, Dict("units" => "degree", "long_name" => "relative azimuth, vSmartMOM convention"))
        put("wavelength", Float64, ("band",), λ_oci, Dict("units" => "nm", "long_name" => "OCI band center"))
        put("wavelength_hres", Float64, ("spectral",), λ_hres, Dict("units" => "nm"))
        put("radiance_hres", Float32, ("spectral", "sif", "view", "rayleigh", "aerosol"), R_hres,
            Dict("units" => unit, "long_name" => "TOA Stokes I with SIF, RT grid"))
        put("radiance_hres_nosif", Float32, ("spectral", "view", "rayleigh", "aerosol"), R0,
            Dict("units" => unit, "long_name" => "TOA Stokes I without SIF, RT grid"))
        put("radiance_clean", Float32, ("band", "sif", "view", "rayleigh", "aerosol"), I_oci,
            Dict("units" => unit, "long_name" => "TOA Stokes I with SIF, OCI-convolved, no noise"))
        put("radiance_nosif", Float32, ("band", "view", "rayleigh", "aerosol"), O0,
            Dict("units" => unit, "long_name" => "TOA Stokes I without SIF, OCI-convolved"))
        put("sif_toa", Float32, ("band", "sif", "view", "rayleigh", "aerosol"), sif_toa,
            Dict("units" => unit, "long_name" => "TOA SIF signal: radiance_clean - radiance_nosif"))
        put("sigma_noise", Float32, ("band", "sif", "view", "rayleigh", "aerosol"), σ,
            Dict("units" => unit, "long_name" => "OCI SNR σ for radiance_clean, σ² = c1 + c2·R"))
        put("sif_waterleaving", Float32, ("band", "sif"), sif_wl,
            Dict("units" => unit, "long_name" => "water-leaving SIF on OCI bands"))
        put("rayleigh_od_column", Float32, ("band",), kernel.RSR_out * τ_col_hres,
            Dict("units" => "1", "long_name" => "column Rayleigh optical depth at scale 1, OCI-convolved"))
        put("aod", Float32, ("band", "aerosol"), aod,
            Dict("units" => "1", "long_name" => "GCHPIO column AOD × aerosol_scale (NaN: no aerosol)"))
        ds.attrib["ocean_yaml"] = SCAN_YAML
        ds.attrib["sza"] = Float64(params.sza)
        ds.attrib["surface_albedo"] = Float64(surf0.albedo)
        ds.attrib["sif_library_index"] = SCAN_SIF_INDEX
        ds.attrib["sif_lambda_nm"] = SIF_λ
        ds.attrib["sif_method"] = "linear in SurfaceSIF: R0 + sif_678 * (R(SIF678=1) - R0)"
        ds.attrib["sunlight_on"] = Int8(SUNLIGHT_ON)
        ds.attrib["aer_column"] = AER_COLUMN
        ds.attrib["aerosol_truncation"] = string(typeof(trunc))
        ds.attrib["architecture"] = string(typeof(params.architecture))
    end
    println("Done: $SCAN_OUT_NC")
end

scan_main()
