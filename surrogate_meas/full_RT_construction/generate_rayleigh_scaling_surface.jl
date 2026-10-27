# Same scan as generate_rayleigh_scaling_spectra.jl (Rayleigh scale × aerosol scale × SIF × views),
# but keeps all rt_run outputs and writes the near-surface fields as well:
#
#   R       TOA upwelling radiance in each view direction            → L_toa
#   T       downwelling radiance at the bottom of the atmosphere      → L_down_boa  (view directions)
#   hdr     upwelling radiance just above the surface (reflected      → Lw          (view directions)
#           downwelling field + surface source, i.e. incl. SIF)
#   bhr_dw  downwelling irradiance just above the surface             → Ed
#           (diffuse + direct beam)
#   bhr_uw  upwelling irradiance just above the surface               → Eu
#   Rrs = Lw / Ed   (sr⁻¹); with SIF off and a Lambertian surface it should equal albedo / π
#
# All fields are converted with the same n_to_radiance factor as R. Every SIF intensity follows from
# two RT runs per (Rayleigh, aerosol) pair, as in the original (RT is linear in SurfaceSIF).
# Writes a new NetCDF (default name ends in _grid_surface.nc); existing files are not touched.
#
#   julia --project=surrogate_meas/full_RT_construction/gchpio_env \
#         surrogate_meas/full_RT_construction/generate_rayleigh_scaling_surface.jl
#   ARCH=CPU AER_COLUMN=77 SCAN_OUT_NC=/path/out.nc julia --project=... ...

include("/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/full_RT_construction/generate_test_spec_lambertian_gchpio.jl")

const SCAN_YAML = get(ENV, "SCAN_YAML",
    joinpath(REPO_RT, "..", "configs", "ocean_lambertian_1010_wavnum0p05_NoRayl.yaml"))
const SCAN_SIF_INDEX = parse(Int, get(ENV, "SIF_INDEX", "1"))
const AER_COLUMN = parse(Int, get(ENV, "AER_COLUMN", "100"))         # 0 = no aerosol
const SUNLIGHT_ON = parse(Bool, get(ENV, "SUNLIGHT_ON", "true"))     # false = solar beam off
const SCAN_OUT_NC = get(ENV, "SCAN_OUT_NC", joinpath(REPO_RT, "output_rayleigh_scaling",
    "rt_toa_rayleigh_scaling_" * splitext(basename(SCAN_YAML))[1] *
    "_" * (SUNLIGHT_ON ? "sunlight_on" : "sunlight_off") *
    "_AerColumn_" * string(AER_COLUMN) * "_grid_surface.nc"))

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
    isfile(SCAN_OUT_NC) && error("$SCAN_OUT_NC exists; this script never overwrites, set SCAN_OUT_NC")
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

    # All rt_run outputs for model `m` with water-leaving SIF amplitude `amp`, ascending λ,
    # in radiance / irradiance units (× n_to_radiance). View fields: (spectral, view).
    function fields(m, amp)
        SIF₀ = zeros(n_stokes, length(ν))
        SIF₀[1, :] .= amp .* π .* I_unit ./ n_to_radiance
        sources = SUNLIGHT_ON ? SolarBeam(F₀=F₀) + SurfaceSIF(SIF₀=SIF₀) : SurfaceSIF(SIF₀=SIF₀)
        R, T, ieR, ieT, hdr, bhr_uw, bhr_dw = rt_run(m; sources=sources)
        view_field(X) = reverse(permutedims(Array(X)[:, 1, :] .* n_to_radiance'); dims=1)
        spec_field(x) = reverse(Array(x) .* n_to_radiance)
        return (L_toa = view_field(R), L_down = view_field(T), Lw = view_field(hdr),
                Ed = spec_field(bhr_dw), Eu = spec_field(bhr_uw))
    end

    τ0 = deepcopy(model.optics.τ_rayl)
    n_hres, n_band, n_v = length(ν), length(λ_oci), length(vza)
    n_r, n_a, n_f = length(ray_scales), length(aer_scales), length(sif_vals)
    names_v = (:L_toa, :L_down, :Lw); names_s = (:Ed, :Eu)
    F0v = Dict(k => zeros(n_hres, n_v, n_r, n_a) for k in names_v)   # SIF off
    F1v = Dict(k => zeros(n_hres, n_v, n_r, n_a) for k in names_v)   # SIF with I(678) = 1
    F0s = Dict(k => zeros(n_hres, n_r, n_a) for k in names_s)
    F1s = Dict(k => zeros(n_hres, n_r, n_a) for k in names_s)
    for (ia, a) in enumerate(aer_scales), (ir, s) in enumerate(ray_scales)
        for ib in eachindex(τ0)
            model.optics.τ_rayl[ib] .= s .* τ0[ib]
        end
        m = col === nothing ? model : inject_aerosols(model, params, scale_column(col, a), p_half, trunc)
        t = @elapsed begin
            f0 = fields(m, 0.0); f1 = fields(m, 1.0)
        end
        for k in names_v
            F0v[k][:, :, ir, ia] .= f0[k]; F1v[k][:, :, ir, ia] .= f1[k]
        end
        for k in names_s
            F0s[k][:, ir, ia] .= f0[k]; F1s[k][:, ir, ia] .= f1[k]
        end
        i678 = argmin(abs.(λ_hres .- 678.2))
        println("  aerosol ×$a  Rayleigh ×$s  $(round(t; digits=1)) s   check (SIF off) Lw/Ed at 678 nm = " *
                "$(round(f0.Lw[i678, 1] / f0.Ed[i678]; sigdigits=5)), albedo/π = $(round(surf0.albedo / π; sigdigits=5)), " *
                "Eu/Ed = $(round(f0.Eu[i678] / f0.Ed[i678]; sigdigits=5))")
    end
    for ib in eachindex(τ0)
        model.optics.τ_rayl[ib] .= τ0[ib]
    end

    conv(x) = reshape(kernel.RSR_out * reshape(x, n_hres, :), n_band, size(x)[2:end]...)
    unit_L = "W m-2 sr-1 um-1"; unit_E = "W m-2 um-1"
    mkpath(dirname(SCAN_OUT_NC))
    NCDataset(SCAN_OUT_NC, "c") do ds
        defDim(ds, "band", n_band); defDim(ds, "spectral", n_hres); defDim(ds, "sif", n_f)
        defDim(ds, "view", n_v); defDim(ds, "rayleigh", n_r); defDim(ds, "aerosol", n_a)
        put(name, T, dims, data, attrib) = defVar(ds, name, T.(data), dims; attrib=attrib)
        put("rayleigh_scale", Float64, ("rayleigh",), ray_scales, Dict("long_name" => "factor on the Rayleigh optical depth"))
        put("aerosol_scale", Float64, ("aerosol",), aer_scales, Dict("long_name" => "factor on the GCHPIO aerosol optical depth"))
        put("sif_678", Float64, ("sif",), sif_vals, Dict("units" => unit_L, "long_name" => "water-leaving SIF at $(SIF_λ) nm"))
        put("vza", Float64, ("view",), vza, Dict("units" => "degree"))
        put("vaz", Float64, ("view",), vaz, Dict("units" => "degree", "long_name" => "relative azimuth, vSmartMOM convention"))
        put("wavelength", Float64, ("band",), λ_oci, Dict("units" => "nm", "long_name" => "OCI band center"))
        put("wavelength_hres", Float64, ("spectral",), λ_hres, Dict("units" => "nm"))
        desc = Dict(:L_toa => ("TOA upwelling radiance (R)", unit_L),
                    :L_down => ("downwelling radiance at the surface in the view directions (T)", unit_L),
                    :Lw => ("upwelling radiance just above the surface, incl. SIF (hdr)", unit_L),
                    :Ed => ("downwelling irradiance just above the surface, diffuse + direct (bhr_dw)", unit_E),
                    :Eu => ("upwelling irradiance just above the surface (bhr_uw)", unit_E))
        for k in names_v
            ln, u = desc[k]
            v0 = conv(F0v[k]); v1 = conv(F1v[k])
            withsif = Array{Float64}(undef, n_band, n_f, n_v, n_r, n_a)
            for (j, amp) in enumerate(sif_vals)
                withsif[:, j, :, :, :] .= v0 .+ amp .* (v1 .- v0)
            end
            put("$(k)_nosif", Float32, ("band", "view", "rayleigh", "aerosol"), v0, Dict("units" => u, "long_name" => "$ln, SIF off, OCI-convolved"))
            put("$(k)", Float32, ("band", "sif", "view", "rayleigh", "aerosol"), withsif, Dict("units" => u, "long_name" => "$ln, with SIF, OCI-convolved"))
            put("$(k)_hres_nosif", Float32, ("spectral", "view", "rayleigh", "aerosol"), F0v[k], Dict("units" => u, "long_name" => "$ln, SIF off, RT grid"))
        end
        for k in names_s
            ln, u = desc[k]
            v0 = conv(F0s[k]); v1 = conv(F1s[k])
            withsif = Array{Float64}(undef, n_band, n_f, n_r, n_a)
            for (j, amp) in enumerate(sif_vals)
                withsif[:, j, :, :] .= v0 .+ amp .* (v1 .- v0)
            end
            put("$(k)_nosif", Float32, ("band", "rayleigh", "aerosol"), v0, Dict("units" => u, "long_name" => "$ln, SIF off, OCI-convolved"))
            put("$(k)", Float32, ("band", "sif", "rayleigh", "aerosol"), withsif, Dict("units" => u, "long_name" => "$ln, with SIF, OCI-convolved"))
            put("$(k)_hres_nosif", Float32, ("spectral", "rayleigh", "aerosol"), F0s[k], Dict("units" => u, "long_name" => "$ln, SIF off, RT grid"))
        end
        # Rrs = Lw / Ed on OCI bands (ratio of band-convolved quantities)
        Lw0 = conv(F0v[:Lw]); Lw1 = conv(F1v[:Lw]); Ed0 = conv(F0s[:Ed]); Ed1 = conv(F1s[:Ed])
        Rrs = Array{Float64}(undef, n_band, n_f, n_v, n_r, n_a)
        for (j, amp) in enumerate(sif_vals), iv in 1:n_v
            Rrs[:, j, iv, :, :] .= (Lw0[:, iv, :, :] .+ amp .* (Lw1[:, iv, :, :] .- Lw0[:, iv, :, :])) ./
                                   (Ed0 .+ amp .* (Ed1 .- Ed0))
        end
        put("Rrs_nosif", Float32, ("band", "view", "rayleigh", "aerosol"), Lw0 ./ reshape(Ed0, n_band, 1, n_r, n_a),
            Dict("units" => "sr-1", "long_name" => "Lw / Ed, SIF off, OCI-convolved (Lambertian: albedo/π)"))
        put("Rrs", Float32, ("band", "sif", "view", "rayleigh", "aerosol"), Rrs,
            Dict("units" => "sr-1", "long_name" => "Lw / Ed with SIF, OCI-convolved"))
        put("sif_waterleaving", Float32, ("band", "sif"), (kernel.RSR_out * reverse(I_unit)) .* sif_vals',
            Dict("units" => unit_L, "long_name" => "water-leaving SIF on OCI bands"))
        put("aod", Float32, ("band", "aerosol"), aod_oci .* aer_scales',
            Dict("units" => "1", "long_name" => "GCHPIO column AOD × aerosol_scale (NaN: no aerosol)"))
        ds.attrib["ocean_yaml"] = SCAN_YAML
        ds.attrib["sza"] = Float64(params.sza)
        ds.attrib["surface_albedo"] = Float64(surf0.albedo)
        ds.attrib["sif_library_index"] = SCAN_SIF_INDEX
        ds.attrib["sunlight_on"] = Int8(SUNLIGHT_ON)
        ds.attrib["aer_column"] = AER_COLUMN
        ds.attrib["aerosol_truncation"] = string(typeof(trunc))
        ds.attrib["rt_outputs"] = "R→L_toa, T→L_down, hdr→Lw, bhr_dw→Ed, bhr_uw→Eu; all × n_to_radiance"
    end
    println("Done: $SCAN_OUT_NC")
end

scan_main()
