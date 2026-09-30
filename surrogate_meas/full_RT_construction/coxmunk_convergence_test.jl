# Cox-Munk glint convergence test vs number of streams.
#
# Rayleigh-only atmosphere (no aerosols, no SIF) over a Cox-Munk ocean, in a narrow
# near-IR continuum window where glint dominates. For each `nstreams` the TOA
# reflectance is computed along the principal plane; the glint peak should stop
# changing once the angular resolution (streams + Fourier terms) is sufficient.
#
# Usage (from the repo root):
#   julia -t 8 surrogate_meas/full_RT_construction/coxmunk_convergence_test.jl
# Environment overrides: NSTREAMS="6,11,16", WIND_SPEED=5.0, OUT_DIR=..., OCEAN_YAML=...

using Pkg
Pkg.activate("/home/zhe2/FraLab/vSmartMOM.jl")
using vSmartMOM
using vSmartMOM.SolarModel
using NCDatasets
using Dates

# configs/ is gitignored, so default to the main checkout's copy
const OCEAN_YAML = get(ENV, "OCEAN_YAML",
    "/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/configs/ocean_coxmunk_0912.yaml")
const OUT_DIR = get(ENV, "OUT_DIR", joinpath(@__DIR__, "output_new_CoxMunk"))

const NSTREAMS = parse.(Int, split(get(ENV, "NSTREAMS", "6,11,16,24,32,40,48"), ","))
const WIND_SPEED = parse(Float64, get(ENV, "WIND_SPEED", "5.0"))
const INCLUDE_WHITECAPS = true
const WHITECAP_ALBEDO = 0.22
const SZA = 30.0
# principal plane: forward (vaz=0, glint side in vSmartMOM convention) and backward (vaz=180)
const VZA = vcat(collect(0.0:5.0:70.0), [15.0, 30.0, 45.0, 60.0])
const VAZ = vcat(zeros(15), fill(180.0, 4))
# near-IR continuum window (cm⁻¹, ascending): 870 → 860 nm
const SPEC_BAND = "(1e7/870):1:(1e7/860)"

const h = 6.62607015e-34
const c = 299792458.0

function yaml_with(nstreams::Int)
    txt = read(OCEAN_YAML, String)
    txt = replace(txt, r"nstreams:\s*\d+" => "nstreams:           $(nstreams)")
    txt = replace(txt, r"- \"\(1e7/850\):1:\(1e7/623\)\"" => "- \"$(SPEC_BAND)\"")
    occursin("nstreams:           $(nstreams)", txt) || error("could not set nstreams in YAML")
    occursin(SPEC_BAND, txt) || error("could not set spec_bands in YAML")
    path = tempname() * ".yaml"
    write(path, txt)
    return path
end

function run_one(nstreams::Int)
    params = parameters_from_yaml(yaml_with(nstreams))
    surf0 = only(params.brdf)
    FT = typeof(surf0.wind_speed)
    params.brdf = [vSmartMOM.CoreRT.CoxMunkSurface{FT}(
        wind_speed = FT(WIND_SPEED), n_water = surf0.n_water,
        whitecap_albedo = FT(WHITECAP_ALBEDO), include_whitecaps = INCLUDE_WHITECAPS,
        shadowing = surf0.shadowing)]
    params.scattering_params = nothing        # Rayleigh only
    params.sza = FT(SZA)
    params.vza = FT.(VZA)
    params.vaz = FT.(VAZ)

    model = model_from_parameters(params)
    ν = params.spec_bands[1]
    F_sol = SolarModel.default_solar_spectrum_at_earth(ν)[:, 2]
    F₀ = zeros(params.polarization_type.n, length(ν))
    F₀[1, :] .= F_sol

    t0 = time()
    R, = rt_run(model; sources = SolarBeam(F₀ = F₀))
    dt = time() - t0

    n_to_rad = @. 100 * h * c * ν                       # photons → W m⁻² sr⁻¹ µm⁻¹
    L = R[:, 1, :] .* reshape(n_to_rad, 1, :)            # (view, ν)
    Fsol_W = F_sol .* n_to_rad                            # W m⁻² µm⁻¹
    ρ = π .* L ./ (reshape(Fsol_W, 1, :) .* cosd(SZA))   # TOA reflectance
    return (; L = vec(sum(L; dims = 2) ./ size(L, 2)),
              ρ = vec(sum(ρ; dims = 2) ./ size(ρ, 2)),
              m_max = maximum(model.solver.m_max_bands),
              Nquad = model.quad_points.Nquad, seconds = dt)
end

mkpath(OUT_DIR)
println("Cox-Munk convergence: wind=$(WIND_SPEED) m/s, sza=$(SZA), nstreams=$(NSTREAMS), window $(SPEC_BAND) cm⁻¹")
println("YAML template: $(OCEAN_YAML)")
results = []
ig = findfirst(i -> VZA[i] == SZA && VAZ[i] == 0, eachindex(VZA))
ib = findfirst(i -> VZA[i] == 30.0 && VAZ[i] == 180, eachindex(VZA))
for n in NSTREAMS
    r = run_one(n)
    push!(results, r)
    println("nstreams=$(lpad(n, 3))  m_max=$(lpad(r.m_max, 3))  Nquad=$(lpad(r.Nquad, 3))  " *
            "ρ_glint(vza=30,fwd)=$(round(r.ρ[ig]; digits = 4))  ρ(vza=30,back)=$(round(r.ρ[ib]; digits = 4))  " *
            "t=$(round(r.seconds; digits = 1)) s")
    flush(stdout)
end

tag = replace(string(WIND_SPEED), "." => "p")
out = joinpath(OUT_DIR, "coxmunk_convergence_ws$(tag)_n$(join(NSTREAMS, "-")).nc")
isfile(out) && rm(out)
NCDataset(out, "c") do ds
    defDim(ds, "view", length(VZA))
    defDim(ds, "nstreams", length(NSTREAMS))
    defVar(ds, "nstreams", Int32.(NSTREAMS), ("nstreams",))
    defVar(ds, "vza", VZA, ("view",); attrib = Dict("units" => "degree"))
    defVar(ds, "vaz", VAZ, ("view",); attrib = Dict("units" => "degree",
        "comment" => "vSmartMOM convention: 0 = forward (glint side), 180 = backscatter"))
    defVar(ds, "reflectance", hcat([r.ρ for r in results]...), ("view", "nstreams");
           attrib = Dict("long_name" => "TOA reflectance pi L / (F0 mu0), window mean"))
    defVar(ds, "radiance", hcat([r.L for r in results]...), ("view", "nstreams");
           attrib = Dict("units" => "W m-2 sr-1 um-1", "long_name" => "TOA radiance, window mean"))
    defVar(ds, "m_max", Int32.([r.m_max for r in results]), ("nstreams",))
    defVar(ds, "Nquad", Int32.([r.Nquad for r in results]), ("nstreams",))
    defVar(ds, "runtime_s", [r.seconds for r in results], ("nstreams",))
    ds.attrib["sza"] = SZA
    ds.attrib["wind_speed"] = WIND_SPEED
    ds.attrib["include_whitecaps"] = Int8(INCLUDE_WHITECAPS)
    ds.attrib["whitecap_albedo"] = WHITECAP_ALBEDO
    ds.attrib["spectral_window_cm-1"] = SPEC_BAND
    ds.attrib["atmosphere"] = "Rayleigh + O2/N2 absorption (YAML profile), no aerosols, no SIF"
    ds.attrib["yaml_template"] = OCEAN_YAML
    ds.attrib["created"] = string(now())
end
println("wrote $out")
