# SVD retrieval on the Rayleigh × aerosol × SIF × view grid from generate_rayleigh_scaling_spectra.jl,
# following retrieve_rt_ensemble.jl: fit each spectrum at the grid SZA with the L1B solar irradiance, report
# retrieved water-leaving SIF on the OCI band nearest SIF_REF_NM (678.2 nm). The truth is the point value of
# the grid's SIF shape at that band centre (sif_true); the RSR band average is kept as sif_true_bandavg.
# Every spectrum is fitted clean and with one noise draw (radiance_clean + randn·sigma_noise, NOISE_SEED).
#
#   cd /home/zhe2/FraLab/PACE_redSIF_PACE
#   julia --project=. surrogate_meas/full_RT_construction/output_rayleigh_scaling/retrieve_rt_grid.jl \
#       [svd_nPC15_npoly3.toml svd_nPC15_npoly5.toml ...]
#   GRID_NC=/path/to/*_grid.nc julia --project=. ...
# → rt_retrieval_grid_AerColumn_<k>_<config>.nc next to the grid file

using TOML
using NCDatasets
using Statistics
using LinearAlgebra
using Random
using Printf

include(joinpath(@__DIR__, "..", "..", "build_single_meas.jl"))

const GRID_NC = get(ENV, "GRID_NC", joinpath(@__DIR__,
    "rt_toa_rayleigh_scaling_ocean_lambertian_1010_wavnum0p05_NoRayl_sunlight_on_AerColumn_100_grid.nc"))
const CONFIG_DIR = joinpath(@__DIR__, "..", "..", "configs")
const DEFAULT_CONFIGS = ["svd_nPC15_npoly3.toml", "svd_nPC15_npoly3_O2A.toml",
                         "svd_nPC15_npoly5.toml", "svd_nPC15_npoly5_O2A.toml"]
const SOLAR_L1B = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/sample_granule_20240830T131442_new_chl.nc"
const SIF_REF_NM = 678.2
const SIF_LIB = get(ENV, "SIF_LIB", "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/SIF_singular_vector.jld2")
const DEFAULT_SIF_AMP_NM = 678.0             # grids written before the sif_lambda_nm attribute
const NOISE_SEED = parse(Int, get(ENV, "NOISE_SEED", "20261014"))
const CASES = ("clean", "noisy")

function load_grid(path)
    NCDataset(path) do ds
        (λ = Float64.(ds["wavelength"][:]),
         R_clean = Float64.(ds["radiance_clean"][:, :, :, :, :]),     # band, sif, view, rayleigh, aerosol
         σ = Float64.(ds["sigma_noise"][:, :, :, :, :]),
         sif_toa = Float64.(ds["sif_toa"][:, :, :, :, :]),
         sif_wl = Float64.(ds["sif_waterleaving"][:, :]),             # band, sif
         sif_678 = Float64.(ds["sif_678"][:]),
         vza = Float64.(ds["vza"][:]),
         vaz = Float64.(ds["vaz"][:]),
         rayleigh = Float64.(ds["rayleigh_scale"][:]),
         aerosol = Float64.(ds["aerosol_scale"][:]),
         rayleigh_od = Float64.(ds["rayleigh_od_column"][:]),
         aod = Float64.(ds["aod"][:, :]),                             # band, aerosol
         sza = Float64(ds.attrib["sza"]),
         aer_column = Int(ds.attrib["aer_column"]),
         sif_index = Int(ds.attrib["sif_library_index"]),
         sif_amp_nm = Float64(get(ds.attrib, "sif_lambda_nm", DEFAULT_SIF_AMP_NM)))
    end
end

"Library shape `j` at `λ` relative to its value at `λ_amp`, the generator's normalization wavelength
(same cubic spline as unit_sif_shapes)."
function sif_point_ratio(j, λ, λ_amp)
    wl, shapes = JLD2.jldopen(SIF_LIB) do f
        Float64.(f["SIF_wavelen"]), Float64.(f["SIF_shapes"])
    end
    knots = range(wl[1]; step=wl[2] - wl[1], length=length(wl))
    itp = CubicSplineInterpolation(knots, shapes[:, j]; extrapolation_bc=Line())
    return itp(λ) / itp(λ_amp)
end

function retrieve_config(g, R_obs, toml)
    cfg = TOML.parsefile(joinpath(CONFIG_DIR, toml))
    sp = get(cfg, "spectral", Dict())
    idx = findall(Float64(get(sp, "lambda_min_nm", 640.0)) .< g.λ .< Float64(get(sp, "lambda_max_nm", 756.0)))
    λ = g.λ[idx]
    sh = prepare_svd_retrieval_setup(cfg, λ)
    solar = load_l1b_solar_on_bands(SOLAR_L1B, λ)
    i678 = argmin(abs.(λ .- SIF_REF_NM))
    dims = size(g.R_clean)[2:end]
    ret = Dict(c => fill(NaN, dims) for c in CASES)
    status = Dict(c => fill(Int16(-1), dims) for c in CASES)
    t0 = time()
    for c in CASES, I in CartesianIndices(dims)
        try
            r = run_pseudo_svd_retrieval(R_obs[c][idx, I], λ, sh, solar; sza_deg=g.sza)
            ret[c][I] = vec(sh.sif_basis * r.sif_coeff)[i678]
            status[c][I] = Int16(r.stats.status)
        catch e
            status[c][I] = Int16(4)
            @warn "Retrieval failed" toml c I exception=e
        end
    end
    @printf("%s: %d bands %.1f–%.1f nm, SIF band %.2f nm, %.0f s\n", toml, length(λ), λ[1], λ[end], λ[i678], time() - t0)
    truth = g.sif_678 .* sif_point_ratio(g.sif_index, λ[i678], g.sif_amp_nm)
    return (ret=ret, status=status, band=idx[i678], λ_band=λ[i678], λ_min=λ[1], λ_max=λ[end],
            n_pc=sh.n_pc, n_leg=sh.n_leg, truth=truth)
end

function summarize(g, r, toml)
    truth = r.truth
    for c in CASES
        e = r.ret[c] .- truth
        @printf("  [%s] converged %d/%d  bias %+.4f  RMSE %.4f\n", c, count(==(1), r.status[c]),
                length(r.status[c]), mean(e), sqrt(mean(e .^ 2)))
    end
end

function write_result(g, r, toml)
    tag = replace(splitext(toml)[1], "svd_" => "")
    out = joinpath(dirname(GRID_NC), "rt_retrieval_grid_AerColumn_$(g.aer_column)_$(tag).nc")
    isfile(out) && rm(out)
    unit = "W m-2 sr-1 um-1"
    band_txt = @sprintf("OCI band at %.2f nm", r.λ_band)
    dims4 = ("sif", "view", "rayleigh", "aerosol")
    NCDataset(out, "c") do o
        defDim(o, "sif", length(g.sif_678)); defDim(o, "view", length(g.vza))
        defDim(o, "rayleigh", length(g.rayleigh)); defDim(o, "aerosol", length(g.aerosol))
        defVar(o, "sif_678", g.sif_678, ("sif",); attrib=Dict("units" => unit,
               "long_name" => "grid SIF amplitude at exactly 678.0 nm (high-res)"))
        defVar(o, "vza", g.vza, ("view",); attrib=Dict("units" => "degree"))
        defVar(o, "vaz", g.vaz, ("view",); attrib=Dict("units" => "degree"))
        defVar(o, "rayleigh_scale", g.rayleigh, ("rayleigh",))
        defVar(o, "aerosol_scale", g.aerosol, ("aerosol",))
        defVar(o, "rayleigh_od", g.rayleigh .* g.rayleigh_od[r.band], ("rayleigh",);
               attrib=Dict("long_name" => "column Rayleigh optical depth, $band_txt"))
        defVar(o, "aod", replace(g.aod[r.band, :], NaN => 0.0), ("aerosol",);
               attrib=Dict("long_name" => "aerosol optical depth, $band_txt"))
        defVar(o, "sif_true", r.truth, ("sif",);
               attrib=Dict("units" => unit, "long_name" => "true water-leaving SIF, point value at the centre of the $band_txt"))
        defVar(o, "sif_true_bandavg", g.sif_wl[r.band, :], ("sif",);
               attrib=Dict("units" => unit, "long_name" => "true water-leaving SIF, RSR-averaged over the $band_txt"))
        defVar(o, "sif_toa", g.sif_toa[r.band, :, :, :, :], dims4;
               attrib=Dict("units" => unit, "long_name" => "TOA SIF signal (with − without SIF), $band_txt"))
        for c in CASES
            defVar(o, "sif_ret_$c", r.ret[c], dims4;
                   attrib=Dict("units" => unit, "long_name" => "retrieved water-leaving SIF ($c spectra), $band_txt"))
            defVar(o, "status_$c", r.status[c], dims4; attrib=Dict("long_name" => "retrieval status, 1 = converged"))
        end
        o.attrib["band_nm"] = r.λ_band
        o.attrib["sif_truth"] = "point value of library shape $(g.sif_index) at band_nm, scaled to sif_678 at $(g.sif_amp_nm) nm"
        o.attrib["fit_window_nm"] = @sprintf("%.1f-%.1f", r.λ_min, r.λ_max)
        o.attrib["svd_config"] = toml
        o.attrib["n_pc"] = r.n_pc
        o.attrib["n_legendre"] = r.n_leg
        o.attrib["input"] = GRID_NC
        o.attrib["noise_seed"] = NOISE_SEED
        o.attrib["sza"] = g.sza
        o.attrib["aer_column"] = g.aer_column
    end
    println("  wrote $out")
end

function main(configs)
    g = load_grid(GRID_NC)
    R_obs = Dict("clean" => g.R_clean,
                 "noisy" => g.R_clean .+ randn(MersenneTwister(NOISE_SEED), size(g.R_clean)) .* g.σ)
    for toml in configs
        r = retrieve_config(g, R_obs, toml)
        summarize(g, r, toml)
        write_result(g, r, toml)
    end
end

main(isempty(ARGS) ? DEFAULT_CONFIGS : ARGS)
