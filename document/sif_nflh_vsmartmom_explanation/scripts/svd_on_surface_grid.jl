# SVD retrieval (nPC15_npoly3, nPC15_npoly3_O2A) on the noise-free TOA spectra L_toa of the *_grid_surface.nc files,
# Usage: julia --project=<repo> svd_on_surface_grid.jl <outdir>   (writes <outdir>/svd_retrievals.csv)
# with the L1B solar irradiance (as retrieve_rt_grid.jl) and with the RT's own solar spectrum (vSmartMOM, OCI-convolved).
using TOML, NCDatasets, Statistics, LinearAlgebra, Printf, DelimitedFiles
include("/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/build_single_meas.jl")
const D = "/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/full_RT_construction/output_rayleigh_scaling"
const FILES = ["rt_toa_rayleigh_scaling_ocean_lambertian_1010_wavnum0p05_NoRayl_sunlight_on_AerColumn_100_grid_surface.nc",
               "rt_toa_rayleigh_scaling_ocean_lambertian_1010_wavnum0p05_NoRayl_sza60p5_sunlight_on_AerColumn_100_grid_surface.nc"]
const CFG = "/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/configs"
const SOLAR_L1B = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/sample_granule_20240830T131442_new_chl.nc"
sb = readdlm(joinpath(@__DIR__, "rt_solar_oci.csv"), ',')   # vSmartMOM solar, OCI-convolved: band nm, W m-2 um-1
open(joinpath(ARGS[1], "svd_retrievals.csv"), "w") do io
    println(io, "sza,aerosol_scale,rayleigh_scale,vza,sif_678,config,solar,band_nm,sif_true_band,sif_ret,status")
    for f in FILES
        ds = NCDataset(joinpath(D, f)); λa = Float64.(ds["wavelength"][:]); L = Float64.(ds["L_toa"][:, :, :, :, :])
        sw = Float64.(ds["sif_waterleaving"][:, :]); sif = Float64.(ds["sif_678"][:]); vza = Float64.(ds["vza"][:])
        rs = Float64.(ds["rayleigh_scale"][:]); ae = Float64.(ds["aerosol_scale"][:]); sza = Float64(ds.attrib["sza"]); close(ds)
        for cf in ("svd_nPC15_npoly3.toml", "svd_nPC15_npoly3_O2A.toml")
            cfg = TOML.parsefile(joinpath(CFG, cf))
            idx = findall(cfg["spectral"]["lambda_min_nm"] .< λa .< cfg["spectral"]["lambda_max_nm"]); λ = λa[idx]
            sh = prepare_svd_retrieval_setup(cfg, λ); i678 = argmin(abs.(λ .- 678.2))
            sols = (("L1B", load_l1b_solar_on_bands(SOLAR_L1B, λ)),
                    ("RT", (E = [sb[argmin(abs.(sb[:, 1] .- l)), 2] for l in λ], esd = 1.0, pace_path = "vSmartMOM default solar, OCI-convolved")))
            t0 = time()
            for (slab, sol) in sols, ia in eachindex(ae), ir in eachindex(rs), iv in eachindex(vza), k in eachindex(sif)
                y = L[idx, k, iv, ir, ia]
                ret, st = try
                    r = run_pseudo_svd_retrieval(y, λ, sh, sol; sza_deg=sza); (vec(sh.sif_basis * r.sif_coeff)[i678], r.stats.status)
                catch; (NaN, 4) end
                @printf(io, "%.1f,%g,%g,%.0f,%g,%s,%s,%.2f,%.6f,%.6f,%d\n", sza, ae[ia], rs[ir], vza[iv], sif[k],
                        replace(cf, "svd_" => "", ".toml" => ""), slab, λ[i678], sw[idx[i678], k], ret, st)
            end
            @printf("SZA %.1f %s done in %.0f s\n", sza, cf, time() - t0)
        end
    end
end
