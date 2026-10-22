# Averaging kernel (AK) and posterior covariance for selected pixels of a retrieved granule.
#
# Re-solves each selected pixel with the pipeline's own `_run_one_svd_retrieval!`
# (global_svd_fit_pipeline/svd_retrieval/svd_helpers.jl), which returns
#   averaging_kernel = S_post * H_obs  (= G*K at the final state),  S_posterior,  dof.
# The swath output only stores DOF and diag(S_post), so the full matrices need this re-run.
#
#   julia --project=/home/zhe2/FraLab/PACE_redSIF_PACE -t 4 retrieval_comparison/ak_selected_pixels.jl \
#         [--runs npoly3,newSVD_npoly3] [--granule 20250724T170203] \
#         [--pixels retrieval_comparison/out/ak/selected_pixels.csv] [--config-dir DIR]
#
# Configs are looked up in <repo>/svd_configs/global_fit_new_Sep26, then in the main checkout's
# copy (some are untracked in git), unless --config-dir is given.
#
# The pipeline deletes its interim radiance file after each run, so it is regenerated
# once with the pipeline's own merge function and cached in out/_interim/ (~1 GB, gitignored).
# Per pixel the setup copies `_process_one_pixel_svd!` (svd_swath_parallel.jl).

using Base.Threads
using Dates
using DelimitedFiles
using ForwardDiff
using Glob
using LinearAlgebra
using NCDatasets
using Printf
using ProgressMeter
using TOML

const HERE = @__DIR__
const REPO = abspath(joinpath(HERE, ".."))
include(joinpath(REPO, "global_svd_fit_pipeline", "svd_retrieval", "svd_granule.jl"))

const CFG_DIR = joinpath(REPO, "svd_configs", "global_fit_new_Sep26")
# Some single-granule configs are untracked in git and live only in the main checkout.
const CFG_DIR_FALLBACK = "/home/zhe2/FraLab/PACE_redSIF_PACE/svd_configs/global_fit_new_Sep26"
const RUN_CFG = Dict(
    "npoly2" => "new.singleGranule.npoly2nPC15_O2A.toml",
    "npoly3" => "new.singleGranule.npoly3nPC15_O2A.toml",
    "npoly5" => "new.singleGranule.npoly5nPC15_O2A.toml",
    "newSVD_npoly2" => "new.newSVD.singleGranule.npoly2nPC15_O2A.toml",
    "newSVD_npoly3" => "new.newSVD.singleGranule.npoly3nPC15_O2A.toml",
)
const RUN_OUTDIR = Dict(
    "npoly2" => "new_svd_retrieval_output_npoly2",
    "npoly3" => "new_svd_retrieval_output_O2A",
    "npoly5" => "new_svd_retrieval_output_npoly5",
    "newSVD_npoly2" => "new_svd_retrieval_output_O2A_newSVD_npoly2",
    "newSVD_npoly3" => "new_svd_retrieval_output_O2A_newSVD",
)
const PACE = "/home/zhe2/data/PACE"
const KIWI = "/kiwi-data/Data/satellite/PACE_OCI"

function parse_args(args)
    o = Dict("runs" => "npoly3", "granule" => "20250724T170203", "config-dir" => "",
             "pixels" => joinpath(HERE, "out", "ak", "selected_pixels.csv"))
    i = 1
    while i <= length(args)
        a = args[i]
        startswith(a, "--") && i < length(args) || error("bad argument $a")
        o[a[3:end]] = args[i + 1]
        i += 2
    end
    return o
end

function day_dir(root, gid)
    return joinpath(root, gid[1:4], gid[5:6], gid[7:8])
end

"""Retrieval config dict as the pipeline builds it: every TOML section except orchestration ones."""
function retrieval_cfg(toml_path)
    pipe = TOML.parsefile(toml_path)
    return Dict{String, Any}(k => v for (k, v) in pipe if !(k in ("download", "general", "svd_retrieval")))
end

function ak_one_pixel(sh, y_obs, sza_deg, max_outer_steps)
    solar_eff = @. sh.solar_band_ctx * cosd(sza_deg) / sh.es
    fm, layout = make_svd_forward_model_λ(sh.λ_ctx, solar_eff, sh.PCs, sh.sif_basis;
        n_pc = sh.n_pc, n_legendre = sh.n_legendre, log_transform = sh.log_transform)
    jac_eval = (x -> ForwardDiff.jacobian(fm, x))
    x_a = copy(sh.x0)
    σ_prior = copy(sh.prior_sigma_template)
    if sh.use_leg01 && length(layout.idx_legendre) >= 1
        y0 = fm(sh.x0)
        ratio = y_obs ./ max.(abs.(y0), eps(Float64))
        w = y_obs .- minimum(y_obs)
        w .+= max(maximum(w), 1.0) * 1e-6
        sv = sqrt.(w ./ maximum(w))
        c01 = (sh.A01 .* sv) \ (ratio .* sv)
        leg0 = first(layout.idx_legendre)
        x_a[leg0] = c01[1]
        σ_prior[leg0] = max(abs(c01[1]) * sh.leg01_frac, sh.prior_min_sigma)
        if length(layout.idx_legendre) >= 2
            leg1 = layout.idx_legendre[2]
            x_a[leg1] = c01[2]
            σ_prior[leg1] = max(abs(c01[2]) * sh.leg01_frac, sh.prior_min_sigma)
        end
    end
    x_tmp = zeros(Float64, layout.n_state)
    stats = _run_one_svd_retrieval!(
        x_tmp, fm, jac_eval, layout, y_obs, copy(sh.x0), x_a, σ_prior, sh.lower, sh.upper,
        sh.x_scale_template, sh.use_leg01, sh.A01, sh.prior_min_sigma, sh.leg01_frac,
        sh.use_band_snr, sh.band_snr_coeffs, sh.meas_sigma, sh.lm,
        (dx_rel_tol = sh.conv.dx_rel_tol, rmse_rel_tol = sh.conv.rmse_rel_tol,
         rmse_abs_tol = sh.conv.rmse_abs_tol, enabled = sh.conv.enabled, window = sh.conv.window,
         redchi2_target = sh.conv.redchi2_target, redchi2_abs_tol = sh.conv.redchi2_abs_tol,
         redchi2_rel_tol = sh.conv.redchi2_rel_tol),
        sh.conv.dx_rel_tol_stall, max_outer_steps)
    return (; x = x_tmp, x_a, σ_prior, stats, y_fit = fm(x_tmp), layout)
end

function main()
    opt = parse_args(ARGS)
    gid = opt["granule"]
    runs = String.(split(opt["runs"], ","))
    sel, hdr = readdlm(opt["pixels"], ','; header = true)
    hdr = vec(String.(hdr))
    col(name) = sel[:, findfirst(==(name), hdr)]
    cats = String.(col("category"))
    scans = Int.(col("scan"))
    pixs = Int.(col("pixel"))
    nsel = length(scans)
    cat_names = unique(cats)
    println("$(nsel) selected pixels, runs = $(runs)")

    l1b = joinpath(day_dir(joinpath(KIWI, "L1B_V3"), gid), "PACE_OCI.$(gid).L1B.V3.nc")
    aop = joinpath(day_dir(joinpath(KIWI, "L2_AOP_V3.2"), gid), "PACE_OCI.$(gid).L2.OC_AOP.V3_2.nc")
    bgc = joinpath(day_dir(joinpath(KIWI, "L2_BGC_V3.2"), gid), "PACE_OCI.$(gid).L2.OC_BGC.V3_2.nc")
    interim = joinpath(HERE, "out", "_interim", "interim_$(gid).nc")
    mkpath(dirname(interim))
    if !isfile(interim)
        println("regenerating interim radiances with the pipeline merge -> $interim")
        preprocess_and_merge_in_memory(l1b, aop, bgc, interim)
    end

    outdir = joinpath(HERE, "out", "ak")
    mkpath(outdir)
    for run in runs
        dirs = isempty(opt["config-dir"]) ? (CFG_DIR, CFG_DIR_FALLBACK) : (opt["config-dir"],)
        cands = [joinpath(d, RUN_CFG[run]) for d in dirs]
        i_cfg = findfirst(isfile, cands)
        i_cfg === nothing && error("config $(RUN_CFG[run]) not found in $(collect(dirs)); pass --config-dir")
        cfgpath = cands[i_cfg]
        println("[$run] config: $cfgpath")
        cfg = retrieval_cfg(cfgpath)
        sh = _build_svd_swath_shared(interim, l1b, cfg)
        max_outer = Int(get(get(cfg, "batch_fit", Dict()), "max_outer_steps", 10))
        nst, nb = sh.layout.n_state, length(sh.λ_ctx)
        stored = NCDataset(joinpath(PACE, RUN_OUTDIR[run], gid[1:4], gid[5:6], gid[7:8],
                                    "interim_$(gid)_svd_retrieval.nc"))

        A = zeros(nst, nst, nsel); S = zeros(nst, nst, nsel)
        X = zeros(nst, nsel); XA = zeros(nst, nsel); SP = zeros(nst, nsel)
        YO = zeros(nb, nsel); YF = zeros(nb, nsel)
        dof = zeros(nsel); trA = zeros(nsel); rchi2 = zeros(nsel); conv = zeros(Int8, nsel)
        status = zeros(Int16, nsel); sif_rerun = zeros(nsel); sif_stored = zeros(nsel)
        conv_stored = zeros(Int8, nsel)

        dsi = NCDataset(interim)
        v = dsi["Rtoa_red"]
        ax = _find_axis_indices(v, "red_wavelength")
        y_sorted = zeros(Float64, length(sh.perm))
        for k in 1:nsel
            inds = ntuple(d -> d == ax.i_pix ? pixs[k] : d == ax.i_scan ? scans[k] : Colon(), ndims(v))
            spec = vec(v[inds...])
            _copy_sorted_spectrum!(y_sorted, spec, sh.perm) || error("pixel $(pixs[k]),$(scans[k]) has missing radiances")
            y_obs = y_sorted[sh.slice_sorted]
            r = ak_one_pixel(sh, copy(y_obs), sh.sza[pixs[k], scans[k]], max_outer)
            A[:, :, k] = r.stats.averaging_kernel
            S[:, :, k] = r.stats.S_posterior
            X[:, k] = r.x; XA[:, k] = r.x_a; SP[:, k] = r.σ_prior
            YO[:, k] = y_obs; YF[:, k] = r.y_fit
            dof[k] = r.stats.dof; trA[k] = r.stats.ak_trace; rchi2[k] = r.stats.reduced_chi2
            conv[k] = r.stats.converged; status[k] = r.stats.status
            sif_rerun[k] = r.x[last(r.layout.idx_sif)]
            sif_stored[k] = Float64(stored["sif_ev1"][pixs[k], scans[k]])
            conv_stored[k] = Int8(stored["converged"][pixs[k], scans[k]])
            println(@sprintf("  [%s] %-20s scan %4d pix %4d  conv=%d (stored %d)  DOF_sig=%.2f  A[sif,sif]=%.3f  sif %.4f (stored %.4f)",
                run, cats[k], scans[k], pixs[k], conv[k], conv_stored[k], trA[k],
                A[last(r.layout.idx_sif), last(r.layout.idx_sif), k], sif_rerun[k], sif_stored[k]))
        end
        close(dsi); close(stored)

        outnc = joinpath(outdir, "ak_$(run)_$(gid).nc")
        NCDataset(outnc, "c") do ds
            defDim(ds, "sel", nsel); defDim(ds, "row", nst); defDim(ds, "col", nst)
            defDim(ds, "state", nst); defDim(ds, "band", nb)
            # Julia arrays are written so that Python reads A[sel, i, j] = d x_hat_i / d x_true_j
            defVar(ds, "A", Float64, ("col", "row", "sel"))[:, :, :] = permutedims(A, (2, 1, 3))
            defVar(ds, "S_post", Float64, ("col", "row", "sel"))[:, :, :] = permutedims(S, (2, 1, 3))
            defVar(ds, "x_hat", Float64, ("state", "sel"))[:, :] = X
            defVar(ds, "x_prior", Float64, ("state", "sel"))[:, :] = XA
            defVar(ds, "sigma_prior", Float64, ("state", "sel"))[:, :] = SP
            defVar(ds, "wavelength", Float64, ("band",))[:] = sh.λ_ctx
            defVar(ds, "y_obs", Float64, ("band", "sel"))[:, :] = YO
            defVar(ds, "y_fit", Float64, ("band", "sel"))[:, :] = YF
            defVar(ds, "dof_signal", Float64, ("sel",))[:] = trA
            defVar(ds, "dof_residual", Float64, ("sel",))[:] = dof
            defVar(ds, "reduced_chi2", Float64, ("sel",))[:] = rchi2
            defVar(ds, "converged", Int8, ("sel",))[:] = conv
            defVar(ds, "converged_stored", Int8, ("sel",))[:] = conv_stored
            defVar(ds, "status", Int16, ("sel",))[:] = status
            defVar(ds, "sif_ev1_rerun", Float64, ("sel",))[:] = sif_rerun
            defVar(ds, "sif_ev1_stored", Float64, ("sel",))[:] = sif_stored
            defVar(ds, "scan", Int32, ("sel",))[:] = scans
            defVar(ds, "pixel", Int32, ("sel",))[:] = pixs
            defVar(ds, "category_id", Int8, ("sel",))[:] = [findfirst(==(c), cat_names) - 1 for c in cats]
            defVar(ds, "chl", Float64, ("sel",))[:] = Float64.(col("chl"))
            defVar(ds, "nflh", Float64, ("sel",))[:] = Float64.(col("nflh"))
            ds.attrib["state_names_csv"] = join(sh.state_names, ",")
            ds.attrib["category_names_csv"] = join(cat_names, ",")
            ds.attrib["run"] = run
            ds.attrib["config_file"] = cfgpath
            ds.attrib["granule"] = gid
            ds.attrib["note"] = "A = S_post*H_obs at the final state; A[sel,i,j] = d xhat_i / d x_true_j"
        end
        println("wrote $outnc")
    end
end

main()
