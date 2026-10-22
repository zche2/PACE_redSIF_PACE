# Averaging kernel (AK) and posterior covariance for selected samples of the synthetic
# GCHPIO full-RT ensemble (Lambertian ocean, GCHPIO TOMAS aerosol optics, 0.05 cm-1 grid).
#
# Re-solves each sample with `run_pseudo_svd_retrieval` (surrogate_meas/build_single_meas.jl),
# the same call retrieve_rt_ensemble.jl used, and keeps `stats.averaging_kernel` and
# `stats.S_posterior`. The config of each run is taken from the stored retrieval file's
# `svd_config` attribute, so it is exactly the one that produced it.
#
#   julia --project=/home/zhe2/FraLab/PACE_redSIF_PACE -t 4 retrieval_comparison/ak_gchpio_ensemble.jl \
#         [--runs gchpio_npoly2,gchpio_npoly3,gchpio_newSVD_npoly3] [--n 30] [--seed 0] [--min-total 120]
#         [--samples-from retrieval_comparison/out/ak/ak_gchpio_npoly3_gchpio.nc]
#   --samples-from: reuse the samples/categories of an existing AK file (required for runs with
#   no stored ensemble retrieval, e.g. gchpio_npoly5, which only needs its config).
#   --n: samples per category; if categories run short (few samples converge in all runs), the
#   selection is topped up round-robin from categories with spare candidates to --min-total.
#
# Sample categories (true ensemble state; converged in every run listed in --runs):
#   clean_nadir  AOD678 <= 25th pct and VZA < 25     high_aod  AOD678 >= 90th pct
#   high_vza     VZA >= 50                           low_sif   true SIF678 < 0.05
#   nonconverged not converged in the first run of --runs (contrast case)
# Output: out/ak/ak_<run>_gchpio.nc (same layout as ak_selected_pixels.jl; `scan` holds the
# 1-based sample index, `chl` holds AOD678 (attribute meta_label), plus the SIF truth).

using LinearAlgebra
using NCDatasets
using Printf
using Random
using Statistics
using TOML

const HERE = @__DIR__
const REPO = abspath(joinpath(HERE, ".."))
include(joinpath(REPO, "surrogate_meas", "build_single_meas.jl"))

# Ensemble outputs and surrogate configs are not in git; they live in the main checkout.
const ENS_DIR = "/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/full_RT_construction/output_new_Lambertian_gchpio_finerwavnum0p05"
const RT_NC = joinpath(ENS_DIR, "rt_toa_ensemble_lambertian_gchpio_n1000.nc")
const SOLAR_L1B = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/sample_granule_20240830T131442_new_chl.nc"
const SIF_REF_NM = 678.2
const STORED = Dict(
    "gchpio_npoly2" => "rt_retrieval_1000_O2A_SIF678BandMismatch_nPC15_npoly2.nc",
    "gchpio_npoly3" => "rt_retrieval_1000_O2A_SIF678BandMismatch_nPC15_npoly3.nc",
    "gchpio_newSVD_npoly3" => "rt_retrieval_1000_O2A_ScatteringIncludedSVD_nPC15_npoly3.nc",
    "gchpio_noO2A_npoly3" => "rt_retrieval_1000_SIF678BandMismatch_nPC15_npoly3.nc",
)
# Runs with no stored ensemble retrieval: only the config is needed (the AK re-run is the retrieval).
# They cannot drive the sample selection, so use --samples-from with them.
const SURR_CFG = "/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/configs"
const CONFIG_ONLY = Dict(
    "gchpio_npoly5" => joinpath(SURR_CFG, "svd_nPC15_npoly5_O2A.toml"),
)
const TRUTH_RUN = "gchpio_npoly3"   # every stored file carries the same sif_678_true

function parse_args(args)
    o = Dict("runs" => "gchpio_npoly2,gchpio_npoly3,gchpio_newSVD_npoly3", "n" => "30", "seed" => "0",
             "min-total" => "120")
    i = 1
    while i <= length(args)
        a = args[i]
        startswith(a, "--") && i < length(args) || error("bad argument $a")
        o[a[3:end]] = args[i + 1]
        i += 2
    end
    return o
end

stored_path(run) = joinpath(ENS_DIR, STORED[run])

"""AOD at the SIF band and VZA for every ensemble sample."""
function ens_meta()
    ens = NCDataset(RT_NC)
    w = Float64.(ens["wavelength"][:])
    aod = Float64.(ens["aod"][argmin(abs.(w .- SIF_REF_NM)), :])
    vza = Float64.(ens["vza"][:])
    close(ens)
    return aod, vza
end

function truth_sif678()
    d = NCDataset(stored_path(TRUTH_RUN)); t = Float64.(d["sif_678_true"][:]); close(d)
    return t
end

"""Samples and categories of an existing AK file, so new runs use exactly the same samples."""
function samples_from(path)
    d = NCDataset(path)
    picks = Int.(d["scan"][:])
    names = split(d.attrib["category_names_csv"], ",")
    cats = String[names[i + 1] for i in Int.(d["category_id"][:])]
    close(d)
    return picks, cats
end

"""Prior mean and sigma exactly as run_pseudo_svd_retrieval builds them (Legendre p0/p1 from the spectrum)."""
function prior_for(sh, fm, layout, y_obs)
    x_a = copy(sh.x0)
    σ = copy(sh.prior_sigma)
    if sh.use_leg01 && length(layout.idx_legendre) >= 1
        y0 = fm(sh.x0)
        ratio = y_obs ./ max.(abs.(y0), eps(Float64))
        w = y_obs .- minimum(y_obs)
        w .+= max(maximum(w), 1.0) * 1e-6
        sv = sqrt.(w ./ maximum(w))
        c01 = (sh.A01 .* sv) \ (ratio .* sv)
        leg0 = first(layout.idx_legendre)
        x_a[leg0] = c01[1]
        σ[leg0] = max(abs(c01[1]) * sh.leg01_frac, sh.prior_min_sigma)
        if length(layout.idx_legendre) >= 2
            leg1 = layout.idx_legendre[2]
            x_a[leg1] = c01[2]
            σ[leg1] = max(abs(c01[2]) * sh.leg01_frac, sh.prior_min_sigma)
        end
    end
    return x_a, σ
end

function select_samples(runs, n, seed, min_total)
    ens = NCDataset(RT_NC)
    w = Float64.(ens["wavelength"][:])
    ib = argmin(abs.(w .- SIF_REF_NM))
    aod = Float64.(ens["aod"][ib, :])
    vza = Float64.(ens["vza"][:])
    close(ens)
    bad = filter(r -> !haskey(STORED, r), runs)
    isempty(bad) || error("runs $bad have no stored retrieval; pass --samples-from <existing ak_*_gchpio.nc>")
    conv = [begin
                d = NCDataset(stored_path(r)); s = Int.(d["status"][:]); close(d); s .== 1
            end for r in runs]
    d = NCDataset(stored_path(runs[1])); sif_true = Float64.(d["sif_678_true"][:]); close(d)
    all_conv = reduce(.&, conv)
    q25, q90 = quantile(aod, 0.25), quantile(aod, 0.90)
    cats = [
        ("clean_nadir", all_conv .& (aod .<= q25) .& (vza .< 25)),
        ("high_aod", all_conv .& (aod .>= q90)),
        ("high_vza", all_conv .& (vza .>= 50)),
        ("low_sif", all_conv .& (sif_true .< 0.05)),
        ("nonconverged", .!conv[1]),
    ]
    rng = MersenneTwister(seed)
    picks, names = Int[], String[]
    for (c, m) in cats
        idx = findall(m)
        println(@sprintf("  %-13s candidates %4d", c, length(idx)))
        isempty(idx) && continue
        for k in shuffle(rng, idx)[1:min(n, length(idx))]
            push!(picks, k); push!(names, c)
        end
    end
    # Top up to min_total from categories that still have unpicked candidates (round robin).
    pools = [(c, shuffle(rng, setdiff(findall(m), picks))) for (c, m) in cats]
    while length(picks) < min_total && any(!isempty(p) for (_, p) in pools)
        for (c, p) in pools
            (isempty(p) || length(picks) >= min_total) && continue
            k = popfirst!(p)
            k in picks && continue
            push!(picks, k); push!(names, c)
        end
    end
    order = sortperm(names; by = c -> findfirst(==(c), first.(cats)))
    return picks[order], names[order], aod, vza, sif_true
end

function main()
    opt = parse_args(ARGS)
    runs = String.(split(opt["runs"], ","))
    if !isempty(get(opt, "samples-from", ""))
        picks, cats = samples_from(opt["samples-from"])
        aod, vza = ens_meta()
        sif_true_all = truth_sif678()
        println("samples reused from $(opt["samples-from"])")
    else
        picks, cats, aod, vza, sif_true_all = select_samples(runs, parse(Int, opt["n"]), parse(Int, opt["seed"]),
                                                             parse(Int, opt["min-total"]))
    end
    cat_names = unique(cats)
    nsel = length(picks)
    println("$(nsel) samples selected, runs = $(runs)")
    outdir = joinpath(HERE, "out", "ak")
    mkpath(outdir)

    for run in runs
        if haskey(STORED, run)
            st = NCDataset(stored_path(run))
            cfgpath = String(st.attrib["svd_config"])
            sif_stored_all = Float64.(st["sif_678_ret"][:])
            conv_stored_all = Int8.(Int.(st["status"][:]) .== 1)
            close(st)
        else
            cfgpath = CONFIG_ONLY[run]
            sif_stored_all = fill(NaN, length(sif_true_all))
            conv_stored_all = fill(Int8(-1), length(sif_true_all))   # -1 = no stored retrieval
        end
        cfg = TOML.parsefile(cfgpath)
        sp = get(cfg, "spectral", Dict())
        λ_min = Float64(get(sp, "lambda_min_nm", 640.0)); λ_max = Float64(get(sp, "lambda_max_nm", 756.0))

        ens = NCDataset(RT_NC)
        λ_all = Float64.(ens["wavelength"][:])
        idx = findall(λ_min .< λ_all .< λ_max)
        λ = λ_all[idx]
        Y = Float64.(ens["radiance_noisy"][idx, :])
        sza = Float64.(ens["sza"][:])
        alb = Float64.(ens["surface_albedo"][:])
        close(ens)

        sh = prepare_svd_retrieval_setup(cfg, λ)
        solar = load_l1b_solar_on_bands(SOLAR_L1B, λ)
        i678 = argmin(abs.(λ .- SIF_REF_NM))
        w678 = sh.sif_basis[i678, 1]
        println("[$run] config $(basename(cfgpath)); $(length(λ)) bands in ($λ_min, $λ_max) nm")

        nst, nb = 0, length(λ)
        A = S = X = XA = SP = YF = nothing
        lay = nothing
        res = Dict{String, Vector{Float64}}(k => zeros(nsel) for k in
            ("trA", "dof", "rchi2", "sif_rerun", "sif_stored", "sif678_rerun", "sif678_stored", "sif678_true", "sigma678"))
        conv = zeros(Int8, nsel); conv_st = zeros(Int8, nsel); status = zeros(Int16, nsel)
        for (k, s) in enumerate(picks)
            r = run_pseudo_svd_retrieval(Y[:, s], λ, sh, solar; sza_deg = sza[s])
            if A === nothing
                nst = r.layout.n_state
                A = zeros(nst, nst, nsel); S = zeros(nst, nst, nsel)
                X = zeros(nst, nsel); XA = zeros(nst, nsel); SP = zeros(nst, nsel); YF = zeros(nb, nsel)
                lay = r.layout
            end
            XA[:, k], SP[:, k] = prior_for(sh, r.fm, r.layout, Y[:, s])
            isif = last(r.layout.idx_sif)
            A[:, :, k] = r.stats.averaging_kernel
            S[:, :, k] = r.stats.S_posterior
            X[:, k] = r.x
            YF[:, k] = r.y_fit
            res["trA"][k] = r.stats.ak_trace
            res["dof"][k] = r.stats.dof
            res["rchi2"][k] = r.stats.reduced_chi2
            res["sif_rerun"][k] = r.x[isif]
            res["sif678_rerun"][k] = w678 * r.x[isif]
            res["sif678_stored"][k] = sif_stored_all[s]
            res["sif_stored"][k] = sif_stored_all[s] / w678
            res["sif678_true"][k] = sif_true_all[s]
            res["sigma678"][k] = abs(w678) * sqrt(max(r.stats.S_posterior[isif, isif], 0.0))
            conv[k] = r.stats.converged; status[k] = r.stats.status; conv_st[k] = conv_stored_all[s]
            println(@sprintf("  [%s] %-13s sample %4d  AOD=%.3f VZA=%4.1f  conv=%d (stored %d)  DOF_sig=%.2f  SIF678 ret %.4f (stored %.4f) true %.4f  σ %.4f",
                run, cats[k], s, aod[s], vza[s], conv[k], conv_st[k], res["trA"][k], res["sif678_rerun"][k],
                res["sif678_stored"][k], res["sif678_true"][k], res["sigma678"][k]))
        end

        outnc = joinpath(outdir, "ak_$(run)_gchpio.nc")
        isfile(outnc) && rm(outnc)
        NCDataset(outnc, "c") do ds
            defDim(ds, "sel", nsel); defDim(ds, "row", nst); defDim(ds, "col", nst)
            defDim(ds, "state", nst); defDim(ds, "band", nb)
            defVar(ds, "A", Float64, ("col", "row", "sel"))[:, :, :] = permutedims(A, (2, 1, 3))
            defVar(ds, "S_post", Float64, ("col", "row", "sel"))[:, :, :] = permutedims(S, (2, 1, 3))
            defVar(ds, "x_hat", Float64, ("state", "sel"))[:, :] = X
            defVar(ds, "x_prior", Float64, ("state", "sel"))[:, :] = XA
            defVar(ds, "sigma_prior", Float64, ("state", "sel"))[:, :] = SP
            defVar(ds, "wavelength", Float64, ("band",))[:] = λ
            defVar(ds, "y_obs", Float64, ("band", "sel"))[:, :] = Y[:, picks]
            defVar(ds, "y_fit", Float64, ("band", "sel"))[:, :] = YF
            defVar(ds, "dof_signal", Float64, ("sel",))[:] = res["trA"]
            defVar(ds, "dof_residual", Float64, ("sel",))[:] = res["dof"]
            defVar(ds, "reduced_chi2", Float64, ("sel",))[:] = res["rchi2"]
            defVar(ds, "converged", Int8, ("sel",))[:] = conv
            defVar(ds, "converged_stored", Int8, ("sel",))[:] = conv_st
            defVar(ds, "status", Int16, ("sel",))[:] = status
            defVar(ds, "sif_ev1_rerun", Float64, ("sel",))[:] = res["sif_rerun"]
            defVar(ds, "sif_ev1_stored", Float64, ("sel",))[:] = res["sif_stored"]
            defVar(ds, "sif678_ret", Float64, ("sel",))[:] = res["sif678_rerun"]
            defVar(ds, "sif678_true", Float64, ("sel",))[:] = res["sif678_true"]
            defVar(ds, "sigma678", Float64, ("sel",))[:] = res["sigma678"]
            defVar(ds, "scan", Int32, ("sel",))[:] = picks
            defVar(ds, "pixel", Int32, ("sel",))[:] = zeros(Int32, nsel)
            defVar(ds, "category_id", Int8, ("sel",))[:] = [findfirst(==(c), cat_names) - 1 for c in cats]
            defVar(ds, "chl", Float64, ("sel",))[:] = aod[picks]
            defVar(ds, "nflh", Float64, ("sel",))[:] = fill(NaN, nsel)
            defVar(ds, "vza", Float64, ("sel",))[:] = vza[picks]
            defVar(ds, "sza", Float64, ("sel",))[:] = sza[picks]
            defVar(ds, "surface_albedo", Float64, ("sel",))[:] = alb[picks]
            ds.attrib["state_names_csv"] = join(svd_state_names(lay), ",")
            ds.attrib["category_names_csv"] = join(cat_names, ",")
            ds.attrib["meta_label"] = "AOD678"
            ds.attrib["run"] = run
            ds.attrib["config_file"] = cfgpath
            ds.attrib["granule"] = "gchpio_n1000"
            ds.attrib["input_ensemble"] = RT_NC
            ds.attrib["note"] = "A = S_post*H_obs; A[sel,i,j] = d xhat_i / d x_true_j; scan = 1-based sample index; chl holds AOD678"
        end
        println("wrote $outnc")
    end
end

main()
