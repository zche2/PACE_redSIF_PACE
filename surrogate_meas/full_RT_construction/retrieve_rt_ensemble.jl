# SVD retrieval on completed RT-ensemble spectra.
#
#   julia --project=. surrogate_meas/full_RT_construction/retrieve_rt_ensemble.jl
#   julia --project=. surrogate_meas/full_RT_construction/retrieve_rt_ensemble.jl \
#       surrogate_meas/configs/svd_nPC15_npoly3.toml
#
# n_pc and n_legendre (npoly) are read from [fit.svd] of the TOML (filename
# nPC*/npoly* is only a fallback). Output files are tagged so npoly3 and
# npoly5 do not overwrite each other.
#
# Reads a snapshot of rt_toa_ensemble.nc (the generator may still be writing),
# subsets to the config window, and fits each finished sample at that
# sample's SZA. Plots retrieved TOA and water-leaving SIF next to truth.

ENV["GKSwstype"] = "100"

using TOML
using NCDatasets
using Plots
using Statistics
using LinearAlgebra
using Dates

include(joinpath(@__DIR__, "..", "build_single_meas.jl"))

const DEFAULT_SVD_TOML = joinpath(@__DIR__, "..", "configs", "svd_nPC15_npoly5.toml")
const RT_NC = get(ENV, "RT_NC", joinpath(@__DIR__, "output_test_realRT_noAerosol", "rt_toa_ensemble.nc"))
const OUT_DIR = joinpath(@__DIR__, "output_test_realRT_noAerosol")
const SOLAR_L1B = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/sample_granule_20240830T131442_new_chl.nc"

function resolve_config_path()
    raw = !isempty(ARGS) ? ARGS[1] : get(ENV, "SVD_TOML", DEFAULT_SVD_TOML)
    path = isabspath(raw) ? raw : abspath(raw)
    isfile(path) || error("SVD config not found: $path")
    return path
end

"""n_pc and n_legendre from [fit.svd], with nPC*/npoly* parsed from the filename if missing."""
function detect_svd_ranks(cfg, config_path)
    svd = get(get(cfg, "fit", Dict()), "svd", Dict())
    n_pc = Int(get(svd, "n_pc", 0))
    n_leg = Int(get(svd, "n_legendre", 0))
    m = match(r"nPC(\d+).*npoly(\d+)"i, basename(config_path))
    if m !== nothing
        fn_pc = parse(Int, m.captures[1])
        fn_leg = parse(Int, m.captures[2])
        if n_pc == 0
            n_pc = fn_pc
        elseif n_pc != fn_pc
            @warn "n_pc in TOML disagrees with filename; using TOML" toml=n_pc filename=fn_pc
        end
        if n_leg == 0
            n_leg = fn_leg
        elseif n_leg != fn_leg
            @warn "n_legendre in TOML disagrees with filename; using TOML" toml=n_leg filename=fn_leg
        end
    end
    (n_pc > 0 && n_leg > 0) || error("Could not detect n_pc / n_legendre from $(basename(config_path))")
    return n_pc, n_leg, "nPC$(n_pc)_npoly$(n_leg)"
end

function snapshot_nc(src)
    snap = joinpath(OUT_DIR, "rt_toa_ensemble_snap.nc")
    cp(src, snap; force=true)
    return snap
end

function load_completed(path, λ_min, λ_max)
    ds = NCDataset(path)
    n = Int(get(ds.attrib, "n_completed", 0))
    n > 0 || error("No completed samples in $path")
    λ_all = Float64.(ds["wavelength"][:])
    idx = findall(λ_min .< λ_all .< λ_max)
    isempty(idx) && error("No OCI bands in ($λ_min, $λ_max) nm")
    ens = (
        λ = λ_all[idx],
        R_noisy = Float64.(ds["radiance_noisy"][idx, 1:n]),
        R_clean = Float64.(ds["radiance_clean"][idx, 1:n]),
        SIF = Float64.(ds["sif_waterleaving"][idx, 1:n]),
        sif_678 = Float64.(ds["sif_678"][1:n]),
        sza = Float64.(ds["sza"][1:n]),
        vza = Float64.(ds["vza"][1:n]),
        idx_lib = Int.(ds["sif_library_index"][1:n]),
    )
    close(ds)
    return ens
end

function pick_examples(sza, sif_678; n_show=4)
    n = length(sza)
    n_show = min(n_show, n)
    order = sortperm(sza)
    picks = unique([order[1], order[cld(n, 2)], order[end], argmax(sif_678)])
    return picks[1:min(n_show, length(picks))]
end

function main()
    svd_toml = resolve_config_path()
    cfg = TOML.parsefile(svd_toml)
    n_pc, n_leg, tag = detect_svd_ranks(cfg, svd_toml)
    λ_min = Float64(get(get(cfg, "spectral", Dict()), "lambda_min_nm", 640.0))
    λ_max = Float64(get(get(cfg, "spectral", Dict()), "lambda_max_nm", 756.0))
    println("Config $(basename(svd_toml)): $tag")

    snap = snapshot_nc(RT_NC)
    ens = load_completed(snap, λ_min, λ_max)
    n = length(ens.sza)
    λ = ens.λ
    println("Retrieving $n completed samples on $(length(λ)) bands in ($λ_min, $λ_max) nm")

    sh = prepare_svd_retrieval_setup(cfg, λ)
    solar = load_l1b_solar_on_bands(SOLAR_L1B, λ)
    i678 = argmin(abs.(λ .- 678.0))
    (sh.n_pc == n_pc && sh.n_leg == n_leg) ||
        error("Setup ranks ($(sh.n_pc), $(sh.n_leg)) do not match detected $tag")
    println("  n_pc=$(sh.n_pc)  n_legendre=$(sh.n_leg)  n_state=$(sh.layout.n_state)  solar=$(basename(solar.pace_path))")

    R_fit = fill(NaN, length(λ), n)
    sif_wl = fill(NaN, length(λ), n)
    resid = fill(NaN, length(λ), n)
    sif_678_ret = fill(NaN, n)
    status = fill(Int16(-1), n)

    t0 = time()
    for i in 1:n
        y = ens.R_noisy[:, i]
        try
            ret = run_pseudo_svd_retrieval(y, λ, sh, solar; sza_deg=ens.sza[i])
            R_fit[:, i] .= ret.y_fit
            wl = vec(sh.sif_basis * ret.sif_coeff)
            sif_wl[:, i] .= wl
            resid[:, i] .= ret.resid
            sif_678_ret[i] = wl[i678]
            status[i] = Int16(ret.stats.status)
        catch e
            status[i] = Int16(4)
            @warn "Retrieval failed" sample=i exception=(e, catch_backtrace())
        end
        if i == 1 || i == n || i % 5 == 0
            println("  $i / $n  ($(round(time()-t0, digits=1)) s)")
        end
    end
    n_ok = count(==(Int16(1)), status)
    println("Converged (status=1): $n_ok / $n")

    out_nc = joinpath(OUT_DIR, "rt_retrieval_$(tag).nc")
    isfile(out_nc) && rm(out_nc)
    ds = NCDataset(out_nc, "c")
    defDim(ds, "band", length(λ))
    defDim(ds, "sample", n)
    defVar(ds, "wavelength", λ, ("band",); attrib=Dict("units"=>"nm"))
    defVar(ds, "R_obs", Float32.(ens.R_noisy), ("band", "sample"))
    defVar(ds, "R_clean", Float32.(ens.R_clean), ("band", "sample"))
    defVar(ds, "R_fit", Float32.(R_fit), ("band", "sample"))
    defVar(ds, "SIF_true", Float32.(ens.SIF), ("band", "sample");
           attrib=Dict("long_name"=>"true water-leaving SIF", "units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "sif_wl_ret", Float32.(sif_wl), ("band", "sample");
           attrib=Dict("long_name"=>"retrieved water-leaving SIF", "units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "sif_678_true", Float32.(ens.sif_678), ("sample",))
    defVar(ds, "sif_678_ret", Float32.(sif_678_ret), ("sample",))
    defVar(ds, "sza", Float32.(ens.sza), ("sample",))
    defVar(ds, "status", status, ("sample",))
    ds.attrib["svd_config"] = svd_toml
    ds.attrib["rank_tag"] = tag
    ds.attrib["n_pc"] = sh.n_pc
    ds.attrib["n_legendre"] = sh.n_leg
    ds.attrib["created"] = string(Dates.now())
    close(ds)
    println("Wrote $out_nc")

    picks = pick_examples(ens.sza, ens.sif_678)
    plots = Any[]
    for i in picks
        ok = status[i] == 1
        status_tag = ok ? "" : "  (status=$(status[i]))"
        p_toa = plot(λ, ens.R_noisy[:, i]; label="true TOA (noisy)", color=:gray, lw=1.2)
        plot!(p_toa, λ, ens.R_clean[:, i]; label="true TOA (clean)", color=:black, lw=1.0, ls=:dash)
        plot!(p_toa, λ, R_fit[:, i]; label="retrieved", color=:crimson, lw=1.6)
        xlabel!(p_toa, "Wavelength (nm)")
        ylabel!(p_toa, "Radiance")
        title!(p_toa, "$tag  TOA  SZA=$(round(ens.sza[i], digits=0))° VZA=$(round(ens.vza[i], digits=0))°$status_tag")

        p_sif = plot(λ, ens.SIF[:, i]; label="true water-leaving", color=:black, lw=1.4)
        plot!(p_sif, λ, sif_wl[:, i]; label="retrieved", color=:darkorange, lw=1.6)
        vline!(p_sif, [678]; label="", color=:gray, ls=:dot, lw=1)
        xlabel!(p_sif, "Wavelength (nm)")
        ylabel!(p_sif, "SIF")
        title!(p_sif, "SIF  true@678=$(round(ens.sif_678[i], digits=3))  ret=$(round(sif_678_ret[i], digits=3))")
        push!(plots, p_toa, p_sif)
    end
    fig = plot(plots...; layout=(length(picks), 2), size=(1100, 280 * length(picks)),
               left_margin=4Plots.mm, bottom_margin=3Plots.mm)
    out_png = joinpath(OUT_DIR, "rt_retrieval_toa_sif_$(tag).png")
    savefig(fig, out_png)
    println("Wrote $out_png  ($(filesize(out_png)) bytes)")

    ok = findall(==(Int16(1)), status)
    if !isempty(ok)
        t = ens.sif_678[ok]
        r = sif_678_ret[ok]
        β = hcat(ones(length(t)), t) \ r
        r_hat = β[1] .+ β[2] .* t
        r2 = 1 - sum((r .- r_hat) .^ 2) / sum((r .- mean(r)) .^ 2)
        bias = mean(r .- t)
        rmse_s = sqrt(mean((r .- t) .^ 2))
        lim = extrema(vcat(t, r))
        pad = 0.08 * (lim[2] - lim[1] + eps())
        lim = (lim[1] - pad, lim[2] + pad)
        p = scatter(t, r; ms=5, alpha=0.45, label="samples (n=$(length(ok)))",
                    xlabel="True water-leaving SIF @ 678 nm",
                    ylabel="Retrieved",
                    title="$tag  SIF@678  bias=$(round(bias, digits=3))  RMSE=$(round(rmse_s, digits=3))  R²=$(round(r2, digits=3))",
                    legend=:topleft, size=(640, 560))
        plot!(p, [lim[1], lim[2]], [lim[1], lim[2]]; color=:black, ls=:dash, label="1:1")
        xx = range(lim[1], lim[2]; length=50)
        plot!(p, xx, β[1] .+ β[2] .* xx; color=:crimson, lw=2,
              label="fit: y=$(round(β[1], digits=3))+$(round(β[2], digits=3))x")
        xlims!(p, lim)
        ylims!(p, lim)
        println("$tag SIF@678: bias=$(round(bias, digits=4))  RMSE=$(round(rmse_s, digits=4))  R²=$(round(r2, digits=3))  slope=$(round(β[2], digits=3))  intercept=$(round(β[1], digits=3))")
        scatter_png = joinpath(OUT_DIR, "rt_retrieval_sif678_$(tag).png")
        savefig(p, scatter_png)
        println("Wrote $scatter_png")
    end
end

main()
