# SVD SIF retrieval on the 5x5 matchup boxes, from the PACE TOA radiance and SZA written by
# extract_matchup_boxes.py, with the same forward model and priors as the swath pipeline
# (build_single_meas.jl / svd_helpers.jl):
#   y = E0·cos(SZA)/π/esd · T_up^α · (Legendre ρ)/π + T_up · SIF_basis·c
# E0 = L1B red_solar_irradiance (red_bands_640_780.csv), esd = L1B earth_sun_distance_correction of the
# granule (matchup_summary.csv). Retrieved SIF is water-leaving, at the band nearest 678.2 nm (as the
# pipeline's sif_radiance_678nm); its 1σ comes from the posterior covariance, its AK from A = S_post·H.
#
#   cd /home/zhe2/FraLab/PACE_redSIF_PACE
#   julia --project=. -t 16 <this dir>/retrieve_matchup_sif.jl [config.toml]
# Default config: surrogate_meas/configs/svd_nPC15_npoly3_O2A.toml. PACE_REPO overrides the repository
# whose retrieval code and configs are used (default: the main checkout).
# → matchup_sif_<config>.csv (one row per box pixel) and matchup_sif_<config>_summary.csv (per matchup:
#   centre pixel and mean / std over valid, converged pixels), next to this script.

using CSV
using DataFrames
using Printf
using Statistics
using TOML

const REPO = get(ENV, "PACE_REPO", "/home/zhe2/FraLab/PACE_redSIF_PACE")
include(joinpath(REPO, "surrogate_meas", "build_single_meas.jl"))

const HERE = @__DIR__
const CONFIG = isempty(ARGS) ? joinpath(REPO, "surrogate_meas", "configs", "svd_nPC15_npoly3_O2A.toml") :
               (isabspath(ARGS[1]) ? ARGS[1] : joinpath(REPO, "surrogate_meas", "configs", ARGS[1]))
const SIF_REF_NM = 678.2

function load_inputs()
    box = CSV.read(joinpath(HERE, "matchup_boxes_5x5.csv"), DataFrame)
    bands = CSV.read(joinpath(HERE, "red_bands_640_780.csv"), DataFrame)
    summ = CSV.read(joinpath(HERE, "matchup_summary.csv"), DataFrame)
    esd = Dict(r.matchup_id => r.earth_sun_distance_correction_L1B for r in eachrow(summ))
    return box, bands, summ, esd
end

function main()
    cfg = TOML.parsefile(CONFIG)
    tag = replace(splitext(basename(CONFIG))[1], "svd_" => "")
    box, bands, summ, esd = load_inputs()

    sp = get(cfg, "spectral", Dict())
    λ_all = Float64.(bands.wavelength_nm)
    win = findall(Float64(get(sp, "lambda_min_nm", 640.0)) .< λ_all .< Float64(get(sp, "lambda_max_nm", 756.0)))
    λ = λ_all[win]
    E0 = Float64.(bands.E0_W_m2_um[win])
    cols = [@sprintf("Rtoa_%.2f", w) for w in λ]
    all(c -> c in names(box), cols) || error("Rtoa columns for the config window missing in matchup_boxes_5x5.csv")

    sh = prepare_svd_retrieval_setup(cfg, λ)
    i678 = argmin(abs.(λ .- SIF_REF_NM))
    b678 = vec(sh.sif_basis[i678, :])
    @printf("%s: %d bands %.2f–%.2f nm, %d PCs, Legendre %d, SIF band %.2f nm; %d box pixels, %d threads\n",
            basename(CONFIG), length(λ), λ[1], λ[end], sh.n_pc, sh.n_leg, λ[i678], nrow(box), Threads.nthreads())

    Y = Matrix{Float64}(box[:, cols])                  # pixels × bands
    n = nrow(box)
    fitme = [all(isfinite, view(Y, i, :)) && isfinite(box.sza[i]) for i in 1:n]
    out = Dict(k => fill(NaN, n) for k in ("sif_678", "sif_678_sigma", "sif_678_ak", "sif_toa_678", "alpha",
                                             "reduced_chi2", "rmse", "dof", "ak_trace", "rel_resid_rms"))
    status = fill(Int16(-1), n)                         # -1 = not fitted (missing Rtoa / SZA)
    n_steps = fill(-1, n)
    n_leg_coeff = length(sh.layout.idx_legendre)        # degrees 0..n_legendre
    leg = fill(NaN, n, n_leg_coeff)

    todo = findall(fitme)
    nmax = parse(Int, get(ENV, "MAX_PIXELS", "0"))         # > 0: fit only the first nmax pixels (testing)
    nmax > 0 && (todo = todo[1:min(nmax, end)])
    first_err = String[]
    err_lock = ReentrantLock()
    t0 = time()
    Threads.@threads for i in todo
        try
            solar = (E=E0, esd=Float64(esd[box.matchup_id[i]]), pace_path=String(box.granule[i]))
            y = Y[i, :]
            r = run_pseudo_svd_retrieval(y, λ, sh, solar; sza_deg=Float64(box.sza[i]))
            L = r.layout
            S = r.stats.S_posterior
            A = r.stats.averaging_kernel
            is_ = collect(L.idx_sif)
            out["sif_678"][i] = dot(b678, r.sif_coeff)
            out["sif_678_sigma"][i] = sqrt(max(dot(b678, S[is_, is_] * b678), 0.0))
            out["sif_678_ak"][i] = length(is_) == 1 ? A[is_[1], is_[1]] : NaN
            out["sif_toa_678"][i] = r.sif_toa_fit[i678]
            out["alpha"][i] = alpha_coeff_from_raw(r.x[first(L.idx_alpha)])
            out["reduced_chi2"][i] = r.stats.reduced_chi2
            out["rmse"][i] = r.stats.rmse
            out["dof"][i] = r.stats.dof
            out["ak_trace"][i] = r.stats.ak_trace
            out["rel_resid_rms"][i] = sqrt(mean((r.resid ./ y) .^ 2))
            leg[i, :] .= r.x[L.idx_legendre]
            status[i] = Int16(r.stats.status)
            n_steps[i] = r.stats.n_steps
        catch e
            status[i] = Int16(4)
            lock(err_lock) do
                isempty(first_err) && push!(first_err, sprint(showerror, e, catch_backtrace()))
            end
        end
    end
    isempty(first_err) || @warn "$(count(==(4), status)) retrievals failed; first error:\n$(first(first_err))"
    @printf("fitted %d pixels in %.0f s; status counts: %s\n", length(todo), time() - t0,
            join(["$s => $(count(==(s), status))" for s in sort(unique(status))], ", "))

    keep = [:matchup_id, :granule, :scan, :pixel, :d_scan, :d_pixel, :is_centre, :lat, :lon, :pixel_time_utc,
            :dt_minutes, :sza, :vza, :saa, :vaa, :watermask, :l2_flags, :flag_ok, :geom_ok, :valid,
            :rtoa_complete, :is_dark, :nflh, :nflh_678, :nflh_680, :chlor_a, :aot_865, :angstrom]
    res = select(box, intersect(keep, Symbol.(names(box))))
    res.status = status
    res.converged = status .== 1
    res.n_steps = n_steps
    for k in ("sif_678", "sif_678_sigma", "sif_678_ak", "sif_toa_678", "alpha", "reduced_chi2", "rmse", "dof",
              "ak_trace", "rel_resid_rms")
        res[!, k] = out[k]
    end
    for j in 1:n_leg_coeff
        res[!, "legendre_$(j - 1)"] = leg[:, j]
    end
    pix_out = joinpath(HERE, "matchup_sif_$(tag).csv")
    CSV.write(pix_out, res)

    # per matchup: centre pixel, and mean / std over valid (L2 flags + geometry) and converged pixels
    rows = NamedTuple[]
    for g in groupby(res, :matchup_id)
        c = g[g.is_centre, :]
        u = g[g.valid .& g.converged, :]
        m(v) = isempty(v) ? NaN : mean(v)
        s(v) = length(v) < 2 ? NaN : std(v)
        push!(rows, (matchup_id=g.matchup_id[1], granule=g.granule[1], n_box=nrow(g), n_valid=count(g.valid),
                     n_valid_converged=nrow(u),
                     sif_678_centre=nrow(c) == 1 ? c.sif_678[1] : NaN,
                     sif_678_sigma_centre=nrow(c) == 1 ? c.sif_678_sigma[1] : NaN,
                     converged_centre=nrow(c) == 1 && c.converged[1],
                     sif_678_mean=m(u.sif_678), sif_678_std=s(u.sif_678), sif_678_sigma_mean=m(u.sif_678_sigma),
                     reduced_chi2_mean=m(u.reduced_chi2), nflh_mean=m(u.nflh), nflh_678_mean=m(u.nflh_678),
                     nflh_680_mean=m(u.nflh_680), chlor_a_mean=m(u.chlor_a)))
    end
    sm = DataFrame(rows)
    sm = leftjoin(sm, select(summ, :matchup_id, :data_type, :FIELD_lat, :FIELD_lon, :passes_protocol), on=:matchup_id)
    sort!(sm, :matchup_id)
    sum_out = joinpath(HERE, "matchup_sif_$(tag)_summary.csv")
    CSV.write(sum_out, sm)

    ok = res.valid .& res.converged
    med(v) = (f = filter(isfinite, v); isempty(f) ? NaN : median(f))
    @printf("valid pixels %d, converged %d (%.1f %%); SIF678 median %.3f, median σ %.3f, median AK %.3f, median reduced χ² %.2f [W m-2 sr-1 um-1]\n",
            count(res.valid), count(ok), 100 * count(ok) / max(count(res.valid), 1),
            med(res.sif_678[ok]), med(res.sif_678_sigma[ok]), med(res.sif_678_ak[ok]), med(res.reduced_chi2[ok]))
    println("wrote $pix_out\n      $sum_out")
end

main()
