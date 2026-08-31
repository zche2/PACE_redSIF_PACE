#!/usr/bin/env julia
# Find continuum bands where total transmittance = solar × atmospheric > threshold.
using NCDatasets
using Statistics
using DelimitedFiles
using Interpolations

const TRANS_NC = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc"
const SOLAR_FILE = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/solar_merged_20240731_600_33300_100.out"
const PROFILE = 5511
const T_THRESH = 0.995
const λ_LO = 600.4
const λ_HI = 894.6

"""Load solar Fraunhofer transmittance onto `λ_target` [nm] (same as MWE solar file)."""
function load_solar_on_bands(solar_path::AbstractString, λ_target::AbstractVector{<:Real}; header_lines::Int=3)
    tbl = readdlm(solar_path, Float64, skipstart=header_lines)
    ν_src = tbl[:, 1]
    T_src = tbl[:, 2]
    λ_src = 1e7 ./ ν_src
    p = sortperm(λ_src)
    λ_s = λ_src[p]
    T_s = T_src[p]
    λt = collect(Float64.(λ_target))
    (minimum(λt) >= minimum(λ_s) && maximum(λt) <= maximum(λ_s)) ||
        error("Target λ outside solar file range [$(minimum(λ_s)), $(maximum(λ_s))]")
    itp = LinearInterpolation(λ_s, T_s, extrapolation_bc=Throw())
    return itp.(λt)
end

# ── atmospheric transmittance (library, one-way, band-convolved) ───────────────
ds = Dataset(TRANS_NC)
λ_all = Float64.(ds["band"][:])
Traw = Float64.(ds["transmittance"][:, :])
close(ds)

if size(Traw, 2) == length(λ_all)
    T_atm = Traw
elseif size(Traw, 1) == length(λ_all)
    T_atm = permutedims(Traw)
else
    error("Unexpected transmittance size $(size(Traw)) vs nband=$(length(λ_all))")
end

band_mask = (λ_all .>= λ_LO) .& (λ_all .<= λ_HI)
λ = λ_all[band_mask]
T_atm = T_atm[:, band_mask]

# ── solar Fraunhofer transmittance on the same bands ──────────────────────────
T_solar = load_solar_on_bands(SOLAR_FILE, λ)

# total = solar × atmospheric
T_tot = T_atm .* T_solar'   # (profile, band) .* (1, band) broadcast → need careful
# T_atm is (n_prof, n_band); T_solar is (n_band,)
T_tot = T_atm .* reshape(T_solar, 1, :)

ip = clamp(PROFILE, 1, size(T_tot, 1))
t_prof = T_tot[ip, :]
t_med = vec(median(T_tot; dims=1))
t_min = vec(minimum(T_tot; dims=1))

t_atm_prof = T_atm[ip, :]
t_atm_med = vec(median(T_atm; dims=1))

println("Atm file:  $TRANS_NC")
println("Solar file: $SOLAR_FILE")
println("Definition: T_total = T_solar × T_atm  (threshold = $T_THRESH)")
println("Window: [$λ_LO, $λ_HI] nm → $(length(λ)) bands")
println("Profiles: $(size(T_tot, 1)); using profile $ip")
println("Solar     T range: $(extrema(T_solar))")
println("Atm med   T range: $(extrema(t_atm_med))")
println("Total med T range: $(extrema(t_med))")
println("Profile $ip total T range: $(extrema(t_prof))")

function report(name, tvec; t_atm_ref=t_atm_med)
    ind = findall(tvec .> T_THRESH)
    println("\n=== $name: T_solar×T_atm > $T_THRESH → $(length(ind)) / $(length(tvec)) bands ===")
    for i in ind
        println("  $(round(λ[i]; digits=3)) nm    T_tot=$(round(tvec[i]; digits=6))" *
                "   (T_sol=$(round(T_solar[i]; digits=5)), T_atm=$(round(t_atm_ref[i]; digits=5)))")
    end
    return ind
end

ind_p = report("profile $ip", t_prof; t_atm_ref=t_atm_prof)
ind_m = report("median across profiles", t_med; t_atm_ref=t_atm_med)

function clusters(ind)
    isempty(ind) && (println("\n(no clusters)"); return)
    println("\nClusters (contiguous high-T_tot bands for median):")
    start = ind[1]
    prev = ind[1]
    for i in ind[2:end]
        if i != prev + 1
            println("  $(round(λ[start]; digits=2))–$(round(λ[prev]; digits=2)) nm  ($(prev-start+1) bands)")
            start = i
        end
        prev = i
    end
    println("  $(round(λ[start]; digits=2))–$(round(λ[prev]; digits=2)) nm  ($(prev-start+1) bands)")
end
clusters(ind_m)

const BASELINE_λ_REF = [
    607.99, 610.36, 612.73, 615.14, 617.6, 620.06, 622.53,
    669.52, 670.76, 671.99, 673.24, 674.51, 675.73, 676.96, 678.21, 679.45,
    754.3, 779.33, 867.11, 869.61, 872.13,
]
println("\n=== Original baseline refs vs median T_solar×T_atm ===")
for r in BASELINE_λ_REF
    i = argmin(abs.(λ .- r))
    pass = t_med[i] > T_THRESH ? "PASS" : "fail"
    println("  ref=$(r) → nearest=$(round(λ[i]; digits=3))" *
            "  T_tot_med=$(round(t_med[i]; digits=5))" *
            "  T_tot_prof=$(round(t_prof[i]; digits=5))" *
            "  T_sol=$(round(T_solar[i]; digits=5))" *
            "  [$pass]")
end

# printable Julia vector of continuum bands (median)
println("\n# Suggested BASELINE_λ_REF (median T_solar×T_atm > $T_THRESH):")
println("const BASELINE_λ_REF = [")
for i in ind_m
    println("    $(round(λ[i]; digits=3)),")
end
println("]")
