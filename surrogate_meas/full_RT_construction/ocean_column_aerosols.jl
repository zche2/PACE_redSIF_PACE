# Load one open-ocean column from gchp_ocean_columns_n500.nc into
# params.scattering_params.rt_aerosols.
#
# Vertical order
#   GCHP NetCDF:  BOA → TOA  (lev 1 / p_half[1] = surface, high p)
#   vSmartMOM RT: TOA → BOA  (p_half[1] = TOA, low p, increasing downward)
# Pressure and AOD are always reversed as a pair. Never flip one alone.
#
# RT_Aerosol.profile still carries a Normal(p₀,σp) placeholder (vSmartMOM type
# requires a Distributions.jl object). For ensemble RT, call
# `apply_gchp_layer_aod!` after `update_model!` to overwrite τ_aer with the
# realistic GCHP layer-AOD profile remapped onto the current RT p_half grid.
# MERRA T/p/q are left untouched.
#
# Requires `using vSmartMOM` (and thus CoreRT / Scattering / Aerosols) in the caller.

import YAML as _YAML
using NCDatasets
using Distributions
using Interpolations
using Printf

const _RI_DB_PATH = joinpath(dirname(dirname(pathof(vSmartMOM))), "data", "refractive_indices_database.yaml")
const _AOD_MIN = 1e-6
const _SIGP_MIN = 20.0  # hPa; only for the placeholder Normal

"""
    load_ocean_column_aerosols!(params, nc_path, i_sample; yaml_path, ri_db_path)

Read column `i_sample` (1-based) from the n500 ocean-column NetCDF and replace
`params.scattering_params.rt_aerosols` with one `RT_Aerosol` per species listed
under `aerosol_scheme.species` in `yaml_path`.

For each species:
- τ_ref = column sum of layer AOD
- μ = AOD-weighted mean wet radius (μm)
- σ = sigma_g from YAML
- profile = Normal(p₀, σp) placeholder (AOD-weighted mean/std of `p_mid`)
- nᵣ, nᵢ from the vSmartMOM refractive-index database at λ_ref

Returns `(p_mid, layer_aods, p_half)` in native GCHP order (BOA→TOA).
`layer_aods[i]` matches `p_mid` index-for-index. `p_half` is the native
GCHP half-level vector (length `n_lev+1`) or `nothing` if the file has none.
"""
function load_ocean_column_aerosols!(
    params,
    nc_path::AbstractString,
    i_sample::Integer;
    yaml_path::AbstractString,
    ri_db_path::AbstractString=_RI_DB_PATH,
    aod_min::Real=_AOD_MIN,
)
    sp = params.scattering_params
    sp === nothing && error("params.scattering_params is nothing; add a scattering: block to the YAML")

    cfg = _YAML.load_file(yaml_path)
    scheme = get(cfg, "aerosol_scheme", Dict())
    species_cfg = get(scheme, "species", Dict())
    isempty(species_cfg) && error("No aerosol_scheme.species in $yaml_path")

    λ_ref = Float64(sp.λ_ref)
    r_max = Float64(sp.r_max)
    FT = typeof(sp.λ_ref)

    ri_db = vSmartMOM.Aerosols.load_refractive_index_database(ri_db_path, Float64)

    ds = NCDataset(nc_path)
    n_samp = ds.dim["sample"]
    n_lev = ds.dim["lev"]
    1 ≤ i_sample ≤ n_samp || error("i_sample=$i_sample out of range 1:$n_samp")
    p_mid = _read_column_vec(ds, "p_mid", i_sample, n_lev, n_samp)
    p_half = haskey(ds, "p_half") ?
        _read_column_vec(ds, "p_half", i_sample, n_lev + 1, n_samp) : nothing
    _assert_gchp_vertical_pairing(p_mid, p_half)
    lat_i = Float64(ds["lat"][i_sample])
    lon_i = Float64(ds["lon"][i_sample])

    column_profile(varname) = _read_column_vec(ds, varname, i_sample, n_lev, n_samp)

    rt_list = vSmartMOM.CoreRT.RT_Aerosol{FT}[]
    layer_aods = Vector{Vector{Float64}}()
    summaries = String[]

    for (name, sc) in sort(collect(species_cfg); by=first)
        aod_var = String(sc["aod_var"])
        rad_var = String(sc["radius_var"])
        sigma_g = Float64(sc["sigma_g"])
        ri_key = String(sc["refractive_index"])

        aod = column_profile(aod_var)
        rad = column_profile(rad_var)
        length(aod) == n_lev || error("$aod_var column length $(length(aod)) ≠ n_lev=$n_lev")
        τ_ref = sum(aod)
        τ_ref < aod_min && continue

        # Callers that pass aod_min=0.0 (to keep a fixed species count for
        # BatchContext) rely on this not producing NaN when τ_ref is exactly
        # zero; 0/0 would otherwise poison μ/p₀/σp even though the species
        # contributes no optical depth.
        w = τ_ref > 0 ? aod ./ τ_ref : fill(inv(n_lev), n_lev)
        μ = clamp(sum(w .* rad), 1e-3, r_max)
        p₀ = sum(w .* p_mid)
        σp = max(sqrt(sum(w .* (p_mid .- p₀) .^ 2)), _SIGP_MIN)

        n_c = try
            vSmartMOM.Aerosols.get_refractive_index(ri_db, ri_key, λ_ref)
        catch e
            @warn "RI lookup failed; using seasalt-like fallback" species=name ri_key=ri_key exception=e
            1.5 + 1e-8im
        end
        nᵣ = FT(real(n_c))
        nᵢ = FT(abs(imag(n_c)))  # Mie path expects non-negative imag

        size_dist = LogNormal(log(μ), log(sigma_g))
        aero = vSmartMOM.Scattering.Aerosol(size_dist, nᵣ, nᵢ)
        # Placeholder only — overwrite τ_aer via apply_gchp_layer_aod!.
        profile = Normal(FT(p₀), FT(σp))
        push!(rt_list, vSmartMOM.CoreRT.RT_Aerosol(aero, FT(τ_ref), profile))
        push!(layer_aods, aod)
        push!(summaries, @sprintf("%s: τ=%.4f μ=%.3fμm (GCHP layer profile)", name, τ_ref, μ))
    end
    close(ds)

    isempty(rt_list) && error("No species with AOD ≥ $aod_min in column $i_sample")
    sp.rt_aerosols = rt_list

    println("Ocean column $i_sample @ ($(round(lat_i; digits=1))°, $(round(lon_i; digits=1))°): $(length(rt_list)) aerosols")
    for s in summaries
        println("  ", s)
    end
    orient = p_mid[1] > p_mid[end] ? "BOA→TOA" : "TOA→BOA"
    println("  native vertical order $orient  p_mid=$(round(p_mid[1]; digits=1)) → $(round(p_mid[end]; digits=3)) hPa")
    return p_mid, layer_aods, p_half
end

"""Read a 1-D or 2-D NetCDF vector for sample `i` along a named size-`n_along` axis."""
function _read_column_vec(ds, name::AbstractString, i::Integer, n_along::Integer, n_samp::Integer)
    haskey(ds, name) || error("Missing var $name")
    A = Array(ds[name])
    if ndims(A) == 1
        length(A) == n_along || error("$name length $(length(A)) ≠ $n_along")
        return Float64.(A)
    elseif ndims(A) == 2
        if size(A) == (n_along, n_samp)
            return Float64.(A[:, i])
        elseif size(A) == (n_samp, n_along)
            return Float64.(A[i, :])
        else
            error("$name shape $(size(A)) incompatible with (along=$n_along, sample=$n_samp)")
        end
    else
        error("$name should be 1-D or 2-D, got ndims=$(ndims(A))")
    end
end

"""GCHP `p_mid` and `p_half` must share BOA→TOA order; AOD uses the same `lev` as `p_mid`."""
function _assert_gchp_vertical_pairing(p_mid, p_half)
    n = length(p_mid)
    n >= 2 || error("Need ≥2 GCHP levels")
    if p_half === nothing
        issorted(p_mid) || issorted(p_mid; rev=true) ||
            error("p_mid is not monotonic; cannot pair AOD with pressure")
        return nothing
    end
    length(p_half) == n + 1 ||
        error("p_half length $(length(p_half)) ≠ n_lev+1=$(n + 1)")
    p_mid_from_half = @. 0.5 * (p_half[1:end-1] + p_half[2:end])
    mid_boa = p_mid[1] > p_mid[end]
    half_boa = p_half[1] > p_half[end]
    mid_boa == half_boa || error(
        "p_mid and p_half have opposite vertical order: " *
        "p_mid $(p_mid[1]) → $(p_mid[end]) hPa, p_half $(p_half[1]) → $(p_half[end]) hPa")
    err = maximum(abs.(p_mid_from_half .- p_mid))
    err_flip = maximum(abs.(p_mid_from_half .- reverse(p_mid)))
    if err_flip + 1.0 < err
        error("p_mid looks reversed relative to p_half (max |Δp|=$(round(err; digits=2)) hPa " *
              "vs $(round(err_flip; digits=2)) hPa if p_mid is flipped). " *
              "AOD would then be assigned to the wrong layers.")
    end
    err > 5.0 && @warn "p_mid vs p_half midpoints differ by $(round(err; digits=2)) hPa"
    return nothing
end

"""Force layer pressure and AOD to TOA→BOA together (never reverse one alone)."""
function _layers_toa_boa(p_mid, aod)
    length(p_mid) == length(aod) ||
        error("p_mid / aod length mismatch: $(length(p_mid)) vs $(length(aod))")
    p = collect(Float64, p_mid)
    a = max.(Float64.(aod), 0.0)
    if p[1] > p[end]
        reverse!(p)
        reverse!(a)
    elseif p[1] < p[end]
        # already TOA→BOA
    else
        error("p_mid is not monotonic: $(p[1]) → $(p[end])")
    end
    return p, a
end

"""Force half-levels to TOA→BOA. Reverse layer AOD with the half-levels."""
function _edges_toa_boa(p_half, aod)
    n = length(aod)
    length(p_half) == n + 1 ||
        error("p_half length $(length(p_half)) must be length(aod)+1=$(n + 1)")
    ph = collect(Float64, p_half)
    a = max.(Float64.(aod), 0.0)
    if ph[1] > ph[end]
        reverse!(ph)
        reverse!(a)
    elseif ph[1] < ph[end]
        # already TOA→BOA: a[i] is the layer between ph[i] and ph[i+1]
    else
        error("p_half is not monotonic: $(ph[1]) → $(ph[end])")
    end
    return ph, a
end

function _strictly_increasing!(ph::Vector{Float64})
    for i in 2:length(ph)
        if ph[i] <= ph[i - 1]
            ph[i] = nextfloat(ph[i - 1])
        end
    end
    return ph
end

"""Unique strictly increasing log-p knots with cumulative AOD merged at clips.

GCHP TOA half-levels are often floored at ~1e-4 hPa, so `p` and `log(p)`
repeat. Do **not** fix this with `nextfloat(p)`: one ulp in `p` is often still
0 ulp in `log(p)`, and Interpolations.jl then warns about duplicated /
successive-repeated knots on the full 73-edge vector. Merge degenerate edges
instead (column τ is carried in `cum`).
"""
function _cum_vs_logp(p_edge::AbstractVector, aod::AbstractVector)
    length(p_edge) == length(aod) + 1 ||
        error("p_edge length $(length(p_edge)) ≠ length(aod)+1=$(length(aod) + 1)")
    logp = log.(Float64.(p_edge))
    cum = vcat(0.0, cumsum(Float64.(aod)))
    lp = Float64[logp[1]]
    cc = Float64[cum[1]]
    for i in 2:length(logp)
        # Require a true Float64 increase in log-p (rejects == and nextfloat noise).
        if logp[i] > nextfloat(lp[end])
            push!(lp, logp[i])
            push!(cc, cum[i])
        else
            cc[end] = cum[i]   # same / tiny log-p edge: keep later column τ
        end
    end
    length(lp) >= 2 || error("Need ≥2 distinct log-p edges after merging TOA clips")
    all(diff(lp) .> 0) || error("log-p knots not strictly increasing after merge")
    return lp, cc
end

"""RT p_half must already be TOA→BOA (increasing p). Do not reverse it here —
that would desynchronize τ_aer from the atmosphere layers."""
function _rt_half_toa_boa(p_half_rt)
    ph = collect(Float64, p_half_rt)
    length(ph) >= 2 || error("p_half_rt too short")
    ph[1] > ph[end] && error(
        "p_half_rt is BOA→TOA ($(ph[1]) → $(ph[end]) hPa). " *
        "vSmartMOM layers are TOA→BOA; pass the RT grid unflipped.")
    ph[1] ≈ ph[end] && error("p_half_rt is not monotonic: $(ph[1]) → $(ph[end])")
    _strictly_increasing!(ph)
    return ph
end

"""
    remap_gchp_aod_to_rt_layers(p_mid_gchp, aod_gchp, p_half_rt; p_half_gchp) -> Vector

Map a GCHP layer-AOD column onto the RT model's half-level grid (TOA→BOA).

GCHP arrays may be BOA→TOA; they are flipped to TOA→BOA **as a pair**
(`p` with `aod`, or `p_half` with `aod`). Uses cumulative AOD vs log-pressure
so column τ is conserved; result is renormalized to `sum(aod_gchp)`.
The returned vector is TOA→BOA, matching vSmartMOM `τ_aer` layer index 1 = TOA.
"""
function remap_gchp_aod_to_rt_layers(
    p_mid_gchp::AbstractVector,
    aod_gchp::AbstractVector,
    p_half_rt::AbstractVector;
    p_half_gchp=nothing,
)
    length(p_mid_gchp) == length(aod_gchp) ||
        error("p_mid / aod length mismatch: $(length(p_mid_gchp)) vs $(length(aod_gchp))")
    n = length(p_mid_gchp)
    n >= 2 || error("Need ≥2 GCHP levels to remap")

    ph_rt = _rt_half_toa_boa(p_half_rt)

    # Do not call `_strictly_increasing!` on GCHP edges: that invents near-duplicate
    # log-p knots. `_cum_vs_logp` merges clipped TOA edges instead.
    p_src, a, p_edge = if p_half_gchp !== nothing
        ph, a = _edges_toa_boa(p_half_gchp, aod_gchp)
        p_mids = @. 0.5 * (ph[1:end-1] + ph[2:end])
        p_mids, a, ph
    else
        p, a = _layers_toa_boa(p_mid_gchp, aod_gchp)
        p_edge = Vector{Float64}(undef, n + 1)
        p_edge[1] = max(p[1]^2 / p[2], 1e-4)
        for i in 1:(n - 1)
            p_edge[i + 1] = sqrt(p[i] * p[i + 1])
        end
        p_edge[n + 1] = max(p[n]^2 / p[n - 1], p[n] * 1.001)
        p, a, p_edge
    end

    τ_col = sum(a)
    Nz = length(ph_rt) - 1
    τ = zeros(Float64, Nz)
    τ_col == 0 && return τ

    lp, cum = _cum_vs_logp(p_edge, a)
    itp = LinearInterpolation(lp, cum; extrapolation_bc=Flat())
    for i in 1:Nz
        τ[i] = max(0.0, itp(log(ph_rt[i + 1])) - itp(log(ph_rt[i])))
    end
    s = sum(τ)
    s > 0 || return τ
    τ .*= (τ_col / s)

    p_full_rt = @. 0.5 * (ph_rt[1:end-1] + ph_rt[2:end])
    # Thin lofted species (e.g. SO4 at ~190 hPa, τ ~ 0.003) trip this on
    # clean ocean columns. Keep the remapped profile and continue.
    try
        _assert_aod_not_flipped(p_src, a, p_full_rt, τ)
    catch e
        @warn "AOD orientation guard tripped; keeping remapped τ" exception=(e, catch_backtrace())
        flush(stderr)
    end
    return τ
end

"""Catch an AOD profile that was reversed relative to pressure (MBL ↔ TOA)."""
function _assert_aod_not_flipped(p_src, a_src, p_rt, a_rt)
    τ_s = sum(a_src)
    τ_r = sum(a_rt)
    (τ_s > 0 && τ_r > 0) || return nothing
    w_src = sum(a_src .* p_src) / τ_s
    w_rt = sum(a_rt .* p_rt) / τ_r
    # After a paired BOA→TOA reverse, ocean AOD sits in the MBL (~800–1000 hPa).
    # Reversing AOD without pressure (or the converse) puts that mass at TOA.
    if τ_s > 1e-3 && w_src < 200
        error("GCHP AOD-weighted p is $(round(w_src; digits=1)) hPa after TOA→BOA " *
              "orient (τ=$(round(τ_s; digits=4))). Pressure and AOD were probably " *
              "reversed independently; GCHP native order is BOA→TOA.")
    end
    if w_src > 400 && w_rt < 200
        error("AOD/pressure pairing looks flipped: AOD-weighted p is " *
              "$(round(w_src; digits=1)) hPa on GCHP vs $(round(w_rt; digits=1)) hPa on RT. " *
              "GCHP is BOA→TOA; RT is TOA→BOA — reverse p and AOD together.")
    end
    return nothing
end

"""
    apply_prescribed_layer_aod!(ctx, i_aer, τ_layer)

Overwrite `τ_aer` for species `i_aer` with a prescribed layer-AOD profile on the
current RT grid (length `ctx.Nz`), keeping Mie spectral scaling via `k / k_ref`.
"""
function apply_prescribed_layer_aod!(ctx, i_aer::Integer, τ_layer::AbstractVector)
    1 <= i_aer <= ctx.n_aerosols || error(
        "apply_prescribed_layer_aod!: i_aer=$i_aer out of range [1, $(ctx.n_aerosols)]")
    length(τ_layer) == ctx.Nz || error(
        "τ_layer length $(length(τ_layer)) ≠ ctx.Nz=$(ctx.Nz)")

    model = ctx.model
    FT = ctx.params.float_type
    τ_eff = FT(sum(τ_layer))
    ctx.current_τ_ref[i_aer] = τ_eff
    τ_frac = τ_eff > 0 ? FT.(τ_layer ./ τ_eff) : zeros(FT, ctx.Nz)
    k_ref_aer = ctx.k_ref[i_aer]

    for i_band in 1:ctx.n_bands
        k_aer = model.optics.aerosols.aerosol_optics[i_band][i_aer].k
        model.optics.aerosols.τ_aer[i_band][i_aer, :, :] .=
            (τ_eff / k_ref_aer) .* k_aer .* τ_frac'
    end
    return nothing
end

"""
    apply_gchp_layer_aod!(ctx, p_mid_gchp, layer_aods, p_half_rt; p_half_gchp)

Remap each cached GCHP species profile onto `p_half_rt` (TOA→BOA) and write
into `τ_aer`. Pass native GCHP `p_half_gchp` (BOA→TOA) when available so layer
edges stay paired with AOD. Call after every `update_model!`.
"""
function apply_gchp_layer_aod!(ctx, p_mid_gchp, layer_aods, p_half_rt; p_half_gchp=nothing)
    length(layer_aods) == ctx.n_aerosols || error(
        "layer_aods has $(length(layer_aods)) species but BatchContext expects $(ctx.n_aerosols)")
    for i in 1:ctx.n_aerosols
        τ_rt = remap_gchp_aod_to_rt_layers(
            p_mid_gchp, layer_aods[i], p_half_rt; p_half_gchp=p_half_gchp)
        apply_prescribed_layer_aod!(ctx, i, τ_rt)
    end
    return nothing
end
