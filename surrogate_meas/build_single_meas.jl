#!/usr/bin/env julia
# Build a single surrogate measurement spectrum from one MERRA2 atmospheric profile.
#
# Workflow:
# (A) Two-way Solar x Atmospheric Transmittance
#   1. Randomly pick a profile from the transmittance NetCDF
#   2. Plot height, pressure, temperature, and H2O mixing ratio
#   3. Prescribe a non-χ² layer reflectance distribution (surface-heavy)
#   4. Compute layer-resolved optical depth on a high-resolution spectral grid
#   5. Form weighted two-way optical depths / transmittances
#   6. Include solar Fraunhofer transmittance with two-way paths
#   7. Convolve with the PACE OCI instrument response function
#   8. Plot hi-res and low-res transmittance on the same panels
# (B) Rebuild Continuum radiance (w/ surface reflectance)
#   1. Pick a random pixel of radiance from L1B data
#   2. Pick baseline windows
#   3. Fit a Legendre polynomial (order of 7) to R_L1B at baseline bands
#   4. Rescale the mean radiance to a random number between 15 and 30
#   5. Plot fitted and rescaled radiance (and derived ρ) vs L1B
# (C) Add SIF
#   1. Calculate a unweighted one-way transmittance spectrum by adding up optical depth of each layer and convolve to lres, T1
#   2. Choose a random SIF shape from SIF spectrum library and rescale to a strength between 0 and 0.5 
#   3. Plot (a) the SIF spectrum, SIF x T1 (b) T1
# (D) Rebuild TOA radiance + white noise
#   1. multiply the continuum radiance by the solar x atm transmittance to get background radiance
#   2. add SIF radiance to the background radiance to get TOA radiance
#   3. calculate standard deviation of the white noise based on instrument SNR and radiance
#   4. add white noise assuming Gaussian distribution to the spectra
#   5. Plot clean / noisy TOA (with and without SIF)
# (E) Pseudo retrieval
#   1. Use svd pipeline to do the retrieval, use E from a random L1B file, set sza=0
#   2. Plot (a) true and fitted TOA, (b) true and fitted SIF, (c) residual vs. white noise (d) relative residual

using Random
using Statistics
using LinearAlgebra
using TOML
using NCDatasets
using Dates
# Headless GR backend — without this, savefig can write a 0-byte PNG on SSH/CI.
get!(ENV, "GKSwstype", "100")
using Plots
using JLD2
using ForwardDiff
using SparseArrays
using Interpolations

const SCRIPT_DIR = @__DIR__
const REPO_ROOT = dirname(SCRIPT_DIR)

include(joinpath(REPO_ROOT, "demo_example", "Simple_PACE_xSecFit_MWE_Functions.jl"))
using .SimplePACEXSecFitMWEFunctions

const MWEF = SimplePACEXSecFitMWEFunctions

include(joinpath(REPO_ROOT, "global_svd_fit_pipeline", "svd_retrieval", "svd_helpers.jl"))
include(joinpath(SCRIPT_DIR, "alpha_mapping.jl"))
include(joinpath(SCRIPT_DIR, "noise_degrade.jl"))

# Surrogate retrieval: α_coeff ∈ [1, 20] (global svd_helpers defaults to [1, 11]).
function make_svd_forward_model_λ(
    λ_obs::AbstractVector{Float64},
    solar_eff::AbstractVector{Float64},
    PCs_obs::Matrix{Float64},
    sif_basis_obs::Matrix{Float64};
    n_pc::Int,
    n_legendre::Int,
    log_transform::Bool,
)
    length(λ_obs) == length(solar_eff) || error("λ and solar_eff length mismatch")
    PCs = Float64.(PCs_obs[:, 1:n_pc])
    SIF = Float64.(sif_basis_obs)
    z_obs = _normalized_grid(λ_obs)
    leg_basis = _legendre_design_matrix(z_obs, n_legendre)
    layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_legendre, n_ev = size(SIF, 2))

    function fm_svd(x::AbstractVector)
        length(x) == layout.n_state || error("SVD state length $(length(x)) != $(layout.n_state)")
        c_vec = @view x[layout.idx_pc]
        alpha_coeff = alpha_coeff_from_raw(x[first(layout.idx_alpha)])
        leg_coeff = @view x[layout.idx_legendre]
        sif_coeff = @view x[layout.idx_sif]
        trans_up = log_transform ? exp.(PCs * c_vec) : 1.0 .+ PCs * c_vec
        trans_updown = exp.(alpha_coeff .* log.(max.(trans_up, eps(Float64))))
        rho_obs = leg_basis * leg_coeff
        sif_toa = trans_up .* (SIF * sif_coeff)
        return @.(solar_eff * trans_updown * rho_obs / π + sif_toa)
    end
    return fm_svd, layout
end

# ── user settings ─────────────────────────────────────────────────────────────
const CONFIG_PATH = joinpath(REPO_ROOT, "demo_example", "Simple_PACE_xSecFit_MWE_zcheVer.toml")
const TRANS_NC = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc"
const OUTPUT_DIR = joinpath(SCRIPT_DIR, "output")
const RANDOM_SEED = 42
const SZA_DEG = 30.0
const VZA_DEG = 0.0
const REFLECTANCE_SHAPE = :exponential   # :exponential or :beta (both non-χ²)
const REFLECTANCE_DECAY = 10.0            # larger => more weight near the surface
# Linear LUT interp avoids cubic overshoot → negative xsecs / transmittance > 1
const LUT_INTERPOLATION = "linear"
const SVD_RETRIEVAL_CONFIG = joinpath(REPO_ROOT, "svd_configs", "global_fit_pipeline.st.npoly3.configs.toml")
const SVD_WINTER_NC = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc"
const PSEUDO_RETRIEVAL_SZA = 0.0

# (B) Reflectance rebuild
const LEGENDRE_ORDER = 5
const RADIANCE_MEAN_RANGE = (15.0, 30.0)
# Spectral window for Legendre ρ fit / reconstruction (L1B red subset)
const REFLECTANCE_λ_MIN = 600.0
const REFLECTANCE_λ_MAX = 890.0
# Continuum windows for Legendre ρ fit only (self-defined, no absorption).
# LS is performed exclusively at these refs (all lie in [600, 890] nm).
const BASELINE_λ_REF = [
    602.918,
    612.732,
    617.605,
    620.06,
    622.53,
    639.79,
    662.068,
    665.795,
    668.265,
    669.518,
    670.755,
    673.245,
    675.73,
    678.205,
    679.445,
    680.68,
    681.918,
    745.535,
    748.04,
    754.295,
    773.075,
    774.337,
    776.832,
    779.335,
    781.842,
    859.595,
    862.112,
]

# (C) SIF
const SIF_STRENGTH_RANGE = (0.0, 0.5)

"""
    half_level_pressures_hpa(ps_hpa, ak, bk)

Reconstruct MERRA-2 half-level pressures [hPa] from surface pressure and `ak`/`bk`.
"""
function half_level_pressures_hpa(ps_hpa::Real, ak::AbstractVector, bk::AbstractVector)
    p_half_pa = ak .+ bk .* (Float64(ps_hpa) * 100.0)
    return p_half_pa ./ 100.0
end

"""
    full_level_pressures_hpa(p_half_hpa)

Mid-layer pressures [hPa] from half-level pressures.
"""
function full_level_pressures_hpa(p_half_hpa::AbstractVector)
    return (p_half_hpa[1:end-1] .+ p_half_hpa[2:end]) ./ 2
end

"""
    scale_height_km(p_half_hpa)

Approximate geometric height [km] from half-level pressure using a constant scale height.
"""
function scale_height_km(p_half_hpa::AbstractVector; H_m=8500.0)
    p_surf = p_half_hpa[end]
    z_m = -H_m .* log.(p_half_hpa ./ p_surf)
    return z_m ./ 1000.0
end

"""
    layer_reflectance_weights(n_layers; shape=:exponential, decay=4.0)

Return a normalized layer weight vector peaked near the surface (layer index = n_layers).
Uses exponential or Beta(2, 5) shapes — intentionally not χ².
"""
function layer_reflectance_weights(
    n_layers::Int;
    shape::Symbol=:exponential,
    decay::Float64=4.0,
)
    layers = collect(1:n_layers)
    w = if shape == :exponential
        # index 1 = top, n_layers = bottom; exponential favors the bottom
        exp.(decay .* (layers .- n_layers) ./ n_layers)
    elseif shape == :beta
        # Beta(2, 5) on [0, 1], evaluated from the bottom upward
        t = (layers .- 1) ./ max(n_layers - 1, 1)
        t .= 1.0 .- t   # flip so t=1 at the surface
        # unnormalized Beta pdf: t^(a-1) * (1-t)^(b-1)
        a, b = 2.0, 5.0
        (t .^ (a - 1)) .* ((1.0 .- t) .^ (b - 1))
    else
        error("Unknown reflectance shape: $shape")
    end
    w ./= sum(w)
    return w
end

"""
    layer_optical_depth(
        spectral_axis,
        p_full,
        temp,
        vcd_dry,
        vcd_h2o,
        o2_sitp,
        h2o_sitp;
        vmr_o2=0.21,
    )

Per-layer and cumulative optical depth on the high-resolution spectral grid.
Returns `(τ_layer, τ_cum)` with shapes `(n_layers, n_λ)`.
"""
function layer_optical_depth(
    spectral_axis::AbstractVector{<:Real},
    p_full::AbstractVector{<:Real},
    temp::AbstractVector{<:Real},
    vcd_dry::AbstractVector{<:Real},
    vcd_h2o::AbstractVector{<:Real},
    o2_sitp,
    h2o_sitp;
    vmr_o2::Float64=0.21,
)
    n_layers = length(p_full)
    n_λ = length(spectral_axis)
    _, p_lut, t_lut = o2_sitp.ranges
    p_min, p_max = extrema(p_lut)
    t_min, t_max = extrema(t_lut)

    τ_layer = zeros(Float64, n_layers, n_λ)
    for l in 1:n_layers
        p_l = clamp(p_full[l], p_min, p_max)
        t_l = clamp(temp[l], t_min, t_max)
        xsec_o2 = vec(o2_sitp(spectral_axis, p_l, t_l))
        xsec_h2o = vec(h2o_sitp(spectral_axis, p_l, t_l))
        τ_layer[l, :] .= xsec_o2 .* (vcd_dry[l] * vmr_o2) .+ xsec_h2o .* vcd_h2o[l]
    end
    τ_cum = cumsum(τ_layer, dims=1)
    return τ_layer, τ_cum
end

"""
    weighted_two_way_optical_depth(τ_cum, weights, amf_down, amf_up)

Effective two-way optical depth using a layer reflectance distribution:
down + up legs to/from each layer, weighted by where reflected light originates.
"""
function weighted_two_way_optical_depth(
    τ_cum::AbstractMatrix{<:Real},
    weights::AbstractVector{<:Real},
    amf_down::Real,
    amf_up::Real,
)
    n_layers, n_λ = size(τ_cum)
    length(weights) == n_layers || error("weights length must match n_layers")
    τ_2way = zeros(Float64, n_λ)
    for l in 1:n_layers
        wl = weights[l]
        τ_2way .+= wl .* (amf_down .+ amf_up) .* view(τ_cum, l, :)
    end
    return τ_2way
end

"""
    column_one_way_optical_depth(τ_layer, amf_up)

Unweighted column optical depth (sum over layers) times upward air-mass factor.
"""
function column_one_way_optical_depth(
    τ_layer::AbstractMatrix{<:Real},
    amf_up::Real,
)
    τ_col = vec(sum(τ_layer, dims=1))
    return Float64(amf_up) .* τ_col
end

"""
    load_profile(ds, profile_index)

Read one atmospheric profile from the MERRA2 transmittance NetCDF.
"""
function load_profile(ds, profile_index::Int)
    n_profiles = ds.dim["profile"]
    1 <= profile_index <= n_profiles || error("profile_index=$profile_index out of range 1:$n_profiles")

    temp = Float64.(ds["temperature"][profile_index, :])
    ps_hpa = Float64(ds["pressure"][profile_index])
    q = Float64.(ds["q"][profile_index, :])
    vcd_dry = Float64.(ds["vcd_dry"][profile_index, :])
    vcd_h2o = Float64.(ds["vcd_h2o"][profile_index, :])
    vmr_h2o = haskey(ds, "vmr_h2o_var") ?
        Float64.(ds["vmr_h2o_var"][profile_index, :]) :
        q
    amf = haskey(ds, "AMF") ? Float64(ds["AMF"][profile_index]) : 1.0
    ak = Float64.(ds.attrib["ak"])
    bk = Float64.(ds.attrib["bk"])

    p_half = half_level_pressures_hpa(ps_hpa, ak, bk)
    p_full = full_level_pressures_hpa(p_half)
    height_km = scale_height_km(p_half)

    return (
        profile_index=profile_index,
        temp=temp,
        ps_hpa=ps_hpa,
        q=q,
        vmr_h2o=vmr_h2o,
        vcd_dry=vcd_dry,
        vcd_h2o=vcd_h2o,
        amf=amf,
        p_half=p_half,
        p_full=p_full,
        height_km=height_km,
    )
end

function plot_profile(prof, weights; out_path)
    n_layers = length(prof.temp)
    layer_idx = collect(1:n_layers)
    height_mid = (prof.height_km[1:end-1] .+ prof.height_km[2:end]) ./ 2
    ykwargs = (yflip=true, ylabel="Layer index", yticks=layer_idx)

    p = plot(
        layout=(1, 3),
        size=(1100, 650),
        legend=false,
        plot_title="Profile $(prof.profile_index): ps=$(round(prof.ps_hpa, digits=1)) hPa, AMF=$(round(prof.amf, digits=2))",
    )

    # Pressure (bottom x) and height (top x) vs layer index
    plot!(p[1], prof.p_full, layer_idx;
        xlabel="Pressure [hPa]", title="Pressure / height", ykwargs...)
    pt = twinx(p[1])
    plot!(pt, height_mid, layer_idx;
        color=:dodgerblue, xlabel="Height [km]", legend=false, ykwargs...)

    plot!(p[2], prof.temp, layer_idx;
        xlabel="Temperature [K]", title="Temperature", ykwargs...)
    plot!(p[3], prof.vmr_h2o, layer_idx;
        xlabel="H2O VMR [mol/mol]", title="Water-vapor mixing ratio", ykwargs...)

    weights_path = replace(out_path, ".png" => "_reflectance_weights.png")
    pw = plot(weights, layer_idx;
        xlabel="Weight [-]",
        title="Layer reflectance weights ($(REFLECTANCE_SHAPE))",
        legend=false,
        size=(500, 650),
        ykwargs...)
    savefig(pw, weights_path)
    println("Saved reflectance-weight figure: ", weights_path)

    savefig(p, out_path)
    println("Saved profile figure: ", out_path)
    return p
end

function plot_transmittance(
    λ_hres,
    λ_band,
    T_solar,
    T2_hres,
    T2_band;
    out_path,
    T_solar_band=nothing,
    T2_atm_band=nothing,
)
    xlim = (minimum(λ_hres), maximum(λ_hres))
    p = plot(
        size=(1100, 450),
        legend=:outerright,
        xlabel="Wavelength [nm]",
        ylabel="Transmittance [-]",
        title="Two-way atmospheric × solar transmittance (hi-res + OCI)",
        xlims=xlim,
    )

    plot!(p, λ_hres, T_solar, label="Solar (hires)", color=:goldenrod, lw=1.0, alpha=0.25)
    plot!(p, λ_hres, T2_hres, label="2-way × solar (hires)", color=:crimson, lw=1.0, alpha=0.25)
    if T_solar_band !== nothing
        plot!(p, λ_band, T_solar_band, label="Solar (lres)", color=:darkorange,
            lw=2.0, marker=:diamond, ms=3, markerstrokewidth=0.5)
    end
    if T2_atm_band !== nothing
        plot!(p, λ_band, T2_atm_band, label="2-way atm (lres)", color=:dodgerblue,
            lw=2.0, marker=:utriangle, ms=3, markerstrokewidth=0.5)
    end
    if T_solar_band !== nothing && T2_atm_band !== nothing
        plot!(p, λ_band, T2_atm_band .* T_solar_band,
            label="2-way atm (lres) × solar (lres)", color=:mediumpurple,
            lw=2.0, marker=:square, ms=3, markerstrokewidth=0.5)
    end
    plot!(p, λ_band, T2_band, label="2-way × solar (lres)", color=:darkred,
        lw=2.0, marker=:circle, ms=4, markerstrokewidth=0.5)

    savefig(p, out_path)
    sz = filesize(out_path)
    sz > 1000 || error("savefig wrote empty/invalid PNG ($sz bytes): $out_path")
    println("Saved transmittance figure: ", out_path, " (", sz, " bytes)")
    return p
end

# ── (B) Reflectance rebuild ───────────────────────────────────────────────────

# Spectral helpers `_normalized_grid` / `_legendre_design_matrix` come from svd_helpers.jl.
const center_wavelength_vec = _normalized_grid
legendre_design(z, order) = _legendre_design_matrix(z, order)

"""Nearest-band indices for every `λ_ref` that falls inside `λ` (and optional [lo, hi])."""
function baseline_band_indices(
    λ::AbstractVector{<:Real},
    λ_ref::AbstractVector{<:Real}=BASELINE_λ_REF;
    λ_lo::Union{Nothing, Float64}=nothing,
    λ_hi::Union{Nothing, Float64}=nothing,
)
    isempty(λ_ref) && error("Empty baseline wavelength list")
    data_lo, data_hi = extrema(λ)
    refs = collect(Float64.(λ_ref))
    if λ_lo !== nothing && λ_hi !== nothing
        refs = filter(r -> λ_lo <= r <= λ_hi, refs)
    end
    # keep refs that can be mapped onto available L1B bands
    refs = filter(r -> data_lo - 2 <= r <= data_hi + 2, refs)
    isempty(refs) && error("No BASELINE_λ_REF fall inside λ span [$data_lo, $data_hi]")
    inds = Int[]
    for r in refs
        i = argmin(abs.(λ .- r))
        if abs(λ[i] - r) > 2.0
            @warn "Baseline λ=$r nm maps to nearest band $(λ[i]) nm (Δ=$(round(abs(λ[i]-r), digits=2)) nm)"
        end
        push!(inds, i)
    end
    return unique(inds), refs
end

"""Subset an L1B observation to [λ_lo, λ_hi]."""
function subset_obs_wavelength(obs, λ_lo::Float64, λ_hi::Float64)
    ind = findall(λ_lo .<= obs.λ .<= λ_hi)
    isempty(ind) && error("No L1B bands in [$λ_lo, $λ_hi] nm")
    return (
        λ=obs.λ[ind],
        E=obs.E[ind],
        R=obs.R[ind],
        sza=obs.sza,
        vza=obs.vza,
        pixel=obs.pixel,
        scan=obs.scan,
        pace_path=obs.pace_path,
    )
end

"""
    pick_random_l1b_pixel(pace_path; max_tries=500)

Load one random finite L1B **full red** radiance spectrum (all red bands).
Prefers ocean (`watermask==1`) with moderate mean radiance for stable ρ fits.
Optional `mean_λ_min`/`mean_λ_max` only restrict the mean-radiance QC window.
"""
function pick_random_l1b_pixel(
    pace_path::AbstractString;
    max_tries::Int=500,
    mean_R_lo::Float64=5.0,
    mean_R_hi::Float64=80.0,
    mean_λ_min::Union{Nothing, Float64}=nothing,
    mean_λ_max::Union{Nothing, Float64}=nothing,
)
    ds = Dataset(pace_path)
    λ_all = Float64.(ds["red_wavelength"][:])
    E_all = Float64.(ds["red_solar_irradiance"][:])
    n_pix, n_scan, n_band = size(ds["radiance_red"])
    length(λ_all) == n_band || error("red_wavelength length mismatch with radiance_red")
    has_water = haskey(ds, "watermask")

    mean_ind = if mean_λ_min !== nothing && mean_λ_max !== nothing
        findall(mean_λ_min .< λ_all .< mean_λ_max)
    else
        collect(1:n_band)
    end
    isempty(mean_ind) && error("No bands in mean-radiance QC window")

    pixel = scan = 0
    R = Float64[]
    sza = vza = NaN
    for _ in 1:max_tries
        pixel = rand(1:n_pix)
        scan = rand(1:n_scan)
        if has_water
            w = ds["watermask"][pixel, scan]
            (ismissing(w) || Int(w) != 1) && continue
        end
        R_try = Float64.(coalesce.(ds["radiance_red"][pixel, scan, :], NaN))
        sza_try = Float64(coalesce(ds["solar_zenith"][pixel, scan], NaN))
        vza_try = Float64(coalesce(ds["sensor_zenith"][pixel, scan], NaN))
        mean_R = mean(R_try[mean_ind])
        if all(isfinite, R_try) && isfinite(sza_try) && sza_try < 85.0 &&
           all(R_try .> 0) && mean_R_lo <= mean_R <= mean_R_hi
            R = R_try
            sza = sza_try
            vza = vza_try
            break
        end
    end
    close(ds)
    isempty(R) && error("Failed to find a valid L1B pixel in $max_tries tries")

    return (
        λ=λ_all,
        E=E_all,
        R=R,
        sza=sza,
        vza=vza,
        pixel=pixel,
        scan=scan,
        pace_path=pace_path,
    )
end

"""
    fit_legendre_reflectance(λ, E, R, sza; order=7)

Least-squares fit of radiance continuum
    R(λ) = Σ aⱼ Pⱼ(z),  z = _normalized_grid(λ) ∈ [-1, 1]
directly to **R_L1B** at all `BASELINE_λ_REF` bands in
[REFLECTANCE_λ_MIN, REFLECTANCE_λ_MAX]. Coefficients are then evaluated on the
full `λ` grid. Apparent ρ = R·π/(E·cos(SZA)) is derived afterwards for plotting.
"""
function fit_legendre_reflectance(
    λ::AbstractVector{<:Real},
    E::AbstractVector{<:Real},
    R::AbstractVector{<:Real},
    sza::Real;
    order::Int=LEGENDRE_ORDER,
    λ_ref::AbstractVector{<:Real}=BASELINE_λ_REF,
    λ_lo::Float64=REFLECTANCE_λ_MIN,
    λ_hi::Float64=REFLECTANCE_λ_MAX,
)
    λc = _normalized_grid(λ)   # z ∈ [-1, 1], same as SVD pipeline
    bl_ind, refs_used = baseline_band_indices(λ, λ_ref; λ_lo=λ_lo, λ_hi=λ_hi)
    n_bl = length(bl_ind)
    n_bl >= order + 1 ||
        error("Need ≥ $(order+1) baseline bands for order=$order; got $n_bl from λ_ref in [$λ_lo, $λ_hi]")
    μ = cosd(Float64(sza))
    μ > 0 || error("cos(SZA) must be positive; got SZA=$sza")

    # Fit R_L1B ONLY at baseline bands
    K_bl = _legendre_design_matrix(λc[bl_ind], order)
    coeffs = K_bl \ R[bl_ind]

    # Reconstruct continuum radiance on full spectrum
    K_full = _legendre_design_matrix(λc, order)
    R_fit = K_full * coeffs
    # Derived apparent reflectance (for plotting / downstream use)
    ρ = R_fit .* π ./ (E .* μ)
    return (
        ρ=ρ, R_fit=R_fit, coeffs=coeffs, bl_ind=bl_ind, λc=λc, μ=μ,
        n_baseline=n_bl, refs_used=refs_used,
    )
end

"""
Rescale so mean reconstructed radiance lies in [lo, hi] (uniform random target).
If `λ` and a window are given, the mean is taken only inside that window.
"""
function rescale_mean_radiance(
    ρ::AbstractVector{<:Real},
    R_fit::AbstractVector{<:Real};
    lo::Float64=RADIANCE_MEAN_RANGE[1],
    hi::Float64=RADIANCE_MEAN_RANGE[2],
    λ::Union{Nothing, AbstractVector{<:Real}}=nothing,
    mean_λ_min::Union{Nothing, Float64}=nothing,
    mean_λ_max::Union{Nothing, Float64}=nothing,
)
    if λ !== nothing && mean_λ_min !== nothing && mean_λ_max !== nothing
        ind = findall(mean_λ_min .< λ .< mean_λ_max)
        isempty(ind) && error("No bands in rescale window [$mean_λ_min, $mean_λ_max]")
        mean_fit = mean(R_fit[ind])
    else
        mean_fit = mean(R_fit)
    end
    mean_fit > 0 || error("Fitted radiance mean must be positive (got $mean_fit)")
    target = lo + (hi - lo) * rand()
    scale = target / mean_fit
    return (ρ=ρ .* scale, R=R_fit .* scale, scale=scale, target_mean=target, mean_fit=mean_fit)
end

function plot_reflectance_rebuild(obs, fit, resc; out_path, highlight_λ_min=nothing, highlight_λ_max=nothing)
    λ_lo = minimum(obs.λ)
    λ_hi = maximum(obs.λ)
    xlim = (λ_lo, λ_hi)

    p = plot(
        size=(1100, 450),
        legend=:outerright,
        xlabel="Wavelength [nm]",
        ylabel="Radiance",
        xlims=xlim,
        title="L1B px=$(obs.pixel), scan=$(obs.scan), SZA=$(round(obs.sza, digits=1))°, " *
              "n_baseline=$(length(fit.bl_ind)), " *
              "mean R: L1B=$(round(mean(obs.R), digits=2)), fit=$(round(resc.mean_fit, digits=2)), " *
              "rescaled→$(round(resc.target_mean, digits=2))",
    )

    plot!(p, obs.λ, obs.R, label="L1B radiance", color=:black, lw=2.0)
    plot!(p, obs.λ, fit.R_fit, label="Fitted radiance (Legendre $LEGENDRE_ORDER on R)",
        color=:dodgerblue, lw=1.5, ls=:dash)
    plot!(p, obs.λ, resc.R, label="Rescaled radiance", color=:crimson, lw=1.5)
    scatter!(p, obs.λ[fit.bl_ind], obs.R[fit.bl_ind],
        label="Baseline ($(length(fit.bl_ind)) refs)", color=:orange, ms=5, markerstrokewidth=0)
    if highlight_λ_min !== nothing && highlight_λ_max !== nothing
        vspan!(p, [highlight_λ_min, highlight_λ_max]; color=:gray, alpha=0.12, label="Working window")
    end

    savefig(p, out_path)
    println("Saved reflectance figure: ", out_path)
    return p
end

"""
    rebuild_reflectance_from_l1b(pace_path; order, out_path)

Full (B) workflow: random L1B pixel → subset to [600, 890] nm →
Legendre fit at **all** `BASELINE_λ_REF` → rescale mean radiance.
"""
function rebuild_reflectance_from_l1b(
    pace_path::AbstractString;
    order::Int=LEGENDRE_ORDER,
    out_path::AbstractString,
    max_fit_tries::Int=40,
    mean_λ_min::Union{Nothing, Float64}=nothing,
    mean_λ_max::Union{Nothing, Float64}=nothing,
    fit_λ_min::Float64=REFLECTANCE_λ_MIN,
    fit_λ_max::Float64=REFLECTANCE_λ_MAX,
)
    local obs, fit, resc
    ok = false
    for attempt in 1:max_fit_tries
        obs_full = pick_random_l1b_pixel(
            pace_path;
            mean_λ_min=mean_λ_min,
            mean_λ_max=mean_λ_max,
        )
        obs = subset_obs_wavelength(obs_full, fit_λ_min, fit_λ_max)
        fit = fit_legendre_reflectance(
            obs.λ, obs.E, obs.R, obs.sza;
            order=order, λ_lo=fit_λ_min, λ_hi=fit_λ_max,
        )
        frac_neg = count(<(0), fit.R_fit) / length(fit.R_fit)
        if mean(fit.R_fit) > 0 && frac_neg < 0.05
            ok = true
            break
        end
        println("  fit attempt $attempt rejected (mean R_fit=$(round(mean(fit.R_fit), digits=2)), " *
                "frac_neg_R=$(round(frac_neg, digits=3)), n_bl=$(length(fit.bl_ind)))")
    end
    ok || error("Could not obtain a stable Legendre reflectance fit in $max_fit_tries tries")

    println("Selected L1B pixel=$(obs.pixel), scan=$(obs.scan), SZA=$(round(obs.sza, digits=2))°, VZA=$(round(obs.vza, digits=2))°")
    println("  fit λ span: [$(round(minimum(obs.λ), digits=2)), $(round(maximum(obs.λ), digits=2))] nm ($(length(obs.λ)) bands)")
    println("  mean L1B radiance: $(round(mean(obs.R), digits=3))")
    println("  baseline windows used: $(fit.n_baseline) / $(length(BASELINE_λ_REF)) refs  →  $(round.(fit.refs_used; digits=2))")
    println("  mean fitted radiance: $(round(mean(fit.R_fit), digits=3))")

    resc = rescale_mean_radiance(
        fit.ρ, fit.R_fit;
        λ=obs.λ, mean_λ_min=mean_λ_min, mean_λ_max=mean_λ_max,
    )
    println("  rescale factor: $(round(resc.scale, digits=4)) → target mean=$(round(resc.target_mean, digits=3))")

    plot_reflectance_rebuild(
        obs, fit, resc;
        out_path=out_path,
        highlight_λ_min=mean_λ_min,
        highlight_λ_max=mean_λ_max,
    )
    return (obs=obs, fit=fit, resc=resc)
end

# ── (C) SIF + one-way T₁ ──────────────────────────────────────────────────────

"""
Load one random SIF shape from `SIF_shapes` in the library JLD2 and map to `λ_dst`
(cubic spline onto OCI bands when the library grid is evenly spaced).
"""
function load_random_sif_shape(
    sif_path::AbstractString,
    λ_dst::AbstractVector{<:Real},
)
    sif = JLD2.load(MWEF.must_exist(sif_path))
    haskey(sif, "SIF_shapes") && haskey(sif, "SIF_wavelen") ||
        error("SIF file must contain SIF_shapes and SIF_wavelen")
    shapes = Matrix{Float64}(sif["SIF_shapes"])
    λ_ref = Float64.(sif["SIF_wavelen"])
    size(shapes, 1) == length(λ_ref) ||
        error("SIF_shapes first dimension must match SIF_wavelen length")
    idx = rand(1:size(shapes, 2))
    shape_band = map_sif_shape_to_bands(λ_ref, shapes[:, idx], λ_dst)
    return (
        λ=collect(Float64.(λ_dst)),
        shape=shape_band,
        library_index=idx,
    )
end

"""Normalize shape to unit peak, then scale peak to `strength` ∈ [lo, hi]."""
function rescale_sif_strength(
    shape::AbstractVector{<:Real};
    lo::Float64=SIF_STRENGTH_RANGE[1],
    hi::Float64=SIF_STRENGTH_RANGE[2],
)
    s = Float64.(shape)
    peak = maximum(abs.(s))
    peak > 0 || error("SIF shape is all zero")
    strength = lo + (hi - lo) * rand()
    return (SIF=(s ./ peak) .* strength, strength=strength, peak_ref=peak)
end

"""
Build unweighted one-way atmospheric transmittance on hi-res and OCI bands (no solar).
"""
function build_one_way_transmittance(
    τ_layer::AbstractMatrix{<:Real},
    amf_up::Real,
    K::AbstractMatrix{<:Real},
)
    τ_1way = column_one_way_optical_depth(τ_layer, amf_up)
    T1_atm = exp.(-τ_1way)
    T1_band = vec(K * T1_atm)
    return (τ_1way=τ_1way, T1_atm=T1_atm, T1_band=T1_band)
end

"""Shared wavelength x-limits for stacked / side-by-side spectral panels."""
function _spectral_xlim(λ::AbstractVector{<:Real})
    λf = Float64.(λ)
    return (minimum(λf), maximum(λf))
end

const _SUBPLOT_X_MARGINS = (
    left_margin=18Plots.mm,
    right_margin=12Plots.mm,
    top_margin=2Plots.mm,
)

function _apply_shared_x!(p, n::Int, xlim::Tuple{<:Real, <:Real})
    m = _SUBPLOT_X_MARGINS
    for i in 1:n
        plot!(p[i];
            xlims=xlim,
            xticks=640:20:760,
            left_margin=m.left_margin,
            right_margin=m.right_margin,
            top_margin=m.top_margin,
        )
    end
    return p
end

"""Two-panel figure: (a) SIF and SIF×T₁, (b) T₁."""
function plot_sif_and_t1(sif, T1; out_path::AbstractString)
    R_sif_toa = sif.SIF .* T1.T1_band
    xlim = _spectral_xlim(sif.λ)
    p = plot(
        layout=(1, 2),
        size=(1100, 420),
        legend=:outerright,
        link=:x,
        plot_title="SIF shape #$(sif.library_index), strength=$(round(sif.strength, digits=3))",
    )
    plot!(p[1], sif.λ, sif.SIF, label="SIF (scaled)", color=:forestgreen, lw=2.0)
    plot!(p[1], sif.λ, R_sif_toa, label="SIF × T₁", color=:darkgreen, lw=1.5, ls=:dash)
    plot!(p[1], xlabel="Wavelength [nm]", ylabel="Radiance", title="SIF spectrum")
    plot!(p[2], sif.λ, T1.T1_band, label="T₁ (1-way atm, lres)", color=:dodgerblue, lw=2.0)
    plot!(p[2], xlabel="Wavelength [nm]", ylabel="Transmittance [-]", title="One-way transmittance")
    _apply_shared_x!(p, 2, xlim)

    savefig(p, out_path)
    println("Saved SIF figure: ", out_path)
    return p
end

"""
    add_sif_component(τ_layer, amf_up, K, sif_path, λ_band; out_path)

(C) workflow: unweighted T₁, random SIF shape rescaled to [0, 0.5], plot.
"""
function add_sif_component(
    τ_layer::AbstractMatrix{<:Real},
    amf_up::Real,
    K::AbstractMatrix{<:Real},
    sif_path::AbstractString,
    λ_band::AbstractVector{<:Real};
    out_path::AbstractString,
)
    T1 = build_one_way_transmittance(τ_layer, amf_up, K)
    picked = load_random_sif_shape(sif_path, λ_band)
    scaled = rescale_sif_strength(picked.shape)
    R_sif_toa = scaled.SIF .* T1.T1_band
    sif = (
        λ=picked.λ,
        SIF=scaled.SIF,
        R_sif_toa=R_sif_toa,
        strength=scaled.strength,
        library_index=picked.library_index,
        peak_ref=scaled.peak_ref,
    )

    println("  SIF library index: $(sif.library_index)")
    println("  SIF strength (peak): $(round(sif.strength, digits=4))")
    println("  mean SIF: $(round(mean(sif.SIF), digits=4))")
    println("  mean SIF×T₁: $(round(mean(sif.R_sif_toa), digits=4))")
    println("  mean T₁: $(round(mean(T1.T1_band), digits=4))")

    plot_sif_and_t1(sif, T1; out_path=out_path)
    return (sif=sif, T1=T1)
end

# ── (D) TOA radiance + SNR white noise ────────────────────────────────────────

"""
Map SIF library shape `y_src(λ_src)` onto OCI band centers `λ_dst`.

Uses cubic spline interpolation when `λ_src` is evenly spaced (as for
`SIF_wavelen` in `SIF_singular_vector.jld2`), otherwise linear interpolation.
Values outside the source wavelength range are set to zero.
"""
function map_sif_shape_to_bands(
    λ_src::AbstractVector{<:Real},
    y_src::AbstractVector{<:Real},
    λ_dst::AbstractVector{<:Real},
)
    length(λ_src) == length(y_src) || error("λ_src / y_src length mismatch")
    λs = collect(Float64.(λ_src))
    ys = Float64.(y_src)
    λt = collect(Float64.(λ_dst))
    dλ = diff(λs)
    step = dλ[1]
    tol = max(1e-10, abs(step) * 1e-8)
    use_cubic = all(abs.(dλ .- step) .<= tol)
    if use_cubic
        λ_knots = range(λs[1], step=step, length=length(λs))
        itp = CubicSplineInterpolation(λ_knots, ys; extrapolation_bc=Line())
        out = itp.(λt)
        out[(λt .< λs[1]) .| (λt .> λs[end])] .= 0.0
        return out
    else
        itp = LinearInterpolation(λs, ys; extrapolation_bc=0.0)
        return itp.(λt)
    end
end

"""
Map `y_src(λ_src)` onto `λ_dst` by nearest-band sampling (OCI band grids align closely).
"""
function map_spectrum_to_bands(
    λ_src::AbstractVector{<:Real},
    y_src::AbstractVector{<:Real},
    λ_dst::AbstractVector{<:Real},
)
    length(λ_src) == length(y_src) || error("λ_src / y_src length mismatch")
    λs = Float64.(λ_src)
    ys = Float64.(y_src)
    out = similar(Float64.(λ_dst))
    for (k, λk) in enumerate(λ_dst)
        i = argmin(abs.(λs .- Float64(λk)))
        out[k] = ys[i]
    end
    return out
end

"""
    noise_std_from_snr(R, c1, c2) -> σ

PACE OCI SNR model:  σ² = c1 + c2 · R  (R clipped at 0).
"""
function noise_std_from_snr(
    R::AbstractVector{<:Real},
    c1::AbstractVector{<:Real},
    c2::AbstractVector{<:Real},
)
    n = length(R)
    length(c1) == n || error("c1 length $(length(c1)) ≠ spectrum length $n")
    length(c2) == n || error("c2 length $(length(c2)) ≠ spectrum length $n")
    return sqrt.(Float64.(c1) .+ Float64.(c2) .* max.(Float64.(R), 0.0))
end

"""
    build_toa_with_noise(λ, R_cont, T2; c1, c2, R_sif_toa=nothing)

Background  R_bg = R_cont · T₂;  clean TOA  R_toa = R_bg + R_sif_toa.
Gaussian noise uses σ = √(c1 + c2 · R) evaluated on the spectrum being noised.
"""
function build_toa_with_noise(
    λ::AbstractVector{<:Real},
    R_cont::AbstractVector{<:Real},
    T2::AbstractVector{<:Real};
    c1::AbstractVector{<:Real},
    c2::AbstractVector{<:Real},
    R_sif_toa::Union{Nothing, AbstractVector{<:Real}}=nothing,
)
    length(λ) == length(R_cont) == length(T2) ||
        error("λ / R_cont / T2 length mismatch")
    R_cont = collect(Float64.(R_cont))
    T2 = collect(Float64.(T2))
    R_bg = R_cont .* T2
    R_sif = R_sif_toa === nothing ? zeros(length(R_bg)) : collect(Float64.(R_sif_toa))
    length(R_sif) == length(R_bg) || error("R_sif_toa length mismatch")
    R_toa = R_bg .+ R_sif

    σ_sif = noise_std_from_snr(R_toa, c1, c2)
    noise_sif = randn(length(R_toa)) .* σ_sif
    R_noisy_sif = R_toa .+ noise_sif

    σ_bg = noise_std_from_snr(R_bg, c1, c2)
    noise_bg = randn(length(R_bg)) .* σ_bg
    R_noisy_bg = R_bg .+ noise_bg

    snr_sif = @. ifelse(σ_sif > 0, R_toa / σ_sif, NaN)
    return (
        λ=collect(Float64.(λ)),
        R_cont=R_cont,
        T2=T2,
        R_bg=R_bg,
        R_sif_toa=R_sif,
        R_toa=R_toa,
        R_noisy=R_noisy_sif,
        R_noisy_bg=R_noisy_bg,
        σ=σ_sif,
        σ_bg=σ_bg,
        noise=noise_sif,
        noise_bg=noise_bg,
        snr=snr_sif,
        c1=collect(Float64.(c1)),
        c2=collect(Float64.(c2)),
    )
end

"""Plot background and TOA spectra, clean and noisy, with and without SIF."""
function plot_toa_surrogate(toa; out_path::AbstractString)
    p = plot(
        size=(1100, 450),
        legend=:outerright,
        xlabel="Wavelength [nm]",
        ylabel="Radiance",
        title="Surrogate TOA  |  mean bg=$(round(mean(toa.R_bg), digits=2)), " *
              "mean TOA=$(round(mean(toa.R_toa), digits=2)), " *
              "median SNR=$(round(median(filter(isfinite, toa.snr)), digits=1))",
    )
    plot!(p, toa.λ, toa.R_bg, label="Background (cont × T₂, no SIF)", color=:gray, lw=1.5, ls=:dot)
    plot!(p, toa.λ, toa.R_toa, label="Clean TOA (+ SIF×T₁)", color=:dodgerblue, lw=2.0)
    plot!(p, toa.λ, toa.R_noisy, label="Noisy TOA (+ SIF)", color=:crimson, lw=1.2, alpha=0.85)
    plot!(p, toa.λ, toa.R_noisy_bg, label="Noisy background (no SIF)", color=:orange, lw=1.0, ls=:dash, alpha=0.8)
    if any(!iszero, toa.R_sif_toa)
        plot!(p, toa.λ, toa.R_sif_toa, label="SIF×T₁ component", color=:forestgreen, lw=1.2, ls=:dashdot)
    end

    savefig(p, out_path)
    println("Saved TOA figure: ", out_path)
    return p
end

"""
    rebuild_toa_with_noise(λ_band, T2_band, refl, sif; band_snr_coeffs, out_path)

(D) workflow: background = continuum × T₂, add SIF×T₁, SNR noise, plot.
"""
function rebuild_toa_with_noise(
    λ_band::AbstractVector{<:Real},
    T2_band::AbstractVector{<:Real},
    refl,
    sif;
    band_snr_coeffs,
    out_path::AbstractString,
)
    band_snr_coeffs === nothing &&
        error("band_snr_coeffs is nothing; enable [kernel] use_band_snr and a valid pace_snr_file")
    c1 = Float64.(band_snr_coeffs["c1"])
    c2 = Float64.(band_snr_coeffs["c2"])
    length(c1) == length(λ_band) ||
        error("SNR c1 length $(length(c1)) ≠ λ_band length $(length(λ_band))")

    R_cont_band = map_spectrum_to_bands(refl.obs.λ, refl.resc.R, λ_band)
    toa = build_toa_with_noise(
        λ_band, R_cont_band, T2_band;
        c1=c1, c2=c2, R_sif_toa=sif.sif.R_sif_toa,
    )

    println("  bands: $(length(toa.λ))  λ ∈ [$(round(minimum(toa.λ), digits=2)), $(round(maximum(toa.λ), digits=2))] nm")
    println("  mean continuum (mapped): $(round(mean(toa.R_cont), digits=3))")
    println("  mean T₂: $(round(mean(toa.T2), digits=4))")
    println("  mean background: $(round(mean(toa.R_bg), digits=3))")
    println("  mean SIF×T₁: $(round(mean(toa.R_sif_toa), digits=4))")
    println("  mean clean TOA: $(round(mean(toa.R_toa), digits=3))")
    println("  mean noisy TOA: $(round(mean(toa.R_noisy), digits=3))")
    println("  σ range: [$(round(minimum(toa.σ), digits=4)), $(round(maximum(toa.σ), digits=4))]")
    println("  median SNR: $(round(median(filter(isfinite, toa.snr)), digits=2))")

    plot_toa_surrogate(toa; out_path=out_path)
    return toa
end

# ── (E) Pseudo SVD retrieval ──────────────────────────────────────────────────

"""Pick a random L1B NetCDF from `base_dir` (prefers filenames containing 'granule')."""
function pick_random_l1b_file(base_dir::AbstractString)
    all_nc = filter(f -> endswith(lowercase(f), ".nc"), readdir(base_dir, join=true))
    granules = filter(f -> occursin("granule", lowercase(f)), all_nc)
    candidates = isempty(granules) ? all_nc : granules
    isempty(candidates) && error("No L1B NetCDF files in $base_dir")
    return rand(candidates)
end

"""Read `red_solar_irradiance` and earth–sun correction; map E onto retrieval bands."""
function load_l1b_solar_on_bands(pace_path::AbstractString, λ_target::AbstractVector{<:Real})
    isfile(pace_path) || error("L1B not found: $pace_path")
    ds = Dataset(pace_path)
    λ_all = Float64.(ds["red_wavelength"][:])
    E_all = Float64.(ds["red_solar_irradiance"][:])
    esd = Float64(get(ds.attrib, "earth_sun_distance_correction", 1.0))
    close(ds)
    E_band = map_spectrum_to_bands(λ_all, E_all, λ_target)
    return (E=E_band, esd=esd, pace_path=pace_path)
end

"""
Build SVD retrieval priors / SNR / basis on the observation band grid `λ_ctx`.
Uses [data], [spectral], [fit], [fit.svd], [batch_fit] from the pipeline TOML.
"""
function prepare_svd_retrieval_setup(retrieval_cfg::Dict{String, Any}, λ_ctx::AbstractVector{Float64})
    spectral_cfg = get(retrieval_cfg, "spectral", Dict{String, Any}())
    fit_cfg = get(retrieval_cfg, "fit", Dict{String, Any}())
    svd_cfg = get(fit_cfg, "svd", Dict{String, Any}())
    data_cfg = get(retrieval_cfg, "data", Dict{String, Any}())
    batch_cfg = get(retrieval_cfg, "batch_fit", Dict{String, Any}())
    isempty(svd_cfg) && error("Retrieval config missing [fit.svd]")

    λ_min = Float64(get(spectral_cfg, "lambda_min_nm", 640.0))
    λ_max = Float64(get(spectral_cfg, "lambda_max_nm", 756.0))
    base_dir = String(get(data_cfg, "base_dir", joinpath(REPO_ROOT, "demo_example")))

    summer_rel = String(get(data_cfg, "summer_nc", TRANS_NC))
    winter_rel = String(get(data_cfg, "winter_nc", SVD_WINTER_NC))
    summer_nc = isabspath(summer_rel) ? summer_rel : joinpath(base_dir, summer_rel)
    winter_nc = isabspath(winter_rel) ? winter_rel : joinpath(base_dir, winter_rel)

    sif_rel = get(data_cfg, "sif_file", "SIF_singular_vector.jld2")
    sif_path = isabspath(String(sif_rel)) ? String(sif_rel) : joinpath(base_dir, String(sif_rel))
    isfile(sif_path) || error("SIF file not found: $sif_path")

    n_pc = Int(get(svd_cfg, "n_pc", 5))
    n_leg = Int(get(svd_cfg, "n_legendre", 3))
    log_trans = Bool(get(svd_cfg, "svd_log_transform", false))
    sif_nev = Int(get(spectral_cfg, "sif_nev", 1))
    normalize_sif = Bool(get(spectral_cfg, "normalize_sif_first_ev", true))

    svd_basis = load_svd_basis(summer_nc, winter_nc, λ_ctx; λ_min=λ_min, λ_max=λ_max, n_pc=n_pc, log_transform=log_trans)
    PCs = Float64.(svd_basis.PCs[:, 1:n_pc])
    sif_basis = MWEF.load_sif_basis(sif_path, λ_ctx; nEV=sif_nev, normalize=normalize_sif)
    layout = svd_state_layout(; n_pc=n_pc, n_legendre=n_leg, n_ev=size(sif_basis, 2))

    prior_sigma_default = Float64(get(fit_cfg, "prior_sigma_default", 1e12))
    prior_min_sigma = Float64(get(fit_cfg, "prior_min_sigma", 1e-3))
    sif_sigma = Float64(get(fit_cfg, "sif_sigma", 1e12))
    alpha_mean = Float64(get(svd_cfg, "alpha_prior_mean", 1.0))
    alpha_sigma = Float64(get(svd_cfg, "alpha_prior_sigma", 0.3))
    use_leg01 = Bool(get(svd_cfg, "use_legendre01_prior", true))
    leg01_frac = Float64(get(svd_cfg, "legendre01_prior_sigma_fraction", 1.0))
    use_leghig = Bool(get(svd_cfg, "use_legendre_higher_prior", true))
    leg_higher_sigma = Float64(get(svd_cfg, "legendre_higher_sigma", 1.0))
    pc_prior_mode = String(get(svd_cfg, "pc_prior_mode", "loading_variance"))
    pc_sigma_scale = Float64(get(svd_cfg, "pc_prior_sigma_scale", 1.0))

    x0 = zeros(Float64, layout.n_state)
    x0[first(layout.idx_alpha)] = alpha_mean
    x0[first(layout.idx_legendre)] = 1.0
    prior_sigma = fill(prior_sigma_default, layout.n_state)
    prior_sigma[first(layout.idx_alpha)] = max(alpha_sigma, prior_min_sigma)
    if pc_prior_mode == "loading_variance"
        n_prof = svd_basis.n_profiles
        for k in 1:n_pc
            prior_sigma[k] = max(svd_basis.S[k] / sqrt(Float64(n_prof)) * pc_sigma_scale, prior_min_sigma)
        end
    end
    if use_leghig && length(layout.idx_legendre) >= 3
        for j in 3:length(layout.idx_legendre)
            prior_sigma[layout.idx_legendre[j]] = max(leg_higher_sigma, prior_min_sigma)
        end
    end
    prior_sigma[layout.idx_sif] .= max(sif_sigma, prior_min_sigma)

    z = _normalized_grid(λ_ctx)
    A01 = hcat(ones(length(z)), z)
    lower = fill(-Inf, layout.n_state)
    upper = fill(Inf, layout.n_state)
    x_scale = ones(Float64, layout.n_state)
    for k in 1:n_pc
        x_scale[k] = prior_sigma[k]
    end
    x_scale[first(layout.idx_alpha)] = prior_sigma[first(layout.idx_alpha)]

    kernel_cfg = get(retrieval_cfg, "kernel", Dict{String, Any}())
    use_band_snr = haskey(fit_cfg, "use_band_snr") ? Bool(fit_cfg["use_band_snr"]) :
                   Bool(get(kernel_cfg, "use_band_snr", true))
    meas_sigma = Float64(get(fit_cfg, "meas_sigma", 0.01))
    pace_snr_rel = get(data_cfg, "pace_snr_file", "PACE_OCI_L1BLUT_baseline_SNR_1.1.txt")
    pace_snr_path = isabspath(String(pace_snr_rel)) ? String(pace_snr_rel) : joinpath(base_dir, String(pace_snr_rel))
    band_snr_coeffs = use_band_snr ?
        load_pace_band_snr_coeffs(pace_snr_path, λ_ctx; λ_min=λ_min, λ_max=λ_max) :
        nothing
    noise_degrade = parse_noise_degrade(retrieval_cfg)
    if noise_degrade !== nothing && band_snr_coeffs !== nothing
        band_snr_coeffs = copy_snr_coeffs_with_degrade(band_snr_coeffs, λ_ctx, noise_degrade)
        n_deg = length(noise_degrade_band_mask(λ_ctx, noise_degrade))
        println("  noise_degrade retrieval: λ∈[$(noise_degrade.lambda_min_nm), $(noise_degrade.lambda_max_nm)] nm, σ×$(noise_degrade.sigma_factor) ($n_deg bands)")
    end

    lm = (
        lambda0 = Float64(get(fit_cfg, "lm_lambda0", 1.0)),
        lambda_up = Float64(get(fit_cfg, "lm_lambda_up", 5.0)),
        lambda_down = Float64(get(fit_cfg, "lm_lambda_down", 0.7)),
        lambda_min = Float64(get(fit_cfg, "lm_lambda_min", 1e-8)),
        lambda_max = Float64(get(fit_cfg, "lm_lambda_max", 1e8)),
        max_inner = Int(get(fit_cfg, "lm_max_inner", 24)),
    )
    conv = (
        dx_rel_tol = Float64(get(fit_cfg, "conv_dx_rel_tol", 1e-6)),
        rmse_rel_tol = Float64(get(fit_cfg, "conv_rmse_rel_tol", 1e-6)),
        rmse_abs_tol = Float64(get(fit_cfg, "conv_rmse_abs_tol", 1e-6)),
        enabled = Bool(get(fit_cfg, "conv_stall_enable", true)),
        window = Int(get(fit_cfg, "conv_stall_window", 3)),
        redchi2_target = Float64(get(fit_cfg, "conv_stall_redchi2_target", 5.0)),
        redchi2_abs_tol = Float64(get(fit_cfg, "conv_stall_redchi2_abs_tol", 0.1)),
        redchi2_rel_tol = Float64(get(fit_cfg, "conv_stall_redchi2_rel_tol", 0.03)),
        dx_rel_tol_stall = Float64(get(fit_cfg, "conv_stall_dx_rel_tol", 5e-3)),
    )
    max_outer_steps = Int(get(batch_cfg, "max_outer_steps", 15))

    return (
        λ_ctx=λ_ctx,
        PCs=PCs,
        sif_basis=sif_basis,
        layout=layout,
        x0=x0,
        prior_sigma=prior_sigma,
        A01=A01,
        lower=lower,
        upper=upper,
        x_scale=x_scale,
        use_leg01=use_leg01,
        leg01_frac=leg01_frac,
        prior_min_sigma=prior_min_sigma,
        use_band_snr=use_band_snr,
        band_snr_coeffs=band_snr_coeffs,
        meas_sigma=meas_sigma,
        noise_degrade=noise_degrade,
        lm=lm,
        conv=conv,
        max_outer_steps=max_outer_steps,
        n_pc=n_pc,
        n_leg=n_leg,
        log_trans=log_trans,
    )
end

"""Split SVD state into background and SIF TOA radiance (matches `make_svd_forward_model_λ`)."""
function svd_obs_components(
    x::AbstractVector{<:Real},
    λ_obs::AbstractVector{<:Real},
    solar_eff::AbstractVector{<:Real},
    PCs::AbstractMatrix{<:Real},
    sif_basis::AbstractMatrix{<:Real},
    layout;
    log_transform::Bool,
)
    c_vec = @view x[layout.idx_pc]
    alpha_coeff = alpha_coeff_from_raw(x[first(layout.idx_alpha)])
    leg_coeff = @view x[layout.idx_legendre]
    sif_coeff = @view x[layout.idx_sif]
    leg_basis = _legendre_design_matrix(_normalized_grid(λ_obs), layout.n_legendre)
    trans_up = log_transform ? exp.(PCs * c_vec) : 1.0 .+ PCs * c_vec
    trans_updown = exp.(alpha_coeff .* log.(max.(trans_up, eps(Float64))))
    rho_obs = leg_basis * leg_coeff
    sif_toa = trans_up .* (sif_basis * sif_coeff)
    y_refl = @.(solar_eff * trans_updown * rho_obs / π)
    return (y_toa=y_refl .+ sif_toa, y_refl=y_refl, sif_toa=sif_toa, sif_coeff=collect(sif_coeff))
end

"""
Run one-pixel SVD LM retrieval on noisy surrogate TOA.
Uses E from a random L1B file and SZA = 0 for `solar_eff = E / (π · esd)`.
"""
function run_pseudo_svd_retrieval(
    y_obs::AbstractVector{<:Real},
    λ_obs::AbstractVector{<:Real},
    sh,
    solar::NamedTuple;
    sza_deg::Float64=PSEUDO_RETRIEVAL_SZA,
)
    length(y_obs) == length(λ_obs) || error("y_obs / λ_obs length mismatch")
    solar_eff = @. solar.E * cosd(sza_deg) / π / solar.esd

    fm, layout = make_svd_forward_model_λ(
        λ_obs,
        solar_eff,
        sh.PCs,
        sh.sif_basis;
        n_pc=sh.n_pc,
        n_legendre=sh.n_leg,
        log_transform=sh.log_trans,
    )
    jac_eval = (x -> ForwardDiff.jacobian(fm, x))
    x_out = zeros(Float64, layout.n_state)
    x_a = copy(sh.x0)
    σ_prior = copy(sh.prior_sigma)
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

    stats = _run_one_svd_retrieval!(
        x_out,
        fm,
        jac_eval,
        layout,
        collect(Float64.(y_obs)),
        copy(sh.x0),
        x_a,
        σ_prior,
        sh.lower,
        sh.upper,
        sh.x_scale,
        sh.use_leg01,
        sh.A01,
        sh.prior_min_sigma,
        sh.leg01_frac,
        sh.use_band_snr,
        sh.band_snr_coeffs,
        sh.meas_sigma,
        sh.lm,
        sh.conv,
        sh.conv.dx_rel_tol_stall,
        sh.max_outer_steps,
    )

    y_fit = fm(x_out)
    parts = svd_obs_components(
        x_out, λ_obs, solar_eff, sh.PCs, sh.sif_basis, layout;
        log_transform=sh.log_trans,
    )
    resid = collect(Float64.(y_obs)) .- y_fit
    return (
        x=x_out,
        y_fit=y_fit,
        y_refl_fit=parts.y_refl,
        sif_toa_fit=parts.sif_toa,
        sif_coeff=parts.sif_coeff,
        resid=resid,
        solar_eff=solar_eff,
        stats=stats,
        fm=fm,
        layout=layout,
        sza_deg=sza_deg,
        solar_path=solar.pace_path,
    )
end

"""Four-panel QA: TOA, SIF@TOA, absolute residual vs. σ, relative residual (%)."""
function plot_pseudo_retrieval(ret; out_path::AbstractString, noise_σ::AbstractVector{<:Real})
    λ = ret.λ
    y_obs = ret.y_obs
    xlim = _spectral_xlim(λ)
    rel_resid_pct = @. ret.resid / max(abs(y_obs), eps(Float64)) * 100.0
    rel_sigma_pct = @. noise_σ / max(abs(y_obs), eps(Float64)) * 100.0

    m = _SUBPLOT_X_MARGINS
    xkw = (
        xlims=xlim,
        left_margin=m.left_margin,
        right_margin=m.right_margin,
        legend=:none,
    )

    p1 = plot(λ, ret.y_true_toa; xkw..., color=:black, lw=2.0, label="True TOA (clean)",
        xticks=:none, xlabel="", ylabel="Radiance", title="(a) TOA",
        top_margin=m.top_margin, bottom_margin=1Plots.mm, size=(900, 230))
    plot!(p1, λ, y_obs; color=:gray, lw=1.0, ls=:dot, alpha=0.8, label="Observed (+ noise)")
    plot!(p1, λ, ret.y_fit; color=:crimson, lw=1.5, ls=:dash, label="Fitted TOA")

    p2 = plot(λ, ret.sif_true_toa; xkw..., color=:forestgreen, lw=2.0, label="True SIF×T₁",
        xticks=:none, xlabel="", ylabel="Radiance", title="(b) SIF at TOA",
        top_margin=1Plots.mm, bottom_margin=1Plots.mm, size=(900, 230))
    plot!(p2, λ, ret.sif_toa_fit; color=:darkgreen, lw=1.5, ls=:dash, label="Retrieved SIF×T₁")

    p3 = plot(λ, ret.resid; xkw..., color=:navy, lw=1.2, label="Residual (obs − fit)",
        xticks=:none, xlabel="", ylabel="Radiance", title="(c) Residual vs. white-noise σ",
        top_margin=1Plots.mm, bottom_margin=1Plots.mm, size=(900, 230))
    plot!(p3, λ, noise_σ; color=:orange, lw=0.8, ls=:dash, label="+σ noise")
    plot!(p3, λ, .-noise_σ; color=:orange, lw=0.8, ls=:dash, label="−σ noise")

    p4 = plot(λ, rel_resid_pct; xkw..., color=:navy, lw=1.2,
        label="Relative residual (obs − fit) / obs",
        xticks=640:20:760, xlabel="Wavelength [nm]", ylabel="Rel. residual [%]",
        title="(d) Relative residual",
        top_margin=1Plots.mm, bottom_margin=8Plots.mm, size=(900, 260))
    plot!(p4, λ, rel_sigma_pct; color=:orange, lw=0.8, ls=:dash, label="+σ/y_obs (%)")
    plot!(p4, λ, .-rel_sigma_pct; color=:orange, lw=0.8, ls=:dash, label="−σ/y_obs (%)")

    p = plot(p1, p2, p3, p4;
        layout=(4, 1),
        link=:x,
        size=(1100, 1100),
        legend=:outerright,
        plot_title="Pseudo SVD retrieval (SZA=$(ret.sza_deg)°, status=$(ret.stats.status), " *
                    "RMSE=$(round(ret.stats.rmse, digits=3)), rχ²=$(round(ret.stats.reduced_chi2, digits=2))",
    )

    savefig(p, out_path)
    println("Saved pseudo-retrieval figure: ", out_path)
    return p
end

"""
    pseudo_retrieval_from_surrogate(toa, sif, λ_band, base_dir; out_path)

(E) workflow: SVD LM on noisy surrogate TOA with L1B E and SZA = 0.
"""
function pseudo_retrieval_from_surrogate(
    toa,
    sif_bundle,
    λ_band::AbstractVector{<:Real},
    base_dir::AbstractString;
    retrieval_cfg::Dict{String, Any},
    out_path::AbstractString,
)
    sh = prepare_svd_retrieval_setup(retrieval_cfg, collect(Float64.(λ_band)))
    l1b_path = pick_random_l1b_file(base_dir)
    solar = load_l1b_solar_on_bands(l1b_path, λ_band)
    println("  L1B solar file: $(basename(l1b_path))")
    println("  SZA for retrieval: $(PSEUDO_RETRIEVAL_SZA)°")
    println("  earth_sun_distance_correction: $(round(solar.esd, digits=4))")

    ret_core = run_pseudo_svd_retrieval(toa.R_noisy, λ_band, sh, solar)
    ret = merge(ret_core, (
        λ=collect(Float64.(λ_band)),
        y_obs=toa.R_noisy,
        y_true_toa=toa.R_toa,
        sif_true_toa=sif_bundle.sif.R_sif_toa,
    ))

    println("  retrieval status: $(ret.stats.status)  converged=$(ret.stats.converged)  steps=$(ret.stats.n_steps)")
    println("  RMSE=$(round(ret.stats.rmse, digits=4))  reduced χ²=$(round(ret.stats.reduced_chi2, digits=3))  DOF=$(round(ret.stats.dof, digits=1))")
    println("  mean |residual|: $(round(mean(abs.(ret.resid)), digits=4))")

    plot_pseudo_retrieval(ret; out_path=out_path, noise_σ=toa.σ)
    return ret
end

"""Write all plot spectra + profile/weight metadata for one surrogate case."""
function write_surrogate_profile_nc(
    out_path::AbstractString;
    profile_index::Int,
    prof,
    weights::AbstractVector{<:Real},
    reflectance_decay::Float64,
    reflectance_shape::AbstractString,
    amf_down::Float64,
    amf_up::Float64,
    λ_hres::AbstractVector{<:Real},
    λ_band::AbstractVector{<:Real},
    τ_2way::AbstractVector{<:Real},
    T_solar::AbstractVector{<:Real},
    T2_atm::AbstractVector{<:Real},
    T2_hres::AbstractVector{<:Real},
    T_solar_band::AbstractVector{<:Real},
    T2_atm_band::AbstractVector{<:Real},
    T2_band::AbstractVector{<:Real},
    reflectance,
    sif,
    toa,
    retrieval,
)
    n_layer = length(weights)
    n_hres = length(λ_hres)
    n_band = length(λ_band)
    length(prof.p_full) == n_layer || error("p_full / weights length mismatch")
    length(τ_2way) == n_hres || error("τ_2way / λ_hres length mismatch")
    length(T2_band) == n_band || error("T2_band / λ_band length mismatch")

    height_mid = (prof.height_km[1:end-1] .+ prof.height_km[2:end]) ./ 2
    λ_refl = collect(Float64.(reflectance.obs.λ))
    n_refl = length(λ_refl)

    isfile(out_path) && rm(out_path)
    ds = Dataset(out_path, "c")
    defDim(ds, "layer", n_layer)
    defDim(ds, "band", n_band)
    defDim(ds, "hres", n_hres)
    defDim(ds, "refl_band", n_refl)

    defVar(ds, "layer_index", Int32.(1:n_layer), ("layer",))
    defVar(ds, "pressure", Float64.(prof.p_full), ("layer",); attrib=Dict("units"=>"hPa"))
    defVar(ds, "height", Float64.(height_mid), ("layer",); attrib=Dict("units"=>"km"))
    defVar(ds, "temperature", Float64.(prof.temp), ("layer",); attrib=Dict("units"=>"K"))
    defVar(ds, "vmr_h2o", Float64.(prof.vmr_h2o), ("layer",); attrib=Dict("units"=>"mol mol-1"))
    defVar(ds, "reflectance_weight", Float64.(weights), ("layer",);
           attrib=Dict("long_name"=>"normalized layer reflectance weights"))

    defVar(ds, "wavelength_hres", Float64.(λ_hres), ("hres",); attrib=Dict("units"=>"nm"))
    defVar(ds, "tau_2way", Float64.(τ_2way), ("hres",))
    defVar(ds, "T_solar_hres", Float64.(T_solar), ("hres",);
           attrib=Dict("long_name"=>"solar continuum-normalized transmittance on hi-res grid"))
    defVar(ds, "T2_atm_hres", Float64.(T2_atm), ("hres",))
    defVar(ds, "T2_hres", Float64.(T2_hres), ("hres",);
           attrib=Dict("long_name"=>"two-way atm × solar on hi-res grid"))

    defVar(ds, "wavelength", Float64.(λ_band), ("band",); attrib=Dict("units"=>"nm"))
    defVar(ds, "T_solar", Float64.(T_solar_band), ("band",))
    defVar(ds, "T2_atm", Float64.(T2_atm_band), ("band",))
    defVar(ds, "T2", Float64.(T2_band), ("band",);
           attrib=Dict("long_name"=>"two-way atm × solar on OCI bands"))
    defVar(ds, "T1", Float64.(sif.T1.T1_band), ("band",);
           attrib=Dict("long_name"=>"one-way atmospheric transmittance on OCI bands"))
    defVar(ds, "SIF", Float64.(sif.sif.SIF), ("band",);
           attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "SIF_toa", Float64.(sif.sif.R_sif_toa), ("band",);
           attrib=Dict("long_name"=>"SIF × T1", "units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "R_cont", Float64.(toa.R_cont), ("band",);
           attrib=Dict("long_name"=>"rescaled continuum radiance", "units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "R_bg", Float64.(toa.R_bg), ("band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "R_toa_clean", Float64.(toa.R_toa), ("band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "R_toa_noisy", Float64.(toa.R_noisy), ("band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "R_toa_noisy_bg", Float64.(toa.R_noisy_bg), ("band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "noise_sigma", Float64.(toa.σ), ("band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "snr", Float64.(toa.snr), ("band",))

    defVar(ds, "R_fit", Float64.(retrieval.y_fit), ("band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "SIF_toa_fit", Float64.(retrieval.sif_toa_fit), ("band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "residual", Float64.(retrieval.resid), ("band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))

    defVar(ds, "wavelength_refl", λ_refl, ("refl_band",); attrib=Dict("units"=>"nm"))
    defVar(ds, "R_l1b", Float64.(reflectance.obs.R), ("refl_band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "R_legendre", Float64.(reflectance.fit.R_fit), ("refl_band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "R_legendre_rescaled", Float64.(reflectance.resc.R), ("refl_band",); attrib=Dict("units"=>"W m-2 sr-1 um-1"))
    defVar(ds, "rho_legendre_rescaled", Float64.(reflectance.resc.ρ), ("refl_band",))

    ds.attrib["title"] = "Surrogate measurement spectra for plotting (profile $profile_index)"
    ds.attrib["profile_index"] = profile_index
    ds.attrib["ps_hpa"] = Float64(prof.ps_hpa)
    ds.attrib["stored_amf"] = Float64(prof.amf)
    ds.attrib["sza_deg"] = Float64(SZA_DEG)
    ds.attrib["vza_deg"] = Float64(VZA_DEG)
    ds.attrib["amf_down"] = amf_down
    ds.attrib["amf_up"] = amf_up
    ds.attrib["reflectance_shape"] = String(reflectance_shape)
    ds.attrib["reflectance_decay"] = reflectance_decay
    ds.attrib["sif_library_index"] = Int(sif.sif.library_index)
    ds.attrib["sif_strength"] = Float64(sif.sif.strength)
    ds.attrib["l1b_pixel"] = Int(reflectance.obs.pixel)
    ds.attrib["l1b_scan"] = Int(reflectance.obs.scan)
    ds.attrib["retrieval_status"] = Int(retrieval.stats.status)
    ds.attrib["retrieval_converged"] = Int(retrieval.stats.converged)
    ds.attrib["retrieval_rmse"] = Float64(retrieval.stats.rmse)
    ds.attrib["retrieval_reduced_chi2"] = Float64(retrieval.stats.reduced_chi2)
    ds.attrib["retrieval_dof"] = Float64(retrieval.stats.dof)
    ds.attrib["created"] = string(Dates.now())
    close(ds)
    println("Saved spectra NetCDF: ", out_path)
    return out_path
end

function main()
    mkpath(OUTPUT_DIR)
    Random.seed!(RANDOM_SEED)

    ENV["PACE_LUT_INTERPOLATION"] = LUT_INTERPOLATION
    ctx = MWEF.prepare_mwe_inputs(CONFIG_PATH)
    cfg = TOML.parsefile(CONFIG_PATH)
    data_cfg = get(cfg, "data", Dict{String, Any}())
    solar_file = get(data_cfg, "solar_file", "solar_merged_20240731_600_33300_100.out")
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)

    ds = Dataset(TRANS_NC)
    n_profiles = ds.dim["profile"]
    profile_index = rand(1:n_profiles)
    prof = load_profile(ds, profile_index)
    close(ds)

    println("Selected profile $profile_index / $n_profiles")
    println("  Surface pressure: $(prof.ps_hpa) hPa")
    println("  Stored AMF: $(prof.amf)")

    weights = layer_reflectance_weights(
        length(prof.temp);
        shape=REFLECTANCE_SHAPE,
        decay=REFLECTANCE_DECAY,
    )
    println("  Reflectance weights — surface layer: $(round(weights[end], digits=4)), top layer: $(round(weights[1], digits=4))")

    profile_fig = joinpath(OUTPUT_DIR, "profile_$(profile_index).png")
    plot_profile(prof, weights; out_path=profile_fig)

    τ_layer, τ_cum = layer_optical_depth(
        ctx.spectral_axis,
        prof.p_full,
        prof.temp,
        prof.vcd_dry,
        prof.vcd_h2o,
        ctx.o2_sitp,
        ctx.h2o_sitp,
    )

    amf_down = 1.0 / cosd(SZA_DEG)
    amf_up = 1.0 / cosd(VZA_DEG)
    τ_2way = weighted_two_way_optical_depth(τ_cum, weights, amf_down, amf_up)

    T2_atm = exp.(-τ_2way)

    T_solar, _ = MWEF.load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
    T_solar = Float64.(T_solar)
    # Normalize solar to unit continuum peak so it overlays cleanly with transmittance
    T_solar_peak = maximum(T_solar)
    T_solar_peak > 0 || error("Solar spectrum is all zero")
    T_solar = T_solar ./ T_solar_peak
    T2_hres = T2_atm .* T_solar

    K = ctx.kernel_rsr_out
    λ_band = collect(Float64.(ctx.λ))
    T2_band = vec(K * T2_hres)
    T_solar_band = vec(K * T_solar)
    T2_atm_band = vec(K * T2_atm)

    trans_fig = joinpath(OUTPUT_DIR, "transmittance_profile$(profile_index).png")
    plot_transmittance(
        ctx.λ_hres,
        λ_band,
        T_solar,
        T2_hres,
        T2_band;
        out_path=trans_fig,
        T_solar_band=T_solar_band,
        T2_atm_band=T2_atm_band,
    )

    summary_path = joinpath(OUTPUT_DIR, "profile_$(profile_index)_summary.txt")
    open(summary_path, "w") do io
        println(io, "profile_index\t$(profile_index)")
        println(io, "ps_hpa\t$(prof.ps_hpa)")
        println(io, "stored_amf\t$(prof.amf)")
        println(io, "sza_deg\t$(SZA_DEG)")
        println(io, "vza_deg\t$(VZA_DEG)")
        println(io, "amf_down\t$(amf_down)")
        println(io, "amf_up\t$(amf_up)")
        println(io, "reflectance_shape\t$(REFLECTANCE_SHAPE)")
        println(io, "reflectance_decay\t$(REFLECTANCE_DECAY)")
        println(io, "tau_2way_min\t$(minimum(τ_2way))")
        println(io, "tau_2way_max\t$(maximum(τ_2way))")
        println(io, "T2_hres_min\t$(minimum(T2_hres))")
        println(io, "T2_hres_max\t$(maximum(T2_hres))")
    end
    println("Saved summary: ", summary_path)

    # ── (B) Rebuild reflectance from a random L1B pixel ───────────────────────
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    pace_file = get(pace_cfg, "pace_file", "sample_granule_20240830T131442_new_chl.nc")
    pace_path = isabspath(pace_file) ? pace_file : joinpath(ctx.paths.base_dir, pace_file)
    λ_min = Float64(get(get(cfg, "spectral", Dict()), "lambda_min_nm", 640.0))
    λ_max = Float64(get(get(cfg, "spectral", Dict()), "lambda_max_nm", 756.0))

    println("\n=== (B) Rebuild reflectance ===")
    refl_fig = joinpath(OUTPUT_DIR, "reflectance_profile$(profile_index).png")
    refl = rebuild_reflectance_from_l1b(
        pace_path;
        order=LEGENDRE_ORDER,
        out_path=refl_fig,
        mean_λ_min=λ_min,   # only for mean-radiance QC / rescale target
        mean_λ_max=λ_max,
    )

    # ── (C) SIF + one-way T₁ ───────────────────────────────────────────────────
    println("\n=== (C) Add SIF ===")
    sif_fig = joinpath(OUTPUT_DIR, "sif_profile$(profile_index).png")
    sif = add_sif_component(
        τ_layer,
        amf_up,
        K,
        ctx.paths.sif_path,
        λ_band;
        out_path=sif_fig,
    )

    # ── (D) Background + SIF → TOA + SNR white noise ─────────────────────────
    println("\n=== (D) Rebuild TOA + white noise ===")
    toa_fig = joinpath(OUTPUT_DIR, "toa_profile$(profile_index).png")
    toa = rebuild_toa_with_noise(
        λ_band,
        T2_band,
        refl,
        sif;
        band_snr_coeffs=ctx.band_snr_coeffs,
        out_path=toa_fig,
    )

    open(summary_path, "a") do io
        println(io, "sif_library_index\t$(sif.sif.library_index)")
        println(io, "sif_strength\t$(sif.sif.strength)")
        println(io, "sif_mean\t$(mean(sif.sif.SIF))")
        println(io, "sif_toa_mean\t$(mean(sif.sif.R_sif_toa))")
        println(io, "T1_band_mean\t$(mean(sif.T1.T1_band))")
        println(io, "toa_mean_bg\t$(mean(toa.R_bg))")
        println(io, "toa_mean_clean\t$(mean(toa.R_toa))")
        println(io, "toa_mean_noisy\t$(mean(toa.R_noisy))")
        println(io, "toa_mean_noisy_bg\t$(mean(toa.R_noisy_bg))")
        println(io, "toa_median_snr\t$(median(filter(isfinite, toa.snr)))")
        println(io, "toa_sigma_min\t$(minimum(toa.σ))")
        println(io, "toa_sigma_max\t$(maximum(toa.σ))")
        println(io, "l1b_pixel\t$(refl.obs.pixel)")
        println(io, "l1b_scan\t$(refl.obs.scan)")
    end

    # ── (E) Pseudo SVD retrieval ───────────────────────────────────────────────
    println("\n=== (E) Pseudo SVD retrieval ===")
    isfile(SVD_RETRIEVAL_CONFIG) || error("SVD retrieval config not found: $SVD_RETRIEVAL_CONFIG")
    svd_cfg = TOML.parsefile(SVD_RETRIEVAL_CONFIG)
    ret_fig = joinpath(OUTPUT_DIR, "retrieval_profile$(profile_index).png")
    retrieval = pseudo_retrieval_from_surrogate(
        toa,
        sif,
        λ_band,
        ctx.paths.base_dir;
        retrieval_cfg=svd_cfg,
        out_path=ret_fig,
    )

    open(summary_path, "a") do io
        println(io, "retrieval_l1b_solar\t$(basename(retrieval.solar_path))")
        println(io, "retrieval_sza_deg\t$(retrieval.sza_deg)")
        println(io, "retrieval_status\t$(retrieval.stats.status)")
        println(io, "retrieval_converged\t$(retrieval.stats.converged)")
        println(io, "retrieval_rmse\t$(retrieval.stats.rmse)")
        println(io, "retrieval_reduced_chi2\t$(retrieval.stats.reduced_chi2)")
        println(io, "retrieval_dof\t$(retrieval.stats.dof)")
    end

    nc_path = joinpath(OUTPUT_DIR, "profile_$(profile_index)_spectra.nc")
    write_surrogate_profile_nc(
        nc_path;
        profile_index=profile_index,
        prof=prof,
        weights=weights,
        reflectance_decay=REFLECTANCE_DECAY,
        reflectance_shape=String(REFLECTANCE_SHAPE),
        amf_down=amf_down,
        amf_up=amf_up,
        λ_hres=collect(Float64.(ctx.λ_hres)),
        λ_band=λ_band,
        τ_2way=τ_2way,
        T_solar=T_solar,
        T2_atm=T2_atm,
        T2_hres=T2_hres,
        T_solar_band=T_solar_band,
        T2_atm_band=T2_atm_band,
        T2_band=T2_band,
        reflectance=refl,
        sif=sif,
        toa=toa,
        retrieval=retrieval,
    )

    return (
        profile=prof,
        weights=weights,
        λ_hres=ctx.λ_hres,
        λ_band=λ_band,
        τ_layer=τ_layer,
        τ_2way=τ_2way,
        T2_hres=T2_hres,
        T2_band=T2_band,
        T_solar=T_solar,
        T2_atm=T2_atm,
        reflectance=refl,
        sif=sif,
        toa=toa,
        retrieval=retrieval,
        spectra_nc=nc_path,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
