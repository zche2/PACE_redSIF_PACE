# Band-wise measurement uncertainty inflation for surrogate ensemble experiments.
# Config: [fit.noise_degrade] in ENSEMBLE_CONFIG, or ENV DEGRADE_NOISE=1, etc.

"""Parse `[fit.noise_degrade]` (optional ENV overrides). Returns `nothing` if disabled."""
function parse_noise_degrade(cfg::Dict)
    fit = get(cfg, "fit", Dict{String, Any}())
    deg = get(fit, "noise_degrade", Dict{String, Any}())
    enabled = Bool(get(deg, "enabled", false))
    if haskey(ENV, "DEGRADE_NOISE")
        enabled = lowercase(get(ENV, "DEGRADE_NOISE", "0")) in ("1", "true", "yes")
    end
    enabled || return nothing
    λ_lo = Float64(get(deg, "lambda_min_nm", parse(Float64, get(ENV, "DEGRADE_LAMBDA_MIN", "683"))))
    λ_hi = Float64(get(deg, "lambda_max_nm", parse(Float64, get(ENV, "DEGRADE_LAMBDA_MAX", "695"))))
    factor = Float64(get(deg, "sigma_factor", parse(Float64, get(ENV, "DEGRADE_SIGMA_FACTOR", "100"))))
    λ_hi > λ_lo || error("noise_degrade: lambda_max_nm must be > lambda_min_nm")
    factor > 0 || error("noise_degrade: sigma_factor must be > 0")
    return (lambda_min_nm=λ_lo, lambda_max_nm=λ_hi, sigma_factor=factor)
end

function noise_degrade_band_mask(λ::AbstractVector{<:Real}, degrade)
    findall(b -> degrade.lambda_min_nm <= Float64(λ[b]) <= degrade.lambda_max_nm, eachindex(λ))
end

"""Multiply measurement σ on degraded bands by `sigma_factor`."""
function apply_noise_degrade_to_sigma!(
    σ::AbstractVector{<:Real},
    λ::AbstractVector{<:Real},
    degrade,
)
    f = degrade.sigma_factor
    for b in noise_degrade_band_mask(λ, degrade)
        σ[b] *= f
    end
    return σ
end

"""Inflate SNR variance coeffs (c1, c2) by σ_factor² so LM downweights the same bands."""
function apply_noise_degrade_to_snr_coeffs!(
    c1::AbstractVector{<:Real},
    c2::AbstractVector{<:Real},
    λ::AbstractVector{<:Real},
    degrade,
)
    f2 = degrade.sigma_factor^2
    for b in noise_degrade_band_mask(λ, degrade)
        c1[b] *= f2
        c2[b] *= f2
    end
    return c1, c2
end

function copy_snr_coeffs_with_degrade(band_snr_coeffs::Dict, λ::AbstractVector{<:Real}, degrade)
    c1 = copy(Float64.(band_snr_coeffs["c1"]))
    c2 = copy(Float64.(band_snr_coeffs["c2"]))
    apply_noise_degrade_to_snr_coeffs!(c1, c2, λ, degrade)
    return Dict("c1" => c1, "c2" => c2)
end

function _noise_degrade_truth_attrs(degrade)
    if degrade === nothing
        return Dict{String, Any}("noise_degrade" => 0)
    end
    return Dict{String, Any}(
        "noise_degrade" => 1,
        "noise_degrade_lambda_min_nm" => degrade.lambda_min_nm,
        "noise_degrade_lambda_max_nm" => degrade.lambda_max_nm,
        "noise_degrade_sigma_factor" => degrade.sigma_factor,
    )
end

function _truth_noise_degrade_matches(truth_nc::AbstractString, degrade)
    ds = Dataset(truth_nc)
    on = get(ds.attrib, "noise_degrade", 0) == 1
    close(ds)
    if degrade === nothing
        return !on
    end
    return on &&
        Float64(get(ds.attrib, "noise_degrade_lambda_min_nm", NaN)) == degrade.lambda_min_nm &&
        Float64(get(ds.attrib, "noise_degrade_lambda_max_nm", NaN)) == degrade.lambda_max_nm &&
        Float64(get(ds.attrib, "noise_degrade_sigma_factor", NaN)) == degrade.sigma_factor
end
