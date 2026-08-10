# SIF basis loader for SVD retrieval (no PACE_SIF / xSecFit dependency).

using Interpolations
using JLD2

"""Validate that a file exists and return the original path."""
function must_exist(path::AbstractString)
    isfile(path) || error("Missing required file: $path")
    return path
end

"""
    load_sif_basis(sif_path, λ; nEV=2, normalize=true)

Interpolate `SIF_U[:, 1:nEV]` onto wavelength grid `λ`.
Uses cubic spline interpolation when `SIF_wavelen` is evenly spaced,
and falls back to linear interpolation otherwise.
Returns matrix with size `(length(λ), nEV)`.
"""
function load_sif_basis(
    sif_path::AbstractString,
    λ::AbstractVector{<:Real};
    nEV::Int = 2,
    normalize::Bool = true,
)
    nEV < 1 && error("nEV must be >= 1")
    sif = JLD2.load(must_exist(sif_path))

    sif_u = convert.(Float64, sif["SIF_U"])
    λ_ref = collect(Float64.(sif["SIF_wavelen"]))

    size(sif_u, 1) == length(λ_ref) ||
        error("SIF_U first dimension must match SIF_wavelen length")

    n_available = size(sif_u, 2)
    n_use = min(nEV, n_available)
    n_use < nEV && @warn "Requested nEV=$nEV but only $n_available available; using $n_use"

    dλ = diff(λ_ref)
    step = dλ[1]
    tol = max(1e-10, abs(step) * 1e-8)
    use_cubic = all(abs.(dλ .- step) .<= tol)
    if use_cubic
        λ_knots = range(λ_ref[1], step = step, length = length(λ_ref))
    else
        @warn "SIF_wavelen is not evenly spaced; using linear interpolation."
    end

    basis = zeros(Float64, length(λ), n_use)
    for iev in 1:n_use
        if use_cubic
            itp_ev = CubicSplineInterpolation(λ_knots, sif_u[:, iev]; extrapolation_bc = Line())
            vals = itp_ev.(λ)
            vals[(λ .< λ_ref[1]) .| (λ .> λ_ref[end])] .= 0.0
            basis[:, iev] .= vals
        else
            itp_ev = LinearInterpolation(λ_ref, sif_u[:, iev], extrapolation_bc = 0.0)
            basis[:, iev] .= itp_ev.(λ)
        end
    end

    if normalize
        s = maximum(abs.(basis[:, 1]))
        s > 0 && (basis ./= s)
    end
    return basis
end
