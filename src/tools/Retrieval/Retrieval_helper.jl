
using LinearAlgebra

"""
    center_wavelength(λ)

Normalize wavelength array to [-1, 1] range centered at midpoint.

Returns: λc (normalized wavelengths)
"""
function center_wavelength(
    λ::Vector{Union{Missing, FT}}
) where {FT <: AbstractFloat}
    
    λ_max = ceil(maximum(λ))
    λ_min = floor(minimum(λ))
    range = (λ_max - λ_min) / 2
    λ_middle = (λ_max + λ_min) / 2
    λc = (λ .- λ_middle) ./ range
    
    return λc
end


"""
root_mean_square(y_obs, y_retrieve)

Calculate root mean square error between observed and retrieved values.

Returns: RMS error
"""
function root_mean_square(y_obs, y_retrieve)
    
    n = length(y_obs)
    totSQ = sum((y_obs .- y_retrieve) .^ 2)
    RMS = sqrt(totSQ / n)
    
    return RMS
end


"""
    scale_transmittance(T, ind)

Normalize transmittance spectrum so maximum at baseline indices equals 1.

# Arguments
- `T`: Transmittance spectrum
- `ind`: Baseline wavelength indices

Returns: T_norm (normalized transmittance)
"""
function scale_transmittance(
    T,
    ind::Vector{Int64}
) where {FT <: AbstractFloat}
    
    T_abs = abs.(T)
    bl = maximum(T_abs[ind])
    T_norm = T_abs ./ bl
    
    return T_norm
end

function Jacobian(
        x, 
        model, 
        len ::Int  # length of measured spectrum
    )
	res = DiffResults.JacobianResult(zeros(len), x);
	ForwardDiff.jacobian!(res, model, x);
	K   = DiffResults.jacobian(res);
	val = DiffResults.value(res);
	return K, val
end

function GainMatrix(
        K,
        px :: Pixel  # pixel struct
    )
    inv_Sₑ = inv(px.Sₑ);
    inv_Sₐ = inv(px.Sₐ);
	return inv( K' * inv_Sₑ * K + inv_Sₐ )K' * inv_Sₑ;
end

function MakePriori!(
    px :: Pixel,
    β   :: Vector{FT};
    γ   :: Vector{FT} = [(secd(px.sza) + secd(px.vza)) / secd(px.vza)],  
                # extra terms for T₁ and T₂ conversion
    SIF :: Vector{FT} = [px.nflh],  
                # coeff. for SIF PCs
    ) where {FT <: AbstractFloat}

    # polynomial terms
    K₀  = px.E .* cosd(px.sza) ./ pi .* hcat(collectPl.(px.λc, lmax=px.nPoly)...)';
    G₀  = inv( K₀'K₀ )K₀';
    x₀  = G₀ * px.R_toa; 

    px.xₐ = [x₀... β... γ... SIF...]';

    return nothing
end
