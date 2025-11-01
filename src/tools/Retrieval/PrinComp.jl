using LinearAlgebra, Interpolations
using NMF

function Spectral_SVD(
    profile::Matrix{FT},
    band::Vector{FT};
    λ_min::FT = 620.,
    λ_max::FT = 860.,
    if_log::Bool = false,
    ) where {FT <: AbstractFloat}

    # --- select fitting window ---
    ind    = findall( λ_min .< band .< λ_max );
    subset = profile[:, ind];
    println("the shape of profile matrix is $(size(subset)) - the second dimension should be wavelength!");

    # --- svd ---
    a = if_log ? log.(subset) : subset
    F = svd(a')
    U = F.U       # Left singular vectors
    S = F.S       # Singular values (1D vector)
    Vt = F.Vt     # Right singular vectors 
    S_norm = S ./ sum(S) * 100

    # return as a struct
    return SpectraOfPC(
        band     = band[ind],
        PrinComp = U,
        VarExp   = S_norm,
        Loading  = Vt,
        if_log   = if_log
    )
    
end

function Spectral_NMF(
    profile::Matrix{FT},
    bandᵢₙ::Vector{FT},
    bandₒᵤₜ::Vector{FT};
    rank::Int = 10,
    ) where {FT <: AbstractFloat}

    # --- Interpolation to align with bandₒᵤₜ ---
    profile_new = zeros(FT, size(profile, 1), length(bandₒᵤₜ))  # Specify type FT
    
    for i in 1:size(profile, 1)
        itpᵢ = LinearInterpolation(bandᵢₙ, profile[i, :], extrapolation_bc=0)
        profile_new[i, :] = itpᵢ.(bandₒᵤₜ)
    end

    println("the shape of profile matrix is $(size(profile_new)) - the second dimension should be wavelength!")

    # --- NMF ---
    W, H = NMF.spa(profile_new, rank)
    NMF.solve!(NMF.SPA{FT}(obj=:mse), profile_new, W, H)

    # return as a struct
    return MatrixFactor{FT}(
        band = bandₒᵤₜ,
        PrinComp = H,
        Loading = W
    )
    
end