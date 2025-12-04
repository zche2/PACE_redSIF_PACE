using LinearAlgebra, Interpolations
using NMF

function Spectral_SVD(
    profile::Matrix{FT},     # [wavelength x sample]
    bandᵢₙ::Vector{FT},
    bandₒᵤₜ::Vector{FT};
    if_log::Bool = false,
    ) where {FT <: AbstractFloat}

    n_sample    = size(profile, 2)
    profile_new = zeros(FT, length(bandₒᵤₜ), n_sample)   # [wavelength x sample]

    for i in 1:n_sample
        itpᵢ = LinearInterpolation(bandᵢₙ, profile[:, i], extrapolation_bc=0)
        profile_new[:, i] = itpᵢ.(bandₒᵤₜ)
    end

    print("Spectra interpolated to target bands: from $(length(bandᵢₙ)) to $(length(bandₒᵤₜ)).\n")

    # --- SVD ---
    if if_log
        profile_new = - log.(profile_new)
    end
    
    F        = svd(profile_new);
    PrinComp = F.U;    
    S        = F.S;  
    Vt       = F.Vt; 
    VarExp   = S ./ sum(S) * 100;
    Loading  = diagm(S) * Vt;

    println("Shape of PrinComp: $(size(PrinComp))\n" *
            "Shape of VarExp: $(size(S))\n"*
            "Shape of Loading (S x V'): $(size(Vt))"
    )

    # return as a struct
    return SpectraOfPC(
        band     = bandₒᵤₜ,
        PrinComp = PrinComp,
        VarExp   = VarExp,
        Loading  = Loading,
        if_log   = if_log
    )
end


function Spectral_NMF(
    profile::Matrix{FT},
    bandᵢₙ::Vector{FT},
    bandₒᵤₜ::Vector{FT};
    rank::Int = 10,
    if_log::Bool = false,
    ) where {FT <: AbstractFloat}

    # --- Interpolation to align with bandₒᵤₜ ---
    profile_new = zeros(FT, size(profile, 1), length(bandₒᵤₜ))  # Specify type FT
    
    for i in 1:size(profile, 1)
        itpᵢ = LinearInterpolation(bandᵢₙ, profile[i, :], extrapolation_bc=0)
        profile_new[i, :] = itpᵢ.(bandₒᵤₜ)
    end

    println("the shape of profile matrix is $(size(profile_new)) - the second dimension should be wavelength!")

    # --- NMF ---
    if if_log
        profile_new = - log.(profile_new)
    end

    W, H = NMF.spa(profile_new, rank)
    NMF.solve!(NMF.SPA{FT}(obj=:mse), profile_new, W, H)

    # return as a struct
    return MatrixFactor{FT}(
        band = bandₒᵤₜ,
        PrinComp = H,
        Loading = W
    )
    
end


function Spectral_SVD_discard(
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