using Parameters
using LinearAlgebra

@with_kw struct SpectraOfPC{FT <: AbstractFloat} 
    band::Vector{FT}                            # Wavelengths (1D array)
    PrinComp::Union{LinearAlgebra.Adjoint{Float64, Matrix{Float64}}, Matrix{FT}}                        
                                                # Columns for Principal Components
    VarExp::Vector{FT}                          # Variance explained (1D array)
    Loading::Union{LinearAlgebra.Adjoint{Float64, Matrix{Float64}}, Matrix{FT}}    
                                                # Rows for each principal components
    if_log::Bool = false                        # Default value for boolean marker
end

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

function center_wavelength(
        λ::Vector{Union{Missing, FT}}
    ) where {FT <: AbstractFloat}
	# get the range and medium of λ and center it to [0,1]
	λ_max = ceil(maximum(λ));
	λ_min = floor(minimum(λ));
	range = (λ_max - λ_min) / 2;
	λ_middle = (λ_max + λ_min) / 2;
	λc       = (λ .- λ_middle) ./ range;
	return λc
end

function root_mean_square(
        y_obs, 
        y_retrieve
    )
    n     = length(y_obs)
	totSQ = sum( ( y_obs .- y_retrieve ) .^ 2 );
	RMS   = sqrt( totSQ / n );
	return RMS
end

function scale_transmittance(
    T,   # ::Vector{FT},
    ind::Vector{Int64},   # baseline indices
) where {FT <: AbstractFloat}
    # rescale transmittance term to keep the maximum value at 1
    T_abs = abs.(T)
    # find max
    bl    = maximum(T_abs[ind]);
    # force the mean val to be 1
    T_norm = T_abs ./ bl
    return T_norm
end