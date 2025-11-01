using LinearAlgebra

Base.@kwdef struct SpectraOfPC{FT <: AbstractFloat} 
    band::Vector{FT}                            # Wavelengths (1D array)
    PrinComp::Union{LinearAlgebra.Adjoint{Float64, Matrix{Float64}}, Matrix{FT}}                        
                                                # Columns for Principal Components
    VarExp::Vector{FT}                          # Variance explained (1D array)
    Loading::Union{LinearAlgebra.Adjoint{Float64, Matrix{Float64}}, Matrix{FT}}    
                                                # Rows for each principal components
    if_log::Bool = false                        # Default value for boolean marker
end

Base.@kwdef struct MatrixFactor{FT <: AbstractFloat} 
    band::Vector{FT}                            # Wavelengths (1D array)
    PrinComp::Matrix{FT}                        # Spectral components (rank × wavelength)
    Loading::Matrix{FT}                         # Loading coefficients (samples × rank)
end

Base.@kwdef struct RetrievalParams
	# specific to measurement
	λ
	λc
	λ_bl_ind
	E
	c₁
	c₂

	# Forward model settings
	forward_model
	nPoly::Int
	nPC::Int
	nSIF::Int
	Sₐ
	βₐ
	PrinComp
	SIFComp

	# Iteration settings
	iteration_method = LM_Iteration!
	nIter::Int = 25
	thr_Converge::Float64 = 1e-6
end

mutable struct Pixel
	# universal for the granule
	λ      # fitting window
	E      # observed extraterrestrial irradiance
	nPoly  # degree of Legendre Polynomial
	nPC    # number of transmittance basis function included
	nSIF   # number of SIF PCs
	"a priori matrix (variance=standard deviation**2)"
	Sₐ    
	"PCs specified, = HighResSVD.PrinComp[:, 1:nPC]"
	trans_mat 
	"SIF shape specified"
	SIF_shape
	"centered wavelength for fast computation of Legendre Polys, = center_wavelength(λ)"
	λc      
	"baseline band for scaling transmittance, = find_baseline_band(λ)"
	λ_bl_ind

	# pixel L1B measurement & set up
	"TOA radiance"
	R_toa
	"solar zenith angle"
	sza
	"viewing zenith angle"
	vza
	"normalized fluorescence height (nFLH)"
	nflh
	"chlor_a concentration"
	chlor_a
	"measurement error (variance)"
	Sₑ
	"flag: 0 - not doing retrieval due to bad input data, refer to l2flag / nflh"
	flag   # 🔴 not necessarily need?
	"a priori estimation"
	xₐ
	
	# retrieval
	"retrieved state vector"
	x
	"modelled radiance"
	y
	"convergence flag"
	ΔRMSE
	"iteration label"
	iter_label

	# Inner constructer
	function Pixel()
        new()
    end
    """
    use @assert all for bounding checking
    specify type when doing outer constructor
    string representation using base show
    """
end
