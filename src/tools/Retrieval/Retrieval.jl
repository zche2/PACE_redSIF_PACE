module Retrieval

using Parameters

include("types.jl")
include("Forward_model.jl")
include("PrinComp.jl")
include("Retrieval_helper.jl")
include("Iteration.jl")

# export types
export SpectraOfPC, MatrixFactor
export RetrievalParams_PCFit, Pixel_PCFit         # PC fit
export RetrievalParams_xSecFit, Pixel_xSecFit     # cross-section fit

# export forward model (can choose which to export)
export forward_model

# export decomposition methods
export Spectral_SVD, Spectral_NMF

# export other tools
export center_wavelength, root_mean_square, scale_transmittance, 
       Jacobian, GainMatrix, MakePriori!
export Retrieval_for_Pixel

# export iteration
export GN_Iteration!, LM_Iteration!

end # module end