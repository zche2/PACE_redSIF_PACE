module Retrieval

using Parameters

include("types.jl")
include("PrinComp.jl")
include("Retrieval_helper.jl")
include("Iteration.jl")

# export types
export SpectraOfPC, MatrixFactor, Pixel

# export decomposition methods
export Spectral_SVD, Spectral_NMF

# export other tools
export center_wavelength, root_mean_square, scale_transmittance, Jacobian, GainMatrix, MakePriori!

# export iteration
export Basic_Iteration!, LM_Iteration!

end # module end