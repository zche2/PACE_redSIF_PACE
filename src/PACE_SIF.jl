module PACE_SIF

# Include toolboxes (each should be a module)
include("./tools/Instrument.jl")
include("./tools/Atmosphere.jl")
include("./tools/Retrieval/Retrieval.jl")

# Make submodules accessible
using .Instrument
using .Retrieval
using .Atmosphere

# Export submodules so users can access them
export Instrument, Atmosphere, Retrieval

# Optionally re-export specific functions/types
# Instrument
export KernelInstrument, conv_matx, interpolate_RSR, read_rescale, λ_to_ν, ν_to_λ

# Atmosphere
export layer_VCD

# Retrieval
# export types
export SpectraOfPC, MatrixFactor
export RetrievalParams_PCFit, Pixel_PCFit         # PC fit
export RetrievalParams_xSecFit, Pixel_xSecFit     # cross-section fit

# export forward model (can choose which to export)
export forward_model, compute_transmittance

# export decomposition methods
export Spectral_SVD, Spectral_NMF

# export other tools
export center_wavelength, root_mean_square, scale_transmittance, 
       Jacobian, GainMatrix, MakePriori!
export Retrieval_for_Pixel

# export iteration
export GN_Iteration!, LM_Iteration!

end  # module PACE_SIF