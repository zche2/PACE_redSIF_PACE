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
export KernelInstrument, conv_matx, interpolate_RSR
export layer_VCD
export Pixel, RetrievalParams, forward_model, Retrieval_for_Pixel
export GN_Iteration!, LM_Iteration!
export Spectral_SVD, Spectral_NMF, SpectraOfPC, MatrixFactor
export center_wavelength, root_mean_square, scale_transmittance

end  # module PACE_SIF