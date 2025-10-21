module PACE_SIF

# Include toolboxes (each should be a module)
include("./tools/Instrument.jl")
include("./tools/Atmosphere.jl")
include("./tools/Retrieval.jl")

# Make submodules accessible
using .Instrument
using .Retrieval
using .Atmosphere

# Export submodules so users can access them
export Instrument, Atmosphere, Retrieval

# Optionally re-export specific functions/types
export KernelInstrument, conv_matx, interpolate_RSR
export Spectral_SVD, Spectral_NMF, SpectraOfPC, MatrixFactor
export center_wavelength, root_mean_square, scale_transmittance
export layer_VCD

end  # module PACE_SIF