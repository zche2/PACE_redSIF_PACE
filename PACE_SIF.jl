module PACE_SIF

    # use toolboxes
    include("./tools/Instrument.jl")
    include("./tools/Atmosphere.jl")
    include("./tools/Retrieval.jl")

    # Re-export everything from submodules
    using .Retrieval
    export Retrieval
    
end