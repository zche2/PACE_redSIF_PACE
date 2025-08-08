using Interpolations
using Parameters

function interpolate_RSR(band, wvlen, RSR, wvlen_out)
    # interpolate RSR to this resolution
    # RSR matrix: m x n
    # m - number of bands; n - number of wavelength in reference spectrum
    RSR_out = zeros(Float64, length(band), length(wvlen_out));
    for i=1:length(band)
        # interpolator: knots & value
        interp_linear = 
            LinearInterpolation(wvlen, RSR[:,i], extrapolation_bc=Line());
        # evaluate at wvlen_out
        RSR_out[i,:]  = interp_linear(wvlen_out);
    end
    # normalize
    RSR_out_norm = RSR_out ./ sum(RSR_out, dims=2);
    return RSR_out_norm
end;

    
@with_kw struct KernelInstrument
    band    # ::Array{Any,1}
    wvlen   # ::Array{FT,1}
    RSR     # ::Matrix{FT}
    wvlen_out  # ::Array{FT,1}
    RSR_out = interpolate_RSR(band, wvlen, RSR, wvlen_out)  # ::Matrix{FT} 
end;


function conv_matx(
    m::KernelInstrument,
    spectrum
    )::Vector{Float64}       # specify the return type
    # check the length of spectrum, make sure the dimensions are matching
    try 
        spec_conv = m.RSR_out * spectrum
        return spec_conv
    catch err
        println("error emerges when convolving the spectrum: ", err)
    end;
end;

