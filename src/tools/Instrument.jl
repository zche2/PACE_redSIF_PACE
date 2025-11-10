module Instrument

using Interpolations
using Parameters
using vSmartMOM, vSmartMOM.Absorption

export KernelInstrument, conv_matx, interpolate_RSR, read_rescale
  
struct KernelInstrument
    band          # central wavelength of each band
    wvlen         # wavelength grid of reference spectrum (OCI)
    RSR           # spectral response function matrix
    wvlen_out     # output wavelength grid (spectrum grid), high resolution
    ν_step        # wavenumber step
    ν_grid        # wavenumber grid
    RSR_out       # interpolated RSR
    
    function KernelInstrument(band, wvlen, RSR, wvlen_out, ν_step)
        RSR_out = interpolate_RSR(band, wvlen, RSR, wvlen_out);
        λ_max = maximum(wvlen_out);
        λ_min = minimum(wvlen_out);
        ν_min = λ_to_ν(λ_max);
        ν_max = λ_to_ν(λ_min);
        ν_grid = ν_min:ν_step:ν_max;
        new(band, wvlen, RSR, wvlen_out, ν_step, ν_grid, RSR_out)
    end
end

# @with_kw struct KernelInstrument
#     band
#     wvlen
#     RSR
#     wvlen_out
#     RSR_out = interpolate_RSR(band, wvlen, RSR, wvlen_out)
#     wavnum_out = @. 1 / (wvlen_out * 1e-7)  # use broadcast macro
# end


function λ_to_ν(λ)
    # convert wavelength (in nm) to wavenumber (in cm-1)
    ν = 1 /(λ * 1E-7);
    return ν
end;


function ν_to_λ(ν)
    # convert wavelength (in nm) to wavenumber (in cm-1)
    λ = 1 /(ν * 1E7);
    return λ
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


function read_rescale(itp_filename::String)
	model = load_interpolation_model(itp_filename);
	ν_grid = model.ν_grid;
	p_grid = model.p_grid;
	t_grid = model.t_grid;
	itp    = model.itp;
	sitp = scale(itp, ν_grid, p_grid, t_grid);
	println("scaled! $itp_filename")
	return sitp
end;


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

end