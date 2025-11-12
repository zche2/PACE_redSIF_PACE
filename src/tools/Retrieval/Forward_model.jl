using LegendrePolynomials

function forward_model(
        x,
        px;
        return_components::Bool=false
    )
    if px isa Pixel_xSecFit
		return forward_model_xSecFit(
            x,
            px;
            return_components=return_components
        )
	elseif px isa Pixel_PCFit
		return forward_model_PCfit(
            x,
            px;
            return_components=return_components
        )
	else
		error("Unknown retrieval method type.")
	end

end

"""
    forward_model_xSecFit(x, px::Pixel; return_components::Bool=false)

Notes:
- fit one- and two-way transmittance respectively, assuming completely independent.
"""

function forward_model_xSecFit(
        x,
        px :: Pixel_xSecFit;       # Pixel struct
        return_components::Bool=false
    )

    # reflectance
    xᵨ    = x[1 : px.nPoly+1]
    v     = collectPl.(px.λc, lmax=px.nPoly);
    ρ     = hcat(v...)' * xᵨ;

    # T↑ transmittance for SIF
    x₁    = x[(px.nPoly+2):(px.nPoly+7)];
    T₁    = compute_transmittance(
                x₁,
                px.InstrumentKernel,
                px.o2_sitp,
                px.h2o_sitp,
            );

    # T↓↑ transmittance for reflected radiance
    x₂    = x[(px.nPoly+8):(px.nPoly+13)];
    T₂    = compute_transmittance(
                x₂,
                px.InstrumentKernel,
                px.o2_sitp,
                px.h2o_sitp,
            );

    # SIF magnitude
    xₛ     = x[end - px.nSIF + 1 : end];
    SIF    = px.SIF_shape * xₛ;

    # TOA radiance
    rad   = @. px.E * cosd(px.sza) / π * T₂ * ρ + SIF * T₁;
    
    if return_components
        return rad, ρ, T₁, T₂, SIF
    else
        return rad
    end

end


function forward_model_PCfit(
        x,
        px :: Pixel_PCFit;       # Pixel struct
        return_components::Bool=false
    )

    # reflectance
    v     = collectPl.(px.λc, lmax=px.nPoly);
    ρ     = hcat(v...)' * x[1 : px.nPoly+1];

    # T↑ transmittance for SIF
    T₁    = (px.trans_mat * x[(px.nPoly+2):(px.nPoly+px.nPC+1)]);

    # T↓↑ transmittance for SIF
    smooth_x = 10. / (1 + exp( -x[px.nPoly+px.nPC+2]) ) + 1.;
    T₂       = @. exp( smooth_x * log(T₁) );

    # SIF magnitude
    SIF   = px.SIF_shape * x[px.nPoly+px.nPC+px.nSIF+2];

    # TOA radiance
    rad   = @. px.E * cosd(px.sza) / π * T₂ * ρ + SIF * T₁;
    
    if return_components
        return rad, ρ, T₁, T₂, SIF
    else
        return rad
    end
end

"""
    compute_transmittance(P_O2, VCD_O2, T_O2, P_H2O, VCD_H2O, T_H2O, MyKernel, o2_sitp, h2o_sitp, ν_grid)

	Calculate convolved transmittance from O2 and H2O atmospheric parameters.
	
	# Arguments
	- `P_O2`: O2 pressure [hPa]
	- `VCD_O2`: O2 vertical column density [molec/cm²]
	- `T_O2`: O2 temperature [K]
	- `P_H2O`: H2O pressure [hPa]
	- `VCD_H2O`: H2O vertical column density [molec/cm²]
	- `T_H2O`: H2O temperature [K]
	- `MyKernel`: Instrument convolution kernel matrix
	- `o2_sitp`: O2 cross-section interpolator (ν, P, T)
	- `h2o_sitp`: H2O cross-section interpolator (ν, P, T)
	- `ν_grid`: Wavenumber grid [cm⁻¹]
	
	# Returns
	- `T_conv`: Convolved transmittance spectrum in wavelength space
"""
function compute_transmittance(
        x,
        MyKernel,
        o2_sitp,
        h2o_sitp,
    )
    # Unpack atmospheric parameters
    P_O2     = x[1]        # O2 pressure [hPa]
    T_O2     = x[2]        # O2 temperature [K]
    VCD_O2   = 10^(x[3])   # O2 vertical column density [molec/cm²]
    P_H2O    = x[4]        # H2O pressure [hPa]
    T_H2O    = x[5]        # H2O temperature [K]
    VCD_H2O  = 10^(x[6])   # H2O vertical column density [molec/cm²]
    
    # detect whether p and T are out of bound
    P_O2   = clamp(P_O2,  0.01, 1040.0);
    P_H2O  = clamp(P_H2O, 0.01, 1040.0);
    T_O2   = clamp(T_O2,  150.0, 330.0);
    T_H2O  = clamp(T_H2O, 150.0, 330.0);

    # Get high-resolution cross-sections
    σ_O2   = o2_sitp(MyKernel.ν_grid, P_O2, T_O2)     # O2 cross-section [cm²/molec]
    σ_H2O  = h2o_sitp(MyKernel.ν_grid, P_H2O, T_H2O)  # H2O cross-section
    
    # Calculate optical depths
    τ_O2    = σ_O2 .* VCD_O2
    τ_H2O   = σ_H2O .* VCD_H2O
    τ_total = τ_O2 .+ τ_H2O
    
    # High-resolution transmittance
    T_highres       = exp.(-τ_total)
    T_highres_wvlen = reverse(T_highres)
    
    # Convolve to instrument resolution
    T_conv          = MyKernel.RSR_out * T_highres_wvlen
    
    return T_conv
end
