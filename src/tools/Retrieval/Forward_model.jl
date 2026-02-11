using LegendrePolynomials


"""
    forward_model_ToBeUpdated(x, px::Pixel; return_components::Bool=false)
    This is incomplete! 🟢 Please update the function name and docstring accordingly.
"""
function forward_model_ToBeUpdated(
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

end;


function forward_model(
        x,
        px :: Pixel_xSecFit,
        params :: RetrievalParams_xSecFit;
        return_components::Bool=false
    )
    # unpack params
    o2_sitp = params.o2_sitp;
    h2o_sitp = params.h2o_sitp;
    InstrumentKernel = params.InstrumentKernel;

    # reflectance
    t_rho = @elapsed begin
        xᵨ    = x[1 : px.nPoly+1]
        v     = collectPl.(px.λc, lmax=px.nPoly);
        ρ     = hcat(v...)' * xᵨ;
    end
    println("  [FwdModel] Reflectance: $(round(t_rho*1e6, digits=1)) μs")

    # T↑ transmittance for SIF
    t_t1 = @elapsed begin
        x₁    = x[(px.nPoly+2):(px.nPoly+7)];
        T₁    = compute_transmittance(
                    x₁,
                    InstrumentKernel,
                    o2_sitp,
                    h2o_sitp,
                );
    end
    println("  [FwdModel] T₁ transmittance: $(round(t_t1*1e3, digits=2)) ms")

    # T↓↑ transmittance for reflected radiance
    t_t2 = @elapsed begin
        x₂    = x[(px.nPoly+8):(px.nPoly+13)];
        T₂    = compute_transmittance(
                    x₂,
                    InstrumentKernel,
                    o2_sitp,
                    h2o_sitp,
                );
    end
    println("  [FwdModel] T₂ transmittance: $(round(t_t2*1e3, digits=2)) ms")

    # SIF magnitude
    t_sif = @elapsed begin
        xₛ   = x[(end - px.nSIF + 1) : end];
        SIF  = px.SIF_shape * xₛ;
    end
    println("  [FwdModel] SIF: $(round(t_sif*1e6, digits=1)) μs")

    # TOA radiance
    t_toa = @elapsed begin
        rad  = @. px.E * cosd(px.sza) / π * T₂ * ρ + SIF * T₁;
    end
    println("  [FwdModel] TOA radiance: $(round(t_toa*1e6, digits=1)) μs")
    println("  [FwdModel] Total: $(round((t_rho+t_t1+t_t2+t_sif+t_toa)*1e3, digits=2)) ms")
    
    if return_components
        return rad, ρ, T₁, T₂, SIF
    else
        return rad
    end

end

"""
    forward_model_xSecFit(x, px::Pixel; return_components::Bool=false)

Notes:
- fit one- and two-way transmittance respectively, assuming completely independent.
"""

function forward_model_xSecFit_old(
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
    @. T₂    = sexp( smooth_x * log(T₁) );

    # SIF magnitude
    SIF    = px.SIF_shape * x[px.nPoly+px.nPC+px.nSIF+2];

    # TOA radiance
    @. rad = px.E * cosd(px.sza) / π * T₂ * ρ + SIF * T₁;
    
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
    t_unpack = @elapsed begin
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
    end
    println("    [CompTrans] Unpack/clamp: $(round(t_unpack*1e6, digits=1)) μs")

    # Get high-resolution cross-sections
    t_xsec = @elapsed begin
        σ_O2   = o2_sitp(MyKernel.ν_grid, P_O2, T_O2)     # O2 cross-section [cm²/molec]
        σ_H2O  = h2o_sitp(MyKernel.ν_grid, P_H2O, T_H2O)  # H2O cross-section
    end
    println("    [CompTrans] Cross-section interp: $(round(t_xsec*1e3, digits=2)) ms")
    
    # Calculate optical depths
    t_tau = @elapsed begin
        τ_O2    = σ_O2 .* VCD_O2
        τ_H2O   = σ_H2O .* VCD_H2O
        τ_total = τ_O2 .+ τ_H2O
    end
    println("    [CompTrans] Optical depth: $(round(t_tau*1e6, digits=1)) μs")
    
    # High-resolution transmittance
    t_trans = @elapsed begin
        T_highres       = exp.(-τ_total)

        # Align to kernel wavelength ordering.
        # Kernel spectral axis can be native wavenumber or wavelength, depending on LUT build.
        spec_grid = collect(MyKernel.ν_grid)
        λ_from_spec = maximum(spec_grid) > 3000.0 ? ν_to_λ.(spec_grid) : spec_grid
        wvlen_out = collect(MyKernel.wvlen_out)

        err_direct = abs(λ_from_spec[1] - wvlen_out[1]) + abs(λ_from_spec[end] - wvlen_out[end])
        err_reverse = abs(λ_from_spec[end] - wvlen_out[1]) + abs(λ_from_spec[1] - wvlen_out[end])
        T_highres_wvlen = err_direct <= err_reverse ? T_highres : reverse(T_highres)
    end
    println("    [CompTrans] Transmittance calc: $(round(t_trans*1e6, digits=1)) μs")
    
    # Convolve to instrument resolution
    t_conv = @elapsed begin
        T_conv          = MyKernel.RSR_out * T_highres_wvlen
    end
    println("    [CompTrans] Convolution: $(round(t_conv*1e6, digits=1)) μs")
    println("    [CompTrans] Total: $(round((t_unpack+t_xsec+t_tau+t_trans+t_conv)*1e3, digits=2)) ms")
    
    return T_conv
end
