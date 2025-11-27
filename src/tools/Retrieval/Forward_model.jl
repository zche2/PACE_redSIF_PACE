using LegendrePolynomials

function forward_model(
        x,
        px :: Pixel;           # Pixel struct
        if_log ::Bool=false,   # log transformation if doing log-SVD to transmittance
        return_components::Bool=false
    )

    # reflectance
    v     = collectPl.(px.λc, lmax=px.nPoly);
    ρ     = hcat(v...)' * x[1 : px.nPoly+1];

    # T↑ and T↓↑ transmittance for SIF
    smooth_x  = 10. / (1 + exp( -x[px.nPoly+px.nPC+2]) ) + 1.;
                    # tanh(-x[px.nPoly+px.nPC+2]) + 2.;

    if if_log
        logT₁     = (px.trans_mat * x[(px.nPoly+2):(px.nPoly+px.nPC+1)]);
        T₁        = @. exp(logT₁);
        T₂        = @. exp( smooth_x * logT₁ );
    else
        T₁       = (px.trans_mat * x[(px.nPoly+2):(px.nPoly+px.nPC+1)]);
        T₂       = @. exp( smooth_x * log(T₁) );
    end

    # SIF magnitude
    SIF   = px.SIF_shape * x[(px.nPoly+px.nPC+3):(px.nPoly+px.nPC+px.nSIF+2)];

    # TOA radiance
    rad   = @. px.E * cosd(px.sza) / π * T₂ * ρ + SIF * T₁;
    
    if return_components
        return rad, ρ, T₁, T₂, SIF
    else
        return rad
    end
end
