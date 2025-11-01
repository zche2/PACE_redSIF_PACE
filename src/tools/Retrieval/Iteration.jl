"""
    GN_Iteration!(px, model; nIter, thr_Converge)

Gauss-Newton iterative retrieval for retrieval.
"""
function GN_Iteration!(
        px :: Pixel,
        model;
        nIter = 20,
        thr_Converge = 1e-8,
    )

    # initial
    xₐ  = px.xₐ;   # priori estimation
    xₙ  = px.x;
    len = length(px.R_toa);
    Kₙ, px.y = Jacobian(xₙ, x -> model(x, px), len);
    RMSE₀    = 1e20; 
    RMSE₁    = root_mean_square(px.R_toa, px.y);
    ΔRMSE    = RMSE₁ - RMSE₀

    # loop
    while ( abs(ΔRMSE) > thr_Converge ) & ( px.iter_label < nIter )
        # k += 1
        # get Gain matrix
        Gₙ     = GainMatrix(Kₙ, px);
        # retrieval
        xₙ₊₁   = xₐ .+ Gₙ * (px.R_toa .- px.y .+ Kₙ * ( px.x .- xₐ ) );
        # update x and y
        px.x   = xₙ₊₁;
        Kₙ₊₁, yₙ₊₁ = Jacobian(xₙ₊₁, x -> model(x, px), len);
        px.y   = yₙ₊₁;
        Kₙ     = Kₙ₊₁;
        # iter ++
        px.iter_label += 1;
        # test convergence
        RMSE₀  = RMSE₁;
        RMSE₁  = root_mean_square(px.R_toa, px.y);
        ΔRMSE  = RMSE₁ - RMSE₀;
        px.ΔRMSE = ΔRMSE;
    end

    return nothing
end

"""
    LM_Iteration!(px::Pixel, model; nIter::Int=20, thr_Converge::Float64=1e-8,
                  γ_init::Float64=10.0, γ⁺::Float64=2.0, γ⁻::Float64=1.5)
Levenberg-Marquardt optimization for retrieval.
"""
function LM_Iteration!(
        px::Pixel,
        model;
        nIter::Int = 20,
        thr_Converge::Float64 = 1e-8,
        γ_init::Float64 = 10.0,
        γ⁺::Float64 = 10.0,          # Factor to increase γ when step fails
        γ⁻::Float64 = 2.0            # Factor to reduce γ when step accepted.
    )
    
    # Initialize
    Sa_inv = inv(px.Sₐ)
    Se_inv = inv(px.Sₑ)
    len    = length(px.R_toa);
    γ      = γ_init;
    
    Kₙ, px.y = Jacobian(px.x, x -> model(x, px), len)
    RMSE₀ = 1e20;
    RMSE₁ = root_mean_square(px.R_toa, px.y);
    ΔRMSE = RMSE₁ - RMSE₀;
    
    # Iteration loop
    while ( abs(ΔRMSE) > thr_Converge ) && ( px.iter_label < nIter )
        
        px.iter_label += 1;

        # Update x
        Δy = Kₙ' * Se_inv * (px.R_toa .- px.y) - Sa_inv * (px.x .- xₐ)
        Gₙ = (1 + γ) * Sa_inv + Kₙ' * Se_inv * Kₙ
        Δx = inv(Gₙ) * Δy
        
        x_trial = px.x .+ Δx
        
        # Evaluate at trial point
        K_trial, y_trial = Jacobian(x_trial, x -> model(x, px), len)
        RMSE_trial       = root_mean_square(px.R_toa, y_trial)
        
        # Check if step improves the fit
        if RMSE_trial < RMSE₁
            # Accept step
            px.x = x_trial
            px.y = y_trial
            Kₙ   = K_trial
            # Update RMSE
            RMSE₀ = RMSE₁
            RMSE₁ = RMSE_trial
            ΔRMSE = RMSE₁ - RMSE₀
            px.ΔRMSE = ΔRMSE
            # Decrease damping (success)
            γ /= γ⁻
            
        else
            # Reject, increase damping
            px.iter_label -= 1
            γ *= γ⁺
        end
    end
    
    return px
end
