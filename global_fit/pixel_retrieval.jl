# =========================================================================================================
# This script does single pixel retrieval with designated parameters. 
# It returns
# (1) econstructed SIF; (2) state vector; (3) residual; (4) convergence status.
# ！for now， we only return the state vector
# all the other components can be recovered from the state vector and the forward model.    
# =========================================================================================================

using PACE_SIF

function process_all_pixels(
    R_toa,      # 3D array: (n_bands, n_pixels, n_scans)
    sza,        # 2D array: (n_pixels, n_scans)
    vza,        # 2D array: (n_pixels, n_scans)
    nflh,       # 2D array: (n_pixels, n_scans)
    chlor_a,    # 2D array: (n_pixels, n_scans)
    flag,       # 2D array: (n_pixels, n_scans)
    params::RetrievalParams
    )

    n_pixels, n_scans, _ = size(R_toa)
    total_pixels = n_pixels * n_scans
    println("Total pixels to process: $total_pixels")

    # elements in state vector
    # recall that state vector = [x₀, β, γ, SIF_terms]
    n_state = 1 + params.nPC + params.nPoly + 1 + params.nSIF
    
    # Preallocate results
    results = Array{Union{Missing, Float64}}(undef, n_pixels, n_scans, n_state)
    
    # a shared flag to stop all threads once error are detected
    stop_flag = Atomic{Bool}(false)
    
    @threads for idx in 1:total_pixels
        if stop_flag[]
            break
        end

        # Convert linear index to 2D indices
        j = ((idx - 1) % n_pixels) + 1      # pixel index
        i = ((idx - 1) ÷ n_pixels) + 1      # scan index

        try
            results[j, i, :] = retrieve_pixel(
                R_toa[:, j, i],    # Extract spectrum for this pixel (all bands)
                sza[j, i],
                vza[j, i],
                nflh[j, i],
                chlor_a[j, i],
                flag[j, i],
                params
            )
        catch e
            @warn "Failed at pixel ($j, $i): $e"
        end
        
        # Progress reporting
        if idx % 50 == 0
            atomic_cas!(stop_flag, false, true)
            println("Processed $idx / $total_pixels pixels")
        end
    end

    return results
end

function retrieve_pixel(
        # "L1B pixel-by-pixel vals"
        R_toa, sza, vza, nflh, chlor_a, flag,
        # params fixed for the retrieval scheme
        params :: RetrievalParams,
    )::Union{Missing, Vector{Float64}}

    # preprocess: if the flag is false, not doing the retrieval
    if !flag
        return missing
    end

    MyPixel       = Pixel();
    MyModel       = params.forward_model
    MyIter        = params.iteration_method
    nIter         = params.nIter
    thr_Converge  = params.thr_Converge
    βₐ            = params.βₐ
    c₁            = params.c₁
    c₂            = params.c₂

    # Step1: construct struct
    MyPixel.λ  = params.λ;
    MyPixel.λc = params.λc;
    MyPixel.λ_bl_ind = params.λ_bl_ind;
    MyPixel.E     = params.E;
    MyPixel.nPoly = params.nPoly;
    MyPixel.nPC   = params.nPC;
    MyPixel.nSIF  = params.nSIF;
    MyPixel.Sₐ    = params.Sₐ;
    MyPixel.trans_mat = params.PrinComp[:, 1:MyPixel.nPC];
    MyPixel.SIF_shape = params.SIFComp[:, 1:MyPixel.nSIF];
    MyPixel.Ŝ     = Matrix{Float64}(undef, size(params.Sₐ));

    MyPixel.R_toa = R_toa;
    MyPixel.sza   = sza;
    MyPixel.vza   = vza;
    MyPixel.nflh  = nflh;
    MyPixel.chlor_a = chlor_a;
    noise         = sqrt.( c₁ .+ c₂ .* MyPixel.R_toa);
    MyPixel.Sₑ    = Diagonal(noise.^2);
    MyPixel.flag  = flag; 

    # set-up
    MakePriori!(MyPixel, βₐ, MyPixel.nSIF);
    MyPixel.x  = MyPixel.xₐ;
    MyPixel.y  = MyModel(MyPixel.x, MyPixel);
    MyPixel.iter_label = 0;

    # Step2: iteration
    try
        MyIter(
            MyPixel, 
            MyModel,
            nIter=nIter,
            thr_Converge=thr_Converge
        )
    catch e
        return missing
    end

    # Step3: return
    # return if converge
    if abs(MyPixel.ΔRMSE) < thr_Converge
        return MyPixel.x[:]
    else
        return missing
    end
end