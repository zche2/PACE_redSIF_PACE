using ForwardDiff, DiffResults, LinearAlgebra, Statistics

"""
    center_wavelength(λ)

Normalize wavelength array to [-1, 1] range centered at midpoint.

Returns: λc (normalized wavelengths)
"""
function center_wavelength(
        λ::Vector{Union{Missing, FT}}
    ) where {FT <: AbstractFloat}
    
    λ_max = ceil(maximum(λ))
    λ_min = floor(minimum(λ))
    range = (λ_max - λ_min) / 2
    λ_middle = (λ_max + λ_min) / 2
    λc = (λ .- λ_middle) ./ range
    
    return λc
end


"""
root_mean_square(y_obs, y_retrieve)

Calculate root mean square error between observed and retrieved values.

Returns: RMS error
"""
function root_mean_square(y_obs, y_retrieve)
    
    n = length(y_obs)
    totSQ = sum((y_obs .- y_retrieve) .^ 2)
    RMS = sqrt(totSQ / n)
    
    return RMS
end


"""
    scale_transmittance(T, ind)

Normalize transmittance spectrum so maximum at baseline indices equals 1.

# Arguments
- `T`: Transmittance spectrum
- `ind`: Baseline wavelength indices

Returns: T_norm (normalized transmittance)
"""
function scale_transmittance(
        T,
        ind::Vector{Int64}
    ) where {FT <: AbstractFloat}
    
    T_abs = abs.(T)
    bl = maximum(T_abs[ind])
    T_norm = T_abs ./ bl
    
    return T_norm
end

function Jacobian(
        x, 
        model, 
        len ::Int  # length of measured spectrum
    )
	res = DiffResults.JacobianResult(zeros(len), x);
	ForwardDiff.jacobian!(res, model, x);
	K   = DiffResults.jacobian(res);
	val = DiffResults.value(res);
	return K, val
end

function GainMatrix(
        K,
        px :: Pixel  # pixel struct
    )
    inv_Sₑ = inv(px.Sₑ);
    inv_Sₐ = inv(px.Sₐ);
	return inv( K' * inv_Sₑ * K + inv_Sₐ )K' * inv_Sₑ;
end

function MakePriori!(
    px  :: Pixel,
    β   :: Vector{FT};
    γ   :: Vector{FT} = [(secd(px.sza) + secd(px.vza)) / secd(px.vza)],  
                # extra terms for T₁ and T₂ conversion
    SIF :: Vector{FT} = [px.nflh],  
                # coeff. for SIF PCs
    ) where {FT <: AbstractFloat}

    # polynomial terms
    K₀  = px.E .* cosd(px.sza) ./ pi .* hcat(collectPl.(px.λc, lmax=px.nPoly)...)';
    G₀  = inv( K₀'K₀ )K₀';
    x₀  = G₀ * px.R_toa; 

    px.xₐ = [x₀... β... γ... SIF...]';

    return nothing
end

function Retrieval_for_Pixel(
		# "L1B pixel-by-pixel vals"
		R_toa, sza, vza, nflh, chlor_a, flag,
		# params fixed for the retrieval scheme
		params,
	)

	if params isa RetrievalParams_xSecFit
		return Retrieval_for_Pixel_xSecFit(
			R_toa, sza, vza, nflh, chlor_a, flag,
			params,
		)
	elseif params isa RetrievalParams_PCFit
		return Retrieval_for_Pixel_PCFit(
			R_toa, sza, vza, nflh, chlor_a, flag,
			params,
		)
	else
		error("Unknown retrieval parameters type.")
	end

end

function Retrieval_for_Pixel_xSecFit(
		# "L1B pixel-by-pixel vals"
		R_toa, sza, vza, nflh, chlor_a, flag,
		# params fixed for the retrieval scheme
		params :: RetrievalParams_xSecFit,
	)
	# preprocess: if the flag is false, not doing the retrieval
	if ismissing(flag)
		return missing
	end

	MyPixel       = Pixel_xSecFit();
	MyModel       = params.forward_model
    MyIter        = params.iteration_method
	thr_Converge  = params.thr_Converge
	βₐ            = params.βₐ
	γₐ 		      = params.γₐ
    c₁            = params.c₁
    c₂            = params.c₂
	"for cross section fit, βₐ and γₐ represent one- and two-way transmittance priors"

	# Step1: construct struct
	MyPixel.λ  = params.λ;
	MyPixel.λc = params.λc;
	MyPixel.E       = params.E;
	MyPixel.nPoly   = params.nPoly;
	MyPixel.nLayer  = params.nLayer;
	MyPixel.nSIF    = params.nSIF;
	MyPixel.Sₐ      = params.Sₐ;
	MyPixel.SIF_shape = params.SIFComp[:, 1:MyPixel.nSIF];
	MyPixel.o2_sitp   = params.o2_sitp;
	MyPixel.h2o_sitp  = params.h2o_sitp;
	MyPixel.InstrumentKernel = params.InstrumentKernel;

	MyPixel.R_toa = R_toa;
	MyPixel.sza   = sza;
	MyPixel.vza   = vza;
	MyPixel.nflh  = nflh;
	MyPixel.chlor_a = chlor_a;
	noise         = sqrt.( c₁ .+ c₂ .* MyPixel.R_toa);
	MyPixel.Sₑ    = Diagonal(noise.^2);
	MyPixel.flag  = flag; 

	# set-up
	MakePriori!(MyPixel, βₐ, γ=γₐ);
	MyPixel.x  = copy(MyPixel.xₐ);
	MyPixel.y  = MyModel(MyPixel.x, MyPixel);
	MyPixel.iter_label = 0;
	MyPixel.ΔRMSE      = Inf;
	
	# Step2: iteration
	try
		MyIter(
			MyPixel, 
			MyModel
		)
	catch e
		println(e)
		return missing
	end

	# Step3: return
	# return if converge
	if abs(MyPixel.ΔRMSE) < thr_Converge
		# println("successfully retrieved")
		return MyPixel
	else
		return missing
	end

end

function Retrieval_for_Pixel_PCFit(
        # "L1B pixel-by-pixel vals"
		R_toa, sza, vza, nflh, chlor_a, flag,
        # params fixed for the retrieval scheme
		params :: RetrievalParams_PCFit,
    )

    # preprocess: if the flag is false, not doing the retrieval
	if ismissing(flag)
		return missing
	end
	
	MyPixel       = Pixel_PCFit();
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

	MyPixel.R_toa = R_toa;
	MyPixel.sza   = sza;
	MyPixel.vza   = vza;
	MyPixel.nflh  = nflh;
	MyPixel.chlor_a = chlor_a;
	noise         = sqrt.( c₁ .+ c₂ .* MyPixel.R_toa);
	MyPixel.Sₑ    = Diagonal(noise.^2);
	MyPixel.flag  = flag; 
	
	# set-up
	MakePriori!(MyPixel, βₐ);
	MyPixel.x  = copy(MyPixel.xₐ);
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
		println(e)
		return missing
	end

	# Step3: return
	# return if converge
	if abs(MyPixel.ΔRMSE) < thr_Converge
		# println("successfully retrieved")
		return MyPixel
	else
		return missing
	end

end