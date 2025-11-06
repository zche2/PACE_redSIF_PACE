### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ cc5c4ce8-8805-45ff-8203-00b18c49875b
begin
	import Pkg; 
	Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE");
end

# ╔═╡ ff37d1f2-08f6-4bd2-9611-cceada1c0e4b
using JLD2, Interpolations, Revise

# ╔═╡ 9bff5713-ab8d-4554-ae7b-b6f8793744d6
# ╠═╡ disabled = true
#=╠═╡
using Base.Threads
  ╠═╡ =#

# ╔═╡ 75a354c8-5b88-459a-9aa2-189229a8a532
using ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics

# ╔═╡ db2ae205-4d00-4135-8b76-a3ab99a5e4e9
using Polynomials, Random

# ╔═╡ 6c2b6b6e-175c-4678-92b5-87b033fa2fa8
using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools

# ╔═╡ 72e9abdf-838b-4f6f-8ec6-7d0253859b78
using PACE_SIF

# ╔═╡ 02bfa302-b66d-11f0-2f66-5d3686c10c23
md"""
## Create pseudo measurement
Created on 2025-10-31 🎃
"""

# ╔═╡ 2b4d2c63-ee1c-4b52-8aa5-eb23d70e2d18
#=╠═╡
println("Running with $(Threads.nthreads()) threads")
  ╠═╡ =#

# ╔═╡ 49880075-3490-4a9d-81a7-57109a7602f5
# wavelenth
λ_min = 610.; λ_max = 820.;

# ╔═╡ 52ff6a31-ecc2-4c47-ab08-ebfb2ead5f9b
md"""
##### Transmittance
"""

# ╔═╡ 31cc57c1-aed8-4de0-86bb-9c0987d0e382
begin
	# MERRA2 generated
	summer = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc");
	winter = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc");
	println("Opened datasets.")
	
	trans = cat(summer["transmittance"][:,:], winter["transmittance"][:,:], dims=1);
	println("\nConcatenated!")

	bands  = summer["band"][:];
	
	close(summer);
	close(winter);
	
	# PACE data
	oci = Dataset(
	"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/sample_granule_20240830T131442_new_chl.nc");
	red_band = oci["red_wavelength"][:];
	nflh     = oci["nflh"][:, :];
	vza      = oci["sensor_zenith"][:, :];
	sza      = oci["solar_zenith"][:, :];
	chlor_a  = oci["chlor_a"][:, :];
	println("\nRead in PACE Dataset")

	# select band (continuum spectrum)
	ind      = findall( λ_min .< red_band .< λ_max );
	E        = oci["red_solar_irradiance"][ind];
	R_toa    = oci["radiance_red"][:, :, ind];
	oci_band = red_band[ind];
	println("\nBand selected tp [$(λ_min), $(λ_max)]")
end

# ╔═╡ d0cd4fd3-a740-493d-a29d-99a0b9403f58
begin
    # transmittance spectra should also be aligned with SIF
    
    trans_new = zeros(size(trans, 1), length(oci_band))
    
    for i in 1:size(trans, 1)
        # Create interpolator for this row
        itp_row = LinearInterpolation(bands, trans[i, :], extrapolation_bc=0)
        
        # Evaluate at OCI bands
        trans_new[i, :] = itp_row.(oci_band)
    end
    
    println("Transmittance interpolated to OCI bands")
end

# ╔═╡ ef0922d4-2572-4327-9499-cb828e27e5e4
begin
	plot(
		oci_band,
		trans_new[1:400:end,:]',
		label="",
		color=:silver,
		size=(600, 250),
		title="number of transmittance spectra: $(size(trans_new,1))",
		titlefontsize=8,
	)
	plot!(oci_band, mean(trans_new', dims=2), lw=2.5, label="")
end

# ╔═╡ 83f4ae07-fbfb-46e0-9806-1a89ff99527a
md"""
##### SIF
"""

# ╔═╡ ca262d91-7318-493b-ad0d-44e34ce3961c
begin
	# load SIF
	SIF_shape_dict = JLD2.load("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/SIF_singular_vector.jld2")
	
	# itp.
	# interpolation in the first dimension and no interp. in the second
	itp    = interpolate(SIF_shape_dict["SIF_shapes"], (BSpline(Linear()), NoInterp()));
	
	# scale
	range₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end];
	range₂ = 1:size(itp, 2);
	sitp   = scale(itp, range₁, range₂);

	# set extrapolation filling value = 0
	setp0  = extrapolate(sitp, 0)

	# interpolation
	SIF_new = reduce(hcat, [setp0.(oci_band, i) for i in range₂]); 

	# scaled
	scale_factor_SIF = 20;
	SIF_new *= scale_factor_SIF;
	
	println("SIF shape interpolated & scaled by a factor of $scale_factor_SIF")
	# vis. 
	plot(
		 oci_band,
		 SIF_new[:, 1:5:end], 
		 label="",
		 size=(600, 250),
		 title="wl SIF radiance, scaled already\n number of shapes: $(size(itp, 2))",
		 titlefontsize=8
	)
end

# ╔═╡ 2a0ac581-5e43-44ef-a8d4-a30ccf4d1c0f
md"""
##### Reflctance
Get some approximation of it by fitting some non-absorption windows
"""

# ╔═╡ 781e7324-2768-4f17-a0cc-96f0bf5d078e
begin
	# wavelength
	λ = oci_band;
	
	# get center wavelength
	λc = center_wavelength(oci_band);

	# find transmittance baseline
	λ_bl_ref = [607.99, 610.36, 612.73, 615.14, 617.6, 620.06, 622.53, 669.52, 670.76, 671.99, 673.24, 674.51, 675.73, 676.96, 678.21, 679.45, 754.3, 779.33, 867.11, 869.61, 872.13];
	bl_ind = map(λ_bl_ref -> argmin(abs.(λ .- λ_bl_ref)), λ_bl_ref);
	
	# Extract rad irradiance where nFLH is detectable but negligible
	valid_mask  = findall(coalesce.(nflh .< 0.005, false))
	R_noSIF     = R_toa[valid_mask,:];
	sza_noSIF   = sza[valid_mask];
	vza_noSIF   = vza[valid_mask];
	R_baseline  = R_noSIF[:,bl_ind];
	E_baseline  = E[bl_ind];
	λ_baseline  = λ[bl_ind];
	λc_baseline = λc[bl_ind];

	# fit
	order = 6;
	n_pixels = size(R_baseline, 1);
	K₀       = hcat(collectPl.(λc_baseline, lmax=order)...)';
	K₀_recon = hcat(collectPl.(λc, lmax=order)...)';
	
	# Preallocate for fitted values and coefficients
	R_fitted = zeros(n_pixels, length(λc))
	coeffs_record = zeros(n_pixels, order + 1)  # Store coefficients for each pixel
	
	# Fit and record coefficients
	for i in 1:n_pixels
		# Fit
		coeff₀ = inv( K₀'K₀ )K₀' * ( R_baseline[i, :] .* pi ./ ( E_baseline .* cosd(sza_noSIF[i]) ) );

	    # Evaluate at all wavelengths
	    R_fitted[i, :] = K₀_recon * coeff₀ .* ( E .* cosd(sza_noSIF[i]) ./ pi);
	    
	    # Store coefficients
	    coeffs_record[i, :] = coeff₀;
	end
	println("Reflectance obtained!")
	
	# vis.
	plot(
		λ,
		R_noSIF[1:250:end,:]',
		label="",
		size=(600, 250),
		title="Background pixel (nFLH<.005) TOA radiance and Polynomial=$order fitted\n number of px: $n_pixels",
		titlefontsize=8
	)
	plot!(
		λ,
		R_fitted[1:250:end,:]',
		ls=:dash,
		label="",
	)
end

# ╔═╡ 7ce2ad51-4d30-4008-b4b2-94e90330df65
md"""
##### Pseudo measurement - add noise
"""

# ╔═╡ fceee754-f52a-4725-b860-95a4c95129bf
begin
	filename = raw"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt";
	lines = readlines(filename);
	end_header_index = findfirst(x -> x == "/end_header", lines);
	data  = readdlm(filename, String, skipstart=end_header_index);

	FPA   = data[:, 1];                   # 1st column: band
	wvlen = parse.(Float64, data[:, 2]);  # 2nd column: center wavelength

	wv_val  = (λ_min .< wvlen .< λ_max);
	snr_ind = findall((FPA .== "Red") .& wv_val);

	# get c1 and c2 at that range
	c1    = parse.(Float64, data[snr_ind, 4]);  # 4th column: c1
	c2    = parse.(Float64, data[snr_ind, 5]);  # 5th column: c2
	println("SNR retrieved")
end

# ╔═╡ 91f48071-853a-4701-a422-896a350482d0
begin
	nᵨ = n_pixels;   # also select SZA from
	nₛ = size(SIF_new, 2);
	nₜ = size(trans_new, 1);

	# generate this amount of samples
	n_sample = 5000;
	Random.seed!(52)
	indᵨ     = rand(1:nᵨ, n_sample);
	indₛ     = rand(1:nₛ, n_sample);
	indₜ₁    = rand(1:nₜ, n_sample);
	indₜ₂    = rand(1:nₜ, n_sample);
	ind_sza  = rand(1:nᵨ, n_sample);   # sza
	ind_vza  = rand(1:nᵨ, n_sample);   # vza
end

# ╔═╡ b42dc5f2-a93f-4254-8601-ea829697d180
begin
	# Preallocate storage for each component
	len_λ   = length(λ);
	ρ_all   = zeros(n_sample, len_λ);
	μ₁_all  = zeros(n_sample);
	μ₂_all  = zeros(n_sample);
	T₁_all  = zeros(n_sample, len_λ);
	T₂_all  = zeros(n_sample, len_λ);
	SIF_all = zeros(n_sample, len_λ);
	pseudo_obs_all = zeros(n_sample, len_λ);
	
	# Loop over samples and store each component
	for i in 1:n_sample
	    # ----- rho -----
	    ρ_all[i, :] = K₀_recon * coeffs_record[indᵨ[i], :];
	    
	    # ----- cos(sza) and cos(vza) -----
	    μ₁_all[i] = cosd(sza_noSIF[ind_sza[i]]);
		μ₂_all[i] = cosd(vza_noSIF[ind_vza[i]]);
	    
	    # ----- Transmittance -----
		σ₁ = @. - 1 / μ₁_all[i] * log( trans_new[indₜ₁[i], :] );
		σ₂ = @. - 1 / μ₂_all[i] * log( trans_new[indₜ₂[i], :] );
	    T₁_all[i, :] = @. exp( - σ₁ );
	    T₂_all[i, :] = @. exp( - σ₁ - σ₂ );
	    
	    # ----- water-leaving SIF -----
	    SIF_all[i, :] = SIF_new[:, indₛ[i]];
	    
	    # ----- TOA -----
	    pseudo_obs_all[i, :] = @. E / pi * μ₁_all[i] * ρ_all[i, :] * T₂_all[i, :] + SIF_all[i, :] * T₁_all[i, :];

		# ----- noise -----
		stds = sqrt.(c1 .+ c2 .* pseudo_obs_all[i, :]);
		pseudo_obs_all[i, :] += randn(len_λ) .* stds;
	    
	    if i % 1000 == 0
	        println("Processed $i / $n_sample samples")
	    end
	end
end

# ╔═╡ b1cc4469-ce38-4dbb-908e-7cc1f0f6a630
begin
	Δn = 400;
	
	# Create 5x1 subplot layout
	p = plot(layout=(5, 1), size=(800, 1000), legend=false)
	
	# Plot each variable
	plot!(p[1], λ, ρ_all[1:Δn:end,:]', title="Reflectance", ylabel="ρ", xlabel="")
	plot!(p[2], λ, T₁_all[1:Δn:end,:]', title="Transmittance (one-way)", ylabel="T", xlabel="")
	plot!(p[3], λ, T₂_all[1:Δn:end,:]', title="Transmittance (two-way)", ylabel="T", xlabel="")
	plot!(p[4], λ, SIF_all[1:Δn:end,:]', title="SIF", ylabel="SIF", xlabel="")
	plot!(p[5], λ, pseudo_obs_all[1:Δn:end,:]', title="TOA", ylabel="Radiance", xlabel="Wavelength (nm)")

end

# ╔═╡ c90eda45-e488-4320-9e17-296da1197e93
md"""
## Try one-pixel retrieval
Go bruins! 🦫
"""

# ╔═╡ 0a000d90-d661-4ebe-bc08-2c8b2a3abc0f
begin
	rank       = 15;
	# NMF
	HighResNMF = Spectral_NMF(
		trans, 
		bands,
		Float64.(collect(skipmissing(oci_band))); 
		rank=rank
	);
	# W and H
	λ₀ = HighResNMF.band;
	W₀ = HighResNMF.Loading;
	H₀ = HighResNMF.PrinComp;
	
	# matrics
    mean_val  = [round(mean(W₀[:, i]), digits=2) for i in 1:rank];
    max_val   = [round(maximum(W₀[:, i]), digits=2) for i in 1:rank];
    min_val   = [round(minimum(W₀[:, i]), digits=2) for i in 1:rank];

	# s.d. for the loading term
	loading_ave = [mean(W₀[:, i]) for i in 1:rank];
	@show loading_sd  = [var(W₀[:, i]) for i in 1:rank];
end

# ╔═╡ 09d8afad-d4f4-4625-9940-f1c82ac957aa
begin
	# set up retrieval scheme
	n     = 10;
	nPC   = rank;
	nSIF  = 1;

	Sₐ   = I(n+nPC+nSIF+2) .+ 0.;
	# update diagonal term
	for i=1:(n+1)
	    Sₐ[i,i] = 1e10;
		# large variance applies no constrain to these polynomial term
	end
	
	# \beta
	for i=(n+2):(n+nPC+1)
		Sₐ[i,i]  = loading_sd[i - (n+1)];
	end
	
	# \gamma
	Sₐ[n+nPC+2, n+nPC+2] = 2;
	# SIF magnitude
	Sₐ[end, end] = 1;
	println("Diagonal terms are: $(diag(Sₐ))")
	
end

# ╔═╡ a3dfbedc-7c7a-4559-bbf3-0252857465b4
begin
	# SIF U
	# interpolation in the first dimension and no interp. in the second
	itp₂    = interpolate(SIF_shape_dict["SIF_U"], (BSpline(Linear()), NoInterp()));
	
	# scale
	r₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end];
	r₂ = 1:size(itp₂, 2);
	sitp₂   = scale(itp₂, r₁, r₂);

	# set extrapolation filling value = 0
	setp0₂  = extrapolate(sitp₂, 0)

	# interpolation
	SIF_PC  = reduce(hcat, [setp0₂.(λ, i) for i in range₂]); 
	
	println("SIF shape interpolated")
end

# ╔═╡ 339dea16-d9b2-4f9f-b784-ce0e3849ef3b
# Create the retrieval parameters
params = RetrievalParams(
    # Measurement specific
    λ  = oci_band,                   # Wavelength array
    λc = λc, 						 # Centered wavelength
    λ_bl_ind = bl_ind,               # Baseline band indices
    E        = E,                    # Solar irradiance
	c₁       = c1, 					 # PACE SNR 
	c₂       = c2, 			       	 # PACE SNR
    
    # Forward model settings
    forward_model = forward_model,
    nPoly = n,                       # Degree of Legendre polynomial
    nPC   = nPC,                     # Number of transmittance PCs
    nSIF  = nSIF,                    # Number of SIF PCs
    Sₐ = Sₐ,   					     # Prior covariance
    βₐ = loading_ave,                # Prior state
    PrinComp = HighResNMF.PrinComp', # Principal components
    SIFComp  = SIF_PC,       # SIF components
    
    # Iteration settings (optional, have defaults)
    iteration_method = LM_Iteration!,
    nIter = 25,
    thr_Converge = 1e-6
)

# ╔═╡ 838bd7b3-7a74-411c-b7d9-f003353d4a1f
begin
	k = 4000;
	# Single pixel ref.
	MyPixel = Retrieval_for_Pixel(
		pseudo_obs_all[k,:],
		sza_noSIF[ind_sza[k]],
		vza_noSIF[ind_vza[k]],
		maximum(SIF_new[:, indₛ[k]]),
		1.0,
		1.0,
		params
	)
end

# ╔═╡ 164ff48c-89a0-45c0-8d3c-af86262a5267
md"""
##### Reconstruct
"""

# ╔═╡ ad35e7ab-5158-4d3d-aa80-a54301b1eedd
begin
	_, ρ, T₁, T₂, SIF = forward_model(MyPixel.x, MyPixel, return_components=true)
	
    # Plot TOA radiance, fitted baseline, and residuals
    p1 = plot(layout=(6,1), size=(800, 1200), legend=true)
    
    # Sample indices
    sample_idx = k;
    
    # 1. ρ (Reflectance)
    plot!(p1[1], λ, ρ_all[sample_idx, :], title="ρ (degree=$n)", ylabel="ρ", lw=1., label="truth")
	plot!(p1[1], λ, ρ, lw=1.5, label="fit")
    
    # 2. T₁ (One-way Transmittance)
    plot!(p1[2], λ, T₁_all[sample_idx, :], title="T₁ (nPC=$nPC)", ylabel="T₁", lw=1., label="truth")
	plot!(p1[2], λ, T₁, lw=1.5, label="fit")
    
    # 3. T₂ (Two-way Transmittance)
    plot!(p1[3], λ, T₂_all[sample_idx, :], title="T₂", ylabel="T₂", lw=1., label="truth")
	plot!(p1[3], λ, T₂, lw=1.5, label="fit")
	
    # 4. SIF
    plot!(p1[4], λ, SIF_all[sample_idx, :], title="SIF", ylabel="SIF", lw=1., label="truth")
	plot!(p1[4], λ, SIF, lw=1.5, label="fit")
	
    # 5. Residual (Observed - Fitted)
    residual = @. MyPixel.y - MyPixel.R_toa;
    plot!(p1[5], λ, residual, title="Residual", label="Fit - Obs",
          ylabel="Residual", lw=1.5)
	
	# 6. SIF
    plot!(p1[6], λ, MyPixel.R_toa, title="TOA radiance", ylabel="radiance", lw=1.5, label="truth")
	plot!(p1[6], λ, MyPixel.y, xlabel="Wavelength [nm]", lw=1., label="fit")

	# title
    plot!(p1, titlefontsize=9)
    p1
end

# ╔═╡ 6df90f60-e318-4196-aefc-822ac5e5d551
md"""
## For all pseudo measurements
Excited to unveil 🎶

"""

# ╔═╡ 9f18e909-5a77-44ac-a510-ac4191282bcd
#=╠═╡
println("Running with $(Threads.nthreads()) threads")
  ╠═╡ =#

# ╔═╡ a82609f4-b8ab-4601-ab91-df6380764195
# ╠═╡ disabled = true
#=╠═╡
@threads for i in 1:n_sample
    try
        Retrieval_all[i] = Retrieval_for_Pixel(
            pseudo_obs_all[i, :],
            sza_pxs[i],
            vza_pxs[i],
            SIF_pxs[i],
            1.0,
            1.0,
            params
        )
    catch e
        @warn "Pixel $i failed" exception=e
        Retrieval_all[i] = missing
    end
end
  ╠═╡ =#

# ╔═╡ 064feee6-92bd-4fcd-852b-0a356b99fdc4
# ╠═╡ disabled = true
#=╠═╡
begin
	# arrs.
	sza_pxs = sza_noSIF[ind_sza];
	vza_pxs = vza_noSIF[ind_vza];
	SIF_pxs = maximum(SIF_new[:, indₛ], dims=1);
	# Multi pixel
	Retrieval_all = Retrieval_for_Pixel.(
						eachslice(pseudo_obs_all, dims=1),
						sza_pxs,
						vza_pxs,
						SIF_pxs,
						1.0,
						1.0,
						Ref(params)
					)
end
  ╠═╡ =#

# ╔═╡ a4588031-e1f2-4d40-9735-41ca7a78c7f7
#=╠═╡
Retrieval_all = Vector{Union{Pixel, Missing}}(undef, n_sample)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─02bfa302-b66d-11f0-2f66-5d3686c10c23
# ╠═9bff5713-ab8d-4554-ae7b-b6f8793744d6
# ╠═2b4d2c63-ee1c-4b52-8aa5-eb23d70e2d18
# ╠═cc5c4ce8-8805-45ff-8203-00b18c49875b
# ╠═ff37d1f2-08f6-4bd2-9611-cceada1c0e4b
# ╠═75a354c8-5b88-459a-9aa2-189229a8a532
# ╠═db2ae205-4d00-4135-8b76-a3ab99a5e4e9
# ╠═6c2b6b6e-175c-4678-92b5-87b033fa2fa8
# ╠═72e9abdf-838b-4f6f-8ec6-7d0253859b78
# ╠═49880075-3490-4a9d-81a7-57109a7602f5
# ╟─52ff6a31-ecc2-4c47-ab08-ebfb2ead5f9b
# ╠═31cc57c1-aed8-4de0-86bb-9c0987d0e382
# ╟─d0cd4fd3-a740-493d-a29d-99a0b9403f58
# ╟─ef0922d4-2572-4327-9499-cb828e27e5e4
# ╟─83f4ae07-fbfb-46e0-9806-1a89ff99527a
# ╟─ca262d91-7318-493b-ad0d-44e34ce3961c
# ╟─2a0ac581-5e43-44ef-a8d4-a30ccf4d1c0f
# ╠═781e7324-2768-4f17-a0cc-96f0bf5d078e
# ╟─7ce2ad51-4d30-4008-b4b2-94e90330df65
# ╠═fceee754-f52a-4725-b860-95a4c95129bf
# ╠═91f48071-853a-4701-a422-896a350482d0
# ╠═b42dc5f2-a93f-4254-8601-ea829697d180
# ╠═b1cc4469-ce38-4dbb-908e-7cc1f0f6a630
# ╟─c90eda45-e488-4320-9e17-296da1197e93
# ╠═0a000d90-d661-4ebe-bc08-2c8b2a3abc0f
# ╠═09d8afad-d4f4-4625-9940-f1c82ac957aa
# ╠═a3dfbedc-7c7a-4559-bbf3-0252857465b4
# ╠═339dea16-d9b2-4f9f-b784-ce0e3849ef3b
# ╠═838bd7b3-7a74-411c-b7d9-f003353d4a1f
# ╟─164ff48c-89a0-45c0-8d3c-af86262a5267
# ╟─ad35e7ab-5158-4d3d-aa80-a54301b1eedd
# ╟─6df90f60-e318-4196-aefc-822ac5e5d551
# ╠═a4588031-e1f2-4d40-9735-41ca7a78c7f7
# ╠═9f18e909-5a77-44ac-a510-ac4191282bcd
# ╠═a82609f4-b8ab-4601-ab91-df6380764195
# ╠═064feee6-92bd-4fcd-852b-0a356b99fdc4
