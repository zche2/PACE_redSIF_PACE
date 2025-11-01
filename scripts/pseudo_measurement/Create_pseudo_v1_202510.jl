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
	indₐ     = rand(1:nᵨ, n_sample);   # sza
end

# ╔═╡ b42dc5f2-a93f-4254-8601-ea829697d180
begin
	# Preallocate storage for each component
	len_λ   = length(λ);
	ρ_all   = zeros(n_sample, len_λ);
	μ_all   = zeros(n_sample);
	T₁_all  = zeros(n_sample, len_λ);
	T₂_all  = zeros(n_sample, len_λ);
	SIF_all = zeros(n_sample, len_λ);
	pseudo_obs_all = zeros(n_sample, len_λ);
	
	# Loop over samples and store each component
	for i in 1:n_sample
	    # ----- rho -----
	    ρ_all[i, :] = K₀_recon * coeffs_record[indᵨ[i], :];
	    
	    # ----- cos(sza) -----
	    μ_all[i] = cosd(sza_noSIF[indₐ[i]]);
	    
	    # ----- Transmittance -----
	    T₁_all[i, :] = trans_new[indₜ₁[i], :];
	    T₂_all[i, :] = @. T₁_all[i, :] * trans_new[indₜ₂[i], :];
	    
	    # ----- water-leaving SIF -----
	    SIF_all[i, :] = SIF_new[:, indₛ[i]];
	    
	    # ----- TOA -----
	    pseudo_obs_all[i, :] = @. E / pi * μ_all[i] * ρ_all[i, :] * T₂_all[i, :] + SIF_all[i, :] * T₁_all[i, :];

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

# ╔═╡ Cell order:
# ╟─02bfa302-b66d-11f0-2f66-5d3686c10c23
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
# ╟─781e7324-2768-4f17-a0cc-96f0bf5d078e
# ╟─7ce2ad51-4d30-4008-b4b2-94e90330df65
# ╠═fceee754-f52a-4725-b860-95a4c95129bf
# ╠═91f48071-853a-4701-a422-896a350482d0
# ╠═b42dc5f2-a93f-4254-8601-ea829697d180
# ╟─b1cc4469-ce38-4dbb-908e-7cc1f0f6a630
