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
## Create pseudo measurement - Understand transmittance construction
Created on 2025-11-13
"""

# ╔═╡ 49880075-3490-4a9d-81a7-57109a7602f5
# wavelenth
λ_min = 620.; λ_max = 860.;

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
	SIF_shape_dict = JLD2.load("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/SIF_singular_vector.jld2")
	
	# Create interpolator
	itp = interpolate(SIF_shape_dict["SIF_shapes"], (BSpline(Linear()), NoInterp()))
	range₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end]
	range₂ = 1:size(itp, 2)
	sitp = scale(itp, range₁, range₂)
	setp0 = extrapolate(sitp, 0)
	
	# Interpolate to OCI bands
	SIF_new = reduce(hcat, [setp0.(oci_band, i) for i in range₂])
	# Scale SIF
	scale_factor_SIF = 20
	SIF_new *= scale_factor_SIF
	println("SIF shapes scaled by factor $scale_factor_SIF: $(size(SIF_new, 2)) shapes")
	# interpolation in the first dimension and no interp. in the second
	itp₂    = interpolate(SIF_shape_dict["SIF_U"], (BSpline(Linear()), NoInterp()));
	# scale
	r₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end];
	r₂ = 1:size(itp₂, 2);
	sitp₂   = scale(itp₂, r₁, r₂);
	# set extrapolation filling value = 0
	setp0₂  = extrapolate(sitp₂, 0)
	# interpolation
	SIF_PC      = reduce(hcat, [setp0₂.(oci_band, i) for i in range₂]); 
	SIF_PC1_max = maximum(SIF_PC[:,1]);
	# scale SIF principle components
	SIF_PC      /= SIF_PC1_max;
	println("SIF shape interpolated")

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
	order = 4;
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
	filename = raw"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt";
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
	Random.seed!(512)
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
	    μ₁_all[i] = 1.;   # cosd(sza_noSIF[ind_sza[i]]);
		μ₂_all[i] = 1.;   # cosd(vza_noSIF[ind_vza[i]]);
	    
	    # ----- Transmittance -----
		σ₁ = @. - 1 / μ₁_all[i] * log( trans_new[indₜ₁[i], :] );
		σ₂ = @. - 1 / μ₂_all[i] * log( trans_new[indₜ₂[i], :] );
	    T₁_all[i, :] .= exp.(- σ₁);
	    T₂_all[i, :] .= exp.(- σ₁ - σ₂);
	    
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
	Δn = 200;
	
	# Create 5x1 subplot layout
	p = plot(layout=(5, 1), size=(800, 1000), legend=false)
	
	# Plot each variable
	plot!(p[1], λ, ρ_all[1:Δn:end,:]', title="Reflectance", ylabel="ρ", xlabel="")
	plot!(p[2], λ, T₁_all[1:Δn:end,:]', title="Transmittance (one-way)", ylabel="T", xlabel="")
	plot!(p[3], λ, T₂_all[1:Δn:end,:]', title="Transmittance (two-way)", ylabel="T", xlabel="")
	plot!(p[4], λ, SIF_all[1:Δn:end,:]', title="SIF", ylabel="SIF", xlabel="")
	plot!(p[5], λ, pseudo_obs_all[1:Δn:end,:]', title="TOA", ylabel="Radiance", xlabel="Wavelength (nm)")

end

# ╔═╡ c9e0820b-1052-4b3b-94e2-918fc956f9da
md"""
> ### Visualize Retrieval
"""

# ╔═╡ 8b914198-b79e-4478-a8da-6a2c43312790
begin
	println("Load LUT for cross sections...")
	
	o2_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_O2.jld2";
	o2_sitp = read_rescale(o2_jld2);
	h2o_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_H2O.jld2";
	h2o_sitp = read_rescale(h2o_jld2);
	metadata = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01.log";
	ν_grid_o2, p_grid_hPa, t_grid = o2_sitp.ranges;
	println("LUT loaded.")

	kernel_dir = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/"
	@load kernel_dir*"KernelInstrument.jld2" MyKernel
	println("Kernel instrument loaded.")
end

# ╔═╡ 65ebbfd1-ea63-4447-88e7-5d32b7437ad4
begin
	# parameters
	println("Setting up retrieval parameters...")
	nPoly  = 6
	nLayer = 1
	nSIF   = 1
	
	Sₐ   = I(nPoly+nSIF+13) .* 0.;
	# update diagonal term
	for i=1:(nPoly+1)
	    Sₐ[i,i] = 1e10;
	    # large variance applies no constrain to these polynomial term
	end
	# SIF magnitude
	Sₐ[end, end] = 1;
	# pressure terms
	index_pressure = [1, 4, 7, 10] .+ nPoly .+ 1;
	Sₐ[CartesianIndex.(index_pressure, index_pressure)] .= 1e4
	# temperature terms
	index_temperature = [2, 5, 8, 11] .+ nPoly .+ 1;
	Sₐ[CartesianIndex.(index_temperature, index_temperature)] .= 1e3;
	# vcd for O2 (log10 transformed)
	index_vcd_o2 = [3, 9] .+ nPoly .+ 1;
	Sₐ[CartesianIndex.(index_vcd_o2, index_vcd_o2)]   .= 0.01;   # variance for O2 statistically is ~0.0025002
	# vcd for H2O (log10 transformed)
	index_vcd_h2o = [6, 12] .+ nPoly .+ 1;
	Sₐ[CartesianIndex.(index_vcd_h2o, index_vcd_h2o)] .= 0.5;    # variance for H2O statistically is ~0.381829
	
	println("Sₐ =====> Diagonal terms updated, with diagonal terms: $(diag(Sₐ))")
	
	# βₐ and γₐ
	βₐ = [800.0; 250.0; 24.6; 800.0; 290.0; 22.6];              # see demo for the stats
	γₐ = βₐ;                       
	println("βₐ and γₐ ======> Set to climatology (?) values")

end

# ╔═╡ 6b25a1ec-1c83-4251-9f08-19213e0f16a4
begin
	# iteration method
	thr_Converge = 1e-4;
	MyIter = (px, MyModel) -> LM_Iteration!(
	    px, 
	    MyModel;
	    thr_Converge = thr_Converge,
	    nIter = 20,
	    γ_init = 1000.0,
	    γ⁺ = 10.0,
	    γ⁻ = 2.0,
	    max_runtime = 400.0
	);
	
	params = RetrievalParams_xSecFit(
	    λ  = λ,
	    λc = λc,
	    E  = E,
	    c₁ = c1,
	    c₂ = c2,
	    o2_sitp  = o2_sitp,
	    h2o_sitp = h2o_sitp,
	    InstrumentKernel = MyKernel,
	    forward_model    = forward_model,
	    nPoly  = nPoly,
	    nLayer = nLayer,
	    nSIF   = nSIF,
	    Sₐ = Sₐ,
	    βₐ = βₐ,
	    γₐ = γₐ,
	    SIFComp  = SIF_PC,      
	    iteration_method = MyIter,
	    thr_Converge = thr_Converge
	);
	println("Retrieval parameters set.")
end

# ╔═╡ 129372a1-cee7-4dde-bacc-7bd1f6782551
begin
	# single pixel retrieval
	i = 501;
	SinglePixel = Retrieval_for_Pixel(
	    pseudo_obs_all[i, :],
	    0.,  # sza_noSIF[ind_sza[i]],
	    0.,  # vza_noSIF[ind_sza[i]],
	    maximum(SIF_all[i, :]),
	    1.0,
	    1.0,
	    params
	);
end

# ╔═╡ 9170dcb7-6532-4088-8cfd-b66374b88a8d
md"""
> ### Reconstruct
"""

# ╔═╡ 17835714-a99f-4b0f-889c-dcb74f38647e
MyReconModel = (x, px) -> forward_model(
	x,
	px, 
	params;
	return_components=true
);


# ╔═╡ 6cb86222-0595-42eb-aa3a-c5b249fade32
SinglePixel

# ╔═╡ 09902dcb-53de-4d32-9ad1-c1242fad4401
_, ρ, T₁, T₂, SIF = MyReconModel(SinglePixel.x, SinglePixel)

# ╔═╡ c94643a2-4d4d-4028-b7d3-fb8433a5272a
begin
	# manual recon.
		
	# reflectance
    xᵨ    = SinglePixel.x[1 : SinglePixel.nPoly+1]
    v     = collectPl.(SinglePixel.λc, lmax=SinglePixel.nPoly);
    thisρ = hcat(v...)' * xᵨ;
	
end

# ╔═╡ 52f8592d-5698-4938-bcf9-95545fa692f5
begin
	plot(
		λ, ρ_all[i,:], size=(800, 300),
		xticks = (620:10:860, string.(620:10:860)),
	)
	plot!(λ, ρ)
	plot!(λ, thisρ)
end

# ╔═╡ 6100cc31-10bc-4dfe-b76f-96dfcb516f7a
plot(
	λ, [T₂_all[i,:] T₂], size=(800, 300),
	xticks = (620:10:860, string.(620:10:860)),
	title="T₂"
)

# ╔═╡ ae030ca8-d93b-4b40-8a68-ecfd9e4fdb12
plot(
	λ, [T₁_all[i,:] T₁], size=(800, 300),
	xticks = (620:10:860, string.(620:10:860)),
	title="T₁"
)

# ╔═╡ c96faaf7-4de2-40b9-9fba-81f33efa6393
plot(
	λ, [SIF_all[i,:] SIF], size=(800, 300),
	xticks = (620:10:860, string.(620:10:860)),
	title="SIF"
)

# ╔═╡ 9075df12-31b1-41c1-a94b-1d9378c110cc
plot(
	λ, [SinglePixel.R_toa SinglePixel.y], size=(800, 300),
	xticks = (620:10:860, string.(620:10:860)),
	title="TOA radiance"
)

# ╔═╡ 3512c031-896f-4ad4-b203-cf505ffe6f62
plot(
	λ, (SinglePixel.R_toa .- SinglePixel.y), size=(800, 300),
	xticks = (620:10:860, string.(620:10:860)),
	title="residual"
)

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
# ╟─31cc57c1-aed8-4de0-86bb-9c0987d0e382
# ╟─d0cd4fd3-a740-493d-a29d-99a0b9403f58
# ╟─ef0922d4-2572-4327-9499-cb828e27e5e4
# ╟─83f4ae07-fbfb-46e0-9806-1a89ff99527a
# ╠═ca262d91-7318-493b-ad0d-44e34ce3961c
# ╟─2a0ac581-5e43-44ef-a8d4-a30ccf4d1c0f
# ╟─781e7324-2768-4f17-a0cc-96f0bf5d078e
# ╟─7ce2ad51-4d30-4008-b4b2-94e90330df65
# ╠═fceee754-f52a-4725-b860-95a4c95129bf
# ╠═91f48071-853a-4701-a422-896a350482d0
# ╠═b42dc5f2-a93f-4254-8601-ea829697d180
# ╟─b1cc4469-ce38-4dbb-908e-7cc1f0f6a630
# ╟─c9e0820b-1052-4b3b-94e2-918fc956f9da
# ╠═8b914198-b79e-4478-a8da-6a2c43312790
# ╠═65ebbfd1-ea63-4447-88e7-5d32b7437ad4
# ╠═6b25a1ec-1c83-4251-9f08-19213e0f16a4
# ╠═129372a1-cee7-4dde-bacc-7bd1f6782551
# ╟─9170dcb7-6532-4088-8cfd-b66374b88a8d
# ╠═17835714-a99f-4b0f-889c-dcb74f38647e
# ╠═6cb86222-0595-42eb-aa3a-c5b249fade32
# ╠═09902dcb-53de-4d32-9ad1-c1242fad4401
# ╠═c94643a2-4d4d-4028-b7d3-fb8433a5272a
# ╟─52f8592d-5698-4938-bcf9-95545fa692f5
# ╟─6100cc31-10bc-4dfe-b76f-96dfcb516f7a
# ╟─ae030ca8-d93b-4b40-8a68-ecfd9e4fdb12
# ╟─c96faaf7-4de2-40b9-9fba-81f33efa6393
# ╟─9075df12-31b1-41c1-a94b-1d9378c110cc
# ╟─3512c031-896f-4ad4-b203-cf505ffe6f62
