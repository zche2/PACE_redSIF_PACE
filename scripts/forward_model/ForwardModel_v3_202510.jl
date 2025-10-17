### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ ae9f792f-5175-4b47-babe-8ee7e50cebe0
import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE");

# ╔═╡ 669f0127-0cfa-42e1-9325-4c309fb225a4
using Polynomials, ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics

# ╔═╡ ca59d107-6c89-46c2-bca3-2d3d604c6d27
using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools

# ╔═╡ e19451d2-ec27-4b5e-a654-94abe2387e16
using JLD2, Interpolations

# ╔═╡ 7c1796ad-0654-4642-b615-035099743c27
include("/home/zhe2/FraLab/PACE_redSIF_PACE/PACE_SIF.jl")

# ╔═╡ bcf01e82-a621-11f0-24e9-db4a3bfd47b5
md"""
> ## Forward model v3
- Fit trans. separately （subtract by 1 / log transformed before performing SVD）
- Fit spectral shape of SIF
- G-N iteration scheme
- Low-degree polynomial: conceptually reproduce nFLH 
🔵 Whole range poly-fit / fit to several baseline wavelengths? - a priori estimation based ONLY on baseline wavelengths and set cov. diag. term to VERY small value.
"""

# ╔═╡ d3e33a05-871c-49a0-a91c-c5df2028d96f
md"""
### Fitting window & data
----
- PACE TOA rad.
- trans. from MERRA2
- SNR
- SIF (has to be interpolated to PACE wv.)
"""

# ╔═╡ a9817961-e3ab-42f0-8b7c-580d4a39bbdd
λ_min = 620.; λ_max = 755.;

# ╔═╡ 3060b59c-55b0-4e35-ac8c-5d7e2279fbd2
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
	
	# SVD
	HighResSVD = PACE_SIF.Spectral_SVD(trans .- 1., bands, λ_min=λ_min, λ_max=λ_max);

	# PACE data
	oci = Dataset(
	"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/sample_granule_20250501T183011_new_chl.nc");
	red_band = oci["red_wavelength"][:];
	nflh     = oci["nflh"][:, :];
	vza      = oci["sensor_zenith"][:, :];
	sza      = oci["solar_zenith"][:, :];
	println("\nRead in PACE Dataset")

	# select band (continuum spectrum)
	ind      = findall( λ_min .< red_band .< λ_max );
	E        = oci["red_solar_irradiance"][ind];
	R_toa    = oci["radiance_red"][:, :, ind];
	oci_band = red_band[ind];
	println("\nBand selected: $oci_band")
end

# ╔═╡ 22c13315-e686-4632-b24c-d6b35011f19d
begin
	filename = raw"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt";
	lines = readlines(filename);
	end_header_index = findfirst(x -> x == "/end_header", lines);
	data  = readdlm(filename, String, skipstart=end_header_index);

	FPA   = data[:, 1];                   # 1st column: band
	wvlen = parse.(Float64, data[:, 2]);  # 2nd column: center wavelength
	c1    = parse.(Float64, data[:, 4]);  # 4th column: c1
	c2    = parse.(Float64, data[:, 5]);  # 5th column: c2

	wv_val  = (λ_min .< wvlen .< λ_max);
	snr_ind = findall((FPA .== "Red") .& wv_val);
	println("SNR retrieved")
end

# ╔═╡ deb62e05-2d00-4704-8288-be7e292e5a5b
begin
	# load SIF
	SIF_shape_dict = JLD2.load("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/SIF_singular_vector.jld2")
	println("SIF data loaded")
end

# ╔═╡ e001c698-ddd3-4a61-9205-25c2e523e28d
begin
	# interpolation in the first dimension and no interp. in the second
	itp    = interpolate(SIF_shape_dict["SIF_U"], (BSpline(Linear()), NoInterp()));
	
	# scale
	range₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end];
	range₂ = 1:size(itp, 2);
	sitp   = scale(itp, range₁, range₂);

	# set extrapolation filling value = 0
	setp0  = extrapolate(sitp, 0)
	println("SIF shape interpolated")
end

# ╔═╡ b6c3c640-e5fe-4322-8dc6-0e4a731704e7
SIF_shape_dict

# ╔═╡ 2bbbc4d9-c869-4dbf-8611-023d62ec24b0
begin
	# SIF_new = [sitp.(λᵢ, range₂) for λᵢ in oci_band]; 
	SIF_new = reduce(hcat, [setp0.(oci_band, i) for i in range₂]); 
	# original
	plot(range₁, SIF_shape_dict["SIF_U"][:,1], label="SIF₁", linewidth=3.)
	plot!(range₁, SIF_shape_dict["SIF_U"][:,2], label="SIF₂", linewidth=3.)
	# interpolated
	scatter!(oci_band, SIF_new[:,1], label="SIF₁ - interpolated", markersize=3.)
	scatter!(oci_band, SIF_new[:,2], label="SIF₂ - interpolated", markersize=3.)
	# xlims
	xlims!(λ_min, λ_max)
end

# ╔═╡ 3ef82576-1699-4d31-8759-2462d2baf6de
begin
	# find transmittance baseline
	bl_wvlen = [607.99, 610.36, 612.73, 615.14, 617.6, 620.06, 622.53, 669.52, 670.76, 671.99, 673.24, 674.51, 675.73, 676.96, 678.21, 679.45, 754.3, 779.33, 867.11, 869.61, 872.13];
	
	# [668.265, 669.518, 670.755, 671.99, 673.245, 674.505, 675.73, 676.962, 678.205, 679.445, 680.68, 751.79, 753.04, 754.295, 776.832, 779.335, 867.115, 869.615, 872.13];

	bl_ind = map(bl_wvlen -> argmin(abs.(oci_band .- bl_wvlen)), bl_wvlen);
	# bl_ind   = collect(1:length(oci_band));
	# println(bl_ind)

	# plot
	temp_ind = findall( λ_min .< bands .< λ_max );
	plot(oci_band, trans[1:400:end, temp_ind]',
		label="",
		size=(600, 250)
	);
	scatter!(oci_band[bl_ind], trans[1, temp_ind[bl_ind]],
		label="baseline pts",
		markersize=2,
	)
	xlims!(λ_min, λ_max)
end

# ╔═╡ e1fe6027-a214-43f9-98cc-e8c431b1b818
findall(coalesce.((nflh .> 0.3) .& (nflh .< 0.4), false))

# ╔═╡ 6f3d41e1-52cd-4a99-b907-47733be3e7b6
md"""
### Forward model 1: Polyfit + transmittance
----
- Sₑ and Sₐ
- Forward model looks like:

$$\rho_{s}=\sum{a_jP_j},\ T(\lambda)=\sum{\beta_i P_i}$$

$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)T(\lambda)}{\pi}$$

	p.s. Starting with an indent shows like this!

- Iterations & reconstruct
"""

# ╔═╡ 5115015a-4479-48aa-80b4-509deb9d038e
begin
	pixel  = 916; scan = 38;
	R_px   = R_toa[pixel, scan, :];
	sza_px = sza[pixel, scan];
	vza_px = vza[pixel, scan];
	n      = 3;
	nPC    = 8;
	# Se
	noise = sqrt.( c1[snr_ind] .+ c2[snr_ind] .* R_px);
	Se    = Diagonal(noise.^2);
	println("measurement error set according to the pixel!")
end

# ╔═╡ 0e24b3a7-cebe-4a10-a37f-f92d37744508
md"""
##### Towards priori s.d. in nPC
the s.d. of 1st PC is usually not that large.
"""

# ╔═╡ 3bcdf18f-17a3-4664-b91e-0c6604825c5d
begin
	println("priori error set!")
	# s.d. of the total loading
	loading_sd = std(HighResSVD.VarExp .* HighResSVD.Loading, dims=2);
	# priori cov
	Sₐ   = I(n+nPC+1) .+ 0.;
	# uodate diagonal term
	for i=1:(n+1)
	    Sₐ[i,i] = 1e-11;
		# large variance applies no constrain to these polynomial term
	end
	for i=(n+2):(n+nPC+1)
	    # Sₐ[i,i] = rel_error .* HighResSVD.VarExp[i - (n+1)];  # rel_error = 0.001
		Sₐ[i,i] = loading_sd[i - (n+1)];
	end
end

# ╔═╡ 4cb8bb3e-2b0a-4cb9-a37f-332f8c9da26d
# center wavelength `oci_band`
λc = PACE_SIF.center_wavelength(oci_band);

# ╔═╡ f6762889-f05d-4ff2-b77e-1eb3ac2d3a91
md"""
##### Struct & more functions
"""

# ╔═╡ 8ce8f9ba-718c-4e60-aaa1-7aba7fe374f4
mutable struct Pixel
	# universal for the granule
	λ      # fitting window
	E      # observed extraterrestrial irradiance
	nPoly  # degree of Legendre Polynomial
	nPC    # number of transmittance basis function included
	nSIF   # number of SIF PCs
	"a priori matrix"
	Sa    
	"PCs specified, = HighResSVD.PrinComp[:, 1:nPC]"
	trans_mat 
	"SIF shape specified"
	SIF_shape
	"centered wavelength for fast computation of Legendre Polys, = center_wavelength(λ)"
	λc        
	"baseline band for scaling transmittance, = find_baseline_band(λ)"
	λ_bl_ind

	# pixel L1B measurement & set up
	"TOA radiance"
	R_toa
	"solar zenith angle"
	sza
	"viewing zenith angle"
	vza
	"measurement error"
	Se
	"flag: 0 - not doing retrieval due to bad input data, refer to l2flag / nflh"
	flag   # 🔴 not necessarily need?
	"a priori estimation"
	xₐ
	
	# retrieval
	"retrieved state vector"
	x
	"modelled radiance"
	y
	"convergence flag"
	ΔRMSE
	"iteration label"
	iter_label

	# Inner constructer
	function Pixel()
        new()
    end
end

# ╔═╡ caa101ef-a81c-4b43-96c0-a2a6810ea756
function forward_model1(
		x,
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	rho   = hcat(v...)' * x[1 : px.nPoly+1];
	
	# transmittance
	T      = px.trans_mat * x[(px.nPoly+2):(px.nPoly+px.nPC+1)];
	# println(typeof(T))
	T_norm = @. T + 1.0 # PACE_SIF.scale_transmittance(T, px.λ_bl_ind);

	# one way vs. two way
	# T2_norm = two_way_trans(T_norm, px.sza, px.vza);
	# SIF magnitude
	# SIF    = @. x[px.nPoly+px.nPC+px.nSIF+1] * SIF_shape(px.λ);
	
	# TOA radiance
	rad    = @. px.E * cosd(px.sza) / pi * T_norm * rho;
	return rad
end

# ╔═╡ c458c47e-e272-4ee9-8936-1ae11a93998c
function Jacobian(x, model; len=length(oci_band))
	res = DiffResults.JacobianResult(zeros(len), x);
	ForwardDiff.jacobian!(res, model, x);
	K   = DiffResults.jacobian(res);
	val = DiffResults.value(res);
	return K, val
end

# ╔═╡ 7f33d270-b1be-49bf-9c66-57598a891901
function GainMatrix(K; Se=Se, Sa=Sa)
	return inv( K'inv(Se)K + inv(Sa) )K'inv(Se)
end

# ╔═╡ e55eacaa-beb5-439f-bcf2-fc647f69058b
function GN_Interation!(
			px :: Pixel;
			model = forward_model1,
			nIter = 20,
			thr_Converge = 1e-8,
		)
	
	# initial
	xₐ = px.xₐ;   # priori estimation
	xₙ = px.x;
	Kₙ, _ = Jacobian(xₙ, x -> model(x, px));
	# k     = px.iter_label;       # number of iterations
	RMSE₀ = 1e20; 
	RMSE₁ = PACE_SIF.root_mean_square(px.R_toa, px.y);
	ΔRMSE = RMSE₁ - RMSE₀

	# loop
	while ( abs(ΔRMSE) > thr_Converge ) & ( px.iter_label < nIter )
		# k += 1
		# get Gain matrix
		Gₙ     = GainMatrix(Kₙ, Se=px.Se, Sa=px.Sa);
		# retrieval
		xₙ₊₁   = xₐ .+ Gₙ * (px.R_toa .- px.y .+ Kₙ * ( px.x .- xₐ ) );
		# update x and y
		px.x   = xₙ₊₁;
		Kₙ₊₁, yₙ₊₁ = Jacobian(xₙ₊₁, x -> model(x, px));
		px.y   = yₙ₊₁;
		# @show size(Kₙ₊₁)
		Kₙ     = Kₙ₊₁;
		# iter ++
		px.iter_label += 1;
		# test convergence
		RMSE₀  = RMSE₁;
		RMSE₁  = PACE_SIF.root_mean_square(px.R_toa, px.y);
		ΔRMSE  = RMSE₁ - RMSE₀;
		px.ΔRMSE = ΔRMSE;
	end
	
	return nothing
end

# ╔═╡ 8f13bb43-45f9-4e93-8fde-1c23b78964e6
function reconstruct1(
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	rho   = hcat(v...)' * px.x[1 : px.nPoly+1];
	
	# transmittance
	T      = px.trans_mat * px.x[(px.nPoly+2):(px.nPoly+px.nPC+1)];
	T_norm = @. T + 1.0; # PACE_SIF.scale_transmittance(T, px.λ_bl_ind);
	
	return rho, T_norm
end

# ╔═╡ 9e09c9bc-691e-4d53-ad59-1f6fb0a32fff
md"""
##### Retrieval
"""

# ╔═╡ 65dc5736-7ffb-4d44-a921-3c541647f196
begin
	MyPixel   = Pixel();

	# Step1: construct struct
	MyPixel.λ  = oci_band;
	MyPixel.λc = λc;
	MyPixel.λ_bl_ind = bl_ind;
	MyPixel.E     = E;
	MyPixel.nPoly = n;
	MyPixel.nPC   = nPC;
	MyPixel.nSIF  = 0;
	MyPixel.Sa    = Sₐ;
	MyPixel.trans_mat = HighResSVD.PrinComp[:, 1:MyPixel.nPC];

	MyPixel.R_toa = R_px;
	MyPixel.sza   = sza_px;
	MyPixel.vza   = vza_px;
	MyPixel.Se    = Se;
	MyPixel.flag  = 1.;       # quite temporary..

	# a priori estimation
	K₀  = 
		MyPixel.E[MyPixel.λ_bl_ind] .* cosd(MyPixel.sza) ./ pi .* hcat(collectPl.(MyPixel.λc[MyPixel.λ_bl_ind], lmax=MyPixel.nPoly)...)';
	# K₀  = 
	# 	MyPixel.E .* cosd(MyPixel.sza) ./ pi .* hcat(collectPl.(MyPixel.λc, lmax=MyPixel.nPoly)...)';
	G₀  = inv( K₀'K₀ )K₀';
	x₀  = G₀ * MyPixel.R_toa[MyPixel.λ_bl_ind];
	tmp = zeros(MyPixel.nPC + MyPixel.nSIF - 2) .+ .001;
	println(x₀)

	MyPixel.xₐ = [x₀... 0. 0.0 tmp...]';

	# set-up
	MyPixel.x  = MyPixel.xₐ;
	MyPixel.y  = forward_model1(MyPixel.x, MyPixel);
	MyPixel.iter_label = 0;

	# Step2: iteration
	GN_Interation!(MyPixel)
end

# ╔═╡ 38daa8c9-2bde-49a1-b446-5aae7a480e4c
# measurement error vs. priori error: which one matters more?
begin
	K, _ = Jacobian(MyPixel.x, x -> forward_model1(x, MyPixel));
	println(diag(K'inv(Se)K), diag(inv(Sₐ)))
	inv( K'inv(Se)K + inv(Sₐ) )
end

# ╔═╡ 957559fd-4e88-468e-9b64-a9cf5964d5bf
rho, T = reconstruct1(MyPixel);

# ╔═╡ c84dc329-0e51-4eb1-9f58-6cae728dc230
begin
	TheTitle = "pixel=$pixel, scan=$scan, nFLH=$(round(nflh[pixel, scan], digits=2))\n nPoly=$(MyPixel.nPoly), nPC=$(MyPixel.nPC)"
	
	plot(
		MyPixel.λ, MyPixel.R_toa,
		label="Observations", linewidth=2,
		color=:grey,
		xlabel="Wavelength [nm]",
		ylabel="TOA radiance \n [W/m²/µm/sr]",
		xlabelfontsize=8,
		ylabelfontsize=8,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title=TheTitle,
		titlefontsize=10,
		size=(600, 300)
		);
	plot!(
		MyPixel.λ, forward_model1(MyPixel.xₐ, MyPixel), label="initial guess"); 
	plot!(
		MyPixel.λ, MyPixel.y, label="@ convergence", linewidth=1);
	plot!(
		MyPixel.λ, MyPixel.E .* rho ./ π .* cosd(MyPixel.sza) , label="ρ", linewidth=1);
	scatter!(oci_band[bl_ind], MyPixel.R_toa[bl_ind], markersize=1.5);
end

# ╔═╡ 2c95ba4c-c0a5-491e-967b-d896814b91f7
begin	
	r1 = plot(
		MyPixel.λ, MyPixel.R_toa .- MyPixel.y, label="Residual (W/m²/µm/sr)", linewidth=1.5, color=:grey);
	title!(TheTitle, titlefontsize=10)
	r2 = plot(MyPixel.λ, (MyPixel.R_toa .- MyPixel.y)./MyPixel.R_toa * 100, 		   label="Relative Residual (%)", linewidth=1.5,
		color=:grey);
	plot(r1, r2, layout=(2,1), size=(600, 400))
end

# ╔═╡ 9c28e65b-9522-47b3-85c0-c7c47852cf72
begin
	rho_fig   = plot(MyPixel.λ, rho, label="surface reflectance",
					title=TheTitle,
					titlefontsize=10,
					)
	trans_fig = plot(MyPixel.λ, T, label="T↑↓")
	plot(rho_fig, trans_fig,
		 layout=(2,1),
		 size=(600, 450)
	)
end

# ╔═╡ 30b65205-de3b-4729-9287-40a9022bb9f2
md"""
### Apply to more pixels to see if there's converging patterns
----
of TOA rad, reflectance, residuals...
"""

# ╔═╡ 6eee7eae-dcc3-4eaf-902a-5436ec04677b
function woSIF_Retrieval(
		# "L1B pixel-by-pixel vals"
		R_toa ,
		sza, vza, flag,   # quite temporary..
		params,
	)
	
	# preprocess: if the flag is false, not doing the retrieval
	if ismissing(flag)
		return missing
	end
	
	MyPixel       = Pixel();
	forward_model = params.forward_model
	nIter         = params.nIter
	thr_Converge  = params.thr_Converge

	# Step1: construct struct
	MyPixel.λ  = params.λ;
	MyPixel.λc = params.λc;
	MyPixel.λ_bl_ind = params.λ_bl_ind;
	MyPixel.E     = params.E;
	MyPixel.nPoly = params.nPoly;
	MyPixel.nPC   = params.nPC;
	MyPixel.nSIF  = params.nSIF;
	MyPixel.Sa    = params.Sₐ;
	MyPixel.trans_mat = HighResSVD.PrinComp[:, 1:MyPixel.nPC];

	MyPixel.R_toa = R_toa;
	MyPixel.sza   = sza;
	MyPixel.vza   = vza;
	noise         = sqrt.( c1[snr_ind] .+ c2[snr_ind] .* MyPixel.R_toa);
	MyPixel.Se    = Diagonal(noise.^2);
	MyPixel.flag  = flag; 
	
	# a priori estimation
	K₀  = 
		MyPixel.E[MyPixel.λ_bl_ind] .* cosd(MyPixel.sza) ./ pi .* hcat(collectPl.(MyPixel.λc[MyPixel.λ_bl_ind], lmax=MyPixel.nPoly)...)';
	# K₀  = 
	# 	MyPixel.E .* cosd(MyPixel.sza) ./ pi .* hcat(collectPl.(MyPixel.λc, lmax=MyPixel.nPoly)...)';
	G₀  = inv( K₀'K₀ )K₀';
	x₀  = G₀ * MyPixel.R_toa[MyPixel.λ_bl_ind];
	# x₀  = G₀ * MyPixel.R_toa;
	tmp = zeros(MyPixel.nPC + MyPixel.nSIF - 2) .+ .001;
	MyPixel.xₐ = [x₀... -1. 0.1 tmp...]';

	# set-up
	MyPixel.x  = [x₀... -1. 0.1 tmp...]';
	MyPixel.y  = forward_model(MyPixel.x, MyPixel);
	MyPixel.iter_label = 0;
	
	# Step2: iteration
	GN_Interation!(
		MyPixel, 
	    model=forward_model,
	    nIter=nIter,
	    thr_Converge=thr_Converge
	)

	# Step3: return
	# return if converge
	if abs(MyPixel.ΔRMSE) < thr_Converge
		# println("successfully retrieved")
		return MyPixel
	else
		return missing
	end
end

# ╔═╡ 6ca8f4f5-c063-457d-8dcb-6b7e66b547eb
begin
	# fix the params
	params = (
		λ = oci_band, 
		E = E,
		λ_bl_ind = bl_ind,
		λc    = λc,
		nPoly = n,
		nPC   = nPC, 
		nSIF  = 0,
		nIter = 20, 
		Sₐ    = Sₐ,
		thr_Converge = 1e-8,
		forward_model = forward_model1,
	)
end

# ╔═╡ 18505102-dd32-4e5a-bc2e-80f77a0e30be
begin
	SIF_index    = findall(coalesce.((nflh .> 0.2) .& (nflh .< 0.4), false));
	number_of_px = size(SIF_index)[1];
	SIF_683 = woSIF_Retrieval.(
	    eachslice(R_toa[SIF_index, :], dims=1),  # rather than dims=(1,2)
		sza[SIF_index],
		vza[SIF_index],
		nflh[SIF_index], 
	    Ref(params) 
	)
end

# ╔═╡ 373cf733-983e-4d68-82f6-473295c081ba
begin
	# make some plot
	Δn = 100;
	p = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="TOA radiance \n [W/m²/µm/sr]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)
	for i in 1:Δn:number_of_px
		plot!(p, oci_band, SIF_683[i].R_toa)
	end
	p
end

# ╔═╡ 451294b0-3451-4034-bd86-22254bcbb561
begin
	# make some plot
	p_rho = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="reflectance [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)
	
	for i in 1:Δn:number_of_px
		# reconstruct
		rhoᵢ, _ = reconstruct1(SIF_683[i]);
		plot!(p_rho, oci_band, rhoᵢ)
	end
	p_rho
end

# ╔═╡ 226ee786-1aad-4576-8797-0a340e2a5142
begin
	# make transmittance
	p_trans = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="transmittance [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)
	
	for i in 1:Δn:number_of_px
		# reconstruct
		_, Tᵢ = reconstruct1(SIF_683[i]);
		plot!(p_trans, oci_band, Tᵢ)
	end
	p_trans
end

# ╔═╡ 269abad8-70df-433d-805b-64d762dacadc
begin
	# residual
	p_resd = plot(
		size=(800, 400), 
		legendcolumns=4,
		xlabel="[nm]",
		ylabel="Residual",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)
	
	for i in 1:Δn:number_of_px
		# reconstruct
		resdᵢ = SIF_683[i].y .- SIF_683[i].R_toa
		plot!(p_resd, oci_band, resdᵢ,
			  label="$(round(SIF_683[i].flag, digits=2))",
		)
	end
	p_resd
end

# ╔═╡ ef592e81-291c-4ffc-8f8a-b30f3a39dc1c
md"""
### Forward model 2: Polyfit + transmittance + SIF shape fixed
----
- Fitting $T_{\downarrow\uparrow}$ and $T_{\uparrow}$ separately:

$$\rho_{s}(\lambda)=\sum{a_jP_j},\ T_{\downarrow\uparrow}(\lambda)=\sum{\beta_i P_i}, \ T_{\uparrow}(\lambda)=\sum{\gamma_i P_i}$$

$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)T_{\downarrow\uparrow}(\lambda)}{\pi} + SIF(\lambda)T_{\uparrow}(\lambda)$$

"""

# ╔═╡ f76d7c40-1ff4-4e39-bbae-cd7dd6834b7a
begin
	println("priori error set! (nPC×2)")
	# priori cov
	Sₐ₂   = I(n+nPC*2+2) .+ 0.;
	# update diagonal term
	for i=1:(n+1)
	    Sₐ₂[i,i]  = 1e10;
		# large variance applies no constrain to these polynomial term
	end
	# \beta
	for i=(n+2):(n+nPC+1)
	    # Sₐ[i,i] = rel_error .* HighResSVD.VarExp[i - (n+1)];  # rel_error = 0.001
		Sₐ₂[i,i]  = loading_sd[i - (n+1)];
	end
	# \gamma
	for i=(n+nPC+2):(n+nPC*2+1)
		Sₐ₂[i,i]  = loading_sd[i - (n+nPC+1)];
	end
	# SIF magnitude
	Sₐ₂[end, end] = 1;
	println("Diagonal terms are: $(diag(Sₐ₂))")
end

# ╔═╡ 2317ce66-50e6-4ae2-9796-27e0f1cd1534
function forward_model2(
		x,
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	ρ     = hcat(v...)' * x[1 : px.nPoly+1];

	# T↓↑ transmittance for solar irradiance
	T2      = px.trans_mat * x[(px.nPoly+2):(px.nPoly+px.nPC+1)];
	T2_norm = @. T2 + 1.0;   # PACE_SIF.scale_transmittance(T2, px.λ_bl_ind);

	# # T↑ transmittance for SIF
	T1       = px.trans_mat * x[(px.nPoly+px.nPC+2):(px.nPoly+px.nPC*2+1)];
	T1_norm  = @. T1 + 1.0;  # PACE_SIF.scale_transmittance(T1, px.λ_bl_ind);

	# # what if T↑ is just T2 scaled by some factors?
	# T1_norm = exp.( x[px.nPoly+px.nPC+2] * log.(T2_norm) )
	
	# SIF magnitude
	SIF      = px.SIF_shape * x[px.nPoly+px.nPC*2+px.nSIF+1];
	
	# TOA radiance
	rad      = @. px.E * cosd(px.sza) / π * T2_norm * ρ + SIF * T1_norm;
	return rad
end

# ╔═╡ c2ae0783-0623-45f1-a009-c824148cd59d
begin
	# define new parameter set
	params₂ = (
		λ = oci_band, 
		E = E,
		λ_bl_ind = bl_ind,
		λc    = λc,
		nPoly = n,
		nPC   = nPC, 
		nSIF  = 1,
		nIter = 20, 
		Sₐ    = Sₐ₂,
		thr_Converge = 1e-8,
		forward_model = forward_model2,
	);
end

# ╔═╡ 64eb2c06-d043-45dc-b60f-b44a3e777db8
function redSIF_Retrieval(
		# "L1B pixel-by-pixel vals"
		R_toa ,
		sza, vza, flag,   # quite temporary..
		params,
	)
	
	# preprocess: if the flag is false, not doing the retrieval
	if ismissing(flag)
		return missing
	end
	
	MyPixel       = Pixel();
	forward_model = params.forward_model
	nIter         = params.nIter
	thr_Converge  = params.thr_Converge

	# Step1: construct struct
	MyPixel.λ  = params.λ;
	MyPixel.λc = params.λc;
	MyPixel.λ_bl_ind = params.λ_bl_ind;
	MyPixel.E     = params.E;
	MyPixel.nPoly = params.nPoly;
	MyPixel.nPC   = params.nPC;
	MyPixel.nSIF  = params.nSIF;
	MyPixel.Sa    = params.Sₐ;
	MyPixel.trans_mat = HighResSVD.PrinComp[:, 1:MyPixel.nPC];
	MyPixel.SIF_shape = SIF_new[:, 1:MyPixel.nSIF]

	MyPixel.R_toa = R_toa;
	MyPixel.sza   = sza;
	MyPixel.vza   = vza;
	noise         = sqrt.( c1[snr_ind] .+ c2[snr_ind] .* MyPixel.R_toa);
	MyPixel.Se    = Diagonal(noise.^2);
	MyPixel.flag  = flag; 

	# a priori estimation
	K₀  = 
		MyPixel.E[MyPixel.λ_bl_ind] .* cosd(MyPixel.sza) ./ pi .* hcat(collectPl.(MyPixel.λc[MyPixel.λ_bl_ind], lmax=MyPixel.nPoly)...)';
	# K₀  = 
	# 	MyPixel.E .* cosd(MyPixel.sza) ./ pi .* hcat(collectPl.(MyPixel.λc, lmax=MyPixel.nPoly)...)';
	G₀  = inv( K₀'K₀ )K₀';
	x₀  = G₀ * MyPixel.R_toa[MyPixel.λ_bl_ind];
	tmp_trans = zeros(MyPixel.nPC - 2) .+ .001;
	tmp_sif   = zeros(MyPixel.nSIF) .+ .001
	MyPixel.xₐ = [x₀... -1. 0.1 tmp_trans... -1. 0.1 tmp_trans... tmp_sif...]';
	
	# set-up
	MyPixel.x  = MyPixel.xₐ;
	MyPixel.y  = forward_model(MyPixel.x, MyPixel);
	MyPixel.iter_label = 0;
	
	# Step2: iteration
	GN_Interation!(
		MyPixel, 
	    model=forward_model,
	    nIter=nIter,
	    thr_Converge=thr_Converge
	)

	# Step3: return
	# return if converge
	if abs(MyPixel.ΔRMSE) < thr_Converge
		# println("successfully retrieved")
		return MyPixel
	else
		return missing
	end
end

# ╔═╡ 16af025b-e817-402e-8441-24e0e94761cf
Retrieval2 = redSIF_Retrieval.(
	eachslice(R_toa[SIF_index, :], dims=1),  # rather than dims=(1,2)
	sza[SIF_index],
	vza[SIF_index],
	nflh[SIF_index], 
	Ref(params₂) 
)

# ╔═╡ 0fac32c1-35bf-4fc5-a62c-68d904bedc73
function reconstruct2(
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	ρ     = hcat(v...)' * px.x[1 : px.nPoly+1];
	
	# T↓↑ transmittance for solar irradiance
	T2      = px.trans_mat * px.x[(px.nPoly+2):(px.nPoly+px.nPC+1)];
	T2_norm = @. T2 + 1.0;   # PACE_SIF.scale_transmittance(T2, px.λ_bl_ind);

	# T↑ transmittance for SIF
	T1       = px.trans_mat * px.x[(px.nPoly+px.nPC+2):(px.nPoly+px.nPC*2+1)];
	T1_norm  = @. T1 + 1.0;  # PACE_SIF.scale_transmittance(T1, px.λ_bl_ind);

	# # what if T↑ is just T2 scaled by some factors?
	# T1_norm = exp.( px.x[px.nPoly+px.nPC+2] * log.(T2_norm) )
	
	# SIF magnitude
	SIF      = px.SIF_shape * px.x[px.nPoly+px.nPC*2+px.nSIF+1];
	
	return ρ, T2_norm, T1_norm, SIF
end

# ╔═╡ 43e06ecb-2e12-4698-bffc-9fd18dde13ea
begin
	# make transmittance
	p_rho₂ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="two-way transmittance [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)
	
	p_trans₂ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="two-way transmittance [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)

	p_trans₁ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="one-way transmittance [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)

	p_SIF₂ = plot(
		size=(800, 400), 
		legendcolumns=4,
		xlabel="[nm]",
		ylabel="SIF",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		label=:outerbottom,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)
	
	for i in 1:(Δn):number_of_px
		if ismissing(Retrieval2[i])
			continue
		end
		# reconstruct
		ρⱼ, T₂ⱼ, T₁ⱼ, SIFⱼ = reconstruct2(Retrieval2[i]);
		plot!(
			p_rho₂, oci_band, ρⱼ, label="$(round(Retrieval2[i].flag, digits=2))")
		plot!(
			p_trans₂, oci_band, T₂ⱼ, label="$(round(Retrieval2[i].flag, digits=2))")
		plot!(
			p_trans₁, oci_band, T₁ⱼ, label="$(round(Retrieval2[i].flag, digits=2))")
		plot!(
			p_SIF₂, oci_band, SIFⱼ .* T₁ⱼ,
			label="",
			linestyle=:dash,
		)
		plot!(
			p_SIF₂, oci_band, SIFⱼ, label="")
	end
end

# ╔═╡ f27387bf-ceb1-4d57-b8f3-3597044962e3
p_rho₂

# ╔═╡ 1ec5f75d-e017-436c-8f53-1ad9469fbd99
p_trans₂

# ╔═╡ 873ebdde-8260-4fb2-8c02-88a0b722cf44
begin
	# residual
	p_resd₂ = plot(
		size=(800, 400), 
		legendcolumns=4,
		xlabel="[nm]",
		ylabel="Residual",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)
	
	for i in 1:Δn:number_of_px
		if ismissing(Retrieval2[i])
			continue
		end
		# get spectral-wise residual
		resdᵢ = Retrieval2[i].y .- Retrieval2[i].R_toa
		plot!(p_resd₂, oci_band, resdᵢ,
			  label="$(round(Retrieval2[i].flag, digits=2))",
		)
	end
	p_resd₂
end

# ╔═╡ 53af5378-fed0-45ce-98ae-a2b5c32c0421
p_trans₁

# ╔═╡ 906167fc-9516-499a-8fb5-229b774dc549
p_SIF₂

# ╔═╡ ba6a2fc3-6ea6-41af-82f3-bddbfe3dd0ce
begin
	# compare background radiance term
	p_trans_rho = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="ρ × T",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ρ vs. ρ×T"
	)
	
	for i in 1:Δn:number_of_px
		if ismissing(Retrieval2[i])
			continue
		end
		# reconstruct
		ρᵢ, Tᵢ = reconstruct1(SIF_683[i]);
		ρⱼ, T₂ⱼ, _, _ = reconstruct2(Retrieval2[i]);
		plot!(
			p_trans_rho, oci_band, ρⱼ .* T₂ⱼ, 
			linestyle=:dash
		)
		plot!(p_trans_rho, oci_band, ρⱼ)
	end
	p_trans_rho
end

# ╔═╡ 16a28e7f-f3f1-4496-ad11-efc6af48992e
begin
	# compare background radiance term
	p_Rtoa_rho = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="Radiance",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ρ vs. R_toa"
	)

	v = collectPl.(λc, lmax=n);
	
	for i in 1:(5*Δn):number_of_px
		if ismissing(Retrieval2[i])
			continue
		end
		# reconstruct
		ρⱼ, T₂ⱼ, _, _ = reconstruct2(Retrieval2[i]);
		
		# reconstructed by the priori
		cᵢ = Retrieval2[i].xₐ[1:params.nPoly+1];
		ρ₀ = hcat(v...)' * cᵢ;
		
		# plot
		# plot!(
		# 	p_Rtoa_rho, oci_band, ρⱼ .* E ./ π .* cosd(Retrieval2[i].sza), 
		# 	linestyle=:dash
		# )
		plot!(
			p_Rtoa_rho, oci_band, ρ₀ .* E ./ π .* cosd(Retrieval2[i].sza), 
			linestyle=:dashdot
		)
		plot!(p_Rtoa_rho, oci_band, Retrieval2[i].R_toa)
	end
	p_Rtoa_rho
end

# ╔═╡ cf673611-4324-464f-aec8-b56aaf97afa7
md"""
### Forward model 3: Polyfit + transmittance + SIF shape fixed
----
Joint fit of $T_{\downarrow\uparrow}$ and $T_{\uparrow}$ by tuning a factor (not neccesarily SVF as before):

$$\rho_{s}(\lambda)=\sum{a_jP_j}, \ T_{\uparrow}(\lambda)=\sum{\beta_i P_i}$$

$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)T_{\downarrow\uparrow}(\lambda)}{\pi} + SIF(\lambda)T_{\uparrow}(\lambda)$$

$$T_{\downarrow\uparrow}(\lambda)=exp(\gamma \times ln(T_{\uparrow}(\lambda)))$$
Where γ is correction factor accounting for （1) light path (VZA and SZA) and （2) upper atmosphere reflectance.
"""

# ╔═╡ 3509ff7b-9a8d-413f-b5d3-96d83a936a32
begin
	println("priori error set! (nPC×2)")
	# priori cov: 3 for 1 constant term, 1 SIF term, and 1 factor
	Sₐ₃   = I(n+nPC+3) .+ 0.;
	# update diagonal term
	for i=1:(n+1)
	    Sₐ₃[i,i] = 1e10;
		# large variance applies no constrain to these polynomial term
	end
	# \beta
	for i=(n+2):(n+nPC+1)
		Sₐ₃[i,i]  = loading_sd[i - (n+1)];
	end
	# \gamma
	Sₐ₃[n+nPC+2, n+nPC+2] = 2;
	# SIF magnitude
	Sₐ₃[end, end] = 1;
	println("Diagonal terms are: $(diag(Sₐ₃))")

	# add a priori estimation of \beta, \gamma, and SIF
	tmp₃ = zeros(nPC-2) .+ .001;
	tmp₃ = [-1.0 0. tmp₃... 1. 0.1]'
	println(tmp₃)
end

# ╔═╡ 798ad5dd-94e9-4f8d-a688-401a8b9e0d1a
function forward_model3(
		x,
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	ρ     = hcat(v...)' * x[1 : px.nPoly+1];

	# T↑ transmittance for SIF
	T₁    = (px.trans_mat * x[(px.nPoly+2):(px.nPoly+px.nPC+1)]) .+ 1.0;

	# T↓↑ transmittance for SIF
	T₂    = @. exp( x[px.nPoly+px.nPC+2] * log(T₁) );

	# SIF magnitude
	SIF   = px.SIF_shape * x[px.nPoly+px.nPC+px.nSIF+2];
	
	# TOA radiance
	rad   = @. px.E * cosd(px.sza) / π * T₂ * ρ + SIF * T₁;
	return rad
end

# ╔═╡ 267ed111-0136-46af-af85-742f303b17bf
begin
	# define new parameter set
	params₃ = (
		λ = oci_band, 
		E = E,
		λ_bl_ind = bl_ind,
		λc    = λc,
		nPoly = n,
		nPC   = nPC, 
		nSIF  = 1,
		nIter = 20, 
		Sₐ    = Sₐ₃,
	    xₐ    = tmp₃,
		thr_Converge = 1e-8,
		forward_model = forward_model3,
	);
end

# ╔═╡ c8e554dc-3621-4a5d-9d75-394ee1c1a56b
function Retrieval(
		# "L1B pixel-by-pixel vals"
		R_toa, sza, vza, flag,
		params,
	)
	
	# preprocess: if the flag is false, not doing the retrieval
	if ismissing(flag)
		return missing
	end
	
	MyPixel       = Pixel();
	forward_model = params.forward_model
	nIter         = params.nIter
	thr_Converge  = params.thr_Converge

	# Step1: construct struct
	MyPixel.λ  = params.λ;
	MyPixel.λc = params.λc;
	MyPixel.λ_bl_ind = params.λ_bl_ind;
	MyPixel.E     = params.E;
	MyPixel.nPoly = params.nPoly;
	MyPixel.nPC   = params.nPC;
	MyPixel.nSIF  = params.nSIF;
	MyPixel.Sa    = params.Sₐ;
	MyPixel.trans_mat = HighResSVD.PrinComp[:, 1:MyPixel.nPC];
	MyPixel.SIF_shape = SIF_new[:, 1:MyPixel.nSIF]

	MyPixel.R_toa = R_toa;
	MyPixel.sza   = sza;
	MyPixel.vza   = vza;
	noise         = sqrt.( c1[snr_ind] .+ c2[snr_ind] .* MyPixel.R_toa);
	MyPixel.Se    = Diagonal(noise.^2);
	MyPixel.flag  = flag; 

	# a priori estimation
	K₀         = 
		MyPixel.E[MyPixel.λ_bl_ind] .* cosd(MyPixel.sza) ./ pi .* hcat(collectPl.(MyPixel.λc[MyPixel.λ_bl_ind], lmax=MyPixel.nPoly)...)';
	G₀         = inv( K₀'K₀ )K₀';
	x₀         = G₀ * MyPixel.R_toa[MyPixel.λ_bl_ind];
	MyPixel.xₐ = [x₀... params.xₐ...]';
	
	# set-up
	MyPixel.x  = MyPixel.xₐ;
	MyPixel.y  = forward_model(MyPixel.x, MyPixel);
	MyPixel.iter_label = 0;
	
	# Step2: iteration
	try
		GN_Interation!(
			MyPixel, 
		    model=forward_model,
		    nIter=nIter,
		    thr_Converge=thr_Converge
		)
	catch e
		# println("Catch error: $e")
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

# ╔═╡ 9952e85a-62f0-4303-b239-49fbcb4db299
Retrieval3 = Retrieval.(
	eachslice(R_toa[SIF_index, :], dims=1),
	sza[SIF_index],
	vza[SIF_index],
	nflh[SIF_index], 
	Ref(params₃) 
)

# ╔═╡ 447c4c0e-d254-4af7-b330-d50833a0e00b
sum(ismissing.(Retrieval3))

# ╔═╡ 83399515-5cbf-4438-a957-5b27ac95dfee
function reconstruct3(
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	ρ     = hcat(v...)' * px.x[1 : px.nPoly+1];
	
	# T↑ transmittance for SIF
	T₁    = (px.trans_mat * px.x[(px.nPoly+2):(px.nPoly+px.nPC+1)]) .+ 1.0;

	# T↓↑ transmittance for SIF
	T₂    = @. exp( px.x[px.nPoly+px.nPC+2] * log(T₁) );

	# SIF magnitude
	SIF   = px.SIF_shape * px.x[px.nPoly+px.nPC+px.nSIF+2];
	
	return ρ, T₂, T₁, SIF
end

# ╔═╡ 85c189c1-9a74-4a6f-85d7-4634fccbfa7e
begin
	# make transmittance
	p_rho₃ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="two-way transmittance [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)
	
	p_trans₃₂ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="two-way transmittance [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)

	p_trans₃₁ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="one-way transmittance [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)

	p_SIF₃ = plot(
		size=(800, 400), 
		legendcolumns=4,
		xlabel="[nm]",
		ylabel="SIF",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		label=:outerbottom,
		title="ensemble of retrieval nFLH=[0.5, 0.7]"
	)
	
	for i in 1:(Δn):number_of_px
		if ismissing(Retrieval3[i])
			continue
		end
		# reconstruct
		ρⱼ, T₃₂ⱼ, T₃₁ⱼ, SIFⱼ = reconstruct3(Retrieval3[i]);
		plot!(
			p_rho₃, oci_band, ρⱼ,
			label="$(round(Retrieval3[i].flag, digits=2))")
		plot!(
			p_trans₃₂, oci_band, T₃₂ⱼ, 
			label="$(round(Retrieval3[i].flag, digits=2))")
		plot!(
			p_trans₃₁, oci_band, T₃₁ⱼ,
			label="$(round(Retrieval3[i].flag, digits=2))")
		plot!(
			p_SIF₃, oci_band, SIFⱼ .* T₃₁ⱼ,
			label="$(round(Retrieval3[i].flag, digits=2))",
			linestyle=:dash,
		)
		plot!(
			p_SIF₃, oci_band, SIFⱼ, label="ori. $(round(Retrieval3[i].flag, digits=2))")
	end
end

# ╔═╡ af488e29-a283-45f9-9938-c3e2b54db23d
p_trans₃₂

# ╔═╡ f09a2337-cf8b-4f2f-9e58-06a1eb312749
p_trans₃₁

# ╔═╡ 3155513a-c84d-43c7-99a9-881174be68e2
p_SIF₃

# ╔═╡ 496c7daf-2fc2-48ef-a88a-e4b8dc0460dc
# residual


# ╔═╡ ae8a34be-a19d-4780-86b7-546944e0f413
md"""
### Where are the nFLH fitting baseline wavelengths?
----
Should be baseline after atmospheric correction, can I apply it?
"""

# ╔═╡ f91b3e1e-e4b0-4556-b8ca-a261c8376c51
begin
	wlrad_bl_wvlen = [649.599976, 650.900024, 652.099976, 653.299988,
	                    654.599976, 655.799988, 657.099976, 658.299988,
	                    659.599976, 710.500000, 711.799988, 713.000000,
	                    714.299988, 716.799988, 719.200012];
	wl_bl_ind      = 
		map(wlrad_bl_wvlen -> argmin(abs.(oci_band .- wlrad_bl_wvlen)), wlrad_bl_wvlen);

	# mark on TOA radiance 
	plot(oci_band, R_px, size=(600, 200), label="TOA rad.")
	scatter!(oci_band[wl_bl_ind], R_px[wl_bl_ind], lw=1., label="nFLH baseline")
	scatter!(oci_band[bl_ind], R_px[bl_ind], lw=1., label="self-defined baseline (no abs.)")
end

# ╔═╡ 9eb3a674-f162-469f-901b-c306b0d9f831
plot(oci_band, MyPixel.E, size=(600, 200))

# ╔═╡ Cell order:
# ╟─bcf01e82-a621-11f0-24e9-db4a3bfd47b5
# ╠═ae9f792f-5175-4b47-babe-8ee7e50cebe0
# ╠═669f0127-0cfa-42e1-9325-4c309fb225a4
# ╠═ca59d107-6c89-46c2-bca3-2d3d604c6d27
# ╠═e19451d2-ec27-4b5e-a654-94abe2387e16
# ╠═7c1796ad-0654-4642-b615-035099743c27
# ╟─d3e33a05-871c-49a0-a91c-c5df2028d96f
# ╠═a9817961-e3ab-42f0-8b7c-580d4a39bbdd
# ╠═3060b59c-55b0-4e35-ac8c-5d7e2279fbd2
# ╟─22c13315-e686-4632-b24c-d6b35011f19d
# ╟─deb62e05-2d00-4704-8288-be7e292e5a5b
# ╠═e001c698-ddd3-4a61-9205-25c2e523e28d
# ╠═b6c3c640-e5fe-4322-8dc6-0e4a731704e7
# ╠═2bbbc4d9-c869-4dbf-8611-023d62ec24b0
# ╠═3ef82576-1699-4d31-8759-2462d2baf6de
# ╠═e1fe6027-a214-43f9-98cc-e8c431b1b818
# ╟─6f3d41e1-52cd-4a99-b907-47733be3e7b6
# ╠═5115015a-4479-48aa-80b4-509deb9d038e
# ╟─0e24b3a7-cebe-4a10-a37f-f92d37744508
# ╠═3bcdf18f-17a3-4664-b91e-0c6604825c5d
# ╠═4cb8bb3e-2b0a-4cb9-a37f-332f8c9da26d
# ╟─f6762889-f05d-4ff2-b77e-1eb3ac2d3a91
# ╠═8ce8f9ba-718c-4e60-aaa1-7aba7fe374f4
# ╠═caa101ef-a81c-4b43-96c0-a2a6810ea756
# ╠═c458c47e-e272-4ee9-8936-1ae11a93998c
# ╠═7f33d270-b1be-49bf-9c66-57598a891901
# ╠═e55eacaa-beb5-439f-bcf2-fc647f69058b
# ╠═8f13bb43-45f9-4e93-8fde-1c23b78964e6
# ╟─9e09c9bc-691e-4d53-ad59-1f6fb0a32fff
# ╠═65dc5736-7ffb-4d44-a921-3c541647f196
# ╠═38daa8c9-2bde-49a1-b446-5aae7a480e4c
# ╠═957559fd-4e88-468e-9b64-a9cf5964d5bf
# ╠═c84dc329-0e51-4eb1-9f58-6cae728dc230
# ╟─2c95ba4c-c0a5-491e-967b-d896814b91f7
# ╟─9c28e65b-9522-47b3-85c0-c7c47852cf72
# ╟─30b65205-de3b-4729-9287-40a9022bb9f2
# ╠═6eee7eae-dcc3-4eaf-902a-5436ec04677b
# ╠═6ca8f4f5-c063-457d-8dcb-6b7e66b547eb
# ╠═18505102-dd32-4e5a-bc2e-80f77a0e30be
# ╟─373cf733-983e-4d68-82f6-473295c081ba
# ╟─451294b0-3451-4034-bd86-22254bcbb561
# ╟─226ee786-1aad-4576-8797-0a340e2a5142
# ╟─269abad8-70df-433d-805b-64d762dacadc
# ╟─ef592e81-291c-4ffc-8f8a-b30f3a39dc1c
# ╠═f76d7c40-1ff4-4e39-bbae-cd7dd6834b7a
# ╠═2317ce66-50e6-4ae2-9796-27e0f1cd1534
# ╠═c2ae0783-0623-45f1-a009-c824148cd59d
# ╠═64eb2c06-d043-45dc-b60f-b44a3e777db8
# ╠═16af025b-e817-402e-8441-24e0e94761cf
# ╠═0fac32c1-35bf-4fc5-a62c-68d904bedc73
# ╠═43e06ecb-2e12-4698-bffc-9fd18dde13ea
# ╟─f27387bf-ceb1-4d57-b8f3-3597044962e3
# ╟─1ec5f75d-e017-436c-8f53-1ad9469fbd99
# ╟─873ebdde-8260-4fb2-8c02-88a0b722cf44
# ╟─53af5378-fed0-45ce-98ae-a2b5c32c0421
# ╟─906167fc-9516-499a-8fb5-229b774dc549
# ╟─ba6a2fc3-6ea6-41af-82f3-bddbfe3dd0ce
# ╟─16a28e7f-f3f1-4496-ad11-efc6af48992e
# ╟─cf673611-4324-464f-aec8-b56aaf97afa7
# ╠═3509ff7b-9a8d-413f-b5d3-96d83a936a32
# ╠═267ed111-0136-46af-af85-742f303b17bf
# ╠═798ad5dd-94e9-4f8d-a688-401a8b9e0d1a
# ╠═c8e554dc-3621-4a5d-9d75-394ee1c1a56b
# ╠═9952e85a-62f0-4303-b239-49fbcb4db299
# ╠═447c4c0e-d254-4af7-b330-d50833a0e00b
# ╠═83399515-5cbf-4438-a957-5b27ac95dfee
# ╠═85c189c1-9a74-4a6f-85d7-4634fccbfa7e
# ╟─af488e29-a283-45f9-9938-c3e2b54db23d
# ╟─f09a2337-cf8b-4f2f-9e58-06a1eb312749
# ╟─3155513a-c84d-43c7-99a9-881174be68e2
# ╠═496c7daf-2fc2-48ef-a88a-e4b8dc0460dc
# ╟─ae8a34be-a19d-4780-86b7-546944e0f413
# ╠═f91b3e1e-e4b0-4556-b8ca-a261c8376c51
# ╠═9eb3a674-f162-469f-901b-c306b0d9f831
