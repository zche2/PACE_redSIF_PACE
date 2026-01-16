### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ b7b035cc-f0c5-463c-90cd-9088a9e23c3a
begin
	import Pkg
	Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")
	using PACE_SIF
end

# ╔═╡ 7a953e1f-de4e-4789-8465-66fb230b1200
using Polynomials, LegendrePolynomials, ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics, Parameters

# ╔═╡ 549101b7-2dc3-4ac1-9658-52e7a8c159d3
md"""
> ##### Prepare data
- defining fitting window
- transmittance spectra
- OC1 L1B data + L2 AOP nFLH
- PACE OCI SNR
"""

# ╔═╡ 843bb5a8-20f2-439a-a6b5-f9b0424f5254
λ_min = 668.2; λ_max = 780.;

# ╔═╡ 31ed06c7-0eee-49e8-b1af-1f3de3cefb99
begin
	summer = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc");
	winter = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc");
	println("Opened datasets.")
	
	trans = cat(summer["transmittance"][:,:], winter["transmittance"][:,:], dims=1);
	println("\nConcatenated!")

	bands  = summer["band"][:];
	
	close(summer);
	close(winter);
	
	# SVD
	HighResSVD = Spectral_SVD(trans, bands, λ_min=λ_min, λ_max=λ_max);
end

# ╔═╡ e0c33dc7-c29b-44df-bae4-d297cf7b0bbe
begin
	oci = Dataset(
		"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/sample_granule_20250501T183011_new_chl.nc");
	red_band = oci["red_wavelength"][:];
	nflh     = oci["nflh"][:, :];
	vza      = oci["sensor_zenith"][:, :];
	sza      = oci["solar_zenith"][:, :];
	println("Read in Dataset")

	# select band (continuum spectrum)
	ind      = findall( λ_min .< red_band .< λ_max );
	E        = oci["red_solar_irradiance"][ind];
	R_toa    = oci["radiance_red"][:, :, ind];
	oci_band = red_band[ind];
end

# ╔═╡ f6876469-d979-4a66-bac2-657d74df08d5
begin
	filename = raw"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt";
	lines = readlines(filename);
	end_header_index = findfirst(x -> x == "/end_header", lines);
	data = readdlm(filename, String, skipstart=end_header_index);

	FPA   = data[:, 1];                   # 1st column: band
	wvlen = parse.(Float64, data[:, 2]);  # 2nd column: center wavelength
	c1    = parse.(Float64, data[:, 4]);  # 4th column: c1
	c2    = parse.(Float64, data[:, 5]);  # 5th column: c2

	wv_val  = (λ_min .< wvlen .< λ_max);
	snr_ind = findall((FPA .== "Red") .& wv_val);
end

# ╔═╡ 39828d72-132b-4049-8911-4aca71cc157d
findall(coalesce.((nflh .> 4.2) .& (nflh .< 4.6), false))

# ╔═╡ 59eb9f37-90a2-4e2a-a8f3-7bbe8afd03ec
md"""
> ##### Toolkits
"""

# ╔═╡ 59fb46d5-67ee-4744-96e6-febd4f7e06c2
function center_wavelength(λ)
	# get the range and medium of λ and center it to [0,1]
	λ_max = ceil(maximum(λ));
	λ_min = floor(minimum(λ));
	range = (λ_max - λ_min) / 2;
	λ_middle = (λ_max + λ_min) / 2;
	# λ_median = median(λ);
	λc    = (λ .- λ_middle) ./ range;
	return λc
end

# ╔═╡ 6b341b63-fcef-4cda-a014-f359c31e1637
function find_baseline_band(
		λ;
		bl_wvlen = [668.265, 669.518, 670.755, 671.99, 673.245, 674.505, 675.73, 676.962, 678.205, 679.445, 680.68, 751.79, 753.04, 754.295, 776.832, 779.335, 867.115, 869.615, 872.13]
	)
	λ_bl_ind = map(bl_wvlen -> argmin(abs.(λ .- bl_wvlen)), bl_wvlen);
	return λ_bl_ind
end

# ╔═╡ 66f97376-3275-4b85-b0ea-3a1775169e33
function two_way_trans(T, sza, vza)
	svf = (secd(sza)+secd(vza)) / secd(vza);
	T2  = exp.( svf .* log.(T));
	return T2
end

# ╔═╡ b317ae66-d58e-409c-8a47-73e8fea4d579
function SIF_shape(λ; λ₀=683., σ=5.)
	return exp.( - ( λ .- λ₀ ).^2 ./ ( 2 * σ^2 ) )
end

# ╔═╡ 81e8620d-8895-4d98-8dee-f5b45c3b2663
function scale_transmittance(T, λ_bl_ind)
	# find max
	T_abs = abs.(T)
	bl_max = maximum(T_abs[λ_bl_ind]);
	# force the mean val to be 1
	T_norm = T_abs ./ bl_max
	return T_norm
end

# ╔═╡ 308ec29c-dfee-4eab-bee8-daceb2d778f3
function Jacobian(x, model; len=length(oci_band))
	res = DiffResults.JacobianResult(zeros(len), x);
	ForwardDiff.jacobian!(res, model, x);
	K   = DiffResults.jacobian(res);
	val = DiffResults.value(res);
	return K, val
end

# ╔═╡ 54d08486-e9b7-4c4b-9559-eacbfa85dd85
function GainMatrix(K; Se=Se, Sa=Sa)
	return inv( K'inv(Se)K + inv(Sa) )K'inv(Se)
end

# ╔═╡ 85a6ce19-a80d-49e9-80da-febf82949b0d
function root_mean_square(y_obs, y_retrieve; n=length(oci_band))
	totSQ = sum( ( y_obs .- y_retrieve ) .^ 2 );
	RMS   = sqrt( totSQ / n );
	return RMS
end

# ╔═╡ e4a2cd9d-1d1f-4df1-8b18-d2f51e98a62d
function construct_Se(c1, c2, R_toa)
	noise   = sqrt.( c1 .+ c2 .* R_toa);
	Se      = Diagonal(noise.^2);
	return Se
end

# ╔═╡ 830f3ee5-828d-44ae-af5a-4a0ff6dfe5e8
function construct_Sa(nPoly, nPC, nSIF; rel_error=.001, VarExp=HighResSVD.VarExp)
	Sa  = zeros(nPoly+nPC+nSIF+1, nPoly+nPC+nSIF+1);
	# Poly
	for i=1:(nPoly+1)
	    Sa[i,i] = 1e20;     
	end
	# PC terms
	for i=(nPoly+2):(nPoly+nPC+1)
	    Sa[i,i] = rel_error .* VarExp[i - (nPoly+1)];
	end
	# SIF terms
	Sa[end-nSIF+1, end-nSIF+1] = 1;
	return Sa
end

# ╔═╡ 25376370-bc83-44d8-a6a2-4dc7db301038
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

# ╔═╡ e19a673b-a21d-4090-a825-fac64bd63483
function forward_model(
		x,
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	rho   = hcat(v...)' * x[1 : px.nPoly+1];
	
	# transmittance
	T      = px.trans_mat * x[(px.nPoly+2):(px.nPoly+px.nPC+1)];
	T_norm = scale_transmittance(T, px.λ_bl_ind);

	# one way vs. two way
	T2_norm = two_way_trans(T_norm, px.sza, px.vza);
	
	# SIF magnitude
	SIF    = @. x[px.nPoly+px.nPC+px.nSIF+1] * SIF_shape(px.λ);
	# TOA radiance
	rad    = @. px.E * cosd(px.sza) / pi * T2_norm * rho + SIF * T_norm;
	return rad
end

# ╔═╡ 8c3d6809-f3d3-4468-9501-f0e9aa482005
function reconstruct(
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	rho   = hcat(v...)' * px.x[1 : px.nPoly+1];
	
	# transmittance
	T      = px.trans_mat * px.x[(px.nPoly+2):(px.nPoly+px.nPC+1)];
	T_norm = scale_transmittance(T, px.λ_bl_ind);

	# one way vs. two way
	T2_norm = two_way_trans(T_norm, px.sza, px.vza);
	
	# SIF magnitude
	SIF    = @. px.x[px.nPoly+px.nPC+px.nSIF+1] * SIF_shape(px.λ);
	
	return rho, T_norm, T2_norm, SIF
end

# ╔═╡ 6eaf833b-336e-4405-b0e3-b6809fdf0a2a
function GN_Interation!(
			px :: Pixel;
			model = forward_model,
			nIter = 20,
			thr_Converge = 1e-8,
		)
	
	# initial
	xₐ = px.xₐ;   # priori estimation
	xₙ = px.x;
	Kₙ, _ = Jacobian(xₙ, x -> model(x, px));
	# k     = px.iter_label;       # number of iterations
	RMSE₀ = 1e20; 
	RMSE₁ = root_mean_square(px.R_toa, px.y);
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
		RMSE₁  = root_mean_square(px.R_toa, px.y);
		ΔRMSE = RMSE₁ - RMSE₀;
		px.ΔRMSE = ΔRMSE;
	end
	
	return nothing
end

# ╔═╡ 07136a7d-2786-4d49-b81f-231c2fe9e592
begin
	pixel = 504; scan = 820;
	@show "lat: ", oci["latitude"][pixel, scan];
	@show "lon: ", oci["longitude"][pixel, scan];
	
	MyPixel   = Pixel();

	# Step1: construct struct
	MyPixel.λ  = oci_band;
	MyPixel.λc = center_wavelength(MyPixel.λ);
	MyPixel.λ_bl_ind = find_baseline_band(MyPixel.λ);
	MyPixel.E     = E;
	MyPixel.nPoly = 3;
	MyPixel.nPC   = 25;
	MyPixel.nSIF  = 1;
	MyPixel.Sa    = construct_Sa(MyPixel.nPoly, MyPixel.nPC, MyPixel.nSIF);
	MyPixel.trans_mat = HighResSVD.PrinComp[:, 1:MyPixel.nPC];

	MyPixel.R_toa = R_toa[pixel, scan, :];
	MyPixel.sza   = sza[pixel, scan];
	MyPixel.vza   = vza[pixel, scan];
	MyPixel.Se    = construct_Se(c1[snr_ind], c2[snr_ind], MyPixel.R_toa);
	MyPixel.flag  = nflh[pixel, scan];   # quite temporary..
	
	# a priori estimation
	K₀  = 
		MyPixel.E .* cosd(MyPixel.sza) ./ pi .* hcat(collectPl.(MyPixel.λc, lmax=MyPixel.nPoly)...)';
	G₀  = inv( K₀'K₀ )K₀';
	x₀  = G₀ * MyPixel.R_toa;
	tmp = zeros(MyPixel.nPC + MyPixel.nSIF - 2) .+ .001;
	MyPixel.xₐ = [x₀... -10. 0.1 tmp...]';

	# set-up
	MyPixel.x  = [x₀... -10. 0.1 tmp...]';
	MyPixel.y  = forward_model(MyPixel.x, MyPixel);
	MyPixel.iter_label = 0;

	# Step2: iteration
	GN_Interation!(MyPixel)
end

# ╔═╡ cd64f729-1462-4044-ae65-7a73e4c3c206
# reconstruct
rho, T1, T2, SIF_px = reconstruct(MyPixel);

# ╔═╡ 24b2a34e-70df-448d-a9e8-19857e20ccf7
begin
	# averaging kernel
	# Jacobian
	K̂, _ = Jacobian(MyPixel.x, x -> forward_model(x, MyPixel));
	Ĝ    = GainMatrix(K̂, Se=MyPixel.Se, Sa=MyPixel.Sa);
	Â    = Ĝ * K̂;
	@show tr(Â)
	heatmap(Â,
		c=:greys,      # <--- Change the colormap here
		clims=(0, 1),  # <--- vmin and vmax
		# aspect_ratio=:equal
		size=(450, 420),
		dpi=400
	)
	xlabel!("Truth")
	ylabel!("Retrieval")
	# savefig("Averaging_kernel_2.pdf")
end

# ╔═╡ 92ff1fa5-daa6-48d6-b4e8-67b382d0dc68
md"""
> ##### Pixel-by-pixel retrieval
"""

# ╔═╡ 03630196-42dd-4b8c-976c-9e020c7306b5
function SIF_Retrieval(
		# "L1B pixel-by-pixel vals"
		R_toa ,
		sza, vza, flag,   # quite temporary..
		params,
		# # "Window and irradiance"
		# λ, E, 
		# # "forward model"
		# forward_model,
		# # "retrieval params"
		# nPoly, nPC, nSIF, nIter, thr_Converge
		# c1, c2
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
	MyPixel.λc = center_wavelength(MyPixel.λ);
	MyPixel.λ_bl_ind = find_baseline_band(MyPixel.λ);
	MyPixel.E     = params.E;
	MyPixel.nPoly = params.nPoly;
	MyPixel.nPC   = params.nPC;
	MyPixel.nSIF  = params.nSIF;
	MyPixel.Sa    = construct_Sa(MyPixel.nPoly, MyPixel.nPC, MyPixel.nSIF);
	MyPixel.trans_mat = HighResSVD.PrinComp[:, 1:MyPixel.nPC];

	MyPixel.R_toa = R_toa;
	MyPixel.sza   = sza;
	MyPixel.vza   = vza;
	MyPixel.Se    = construct_Se(c1[snr_ind], c2[snr_ind], MyPixel.R_toa);
	MyPixel.flag  = flag; 
	# a priori estimation
	K₀  = 
		MyPixel.E .* cosd(MyPixel.sza) ./ pi .* hcat(collectPl.(MyPixel.λc, lmax=MyPixel.nPoly)...)';
	G₀  = inv( K₀'K₀ )K₀';
	x₀  = G₀ * MyPixel.R_toa;
	tmp = zeros(MyPixel.nPC + MyPixel.nSIF - 2) .+ .001;
	MyPixel.xₐ = [x₀... -10. 0.1 tmp...]';

	# set-up
	MyPixel.x  = [x₀... -10. 0.1 tmp...]';
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
		return MyPixel.x[end-MyPixel.nSIF+1]
	else
		return missing
	end
end

# ╔═╡ 2aa42e02-59fa-4fb6-b893-a64c14725b04
begin
	# fix the params
	params = (
		λ = oci_band, 
		E = E,
		forward_model = forward_model,
		nPoly = 3,
		nPC   = 25, 
		nSIF  = 1,
		nIter = 20, 
		thr_Converge = 1e-8
	)
	
end

# ╔═╡ 7635c5d4-68c6-4b40-9f8f-8c739c1ba284
# try single pixel
begin
	@show SIF_Retrieval(
	    R_toa[pixel, scan, :],
		sza[pixel, scan],
		vza[pixel, scan],
		nflh[pixel, scan], 
	    params
	)
end

# ╔═╡ 9399daaf-7526-4d37-9b53-16fd90f231a7
SIF_683 = SIF_Retrieval.(
    eachslice(R_toa, dims=(1,2)),
	sza,
	vza,
	nflh, 
    Ref(params) 
)


# ╔═╡ 46735ae0-c9b9-4609-ad14-6a9f3e6c9a53
rho_red  = mean(oci["rhot_red"][:, :, ind], dims=3);

# ╔═╡ 3d3566e8-1b57-4f91-b5a3-f5a5f532ee87
chlor_a  = coalesce.(log.(oci["chlor_a"]), -6);

# ╔═╡ 7498f09a-37b1-4361-a87e-ca1b0c9021d5
begin
	@show chlor_a[pixel, scan]
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
		);
	plot!(
		MyPixel.λ, forward_model(MyPixel.xₐ, MyPixel), label="initial guess"); 
	plot!(
		MyPixel.λ, MyPixel.y, label="@ convergence", size=(600, 200), linewidth=1);
	
end

# ╔═╡ cb511427-f5dc-4157-9b50-70d5fdd4ba45
begin
	rho_fig   = plot(MyPixel.λ, rho, label="surface reflectance",
					title=TheTitle,
					titlefontsize=10,
					)
	trans_fig = plot(MyPixel.λ, T1, label="T↑")
	plot!(trans_fig,MyPixel.λ,  T2, label="T↓↑")
	SIF_fig   = plot(MyPixel.λ, SIF_px,
					label="SIF₆₈₃=$(round(MyPixel.x[end], digits=3))",
					xlabel="Wavelength [nm]",
					xlabelfontsize=10,
					)
	plot(rho_fig, trans_fig, SIF_fig,
		 layout=(3,1),
		 size=(600, 450)
	)
end

# ╔═╡ e36b9dda-bf56-4af3-ba0f-b9e49730a959
begin	
	r1 = plot(
		MyPixel.λ, MyPixel.R_toa .- MyPixel.y, label="Residual (W/m²/µm/sr)", linewidth=1.5, color=:grey);
	title!(TheTitle, titlefontsize=10)
	r2 = plot(MyPixel.λ, (MyPixel.R_toa .- MyPixel.y)./MyPixel.R_toa * 100, 		   label="Relative Residual (%)", linewidth=1.5,
		color=:grey);
	plot(r1, r2, layout=(2,1), size=(600, 400))
end

# ╔═╡ 9fdf7797-fcc2-4bf1-b209-df6921b15e55
begin
	scatter(
		vec(SIF_683), vec(nflh),
		markerstrokewidth=0,  
		markersize=3,    
		markeralpha=.3,
		marker_z=vec(chlor_a),
		xlabel="SIF₆₈₃ (W/m²/µm/sr)",
	    ylabel="nFLH (W/m²/µm/sr)",
		label = "",
		colorbar_title="Chl concentration (mg/m³)",
		clims = (-1, 2.5),
		size  = (500, 400),
		aspect_ratio=:equal,
		xticks=-.4:1.:5.6,             
	    yticks=-.4:1.:5.6,
		xlim=(-0.4, 2.6),
		ylim=(-0.4, 2.6),
		dpi=400,
	    # title="20250501T183011 L1B",
	)
	plot!(-.3:4.6, -.3:4.6,
		color=:silver,
		linestyle=:dot,
		linewidth=2,
		label="1:1 Line"
	)
	# savefig("20250501T183011_SIF683_20250904.png")
end

# ╔═╡ Cell order:
# ╠═b7b035cc-f0c5-463c-90cd-9088a9e23c3a
# ╠═7a953e1f-de4e-4789-8465-66fb230b1200
# ╟─549101b7-2dc3-4ac1-9658-52e7a8c159d3
# ╠═843bb5a8-20f2-439a-a6b5-f9b0424f5254
# ╟─31ed06c7-0eee-49e8-b1af-1f3de3cefb99
# ╠═e0c33dc7-c29b-44df-bae4-d297cf7b0bbe
# ╠═f6876469-d979-4a66-bac2-657d74df08d5
# ╠═39828d72-132b-4049-8911-4aca71cc157d
# ╟─59eb9f37-90a2-4e2a-a8f3-7bbe8afd03ec
# ╟─59fb46d5-67ee-4744-96e6-febd4f7e06c2
# ╠═6b341b63-fcef-4cda-a014-f359c31e1637
# ╟─66f97376-3275-4b85-b0ea-3a1775169e33
# ╟─b317ae66-d58e-409c-8a47-73e8fea4d579
# ╟─81e8620d-8895-4d98-8dee-f5b45c3b2663
# ╠═308ec29c-dfee-4eab-bee8-daceb2d778f3
# ╠═54d08486-e9b7-4c4b-9559-eacbfa85dd85
# ╟─85a6ce19-a80d-49e9-80da-febf82949b0d
# ╟─e4a2cd9d-1d1f-4df1-8b18-d2f51e98a62d
# ╟─830f3ee5-828d-44ae-af5a-4a0ff6dfe5e8
# ╠═25376370-bc83-44d8-a6a2-4dc7db301038
# ╠═e19a673b-a21d-4090-a825-fac64bd63483
# ╠═8c3d6809-f3d3-4468-9501-f0e9aa482005
# ╠═6eaf833b-336e-4405-b0e3-b6809fdf0a2a
# ╠═07136a7d-2786-4d49-b81f-231c2fe9e592
# ╠═7498f09a-37b1-4361-a87e-ca1b0c9021d5
# ╠═cd64f729-1462-4044-ae65-7a73e4c3c206
# ╟─cb511427-f5dc-4157-9b50-70d5fdd4ba45
# ╟─e36b9dda-bf56-4af3-ba0f-b9e49730a959
# ╠═24b2a34e-70df-448d-a9e8-19857e20ccf7
# ╟─92ff1fa5-daa6-48d6-b4e8-67b382d0dc68
# ╠═03630196-42dd-4b8c-976c-9e020c7306b5
# ╠═2aa42e02-59fa-4fb6-b893-a64c14725b04
# ╠═7635c5d4-68c6-4b40-9f8f-8c739c1ba284
# ╠═9399daaf-7526-4d37-9b53-16fd90f231a7
# ╠═46735ae0-c9b9-4609-ad14-6a9f3e6c9a53
# ╠═3d3566e8-1b57-4f91-b5a3-f5a5f532ee87
# ╠═9fdf7797-fcc2-4bf1-b209-df6921b15e55
