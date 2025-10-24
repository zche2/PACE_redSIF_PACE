### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ bb10dd17-e8be-4bf1-9caa-c2f645d07046
begin
	import Pkg; 
	Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE");
	# Pkg.develop(path="/home/zhe2/FraLab/PACE_redSIF_PACE")
end

# ╔═╡ 489b1161-e108-4506-915b-07143c546d70
using JLD2, Interpolations, Revise

# ╔═╡ d8546b0b-d901-4520-96b9-e3601f0ab70f
using Polynomials, ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics

# ╔═╡ 3b32eb28-b9ee-4103-bdaf-f65e0cd89760
using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools

# ╔═╡ 293e259b-af97-476b-a8be-103c731e6d1a
using PACE_SIF

# ╔═╡ 64ad88f8-af61-11f0-230d-a525e3af3a49
md"""
> ## Forward model v4
- Apply NMF rather than SVD to transmittance spectra
🟡 To be updated: the second PC is almost fixed (should not be a free parameter).
"""

# ╔═╡ 584e1386-e522-4774-adee-f0d54f55db78
md"""
### Spin up
---
"""

# ╔═╡ 5afd8471-0ec5-454a-aeb3-1921a60bd48f
# wavelenth
λ_min = 610.; λ_max = 820.;

# ╔═╡ 8c4f5fa0-0929-4a55-8e52-bbcd4adfcbb1
# wavenumber
n₁ = 1e7 / λ_min; n₂ = 1e7 / λ_max; println("wavenumber is $n₂ to $n₁")

# ╔═╡ 3e56380e-6ce7-4c34-a9c5-fa75a06a6b6a
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
	println("\nBand selected: $oci_band")
end

# ╔═╡ 41ee0706-efe3-4c5a-859e-dd8f0390732b
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

# ╔═╡ 5e3d0da7-1910-4879-b2e0-f4a1ea32f5e2
begin
	# load SIF
	SIF_shape_dict = JLD2.load("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/SIF_singular_vector.jld2")
	println("SIF data loaded")
end

# ╔═╡ 9cb2e8d3-2fdb-4a6d-8f81-1e7826e6cc8e
begin
	# interpolation in the first dimension and no interp. in the second
	itp    = interpolate(SIF_shape_dict["SIF_U"], (BSpline(Linear()), NoInterp()));
	
	# scale
	range₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end];
	range₂ = 1:size(itp, 2);
	sitp   = scale(itp, range₁, range₂);

	# set extrapolation filling value = 0
	setp0  = extrapolate(sitp, 0)

	# interpolation
	SIF_new = reduce(hcat, [setp0.(oci_band, i) for i in range₂]); 
	
	println("SIF shape interpolated")
end

# ╔═╡ b8ba6ff5-0d91-4201-b38e-8afd6c60fdfa
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

# ╔═╡ 7c4fcde6-0e72-4229-a60c-969eda6050b1
md"""
### NMF
---
"""

# ╔═╡ 41cc9e03-be5d-49f9-8769-2e1187ee900c
begin
	rank       = 10;
	# NMF
	HighResNMF = Spectral_NMF(
		trans, 
		bands,
		Float64.(collect(skipmissing(oci_band))); 
		rank=rank
	);
end

# ╔═╡ cf4159bd-23b1-4ffd-be6b-02ec2aa6e34f
begin
	# W and H
	λ₀ = HighResNMF.band;
	W₀ = HighResNMF.Loading;
	H₀ = HighResNMF.PrinComp;
	
	# matrics
    mean_val  = [round(mean(W₀[:, i]), digits=2) for i in 1:rank];
    max_val   = [round(maximum(W₀[:, i]), digits=2) for i in 1:rank];
    min_val   = [round(minimum(W₀[:, i]), digits=2) for i in 1:rank];
    mean_spec = mean(trans[:, ind], dims=1);

    # Create a plot with k panels (one for each row)
    plot(
        [
			begin
            p = plot(oci_band, H₀[i, :], label="",title="$(mean_val[i]) ($(min_val[i]), $(max_val[i]))", lw=2.)
            plot!(p, oci_band, vec(mean_spec), color=:silver, label="", lw=2., alpha=.3)
            p
        end for i in 1:rank]..., 
        layout=(rank÷2, 2), 
        size=(1200, 800),  
    )
end

# ╔═╡ 31069e84-bde9-4377-a75e-2cc4e9c9fc49
md"""
### Adapted from v3 - Forward model 4
---
Joint fit of $T_{\downarrow\uparrow}$ and $T_{\uparrow}$ by tuning a factor: 

$$\rho_{s}(\lambda)=\sum{a_jP_j}, \ T_{\uparrow}(\lambda)=\sum{\beta_i P_i}$$

$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)T_{\downarrow\uparrow}(\lambda)}{\pi} + SIF(\lambda)T_{\uparrow}(\lambda)$$

$$T_{\downarrow\uparrow}(\lambda)=exp(\gamma \times ln(T_{\uparrow}(\lambda)))$$
Where γ is correction factor accounting for （1) light path (VZA and SZA) and （2) upper atmosphere reflectance.
"""

# ╔═╡ 08343917-236b-4bfb-be42-9b560139906c
begin
	n    = 2;
	nPC  = rank;
	nSIF = 1;

	# s.d. for the loading term
	loading_sd  = [std(W₀[:, i]) for i in 1:rank];
	loading_ave = [mean(W₀[:, i]) for i in 1:rank];
	
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

	# force Sₐ[nPoly+3, nPoly+3]=1e-10 (almost certain?)
	Sₐ[n+3, n+3] = 1e-20;
	
	# \gamma
	Sₐ[n+nPC+2, n+nPC+2] = 2;
	# SIF magnitude
	Sₐ[end, end] = 1;
	println("Diagonal terms are: $(diag(Sₐ))")

end

# ╔═╡ 403f6061-e5ca-40ba-9a64-afc0fe018bba
# center wavelength `oci_band`
λc = center_wavelength(oci_band);

# ╔═╡ c437a631-ab57-4b27-9550-8aa9c7e4b27e
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
	"normalized fluorescence height (nFLH)"
	nflh
	"chlor_a concentration"
	chlor_a
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

# ╔═╡ 4fc5d644-d346-42d8-88e0-37ac0cd46ecb
function Jacobian(x, model; len=length(oci_band))
	res = DiffResults.JacobianResult(zeros(len), x);
	ForwardDiff.jacobian!(res, model, x);
	K   = DiffResults.jacobian(res);
	val = DiffResults.value(res);
	return K, val
end

# ╔═╡ 4ac45a58-1988-43e6-93af-8a27a3d80a72
function GainMatrix(K; Se=Se, Sa=Sa)
	return inv( K'inv(Se)K + inv(Sa) )K'inv(Se)
end

# ╔═╡ e1ade8ed-21b6-467d-a291-02ae3fe84d5e
function GN_Interation!(
			px :: Pixel,
			model;
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
		RMSE₁  = root_mean_square(px.R_toa, px.y);
		ΔRMSE  = RMSE₁ - RMSE₀;
		px.ΔRMSE = ΔRMSE;
	end
	
	return nothing
end

# ╔═╡ f0b8aed0-1455-4cd6-b42c-697287943e59
begin
	# sigmoid apporximation function, bounded (for now) by [1,2]
	sigm(x) = 1. / (1 + exp(-x)) + 1.;
end

# ╔═╡ cb604d64-a31f-4f97-bf48-a816d8abf83b
function MakePriori!(
		px :: Pixel,
		β :: Vector{FT}
	) where {FT <: AbstractFloat}
	# polynomial terms
	K₀  = px.E .* cosd(px.sza) ./ pi .* hcat(collectPl.(px.λc, lmax=px.nPoly)...)';
	G₀  = inv( K₀'K₀ )K₀';
	x₀  = G₀ * px.R_toa; 
	# gamma
	γ   = (secd(px.sza) + secd(px.vza)) / secd(px.vza);
	# SIF
	SIF   = px.nflh       # nflh as an approximation
	px.xₐ = [x₀... β... γ SIF]';

	return nothing
end

# ╔═╡ 509c5b50-098e-4d66-878a-8cbc0de50ee9
function forward_model4(
		x,
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	ρ     = hcat(v...)' * x[1 : px.nPoly+1];

	# T↑ transmittance for SIF
	T₁    = (px.trans_mat * x[(px.nPoly+2):(px.nPoly+px.nPC+1)]);

	# T↓↑ transmittance for SIF
	smooth_x = sigm(x[px.nPoly+px.nPC+2]);
	T₂       = @. exp( smooth_x * log(T₁) );

	# SIF magnitude
	SIF   = px.SIF_shape * x[px.nPoly+px.nPC+px.nSIF+2];
	
	# TOA radiance
	rad   = @. px.E * cosd(px.sza) / π * T₂ * ρ + SIF * T₁;
	return rad
end

# ╔═╡ 1fc5b94b-6246-4cd0-881b-b010cef82a0c
function Retrieval4(
		# "L1B pixel-by-pixel vals"
		R_toa, sza, vza, nflh, chlor_a, flag,
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
	βₐ            = params.βₐ

	# Step1: construct struct
	MyPixel.λ  = params.λ;
	MyPixel.λc = params.λc;
	MyPixel.λ_bl_ind = params.λ_bl_ind;
	MyPixel.E     = params.E;
	MyPixel.nPoly = params.nPoly;
	MyPixel.nPC   = params.nPC;
	MyPixel.nSIF  = params.nSIF;
	MyPixel.Sa    = params.Sₐ;
	MyPixel.trans_mat = params.PrinComp[:, 1:MyPixel.nPC];
	MyPixel.SIF_shape = SIF_new[:, 1:MyPixel.nSIF];

	MyPixel.R_toa = R_toa;
	MyPixel.sza   = sza;
	MyPixel.vza   = vza;
	MyPixel.nflh  = nflh;
	MyPixel.chlor_a = chlor_a;
	noise         = sqrt.( c1[snr_ind] .+ c2[snr_ind] .* MyPixel.R_toa);
	MyPixel.Se    = Diagonal(noise.^2);
	MyPixel.flag  = flag; 
	
	# set-up
	MakePriori!(MyPixel, βₐ);
	MyPixel.x  = MyPixel.xₐ;
	MyPixel.y  = forward_model(MyPixel.x, MyPixel);
	MyPixel.iter_label = 0;
	
	# Step2: iteration
	try
		GN_Interation!(
			MyPixel, 
			forward_model,
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

# ╔═╡ 60c4eacc-6db1-43ad-945c-6af2147101f5
md"""
#### Select pixels & Retrieval
---
"""

# ╔═╡ f6155dde-0a72-4750-a82f-b045b8ada62b
begin
	nFLH_min = 0.02;
	nFLH_max = 0.8;
	SIF_index_all = findall(
		coalesce.((nflh .> nFLH_min) .& (nflh .< nFLH_max), false)
	);
	SIF_index     = SIF_index_all[1:20:end];
	println("# of pixels included: $(length(SIF_index))")
end

# ╔═╡ 8f3b7a74-fc81-47ba-aac0-c6179fa02e4d
begin
	# define new parameter set
	params₄ = (
		λ = oci_band, 
		E = E,
		λ_bl_ind = missing,
		λc    = λc,
		nPoly = n,
		nPC   = nPC, 
		nSIF  = 1,
		nIter = 20, 
		Sₐ    = Sₐ,
		βₐ    = loading_ave,
		PrinComp      = HighResNMF.PrinComp',
		thr_Converge  = 1e-5,
		forward_model = forward_model4,
	);
end

# ╔═╡ 06a27f81-3d65-4e39-a627-4f21b1c9830a
MyRetrieval = Retrieval4.(
	eachslice(R_toa[SIF_index, :], dims=1),
	sza[SIF_index],
	vza[SIF_index],
	nflh[SIF_index],      # nflh
	chlor_a[SIF_index],   # chlor_a
	nflh[SIF_index],      # flag
	Ref(params₄) 
)

# ╔═╡ 73286e0b-0356-4a79-8ef0-8e6d05d6cfe7
sum(ismissing.(MyRetrieval))

# ╔═╡ 09b0e9ae-123a-4fb9-a1cc-4efa9d6a07e4
md"""
#### Reconstruct
---
"""

# ╔═╡ c18b7fb5-c0c4-4c4e-bc69-2fe5fc0442ab
function reconstruct4(
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	ρ     = hcat(v...)' * px.x[1 : px.nPoly+1];

	# T↑ transmittance for SIF
	T₁    = (px.trans_mat * px.x[(px.nPoly+2):(px.nPoly+px.nPC+1)]);

	# T↓↑ transmittance for SIF
	smooth_x = sigm(px.x[px.nPoly+px.nPC+2]);
	T₂       = @. exp( smooth_x * log(T₁) );

	# SIF magnitude
	SIF   = px.SIF_shape * px.x[px.nPoly+px.nPC+px.nSIF+2];
	
	return ρ, T₂, T₁, SIF
end

# ╔═╡ 68ee1b3b-e624-431c-8908-142a6da86074
begin
	number_of_px = size(SIF_index)[1];
	println("make plots for `MyRetrieval`.")
	
	# make transmittance
	p_rho₄ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="ρ [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
	)
	
	p_trans₄₂ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="two-way transmittance [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
	)

	p_trans₄₁ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="one-way transmittance [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
	)

	p_SIF₄ = plot(
		size=(800, 400), 
		legendcolumns=4,
		xlabel="[nm]",
		ylabel="SIF",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		label=:outerbottom,
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
	)

	# null array to store that
	resd₄      = [];
	nflh_px    = [];
	SIF₆₈₃_px₄ = [];
	chl_px     = [];

	# some wavelength range
	nflh_bl_ind = argmin(abs.(oci_band .- 678.2));
	resd_bl_ind = findall(λ_min .< oci_band .< λ_max);  # 650.
	
	for i in 1:number_of_px
		if ismissing(MyRetrieval[i])
			push!(resd₄, missing);
			push!(nflh_px, missing);
			push!(SIF₆₈₃_px₄, missing);
			push!(chl_px, missing);
			continue
		end
		# reconstruct
		ρ₄ⱼ, T₄₂ⱼ, T₄₁ⱼ, SIF₄ⱼ = reconstruct4(MyRetrieval[i]);

		# push
		resdᵢ = MyRetrieval[i].y .- MyRetrieval[i].R_toa;
		nflhᵢ = MyRetrieval[i].nflh;
		ā     = norm(resdᵢ[resd_bl_ind], 2);
		push!(resd₄, ā);
		push!(nflh_px, nflhᵢ);
		push!(SIF₆₈₃_px₄, SIF₄ⱼ[nflh_bl_ind]);
		push!(chl_px, MyRetrieval[i].chlor_a);
	
		if i%100==0
			# plot
			plot!(
				p_rho₄, oci_band, ρ₄ⱼ,
				label="$(round(MyRetrieval[i].flag, digits=2))")
			plot!(
				p_trans₄₂, oci_band, T₄₂ⱼ, 
				label="")
			plot!(
				p_trans₄₁, oci_band, T₄₁ⱼ,
				label="$(round(sigm(MyRetrieval[i].x[MyRetrieval[i].nPoly+MyRetrieval[i].nPC+2]), digits=2))"*" $( round((secd(MyRetrieval[i].sza) + secd(MyRetrieval[i].vza)) / secd(MyRetrieval[i].vza) , digits=2))")
			plot!(
				p_SIF₄, oci_band, SIF₄ⱼ .* T₄₁ⱼ,
				label="",
				linestyle=:dash)
			plot!(
				p_SIF₄, oci_band, SIF₄ⱼ, label="")
		end
	
	end
end

# ╔═╡ efe7e97e-a1b8-4c05-bd9d-5222f4678d5c
p_trans₄₂

# ╔═╡ e6f7310a-0e11-4e67-933a-f32292f45047
p_trans₄₁

# ╔═╡ b678620c-4a20-43f9-b22e-7ffd199bde77
p_SIF₄

# ╔═╡ 43ad9bd4-dafd-47a1-99c5-d078e3381c4f
begin
	# residual
	p_resd₄ = plot(
		size=(800, 400), 
		legendcolumns=4,
		xlabel="[nm]",
		ylabel="Residual",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
	)

	for i in 1:50:number_of_px
		if ismissing(MyRetrieval[i])
			continue
		end
		# get spectral-wise residual
		resdᵢ = MyRetrieval[i].y .- MyRetrieval[i].R_toa;
		
		plot!(
			p_resd₄, oci_band, resdᵢ, label="",
		)
	end
	p_resd₄
end

# ╔═╡ d297a09a-5b97-43ef-ae01-367444ad87fe
begin
	p_hist = histogram2d(
		   SIF₆₈₃_px₄, nflh_px, bins=100,
           xlabel="SIF₆₈₃ (W/m²/µm/sr) - alg4 (NMF)",
	       ylabel="nFLH (W/m²/µm/sr)",
           title="NMF Retrieval vs. nFLH",
           colorbar_title="Count",
		   # xlim=( 0.0, 0.25 ),
		   # ylim=( 0.05, 0.25 ),
           color=:viridis)
end

# ╔═╡ 643d1724-3d05-467b-8cf2-2ac0bd97dd13
histogram2d(
   resd₄, log.(chl_px), bins=100,
   xlabel="Residual (l₂-norm [620 nm, 650 nm])",
   ylabel="Chl-a concentration",
   title="chlor_a vs. Residual",
   colorbar_title="Count",
   # xlim=( 0.0, 0.25 ),
   # ylim=( 0., 1. ),
   color=:viridis
)

# ╔═╡ 38afdb15-c285-46c7-b40c-0d320fd500d5
md"""
### Compare with forward model 3
---
"""

# ╔═╡ d9af203c-fb1f-49f9-8a32-386689c21245
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
	smooth_x = sigm(x[px.nPoly+px.nPC+2]);
	T₂    = @. exp( smooth_x * log(T₁) );

	# SIF magnitude
	SIF   = px.SIF_shape * x[px.nPoly+px.nPC+px.nSIF+2];
	
	# TOA radiance
	rad   = @. px.E * cosd(px.sza) / π * T₂ * ρ + SIF * T₁;
	return rad
end

# ╔═╡ 875b3251-edbc-4b38-be48-53017ac3df0f
function reconstruct3(
		px :: Pixel,       # Pixel struct
	)

	# reflectance
	v     = collectPl.(px.λc, lmax=px.nPoly);
	ρ     = hcat(v...)' * px.x[1 : px.nPoly+1];
	
	# T↑ transmittance for SIF
	T₁    = (px.trans_mat * px.x[(px.nPoly+2):(px.nPoly+px.nPC+1)]) .+ 1.0;

	# T↓↑ transmittance for SIF
	T₂    = @. exp( sigm(px.x[px.nPoly+px.nPC+2]) * log(T₁) );

	# SIF magnitude
	SIF   = px.SIF_shape * px.x[px.nPoly+px.nPC+px.nSIF+2];
	
	return ρ, T₂, T₁, SIF
end

# ╔═╡ 32d96696-6758-4268-b1d6-c75da709aba9
# SVD
HighResSVD = Spectral_SVD(trans .- 1., bands, λ_min=λ_min, λ_max=λ_max);

# ╔═╡ b94a2d76-f930-45b6-b496-c39013f6ab18
begin
	# s.d. of the total loading
	SVDloading_sd = std(HighResSVD.VarExp .* HighResSVD.Loading, dims=2);
	# priori cov: 3 for 1 constant term, 1 SIF term, and 1 factor
	Sₐ₃   = I(n+nPC+3) .+ 0.;
	# update diagonal term
	for i=1:(n+1)
	    Sₐ₃[i,i] = 1e10;
		# large variance applies no constrain to these polynomial term
	end
	# \beta
	for i=(n+2):(n+nPC+1)
		Sₐ₃[i,i]  = SVDloading_sd[i - (n+1)];
	end
	# \gamma
	Sₐ₃[n+nPC+2, n+nPC+2] = 2;
	# SIF magnitude
	Sₐ₃[end, end] = 1;
	println("Diagonal terms are: $(diag(Sₐ₃))")
end

# ╔═╡ 2fe0df06-4dcd-4a2a-bfe1-846923c60b80
begin
	# define new parameter set
	params₃ = (
		λ = oci_band, 
		E = E,
		λ_bl_ind = missing,
		λc    = λc,
		nPoly = n,
		nPC   = nPC, 
		nSIF  = 1,
		nIter = 20, 
		Sₐ    = Sₐ₃,
		βₐ    = mean(HighResSVD.VarExp .* HighResSVD.Loading, dims=2)[1:nPC],
		PrinComp      = HighResSVD.PrinComp[:, 1:nPC],
		thr_Converge  = 1e-6,
		forward_model = forward_model3,
	);
	
	SVDRetrieval = Retrieval4.(
		eachslice(R_toa[SIF_index, :], dims=1),
		sza[SIF_index],
		vza[SIF_index],
		nflh[SIF_index],      # nflh
		chlor_a[SIF_index],   # chlor_a
		nflh[SIF_index],      # flag
		Ref(params₃) 
	)
end

# ╔═╡ f3958efb-ca4a-4775-a5c3-aac8d939c331
sum(ismissing.(SVDRetrieval))

# ╔═╡ 5d614c07-3b45-458e-abc8-97efb9e6cf0e
begin
	println("make plots for `SVD-Retrieval`.")
	
	# make transmittance
	p_rho₃ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="ρ [-]",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
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
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
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
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
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
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
	)

	# null array to store that
	resd₃      = [];
	nflh_px₃   = [];
	SIF₆₈₃_px₃ = [];
	chl_px₃    = [];

	for i in 1:number_of_px
		if ismissing(SVDRetrieval[i])
			push!(resd₃, missing);
			push!(nflh_px₃, missing);
			push!(SIF₆₈₃_px₃, missing);
			push!(chl_px₃, missing);
			continue
		end
		# reconstruct
		ρ₃ⱼ, T₃₂ⱼ, T₃₁ⱼ, SIF₃ⱼ = reconstruct3(SVDRetrieval[i]);

		# push
		resdᵢ = SVDRetrieval[i].y .- SVDRetrieval[i].R_toa;
		nflhᵢ = SVDRetrieval[i].nflh;
		ā     = norm(resdᵢ[resd_bl_ind], 2);
		push!(resd₃, ā);
		push!(nflh_px₃, nflhᵢ);
		push!(SIF₆₈₃_px₃, SIF₃ⱼ[nflh_bl_ind]);
		push!(chl_px₃, SVDRetrieval[i].chlor_a);
	
		if i%100==0
			# plot
			plot!(
				p_rho₃, oci_band, ρ₃ⱼ,
				label="$(round(SVDRetrieval[i].flag, digits=2))")
			plot!(
				p_trans₃₂, oci_band, T₃₂ⱼ, 
				label="")
			plot!(
				p_trans₃₁, oci_band, T₃₁ⱼ,
				label="$(round(sigm(SVDRetrieval[i].x[SVDRetrieval[i].nPoly+SVDRetrieval[i].nPC+2]), digits=2))"*" $( round((secd(SVDRetrieval[i].sza) + secd(SVDRetrieval[i].vza)) / secd(SVDRetrieval[i].vza) , digits=2))")
			plot!(
				p_SIF₃, oci_band, SIF₃ⱼ .* T₃₁ⱼ,
				label="",
				linestyle=:dash)
			plot!(
				p_SIF₃, oci_band, SIF₃ⱼ, label="")
		end
	
	end
end

# ╔═╡ 73b635a2-1cc8-4fe3-98bf-5cbdd4ebf160
p_trans₃₂

# ╔═╡ aafef036-e6a4-471d-8fc2-14f7798dfc23
p_trans₃₁

# ╔═╡ 9df9c1ee-0ef9-4a9f-9284-e5a85ebdc557
p_SIF₃

# ╔═╡ c590844b-1529-4ce4-ad2d-1a6a094abe96
begin
	# residual
	p_resd₃ = plot(
		size=(800, 400), 
		legendcolumns=4,
		xlabel="[nm]",
		ylabel="Residual",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
	)

	for i in 1:20:number_of_px
		if ismissing(SVDRetrieval[i])
			continue
		end
		# get spectral-wise residual
		resdᵢ = SVDRetrieval[i].y .- SVDRetrieval[i].R_toa;
		
		plot!(
			p_resd₃, oci_band, resdᵢ, label="",
		)
	end
	p_resd₃
end

# ╔═╡ a0005415-37d5-4533-8ac0-68ee21373997
begin
	p_hist₃ = histogram2d(
		   SIF₆₈₃_px₃, nflh_px, bins=90,
           xlabel="SIF₆₈₃ (W/m²/µm/sr) - alg3 (SVD)",
	       ylabel="nFLH (W/m²/µm/sr)",
           title="SVD retrieval vs. nFLH",
           colorbar_title="Count",
		   # xlim=( 0.0, 0.25 ),
		   # ylim=( 0.05, 0.25 ),
           color=:viridis)
end

# ╔═╡ 1c7b4659-9d4f-41b8-8e66-1a8730a255dd
histogram2d(
   SIF₆₈₃_px₃, SIF₆₈₃_px₄, bins=150,
   xlabel="SIF₆₈₃ (W/m²/µm/sr) - alg3 (SVD)",
   ylabel="SIF₆₈₃ (W/m²/µm/sr) - alg4 (NMF)",
   title="NMF vs. SVD",
   colorbar_title="Count",
   color=:viridis
)

# ╔═╡ 24db2139-9fbe-45e2-b835-302cdfc587b9
histogram2d(
   resd₃, resd₄, bins=150,
   xlabel="Residual (l2-norm) - alg3 (SVD)",
   ylabel="Residual (l2-norm) - alg4 (NMF)",
   title="NMF vs. SVD Residual @ [620 nm, 650 nm]",
   colorbar_title="Count",
   color=:viridis
)

# ╔═╡ bbaabf07-5c02-4c26-9dfc-98671d8d181c
begin
	# make transmittance
	pᵥ = plot(
		size=(800, 400), 
		legendcolumns=3,
		xlabel="[nm]",
		ylabel="Transmittance",
		xlabelfontsize=10,
		ylabelfontsize=10,
		left_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		title="ensemble of retrieval nFLH=[$nFLH_min, $nFLH_max]"
	)

    colors = palette(:tab10)  # Use a specific color palette
    color_idx = 1
    
    for i in 1:number_of_px
        if ismissing(MyRetrieval[i]) || ismissing(SVDRetrieval[i])
            continue
        end
        
        # reconstruct
        _, T₄₂ⱼ, T₄₁ⱼ, _ = reconstruct4(MyRetrieval[i]);
        _, T₃₂ⱼ, T₃₁ⱼ, _ = reconstruct3(SVDRetrieval[i]);

        if (i%50 == 0) # && (0.10 < MyRetrieval[i].nflh < 0.25)
            # plot with same color
            plot!(pᵥ, oci_band, T₄₁ⱼ, 
				label="$(MyRetrieval[i].nflh)",
				color=colors[color_idx]
			)
            plot!(pᵥ, oci_band, T₃₁ⱼ, label="", 
				linestyle=:dash,
				color=colors[color_idx]
			)
            
            # Cycle colors (wrap around if needed)
            color_idx = mod1(color_idx + 1, length(colors))
        end
    end
	pᵥ
end

# ╔═╡ 7d2a01d0-bbb5-4025-9151-f6fb99f709db
md"""
### Some random plot
---
"""

# ╔═╡ 7a204dd7-85c6-419f-ae27-4b2d78558082
plot(oci_band, E, size=(800, 400))

# ╔═╡ 1e110787-844a-43d8-816f-bfeeeb182753
plot(oci_band, MyRetrieval[10].y, size=(800, 400)); plot!(oci_band, MyRetrieval[10].R_toa, size=(800, 400))

# ╔═╡ 45be8693-430b-470e-971c-4fef038bebf7
begin
	k     = 10;
	y_obs = MyRetrieval[k].R_toa;
	y_fit = MyRetrieval[k].y;
	plot(oci_band, y_fit ./ y_obs, size=(800, 400));
	# plot(oci_band, y_fit .- y_obs, size=(800, 400))
end

# ╔═╡ Cell order:
# ╟─64ad88f8-af61-11f0-230d-a525e3af3a49
# ╠═bb10dd17-e8be-4bf1-9caa-c2f645d07046
# ╠═d8546b0b-d901-4520-96b9-e3601f0ab70f
# ╠═3b32eb28-b9ee-4103-bdaf-f65e0cd89760
# ╠═489b1161-e108-4506-915b-07143c546d70
# ╠═293e259b-af97-476b-a8be-103c731e6d1a
# ╟─584e1386-e522-4774-adee-f0d54f55db78
# ╠═5afd8471-0ec5-454a-aeb3-1921a60bd48f
# ╠═8c4f5fa0-0929-4a55-8e52-bbcd4adfcbb1
# ╠═3e56380e-6ce7-4c34-a9c5-fa75a06a6b6a
# ╟─41ee0706-efe3-4c5a-859e-dd8f0390732b
# ╟─5e3d0da7-1910-4879-b2e0-f4a1ea32f5e2
# ╠═9cb2e8d3-2fdb-4a6d-8f81-1e7826e6cc8e
# ╠═b8ba6ff5-0d91-4201-b38e-8afd6c60fdfa
# ╟─7c4fcde6-0e72-4229-a60c-969eda6050b1
# ╠═41cc9e03-be5d-49f9-8769-2e1187ee900c
# ╠═cf4159bd-23b1-4ffd-be6b-02ec2aa6e34f
# ╟─31069e84-bde9-4377-a75e-2cc4e9c9fc49
# ╠═08343917-236b-4bfb-be42-9b560139906c
# ╠═403f6061-e5ca-40ba-9a64-afc0fe018bba
# ╠═c437a631-ab57-4b27-9550-8aa9c7e4b27e
# ╟─4fc5d644-d346-42d8-88e0-37ac0cd46ecb
# ╟─4ac45a58-1988-43e6-93af-8a27a3d80a72
# ╠═e1ade8ed-21b6-467d-a291-02ae3fe84d5e
# ╠═f0b8aed0-1455-4cd6-b42c-697287943e59
# ╠═cb604d64-a31f-4f97-bf48-a816d8abf83b
# ╠═509c5b50-098e-4d66-878a-8cbc0de50ee9
# ╠═1fc5b94b-6246-4cd0-881b-b010cef82a0c
# ╟─60c4eacc-6db1-43ad-945c-6af2147101f5
# ╠═f6155dde-0a72-4750-a82f-b045b8ada62b
# ╠═8f3b7a74-fc81-47ba-aac0-c6179fa02e4d
# ╠═06a27f81-3d65-4e39-a627-4f21b1c9830a
# ╠═73286e0b-0356-4a79-8ef0-8e6d05d6cfe7
# ╟─09b0e9ae-123a-4fb9-a1cc-4efa9d6a07e4
# ╠═c18b7fb5-c0c4-4c4e-bc69-2fe5fc0442ab
# ╠═68ee1b3b-e624-431c-8908-142a6da86074
# ╟─efe7e97e-a1b8-4c05-bd9d-5222f4678d5c
# ╟─e6f7310a-0e11-4e67-933a-f32292f45047
# ╟─b678620c-4a20-43f9-b22e-7ffd199bde77
# ╟─43ad9bd4-dafd-47a1-99c5-d078e3381c4f
# ╟─d297a09a-5b97-43ef-ae01-367444ad87fe
# ╟─643d1724-3d05-467b-8cf2-2ac0bd97dd13
# ╟─38afdb15-c285-46c7-b40c-0d320fd500d5
# ╟─d9af203c-fb1f-49f9-8a32-386689c21245
# ╟─875b3251-edbc-4b38-be48-53017ac3df0f
# ╠═32d96696-6758-4268-b1d6-c75da709aba9
# ╠═b94a2d76-f930-45b6-b496-c39013f6ab18
# ╠═2fe0df06-4dcd-4a2a-bfe1-846923c60b80
# ╠═f3958efb-ca4a-4775-a5c3-aac8d939c331
# ╠═5d614c07-3b45-458e-abc8-97efb9e6cf0e
# ╟─73b635a2-1cc8-4fe3-98bf-5cbdd4ebf160
# ╟─aafef036-e6a4-471d-8fc2-14f7798dfc23
# ╟─9df9c1ee-0ef9-4a9f-9284-e5a85ebdc557
# ╟─c590844b-1529-4ce4-ad2d-1a6a094abe96
# ╟─a0005415-37d5-4533-8ac0-68ee21373997
# ╟─1c7b4659-9d4f-41b8-8e66-1a8730a255dd
# ╟─24db2139-9fbe-45e2-b835-302cdfc587b9
# ╟─bbaabf07-5c02-4c26-9dfc-98671d8d181c
# ╟─7d2a01d0-bbb5-4025-9151-f6fb99f709db
# ╠═7a204dd7-85c6-419f-ae27-4b2d78558082
# ╠═1e110787-844a-43d8-816f-bfeeeb182753
# ╠═45be8693-430b-470e-971c-4fef038bebf7
