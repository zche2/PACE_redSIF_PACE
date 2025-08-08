### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 99ff878a-6e71-11f0-17ed-b7ac188e90b8
import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE");

# ╔═╡ 04c805e7-45b5-4878-b288-0cf1d02d31fc
using Polynomials, ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics

# ╔═╡ 0b112f15-6cc7-4f02-849e-e0ef8a71b639
using LegendrePolynomials

# ╔═╡ 922ddadd-a129-406d-9de6-892899786e73
using Parameters

# ╔═╡ 0ec3629f-0278-42b1-8ab8-f399d4d4f216
include("../PACE_SIF.jl")

# ╔═╡ 05837924-482b-4564-a770-3544f736889b
md"""
> #### Load transmittance spectra and do SVD
"""

# ╔═╡ 379babe3-7d99-431b-b5db-499ee9b5b406
begin
	summer = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel.nc");
	winter = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel.nc");
	println("Opened datasets.")

	
	temp  = cat(summer["temperature"][:,:], winter["temperature"][:,:], dims=1);
	psurf = cat(summer["pressure"][:], winter["pressure"][:], dims=1);
	q     = cat(summer["q"][:,:], winter["q"][:,:], dims=1);
	AMF   = cat(summer["AMF"][:], winter["AMF"][:], dims=1);
	trans = cat(summer["transmittance"][:,:], winter["transmittance"][:,:], dims=1);
	println("\nConcatenated!")

	bands  = summer["band"][:];
	
end

# ╔═╡ 6f24e4fe-94b5-45bd-bf46-a98a0fdbaf48
begin
	λ_min = 620.;
	λ_max = 860.;
	# get principal components, variance explained by each component (normalized to 100%), and spatial loading
	HighResSVD = PACE_SIF.Spectral_SVD(trans, bands, λ_min=λ_min, λ_max=λ_max);
end

# ╔═╡ 401b62ff-9966-40b7-ac5d-ed5d704ddda3
mean(HighResSVD.PrinComp[:,1:4], dims=1)

# ╔═╡ 0ccf2db1-9080-4d29-bfc7-11dffa706f62
md"""
> ##### sample data from PACE-OCI
"""

# ╔═╡ a42bd26f-46d5-44a4-81d8-7788899b95bc
begin
	oci = Dataset(
		"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample_swath_20250130T20.nc");
	# findall(coalesce.(nflh .> 0.4, false))
	pixel  = 769;  # cross-track
	scan   = 1551;
	E      = oci["red_solar_irradiance"][:];
	red_band = oci["red_bands"][:];
	cloud    = oci["cloud_flag_dilated"][:, :];
	nflh     = oci["nflh"][:, :];
	# select the pixel
	ind      = findall( λ_min .< red_band .< λ_max );
	E        = oci["red_solar_irradiance"][ind];
	oci_band = red_band[ind];
	rhot     = oci["rhot_red"][pixel, scan, ind];
	R_TOA    = oci["radiance_red"][pixel, scan, ind];
	vza      = oci["sensor_zenith"][pixel, scan];
	sza      = oci["solar_zenith"][pixel, scan];
	println("Data read!")
	@show oci["chlor_a"][pixel, scan];
	@show nflh[pixel, scan];
	
end

# ╔═╡ acacde64-9957-409d-ae67-428d13428e9d
begin
	# the PCs look like:
	plot(oci_band, HighResSVD.PrinComp[:,1:4], size=(500, 200))
	title!("eigen vectors")
end

# ╔═╡ 063343a5-5879-4cb7-91ad-5068fe0b33d2
md"""
> ##### SNR -> measurement covariance matrix $S_{\epsilon}$
"""

# ╔═╡ 466c6800-dd8d-4b11-b33b-bff17dfcf387
begin
	filename = raw"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt";
	lines = readlines(filename);
	end_header_index = findfirst(x -> x == "/end_header", lines);
	data = readdlm(filename, String, skipstart=end_header_index);

	FPA   = data[:, 1];                   # 1st column: band
	wvlen = parse.(Float64, data[:, 2]);  # 2nd column: center wavelength
	c1    = parse.(Float64, data[:, 4]);  # 4th column: c1
	c2    = parse.(Float64, data[:, 5]);  # 5th column: c2

	snr_ind = findall((FPA .== "Red") .& (λ_min .< wvlen .< λ_max));
	# see instruction in .txt file
	noise   = sqrt.( c1[snr_ind] .+ c2[snr_ind] .* R_TOA);
	Se      = Diagonal(noise.^2);
end

# ╔═╡ f80f7a81-000a-4784-9d10-713406303102
plot(oci_band, R_TOA, size=(500, 200), label="obs"); ylabel!("W/m2/µm/sr")

# ╔═╡ 40796346-a6e5-49e4-a619-b28e5fda2521
md"""
> ##### Start with polynomial fit
$$\rho_{s}=\sum{a_jP_j}$$
$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)}{\pi}$$
Where $P_j$ are Legendre polynomials. This is a linear problem with Jacobian equal to $P_j$. $\pi$ converts flux to radiance.
"""

# ╔═╡ e61f28fc-25bd-41c8-9290-a849e43d5776
function PolyFit(x, λ;
	    		 degree::Int = 3,
				 λ_min=λ_min,
				 λ_max=λ_max,
		)
		# adjust to [-1,1]
		λm    = mean(λ);
		range = λ_max - λ_min;
		λc    = (λ .- λm) ./ range;
		# println(λc)
		# collection of Legendre polys.
		v = collectPl.(λc, lmax=degree);
		myFit = x * hcat(v...);
		return myFit, hcat(v...)
end

# ╔═╡ d42e7040-06ec-4b94-bece-78ab4728cfc8
begin
	n     = 10;
	λm    = mean(oci_band);
	range = λ_max - λ_min;
	λc    = (oci_band .- λm) ./ range;
	K     = E .* cosd(sza) ./ pi .* hcat(collectPl.(λc, lmax=n)...)';
	G     = inv( K'inv(Se)K )K'inv(Se);
	x̂     = G * R_TOA;
	ŷ     = K * x̂;
end

# ╔═╡ 9075c8bf-efaf-4836-9acc-f8de96fb840d
1e7 / 900

# ╔═╡ 6feb3c6c-cd8c-469a-bb41-cde87a967f04
println("state vectors: $x̂")

# ╔═╡ c7845a6b-c7de-4bee-84b7-8bd922dccae4
begin
	plot(oci_band, R_TOA, size=(500, 200), label="obs.")
	plot!(oci_band, ŷ, label="Legendre polynomial fit n=$n")
	title!("Retrieved and observed TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ 073ae8ca-e8f3-41e1-b460-cf00d51cf2fb
begin
	plot(oci_band,
		R_TOA ./ (E .* cosd(sza) ) .* pi,
		size=(500, 200),
		label="obs")
	plot!(oci_band, ŷ ./ (E .* cosd(sza) ) .* pi, label="poly.fit ρₜ")
	title!("Total reflectance", titlefont=10)
end

# ╔═╡ 434ee765-087e-456a-9696-2ba87fa3c5f3
md"""
> ##### Add transmittance
$$\rho_{s}=\sum{a_jP_j}$$
$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)T(\lambda)}{\pi}$$
 $T(\lambda)$ is set to have a maximum of 1.
"""

# ╔═╡ 3cb579f3-c9c4-48b3-997d-967f4e1df546
begin
	# define priori error matrix
	nPC = 7; # use 2 eigen vectors
	# priori cov
	Sa  = zeros(n+nPC+1, n+nPC+1);
	# uodate diagonal term
	for i=1:(n+1)
	    Sa[i,i] = 1e20;     
		# large variance applies no constrain to these polynomial term
	end
	rel_error   = .001
	for i=(n+2):(n+nPC+1)
	    Sa[i,i] = rel_error .* HighResSVD.VarExp[i - (n+1)];
	end
	Sa
end

# ╔═╡ c4d3782c-f85d-492e-a805-61d6f98fb657
function forward_model1(
			x; 
			λ = oci_band,     # wavelength range
			nPoly::Int = n,   # degree of polynomials
			nPC::Int   = nPC,   # number of eigen vectors used
			trans_mat  = HighResSVD.PrinComp[:, 1:nPC],
			sza        = sza,
			vza        = vza,
			E          = E,
		)
	
	# adjust to [-1,1]
	λm    = mean(λ);
	range = λ_max - λ_min;
	λc    = (λ .- λm) ./ range;
	v     = collectPl.(λc, lmax=nPoly);
	# reflectance
	rho   = hcat(v...)' * x[1 : nPoly+1];
	# transmittance
	T     = trans_mat * x[(nPoly+2):(nPoly+nPC+1)];
	T_min = minimum(T);
	T_max = maximum(T);
	factor = maximum(abs.([T_min, T_max]))
	T_norm = abs.(T) / factor
	# TOA radiance
	rad    = E .* cosd(sza) ./ pi .* T_norm .* rho
	return rad
end


# ╔═╡ d5cfaed6-0063-4649-83da-a64727487741
begin
	xa = [x̂... -6. .05 .001 .001 .001 .001 .001]';
	rad = forward_model1(xa, nPoly=n);
	plot(oci_band, R_TOA, size=(500, 200), label="obs.")
	plot!(oci_band, rad, label="initial guess, n=$n")
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ d0cbf663-ac73-413a-951c-f99bf8d2cd8d
md"""
🟢 defining functions to calculate Jacobian, gain matrix and do the iteration
"""

# ╔═╡ dd2cd8cb-ed9e-4b6b-af99-59fe26809d39
function Jacobian(x, model; len=length(oci_band))
	res = DiffResults.JacobianResult(zeros(len), x);
	ForwardDiff.jacobian!(res, model, x);
	K   = DiffResults.jacobian(res);
	val = DiffResults.value(res);
	return K, val
end

# ╔═╡ 7dcb675f-fd35-46ed-ba58-82df3d68627b
function GainMatrix(K, Se=Se, Sa=Sa)
	return inv( K'inv(Se)K + inv(Sa) )K'inv(Se)
end

# ╔═╡ 3923d033-4639-43a3-a693-8d77c04dd186
@with_kw struct retrieval
    x    # state vector
	y_x  # value evaluated at x
	K    # Jacobian
	G = GainMatrix(K) # Gain matrix
	A = G*K                   # averaging kernel
end;

# ╔═╡ a512b192-5d5e-4688-8b84-f0bc27aa18e7
function iter(
		m ::retrieval,   # retrieval of last iteration
		xa,              # priori
	    rad;   		     # measurements
		model = x -> forward_model1(x)
	)
	# get results from last iteration x̂ₙ, note that K and y are evaluated at x̂ₙ
	xn     = m.x
	Kn     = m.K
	yn     = m.y_x
	G      = m.G
	x̂      = xa .+ G * (rad .- yn .+ Kn * (xn .- xa));
	K_n1, y_n1 = Jacobian(x̂, model);
	# update 
	m_new  = retrieval(
		x   = x̂,
		y_x = y_n1,
		K   = K_n1,
	)
	return m_new
end

# ╔═╡ 93c48028-a4bb-4d6d-9bc4-85749a675793
md"""
> ##### first iter.
"""

# ╔═╡ bdcc5bf7-7ab0-43a2-8710-09b4b4366b1a
begin
	# start from xa
	Ka, ya = Jacobian(xa, x -> forward_model1(x))
	ma = retrieval(x=xa, y_x=ya, K=Ka)
end

# ╔═╡ b621fa58-9f13-48a2-9144-b3a3cb5292ac
begin
	# 1st iteration
	m1 = iter(ma, xa, R_TOA);
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=4, linestyle=:dash, color=:black)
	plot!(oci_band, ma.y_x, label="initial guess, n=$n", linewidth=2)
	plot!(oci_band, m1.y_x, label="iter#1, n=$n", linewidth=2)
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ 7b0a281d-daaa-4aaa-a001-12be469225f9
md"""
🟢 recover transmittance
"""

# ╔═╡ c17a958d-fec3-445a-ba1f-59f65ad63af6
begin
	T1     = HighResSVD.PrinComp[:, 1:nPC] * m1.x[(n+2):(n+nPC+1)];
	T1_min = minimum(T1);
	T1_max = maximum(T1);
	factor1 = maximum(abs.([T1_min, T1_max]))
	T1_norm = abs.(T1) / factor1
	# @show T1_norm

	plot(oci_band, T1_norm, size=(500, 200), label="iter#1")
	
end

# ╔═╡ 97672495-e0b1-4952-9f84-a26e926c7235
md"""
🟣 reflectance
"""

# ╔═╡ 3d80255b-7409-4d8b-9fb7-b05ed286b18a
begin
	v1   = collectPl.(λc, lmax=n);
	rho1 = hcat(v1...)' * m1.x[1 : n+1];
	plot(oci_band, ŷ ./ (E .* cosd(sza) ) .* pi,
		label="linear fit ρₜ",
		size=(500, 200))
	plot!(oci_band, rho1, label="iter#1 ρₜ")
	title!("Total reflectance", titlefont=10)
end

# ╔═╡ 556e3e8b-aae5-4462-9aab-1f5c3f90c5a4
md"""
> ##### 2nd iter.
"""

# ╔═╡ abb9b4e8-9c9c-4d82-8190-06ededcbfd52
begin
	# 2nd iteration
	m2 = iter(m1, xa, R_TOA);
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=4, linestyle=:dash, color=:black)
	plot!(oci_band, ma.y_x, label="initial guess, n=$n", linewidth=2)
	plot!(oci_band, m1.y_x, label="iter#1, n=$n", linewidth=2)
	plot!(oci_band, m2.y_x, label="iter#2, n=$n", linewidth=2)
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ 2a4b61f9-328a-4e92-ae84-58bdda55dc74
begin
	T2     = HighResSVD.PrinComp[:, 1:nPC] * m2.x[(n+2):(n+nPC+1)];
	T2_min = minimum(T2);
	T2_max = maximum(T2);
	factor2 = maximum(abs.([T2_min, T2_max]))
	T2_norm = abs.(T2) / factor2
	# @show T1_norm
	
	plot(oci_band, T1_norm, size=(500, 200), label="iter#1")
	plot!(oci_band, T2_norm, label="iter#2")
	title!("transmittance")
end

# ╔═╡ 33a4a5b0-ae07-4536-9c45-a2043d136f9f
begin
	v2   = collectPl.(λc, lmax=n);
	rho2 = hcat(v2...)' * m2.x[1 : n+1];
	plot(oci_band, ŷ ./ (E .* cosd(sza) ) .* pi,
		label="linear fit ρₜ",
		size=(500, 200))
	plot!(oci_band, rho1, label="iter#1 ρₜ", lw=2)
	plot!(oci_band, rho2, label="iter#2 ρₜ")
	title!("Total reflectance", titlefont=10)
end

# ╔═╡ 93b65f52-c5a5-4580-a64b-5a50a44208af
begin
	# resildual
	plot(oci_band, R_TOA .- m1.y_x, size=(500, 200), label="residual, iter#1")
	plot!(oci_band, R_TOA .- m2.y_x, label="residual, iter#2")
	title!("Residuals", titlefont=10)
end

# ╔═╡ 26aab1d6-706f-4b9d-9447-883528f9ed09
md"""
> ##### 3rd iter.
"""

# ╔═╡ 6adb8f15-c885-4f01-b9e9-34725602ac1e
m3 = iter(m2, xa, R_TOA);

# ╔═╡ 81006ea3-792a-4a6e-b7af-88069e2a87bd
begin
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=4, linestyle=:dash, color=:black)
	plot!(oci_band, ma.y_x, label="initial guess, n=$n", linewidth=2)
	plot!(oci_band, m1.y_x, label="iter#1, n=$n", linewidth=2)
	plot!(oci_band, m2.y_x, label="iter#2, n=$n", linewidth=2)
	plot!(oci_band, m3.y_x, label="iter#3, n=$n", linewidth=2)
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ c5c8d823-e757-46ac-b315-4c166ec12903
begin
	# resildual
	plot(oci_band, R_TOA .- m1.y_x, size=(500, 200), label="residual, iter#1")
	plot!(oci_band, R_TOA .- m2.y_x, label="residual, iter#2")
	plot!(oci_band, R_TOA .- m3.y_x, label="residual, iter#3")
	title!("Residuals", titlefont=10)
end

# ╔═╡ Cell order:
# ╠═99ff878a-6e71-11f0-17ed-b7ac188e90b8
# ╠═04c805e7-45b5-4878-b288-0cf1d02d31fc
# ╠═0b112f15-6cc7-4f02-849e-e0ef8a71b639
# ╠═922ddadd-a129-406d-9de6-892899786e73
# ╠═0ec3629f-0278-42b1-8ab8-f399d4d4f216
# ╟─05837924-482b-4564-a770-3544f736889b
# ╠═379babe3-7d99-431b-b5db-499ee9b5b406
# ╠═6f24e4fe-94b5-45bd-bf46-a98a0fdbaf48
# ╟─acacde64-9957-409d-ae67-428d13428e9d
# ╠═401b62ff-9966-40b7-ac5d-ed5d704ddda3
# ╟─0ccf2db1-9080-4d29-bfc7-11dffa706f62
# ╠═a42bd26f-46d5-44a4-81d8-7788899b95bc
# ╟─063343a5-5879-4cb7-91ad-5068fe0b33d2
# ╠═466c6800-dd8d-4b11-b33b-bff17dfcf387
# ╠═f80f7a81-000a-4784-9d10-713406303102
# ╟─40796346-a6e5-49e4-a619-b28e5fda2521
# ╠═e61f28fc-25bd-41c8-9290-a849e43d5776
# ╠═d42e7040-06ec-4b94-bece-78ab4728cfc8
# ╠═9075c8bf-efaf-4836-9acc-f8de96fb840d
# ╠═6feb3c6c-cd8c-469a-bb41-cde87a967f04
# ╟─c7845a6b-c7de-4bee-84b7-8bd922dccae4
# ╟─073ae8ca-e8f3-41e1-b460-cf00d51cf2fb
# ╟─434ee765-087e-456a-9696-2ba87fa3c5f3
# ╠═c4d3782c-f85d-492e-a805-61d6f98fb657
# ╠═d5cfaed6-0063-4649-83da-a64727487741
# ╠═3cb579f3-c9c4-48b3-997d-967f4e1df546
# ╟─d0cbf663-ac73-413a-951c-f99bf8d2cd8d
# ╠═dd2cd8cb-ed9e-4b6b-af99-59fe26809d39
# ╠═7dcb675f-fd35-46ed-ba58-82df3d68627b
# ╠═3923d033-4639-43a3-a693-8d77c04dd186
# ╠═a512b192-5d5e-4688-8b84-f0bc27aa18e7
# ╟─93c48028-a4bb-4d6d-9bc4-85749a675793
# ╠═bdcc5bf7-7ab0-43a2-8710-09b4b4366b1a
# ╠═b621fa58-9f13-48a2-9144-b3a3cb5292ac
# ╟─7b0a281d-daaa-4aaa-a001-12be469225f9
# ╟─c17a958d-fec3-445a-ba1f-59f65ad63af6
# ╟─97672495-e0b1-4952-9f84-a26e926c7235
# ╟─3d80255b-7409-4d8b-9fb7-b05ed286b18a
# ╟─556e3e8b-aae5-4462-9aab-1f5c3f90c5a4
# ╠═abb9b4e8-9c9c-4d82-8190-06ededcbfd52
# ╟─2a4b61f9-328a-4e92-ae84-58bdda55dc74
# ╟─33a4a5b0-ae07-4536-9c45-a2043d136f9f
# ╠═93b65f52-c5a5-4580-a64b-5a50a44208af
# ╟─26aab1d6-706f-4b9d-9447-883528f9ed09
# ╠═6adb8f15-c885-4f01-b9e9-34725602ac1e
# ╟─81006ea3-792a-4a6e-b7af-88069e2a87bd
# ╠═c5c8d823-e757-46ac-b315-4c166ec12903
