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

# ╔═╡ 3a475f5d-b7f2-4dad-9335-82b6bf6e368b
using NonlinearSolve, BenchmarkTools

# ╔═╡ 0ec3629f-0278-42b1-8ab8-f399d4d4f216
include("/home/zhe2/FraLab/PACE_redSIF_PACE/PACE_SIF.jl")

# ╔═╡ 857eaa38-cc95-42d7-82f6-853ffa39dfe6
md"""
> ## Two solutions to Fix Non-convergence in V2.1
- "Baseline fit": bypass O2 B-band and SIF
- Scale transmittance according to baseline region: Avoid absorption feature
"""

# ╔═╡ 05837924-482b-4564-a770-3544f736889b
md"""
> #### Load transmittance spectra and do SVD
"""

# ╔═╡ 379babe3-7d99-431b-b5db-499ee9b5b406
begin
	summer = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc");
	winter = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc");
	println("Opened datasets.")

	
	temp  = cat(summer["temperature"][:,:], winter["temperature"][:,:], dims=1);
	psurf = cat(summer["pressure"][:], winter["pressure"][:], dims=1);
	q     = cat(summer["q"][:,:], winter["q"][:,:], dims=1);
	AMF   = cat(summer["AMF"][:], winter["AMF"][:], dims=1);
	trans = cat(summer["transmittance"][:,:], winter["transmittance"][:,:], dims=1);
	println("\nConcatenated!")

	bands  = summer["band"][:];

	close(summer);
	close(winter);
end

# ╔═╡ 3ac1d3eb-a22b-441c-8343-062f1d733779
bands

# ╔═╡ 6f24e4fe-94b5-45bd-bf46-a98a0fdbaf48
begin
	λ_min = 610.;
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
		"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample_granule_20250501T183011_new.nc");
	pixel  = 456;  # cross-track
	scan   = 901;
	red_band = oci["red_wavelength"][:];
	# red_band = oci["red_bands"][:];
	# cloud    = oci["cloud_flag_dilated"][:, :];
	nflh     = oci["nflh"][:, :];
end

# ╔═╡ cc1acba5-d114-4579-a64f-8546c2df40b1
begin
	# select band (continuum spectrum)
	ind      = findall( λ_min .< red_band .< λ_max );
	E        = oci["red_solar_irradiance"][ind];
	oci_band = red_band[ind];
	rhot     = oci["rhot_red"][pixel, scan, ind];
	R_TOA    = oci["radiance_red"][pixel, scan, ind];
	vza      = oci["sensor_zenith"][pixel, scan];
	sza      = oci["solar_zenith"][pixel, scan];
	println("Data read!")
	# @show oci["chlor_a"][pixel, scan];
	@show nflh[pixel, scan];

	# close
	close(oci)
end

# ╔═╡ 672286a7-5b44-49f3-8098-9371f5928826
begin
	# select fitting band
	λ_left  = 650.; # 670.;
	λ_right = 640.; # 710.;
	ind_fit = findall( ( oci_band .< λ_left) .| (λ_right .< oci_band) );
	oci_band_fit = oci_band[ind_fit];
	R_TOA_fit    = R_TOA[ind_fit];
end

# ╔═╡ f80f7a81-000a-4784-9d10-713406303102
begin
	p1 = plot(oci_band, R_TOA, size=(500, 300), label="obs");
	scatter!(p1, oci_band_fit, R_TOA_fit, label="fitting band (TBD)", markersize=1.5);

	p2 = plot(oci_band, E, size=(500, 300), label="obs");
	scatter!(p2, oci_band_fit, E[ind_fit], label="solar irr.", markersize=1.5)

	plot(p1, p2, layout=(2, 1))
	ylabel!("W/m2/µm/sr")
end

# ╔═╡ acacde64-9957-409d-ae67-428d13428e9d
begin
	# the PCs look like:
	plot(oci_band, HighResSVD.PrinComp[:,1:4], size=(500, 200))
	title!("eigen vectors")
end

# ╔═╡ 0d68673e-5d07-4703-96f6-d1a4ef919f0e
findall(coalesce.(nflh .> 1., false))

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

	wv_val  = (λ_min .< wvlen .< λ_left) .| (λ_right .< wvlen .< λ_max);
	snr_ind = findall((FPA .== "Red") .& wv_val);

	# to make sure I use the right wvlen
	@show wvlen[snr_ind]
	
	# see instruction in .txt file
	noise   = sqrt.( c1[snr_ind] .+ c2[snr_ind] .* R_TOA_fit);
	Se      = Diagonal(noise.^2);
end

# ╔═╡ 3c5964a2-c1e5-4dff-9932-5db894771191
begin
	# compare to make sure wavelengths are consistent
	scatter(oci_band_fit, wvlen[snr_ind], label="fit wv", size=(300, 300))
	plot!([610, 850], [610, 850], label="1:1")
end

# ╔═╡ 434ee765-087e-456a-9696-2ba87fa3c5f3
md"""
> ##### Start with polynomial fit +Transmittance
$$\rho_{s}=\sum{a_jP_j}$$
$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)T(\lambda)}{\pi}$$
 $T(\lambda)$ is set to have a maximum of 1.
"""

# ╔═╡ 40253fb3-981f-4c2d-9f43-ce1c802fc6ef
md"""
🟢 For polynomial term, the argument needs to satisfy -1 <= x <= 1.
"""

# ╔═╡ 9dcc1616-91d6-45d8-9873-2c449b6e321e
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

# ╔═╡ b2dbaa15-ac40-40cd-b022-d49193febaa9
bl_wvlen=[649.599976, 650.900024, 652.099976, 653.299988, 654.599976, 655.799988, 657.099976, 658.299988, 659.599976, 710.500000, 711.799988, 713.000000, 714.299988, 716.799988, 719.200012]

# ╔═╡ ab74fe0c-cfa8-45fc-b4fd-8fea3f93c51b
function scale_transmittance(T, λ; 
							 bl_wvlen=[649.599976, 650.900024, 652.099976, 653.299988, 654.599976, 655.799988, 657.099976, 658.299988, 659.599976, 710.500000, 711.799988, 713.000000, 714.299988, 716.799988, 719.200012]
		)
	# Use a one-liner to find the indices
	closest_ind = map(bl_wvlen -> argmin(abs.(λ .- bl_wvlen)), bl_wvlen);
	# find max
	bl_mean = mean(T[closest_ind]);
	# force the mean val to be 1
	T_norm = T ./ bl_mean
	return T_norm
end

# ╔═╡ e59a4998-c578-42c3-b4e8-61585544f69b
begin
	# number of polynormial terms and PCs used
	n = 5; nPC = 10;
	# inital guess 
	λc    = center_wavelength(oci_band_fit)
	K     = E[ind_fit] .* cosd(sza) ./ pi .* hcat(collectPl.(λc, lmax=n)...)';
	G     = inv( K'inv(Se)K )K'inv(Se);
	x̂     = G * R_TOA_fit;
	ŷ     = K * x̂;
end

# ╔═╡ 3948942d-e754-445a-aa4f-e7dc79537822
begin
	# define priori error matrix
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
	λc    = center_wavelength(λ)
	v     = collectPl.(λc, lmax=nPoly);
	# reflectance
	rho   = hcat(v...)' * x[1 : nPoly+1];
	# transmittance
	T      = trans_mat * x[(nPoly+2):(nPoly+nPC+1)];
	T_norm = scale_transmittance(T, λ)
	# TOA radiance
	rad    = E .* cosd(sza) ./ pi .* T_norm .* rho
	return rad
end


# ╔═╡ 407cf364-32c6-4d9f-9596-6a04bbd5a588
begin
	# wrap the func
	forward_model1_wrap = x -> forward_model1(x,
		λ = oci_band_fit, 
		nPoly = n, 
		trans_mat = HighResSVD.PrinComp[ind_fit, 1:nPC],
		E = E[ind_fit],
	)
end

# ╔═╡ d5cfaed6-0063-4649-83da-a64727487741
begin
	tmp  = zeros(nPC-2) .+ .001;
	xa   = [x̂... -6. .05 tmp...]';
	rad  = forward_model1_wrap(xa);
	plot(oci_band, R_TOA, size=(500, 200), label="obs.")
	scatter!(oci_band_fit, rad, 
			 label="initial guess, n=$n, nPC=$nPC",
		     markersize=1.5
	)
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
function GainMatrix(K; Se=Se, Sa=Sa)
	return inv( K'inv(Se)K + inv(Sa) )K'inv(Se)
end

# ╔═╡ 3923d033-4639-43a3-a693-8d77c04dd186
@with_kw struct retrieval
    x    # state vector
	y_x  # value evaluated at x
	K    # Jacobian
	G = GainMatrix(K)         # Gain matrix
	A = G*K                   # averaging kernel
end;

# ╔═╡ a512b192-5d5e-4688-8b84-f0bc27aa18e7
function iter(
		m ::retrieval,   # retrieval of last iteration
		xa,              # priori
	    rad; 
		Sa    = Sa,      # measurements
		model = forward_model1_wrap
	)
	
	# get results from last iteration x̂ₙ, note that K and y are evaluated at x̂ₙ
	xn     = m.x
	Kn     = m.K
	yn     = m.y_x
	G      = m.G
	x̂      = xa .+ G * (rad .- yn .+ Kn * (xn .- xa));
	K_n1, y_n1 = Jacobian(x̂, model, len=length(rad));

	# update 
	m_new  = retrieval(
		x   = x̂,
		y_x = y_n1,
		K   = K_n1,
		G   = GainMatrix(K_n1, Sa=Sa)
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
	Ka, ya = Jacobian(xa, forward_model1_wrap, len=length(oci_band_fit))
	ma = retrieval(x=xa, y_x=ya, K=Ka)
end

# ╔═╡ b621fa58-9f13-48a2-9144-b3a3cb5292ac
begin
	# 1st iteration
	m1 = iter(ma, xa, R_TOA_fit);
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=4, linestyle=:dash, color=:black)
	scatter!(oci_band_fit, ma.y_x, label="initial guess, n=$n", markersize=1.5)
	scatter!(oci_band_fit, m1.y_x, label="iter#1, n=$n", markersize=1.5)
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ e9bc8ce0-14a1-4cbe-9df0-c5b5098ecede
md"""
> ##### The trick here is to believe (close-to) linear fit out side of SIF/O$_2$ B-Band, and recover the O$_2$ absorption within the band
"""

# ╔═╡ 7b0a281d-daaa-4aaa-a001-12be469225f9
md"""
🟢 recover transmittance
"""

# ╔═╡ c17a958d-fec3-445a-ba1f-59f65ad63af6
begin
	T1     = HighResSVD.PrinComp[:, 1:nPC] * m1.x[(n+2):(n+nPC+1)];
	T1_norm = scale_transmittance(T1, oci_band_fit)
	# @show T1_norm

	plot(oci_band, T1_norm, size=(500, 200), label="iter#1")
	
end

# ╔═╡ 97672495-e0b1-4952-9f84-a26e926c7235
md"""
🟣 reflectance
"""

# ╔═╡ 3d80255b-7409-4d8b-9fb7-b05ed286b18a
begin
	v1   = collectPl.(center_wavelength(oci_band), lmax=n);
	rho1 = hcat(v1...)' * m1.x[1 : n+1];
	plot(oci_band, rho1, label="iter#1 ρₜ", size=(500, 200))
	title!("Total reflectance", titlefont=10)
end

# ╔═╡ 0d03bf4e-64d7-4f52-b52c-8ad17a93157c
# fit the full spectral range and get the residual
spectra1 = forward_model1(m1.x); residual1 = R_TOA .- spectra1;

# ╔═╡ aed98d25-2b7a-4755-bb5a-3acbd1bae0a4
begin
	# fit the full spectral range and get the residual
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=4,
		linestyle=:dash, color=:black)
	plot!(oci_band, spectra1, label="iter#1, n=$n")
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ f8320953-b1b2-4954-8215-2fa6f27cb87e
begin
	plot(oci_band, residual1, label="iter#1, n=$n", size=(500, 200))
	title!("Residual (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ 556e3e8b-aae5-4462-9aab-1f5c3f90c5a4
md"""
> ##### 2nd iter.
"""

# ╔═╡ d8810003-05c7-495f-b4d1-77a057698d2e
m1.x

# ╔═╡ abb9b4e8-9c9c-4d82-8190-06ededcbfd52
begin
	# 2nd iteration
	m2 = iter(m1, xa, R_TOA_fit);
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=4, linestyle=:dash, color=:black)
	scatter!(oci_band_fit, ma.y_x, label="initial guess, n=$n", markersize=1.5)
	scatter!(oci_band_fit, m1.y_x, label="iter#1, n=$n", markersize=1.5)
	scatter!(oci_band_fit, m2.y_x, label="iter#2, n=$n", markersize=1.5)
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ c9d962f6-8722-4faa-b34f-092de7a76bcf
m2.x

# ╔═╡ 2a4b61f9-328a-4e92-ae84-58bdda55dc74
begin
	T2     = HighResSVD.PrinComp[:, 1:nPC] * m2.x[(n+2):(n+nPC+1)];
	# T2_min = minimum(T2);
	# T2_max = maximum(T2);
	# factor2 = maximum(abs.([T2_min, T2_max]))
	T2_norm = scale_transmittance(T2, oci_band_fit)
	# @show T1_norm
	
	plot(oci_band, T1_norm, size=(500, 200), label="iter#1")
	plot!(oci_band, T2_norm, label="iter#2")
	title!("transmittance")
end

# ╔═╡ 33a4a5b0-ae07-4536-9c45-a2043d136f9f
begin
	v2   = collectPl.(center_wavelength(oci_band), lmax=n);
	rho2 = hcat(v2...)' * m2.x[1 : n+1];
	plot(oci_band, rho1, label="iter#1 ρₜ", lw=2, size=(500, 200))
	plot!(oci_band, rho2, label="iter#2 ρₜ")
	title!("Total reflectance", titlefont=10)
end

# ╔═╡ 93b65f52-c5a5-4580-a64b-5a50a44208af
begin
	spectra2  = forward_model1(m2.x);
	residual2 = R_TOA .- spectra2;
	# resildual
	plot(oci_band, residual1, size=(500, 200), label="residual, iter#1")
	plot!(oci_band, residual2, label="residual, iter#2")
	title!("Residuals", titlefont=10)
end

# ╔═╡ 0ff68f72-cc18-415b-a7b9-b94d49ee74dd
#=╠═╡
begin
	plot(oci_band, ma_new.x[end-nSIF+1] .* SIF_shape(oci_band, λ₀=ma_new.x[end-nSIF+2], σ=ma_new.x[end-nSIF+3]), label="initial guess (nFLH)", size=(500, 150))
	plot!(oci_band, m1_new.x[end-nSIF+1] .* SIF_shape(oci_band, λ₀=m1_new.x[end-nSIF+2], σ=m1_new.x[end-nSIF+3]), label="iter#1")
	plot!(oci_band, m2_new.x[end-nSIF+1] .* SIF_shape(oci_band, λ₀=m2_new.x[end-nSIF+2], σ=m2_new.x[end-nSIF+3]), label="iter#2")
	plot!(oci_band, m3_new.x[end-nSIF+1] .* SIF_shape(oci_band, λ₀=m3_new.x[end-nSIF+2], σ=m3_new.x[end-nSIF+3]), label="iter#3")
	plot!(oci_band, m4_new.x[end-nSIF+1] .* SIF_shape(oci_band, λ₀=m4_new.x[end-nSIF+2], σ=m4_new.x[end-nSIF+3]), label="iter#4")
	title!("retrieved SIF")
end
  ╠═╡ =#

# ╔═╡ 0bf97b73-04c9-4eb5-906a-23827a2c5f3a
md"""
> ##### explicitly define loss function
"""

# ╔═╡ 98bbf74e-6d47-4f25-b060-3f3c6d289a1a
# ╠═╡ disabled = true
#=╠═╡
function loss_function(x, p)
	# p is params for forward model and error matrix
	λ     = p.λ;       # wavelength range
	nPoly = p.nPoly;   # degree of polynomials
	nPC   = p.nPC;     # number of eigen vectors used
	nSIF  = p.nSIF;
	trans_mat = p.trans_mat;
	E     = p.E;
	xa    = p.xa;      # priori
	ŷ     = p.ŷ;
	sza   = p.sza;
	vza   = p.vza;
	Se    = p.Se;
	Sa    = p.Sa;

	# evaluate forward model @ x
	y = forward_model2(x, λ=oci_band, nPoly=nPoly, nPC=nPC,
						nSIF=nSIF, trans_mat=trans_mat, sza=sza,
						vza=vza, E=E);
	# J w.r.t. x
	J = sum(abs.(ŷ .- y))
	# L = (ŷ .- y)' * inv(Se) * (ŷ .- y) ;
	# R = (x .- xa)' * inv(Sa) * (x .- xa) ;
	# J = L .+ R ;
	return J
end
  ╠═╡ =#

# ╔═╡ e648c84e-2ecb-4fe2-997b-7c96cc4a1940
#=╠═╡
function loss_function(x, p)
	# p is params for forward model and error matrix
	λ     = p.λ;       # wavelength range
	nPoly = p.nPoly;   # degree of polynomials
	nPC   = p.nPC;     # number of eigen vectors used
	nSIF  = p.nSIF;
	trans_mat = p.trans_mat;
	E     = p.E;
	xa    = p.xa;      # priori
	ŷ     = p.ŷ;
	sza   = p.sza;
	vza   = p.vza;
	Se    = p.Se;
	Sa    = p.Sa;

	# evaluate forward model @ x
	y = forward_model2(x, λ=oci_band, nPoly=nPoly, nPC=nPC,
						nSIF=nSIF, trans_mat=trans_mat, sza=sza,
						vza=vza, E=E);
	# J w.r.t. x
	L = (ŷ .- y)' * inv(Se) * (ŷ .- y) ;
	R = (x .- xa)' * inv(Sa) * (x .- xa) ;
	J = L .+ R ;
	return J
end
  ╠═╡ =#

# ╔═╡ 0aad0a27-d51e-4da5-b11b-d0c04859af73
#=╠═╡
begin
	# define a problem
	param = (
		λ     = oci_band,
		nPoly = n,
		nPC   = nPC,
		nSIF  = nSIF,
		trans_mat = HighResSVD.PrinComp[:, 1:nPC],
		E     = E,
		xa    = ma_new.x,
		ŷ     = R_TOA,
		sza   = sza,
		vza   = vza,
		Se    = Se,
		Sa    = Sa_new
	)

	x0 = m4_new.x;
	# define non linear prob
	prob = NonlinearLeastSquaresProblem(loss_function, x0, param);
end
  ╠═╡ =#

# ╔═╡ 60a70269-81a9-4f93-9155-f2d769432ddc
#=╠═╡
@time sol_gn = solve(prob, NewtonRaphson(), store_trace = Val(true))
  ╠═╡ =#

# ╔═╡ 9ecaf87a-22c0-45d6-b6d8-93a9bb74e15d
#=╠═╡
@time sol_lm = solve(
					 prob, 
					 LevenbergMarquardt(),
					 store_trace = Val(true),
				)
  ╠═╡ =#

# ╔═╡ 4ab52ce5-f5ea-4f76-a922-228c28a67005
#=╠═╡
u_gn = sol_gn.u; u_lm = sol_lm.u;
  ╠═╡ =#

# ╔═╡ 1945df2f-7f87-4fb6-ad5e-349c3008e4ee
#=╠═╡
begin
	# residual
	K_lm, y_lm = Jacobian(u_lm, x -> forward_model2(x));
	K_gn, y_gn = Jacobian(u_gn, x -> forward_model2(x));

	# evaluate averaging kernel
	G_lm = GainMatrix(K_lm, Se=Se, Sa=Sa_new);
	G_gn = GainMatrix(K_gn, Se=Se, Sa=Sa_new);

	# Averaging kernel
	A_lm = G_lm * K_lm;
	A_gn = G_gn * K_gn;
	
end
  ╠═╡ =#

# ╔═╡ 8740635d-4c4a-4fdf-a487-1ae2b158ff96
#=╠═╡
heatmap(A_lm, size=(450, 400)); title!("G-N")
  ╠═╡ =#

# ╔═╡ 825c5181-e807-43e7-a086-ce7abda4999d
#=╠═╡
tr(A_lm)
  ╠═╡ =#

# ╔═╡ 755931f2-2c51-4cf9-ae12-cba4add9c7be
#=╠═╡
begin
	# Create individual plots
	p1__ = plot(oci_band, K_gn[:,1], title="Jacobian @ the priori", label="refl. const. Jac.")
	plot!(oci_band, K_gn[:,2], title="Jacobian @ the priori", label="refl. 1st Jac.")
	p2__ = plot(oci_band, K_gn[:,n+2], label="PC#1 Jac.")
	p3__ = plot(oci_band, K_gn[:,n+3], label="PC#2 Jac.")
	plot!(p3__, oci_band, K_gn[:,n+4], label="PC#3 Jac.")
	p4__ = plot(oci_band, K_gn[:,n+5], label="PC#4 Jac.")
	p5__ = plot(oci_band, K_gn[:,end-nSIF+1], label="SIF mag. Jac.")
	plot!(p5__, oci_band, K_gn[:,end-nSIF+2], label="SIF λ₀ Jac.")
	plot!(p5__, oci_band, K_gn[:,end-nSIF+3], label="SIF σ Jac.")
	
	# Combine the plots in a 1x4 gri6
	plot(p1, p2, p3, p4, p5, layout=(5,1), size=(600, 600))
	# xlims!(640, 730)
end
  ╠═╡ =#

# ╔═╡ ce85ad78-1c64-4822-84cc-e9f748105145
#=╠═╡
begin
	# vis
	plot(oci_band, R_TOA .- m4_new.y_x, size=(700, 200), label="residual, iter#4, w/ SIF fit")
	plot!(oci_band, R_TOA .- m4.y_x, label="residual, iter#4, w/o SIF fit")
	plot!(oci_band, R_TOA .- y_gn, label="G-N")
	plot!(oci_band, R_TOA .- y_lm, label="L-M")
	title!("Residuals", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 12f0d3f8-a95e-46e4-adfa-46e41746a284
#=╠═╡
begin
	rho_gn = hcat(v1...)' * u_gn[1 : n+1];
	rho_lm = hcat(v1...)' * u_lm[1 : n+1];
	
	plot(oci_band, rho0_new, label="a priori", size=(1000, 300))
	plot!(oci_band, rho1_new, label="iter#1 ρₜ")
	plot!(oci_band, rho4_new, label="iter#4 ρₜ")

	plot!(oci_band, rho_gn, label="G-N")
	plot!(oci_band, rho_lm, label="L-M")
	
	title!("Total reflectance", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 96f6e86b-acd2-4e13-b5f6-de049dd6b43d
#=╠═╡
begin
	plot(oci_band, ma_new.x[end-nSIF+1] .* SIF_shape(oci_band, λ₀=ma_new.x[end-nSIF+2], σ=ma_new.x[end-nSIF+3]), label="initial guess (nFLH)", size=(1000, 300))
	plot!(oci_band, m1_new.x[end-nSIF+1] .* SIF_shape(oci_band, λ₀=m1_new.x[end-nSIF+2], σ=m1_new.x[end-nSIF+3]), label="iter#1")
	plot!(oci_band, m2_new.x[end-nSIF+1] .* SIF_shape(oci_band, λ₀=m2_new.x[end-nSIF+2], σ=m2_new.x[end-nSIF+3]), label="iter#2")
	plot!(oci_band, m3_new.x[end-nSIF+1] .* SIF_shape(oci_band, λ₀=m3_new.x[end-nSIF+2], σ=m3_new.x[end-nSIF+3]), label="iter#3")
	plot!(oci_band, m4_new.x[end-nSIF+1] .* SIF_shape(oci_band, λ₀=m4_new.x[end-nSIF+2], σ=m4_new.x[end-nSIF+3]), label="iter#4")

	plot!(oci_band, u_gn[end-nSIF+1] .* SIF_shape(oci_band, λ₀=u_gn[end-nSIF+2], σ=u_gn[end-nSIF+3]), label="GaussianNewton\n (Take the forth iteration as x0)")

	plot!(oci_band, u_lm[end-nSIF+1] .* SIF_shape(oci_band, λ₀=u_lm[end-nSIF+2], σ=u_lm[end-nSIF+3]), label="Levenberg-Marquardt")
	
	title!("retrieved SIF")
end
  ╠═╡ =#

# ╔═╡ 648b4340-c819-44dd-b666-714db5e5c62a
#=╠═╡
begin
	T_gn = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * u_gn[(n+2):(n+nPC+1)]);
	T_lm = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * u_lm[(n+2):(n+nPC+1)]);
	
	plot(oci_band, T0_new, size=(1000, 300), label="inital guess", dpi=300)
	plot!(oci_band, T1_new, label="iter#1")
	plot!(oci_band, T4_new, label="iter#4")
	plot!(oci_band, T_gn, label="G-N")
	plot!(oci_band, T_lm, label="L-M")
	title!("transmittance")
end
  ╠═╡ =#

# ╔═╡ 2c3a3310-7355-4087-a566-271768b33bd6
#=╠═╡
begin
	# vis. reduction in loss function
	fnorms_gn = [h.fnorm for h in sol_gn.trace.history];
	fnorms_lm = [h.fnorm for h in sol_lm.trace.history];
	p1_ = plot(fnorms_gn, yaxis=:log, label="G-N")
	title!("Reduction in Residual (L2-norm)")
	p2_ = plot(fnorms_lm, label="L-M")
	plot(p1_, p2_, layout=(2,1), size=(600, 300))
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─857eaa38-cc95-42d7-82f6-853ffa39dfe6
# ╠═99ff878a-6e71-11f0-17ed-b7ac188e90b8
# ╠═04c805e7-45b5-4878-b288-0cf1d02d31fc
# ╠═0b112f15-6cc7-4f02-849e-e0ef8a71b639
# ╠═922ddadd-a129-406d-9de6-892899786e73
# ╠═3a475f5d-b7f2-4dad-9335-82b6bf6e368b
# ╠═0ec3629f-0278-42b1-8ab8-f399d4d4f216
# ╟─05837924-482b-4564-a770-3544f736889b
# ╠═379babe3-7d99-431b-b5db-499ee9b5b406
# ╠═3ac1d3eb-a22b-441c-8343-062f1d733779
# ╠═6f24e4fe-94b5-45bd-bf46-a98a0fdbaf48
# ╠═401b62ff-9966-40b7-ac5d-ed5d704ddda3
# ╟─0ccf2db1-9080-4d29-bfc7-11dffa706f62
# ╠═a42bd26f-46d5-44a4-81d8-7788899b95bc
# ╠═cc1acba5-d114-4579-a64f-8546c2df40b1
# ╠═672286a7-5b44-49f3-8098-9371f5928826
# ╟─f80f7a81-000a-4784-9d10-713406303102
# ╟─acacde64-9957-409d-ae67-428d13428e9d
# ╠═0d68673e-5d07-4703-96f6-d1a4ef919f0e
# ╟─063343a5-5879-4cb7-91ad-5068fe0b33d2
# ╠═466c6800-dd8d-4b11-b33b-bff17dfcf387
# ╟─3c5964a2-c1e5-4dff-9932-5db894771191
# ╟─434ee765-087e-456a-9696-2ba87fa3c5f3
# ╟─40253fb3-981f-4c2d-9f43-ce1c802fc6ef
# ╠═9dcc1616-91d6-45d8-9873-2c449b6e321e
# ╠═b2dbaa15-ac40-40cd-b022-d49193febaa9
# ╠═ab74fe0c-cfa8-45fc-b4fd-8fea3f93c51b
# ╠═e59a4998-c578-42c3-b4e8-61585544f69b
# ╠═3948942d-e754-445a-aa4f-e7dc79537822
# ╠═c4d3782c-f85d-492e-a805-61d6f98fb657
# ╠═407cf364-32c6-4d9f-9596-6a04bbd5a588
# ╠═d5cfaed6-0063-4649-83da-a64727487741
# ╟─d0cbf663-ac73-413a-951c-f99bf8d2cd8d
# ╠═dd2cd8cb-ed9e-4b6b-af99-59fe26809d39
# ╠═7dcb675f-fd35-46ed-ba58-82df3d68627b
# ╠═3923d033-4639-43a3-a693-8d77c04dd186
# ╠═a512b192-5d5e-4688-8b84-f0bc27aa18e7
# ╟─93c48028-a4bb-4d6d-9bc4-85749a675793
# ╠═bdcc5bf7-7ab0-43a2-8710-09b4b4366b1a
# ╠═b621fa58-9f13-48a2-9144-b3a3cb5292ac
# ╟─e9bc8ce0-14a1-4cbe-9df0-c5b5098ecede
# ╟─7b0a281d-daaa-4aaa-a001-12be469225f9
# ╠═c17a958d-fec3-445a-ba1f-59f65ad63af6
# ╟─97672495-e0b1-4952-9f84-a26e926c7235
# ╟─3d80255b-7409-4d8b-9fb7-b05ed286b18a
# ╠═0d03bf4e-64d7-4f52-b52c-8ad17a93157c
# ╠═aed98d25-2b7a-4755-bb5a-3acbd1bae0a4
# ╠═f8320953-b1b2-4954-8215-2fa6f27cb87e
# ╟─556e3e8b-aae5-4462-9aab-1f5c3f90c5a4
# ╠═d8810003-05c7-495f-b4d1-77a057698d2e
# ╠═c9d962f6-8722-4faa-b34f-092de7a76bcf
# ╠═abb9b4e8-9c9c-4d82-8190-06ededcbfd52
# ╠═2a4b61f9-328a-4e92-ae84-58bdda55dc74
# ╠═33a4a5b0-ae07-4536-9c45-a2043d136f9f
# ╠═93b65f52-c5a5-4580-a64b-5a50a44208af
# ╟─0ff68f72-cc18-415b-a7b9-b94d49ee74dd
# ╟─0bf97b73-04c9-4eb5-906a-23827a2c5f3a
# ╠═98bbf74e-6d47-4f25-b060-3f3c6d289a1a
# ╠═e648c84e-2ecb-4fe2-997b-7c96cc4a1940
# ╠═0aad0a27-d51e-4da5-b11b-d0c04859af73
# ╠═60a70269-81a9-4f93-9155-f2d769432ddc
# ╠═9ecaf87a-22c0-45d6-b6d8-93a9bb74e15d
# ╠═4ab52ce5-f5ea-4f76-a922-228c28a67005
# ╠═1945df2f-7f87-4fb6-ad5e-349c3008e4ee
# ╠═8740635d-4c4a-4fdf-a487-1ae2b158ff96
# ╠═825c5181-e807-43e7-a086-ce7abda4999d
# ╠═755931f2-2c51-4cf9-ae12-cba4add9c7be
# ╟─ce85ad78-1c64-4822-84cc-e9f748105145
# ╠═12f0d3f8-a95e-46e4-adfa-46e41746a284
# ╟─96f6e86b-acd2-4e13-b5f6-de049dd6b43d
# ╟─648b4340-c819-44dd-b666-714db5e5c62a
# ╠═2c3a3310-7355-4087-a566-271768b33bd6
