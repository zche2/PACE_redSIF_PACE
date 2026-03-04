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
	
end

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

# ╔═╡ acacde64-9957-409d-ae67-428d13428e9d
#=╠═╡
begin
	# the PCs look like:
	plot(oci_band, HighResSVD.PrinComp[:,1:4], size=(500, 200))
	title!("eigen vectors")
end
  ╠═╡ =#

# ╔═╡ 0d68673e-5d07-4703-96f6-d1a4ef919f0e
#=╠═╡
findall(coalesce.(nflh .> .6, false))
  ╠═╡ =#

# ╔═╡ 063343a5-5879-4cb7-91ad-5068fe0b33d2
md"""
> ##### SNR -> measurement covariance matrix $S_{\epsilon}$
"""

# ╔═╡ 466c6800-dd8d-4b11-b33b-bff17dfcf387
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ f80f7a81-000a-4784-9d10-713406303102
#=╠═╡
plot(oci_band, R_TOA, size=(500, 200), label="obs"); ylabel!("W/m2/µm/sr")
  ╠═╡ =#

# ╔═╡ 434ee765-087e-456a-9696-2ba87fa3c5f3
md"""
> ##### Start with polynomial fit +Transmittance
$$\rho_{s}=\sum{a_jP_j}$$
$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)T(\lambda)}{\pi}$$
 $T(\lambda)$ is set to have a maximum of 1.
"""

# ╔═╡ 810f2417-3f25-4fab-a88f-a99642a1e2c6
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

# ╔═╡ e59a4998-c578-42c3-b4e8-61585544f69b
#=╠═╡
begin
	# number of polynormial terms and PCs used
	n = 5; nPC = 25;
	# inital guess 
	λm    = mean(oci_band);
	range = λ_max - λ_min;
	λc    = (oci_band .- λm) ./ range;
	K     = E .* cosd(sza) ./ pi .* hcat(collectPl.(λc, lmax=n)...)';
	G     = inv( K'inv(Se)K )K'inv(Se);
	x̂     = G * R_TOA;
	ŷ     = K * x̂;
end
  ╠═╡ =#

# ╔═╡ c4d3782c-f85d-492e-a805-61d6f98fb657
#=╠═╡
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
	λc    = center_wavelength(λ);
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

  ╠═╡ =#

# ╔═╡ d5cfaed6-0063-4649-83da-a64727487741
#=╠═╡
begin
	tmp  = zeros(nPC-2) .+ .001;
	xa   = [x̂... -6. .05 tmp...]';
	rad  = forward_model1(xa, nPoly=n);
	plot(oci_band, R_TOA, size=(500, 200), label="obs.")
	plot!(oci_band, rad, label="initial guess, n=$n")
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 3cb579f3-c9c4-48b3-997d-967f4e1df546
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ d0cbf663-ac73-413a-951c-f99bf8d2cd8d
md"""
🟢 defining functions to calculate Jacobian, gain matrix and do the iteration
"""

# ╔═╡ dd2cd8cb-ed9e-4b6b-af99-59fe26809d39
#=╠═╡
function Jacobian(x, model; len=length(oci_band))
	res = DiffResults.JacobianResult(zeros(len), x);
	ForwardDiff.jacobian!(res, model, x);
	K   = DiffResults.jacobian(res);
	val = DiffResults.value(res);
	return K, val
end
  ╠═╡ =#

# ╔═╡ 7dcb675f-fd35-46ed-ba58-82df3d68627b
#=╠═╡
function GainMatrix(K; Se=Se, Sa=Sa)
	return inv( K'inv(Se)K + inv(Sa) )K'inv(Se)
end
  ╠═╡ =#

# ╔═╡ 3923d033-4639-43a3-a693-8d77c04dd186
@with_kw struct retrieval
    x    # state vector
	y_x  # value evaluated at x
	K    # Jacobian
	G = GainMatrix(K)         # Gain matrix
	A = G*K                   # averaging kernel
end;

# ╔═╡ a512b192-5d5e-4688-8b84-f0bc27aa18e7
#=╠═╡
function iter(
		m ::retrieval,   # retrieval of last iteration
		xa,              # priori
	    rad; 
		Sa    = Sa,      # measurements
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
		G   = GainMatrix(K_n1, Sa=Sa)
	)
	return m_new
end
  ╠═╡ =#

# ╔═╡ 93c48028-a4bb-4d6d-9bc4-85749a675793
md"""
> ##### first iter.
"""

# ╔═╡ bdcc5bf7-7ab0-43a2-8710-09b4b4366b1a
#=╠═╡
begin
	# start from xa
	Ka, ya = Jacobian(xa, x -> forward_model1(x))
	ma = retrieval(x=xa, y_x=ya, K=Ka)
end
  ╠═╡ =#

# ╔═╡ b621fa58-9f13-48a2-9144-b3a3cb5292ac
#=╠═╡
begin
	# 1st iteration
	m1 = iter(ma, xa, R_TOA, Sa=Sa);
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=4, linestyle=:dash, color=:black)
	plot!(oci_band, ma.y_x, label="initial guess, n=$n", linewidth=2)
	plot!(oci_band, m1.y_x, label="iter#1, n=$n", linewidth=2)
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 7b0a281d-daaa-4aaa-a001-12be469225f9
md"""
🟢 recover transmittance
"""

# ╔═╡ c17a958d-fec3-445a-ba1f-59f65ad63af6
#=╠═╡
begin
	T1     = HighResSVD.PrinComp[:, 1:nPC] * m1.x[(n+2):(n+nPC+1)];
	T1_min = minimum(T1);
	T1_max = maximum(T1);
	factor1 = maximum(abs.([T1_min, T1_max]))
	T1_norm = abs.(T1) / factor1
	# @show T1_norm

	plot(oci_band, T1_norm, size=(500, 200), label="iter#1")
	
end
  ╠═╡ =#

# ╔═╡ 97672495-e0b1-4952-9f84-a26e926c7235
md"""
🟣 reflectance
"""

# ╔═╡ 3d80255b-7409-4d8b-9fb7-b05ed286b18a
#=╠═╡
begin
	v1   = collectPl.(λc, lmax=n);
	rho1 = hcat(v1...)' * m1.x[1 : n+1];
	plot(oci_band, ŷ ./ (E .* cosd(sza) ) .* pi,
		label="linear fit ρₜ",
		size=(500, 200))
	plot!(oci_band, rho1, label="iter#1 ρₜ")
	title!("Total reflectance", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 556e3e8b-aae5-4462-9aab-1f5c3f90c5a4
md"""
> ##### 2nd iter.
"""

# ╔═╡ abb9b4e8-9c9c-4d82-8190-06ededcbfd52
#=╠═╡
begin
	# 2nd iteration
	m2 = iter(m1, xa, R_TOA);
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=4, linestyle=:dash, color=:black)
	plot!(oci_band, ma.y_x, label="initial guess, n=$n", linewidth=2)
	plot!(oci_band, m1.y_x, label="iter#1, n=$n", linewidth=2)
	plot!(oci_band, m2.y_x, label="iter#2, n=$n", linewidth=2)
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 2a4b61f9-328a-4e92-ae84-58bdda55dc74
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 33a4a5b0-ae07-4536-9c45-a2043d136f9f
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 93b65f52-c5a5-4580-a64b-5a50a44208af
#=╠═╡
begin
	# resildual
	plot(oci_band, R_TOA .- m1.y_x, size=(500, 200), label="residual, iter#1")
	plot!(oci_band, R_TOA .- m2.y_x, label="residual, iter#2")
	title!("Residuals", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 4872795f-8afd-41fa-abbf-ebf2cba48bb0
md"""
> ##### More iterations
"""

# ╔═╡ d88a1ffb-04ac-41ab-bcb4-a039e9516f03
#=╠═╡
m3 = iter(m2, xa, R_TOA);
  ╠═╡ =#

# ╔═╡ 79a5ff51-b649-4d56-80e4-95ec9470fced
#=╠═╡
m4 = iter(m3, xa, R_TOA);
  ╠═╡ =#

# ╔═╡ 5587f936-e9c0-4b50-b2cb-f1dc6ca93eba
md"""
> ##### Add SIF with designated shape and peak wavelength
$$\rho_{s}=\sum{a_jP_j}$$
$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)T_{\downarrow\uparrow}(\lambda)}{\pi} + SIF(\lambda)T_{\uparrow}(\lambda)$$
 $T_{\uparrow}(\lambda)$ is set to have a maximum of 1, and relates with $T_{\downarrow\uparrow}(\lambda)$ by:

$$T_{\downarrow\uparrow}(\lambda)=exp(SVF \times ln(T_{\uparrow}(\lambda)))$$
Where SVF is some solar/viewing zenith angle correction factor:

$$SVF=\frac{Sec(SZA)+Sec(VZA)}{Sec(VZA)}$$

Assuming a Gaussian shape of SIF emisison:

$$SIF(\lambda)=Aexp(-\frac{(\lambda-\lambda_0)^2}{2\sigma^2})$$
Where $\lambda_0=683$ is the peak wavelength and $\sigma=9$ (to be tuned).

🔴 Q: Shall I give A more constraint? e.g., non-negative? - probably not, it is possible to have negative nFLH detected during some cyanbacteria bloom.
"""

# ╔═╡ 42cd6dc2-4c8a-4862-b45d-242f51ae9bfb
md"""
> ##### Check one-way / two-way transmittance
"""

# ╔═╡ eaacce7b-c6ad-4e8f-9768-74c7c88fec1a
function two_way_trans(T, sza, vza)
	svf = (secd(sza)+secd(vza)) / secd(vza);
	T2  = exp.( svf .* log.(T));
	return T2
end

# ╔═╡ 0474ac88-6908-4337-9039-277751c6bd75
function SIF_shape(λ; λ₀=683., σ=5.0)
	return exp.( - ( λ .- λ₀ ).^2 ./ ( 2 * σ^2 ) )
end

# ╔═╡ b64fdb6c-e4d3-4ca8-aea9-02618bcecd02
#=╠═╡
begin
	# one-way to two-way transmittance
	plot(oci_band, T1_norm, size=(500, 200), label="T↑")
	plot!(oci_band, two_way_trans(T1_norm, sza, vza), label="T↓↑")
	title!("transmittance")
end
  ╠═╡ =#

# ╔═╡ e8040280-82ea-4af2-b8dd-1653deb069a5
#=╠═╡
begin
	# plot spectra of SIF
	plot(oci_band, SIF_shape(oci_band), size=(500, 100), label="SIF (normalzied)" )
	xlims!(645, 720)
end
  ╠═╡ =#

# ╔═╡ 22a31224-df90-4484-807b-a3f2d36a178d
md"""
🟢 defining new forward model
"""

# ╔═╡ 5ad3d348-3997-4e11-8221-602f6c6d3b3e
#=╠═╡
begin
	# define priori error matrix
	nSIF    = 3;
	# priori cov
	Sa_new  = zeros(n+nPC+nSIF+1, n+nPC+nSIF+1);
	# uodate diagonal term
	for i=1:(n+1)
	    Sa_new[i,i] = 1e20;     
		# large variance applies no constrain to these polynomial term
	end
	for i=(n+2):(n+nPC+1)
	    Sa_new[i,i] = rel_error .* HighResSVD.VarExp[i - (n+1)];
	end
	# SIF uncertainty
	Sa_new[end-2, end-2] = 1.5;
	Sa_new[end-1, end-1] = 5.;
	Sa_new[end, end] = 10;
	Sa_new
end
  ╠═╡ =#

# ╔═╡ 5aee9d76-54b1-4143-b0df-b6005b1d14e7
#=╠═╡
function forward_model2(
			x; 
			λ = oci_band,     # wavelength range
			nPoly::Int = n,   # degree of polynomials
			nPC::Int   = nPC,   # number of eigen vectors used
			nSIF::Int  = nSIF,
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
	# upward transmittance
	T     = trans_mat * x[(nPoly+2):(nPoly+nPC+1)];
	T_min = minimum(T);
	T_max = maximum(T);
	factor  = maximum(abs.([T_min, T_max]));
	T1_norm = abs.(T) / factor;
	# downward transmittance (no need to normalize)
	T2_norm = two_way_trans(T1_norm, sza, vza);
	# water-leaving SIF
	SIF_w   = x[end-nSIF+1] .* SIF_shape(λ, λ₀=x[end-nSIF+2], σ=x[end-nSIF+3]);
	# TOA radiance
	rad    = E .* cosd(sza) ./ pi .* T2_norm .* rho + SIF_w .* T1_norm;
	return rad
end

  ╠═╡ =#

# ╔═╡ 31f80fa5-8229-4c4b-bb26-b27bdce65c2f
#=╠═╡
begin
	# set inital guess
	x_for_PCs    = zeros(nPC);
	x_for_PCs[1] = -6.;
	x_for_PCs[2] = .1;
	nflh_val = nflh[pixel, scan]
	# @show xa_new   = [x̂... x_for_PCs... nflh_val 683. 5.]';
	@show xa_new   = [x̂... x_for_PCs... 0. 683. 5.]';
	# xa_new   = [m2.x... nflh_val]';
	
	# wrap up in a struct
	Ka_new, ya_new = Jacobian(xa_new, x -> forward_model2(x));
	ma_new = retrieval(x=xa_new, y_x=ya_new, K=Ka_new, G=GainMatrix(Ka_new,Sa=Sa_new));
	
end
  ╠═╡ =#

# ╔═╡ 261af164-b584-4aca-b064-714a226b79a8
#=╠═╡
nflh_val
  ╠═╡ =#

# ╔═╡ 9f6bc382-7e26-4832-a57e-7d664aa8c322
#=╠═╡
begin
	plot(oci_band, R_TOA, size=(600, 200), label="obs.")
	plot!(oci_band, ya_new, label="initial guess, n=$n, nPC=$nPC, SIF=$nflh_val")
	# plot!(oci_band, ya, label="another initial guess, n=$n, nPC=$nPC, SIF=0")
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 32c99105-fa5f-422a-bf85-6696d983b0c3
md"""
🟢 visualize Jacobian for each term

Notes: 
- PC#3 the O2 absoprtion is most likely to confound with SIF?
- PC#1 and #2, the shapes look so alike (see the PC )

"""

# ╔═╡ fadd7cf1-7a27-4202-94c0-60f09c9657e7
#=╠═╡
begin
	# Create individual plots
	p1 = plot(oci_band, Ka_new[:,1], title="Jacobian @ the priori", label="refl. const. Jac.")
	plot!(oci_band, Ka_new[:,2], title="Jacobian @ the priori", label="refl. 1st Jac.")
	p2 = plot(oci_band, Ka_new[:,n+2], label="PC#1 Jac.")
	p3 = plot(oci_band, Ka_new[:,n+3], label="PC#2 Jac.")
	plot!(p3, oci_band, Ka_new[:,n+4], label="PC#3 Jac.")
	p4 = plot(oci_band, Ka_new[:,n+5], label="PC#4 Jac.")
	p5 = plot(oci_band, Ka_new[:,end-nSIF+1], label="SIF mag. Jac.")
	plot!(p5, oci_band, Ka_new[:,end-nSIF+2], label="SIF λ₀ Jac.")
	plot!(p5, oci_band, Ka_new[:,end-nSIF+3], label="SIF σ Jac.")
	
	# Combine the plots in a 1x4 gri6
	plot(p1, p2, p3, p4, p5, layout=(5,1), size=(600, 600))
	# xlims!(640, 730)
end
  ╠═╡ =#

# ╔═╡ 2068e446-0265-4b91-b2f0-9e5a2d7e5dec
md"""
> ##### Iterations
"""

# ╔═╡ 369724ff-0c8a-4ce6-9943-fae894146392
#=╠═╡
begin
	# 1st iteration
	m1_new = iter(ma_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=3, linestyle=:dash, color=:blue)
	plot!(oci_band, ma_new.y_x, label="initial guess, SIF=$(ma_new.x[end-nSIF+1])", linewidth=2)
	plot!(oci_band, m1_new.y_x, label="iter#1, SIF=$(m1_new.x[end-nSIF+1])", linewidth=2)
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 05caf90e-0895-4664-ae2e-d8ac88cf2b40
#=╠═╡
begin
	# 2nd iteration
	m2_new = iter(m1_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=3, linestyle=:dash, color=:blue)
	plot!(oci_band, ma_new.y_x, label="initial guess, n=$n", linewidth=1)
	plot!(oci_band, m1_new.y_x, label="iter#1, n=$n", linewidth=1)
	plot!(oci_band, m2_new.y_x, label="iter#2, n=$n", linewidth=1)
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 9d99d5a9-be2c-4a6f-bdba-189f35591a78
#=╠═╡
begin
	# 2nd iteration
	m3_new = iter(m2_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	
	plot(oci_band, R_TOA, size=(500, 200), label="obs.", linewidth=3, linestyle=:dash, color=:blue)
	plot!(oci_band, ma_new.y_x, label="initial guess, n=$n", linewidth=1)
	plot!(oci_band, m1_new.y_x, label="iter#1, n=$n", linewidth=1)
	plot!(oci_band, m2_new.y_x, label="iter#2, n=$n", linewidth=1)
	plot!(oci_band, m3_new.y_x, label="iter#3, n=$n", linewidth=1)
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 1e1881c5-4bea-4159-b750-87e941abec5a
#=╠═╡
begin
	# resildual
	plot(oci_band, R_TOA .- m1_new.y_x, size=(700, 200), label="residual, iter#1")
	plot!(oci_band, R_TOA .- m2_new.y_x, label="residual, iter#2")
	plot!(oci_band, R_TOA .- m3_new.y_x, label="residual, iter#3")
	title!("Residuals", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ c2bcb483-fa85-416b-a0fb-15870d62084e
@show "reconstructed SIF peak emission", m1_new.x[end-nSIF+1], m2_new.x[end-nSIF+1], m3_new.x[end-nSIF+1], m4_new.x[end-nSIF+1]

# ╔═╡ fd12c14c-fd48-411e-bdae-00109afcce8b
md"""
🔴 to make it comparable with nFLH globalwise, normalize it with SZA or VZA?
"""

# ╔═╡ 581fd278-129c-4fae-841c-16458a656f52
#=╠═╡
# more iter
m4_new = iter(m3_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
  ╠═╡ =#

# ╔═╡ 61653ce7-35fe-4ee4-93db-46b0918446c2
#=╠═╡
"reconstructed SIF peak emission (normalized by SZA)", m1_new.x[end-nSIF+1] / cosd(sza), m2_new.x[end-nSIF+1] / cosd(sza), m3_new.x[end-nSIF+1] / cosd(sza), m4_new.x[end-nSIF+1] / cosd(sza)
  ╠═╡ =#

# ╔═╡ ea37b0a0-36c4-4e00-894f-b37cf4b51b8a
#=╠═╡
begin
	# resildual
	plot(oci_band, R_TOA .- m4_new.y_x, size=(700, 200), label="residual, iter#4, w/ SIF fit")
	plot!(oci_band, R_TOA .- m4.y_x, label="residual, iter#4, w/o SIF fit")
	title!("Residuals", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 20c699ed-4455-4a70-ae91-9e04a8e9f365
md"""
> ##### Reconstruct each component 
"""

# ╔═╡ cce6b779-bb3a-4307-86be-2fd6bd9f716f
function scale_transmittance(T)
	T_min = minimum(T);
	T_max = maximum(T);
	factor = maximum(abs.([T_min, T_max]))
	return abs.(T) / factor
end

# ╔═╡ e0195a93-c1a8-47c2-b699-63eb9911f87d
#=╠═╡
begin
	rho0_new = hcat(v1...)' * ma_new.x[1 : n+1];
	rho1_new = hcat(v1...)' * m1_new.x[1 : n+1];
	rho2_new = hcat(v1...)' * m2_new.x[1 : n+1];
	rho4_new = hcat(v1...)' * m4_new.x[1 : n+1];
	
	plot(oci_band, rho0_new, label="a priori", size=(500, 150))
	plot!(oci_band, rho1_new, label="iter#1 ρₜ")
	plot!(oci_band, rho2_new, label="iter#2 ρₜ")
	plot!(oci_band, rho4_new, label="iter#4 ρₜ")
	title!("Total reflectance", titlefont=10)
end
  ╠═╡ =#

# ╔═╡ 4e4054fd-d88f-43aa-a043-3ffd503e2278
#=╠═╡
begin
	T0_new = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * ma_new.x[(n+2):(n+nPC+1)]);
	T1_new = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * m1_new.x[(n+2):(n+nPC+1)]);
	T2_new = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * m2_new.x[(n+2):(n+nPC+1)]);
	T3_new = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * m3_new.x[(n+2):(n+nPC+1)]);
	T4_new = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * m4_new.x[(n+2):(n+nPC+1)]);
	
	plot(oci_band, T0_new, size=(600, 200), label="inital guess", dpi=300)
	plot!(oci_band, T1_new, label="iter#1")
	plot!(oci_band, T2_new, label="iter#2")
	plot!(oci_band, T3_new, label="iter#3")
	plot!(oci_band, T4_new, label="iter#4")
	title!("transmittance")
end
  ╠═╡ =#

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

# ╔═╡ 152d1672-0eca-472a-adae-d73604521b98
#=╠═╡
ma_new.x[end-nSIF+1], m4_new.x[end-nSIF+1]
  ╠═╡ =#

# ╔═╡ 0bf97b73-04c9-4eb5-906a-23827a2c5f3a
md"""
> ##### explicitly define loss function
"""

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
@time sol_gn = solve(prob, NewtonRaphson(), store_trace = Val(true))

# ╔═╡ 9ecaf87a-22c0-45d6-b6d8-93a9bb74e15d
@time sol_lm = solve(
					 prob, 
					 LevenbergMarquardt(),
					 store_trace = Val(true),
				)

# ╔═╡ 4ab52ce5-f5ea-4f76-a922-228c28a67005
u_gn = sol_gn.u; u_lm = sol_lm.u;

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
begin
	# vis. reduction in loss function
	fnorms_gn = [h.fnorm for h in sol_gn.trace.history];
	fnorms_lm = [h.fnorm for h in sol_lm.trace.history];
	p1_ = plot(fnorms_gn, yaxis=:log, label="G-N")
	title!("Reduction in Residual (L2-norm)")
	p2_ = plot(fnorms_lm, label="L-M")
	plot(p1_, p2_, layout=(2,1), size=(600, 300))
end

# ╔═╡ 5c03edd1-b980-4ed8-a94d-2541b5f2cfcd
md"""
> ##### Check nFLH retrieval
"""

# ╔═╡ f1fd5c02-0751-4689-856b-a9f2cbad29c8
baseline_wv = [649.599976, 650.900024, 652.099976, 653.299988, 654.599976, 655.799988, 657.099976, 658.299988, 659.599976, 710.500000, 711.799988, 713.000000, 714.299988, 716.799988, 719.200012]

# ╔═╡ 41465a2f-f3ff-4792-8ac1-9dabce29e09b
#=╠═╡
begin
	oci = Dataset(
		"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample_granule_20250501T183011_new.nc");
	pixel  = 471;  # cross-track
	scan   = 899;
	red_band = oci["red_wavelength"][:];
	# red_band = oci["red_bands"][:];
	# cloud    = oci["cloud_flag_dilated"][:, :];
	nflh     = oci["nflh"][:, :];

	ind      = findall( λ_min .< red_band .< λ_max );

	E        = oci["red_solar_irradiance"][ind];
	oci_band = red_band[ind];
	rhot     = oci["rhot_red"][pixel, scan, ind];
	R_TOA    = oci["radiance_red"][pixel, scan, ind];
	vza      = oci["sensor_zenith"][pixel, scan];
	sza      = oci["solar_zenith"][pixel, scan];

	@show nflh[pixel, scan]
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

# ╔═╡ a42bd26f-46d5-44a4-81d8-7788899b95bc
# ╠═╡ disabled = true
#=╠═╡
begin
	oci = Dataset(
		"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample_swath_20250130T20.nc");
	pixel  = 776;  # cross-track
	scan   = 1525;
	# E      = oci["red_solar_irradiance"][:];
	red_band = oci["red_wavelength"][:];
	# cloud    = oci["cloud_flag_dilated"][:, :];
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
  ╠═╡ =#

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

# ╔═╡ Cell order:
# ╠═99ff878a-6e71-11f0-17ed-b7ac188e90b8
# ╠═04c805e7-45b5-4878-b288-0cf1d02d31fc
# ╠═0b112f15-6cc7-4f02-849e-e0ef8a71b639
# ╠═922ddadd-a129-406d-9de6-892899786e73
# ╠═3a475f5d-b7f2-4dad-9335-82b6bf6e368b
# ╠═0ec3629f-0278-42b1-8ab8-f399d4d4f216
# ╟─05837924-482b-4564-a770-3544f736889b
# ╠═379babe3-7d99-431b-b5db-499ee9b5b406
# ╠═6f24e4fe-94b5-45bd-bf46-a98a0fdbaf48
# ╟─acacde64-9957-409d-ae67-428d13428e9d
# ╠═401b62ff-9966-40b7-ac5d-ed5d704ddda3
# ╟─0ccf2db1-9080-4d29-bfc7-11dffa706f62
# ╠═a42bd26f-46d5-44a4-81d8-7788899b95bc
# ╠═41465a2f-f3ff-4792-8ac1-9dabce29e09b
# ╠═0d68673e-5d07-4703-96f6-d1a4ef919f0e
# ╟─063343a5-5879-4cb7-91ad-5068fe0b33d2
# ╠═466c6800-dd8d-4b11-b33b-bff17dfcf387
# ╠═f80f7a81-000a-4784-9d10-713406303102
# ╟─434ee765-087e-456a-9696-2ba87fa3c5f3
# ╠═810f2417-3f25-4fab-a88f-a99642a1e2c6
# ╠═e59a4998-c578-42c3-b4e8-61585544f69b
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
# ╠═2a4b61f9-328a-4e92-ae84-58bdda55dc74
# ╟─33a4a5b0-ae07-4536-9c45-a2043d136f9f
# ╠═93b65f52-c5a5-4580-a64b-5a50a44208af
# ╟─4872795f-8afd-41fa-abbf-ebf2cba48bb0
# ╠═d88a1ffb-04ac-41ab-bcb4-a039e9516f03
# ╠═79a5ff51-b649-4d56-80e4-95ec9470fced
# ╟─5587f936-e9c0-4b50-b2cb-f1dc6ca93eba
# ╟─42cd6dc2-4c8a-4862-b45d-242f51ae9bfb
# ╠═eaacce7b-c6ad-4e8f-9768-74c7c88fec1a
# ╠═0474ac88-6908-4337-9039-277751c6bd75
# ╟─b64fdb6c-e4d3-4ca8-aea9-02618bcecd02
# ╟─e8040280-82ea-4af2-b8dd-1653deb069a5
# ╟─22a31224-df90-4484-807b-a3f2d36a178d
# ╠═5aee9d76-54b1-4143-b0df-b6005b1d14e7
# ╠═5ad3d348-3997-4e11-8221-602f6c6d3b3e
# ╠═31f80fa5-8229-4c4b-bb26-b27bdce65c2f
# ╠═261af164-b584-4aca-b064-714a226b79a8
# ╟─9f6bc382-7e26-4832-a57e-7d664aa8c322
# ╟─32c99105-fa5f-422a-bf85-6696d983b0c3
# ╠═fadd7cf1-7a27-4202-94c0-60f09c9657e7
# ╟─2068e446-0265-4b91-b2f0-9e5a2d7e5dec
# ╠═369724ff-0c8a-4ce6-9943-fae894146392
# ╟─05caf90e-0895-4664-ae2e-d8ac88cf2b40
# ╟─9d99d5a9-be2c-4a6f-bdba-189f35591a78
# ╠═1e1881c5-4bea-4159-b750-87e941abec5a
# ╠═c2bcb483-fa85-416b-a0fb-15870d62084e
# ╟─fd12c14c-fd48-411e-bdae-00109afcce8b
# ╠═61653ce7-35fe-4ee4-93db-46b0918446c2
# ╠═581fd278-129c-4fae-841c-16458a656f52
# ╟─ea37b0a0-36c4-4e00-894f-b37cf4b51b8a
# ╟─20c699ed-4455-4a70-ae91-9e04a8e9f365
# ╠═cce6b779-bb3a-4307-86be-2fd6bd9f716f
# ╟─e0195a93-c1a8-47c2-b699-63eb9911f87d
# ╟─4e4054fd-d88f-43aa-a043-3ffd503e2278
# ╟─0ff68f72-cc18-415b-a7b9-b94d49ee74dd
# ╠═152d1672-0eca-472a-adae-d73604521b98
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
# ╟─12f0d3f8-a95e-46e4-adfa-46e41746a284
# ╟─96f6e86b-acd2-4e13-b5f6-de049dd6b43d
# ╟─648b4340-c819-44dd-b666-714db5e5c62a
# ╠═2c3a3310-7355-4087-a566-271768b33bd6
# ╟─5c03edd1-b980-4ed8-a94d-2541b5f2cfcd
# ╠═f1fd5c02-0751-4689-856b-a9f2cbad29c8
