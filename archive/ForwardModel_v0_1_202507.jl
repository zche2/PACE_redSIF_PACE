### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 1240caa6-6d5a-11f0-20c1-3942248725aa
import Pkg; Pkg.activate("..");

# ╔═╡ ce694051-4175-4866-b763-10c971239fed
using Polynomials, ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles


# ╔═╡ f91e476e-229e-4193-b2d9-f519d5987dda
using NCDatasets

# ╔═╡ 8ce45896-b10e-47bd-b95b-1c0efd596460
using Statistics

# ╔═╡ 624f5eb1-50fe-4911-b1cb-66239fbafb26
include("../PACE_SIF.jl")

# ╔═╡ 81024f0c-48c3-44ac-bff5-1ccfbbb29188
md"""
> #### Read and calculate SNR
"""

# ╔═╡ f0dc3dc7-a718-4f00-8e82-f451acdeccc2
begin
	filename = raw"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt";
	lines = readlines(filename);
	end_header_index = findfirst(x -> x == "/end_header", lines);
	data = readdlm(filename, String, skipstart=end_header_index);
	# row 2: wavelength, row 4: c1, row 5: c2
	# println(data)
end

# ╔═╡ b85c4fe0-07eb-4b80-bf4a-9a7f5aa161f8
# ╠═╡ disabled = true
#=╠═╡
# function SNR_LUT(
# 	rad::Vector{FT};    # TOA radiance (corresponds to measurements)
# 	λ_min::FT = 314.,          # minimum wavelength
# 	λ_max::FT = 2259.,   		# maximum wavelength
# 	FileDir::String = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt",
# 	band::String = "Red"
# 	) where {FT <: AbstractFloat}

# 	# return signal to noise ratio (SNR) and noise (L * SNR)

# 	# ---- step 1: read files ----
# 	lines = readlines(FileDir)
# 	lines = readlines(filename);
# 	end_header_index = findfirst(x -> x == "/end_header", lines);
# 	data  = readdlm(filename, String, skipstart=end_header_index);
# 	FPA   = data[:, 1];  # 1st column: band
# 	wvlen = parse.(Float64, data[:, 2]);  # 2nd column: center wavelength
# 	c1    = parse.(Float64, data[:, 4]);  # 4th column: c1
# 	c2    = parse.(Float64, data[:, 5]);  # 5th column: c2

# 	# ---- step 2: select wavelength that matches radiance ----
# 	ind = findall(( λ_min .< wvlen .< λ_max) .& (FPA .== band));
# 	println(ind)
# 	nse = sqrt.( c1[ind] .+ c2[ind] .* rad);
# 	SNR = rad ./ nse;

# 	return SNR, nse
# end
  ╠═╡ =#

# ╔═╡ 83972c75-98a7-4909-844d-3fa40ab43dc4
function SNR_LUT(
	rad,
	c1,
	c2
	)

	# return signal to noise ratio (SNR) and noise (L * SNR)
	nse = sqrt.( c1 .+ c2 .* rad);
	SNR = rad ./ nse;

	return SNR, nse
end

# ╔═╡ 0a0a382e-b8e8-4b20-b7e1-1b2d42e2c5d1
begin
	FPA   = data[:, 1];                   # 1st column: band
	wvlen = parse.(Float64, data[:, 2]);  # 2nd column: center wavelength
	c1    = parse.(Float64, data[:, 4]);  # 4th column: c1
	c2    = parse.(Float64, data[:, 5]);  # 5th column: c2
end

# ╔═╡ fd61a65b-7f71-43c1-983b-3a374bdbb853
md"""
> ##### L_typ
"""

# ╔═╡ 843b5d9c-4f3e-4d39-a53c-c484eb31064d
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
	λ_min    = 620.;
	λ_max    = 860.;
	ind      = findall( λ_min .< red_band .< λ_max );
	E        = oci["red_solar_irradiance"][ind];
	@show oci_band = red_band[ind];
	rhot     = oci["rhot_red"][pixel, scan, ind];
	R_TOA    = oci["radiance_red"][pixel, scan, ind];
	vza      = oci["sensor_zenith"][pixel, scan];
	sza      = oci["solar_zenith"][pixel, scan];
	println("Data read!")
	@show oci["chlor_a"][pixel, scan];
	@show nflh[pixel, scan];
	# select the band
	
end

# ╔═╡ 7e11ef43-b446-42e6-aff8-30956d339099
begin
	snr_ind = findall((FPA .== "Red") .& (λ_min .< wvlen .< λ_max));
	# check dimensions
	@show size(oci_band)
	@show size(snr_ind)
	# SNR and noise
	snr, nse = SNR_LUT(R_TOA, c1[snr_ind], c2[snr_ind]);
	snr_bl, nse_bl = SNR_LUT(ones(size(oci_band)) .* 12., c1[snr_ind], c2[snr_ind]);
end

# ╔═╡ 11ae2d74-3db5-4877-b1cf-c65b8e6536ab
Se = Diagonal(nse.^2)

# ╔═╡ 00ce12d4-9b19-4273-be1d-acd5f288a747
begin
	# plot(red_band, E)
	plot(oci_band, R_TOA, size=(500, 200), label="TOA radiance observation")
	xlabel!("wv [nm]")
	ylabel!("TOA Rad [W/m2/µm/sr]")
end

# ╔═╡ 5a50a92f-3719-45cd-a8de-c86222372e93
begin
	plot(oci_band, snr, size=(500, 200), label="TOA radiance observation")
	plot!(oci_band, snr_bl, label="uniform 12 W/m2/µm/sr")
	xlabel!("wv [nm]")
	ylabel!("SNR")
end

# ╔═╡ eae34753-17b1-4b67-b898-b145aa2f5b9d
md"""
> ##### Start with a even simpler model.

$$R_{TOA}=E(\lambda)cos(SZA)\rho_s(\lambda)$$

where $R_{TOA}$ (W/m2/µm) is upwelling flux at TOA and $E(\lambda)$ is extraterrestrial solar irradiance (W/m2/µm).

where $\rho_s$ is polynomial function (or Landre? funciton)
"""

# ╔═╡ 2c965979-8308-4f7b-9f8b-a17737dcba4e
begin
	n  = length(oci_band);
	λm = mean(oci_band);
	λc = oci_band .- λm;   # centered wavelength.
	# y_prime is scaled TOA radiance
	y_prime = R_TOA .* pi
	# define Jacobian matrix
	K = E .* cosd(sza) .* [ones(n) λc λc.^2 λc.^3];
	# suppose it is a linear problem, Sa=0
	G = inv( K'inv(Se)K )K'inv(Se);
	# retrieve
	x̂ = G * y_prime;
end

# ╔═╡ 37688d35-2c40-46b9-9bf1-5bb54aaf2b5a
begin
	# reconstruct
	ŷ = K * x̂;
	# averaging kernel
	A = G * K;
end

# ╔═╡ 1976224e-38d4-4b4f-aad9-976787a2a051
# posterior error
Ŝ = inv(K'inv(Se)K)

# ╔═╡ e6c26de2-d1c6-418a-830e-af9f4eeb6343
begin
	plot(oci_band, R_TOA, size=(500, 200), label="obs.")
	plot!(oci_band, ŷ ./ pi, label="poly.fit")
	title!("Retrieved and observed TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ cf02f0e4-f579-4281-84de-d6df73584432
begin
	plot(oci_band,
		R_TOA ./ (E .* cosd(sza) ) .* pi,
		size=(500, 200),
		label="poly.fit")
	plot!(oci_band, ŷ ./ (E .* cosd(sza) ), label="ρ_t")
	title!("Reflectance normalized by π", titlefont=10)
end

# ╔═╡ a4c7ab11-27bb-4f04-8d83-f730597157d8
md"""
> Averaging kernel represents amount of derived state that comes from radiance (truth):
$$x̂ = (I-A)xₐ + Ax$$

The trace is degree of freedoms?
"""

# ╔═╡ f89e0853-c19d-4247-964e-51833e1b6a2a
begin
	heatmap(A, size=(300, 300))
	title!("Averaging kernel?")
end

# ╔═╡ fcd6ec40-afa7-4a33-a07b-0461ba5cc735
begin
	# mean residuals
	plot(oci_band, ŷ ./ pi .- R_TOA, size=(500, 200), label="obs.")
	title!("retrieval - obs [W/m2/µm/sr]", titlefont=10)
end

# ╔═╡ 543a8d64-afbd-464c-ac47-5abd65b377ae
begin
	# mean residuals
	plot(oci_band,( ŷ ./ pi .- R_TOA ) ./ R_TOA, size=(500, 200), label="obs.")
	title!("(retrieval - obs) / obs", titlefont=10)
	ylabel!("%")
end

# ╔═╡ c8eac9b1-7aea-4516-ae2f-f74caeb8e270
md"""
> ##### simple forward model
add PC1 (still linear)

$$R_{TOA}=E(\lambda)cos(SZA)\rho_s(\lambda)T_2(\lambda)$$
"""

# ╔═╡ ea24ed55-856c-4c40-a4d3-3bf6d05fac4c


# ╔═╡ 477c9072-3c97-4b08-93fc-4edb8f3d188c
md"""
> ##### Add transmittance - fine structures
$$R_{TOA}=E(\lambda)cos(SZA)\rho_s(\lambda)T_2(\lambda)$$
"""

# ╔═╡ 968a26be-c38a-4ad2-9e7e-45755745c196


# ╔═╡ c57c7d56-f418-4c4e-b29f-53e718fb7e34
begin
	# use file `transmittance_winter_FineWvResModel`
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

# ╔═╡ 50265f0b-28d7-427f-aa98-928f72dcf632
begin
	# get PCs
	HighResSVD = PACE_SIF.Spectral_SVD(trans, bands, λ_min=λ_min, λ_max=λ_max)
end

# ╔═╡ e0817a08-18cc-4e1f-9386-8c92990a7a90
begin
	@show deno = maximum(abs.(HighResSVD.PrinComp[:, 1]));
	T2 = HighResSVD.PrinComp[:, 1] / deno;
	plot(oci_band, -T2, size=(500, 200))
	title!("PC1 normalized")
end

# ╔═╡ 7dace1ae-7d77-4f48-b6c2-8ad32faa63a3
begin
	K_ln = E .* cosd(sza) .* T2 .* [ones(n) λc λc.^2 λc.^3];
	G_ln = inv( K_ln'inv(Se)K_ln )K_ln'inv(Se);
	x̂_ln = G_ln * y_prime;
	
	# reconstruct
	ŷ_ln = K_ln * x̂_ln;
	
	plot(oci_band, R_TOA, size=(500, 200), label="obs.")
	plot!(oci_band, ŷ_ln ./ pi, label="poly.fit")
	title!("Retrieved and observed TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ 61b1d592-5758-4d5c-9017-8425d6c06667
begin
	# mean residuals
	plot(oci_band,( ŷ_ln ./ pi .- R_TOA ) ./ R_TOA, size=(500, 200), label="mean residual")
	title!("(retrieval - obs) / obs", titlefont=10)
	ylabel!("%")
end

# ╔═╡ eb91e67b-2e98-48db-8964-369be673907e


# ╔═╡ 9378369f-c840-4aae-8e4c-6780f47f5755


# ╔═╡ e499f423-82d0-407b-a8da-6ae7fd642bf6
function construct_T2(
		trans_mat,
		x,
		sza=sza,
		vza=vza
	)
	T_up = trans_mat * x;
	svc  = (secd(sza) + secd(vza)) / secd(vza); 
	T₂   = exp.( svc .* log.(T_up) );   
	max_trans = maximum(abs.(T₂));
	T₂_norm   = T₂ ./ max_trans;
	return T₂_norm
end

# ╔═╡ b0fd9331-65d9-470d-a199-d4d1a632ecb9
# ╠═╡ disabled = true
#=╠═╡
begin
	function forward_model_v0(
		x;
		λ  = λ,               # fiting window (nm)
		T2 = -T2,
		E  = E,               # extra terrestrial irradiance (W/m2/µm)
		sza = sza,
		vza = vza,
		nPoly::Int     = 3,   # order of polynomial
		)	
		
		# --- surface reflectance ---
		poly  = Polynomial(x[1:nPoly+1]);
		λm    = mean(λ);  # the polynomial is evaluated after centralizing the wavelength
		ρₛ    = poly.(λ .- λm);
		# T2_norm = construct_T2(T2, x[end])
		R_TOA = E .* cosd(sza) .* ρₛ .* T2_norm .* pi
		return R_TOA
	end
		
end
  ╠═╡ =#

# ╔═╡ 534410d8-c02c-4689-bcfa-a8b8ccbeb0fa


# ╔═╡ b35af8fb-4418-4ff5-951f-8f7e3c047978


# ╔═╡ 5ae8d300-677b-40aa-bce0-84ad3bef91e9


# ╔═╡ 2b0fe76c-cae4-4e92-bbd2-7444a6450797


# ╔═╡ 968df5d2-2446-489a-8faf-ffcd9a4f95c1


# ╔═╡ 1cd18440-06fb-49cb-99d1-f499bb1ec237


# ╔═╡ 19c30312-3d43-4320-8423-c75a992e3d11
function normalize_rho(rho)
	max_rho = maximum(abs.(rho));
	return rho ./ max_rho
end

# ╔═╡ 14ffa2ec-d23f-4489-9e8f-f0ab13d175c0
function forward_model(
	x;
	λ  = λ,               # fiting window (nm)
	trans_mat = trans_mat,   
	# SVD.U truncated matrix to represent upward transmittance
	E  = E,               # extra terrestrial irradiance (W/m2/µm)
	nPoly::Int     = 3,   # order of polynomial
	nPC_trans::Int = 2,   # number of PC used
	sza = 45.,        	  # solar zenith angle (˚)
	vza = 30.,            # viewing zenith angle (˚)
	if_log_SVD::Bool = false,
						  # forward model is slightly different if the SVD is conducted in log space
	pi_scale::Bool = true
	)
	#=
	- this function gives the radiance at TOA (unit: W/m2/µm/sr)
	- the column state vector x includs (in order) :
		nPoly: polynomial term
		nPC_trans: one-way transmittance PC
	=#
	
	# --- surface reflectance ---
	poly = Polynomial(x[1:nPoly+1]);
	λm   = mean(λ);  # the polynomial is evaluated after centralizing the wavelength
	ρₛ   = poly.(λ .- λm);

	# --- check the dimension of PC matrix ---
	if size(trans_mat)[2] != nPC_trans
		println("Dimension mismatch: trans_mat ❗️")
	end

	if if_log_SVD
		println("principle components from log-SVD! to be updated 😂")
		return
	end
	
	T_up = trans_mat * x[(nPoly+2):(nPoly+nPC_trans+1)]; # get upward transmittance
	# svc  = (secd(sza) + secd(vza)) / secd(vza);          # solar-viewer correction
	# T₂   = exp.( svc .* log.(T_up) );                    # get two-way transmittance

	# normalize T2 [0-1]
	# println("size of T2", T₂)
	
	# max_trans = maximum(abs.(T₂));
	# println("size of max_trans", size(max_trans))
	# println("maximum transmittance", max_trans)
	# T₂_norm   = T₂ ./ max_trans;
	# rho_norm = normalize_rho(ρₛ);
	T₂_norm  = construct_T2(trans_mat, x[(nPoly+2):(nPoly+nPC_trans+1)]);
	
	
	# TOA upwelling flux
	scale_factor = pi_scale ? pi : 1;
	# F_TOA = E .* cosd(sza) .* ρₛ .* T₂ ./ scale_factor
	R_TOA = E .* cosd(sza) .* ρₛ .* T₂_norm ./ scale_factor
	# R_TOA = E .* cosd(sza) .* abs.(ρₛ) .* T₂_norm ./ scale_factor

	return R_TOA
end

# ╔═╡ 89ee5080-35f5-43fd-951f-19f7f4f0b658
begin
	degree    = 3;
	nPC_trans = 2;
	model_for_jacobian = x -> forward_model(
		x,
		λ=collect(skipmissing(oci_band)),
		trans_mat=HighResSVD.PrinComp[:, 1:nPC_trans],
		E=E,
		nPoly=degree,
		nPC_trans=nPC_trans,
		sza=sza,
		vza=vza,
	)
end

# ╔═╡ 860d4ae3-8426-4c70-aea3-7c33fc2d471f
begin
	Sₐ = zeros(6,6);
	# the polynomial terms are not constraining retrieval at all!
	for i=1:(degree+1)
	    Sₐ[i,i] = 1e2;
	end
	# error in coeffs of PC are proportional to their explained variance? - now just assign a constant
	# HighResSVD.VarExp or (rel_error*HighResSVD.VarExp)^2   
	rel_error   = .001
	for i=(degree+2):(degree+nPC_trans+1)
	    Sₐ[i,i] = rel_error .* HighResSVD.VarExp[i - (degree+1)];
	end
	# error in SIF is determined by the shape
	# Sₐ[end, end] = 1e20;
	@show Sₐ
end

# ╔═╡ 1e86886e-2c60-4039-b657-0930e505eba4
begin
	x₁ = vcat(x̂, -7, .3);
	F₁ = model_for_jacobian(x₁)
	y₁ = R_TOA .- F₁;
	# allocate
	result1 = DiffResults.JacobianResult(zeros(length(oci_band)), x₁);
	@time ForwardDiff.jacobian!(result1, model_for_jacobian, x₁);
	# Jacobian
	K1 = DiffResults.jacobian(result1);
	# Gaining matrix
	G1 = inv( K1'inv(Se)K1 + inv(Sₐ) )K1'inv(Se);     # inv( K1'inv(Se)K1 )K1'inv(Se);
	# 
	@show x̂₁  = x₁ + G1 * y_prime;
end

# ╔═╡ c9d2527e-1f7e-41ad-aabc-c66b1638123a
begin
	plot(oci_band, construct_T2(HighResSVD.PrinComp[:, 1:2], x̂₁[end-1:end]), size=(500, 200))
	plot!(oci_band, construct_T2(HighResSVD.PrinComp[:, 1], x̂₁[end-1]))
	title!("transmittance reconstructed")
end

# ╔═╡ 60a19ce3-d57d-4b88-9f65-244fc9df2e14
begin
	plot(oci_band, model_for_jacobian(x̂₁), size=(500, 200))
	plot!(oci_band, R_TOA, size=(500, 200))
end

# ╔═╡ 518026a7-2743-48b0-b5b6-af588cb94812
begin
	poly1 = Polynomial(x̂₁[1:degree+1]);
	ρs1  = poly1.(λc);
	
	plot(oci_band,
		R_TOA ./ (E .* cosd(sza) ) .* pi,
		size=(500, 200),
		label="obs.")
	plot!(oci_band, ρs1, label="ρ_t")
	title!("Reflectance normalized by π", titlefont=10)
end

# ╔═╡ 98963786-46c9-4b90-955b-8f05a971bed2
md"""
1 iter

"""

# ╔═╡ bf91ec74-4406-4c3f-b773-457abc9f45cb
begin
	F₂ = DiffResults.value(result1);
	# allocate
	result2 = DiffResults.JacobianResult(zeros(length(oci_band)), x̂₁);
	@time ForwardDiff.jacobian!(result2, model_for_jacobian, x̂₁);
	# Jacobian
	K2 = DiffResults.jacobian(result2);
	y₂ = R_TOA .- F₂ .+ K2 * (x̂₁ - x₁);   # R_TOA .- model_for_jacobian(x̂₁);
	# Gaining matrix
	G2 = inv( K2'inv(Se)K2 + inv(Sₐ) )K2'inv(Se);     # inv( K1'inv(Se)K1 )K1'inv(Se);
	# 
	@show x̂₂  = x₁ + G2 * y₂;
end

# ╔═╡ 97ab2d84-1267-4240-9629-8f72f6ae0ee8
begin
	plot(oci_band, model_for_jacobian(x̂₂), size=(500, 200))
	plot!(oci_band, R_TOA, size=(500, 200))
end

# ╔═╡ 28083551-d1bf-4ec7-a5e3-5a50985b717a
begin
	plot(oci_band, construct_T2(HighResSVD.PrinComp[:, 1:2], x̂₂[end-1:end]), size=(500, 200))
	plot!(oci_band, construct_T2(HighResSVD.PrinComp[:, 1], x̂₂[end-1]))
	title!("transmittance reconstructed")
end

# ╔═╡ 00f88fa4-6175-4ce9-92d4-a46360989381
begin
	poly = Polynomial(x̂₂[1:degree+1]);
	ρs2  = poly.(λc);
	
	plot(oci_band,
		R_TOA ./ (E .* cosd(sza) ) .* pi,
		size=(500, 200),
		label="obs.")
	plot!(oci_band, ρs2, label="ρ_t")
	title!("Reflectance normalized by π", titlefont=10)
end

# ╔═╡ b2913214-20db-4069-bfa9-72942054769c
md"""
2 iter
"""

# ╔═╡ 4b04df67-d191-4dc8-a7fd-e10e8774e41c
begin
	F₃ = DiffResults.value(result2);
	# allocate
	result3 = DiffResults.JacobianResult(zeros(length(oci_band)), x̂₂);
	@time ForwardDiff.jacobian!(result3, model_for_jacobian, x̂₂);
	# Jacobian
	K3 = DiffResults.jacobian(result3);
	y₃ = R_TOA .- F₃ .+ K3 * (x̂₂ - x₁);   # R_TOA .- model_for_jacobian(x̂₁);
	# Gaining matrix
	G3 = inv( K3'inv(Se)K3 + inv(Sₐ) )K3'inv(Se);     # inv( K1'inv(Se)K1 )K1'inv(Se);
	# 
	@show x̂₃  = x₁ + G3 * y₃;
end

# ╔═╡ 5c55d2d5-dc76-439b-a59b-d14f6554bee7
begin
	plot(oci_band, model_for_jacobian(x̂₃), size=(500, 200), label="retrieval")
	plot!(oci_band, R_TOA, size=(500, 200), label="obs")
end

# ╔═╡ 860c177b-3acb-4bc4-8337-f2d2f88b5150
begin
	poly3 = Polynomial(x̂₃[1:degree+1]);
	ρs3  = poly3.(λc);
	
	plot(oci_band,
		R_TOA ./ (E .* cosd(sza) ) .* pi,
		size=(500, 200),
		label="obs.")
	plot!(oci_band, ρs3, label="ρ_t")
	title!("Reflectance normalized by π", titlefont=10)
end

# ╔═╡ 7409d914-99f7-43c1-bbc2-b1565426d21a
begin
	plot(oci_band, construct_T2(HighResSVD.PrinComp[:, 1:2], x̂₃[end-1:end]), size=(500, 200), label="2PCs")
	plot!(oci_band, construct_T2(HighResSVD.PrinComp[:, 1], x̂₃[end-1]), label="1PC")
	title!("transmittance reconstructed")
end

# ╔═╡ 5e8eca23-f86e-41da-9f47-3c24ae725b8c
md"""
3 iter
"""

# ╔═╡ cb2c1d1d-3d49-44bb-a9e9-a1b764e98000
begin
	F4 = DiffResults.value(result3);
	# allocate
	result4 = DiffResults.JacobianResult(zeros(length(oci_band)), x̂₃);
	@time ForwardDiff.jacobian!(result4, model_for_jacobian, x̂₃);
	# Jacobian
	K4 = DiffResults.jacobian(result4);
	y₄ = R_TOA .- F4 .+ K4 * (x̂₃ - x₁);   # R_TOA .- model_for_jacobian(x̂₁);
	# Gaining matrix
	G4 = inv( K4'inv(Se)K4 + inv(Sₐ) )K4'inv(Se);     # inv( K1'inv(Se)K1 )K1'inv(Se);
	# 
	@show x̂₄  = x₁ + G4 * y₄;
end

# ╔═╡ ef1f31bc-fbe9-4fb3-b8c2-2742bdfcf444
begin
	plot(oci_band, model_for_jacobian(x̂₄), size=(500, 200), label="retrieval")
	plot!(oci_band, R_TOA, size=(500, 200), label="obs")
end

# ╔═╡ 3725077f-26d3-46e5-8080-056961b51503
begin
	poly4 = Polynomial(x̂₄[1:degree+1]);
	ρs4  = poly4.(λc);
	
	plot(oci_band,
		R_TOA ./ (E .* cosd(sza) ) .* pi,
		size=(500, 200),
		label="obs.")
	plot!(oci_band, ρs4, label="ρ_t")
	title!("Reflectance normalized by π", titlefont=10)
end

# ╔═╡ Cell order:
# ╠═1240caa6-6d5a-11f0-20c1-3942248725aa
# ╠═ce694051-4175-4866-b763-10c971239fed
# ╠═f91e476e-229e-4193-b2d9-f519d5987dda
# ╠═8ce45896-b10e-47bd-b95b-1c0efd596460
# ╟─81024f0c-48c3-44ac-bff5-1ccfbbb29188
# ╠═f0dc3dc7-a718-4f00-8e82-f451acdeccc2
# ╟─b85c4fe0-07eb-4b80-bf4a-9a7f5aa161f8
# ╠═83972c75-98a7-4909-844d-3fa40ab43dc4
# ╠═0a0a382e-b8e8-4b20-b7e1-1b2d42e2c5d1
# ╟─fd61a65b-7f71-43c1-983b-3a374bdbb853
# ╠═843b5d9c-4f3e-4d39-a53c-c484eb31064d
# ╠═7e11ef43-b446-42e6-aff8-30956d339099
# ╠═11ae2d74-3db5-4877-b1cf-c65b8e6536ab
# ╟─00ce12d4-9b19-4273-be1d-acd5f288a747
# ╠═5a50a92f-3719-45cd-a8de-c86222372e93
# ╟─eae34753-17b1-4b67-b898-b145aa2f5b9d
# ╠═2c965979-8308-4f7b-9f8b-a17737dcba4e
# ╠═37688d35-2c40-46b9-9bf1-5bb54aaf2b5a
# ╠═1976224e-38d4-4b4f-aad9-976787a2a051
# ╠═e6c26de2-d1c6-418a-830e-af9f4eeb6343
# ╠═cf02f0e4-f579-4281-84de-d6df73584432
# ╟─a4c7ab11-27bb-4f04-8d83-f730597157d8
# ╟─f89e0853-c19d-4247-964e-51833e1b6a2a
# ╠═fcd6ec40-afa7-4a33-a07b-0461ba5cc735
# ╠═543a8d64-afbd-464c-ac47-5abd65b377ae
# ╟─c8eac9b1-7aea-4516-ae2f-f74caeb8e270
# ╠═7dace1ae-7d77-4f48-b6c2-8ad32faa63a3
# ╠═61b1d592-5758-4d5c-9017-8425d6c06667
# ╠═ea24ed55-856c-4c40-a4d3-3bf6d05fac4c
# ╠═477c9072-3c97-4b08-93fc-4edb8f3d188c
# ╠═968a26be-c38a-4ad2-9e7e-45755745c196
# ╠═624f5eb1-50fe-4911-b1cb-66239fbafb26
# ╠═c57c7d56-f418-4c4e-b29f-53e718fb7e34
# ╠═50265f0b-28d7-427f-aa98-928f72dcf632
# ╠═e0817a08-18cc-4e1f-9386-8c92990a7a90
# ╠═eb91e67b-2e98-48db-8964-369be673907e
# ╠═9378369f-c840-4aae-8e4c-6780f47f5755
# ╠═e499f423-82d0-407b-a8da-6ae7fd642bf6
# ╠═b0fd9331-65d9-470d-a199-d4d1a632ecb9
# ╠═534410d8-c02c-4689-bcfa-a8b8ccbeb0fa
# ╠═b35af8fb-4418-4ff5-951f-8f7e3c047978
# ╠═5ae8d300-677b-40aa-bce0-84ad3bef91e9
# ╠═2b0fe76c-cae4-4e92-bbd2-7444a6450797
# ╠═968df5d2-2446-489a-8faf-ffcd9a4f95c1
# ╠═1cd18440-06fb-49cb-99d1-f499bb1ec237
# ╠═19c30312-3d43-4320-8423-c75a992e3d11
# ╠═14ffa2ec-d23f-4489-9e8f-f0ab13d175c0
# ╠═89ee5080-35f5-43fd-951f-19f7f4f0b658
# ╠═860d4ae3-8426-4c70-aea3-7c33fc2d471f
# ╠═1e86886e-2c60-4039-b657-0930e505eba4
# ╠═c9d2527e-1f7e-41ad-aabc-c66b1638123a
# ╠═60a19ce3-d57d-4b88-9f65-244fc9df2e14
# ╠═518026a7-2743-48b0-b5b6-af588cb94812
# ╟─98963786-46c9-4b90-955b-8f05a971bed2
# ╠═bf91ec74-4406-4c3f-b773-457abc9f45cb
# ╠═97ab2d84-1267-4240-9629-8f72f6ae0ee8
# ╠═28083551-d1bf-4ec7-a5e3-5a50985b717a
# ╠═00f88fa4-6175-4ce9-92d4-a46360989381
# ╟─b2913214-20db-4069-bfa9-72942054769c
# ╠═4b04df67-d191-4dc8-a7fd-e10e8774e41c
# ╠═5c55d2d5-dc76-439b-a59b-d14f6554bee7
# ╠═860c177b-3acb-4bc4-8337-f2d2f88b5150
# ╟─7409d914-99f7-43c1-bbc2-b1565426d21a
# ╟─5e8eca23-f86e-41da-9f47-3c24ae725b8c
# ╠═cb2c1d1d-3d49-44bb-a9e9-a1b764e98000
# ╠═ef1f31bc-fbe9-4fb3-b8c2-2742bdfcf444
# ╠═3725077f-26d3-46e5-8080-056961b51503
