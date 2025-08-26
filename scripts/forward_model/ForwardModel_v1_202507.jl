### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ bbd4b752-5e7e-11f0-0ac9-43833ecaae94
import Pkg; Pkg.activate("/FraLab/PACE_redSIF_PACE");

# ╔═╡ 9f7e209e-454c-4cc6-9d5c-838c373620c1
using Polynomials, ForwardDiff, DiffResults, Plots, LinearAlgebra

# ╔═╡ 4a17175f-d424-4f09-95b9-407f2cbd74f1
using Statistics

# ╔═╡ 7a8de7e1-a009-4cff-8520-f93b7efcecf9
using NCDatasets

# ╔═╡ 78653793-ece0-49b3-ace8-24ac0a84a695
include("../../PACE_SIF.jl")

# ╔═╡ 3452c8b6-884e-41d0-a52a-a033545ed5f5
md"""

### Constructing a Forward Model
See the slide for the formula of the forward model

$$R_{TOA}=E(\lambda)cos(SZA)\rho_s(\lambda)T_2(\lambda) + \pi SIF(\lambda)T_{\uparrow}(\lambda)$$

where $R_{TOA}$ (W/m2/µm) is upwelling flux at TOA and $E(\lambda)$ is extraterrestrial solar irradiance (W/m2/µm).

"""

# ╔═╡ 8e1a0685-bcec-45ac-aab9-87cca07361bc
# function forward_model(
# 	x::AbstractArray{FT};
# 	λ::Vector{FT}  = λ,    # fiting window (nm)
# 	trans_mat::Matrix{FT} = trans_mat,   # SVD.U truncated matrix to represent upward transmittance
# 	E::Vector{FT}  = E,    # extra terrestrial irradiance (W/m2/µm)
# 	nPoly::Int     = 3,   # order of polynomial
# 	nPC_trans::Int = 5,   # number of PC used
# 	sza::FT = 45.,        # solar zenith angle (˚)
# 	vza::FT = 30.,        # viewing zenith angle (˚)
# 	λ₀::FT  = 683.,
# 	σ::FT   = 10, 		  # λ₀ and σ prescribe the Gaussian shape of SIF
# 	if_log_SVD::Bool = false
# 						  # forward model is slightly different if the SVD is conducted in log space
# ) where{FT}
# 	#=
# 	- this function gives TOA reflecting FLUX (unit: W/m2/µm)
# 	- the column state vector x includs (in order) :
# 		nPoly: polynomial term
# 		nPC_trans: one-way transmittance Pc
# 		SIF: peak emission (W/m2/µm/sr), multiplying by pi to be consistent with other terms
# 	=#
	
# 	# --- surface reflectance ---
# 	poly = Polynomial(x[1:nPoly+1]);
# 	λm   = mean(λ);  # the polynomial is evaluated after centralizing the wavelength
# 	ρₛ   = poly.(λ .- λm);

# 	# --- check the dimension of PC matrix ---
# 	if size(trans_mat)[1] != nPC_trans
# 		println("Dimension mismatch: trans_mat")
# 	end

# 	if if_log_SVD
# 		println("principle components from log-SVD! to be updated 😂")
# 		return
# 	end
	
# 	T_up = x[(nPoly+2):(nPoly+nPC_trans+1)] * trans_mat;  # get upward transmittance
# 	svc  =  (secd(sza) + secd(vza)) / secd(vza);          # solar-viewer correction
# 	T₂   = exp.(svc .* log.(T_up));                       # get two-way transmittance
	
# 	# --- shape of SIF --- 
# 	SIF  = x[end] * exp.( - ( λ .- λ₀ ).^2 ./ ( 2 * σ^2 ) );

# 	# TOA upwelling flux
# 	R_TOA = E .* cosd(sza) .* ρₛ .* T₂ .+ pi .* SIF .* T_up

# 	return R_TOA
# end

# ╔═╡ 169c4d69-9b8f-425e-a1b6-9e03c2c49236
md"""
---
### SVD: Get basis functions & compare different spectral resolution
"""

# ╔═╡ f029a56f-f4bd-4fbc-b81b-9b26c5110840
begin
	# use file `transmittance_winter_FineWvResModel`
	summer = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/transmittance_summer_FineWvResModel.nc");
	winter = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/transmittance_winter_FineWvResModel.nc");
	println("Opened datasets.")

	
	temp  = cat(summer["temperature"][:,:], winter["temperature"][:,:], dims=1);
	psurf = cat(summer["pressure"][:], winter["pressure"][:], dims=1);
	q     = cat(summer["q"][:,:], winter["q"][:,:], dims=1);
	AMF   = cat(summer["AMF"][:], winter["AMF"][:], dims=1);
	trans = cat(summer["transmittance"][:,:], winter["transmittance"][:,:], dims=1);
	println("\nConcatenated!")

	bands  = summer["band"][:];
end

# ╔═╡ 5d629d4f-eb2a-42cd-8756-6d38e4079ece
begin
	summer1 = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/transmittance_summer_DefaultModel.nc");
	winter1 = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/transmittance_winter_DefaultModel.nc");
	println("Opened datasets.")

	
	temp1  = cat(summer1["temperature"][:,:], winter1["temperature"][:,:], dims=1);
	psurf1 = cat(summer1["pressure"][:], winter1["pressure"][:], dims=1);
	q1     = cat(summer1["q"][:,:], winter1["q"][:,:], dims=1);
	AMF1   = cat(summer1["AMF"][:], winter1["AMF"][:], dims=1);
	trans1 = cat(summer1["transmittance"][:,:], winter1["transmission"][:,:], dims=1);
	println("\nConcatenated! $(size(trans1)[1]) profiles in total.")

	bands1  = summer1["band"][:];
end

# ╔═╡ 7174ab61-ee4c-45b6-9861-8213ef5f7902
λ_min = 620.;λ_max = 860.;

# ╔═╡ f6d0e8f6-c251-4ac5-a359-fc7bba02c875
begin
	LowResSVD = PACE_SIF.Spectral_SVD(trans, bands, λ_min=λ_min, λ_max=λ_max)
	HighResSVD = PACE_SIF.Spectral_SVD(trans, bands, λ_min=λ_min, λ_max=λ_max)
end

# ╔═╡ 0d7e9b30-ba60-4f91-bd21-e502c88d9b34
LowResSVD.Loading

# ╔═╡ abdbef02-361a-4caf-ba97-529f3294023a
begin
	nPC = 7;
	gr()
	MyLayout = (nPC, 1);

	p = plot(LowResSVD.band, LowResSVD.PrinComp[:,1],
			 label = "PC1 ($(round.(LowResSVD.VarExp[1], digits=3))%)",
			 subplot = 1,
			 legend=:outerright,
			 layout = MyLayout,
			 size=(800, 900),
	)
	plot!(p, HighResSVD.band, HighResSVD.PrinComp[:,1],
		 label = "PC1 ($(round.(HighResSVD.VarExp[1], digits=3))%)",
		 subplot = 1,
		 legend=:outerright,
		 layout = MyLayout,
	)
	
	for i in 2:nPC
		plot!(p, LowResSVD.band, LowResSVD.PrinComp[:,i],
				 label = "PC$i ($(round.(LowResSVD.VarExp[i], digits=3))%)",
				 subplot = i,
				 legend=:outerright,
				 layout = MyLayout,
		)
		plot!(p, HighResSVD.band, HighResSVD.PrinComp[:,i],
				 label = "PC$i ($(round.(HighResSVD.VarExp[i], digits=3))%)",
				 subplot = i,
				 legend=:outerright,
				 layout = MyLayout,
		)
	end
	plot!(plot_title="Use different resolution in simulation\nBlue [∆ν=0.01] vs. Red [∆ν=0.1]", plot_titlefontsize=12)
	
	current(p)
end

# ╔═╡ 30c7a602-bd01-4bd2-bbb2-616c3ba9a9f3
md"""
---
### Read in OCI data -> select some of them as sample data
"""

# ╔═╡ 59289b0a-fb24-4a8a-916c-b7ac78c31eca
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
	# select the band
	
end

# ╔═╡ 8a01a1b2-0444-4619-bb88-42e1d0d77da2
findall(coalesce.(nflh .> 0.4, false))

# ╔═╡ 26e1cafa-b0ba-4a9e-98b3-d21e34884e6e
begin
	# plot(red_band, E)
	plot(oci_band, R_TOA, size=(600, 300))
	xlabel!("wv [nm]")
	ylabel!("TOA Rad [W/m2/µm/sr]")
end

# ╔═╡ 79441619-248b-48bc-83f7-15b879821c8f
md"""
> ##### find proper polynomial coeffs
"""

# ╔═╡ b0fe3954-9a9a-4c47-84f1-f0b9764552b2
begin
	degree = 2;
	ind_fit = findall( (oci_band .< 683) .|| (770 .< oci_band .< 810) .|| (oci_band .> 840) );
		# findall( ( λ_min .< oci_band .< λ_max ) )
		# findall( (oci_band .< 683) .|| (770 .< oci_band .< 810) );
	@show "choose some baseline wavelength", ind_fit
	# centralize
	x_mean = mean(oci_band);
	xs     = collect(skipmissing(oci_band[ind_fit] .- x_mean));
	ys     = collect(skipmissing(rhot[ind_fit]));
	poly_fitted = fit(xs, ys, degree);
	poly_fitted2 = fit(xs .+ x_mean, ys, degree);
	yfit   = poly_fitted.(xs)
	yfit2   = poly_fitted2.(xs .+ x_mean)
	
	plot(oci_band, rhot, size=(600, 300), 
		label="data", linewidth=2)
	# plot!(poly_fitted, extrema(xs)...)
	plot!(xs .+ x_mean, yfit,
		label="centralized fit", linewidth=2)
	plot!(xs .+ x_mean, yfit2,
		label="wavelength as xs", linewidth=2)
	xlabel!("wavelength [nm]")
	ylabel!("rhot")
	title!("degree=$degree")
end

# ╔═╡ ed5ea8ee-ab34-428c-8a60-22947605719b
@show poly_fitted, poly_fitted[:]

# ╔═╡ 8f8440f0-ee10-46d7-a8bb-4fef2d9f31e0
md"""
> ##### Proper strength of transmittance vectors
"""

# ╔═╡ be926da4-e47a-45ec-842a-f7ca9de7f3fa
begin
	len = size(HighResSVD.PrinComp)[1];
	nPC_trans = 10;
	PC_coeff  = diagm(HighResSVD.VarExp[1:nPC_trans]) * HighResSVD.Loading[1:nPC_trans,2800] * 10
end

# ╔═╡ 9d721c7d-a471-4017-8d36-b858893a5d4e
md"""
> ##### Take these as priori, evaluate Jacobians
polynomial, coeff. of PCs, and SIF (nFLH)
"""

# ╔═╡ 0a978b0c-1cf4-4e35-a072-611819b2204c
begin
	x = zeros(degree + nPC_trans + 2);
	x[1:degree+1] = poly_fitted[:];
	x[degree+2 : degree+1+nPC_trans] = PC_coeff;
	x[end] = nflh[pixel, scan];
	println("priori constructed!\n", x)
end

# ╔═╡ 41f16839-c13c-489c-91a4-3a09479a9457
# y and x are only used for type and shape information and are not stored in the returned DiffResult
# some random generated in `result` at this stage, just allocate some storage
result = DiffResults.JacobianResult(zeros(length(oci_band)), x);

# ╔═╡ 98174303-e52c-46ec-8bec-1a8acbf4e76d
function forward_model(
	x;
	λ  = λ,    # fiting window (nm)
	trans_mat = trans_mat,   
	# SVD.U truncated matrix to represent upward transmittance
	E  = E,   # extra terrestrial irradiance (W/m2/µm)
	nPoly::Int     = 3,   # order of polynomial
	nPC_trans::Int = 5,   # number of PC used
	sza = 45.,        # solar zenith angle (˚)
	vza = 30.,        # viewing zenith angle (˚)
	λ₀  = 683.,
	σ   = 10, 		  # λ₀ and σ prescribe the Gaussian shape of SIF
	if_log_SVD::Bool = false
						  # forward model is slightly different if the SVD is conducted in log space
)
	#=
	- this function gives TOA reflecting FLUX (unit: W/m2/µm)
	- the column state vector x includs (in order) :
		nPoly: polynomial term
		nPC_trans: one-way transmittance Pc
		SIF: peak emission (W/m2/µm/sr), multiplying by pi to be consistent with other terms
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
	svc  = (secd(sza) + secd(vza)) / secd(vza);          # solar-viewer correction
	T₂   = exp.( svc .* log.(T_up) );                    # get two-way transmittance
	
	# --- shape of SIF --- 
	SIF  = x[end] * exp.( - ( λ .- λ₀ ).^2 ./ ( 2 * σ^2 ) );

	# TOA upwelling flux
	R_TOA = E .* cosd(sza) .* ρₛ .* T₂ .+ pi .* SIF .* T_up

	return R_TOA
end

# ╔═╡ a8d6634d-2482-48c3-9173-52b198ad3dfc
TryResult = forward_model(
					x,
					λ=collect(skipmissing(oci_band)),
					trans_mat=HighResSVD.PrinComp[:, 1:nPC_trans],
					E=E,
					nPoly=degree,
					nPC_trans=nPC_trans,
					sza=sza,
					vza=vza,
				)

# ╔═╡ abe8c084-b018-4926-b762-b4263c03b56e
begin
	plot(oci_band, TryResult,
		label="priori",
		size=(600, 300)
	)
	plot!(oci_band, R_TOA,
		label="measurement"
	)
	xlabel!("wv [nm]")
	ylabel!("TOA Rad [W/m2/µm/sr]")
end

# ╔═╡ aeb6c5da-3ce0-4092-a502-aadf36519f7c
begin
	# model for Jacobian (partial function to include parameters)
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

# ╔═╡ 9c788de3-982d-4fd3-8cde-d036d17a6720
@time ForwardDiff.jacobian!(result, model_for_jacobian, x);

# ╔═╡ 729ddce4-c5d8-410c-a65c-ad640d3ae1bf
typeof(x)

# ╔═╡ 7050b5e9-fb25-4d2e-98fc-cf354b2dc4fa
begin
	# extract Jacobian
	K = DiffResults.jacobian(result);
	println("size of K: ", size(K))
end

# ╔═╡ 9cbe62b6-6efa-4982-bd1a-e9515aa922d3
F = DiffResults.value(result); # i.e., TryResult (priori)

# ╔═╡ 2aea0e9a-4f6d-470b-b1f6-d16c3f73c034
begin
	label_poly = ["p$i" for i in 0:degree];
	label_PC   = ["PC$j" for j in 1:nPC_trans];
	MyLabel    = reshape(vcat(label_poly, label_PC, "SIF"), 1, :);
	plot(oci_band, K, 
		label=MyLabel
	)
	title!("Jacobian of SIF component (magnitude)")
	xlabel!("[nm]")
end

# ╔═╡ a55d808a-866c-4a4e-a117-2e2cbedb4e32
println("priori: ", x)

# ╔═╡ 30fdfd3c-d501-46ac-94df-3c6f918a43a5
md"""
> ##### Construct error convariance before retrieval
"""

# ╔═╡ 3848c0bf-3b2b-4002-9660-3796e7b1a59b
begin
	# priori
	n_state = length(x);
	Sₐ = zeros(n_state,n_state);
	# the polynomial terms are not constraining retrieval at all!
	for i=1:(degree+1)
	    Sₐ[i,i] = 1e20;
	end
	# error in coeffs of PC are proportional to their explained variance? - now just assign a constant
	# HighResSVD.VarExp or (rel_error*HighResSVD.VarExp)^2   
	rel_error   = 2.5
	for i=(degree+2):(degree+nPC_trans+1)
	    Sₐ[i,i] = rel_error;
	end
	# error in SIF is determined by the shape
	Sₐ[end, end] = 10.;

	# measurement noise
	noise = .01
	Se = Diagonal((ones(length(oci_band)).*noise).^2);

	# @show Sₐ
	# @show Se
end

# ╔═╡ e7aa5697-ad4a-4a28-b0d7-58f367faa201
Sₐ

# ╔═╡ 91514d8e-c770-4f64-bc33-5de86ab7c772
md"""
The time to compute the inverse of a matrix depends more on its dimension than on its structure for a general, dense matrix. However, the structure becomes dominant when it allows for specialized, much faster algorithms.

In Julia, the calculation of the inverse depends heavily on the structure of the matrix, as the system intelligently dispatches to specialized algorithms when a specific structure is identified and represented by a dedicated type. If no special structure is present or known, it defaults to general algorithms whose cost is primarily determined by the matrix's dimension.
"""

# ╔═╡ 6a350d74-d86c-4112-8836-c99bcf67eb8b
# gaining matrix G: this calculation involves 3 inverse
@time G = inv( K'inv(Se)K + inv(Sₐ) )K'inv(Se);

# ╔═╡ fd11e0b8-8bd0-4313-b2f0-64936b614a2c
# gaining matrix G: this calculation involves one inverse
@time G2 = Sₐ * K' * inv(K*Sₐ*K' + Se);

# ╔═╡ 2498e131-7b2a-48ba-a618-cb4f0ae990e5
md"""
🟠 My conclusion: since Se is often a diagonal matrix, this special structure renders a faster calculation of inverse despite larger size.

🟣 See Rodgers book Eqn(3.27) at p56 for two equivalent formula of gaining matrix G.

🔵 Realize calculation using G₂ is potentially unstable, if Sa has elements with orders of magnitude different, which preserves in Sₐ * K' term...
"""

# ╔═╡ 9f99b16d-c6e6-415e-a547-5d6b998ae589
x

# ╔═╡ 6b2b3870-2383-4301-b5ef-7ca9933a1a9d
G2

# ╔═╡ bc3a3207-5ab5-4b8a-a5e1-500dabf77198
md"""
> ##### Retrieval and evaluation
"""

# ╔═╡ a7e14103-51d2-4dd1-a545-ab1d8b973bec
begin
	y_prime = R_TOA .- F;
	x̂₁  = x + G * y_prime;
	x̂₁_G2  = x + G2 * y_prime;
	println("iter=1, posterior=$x̂₁")
	println("SIF=$(x̂₁[end])")
	# F(x̂)
	ŷ₁ = model_for_jacobian(x̂₁);
end

# ╔═╡ 8a81bb41-a313-46ef-a258-255248c125cf
md"""
🔵 Again, Ŝ₂ is unstable (see below), but there's a shortcut calculating posterior uncertainty, shown as Ŝ₃ (it's much faster!!)
"""

# ╔═╡ 4f88ebb9-3684-46b8-b000-fa643406cf86
# error covariance
@time Ŝ₁ = inv(inv(Sₐ) + K'inv(Se)K);

# ╔═╡ d4f378b8-5c06-4b44-9b5a-e23095312560
# complicated expression but might be faster?
@time Ŝ₂ = Sₐ - Sₐ * K' * inv(K*Sₐ*K' + Se) * K * Sₐ;

# ╔═╡ d1bca496-a43f-4a05-a4b2-d2914d8f555b
# averaging kernel
A = G*K;

# ╔═╡ 6fd77459-838a-40f7-aec4-8e8767a35f6e
# the third way to compute it with G / averaging kernel at hand
begin
	# identical matrix
	eye = I(n_state);
	@time Ŝ₃ = (eye - G*K)Sₐ;
	@time Ŝ₄ = (eye - A)Sₐ;
end

# ╔═╡ e77ce994-f43d-4601-855c-ce4e1e659f36
begin
	plot(oci_band, TryResult, label="priori", size=(600, 300))
	plot!(oci_band, model_for_jacobian(x̂₁), label="posterior, iter=1, G");
	# plot!(oci_band, model_for_jacobian(x̂₂), label="posterior, iter=1, G₂");
	plot!(oci_band, R_TOA, label="measurement")
	title!("Retrieval")
end

# ╔═╡ 53ad8dbf-22f1-4ac2-abed-0851d2b82848
begin
	plot(oci_band, model_for_jacobian(x̂₁) .- R_TOA, label="posterior, iter=1, G", size=(600, 300));
	title!("ŷ - y")
end

# ╔═╡ f6644046-ae9e-4f9b-a646-55d05cd82792
begin
	heatmap(A, size=(450, 400))
	title!("Averaging kernel?")
end

# ╔═╡ e2972674-6697-441f-9b5b-f5b30591304e
begin
	heatmap(Ŝ₁, size=(450, 400))
	title!("posterior covariance, Ŝ=(I-A)Sₐ")
end

# ╔═╡ 92852751-2dfb-4247-91ba-5dbcfa6e03a6
md"""
> ##### One more iteration
"""

# ╔═╡ ec94d132-a4ba-4172-86f6-33d4fc618346
begin
	result1 = DiffResults.JacobianResult(zeros(length(oci_band)), x)
	# evaluate Jacobian and x1
	ForwardDiff.jacobian!(result1, model_for_jacobian, x̂₁);
	K₂ = DiffResults.jacobian(result1);
	# F(x) from the last retrieval
	F₂ = DiffResults.value(result1);
	G₂ = inv( K₂'inv(Se)K₂ + inv(Sₐ) )K₂'inv(Se);
	y₂ = R_TOA .- F₂ .- K₂ * (x̂₁ - x);
	x̂₂ = x + G₂ * y₂;
end

# ╔═╡ 3264f809-07d3-490a-8a25-18befa1abdba
begin
	println("iter=2, posterior=$x̂₂")
	println("SIF=$(x̂₂[end])")
end

# ╔═╡ a7eb4318-f0c2-4033-b815-fe2692d83ef9
begin
	plot(oci_band, K₂, 
		label=MyLabel,
		linestyle=:dash
	)
	plot!(oci_band, K, 
		label=MyLabel
	)
	title!("Jacobian")
	xlabel!("[nm]")
	ylims!(-30,30)
end

# ╔═╡ 108fb8e3-bae8-4edb-9c81-0ea9ebcf6d2c
begin
	plot(oci_band, TryResult, label="priori", size=(600, 300))
	plot!(oci_band, R_TOA, label="measurement")
	plot!(oci_band, model_for_jacobian(x̂₁), label="posterior, iter=1");
	plot!(oci_band, model_for_jacobian(x̂₂), label="posterior, iter=2");
	title!("Retrieval")
end

# ╔═╡ 71061990-28b7-4af5-8b34-6bf710206c31
begin
	result2 = DiffResults.JacobianResult(zeros(length(oci_band)), x̂₂)
	# evaluate Jacobian at x2
	ForwardDiff.jacobian!(result2, model_for_jacobian, x̂₂);
	K₃ = DiffResults.jacobian(result2);
	# F(x) from the last retrieval
	F₃ = DiffResults.value(result2);
	G₃ = inv( K₃'inv(Se)K₃ + inv(Sₐ) )K₃'inv(Se);
	y₃ = R_TOA .- F₃ .- K₂ * (x̂₂ - x);
	x̂₃ = x + G₃ * y₃;

	println("iter=3, posterior=$x̂₃")
	println("SIF=$(x̂₃[end])")
end

# ╔═╡ 428e458f-02b3-43ce-986f-adf017737050
begin
	plot(oci_band, TryResult, label="priori", size=(600, 300), linestyle=:dash)
	plot!(oci_band, R_TOA, label="measurement", linewidth=3)
	plot!(oci_band, model_for_jacobian(x̂₁), label="posterior, iter=1");
	plot!(oci_band, model_for_jacobian(x̂₂), label="posterior, iter=2");
	plot!(oci_band, model_for_jacobian(x̂₃), label="posterior, iter=3");
	title!("Retrieval")
end

# ╔═╡ Cell order:
# ╠═bbd4b752-5e7e-11f0-0ac9-43833ecaae94
# ╠═9f7e209e-454c-4cc6-9d5c-838c373620c1
# ╠═78653793-ece0-49b3-ace8-24ac0a84a695
# ╠═4a17175f-d424-4f09-95b9-407f2cbd74f1
# ╠═3452c8b6-884e-41d0-a52a-a033545ed5f5
# ╟─8e1a0685-bcec-45ac-aab9-87cca07361bc
# ╟─169c4d69-9b8f-425e-a1b6-9e03c2c49236
# ╠═7a8de7e1-a009-4cff-8520-f93b7efcecf9
# ╠═f029a56f-f4bd-4fbc-b81b-9b26c5110840
# ╟─5d629d4f-eb2a-42cd-8756-6d38e4079ece
# ╠═7174ab61-ee4c-45b6-9861-8213ef5f7902
# ╠═f6d0e8f6-c251-4ac5-a359-fc7bba02c875
# ╠═0d7e9b30-ba60-4f91-bd21-e502c88d9b34
# ╠═abdbef02-361a-4caf-ba97-529f3294023a
# ╟─30c7a602-bd01-4bd2-bbb2-616c3ba9a9f3
# ╠═59289b0a-fb24-4a8a-916c-b7ac78c31eca
# ╠═8a01a1b2-0444-4619-bb88-42e1d0d77da2
# ╠═26e1cafa-b0ba-4a9e-98b3-d21e34884e6e
# ╟─79441619-248b-48bc-83f7-15b879821c8f
# ╠═b0fe3954-9a9a-4c47-84f1-f0b9764552b2
# ╠═ed5ea8ee-ab34-428c-8a60-22947605719b
# ╟─8f8440f0-ee10-46d7-a8bb-4fef2d9f31e0
# ╠═be926da4-e47a-45ec-842a-f7ca9de7f3fa
# ╟─9d721c7d-a471-4017-8d36-b858893a5d4e
# ╠═0a978b0c-1cf4-4e35-a072-611819b2204c
# ╠═41f16839-c13c-489c-91a4-3a09479a9457
# ╠═98174303-e52c-46ec-8bec-1a8acbf4e76d
# ╠═a8d6634d-2482-48c3-9173-52b198ad3dfc
# ╟─abe8c084-b018-4926-b762-b4263c03b56e
# ╠═aeb6c5da-3ce0-4092-a502-aadf36519f7c
# ╠═9c788de3-982d-4fd3-8cde-d036d17a6720
# ╠═729ddce4-c5d8-410c-a65c-ad640d3ae1bf
# ╠═7050b5e9-fb25-4d2e-98fc-cf354b2dc4fa
# ╠═9cbe62b6-6efa-4982-bd1a-e9515aa922d3
# ╠═2aea0e9a-4f6d-470b-b1f6-d16c3f73c034
# ╠═a55d808a-866c-4a4e-a117-2e2cbedb4e32
# ╟─30fdfd3c-d501-46ac-94df-3c6f918a43a5
# ╠═3848c0bf-3b2b-4002-9660-3796e7b1a59b
# ╠═e7aa5697-ad4a-4a28-b0d7-58f367faa201
# ╟─91514d8e-c770-4f64-bc33-5de86ab7c772
# ╠═6a350d74-d86c-4112-8836-c99bcf67eb8b
# ╠═fd11e0b8-8bd0-4313-b2f0-64936b614a2c
# ╟─2498e131-7b2a-48ba-a618-cb4f0ae990e5
# ╠═9f99b16d-c6e6-415e-a547-5d6b998ae589
# ╠═6b2b3870-2383-4301-b5ef-7ca9933a1a9d
# ╟─bc3a3207-5ab5-4b8a-a5e1-500dabf77198
# ╠═a7e14103-51d2-4dd1-a545-ab1d8b973bec
# ╟─8a81bb41-a313-46ef-a258-255248c125cf
# ╠═4f88ebb9-3684-46b8-b000-fa643406cf86
# ╠═d4f378b8-5c06-4b44-9b5a-e23095312560
# ╠═d1bca496-a43f-4a05-a4b2-d2914d8f555b
# ╠═6fd77459-838a-40f7-aec4-8e8767a35f6e
# ╠═e77ce994-f43d-4601-855c-ce4e1e659f36
# ╠═53ad8dbf-22f1-4ac2-abed-0851d2b82848
# ╠═f6644046-ae9e-4f9b-a646-55d05cd82792
# ╠═e2972674-6697-441f-9b5b-f5b30591304e
# ╟─92852751-2dfb-4247-91ba-5dbcfa6e03a6
# ╠═ec94d132-a4ba-4172-86f6-33d4fc618346
# ╠═3264f809-07d3-490a-8a25-18befa1abdba
# ╟─a7eb4318-f0c2-4033-b815-fe2692d83ef9
# ╠═108fb8e3-bae8-4edb-9c81-0ea9ebcf6d2c
# ╠═71061990-28b7-4af5-8b34-6bf710206c31
# ╠═428e458f-02b3-43ce-986f-adf017737050
