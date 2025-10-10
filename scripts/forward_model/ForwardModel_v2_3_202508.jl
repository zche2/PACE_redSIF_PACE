### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 882b69f8-b302-4149-9161-91eadc453794
import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE");

# ╔═╡ 8e511790-d862-497a-83c4-63daed596bba
using Polynomials, ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics

# ╔═╡ 6cbac1ab-f709-4f7c-8c38-d5eed2ae8bdc
using LegendrePolynomials

# ╔═╡ 9d6fdf32-a70e-41cb-b66e-c44b7954467f
using Parameters

# ╔═╡ 8daf61c0-1ae8-4a53-af13-231c18ffb86f
include("../../PACE_SIF.jl")

# ╔═╡ 783f3b03-4a15-4169-bf1b-1bc01338a853
md"""
> #### Add a Gaussian-shaped Chl absorption (actually transmittance considering backscattering)
![Chl absorption spectra](https://www.oceanopticsbook.info/packages/iws_l2h/conversion/files/a_phyt_species.png)
"""

# ╔═╡ fcfca6d1-715d-4a68-a132-4a6a6dbeb091
md"""
> #### Load transmittance spectra and do SVD
"""

# ╔═╡ 49a138de-74de-4795-b94a-3a0167ebbb33
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

# ╔═╡ 70299897-5cbf-401d-b960-196cbbaeef83
a = mean(trans, dims=1)

# ╔═╡ dad78a91-8378-4ec2-91bf-b36bb9ad1f66
bands[a[:] .> 0.999]

# ╔═╡ bafee797-8ac2-4132-ad6b-b429bf0ffd41
begin
	λ_min = 620.;
	λ_max = 860.;
	# get principal components, variance explained by each component (normalized to 100%), and spatial loading
	HighResSVD = PACE_SIF.Spectral_SVD(trans, bands, λ_min=λ_min, λ_max=λ_max);
end

# ╔═╡ 597be477-0c29-409a-abee-e119a40f8a51
mean(HighResSVD.PrinComp[:,1:4], dims=1)

# ╔═╡ 72e10087-daa4-4d15-b40d-7a149feeefb3
md"""
> ##### sample data from PACE-OCI
"""

# ╔═╡ 94ae14de-c80b-4181-b12b-808ef37e574a
begin
	oci = Dataset(
		"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample_granule_20250501T183011_new.nc");
	pixel  = 244;  # cross-track
	scan   = 16;
	red_band = oci["red_wavelength"][:];
	# red_band = oci["red_bands"][:];
	# cloud    = oci["cloud_flag_dilated"][:, :];
	nflh     = oci["nflh"][:, :];
	println("Read in Dataset")

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

# ╔═╡ edeeeb48-d711-493f-b720-529456d6e40c
begin
	p1 = plot(oci_band, R_TOA, size=(500, 300), label="obs");
	# scatter!(p1, oci_band, R_TOA_fit, label="fitting band (TBD)", markersize=1.5);

	p2 = plot(oci_band, E, size=(500, 300), label="obs");
	# scatter!(p2, oci_band_fit, E[ind_fit], label="solar irr.", markersize=1.5)

	plot(p1, p2, layout=(2, 1))
	ylabel!("W/m2/µm/sr")
end

# ╔═╡ 3f6749aa-bbc2-4437-bc1c-e9f2db087677
begin
	# the PCs look like:
	plot(oci_band, HighResSVD.PrinComp[:,1:3], size=(500, 200))
	title!("eigen vectors")
end

# ╔═╡ 19b01b03-2245-4b97-9319-fba76795c927
HighResSVD.VarExp

# ╔═╡ 0431db0e-81b8-4e03-90ab-d3a0c65ebc59
findall(coalesce.((nflh .> .3) .& (nflh .< .4), false))

# ╔═╡ cf4b298d-d68b-4e5e-810e-f74ad37d02c1
md"""
> ##### Forward model: Start with polynomial fit +Transmittance
$$\rho_{s}=\sum{a_jP_j},\ T(\lambda)=\sum{\beta_i P_i}$$

$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)T(\lambda)}{\pi}$$
 $T(\lambda)$ is set to have a maximum of 1.

"""

# ╔═╡ 577d4731-9120-491a-ba11-75f504fd0c6b
md"""
🟢 $S_{\epsilon}$ and Sa
"""

# ╔═╡ 4a9c1228-5ef9-4a6e-a1a4-922cdb2d9f45
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

	# to make sure I use the right wvlen
	# @show wvlen[snr_ind]
	
	# see instruction in .txt file
	noise   = sqrt.( c1[snr_ind] .+ c2[snr_ind] .* R_TOA);
	Se      = Diagonal(noise.^2);
	println("Get measurement error from PACE SNR file.")
end

# ╔═╡ 24018a70-e6ba-406d-b4db-e848e00aed4a
md"""
🟢 Surface reflectance


For polynomial term, the argument needs to satisfy -1 <= x <= 1.
"""

# ╔═╡ 5aaf56d2-194e-4752-a620-007deafeaa01
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

# ╔═╡ 49b1186d-7fb4-4fa8-843f-57e7c5f1d2c6
λc = center_wavelength(oci_band)

# ╔═╡ 8d980f1d-a327-4bb4-be18-f8830202a4ce
md"""
🟢 Transmittance
"""

# ╔═╡ 56e0071f-196f-4931-9d49-962d901a67e7
function scale_transmittance(T, λ_bl_ind)
	# find max
	T_abs = abs.(T)
	bl_max = maximum(T_abs[λ_bl_ind]);
	# force the mean val to be 1
	T_norm = T_abs ./ bl_max
	return T_norm
end

# ╔═╡ 95e76e2f-9653-40d8-a23b-66b700b02db9
md"""
🟠 to scale the transmittance, I calculated the average of all trans. spectra, choosing bands where the mean trans. consistently > 0.999.

previously, I got from oci_band[sortperm(abs.(m1.K[:,6]))]. 

This was used to rescale the transmittance spectra
"""

# ╔═╡ 19ec88d2-648f-4603-9540-89b2c4a706e4
begin
	bl_wvlen = [668.265, 669.518, 670.755, 671.99, 673.245, 674.505, 675.73, 676.962, 678.205, 679.445, 680.68, 751.79, 753.04, 754.295, 776.832, 779.335, 867.115, 869.615, 872.13];
	# [610.36, 612.732, 615.145, 617.605, 620.06, 622.53, 668.265, 669.518, 670.755, 671.99, 673.245, 674.505, 675.73, 676.962, 678.205, 679.445, 680.68, 751.79, 753.04, 754.295, 776.832, 779.335, 867.115, 869.615, 872.13];
	# this is the least sensitive band
	# [801.883, 734.272, 668.246, 667.007, 736.77, 804.388, 669.486, 731.768, 665.766, 799.374, 847.038]

	# this is baseline band for Rrs
	# bl_wvlen = [649.599976, 650.900024, 652.099976, 653.299988, 654.599976, 655.799988, 657.099976, 658.299988, 659.599976, 710.500000, 711.799988, 713.000000, 714.299988, 716.799988, 719.200012]
	λ_bl_ind = map(bl_wvlen -> argmin(abs.(oci_band .- bl_wvlen)), bl_wvlen);
	oci_band[λ_bl_ind]
end

# ╔═╡ d0816873-8573-4d19-bbc9-ba26ecabf147
md"""
🟢 Jacobian
"""

# ╔═╡ 8799c188-7aa7-4e35-b770-0e89fce160d1
function Jacobian(x, model; len=length(oci_band))
	res = DiffResults.JacobianResult(zeros(len), x);
	ForwardDiff.jacobian!(res, model, x);
	K   = DiffResults.jacobian(res);
	val = DiffResults.value(res);
	return K, val
end

# ╔═╡ 9c791f97-782d-4dd2-8978-ae29398aac51
md"""
🟢 Gain Matrix
"""

# ╔═╡ 4756d0d2-0855-4267-af76-9eec3935602a
md"""
🟢 Forward model
"""

# ╔═╡ b93b300e-bbd5-4374-aacf-2e015a3d32d6
n = 5; nPC = 15;

# ╔═╡ adc4d73e-8b36-48c6-be65-9822e71c9baa
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
	println("Get prior error covariance")
end

# ╔═╡ 2df8cd2d-a9e6-41cc-9499-24f506c2ed03
function GainMatrix(K; Se=Se, Sa=Sa)
	return inv( K'inv(Se)K + inv(Sa) )K'inv(Se)
end

# ╔═╡ 2b1ce2d9-d232-4d29-a87e-596c30d332bb
function forward_model1(
			x; 
			λ = oci_band,     # wavelength range
			nPoly::Int = n,   # degree of polynomials
			nPC::Int   = nPC,   # number of eigen vectors used
			trans_mat  = HighResSVD.PrinComp[:, 1:nPC],
			sza        = sza,
			vza        = vza,
			E          = E,
			λc         = λc,
			λ_bl_ind   = λ_bl_ind
		)
	
	# adjust to [-1,1]
	v     = collectPl.(λc, lmax=nPoly);
	# reflectance
	rho   = hcat(v...)' * x[1 : nPoly+1];
	# transmittance
	T      = trans_mat * x[(nPoly+2):(nPoly+nPC+1)];
	T_norm = scale_transmittance(T, λ_bl_ind);
	# TOA radiance
	rad    = E .* cosd(sza) ./ pi .* T_norm .* rho
	return rad
end

# ╔═╡ f5a3084c-cc48-48f3-99e1-fb568a9aeab6
md"""
> ##### Retrieval: xa + iteration
"""

# ╔═╡ a2d146df-127b-414d-bb3c-2182da643b15
@with_kw struct retrieval
    x    # state vector
	y_x  # value evaluated at x
	K    # Jacobian
	G = GainMatrix(K)         # Gain matrix
	A = G*K                   # averaging kernel
end;

# ╔═╡ dc55c80b-fcb8-46e4-a35a-dd1327ba922f
md"""
🟢 Iteration
"""

# ╔═╡ 020a0a6d-fdfa-4a21-852e-998d95e27300
function iter(
		m ::retrieval,   # retrieval of last iteration
		xa,              # priori
	    R_msr; 
		Sa    = Sa,      # measurements
		model = x -> forward_model1(x)
	)
	
	# get results from last iteration x̂ₙ, note that K and y are evaluated at x̂ₙ
	xn     = m.x
	Kn     = m.K
	yn     = m.y_x
	G      = m.G
	x_n1   = xa .+ G * (R_msr .- yn .+ Kn * (xn .- xa));
	K_n1, y_n1 = Jacobian(x_n1, model);

	# update 
	m_new  = retrieval(
		x   = x_n1,
		y_x = y_n1,
		K   = K_n1,
		G   = GainMatrix(K_n1, Sa=Sa)
	)
	return m_new
end

# ╔═╡ f9f0752c-d19d-4858-aba6-f315336093fa
begin
	K     = E .* cosd(sza) ./ pi .* hcat(collectPl.(λc, lmax=n)...)';
	G     = inv( K'inv(Se)K )K'inv(Se);
	x̂     = G * R_TOA;
	ŷ     = K * x̂;
	tmp   = zeros(nPC-2) .+ .001;
	xa    = [x̂... -10. 0.1 tmp...]';
	rad   = forward_model1(xa, nPoly=n);
	plot(oci_band, R_TOA, size=(500, 200), label="obs.")
	plot!(oci_band, rad, label="initial guess, n=$n, nPC=$nPC")
	title!("TOA radiance (W/m2/µm/sr)", titlefont=10)
end

# ╔═╡ 123f5138-3d38-4764-8c9a-7f06cc2c8d10
begin
	# start from xa
	Ka, ya = Jacobian(xa, x -> forward_model1(x))
	ma = retrieval(x=xa, y_x=ya, K=Ka, G=GainMatrix(Ka, Sa=Sa))
end

# ╔═╡ 2477de81-b60f-4157-817f-881aad0329b5
begin
	T_try = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * ma.x[(n+2):(n+nPC+1)], λ_bl_ind);

	plot(oci_band, T_try, size=(600, 200))
	scatter!(oci_band[λ_bl_ind], T_try[λ_bl_ind], markersize=2.5, label="ref pts to scale the spectra")
end

# ╔═╡ 39beb5af-e208-4adf-ab40-1e21fe474c86
begin
	results_history = []
	push!(results_history, (iteration=0, result_val=ma))
	for i = 1:15
		m  = iter(results_history[i].result_val, xa, R_TOA,
			      Sa=Sa, model=x->forward_model1(x));
		println("Current iteration $(i)")
		# push 
		push!(results_history, (iteration=i, result_val=m))
	end
end

# ╔═╡ f7226689-993f-47dc-afc8-9bd1768aee50
begin
	p11 = plot(oci_band, R_TOA,
		size=(600, 250),
		label="Observed",
		linewidth=3,
		linestyle=:dash,
		color=:black,
		# xlabel="Wavelength [nm]",
		ylabel="W/m²/µm/sr",
		title="Top of Atmosphere Radiance",
		titlefontsize=10,
		xlabelfontsize=10,
		ylabelfontsize=10,
		legend=:outerright,
		# top_margin=5Plots.mm,
		# bottom_margin=5Plots.mm,
		left_margin=5Plots.mm,
		# dpi=400,
	)
	plot!(p11, oci_band, ma.y_x, label="Initial guess", linewidth=2)
	plot!(p11, oci_band, results_history[1].result_val.y_x,
		  label="iter#1", linewidth=1.5, color=:peru)
	plot!(p11, oci_band, results_history[14].result_val.y_x,
		  label="iter#14", linewidth=1.5, color=:blue)
	plot!(p11, oci_band, results_history[15].result_val.y_x,
		  label="iter#15", linewidth=1.5, color=:green)
	
	# title!("Observed & Retrieved Radiance\n degree of polynomial=$n, number of PC=$nPC", titlefont=8)


	# savefig(p, "../../demo_example/Figures/NullRetrieval.png")
end

# ╔═╡ 4a087376-12e9-47c8-8c17-6296a60510d9
begin
	v     = collectPl.(λc, lmax=n);
	rho_a = hcat(v...)' * ma.x[1 : n+1];
	rho_1 = hcat(v...)' * results_history[1].result_val.x[1 : n+1];
	rho_14 = hcat(v...)' * results_history[14].result_val.x[1 : n+1];
	rho_15 = hcat(v...)' * results_history[15].result_val.x[1 : n+1];

	# plot(oci_band, R_TOA, size=(600, 200), label="obs.", linewidth=4, linestyle=:dash, color=:black)
	plot(oci_band, rho_a, label="initial guess", linewidth=1, size=(600, 200))
	plot!(oci_band, rho_1, label="iter#1", linewidth=1)
	plot!(oci_band, rho_14, label="iter#14", linewidth=1)
	plot!(oci_band, rho_15, label="iter#15", linewidth=1)
	title!("Reconstructed Surface Reflectance (W/m2/µm/sr), n=$n, nPC=$nPC", titlefont=10)
end

# ╔═╡ 882e041e-2b15-4784-85be-790785fbfa16
begin
	T1 = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * results_history[1].result_val.x[(n+2):(n+nPC+1)], λ_bl_ind);
	T2 = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * results_history[14].result_val.x[(n+2):(n+nPC+1)], λ_bl_ind);
	T3 = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * results_history[15].result_val.x[(n+2):(n+nPC+1)], λ_bl_ind);

	plot(oci_band, T1, label="iter#1", size=(600, 200))
	plot!(oci_band, T2, label="iter#14")
	plot!(oci_band, T3, label="iter#15")

	title!("Reconstructed Transmittance Spectra")
end

# ╔═╡ dd532004-54dc-45c0-a462-2a17c43d75fc
begin
	# plot(oci_band, R_TOA .- ma.y_x, size=(600, 200), label="initial guess", linewidth=2)
	p12 = plot(oci_band, R_TOA .- results_history[1].result_val.y_x,
		label="iter#1",
		linewidth=1, size=(800, 250),
		xlabel="[nm]",
		ylabel="W/m²/µm/sr",
		title="Residual",
		titlefontsize=10,
		xlabelfontsize=10,
		ylabelfontsize=10,
		# top_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		left_margin=5Plots.mm,
		legend=:outerright,
		color=:peru,
		dpi=400,
	)
	plot!(p12, oci_band, R_TOA .- results_history[12].result_val.y_x, label="iter#12", linewidth=1.5, color=:coral)
	plot!(p12, oci_band, R_TOA .- results_history[14].result_val.y_x, label="iter#14", linewidth=1.5, color=:blue)
	plot!(p12, oci_band, R_TOA .- results_history[15].result_val.y_x, label="iter#15", linewidth=1.5, color=:green)
	# plot!(oci_band, R_TOA .- m10.y_x, label="iter#10", linewidth=1)
	# vline!([683.], label="683 nm")
	# title!("Residual", titlefont=15)

	# savefig("Residual_NullFit.png")
end

# ╔═╡ d91fdc59-af9a-4a49-946a-6a5060536e62
begin
	# plot(oci_band, R_TOA .- ma.y_x, size=(600, 200), label="initial guess", linewidth=2)
	p13 = plot(oci_band, (R_TOA .- results_history[1].result_val.y_x) ./ R_TOA,
		label="iter#1",
		linewidth=1, size=(800, 250),
		xlabel="[nm]",
		ylabel="Rel.",
		title="Residual",
		titlefontsize=10,
		xlabelfontsize=10,
		ylabelfontsize=10,
		# top_margin=5Plots.mm,
		bottom_margin=5Plots.mm,
		left_margin=5Plots.mm,
		legend=:outerright,
		color=:peru,
		dpi=400,
	)
	plot!(p13, oci_band, (R_TOA .- results_history[12].result_val.y_x) ./ R_TOA, label="iter#12", linewidth=1.5, color=:coral)
	plot!(p13, oci_band, (R_TOA .- results_history[14].result_val.y_x) ./ R_TOA, label="iter#14", linewidth=1.5, color=:blue)
	plot!(p13, oci_band, (R_TOA .- results_history[15].result_val.y_x) ./ R_TOA, label="iter#15", linewidth=1.5, color=:green)
	# plot!(oci_band, R_TOA .- m10.y_x, label="iter#10", linewidth=1)
	# vline!([683.], label="683 nm")
	title!(p13, "Rel. Residual", titlefont=15)

	# savefig("Residual_NullFit.png")
end

# ╔═╡ 146a6efb-f8ba-4423-90da-0f37dcb6f6cf
results_history[15].result_val.x

# ╔═╡ 187093c8-fb6f-4afa-8497-13a1173a88db
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
Where $\lambda_0=683$ is the peak wavelength and $\sigma=5$ (to be tuned).
"""

# ╔═╡ 75b831fb-5f82-4c7e-a7fb-daa8ed15d7aa


# ╔═╡ 12c0128b-28e6-4a5e-8da3-da2b7211bba8


# ╔═╡ 39469984-ef88-4a70-bbb6-a665fd51e0d1


# ╔═╡ 2657b89e-7646-47ed-a800-a41f245dfe76


# ╔═╡ 04831ff7-0489-4cd5-aef8-56b01dc41ae5


# ╔═╡ 5c369951-3863-432c-8c97-6bc2937e4cb5


# ╔═╡ e5abf37e-2937-467d-97f9-0d7d55523322


# ╔═╡ 1fb84295-7cac-4103-945f-43d3f1748cf7


# ╔═╡ d34cffcc-24f7-47ab-bea8-cfe645699ad0


# ╔═╡ 6b06645b-9340-423b-8736-83e8c16e9794


# ╔═╡ 51909fd0-cf26-4e8c-929a-2490aeed850b


# ╔═╡ 3311fd3d-4552-45f4-8b02-8fffffbba97d


# ╔═╡ dde7cf1f-0856-4671-bfc0-d277f69f2bfe


# ╔═╡ 5398e782-7046-4904-9186-6b06b9b9101c


# ╔═╡ d578fb7a-7204-4cd7-aa5b-75cf9de635b0


# ╔═╡ 938018dd-9923-4677-a904-b01002587596


# ╔═╡ a2d3ef4e-dc6f-4506-87dd-01af2412b7a6


# ╔═╡ 5644bced-9475-4e94-b287-5cc3b6af0f32


# ╔═╡ 33354c81-129c-41f5-9cef-2bb88239484c


# ╔═╡ 4ff0b225-a432-429c-9b5b-eb5e70120d46


# ╔═╡ 03c43815-b1d0-4b46-bb16-13f89ded8570


# ╔═╡ 53f86a2e-5cb9-4d89-bd1d-fcdf6b3ae671


# ╔═╡ bffaed6b-239d-4789-b211-a2f04a39c922


# ╔═╡ f360e5f9-1949-4d17-85dd-0f329bc89a85


# ╔═╡ f28b764f-7631-41e3-9a14-a6f0af855432


# ╔═╡ 00f7bef7-866d-4e43-97a1-a68539c8d65a


# ╔═╡ ec242bb5-ad9c-47a4-9003-30c4cc66bc31


# ╔═╡ b65a660e-697b-4172-bd5d-a5268a698588


# ╔═╡ 4b2800e4-3a16-4dee-8445-9dd5af3889b6


# ╔═╡ 22153864-ce3d-4e23-9c91-a3605216ecbc


# ╔═╡ 505d93b0-de73-4477-a343-99f04b570d24


# ╔═╡ c11884ce-86db-4dc9-8cf7-f5f6fc62e726


# ╔═╡ 16410d72-27c8-4486-b63e-1672fbf08379


# ╔═╡ 1e186f7e-6180-4d26-9443-b135821418b6


# ╔═╡ f5927565-56bc-4a78-8aee-e8082d5e793f


# ╔═╡ 439dc428-445f-4852-9599-d014db6b7a48


# ╔═╡ 21bfe19b-0ac1-4eb5-bb52-afe13bdc25e1


# ╔═╡ ec9aabfb-55f8-4445-89e2-33c61ad514ce


# ╔═╡ 781270e6-b8fe-4e0e-8003-09dad8a64ffb


# ╔═╡ 171bf4e3-a944-4fee-a309-971219cde2e6


# ╔═╡ f1dd70e1-e973-407a-807a-5a0d4f8d9d31


# ╔═╡ f1cf59d2-c75c-486f-b625-ddf2fa6b178d


# ╔═╡ c6d9dd98-8987-4de8-b2f5-92cc3472e63e


# ╔═╡ 82527c07-a8df-46e0-822b-e5b4b118d1d4


# ╔═╡ b6b5ba2f-e2bb-4517-b229-d8afb94b664c


# ╔═╡ fa52b56b-0ed1-4cce-88f5-e927a22a6740


# ╔═╡ 044c2e94-2a7e-40b9-b165-008171d89d93


# ╔═╡ 2ac3f305-7640-4d1a-a430-82ae20e8eb49


# ╔═╡ c1bad7f3-694b-4935-a8af-f9d26583556a


# ╔═╡ cf71c986-5980-4a2f-90db-acee742df2d4


# ╔═╡ 82fe3334-b060-4c25-a1f2-c488191ed05e


# ╔═╡ 560def2a-2346-46e5-85cc-25630f1d352f


# ╔═╡ 2d1ea09d-82c4-4f5f-b684-a28e3084d6c7


# ╔═╡ c0232ab7-2e13-490c-a534-117548f20003


# ╔═╡ 688980f3-35f2-4dde-8765-3cde4b984e7d


# ╔═╡ e52c9e88-a39d-47a0-8559-521e71dbf7c4


# ╔═╡ 0fcafeb5-774c-4b47-8d12-9f1c1a4c4651


# ╔═╡ 2de4992d-bf0e-4f8d-bbfb-c92b1002e833


# ╔═╡ 58288f8a-de66-4d06-9905-3a86076b6159


# ╔═╡ 6d52f464-b558-4722-934c-62f7dd42ad08


# ╔═╡ 8eecad64-5e90-4a8b-afd4-40d321ea8208


# ╔═╡ 34b5454a-f49c-4f3e-9df2-cbf5f06709c0


# ╔═╡ e76ca676-4ec1-4c3a-9177-596666f269ed


# ╔═╡ 17e16d1d-ca46-4529-9241-0a66fdc52661


# ╔═╡ e15baf76-2be6-4cdc-b68f-ded591feb9b3


# ╔═╡ bb11382f-4b01-4b35-a578-6bc6f105024c


# ╔═╡ 7c144dc1-71e1-4e85-856d-d21395ed8b8e


# ╔═╡ 7176567e-ca03-4ec8-9b81-b8b2449ddf80


# ╔═╡ 68eea47d-9623-4b3a-9de0-bb35934e94d9


# ╔═╡ 95ff596b-25b9-4aa4-aa61-23b72553ac8e


# ╔═╡ 31ca31d7-14eb-4b01-9ebb-b9c135fd7302


# ╔═╡ Cell order:
# ╟─783f3b03-4a15-4169-bf1b-1bc01338a853
# ╠═882b69f8-b302-4149-9161-91eadc453794
# ╠═8e511790-d862-497a-83c4-63daed596bba
# ╠═6cbac1ab-f709-4f7c-8c38-d5eed2ae8bdc
# ╠═9d6fdf32-a70e-41cb-b66e-c44b7954467f
# ╠═8daf61c0-1ae8-4a53-af13-231c18ffb86f
# ╠═fcfca6d1-715d-4a68-a132-4a6a6dbeb091
# ╠═49a138de-74de-4795-b94a-3a0167ebbb33
# ╠═70299897-5cbf-401d-b960-196cbbaeef83
# ╠═dad78a91-8378-4ec2-91bf-b36bb9ad1f66
# ╠═bafee797-8ac2-4132-ad6b-b429bf0ffd41
# ╠═597be477-0c29-409a-abee-e119a40f8a51
# ╟─72e10087-daa4-4d15-b40d-7a149feeefb3
# ╠═94ae14de-c80b-4181-b12b-808ef37e574a
# ╟─edeeeb48-d711-493f-b720-529456d6e40c
# ╠═3f6749aa-bbc2-4437-bc1c-e9f2db087677
# ╠═19b01b03-2245-4b97-9319-fba76795c927
# ╠═0431db0e-81b8-4e03-90ab-d3a0c65ebc59
# ╠═cf4b298d-d68b-4e5e-810e-f74ad37d02c1
# ╟─577d4731-9120-491a-ba11-75f504fd0c6b
# ╠═4a9c1228-5ef9-4a6e-a1a4-922cdb2d9f45
# ╟─24018a70-e6ba-406d-b4db-e848e00aed4a
# ╠═5aaf56d2-194e-4752-a620-007deafeaa01
# ╠═49b1186d-7fb4-4fa8-843f-57e7c5f1d2c6
# ╟─8d980f1d-a327-4bb4-be18-f8830202a4ce
# ╠═56e0071f-196f-4931-9d49-962d901a67e7
# ╟─95e76e2f-9653-40d8-a23b-66b700b02db9
# ╠═19ec88d2-648f-4603-9540-89b2c4a706e4
# ╟─d0816873-8573-4d19-bbc9-ba26ecabf147
# ╠═8799c188-7aa7-4e35-b770-0e89fce160d1
# ╟─9c791f97-782d-4dd2-8978-ae29398aac51
# ╟─4756d0d2-0855-4267-af76-9eec3935602a
# ╠═b93b300e-bbd5-4374-aacf-2e015a3d32d6
# ╠═adc4d73e-8b36-48c6-be65-9822e71c9baa
# ╠═2df8cd2d-a9e6-41cc-9499-24f506c2ed03
# ╠═2b1ce2d9-d232-4d29-a87e-596c30d332bb
# ╟─f5a3084c-cc48-48f3-99e1-fb568a9aeab6
# ╠═a2d146df-127b-414d-bb3c-2182da643b15
# ╟─dc55c80b-fcb8-46e4-a35a-dd1327ba922f
# ╠═020a0a6d-fdfa-4a21-852e-998d95e27300
# ╠═f9f0752c-d19d-4858-aba6-f315336093fa
# ╠═123f5138-3d38-4764-8c9a-7f06cc2c8d10
# ╠═2477de81-b60f-4157-817f-881aad0329b5
# ╠═39beb5af-e208-4adf-ab40-1e21fe474c86
# ╟─f7226689-993f-47dc-afc8-9bd1768aee50
# ╟─4a087376-12e9-47c8-8c17-6296a60510d9
# ╟─882e041e-2b15-4784-85be-790785fbfa16
# ╟─dd532004-54dc-45c0-a462-2a17c43d75fc
# ╠═d91fdc59-af9a-4a49-946a-6a5060536e62
# ╠═146a6efb-f8ba-4423-90da-0f37dcb6f6cf
# ╠═187093c8-fb6f-4afa-8497-13a1173a88db
# ╠═75b831fb-5f82-4c7e-a7fb-daa8ed15d7aa
# ╠═12c0128b-28e6-4a5e-8da3-da2b7211bba8
# ╠═39469984-ef88-4a70-bbb6-a665fd51e0d1
# ╠═2657b89e-7646-47ed-a800-a41f245dfe76
# ╠═04831ff7-0489-4cd5-aef8-56b01dc41ae5
# ╠═5c369951-3863-432c-8c97-6bc2937e4cb5
# ╠═e5abf37e-2937-467d-97f9-0d7d55523322
# ╠═1fb84295-7cac-4103-945f-43d3f1748cf7
# ╠═d34cffcc-24f7-47ab-bea8-cfe645699ad0
# ╠═6b06645b-9340-423b-8736-83e8c16e9794
# ╠═51909fd0-cf26-4e8c-929a-2490aeed850b
# ╠═3311fd3d-4552-45f4-8b02-8fffffbba97d
# ╠═dde7cf1f-0856-4671-bfc0-d277f69f2bfe
# ╠═5398e782-7046-4904-9186-6b06b9b9101c
# ╠═d578fb7a-7204-4cd7-aa5b-75cf9de635b0
# ╠═938018dd-9923-4677-a904-b01002587596
# ╠═a2d3ef4e-dc6f-4506-87dd-01af2412b7a6
# ╠═5644bced-9475-4e94-b287-5cc3b6af0f32
# ╠═33354c81-129c-41f5-9cef-2bb88239484c
# ╠═4ff0b225-a432-429c-9b5b-eb5e70120d46
# ╠═03c43815-b1d0-4b46-bb16-13f89ded8570
# ╠═53f86a2e-5cb9-4d89-bd1d-fcdf6b3ae671
# ╠═bffaed6b-239d-4789-b211-a2f04a39c922
# ╠═f360e5f9-1949-4d17-85dd-0f329bc89a85
# ╠═f28b764f-7631-41e3-9a14-a6f0af855432
# ╠═00f7bef7-866d-4e43-97a1-a68539c8d65a
# ╠═ec242bb5-ad9c-47a4-9003-30c4cc66bc31
# ╠═b65a660e-697b-4172-bd5d-a5268a698588
# ╠═4b2800e4-3a16-4dee-8445-9dd5af3889b6
# ╠═22153864-ce3d-4e23-9c91-a3605216ecbc
# ╠═505d93b0-de73-4477-a343-99f04b570d24
# ╠═c11884ce-86db-4dc9-8cf7-f5f6fc62e726
# ╠═16410d72-27c8-4486-b63e-1672fbf08379
# ╠═1e186f7e-6180-4d26-9443-b135821418b6
# ╠═f5927565-56bc-4a78-8aee-e8082d5e793f
# ╠═439dc428-445f-4852-9599-d014db6b7a48
# ╠═21bfe19b-0ac1-4eb5-bb52-afe13bdc25e1
# ╠═ec9aabfb-55f8-4445-89e2-33c61ad514ce
# ╠═781270e6-b8fe-4e0e-8003-09dad8a64ffb
# ╠═171bf4e3-a944-4fee-a309-971219cde2e6
# ╠═f1dd70e1-e973-407a-807a-5a0d4f8d9d31
# ╠═f1cf59d2-c75c-486f-b625-ddf2fa6b178d
# ╠═c6d9dd98-8987-4de8-b2f5-92cc3472e63e
# ╠═82527c07-a8df-46e0-822b-e5b4b118d1d4
# ╠═b6b5ba2f-e2bb-4517-b229-d8afb94b664c
# ╠═fa52b56b-0ed1-4cce-88f5-e927a22a6740
# ╠═044c2e94-2a7e-40b9-b165-008171d89d93
# ╠═2ac3f305-7640-4d1a-a430-82ae20e8eb49
# ╠═c1bad7f3-694b-4935-a8af-f9d26583556a
# ╠═cf71c986-5980-4a2f-90db-acee742df2d4
# ╠═82fe3334-b060-4c25-a1f2-c488191ed05e
# ╠═560def2a-2346-46e5-85cc-25630f1d352f
# ╠═2d1ea09d-82c4-4f5f-b684-a28e3084d6c7
# ╠═c0232ab7-2e13-490c-a534-117548f20003
# ╠═688980f3-35f2-4dde-8765-3cde4b984e7d
# ╠═e52c9e88-a39d-47a0-8559-521e71dbf7c4
# ╠═0fcafeb5-774c-4b47-8d12-9f1c1a4c4651
# ╠═2de4992d-bf0e-4f8d-bbfb-c92b1002e833
# ╠═58288f8a-de66-4d06-9905-3a86076b6159
# ╠═6d52f464-b558-4722-934c-62f7dd42ad08
# ╠═8eecad64-5e90-4a8b-afd4-40d321ea8208
# ╠═34b5454a-f49c-4f3e-9df2-cbf5f06709c0
# ╠═e76ca676-4ec1-4c3a-9177-596666f269ed
# ╠═17e16d1d-ca46-4529-9241-0a66fdc52661
# ╠═e15baf76-2be6-4cdc-b68f-ded591feb9b3
# ╠═bb11382f-4b01-4b35-a578-6bc6f105024c
# ╠═7c144dc1-71e1-4e85-856d-d21395ed8b8e
# ╠═7176567e-ca03-4ec8-9b81-b8b2449ddf80
# ╠═68eea47d-9623-4b3a-9de0-bb35934e94d9
# ╠═95ff596b-25b9-4aa4-aa61-23b72553ac8e
# ╠═31ca31d7-14eb-4b01-9ebb-b9c135fd7302
