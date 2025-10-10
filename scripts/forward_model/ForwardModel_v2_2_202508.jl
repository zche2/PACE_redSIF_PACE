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

# ╔═╡ a69eee3d-ff9c-4035-9a23-667322fec0d9
a = mean(trans, dims=1)

# ╔═╡ 7f1126f2-8819-4759-b9d1-3baa2ef38dda
bands[a[:] .> 0.999]

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
		"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample_granule_20250501T183011_new.nc");
	pixel  = 517;  # cross-track
	scan   = 807;
	red_band = oci["red_wavelength"][:];
	# red_band = oci["red_bands"][:];
	# cloud    = oci["cloud_flag_dilated"][:, :];
	nflh     = oci["nflh"][:, :];
	println("Read in Dataset")
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

# ╔═╡ f80f7a81-000a-4784-9d10-713406303102
begin
	p1 = plot(oci_band, R_TOA, size=(500, 300), label="obs");
	# scatter!(p1, oci_band, R_TOA_fit, label="fitting band (TBD)", markersize=1.5);

	p2 = plot(oci_band, E, size=(500, 300), label="obs");
	# scatter!(p2, oci_band_fit, E[ind_fit], label="solar irr.", markersize=1.5)

	plot(p1, p2, layout=(2, 1))
	ylabel!("W/m2/µm/sr")
end

# ╔═╡ acacde64-9957-409d-ae67-428d13428e9d
begin
	# the PCs look like:
	plot(oci_band, HighResSVD.PrinComp[:,1:3], size=(500, 200))
	title!("eigen vectors")
end

# ╔═╡ c92b4782-6eb8-42a0-83c2-7e7e1d9544fe
HighResSVD.VarExp

# ╔═╡ 0d68673e-5d07-4703-96f6-d1a4ef919f0e
findall(coalesce.((nflh .> 0.6) .& (nflh .< 0.7), false))

# ╔═╡ 434ee765-087e-456a-9696-2ba87fa3c5f3
md"""
> ##### Forward model: Start with polynomial fit +Transmittance
$$\rho_{s}=\sum{a_jP_j},\ T(\lambda)=\sum{\beta_i P_i}$$

$$R_{TOA}=\frac{E(\lambda)cos(SZA)\rho_s(\lambda)T(\lambda)}{\pi}$$
 $T(\lambda)$ is set to have a maximum of 1.

"""

# ╔═╡ 063343a5-5879-4cb7-91ad-5068fe0b33d2
md"""
🟢 $S_{\epsilon}$ and Sa
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

	wv_val  = (λ_min .< wvlen .< λ_max);
	snr_ind = findall((FPA .== "Red") .& wv_val);

	# to make sure I use the right wvlen
	# @show wvlen[snr_ind]
	
	# see instruction in .txt file
	noise   = sqrt.( c1[snr_ind] .+ c2[snr_ind] .* R_TOA);
	Se      = Diagonal(noise.^2);
	println("Get measurement error from PACE SNR file.")
end

# ╔═╡ 40253fb3-981f-4c2d-9f43-ce1c802fc6ef
md"""
🟢 Surface reflectance


For polynomial term, the argument needs to satisfy -1 <= x <= 1.
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

# ╔═╡ 2d8e7b04-4c7f-475a-95d1-bf60318fd3ed
λc = center_wavelength(oci_band)

# ╔═╡ f8ae9017-49ad-48c8-b361-6161674a3175
md"""
🟢 Transmittance
"""

# ╔═╡ ab74fe0c-cfa8-45fc-b4fd-8fea3f93c51b
function scale_transmittance(T, λ_bl_ind)
	# find max
	T_abs = abs.(T)
	bl_max = maximum(T_abs[λ_bl_ind]);
	# force the mean val to be 1
	T_norm = T_abs ./ bl_max
	return T_norm
end

# ╔═╡ 3ccac42c-4a86-41c1-b7a8-5a2a21209d12
md"""
🟠 I got from oci_band[sortperm(abs.(m1.K[:,6]))]. 

This was used to rescale the transmittance spectra
"""

# ╔═╡ dfb8d9ec-b8e3-49e0-81aa-3b172b1a4fa0
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

# ╔═╡ 6f2670c9-0438-49ed-b75b-6c03b3a2e325
md"""
🟢 Jacobian
"""

# ╔═╡ 17141496-cb08-48a5-ac95-237ff94a51ee
function Jacobian(x, model; len=length(oci_band))
	res = DiffResults.JacobianResult(zeros(len), x);
	ForwardDiff.jacobian!(res, model, x);
	K   = DiffResults.jacobian(res);
	val = DiffResults.value(res);
	return K, val
end

# ╔═╡ 7b675457-7113-4ad6-89fc-0a324b6cfe2b
md"""
🟢 Gain Matrix
"""

# ╔═╡ 00a1a14b-4853-496c-914a-53f5a8196753
md"""
🟢 Forward model
"""

# ╔═╡ f9412a2b-f201-4c9f-8fd0-0b10050d2284
n = 5; nPC = 15; 

# ╔═╡ 37260f9b-de49-4420-b46b-16cc46f10ffc
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

# ╔═╡ 62bc5487-1c31-4c22-917f-884b3d8dda61
function GainMatrix(K; Se=Se, Sa=Sa)
	return inv( K'inv(Se)K + inv(Sa) )K'inv(Se)
end

# ╔═╡ f94e037d-fc09-4e26-9cb3-14a043dd057a
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

# ╔═╡ f57ca10e-bb34-45cd-9c98-d771ff81e6ba
md"""
> ##### Retrieval: xa + iteration
"""

# ╔═╡ 036ddbba-5ea4-49a3-9ae8-c0933914f5c0
@with_kw struct retrieval
    x    # state vector
	y_x  # value evaluated at x
	K    # Jacobian
	G = GainMatrix(K)         # Gain matrix
	A = G*K                   # averaging kernel
end;

# ╔═╡ 2af5dfbd-8711-4fd2-9a3c-b6ef4fb9802e
md"""
🟢 Iteration
"""

# ╔═╡ c733529b-4c40-45b6-99f5-37ccb3da0df5
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

# ╔═╡ ad15e249-fa7c-4f44-8c56-df67c1a35556
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

# ╔═╡ a14d6682-c4ee-4962-bc82-f60f924878b0
begin
	# start from xa
	Ka, ya = Jacobian(xa, x -> forward_model1(x))
	ma = retrieval(x=xa, y_x=ya, K=Ka, G=GainMatrix(Ka, Sa=Sa))
end

# ╔═╡ ca699bd2-0786-4742-ade2-8fcff0184da2
begin
	T_try = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * ma.x[(n+2):(n+nPC+1)], λ_bl_ind);

	plot(oci_band, T_try, size=(600, 200))
	scatter!(oci_band[λ_bl_ind], T_try[λ_bl_ind], markersize=2.5, label="ref pts to scale the spectra")
end

# ╔═╡ ffd59e7b-b8fe-4dbe-b996-91981ea32ed5
begin
	# 1st iteration
	m1 = iter(ma, xa, R_TOA, Sa=Sa);
	# 2nd
	m2 = iter(m1, xa, R_TOA, Sa=Sa);
	# 3rd
	m3 = iter(m2, xa, R_TOA, Sa=Sa);
	# 4th
	m4 = iter(m3, xa, R_TOA, Sa=Sa);
	m5 = iter(m4, xa, R_TOA, Sa=Sa);
	m6 = iter(m5, xa, R_TOA, Sa=Sa);
	m7 = iter(m6, xa, R_TOA, Sa=Sa);
	m8 = iter(m7, xa, R_TOA, Sa=Sa);
	m9 = iter(m8, xa, R_TOA, Sa=Sa);
	m10 = iter(m9, xa, R_TOA, Sa=Sa);
	
end

# ╔═╡ b08e27a2-7164-42d4-b8d5-453403c70e67
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
		dpi=400,
	)
	plot!(p11, oci_band, ma.y_x, label="Initial guess", linewidth=2)
	plot!(p11, oci_band, m1.y_x, label="iter#1", linewidth=1.5, color=:peru)
	# plot!(p11, oci_band, m2.y_x, label="iter#2", linewidth=1)
	plot!(p11, oci_band, m9.y_x, label="iter#9", linewidth=1.5, color=:blue)
	plot!(p11, oci_band, m10.y_x, label="iter#10", linewidth=1.5, color=:green)
	
	# title!("Observed & Retrieved Radiance\n degree of polynomial=$n, number of PC=$nPC", titlefont=8)


	# savefig(p, "../../demo_example/Figures/NullRetrieval.png")
end

# ╔═╡ 8ba4d50a-1430-4e99-bb03-3bc98393a610
begin
	v     = collectPl.(λc, lmax=n);
	rho_a = hcat(v...)' * ma.x[1 : n+1];
	rho_1 = hcat(v...)' * m1.x[1 : n+1];
	rho_9 = hcat(v...)' * m2.x[1 : n+1];
	rho_10 = hcat(v...)' * m3.x[1 : n+1];

	# plot(oci_band, R_TOA, size=(600, 200), label="obs.", linewidth=4, linestyle=:dash, color=:black)
	plot(oci_band, rho_a, label="initial guess", linewidth=1, size=(600, 200))
	plot!(oci_band, rho_1, label="iter#1", linewidth=1)
	plot!(oci_band, rho_9, label="iter#2", linewidth=1)
	plot!(oci_band, rho_10, label="iter#3", linewidth=1)
	title!("Reconstructed Surface Reflectance (W/m2/µm/sr), n=$n, nPC=$nPC", titlefont=10)
end

# ╔═╡ bd04263a-7f2c-4ffd-8b55-d1c1bc21fcfa
begin
	T1 = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * m1.x[(n+2):(n+nPC+1)], λ_bl_ind);
	T2 = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * m9.x[(n+2):(n+nPC+1)], λ_bl_ind);
	T3 = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * m10.x[(n+2):(n+nPC+1)], λ_bl_ind);

	plot(oci_band, T1, label="iter#1", size=(600, 200))
	plot!(oci_band, T2, label="iter#9")
	plot!(oci_band, T3, label="iter#10")

	title!("Reconstructed Transmittance Spectra")
end

# ╔═╡ 2ad08504-ef3d-41fb-9c07-d336b5a12368
begin
	# plot(oci_band, R_TOA .- ma.y_x, size=(600, 200), label="initial guess", linewidth=2)
	p12 = plot(oci_band, R_TOA .- m1.y_x,
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
	plot!(p12, oci_band, R_TOA .- m5.y_x, label="iter#5", linewidth=1.5, color=:blue)
	plot!(p12, oci_band, R_TOA .- m6.y_x, label="iter#6", linewidth=1.5, color=:green)
	# plot!(oci_band, R_TOA .- m10.y_x, label="iter#10", linewidth=1)
	# vline!([683.], label="683 nm")
	# title!("Residual", titlefont=15)

	# savefig("Residual_NullFit.png")
end

# ╔═╡ a198ec29-49d3-43d0-91df-719528d7440d
begin
	plot(p11, p12, layout=(2, 1), link=:x, size=(1000, 500))
	plot!(
		xtickfontsize=10, ytickfontsize=10,
		xlabelfontsize=10, ylabelfontsize=10, 
	)
	# savefig("../../demo_example/Figures/Fit_and_Residual.png")
end

# ╔═╡ 9a9a5aab-1de3-4855-8f91-54caa9a4f8e5
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

# ╔═╡ 0ef0267f-dda5-4404-b8f5-3dc0c6019cd9
md"""
🟢 Two-way vs. one-way transmittance
"""

# ╔═╡ def2cc7c-84ea-4805-8b70-03df91f25480
function two_way_trans(T, sza, vza)
	svf = (secd(sza)+secd(vza)) / secd(vza);
	T2  = exp.( svf .* log.(T));
	return T2
end

# ╔═╡ b28d897a-679d-47d1-903e-cb293a4dd98a
function SIF_shape(λ; λ₀=683., σ=5.0)
	return exp.( - ( λ .- λ₀ ).^2 ./ ( 2 * σ^2 ) )
end

# ╔═╡ 6594e6ed-3b33-4aae-8371-3b8aff20f17f
begin
	# one-way to two-way transmittance
	plot(oci_band, T3, size=(500, 200), label="T↑")
	plot!(oci_band, two_way_trans(T3, sza, vza), label="T↓↑")
	title!("From one-way to two-way transmittance")
end

# ╔═╡ 32f84080-38ae-48ca-a20f-15097dab98e2
begin
	# plot spectra of SIF
	plot(oci_band, SIF_shape(oci_band), size=(500, 100), label="SIF (normalzied)" )
	xlims!(645, 720)
end

# ╔═╡ 9a86432c-250f-4318-bb33-81ecd8817f57
md"""
🟢 Defining a new forward model
"""

# ╔═╡ 97fd81d0-cf4b-4aa7-af14-110cf6e5d8f5
begin
	# define priori error matrix
	nSIF    = 3;
	# priori cov
	Sa_new  = zeros(n+nPC+nSIF+1, n+nPC+nSIF+1);
	# uodate diagonal term
	for i=1:(n+1)
	    Sa_new[i,i] = 1;     
		# large variance applies no constrain to these polynomial term
	end
	for i=(n+2):(n+nPC+1)
	    Sa_new[i,i] = rel_error .* HighResSVD.VarExp[i - (n+1)];
	end
	# SIF uncertainty
	Sa_new[end-2, end-2] = .5;
	Sa_new[end-1, end-1] = 1e-5;
	Sa_new[end, end] = 1e-5;
	Sa_new
end

# ╔═╡ 1b72a193-059d-4de4-a092-08bc16aee1e0
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
			λc         = λc,
			λ_bl_ind   = λ_bl_ind
		)
	
	# adjust to [-1,1]
	v     = collectPl.(λc, lmax=nPoly);
	# reflectance
	rho   = hcat(v...)' * x[1 : nPoly+1];
	# upward transmittance
	T     = trans_mat * x[(nPoly+2):(nPoly+nPC+1)];
	T1_norm = scale_transmittance(T, λ_bl_ind);
	# downward transmittance (no need to normalize)
	T2_norm = two_way_trans(T1_norm, sza, vza);
	# water-leaving SIF
	SIF_w   = x[end-nSIF+1] .* SIF_shape(λ, λ₀=x[end-nSIF+2], σ=x[end-nSIF+3]);
	# TOA radiance
	rad    = E .* cosd(sza) ./ pi .* T2_norm .* rho + SIF_w .* T1_norm;
	return rad
end


# ╔═╡ 80bbbc67-c9ac-4c17-b623-4cbd0d212a08
begin
	# set inital guess
	x_for_PCs    = zeros(nPC);
	x_for_PCs[1] = -6.;
	x_for_PCs[2] = .1;
	nflh_val = nflh[pixel, scan]
	@show xa_new   = [x̂... x_for_PCs... nflh_val 683. 5.]';
	# @show xa_new   = [x̂... x_for_PCs... 0. 683. 5.]';

	# wrap up in a struct
	Ka_new, ya_new = Jacobian(xa_new, x -> forward_model2(x));
	ma_new = retrieval(x=xa_new, y_x=ya_new, K=Ka_new, G=GainMatrix(Ka_new,Sa=Sa_new));
end

# ╔═╡ 465c7802-3da8-4321-9745-a8f587cc2ba4
begin
	# 1st iteration
	m1_new = iter(ma_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	# 2nd
	m2_new = iter(m1_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	# 3rd
	m3_new = iter(m2_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));	# 4th
	m4_new = iter(m3_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	m5_new = iter(m4_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	m6_new = iter(m5_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	m7_new = iter(m6_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	# more
	m8_new = iter(m7_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	m9_new = iter(m8_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
	m10_new = iter(m9_new, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
end

# ╔═╡ 8d151e51-ed65-4bfb-bd89-734d568691a1
# 10 more times
begin
	m0 = m10_new;
	for i = 1:10
		m  = iter(m0, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));
		m0 = m;
		# println("This time: $(m0[end-nSIF+1, end])")
	end
end

# ╔═╡ 4989b4b5-d6cd-4868-b13c-8db55da442a2
# one more
m = iter(m0, xa_new, R_TOA, Sa=Sa_new, model=x->forward_model2(x));

# ╔═╡ 5a2b0a13-399c-452d-8c14-f32044140aaf
begin
	plot(oci_band, R_TOA, size=(600, 200), label="obs.", linewidth=4, linestyle=:dash, color=:black)
	plot!(oci_band, ma_new.y_x, label="initial guess", linewidth=2)
	plot!(oci_band, m1_new.y_x, label="iter#1", linewidth=1)
	plot!(oci_band, m2_new.y_x, label="iter#2", linewidth=1)
	plot!(oci_band, m5_new.y_x, label="iter#5", linewidth=1)
	plot!(oci_band, m6_new.y_x, label="iter#6", linewidth=1)
	title!("TOA radiance (W/m2/µm/sr), n=$n, nPC=$nPC", titlefont=10)
end

# ╔═╡ 8d34205b-83dc-44b5-8818-2c9baec740bc
begin
	T1_new = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * m1_new.x[(n+2):(n+nPC+1)], λ_bl_ind);
	T2_new = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * m4_new.x[(n+2):(n+nPC+1)], λ_bl_ind);
	T3_new = scale_transmittance(HighResSVD.PrinComp[:, 1:nPC] * m6_new.x[(n+2):(n+nPC+1)], λ_bl_ind);

	plot(oci_band, T1_new, label="iter#1", size=(600, 200))
	plot!(oci_band, T2_new, label="iter#4")
	plot!(oci_band, T3_new, label="iter#6")

	title!("Reconstructed Transmittance Spectra")
end

# ╔═╡ d07562e1-9757-4cbe-9a8b-54273fef5b60
begin
	p21 = plot(oci_band, R_TOA,
		size=(800, 250),
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
		dpi=400,
	)
	plot!(p21, oci_band, ma_new.y_x, label="Initial guess", linewidth=2)
	plot!(p21, oci_band, m1_new.y_x, label="iter#1", linewidth=1.5, color=:peru)
	# plot!(p11, oci_band, m2.y_x, label="iter#2", linewidth=1)
	plot!(p21, oci_band, m9_new.y_x, label="iter#9", linewidth=1.5, color=:blue)
	plot!(p21, oci_band, m10_new.y_x, label="iter#10", linewidth=1.5, color=:green)
	plot!(p21, oci_band, m10.y_x, label="w/o SIF", linewidth=1.5,
		linestyle=:dash,
		color=:silver)
	
	# title!("Observed & Retrieved Radiance\n degree of polynomial=$n, number of PC=$nPC", titlefont=8)


	# savefig(p, "../../demo_example/Figures/NullRetrieval.png")
end

# ╔═╡ 2e92b7d1-555e-4a55-a2e8-aca814f065e5
begin
	# plot(oci_band, R_TOA .- ma.y_x, size=(600, 200), label="initial guess", linewidth=2)
	p22 = plot(oci_band, R_TOA .- m1_new.y_x,
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
	plot!(p22, oci_band, R_TOA .- m9_new.y_x, label="iter#9", linewidth=1.5, color=:blue)
	plot!(p22, oci_band, R_TOA .- m10_new.y_x, label="iter#10", linewidth=1.5, color=:green)
	plot!(p22, oci_band, R_TOA .- m10.y_x, label="w/o SIF", linewidth=1.5,
		linestyle=:dash,
		color=:silver)

	# savefig("Residual_NullFit.png")
end

# ╔═╡ 04c2f928-bd8d-4f16-91b6-1ed764e933f3
begin
	plot(p21, p22, layout=(2, 1), link=:x, size=(1000, 500))
	plot!(
		xtickfontsize=10, ytickfontsize=10,
		xlabelfontsize=10, ylabelfontsize=10, 
	)
	# savefig("../../demo_example/Figures/Fit_and_Residual_withSIF_pixel$(pixel)_scan$(scan)_nflh$(nflh[pixel, scan]).png")
end

# ╔═╡ 92de6876-4693-4734-9a62-3ea6b4170428
begin
	xs = [string("x", i) for i = 1:length(xa_new)]
	heatmap(xs, oci_band, m10_new.K,
			c=:RdBu
	)
	xlabel!("State Vector")
	ylabel!("[nm]")
	title!("Jacobian @ 10th iteration")
end

# ╔═╡ df0375d9-e3b4-4bf7-855b-d4417322e42d
begin
	λ    = oci_band;
	SIF1 = ma_new.x[end-nSIF+1] .* SIF_shape(
		λ, λ₀=ma_new.x[end-nSIF+2], σ=ma_new.x[end-nSIF+3]);
	SIF2 = m1_new.x[end-nSIF+1] .* SIF_shape(
		λ, λ₀=m1_new.x[end-nSIF+2], σ=m1_new.x[end-nSIF+3]);
	SIF3 = m9_new.x[end-nSIF+1] .* SIF_shape(
		λ, λ₀=m9_new.x[end-nSIF+2], σ=m9_new.x[end-nSIF+3]);
	SIF4 = m10_new.x[end-nSIF+1] .* SIF_shape(
		λ, λ₀=m10_new.x[end-nSIF+2], σ=m10_new.x[end-nSIF+3]);

	SIF5 = m0.x[end-nSIF+1] .* SIF_shape(
		λ, λ₀=m0.x[end-nSIF+2], σ=m0.x[end-nSIF+3]);

	SIF6 = m.x[end-nSIF+1] .* SIF_shape(
		λ, λ₀=m.x[end-nSIF+2], σ=m.x[end-nSIF+3]);

	plot(oci_band, SIF1, label="initial guess",
		 linewidth=1, size=(800, 200),
		 legend=:outerright,
		 left_margin=5Plots.mm,
		 bottom_margin=5Plots.mm,
		 dpi=400,
	)
	plot!(oci_band, SIF2, label="iter#1", linewidth=1)
	plot!(oci_band, SIF3, label="iter#9", linewidth=1)
	plot!(oci_band, SIF4, label="iter#10", linewidth=1)
	plot!(oci_band, SIF5, label="iter#19", linewidth=1)
	plot!(oci_band, SIF6, label="iter#20", linewidth=1)
	# xlims!(645, 720)
	title!("Retrieved SIF", titlefont=10)
	xlabel!("Wavelength [nm]")
	ylabel!("W/m²/µm/sr")

	# savefig("../../demo_example/Figures/Estimated_SIF_pixel$(pixel)_scan$(scan)_nflh$(nflh[pixel, scan]).png")
end

# ╔═╡ c1c03191-2846-4aac-997e-be4839c3ed0e
md"""
🟢 Averaging kernel
"""

# ╔═╡ 9faeb35f-61c0-48fb-8806-ac5823db24fe
begin
	heatmap(m.A,
			c=:greys, # <--- Change the colormap here
			clims=(0, 3),  # <--- vmin and vmax
	    	# aspect_ratio=:equal
			size=(450, 450),
			dpi=400,
	)
	title!("Averaging kernel")
	# savefig("../../demo_example/Figures/Averaging_kernel.png")
end

# ╔═╡ 572179ce-96d2-4d40-8e2a-8c54a8c15cfa
md"""
average root mean square as a simple test of convergence
"""

# ╔═╡ 3e513d67-5e6e-4b4a-b65c-10114f66205a
function RMS(y_obs, y_retrieval; nband=length(oci_band))
	totRMS = sum((y_obs .- y_retrieval) .^ 2);
	aveRMS = totRMS / nband;
	return aveRMS
end

# ╔═╡ 776b0716-14da-4fb9-b74b-61b4bb45b961
begin
	@show RMS(R_TOA, m0.y_x) - RMS(R_TOA, m.y_x)
end

# ╔═╡ 81cbeffc-0d5a-468d-8b12-fa6d4e4cd5c2
begin
	# measured uncertainty
	Ŝ = inv( m.K' * inv(Se) * m.K + inv(Sa_new) );
	heatmap(log(Ŝ))
end

# ╔═╡ c63ab3c9-0d61-48f3-9e0c-ee4a91aa095c
md"""
> ##### Use Loss Function Method
"""

# ╔═╡ 4e95cfc7-3cf6-4f67-9ed5-375aa25f848b
@time sol_lm = solve(
					 prob, 
					 LevenbergMarquardt(),
					 store_trace = Val(true),
					
				)

# ╔═╡ 26073b53-d3b3-4280-8eaf-428bb98667e7
@time sol_gn = solve(prob, NewtonRaphson(), store_trace = Val(true))

# ╔═╡ 965a7ebc-8cf9-445f-87e5-02e679443af1
begin
	# vis. reduction in loss function
	fnorms_gn = [h.fnorm for h in sol_gn.trace.history];
	fnorms_lm = [h.fnorm for h in sol_lm.trace.history];
	p1_ = plot(fnorms_gn, yaxis=:log, label="G-N")
	title!("Reduction in Residual (L2-norm)")
	p2_ = plot(fnorms_lm, yaxis=:log, label="L-M")
	plot(p1_, p2_, layout=(2,1), size=(600, 300))
end

# ╔═╡ 0239f879-1c2f-4543-8b27-179047f6e685
sol_lm

# ╔═╡ ab01db6f-09cc-4bee-8ff0-60ec1eb1a682
# sol_lm.u
sol_gn.u

# ╔═╡ ea2f7e0d-9bc0-4856-b34b-d1b4a73411d6
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
						vza=vza, E=E, λc=λc, λ_bl_ind=λ_bl_ind);
	# J w.r.t. x
	L = (ŷ .- y)' * inv(Se) * (ŷ .- y) ;
	R = (x .- xa)' * inv(Sa) * (x .- xa) ;
	J = L .+ R ;
	return J
end
  ╠═╡ =#

# ╔═╡ a2c65a3d-c958-4ca1-b41b-96efd420b077
# ╠═╡ disabled = true
#=╠═╡
function loss_function(x, p)
	# p is params for forward model and error matrix
	λ     = p.λ;       # wavelength range
	nPoly = p.nPoly;   # degree of polynomials
	nPC   = p.nPC;     # number of eigen vectors used
	trans_mat = p.trans_mat;
	E     = p.E;
	xa    = p.xa;      # priori
	ŷ     = p.ŷ;
	sza   = p.sza;
	vza   = p.vza;
	Se    = p.Se;
	Sa    = p.Sa;

	# evaluate forward model @ x
	y = forward_model1(x, λ=oci_band, nPoly=nPoly, nPC=nPC,
						trans_mat=trans_mat, sza=sza,
						vza=vza, E=E, λc=λc, λ_bl_ind=λ_bl_ind);
	# J w.r.t. x
	L = (ŷ .- y)' * inv(Se) * (ŷ .- y) ;
	R = (x .- xa)' * inv(Sa) * (x .- xa) ;
	J = L .+ R ;
	return J
end
  ╠═╡ =#

# ╔═╡ 7a1d8cdc-6703-47c3-8846-e64f5f38f603
# ╠═╡ disabled = true
#=╠═╡
begin
	# define a problem
	param = (
		λ     = oci_band,
		nPoly = n,
		nPC   = nPC,
		trans_mat = HighResSVD.PrinComp[:, 1:nPC],
		E     = E,
		xa    = ma.x,
		ŷ     = R_TOA,
		sza   = sza,
		vza   = vza,
		Se    = Se,
		Sa    = Sa
	)

	x0 = m9.x;
	# define non linear prob
	prob = NonlinearLeastSquaresProblem(loss_function, x0, param);
end
  ╠═╡ =#

# ╔═╡ 4d829873-6484-4b1d-b8d5-752000c1c331
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

	x0 = m10_new.x;
	# define non linear prob
	prob = NonlinearLeastSquaresProblem(loss_function, x0, param);
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
# ╠═a69eee3d-ff9c-4035-9a23-667322fec0d9
# ╠═7f1126f2-8819-4759-b9d1-3baa2ef38dda
# ╠═6f24e4fe-94b5-45bd-bf46-a98a0fdbaf48
# ╠═401b62ff-9966-40b7-ac5d-ed5d704ddda3
# ╟─0ccf2db1-9080-4d29-bfc7-11dffa706f62
# ╠═a42bd26f-46d5-44a4-81d8-7788899b95bc
# ╠═cc1acba5-d114-4579-a64f-8546c2df40b1
# ╟─f80f7a81-000a-4784-9d10-713406303102
# ╟─acacde64-9957-409d-ae67-428d13428e9d
# ╠═c92b4782-6eb8-42a0-83c2-7e7e1d9544fe
# ╠═0d68673e-5d07-4703-96f6-d1a4ef919f0e
# ╟─434ee765-087e-456a-9696-2ba87fa3c5f3
# ╟─063343a5-5879-4cb7-91ad-5068fe0b33d2
# ╟─466c6800-dd8d-4b11-b33b-bff17dfcf387
# ╟─37260f9b-de49-4420-b46b-16cc46f10ffc
# ╟─40253fb3-981f-4c2d-9f43-ce1c802fc6ef
# ╟─9dcc1616-91d6-45d8-9873-2c449b6e321e
# ╟─2d8e7b04-4c7f-475a-95d1-bf60318fd3ed
# ╟─f8ae9017-49ad-48c8-b361-6161674a3175
# ╠═ab74fe0c-cfa8-45fc-b4fd-8fea3f93c51b
# ╟─3ccac42c-4a86-41c1-b7a8-5a2a21209d12
# ╠═dfb8d9ec-b8e3-49e0-81aa-3b172b1a4fa0
# ╟─6f2670c9-0438-49ed-b75b-6c03b3a2e325
# ╟─17141496-cb08-48a5-ac95-237ff94a51ee
# ╟─7b675457-7113-4ad6-89fc-0a324b6cfe2b
# ╟─62bc5487-1c31-4c22-917f-884b3d8dda61
# ╟─00a1a14b-4853-496c-914a-53f5a8196753
# ╠═f9412a2b-f201-4c9f-8fd0-0b10050d2284
# ╠═f94e037d-fc09-4e26-9cb3-14a043dd057a
# ╟─f57ca10e-bb34-45cd-9c98-d771ff81e6ba
# ╠═036ddbba-5ea4-49a3-9ae8-c0933914f5c0
# ╟─2af5dfbd-8711-4fd2-9a3c-b6ef4fb9802e
# ╠═c733529b-4c40-45b6-99f5-37ccb3da0df5
# ╠═ad15e249-fa7c-4f44-8c56-df67c1a35556
# ╟─ca699bd2-0786-4742-ade2-8fcff0184da2
# ╠═a14d6682-c4ee-4962-bc82-f60f924878b0
# ╠═ffd59e7b-b8fe-4dbe-b996-91981ea32ed5
# ╠═b08e27a2-7164-42d4-b8d5-453403c70e67
# ╠═8ba4d50a-1430-4e99-bb03-3bc98393a610
# ╠═bd04263a-7f2c-4ffd-8b55-d1c1bc21fcfa
# ╠═2ad08504-ef3d-41fb-9c07-d336b5a12368
# ╠═a198ec29-49d3-43d0-91df-719528d7440d
# ╟─9a9a5aab-1de3-4855-8f91-54caa9a4f8e5
# ╟─0ef0267f-dda5-4404-b8f5-3dc0c6019cd9
# ╟─def2cc7c-84ea-4805-8b70-03df91f25480
# ╟─b28d897a-679d-47d1-903e-cb293a4dd98a
# ╠═6594e6ed-3b33-4aae-8371-3b8aff20f17f
# ╟─32f84080-38ae-48ca-a20f-15097dab98e2
# ╟─9a86432c-250f-4318-bb33-81ecd8817f57
# ╟─1b72a193-059d-4de4-a092-08bc16aee1e0
# ╠═97fd81d0-cf4b-4aa7-af14-110cf6e5d8f5
# ╠═80bbbc67-c9ac-4c17-b623-4cbd0d212a08
# ╠═465c7802-3da8-4321-9745-a8f587cc2ba4
# ╠═8d151e51-ed65-4bfb-bd89-734d568691a1
# ╠═4989b4b5-d6cd-4868-b13c-8db55da442a2
# ╠═5a2b0a13-399c-452d-8c14-f32044140aaf
# ╠═8d34205b-83dc-44b5-8818-2c9baec740bc
# ╠═d07562e1-9757-4cbe-9a8b-54273fef5b60
# ╠═2e92b7d1-555e-4a55-a2e8-aca814f065e5
# ╠═04c2f928-bd8d-4f16-91b6-1ed764e933f3
# ╟─92de6876-4693-4734-9a62-3ea6b4170428
# ╠═df0375d9-e3b4-4bf7-855b-d4417322e42d
# ╟─c1c03191-2846-4aac-997e-be4839c3ed0e
# ╟─9faeb35f-61c0-48fb-8806-ac5823db24fe
# ╟─572179ce-96d2-4d40-8e2a-8c54a8c15cfa
# ╠═3e513d67-5e6e-4b4a-b65c-10114f66205a
# ╠═776b0716-14da-4fb9-b74b-61b4bb45b961
# ╠═81cbeffc-0d5a-468d-8b12-fa6d4e4cd5c2
# ╟─c63ab3c9-0d61-48f3-9e0c-ee4a91aa095c
# ╟─a2c65a3d-c958-4ca1-b41b-96efd420b077
# ╠═7a1d8cdc-6703-47c3-8846-e64f5f38f603
# ╠═ea2f7e0d-9bc0-4856-b34b-d1b4a73411d6
# ╠═4d829873-6484-4b1d-b8d5-752000c1c331
# ╠═4e95cfc7-3cf6-4f67-9ed5-375aa25f848b
# ╠═26073b53-d3b3-4280-8eaf-428bb98667e7
# ╠═965a7ebc-8cf9-445f-87e5-02e679443af1
# ╠═0239f879-1c2f-4543-8b27-179047f6e685
# ╠═ab01db6f-09cc-4bee-8ff0-60ec1eb1a682
