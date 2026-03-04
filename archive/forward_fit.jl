### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ c77b41ac-e0c1-48fd-a987-46f7e9a9a44a
begin
	import Pkg; Pkg.activate("..")
	using Markdown, InteractiveUtils, Plots, NCDatasets
	using vSmartMOM,  vSmartMOM.Absorption
	include("/home/zhe2/FraLab/PACE_redSIF_PACE/tools.jl")
end

# ╔═╡ 87d03766-c016-4851-b01f-c65f5954d2a0
using ProgressMeter

# ╔═╡ 344e838e-1388-44ba-82e5-01849e51a471
using LinearAlgebra

# ╔═╡ 6bd228d0-340f-4a34-9453-38c79d6a57f1
using Distributions

# ╔═╡ dc0f9b7e-5218-11f0-2433-d593f48f483b
md"""
> # Test the simple forward model
> 2025/06/25

Description: generate synthetic transmission based on layered atmosphric profile, test whether principle components obtained from T/p LUT is sufficient in capturing the variation?
"""

# ╔═╡ 07631bb2-e345-4316-9be1-057e005f3ed3
md"""
> #### Read dataset and construct layered cross section
"""

# ╔═╡ 19a6460c-339d-41d9-8f1e-1b4854199b94
begin
	file = "../sample_data/atm_profile_vmrvcd.nc";
	ds   = Dataset(file)

	lat  = ds["latitude"][:]
	lon  = ds["longitude"][:]
	# just choose a profile
	myLat = 34.1377
	myLon = 118.1253
	iLat = argmin(abs.(lat .- myLat))
	iLon = argmin(abs.(lon .- myLon))
	TimeIndex = 1;

	vcd_h2o = ds["vcd_h2o"][iLon, iLat,:,TimeIndex]
	vcd_dry = ds["vcd_dry"][iLon, iLat,:,TimeIndex]
	T    = ds["temperature"][iLon, iLat,:,TimeIndex]
	p    = ds["p"][iLon, iLat,:,TimeIndex]
	close(ds)
end

# ╔═╡ d7aa862d-5473-4dfe-8735-73fea5d846e1
begin
	# make hitran table
	res   = 0.05;
	ν_min = 11627 
	ν_max = 16129
	ν = ν_min:res:ν_max;
	λ = 1 ./ ν .* 1E7;  # cm^-1 to nm
	if_wavenumber = true;

	o2_par  = Absorption.read_hitran(artifact("O2"), mol=7, iso=1, ν_min=ν_min, ν_max=ν_max);
	h2o_par = Absorption.read_hitran(artifact("H2O"), mol=1, iso=1, ν_min=ν_min, ν_max=ν_max);

	o2_voigt   = make_hitran_model(o2_par, Voigt(), wing_cutoff=10, architecture=CPU());
	h2o_voigt   = make_hitran_model(h2o_par, Voigt(), wing_cutoff=10, architecture=CPU());
end

# ╔═╡ b447ba1b-f175-4d0e-9055-8244f77660ee
begin
	# Create matrix of cross sections for each atmospheric layer
	n_layers = 72;
	cs_matrix_h2o = zeros((length(ν),n_layers))
	cs_matrix_o2 = zeros((length(ν),n_layers))
	
	# Loop over each layer 
	@showprogress for i=1:n_layers
	    p_ = p[i] / 100 # in hPa
	    T_ = T[i]
	    cs_matrix_h2o[:,i] = absorption_cross_section(h2o_voigt, ν, p_, T_);
		cs_matrix_o2[:,i] = absorption_cross_section(o2_voigt, ν, p_, T_);
	end
end

# ╔═╡ 56cfc16f-90b3-4fa8-a024-636477ba9405
begin
	plot(λ,cs_matrix_h2o[:,20:20:end])
	xlims!(650, 750)
	xlabel!("λ [nm]")
	ylabel!("σ (cm²/molec)")
	title!("xSection of H2O at different layers\n(every 20th)")
end

# ╔═╡ 69dd2de4-acc5-4667-88e4-dcbb9ba69c45
begin
	SZA   = 20.0;
	# Single air mass factor (1 way path)
	AMF   = 1/cosd(SZA);
	#VCD_h2o = dz*rho_N_h2o
	τ_h2o = cs_matrix_h2o .* vcd_h2o';
	τ_o2  = cs_matrix_o2 .* vcd_dry' .* 0.21;
	trans_h2o = exp.(-AMF*sum(τ_h2o,dims=2));
	trans_o2  = exp.(-AMF*sum(τ_o2, dims=2));
end

# ╔═╡ 33afe984-437a-4a59-9925-01ddac3448f0
begin
	# read data
	filename = "/home/zhe2/FraLab/PACE_redSIF_PACE/sample_data/PACE_OCI_RSRs.nc";
	pace = Dataset(filename, "r");

	wavlen = pace["wavelength"][:];
	RSR = pace["RSR"];
	band = pace["bands"];

	λ_ref = if_wavenumber ? reverse(λ) : λ;
	trans_h2o_ref = if_wavenumber ? reverse(trans_h2o, dims=1) : trans_h2o ;
	trans_o2_ref  = if_wavenumber ? reverse(trans_o2, dims=1) : trans_o2 ;
	trans_tot_ref = trans_h2o_ref .* trans_o2_ref;
	# println(λ_ref)

	ind₁   = findall( λ_ref[1] .< wavlen .< λ_ref[end]);
	ind₂   = findall( λ_ref[1] .< band   .< λ_ref[end]);
	λ_msr  = wavlen[ind₁];
	MyRSR  = RSR[ind₁, ind₂];
	
	MyKernel = KernelInstrument(
		band=band[ind₂],
		wvlen=λ_msr,
		RSR=RSR[ind₁, ind₂],
		wvlen_out=λ_ref
	);
	println(size(MyKernel.RSR_out))
end

# ╔═╡ dba299f3-e61a-4443-9afd-a824454757a2
begin
	# convolve
	trans_h2o_conv = MyKernel.RSR_out * trans_h2o_ref;
	trans_o2_conv  = MyKernel.RSR_out * trans_o2_ref;
	trans_tot_conv = MyKernel.RSR_out * trans_tot_ref;
	trans_conv_tot = trans_o2_conv .* trans_h2o_conv;

	plot(λ, trans_h2o, alpha=0.3, label="H2O")
	plot!(λ, trans_o2, alpha=0.3, label="O2")
	plot!(band[ind₂], trans_h2o_conv, linewidth=1., label="H2O - convolved")
	plot!(band[ind₂], trans_o2_conv, linewidth=1., label="O2 - convolved")
	plot!(band[ind₂], trans_tot_conv, linewidth=2., label="*(T1̇ T2)")
	plot!(band[ind₂], trans_conv_tot, linewidth=2., label="*(T1)̇ *(T2))")

	xlabel!("wv [nm]")
	ylabel!("transmission")
	
end

# ╔═╡ 50060b29-cbbc-4bbe-88c1-dac350110440
begin
	plot(band[ind₂], trans_conv_tot .- trans_tot_conv, alpha=1, label="diff")
	title!("when2convolve?")
end

# ╔═╡ 77f85f11-f493-4e0b-8019-b8459ce188b2
md"""
> #### fit the transmission spectra using several PCs
the target spectrum without noise is the variable "trans_tot_conv".

Forward model:

$$T=exp(∑α_{i}\phi_{i}+β_{j}\zeta_{j})$$
"""

# ╔═╡ bccb5c3e-91e3-477f-8f3a-93431e5dd7c7
begin
	h2o   = Dataset("../sample_data/H2O_transmission.nc", "r");
	o2    = Dataset("../sample_data/O2_transmission.nc", "r");
	h2o_mat = h2o["H2O_trans"][:];
	o2_mat  = o2["O2_trans"][:];
	band_pace = h2o["band"][:];

	λ_min = 620.;
	λ_max = 850.;
	ind = findall( λ_min .< band_pace .< λ_max);
	
	h2o_mat_reshape = log.(reshape(h2o_mat, size(h2o_mat, 1), :));
	o2_mat_reshape  = log.(reshape(o2_mat, size(o2_mat, 1), :));

	svd_h2o = svd(h2o_mat_reshape);
	svd_o2  = svd(o2_mat_reshape);
end

# ╔═╡ a87c6615-6370-4970-af06-ce09e3929f1a
begin
	# construct coeff. matrix
	m = 3;    # m PCs for H2O
	n = 3;    # n PCs for O2
	K = zeros(Float64, (length(band_pace), m+n)) .* NaN;
	K[:, 1:m] = svd_h2o.U[:, 1:m];
	K[:, (m+1):(m+n)] = svd_o2.U[:, 1:n];
end

# ╔═╡ 99d01ec6-3bc5-43c6-a034-cb91be92614f
x = inv(K'K)K'log.(trans_tot_conv);

# ╔═╡ 624e983f-d097-4601-97c8-f492f125fd72
begin
	plot(band[ind₂], trans_tot_conv, alpha=.7, linewidth=2., label="*(T1̇ T2)")
	plot!(band[ind₂], exp.(K * x), alpha=.7, linewidth=2., label="retrieval")

	xlabel!("wv [nm]")
	ylabel!("transmission")
	title!("H2O PC×$m, O2 PC×$n\natm. profile, SZA=$SZA")
end

# ╔═╡ e93e02a1-f3f9-4ede-9ee1-85a6aa940f1d
begin
	plot(band[ind₂], exp.(K * x)-trans_tot_conv, label="fit - true")
	title!("H2O PC×$m, O2 PC×$n")
end

# ╔═╡ 9ce9ed34-efcc-4c35-9694-0bb15866686c
md"""
> #### Transmission of an Gaussian SIF signal
"""

# ╔═╡ a93e5422-b380-4883-93d6-f0ebc7fe21c7
begin
	µ = 680;
	σ = 14;
	SIF_dist = Normal(µ, σ) ;
	
	SIF_Lw  = pdf.(SIF_dist, band[ind₂]) * 20; 
	SIF_toa = trans_tot_conv .* SIF_Lw;
end

# ╔═╡ f48fca65-91dc-433a-9c83-4898f23bef5b
begin
	plot(λ, trans_h2o, alpha=0.3, label="H2O", legend=:outerright)
	plot!(λ, trans_o2, alpha=0.3, label="O2")
	plot!(band[ind₂], trans_h2o_conv, linewidth=1., label="H2O - convolved")
	plot!(band[ind₂], trans_o2_conv, linewidth=1., label="O2 - convolved")
	plot!(band[ind₂], trans_tot_conv, linewidth=2., label="total transmission")
	plot!(band[ind₂], SIF_Lw, linewidth=2., label="water leaving SIF")
	plot!(band[ind₂], SIF_toa, linewidth=1.5, label="TOA signal")

	xlabel!("wv [nm]")
	ylabel!("transmission [unitless] / SIF [a.u.]")
end

# ╔═╡ ea285113-fff7-4251-be33-b0212ad7f3a0
md"""
> #### Forward model 2

Top of atmosphere radiance: 

$$F=exp(∑α_{i}\phi_{i}+β_{j}\zeta_{j}) × exp(∑(x_iSIF_i))$$

after log transformation:

$$logF = ∑(x_iSIF_i) + ∑α_{i}\phi_{i}+β_{j}\zeta_{j}$$

"""


# ╔═╡ d1e347d6-bc69-4ec1-ad3e-dfa3ff793856
md"""
construct SIF PCs in log-space
"""

# ╔═╡ 35cdef23-db34-434d-90d2-72ed10ce6ec8
begin
	µ_arr = 678:.05:688;
	σ_arr = 10:.1:15;
	SIF_lst = [pdf.(Normal(i, j), band[ind₂]) * 20 for i in µ_arr for j in σ_arr];
	SIF_arr = hcat(SIF_lst...)
	println("SIF spectra generated with dims: $(size(SIF_arr))")
	svd_SIF  = svd(log.(SIF_arr));
end

# ╔═╡ 303e6a38-c1c0-4edf-87d9-d0e8537712d8
begin
	k = 3;
	println(svd_SIF.S[1:k] ./ sum(svd_SIF.S) * 100)
end

# ╔═╡ 7b91862a-ccf8-4806-8d88-2d9cca77360e
begin
	plot(band[ind₂], svd_SIF.U[:,1:k], linewidth=1.5)
	title!("first $k PCs for synthetic SIF in log space")
end

# ╔═╡ 0d86fa27-93b4-4079-9022-e0e221ec1a69
md"""
> #### SIF retrieval using linearized forward model (log F)
"""

# ╔═╡ ec54d597-f344-48a4-841c-962264cd2aa9
begin
	# construct coeff. matrix
	m_h2o = 3;    # m PCs for H2O
	n_o2  = 3;    # n PCs for O2
	l_sif = 3;
	K2 = zeros(Float64, (length(band_pace), m_h2o+n_o2+l_sif)) .* NaN;
	K2[:, 1:m_h2o] = svd_h2o.U[:, 1:m_h2o];
	K2[:, (m_h2o+1):(m_h2o+n_o2)] = svd_o2.U[:, 1:n_o2];
	K2[:, (m_h2o+n_o2+1):(m_h2o+n_o2+l_sif)] = svd_SIF.U[:, 1:l_sif];
	# invert
	x2 = inv(K2'K2)K2'log.(SIF_toa);
	@show x2
end

# ╔═╡ ac88b48e-aea0-4192-9743-22133bf74eb9
begin
	plot(band[ind₂], SIF_toa,
		linewidth=2.,
		label="TOA signal", legend=:outerright)
	plot!(band[ind₂],
		exp.(K2 * x2),
		linewidth=2.,
		label="inversion - TOA")
	plot!(band[ind₂], 
		exp.(K2[:, (m_h2o+n_o2+1):(m_h2o+n_o2+l_sif)] * x2[ (m_h2o+n_o2+1):(m_h2o+n_o2+l_sif)]),
		label="inversion - water leaving SIF")

	ylims!(0,.7)
	xlabel!("wv [nm]")
	ylabel!("SIF [-]")
	title!("log-space, no noise\nH2O PC×$m_h2o, O2 PC×$n_o2, SIF PC×$l_sif")
end

# ╔═╡ 9a6ed22a-eb0f-4eab-a044-e02e1a160893
begin
	plot(band[ind₂], exp.(K2 * x2)-SIF_toa,
		label="fit - true",
		legend=:outerright
	)
	title!("log-space, no noise\nH2O PC×$m_h2o, O2 PC×$n_o2, SIF PC×$l_sif")
end

# ╔═╡ Cell order:
# ╟─dc0f9b7e-5218-11f0-2433-d593f48f483b
# ╠═c77b41ac-e0c1-48fd-a987-46f7e9a9a44a
# ╠═87d03766-c016-4851-b01f-c65f5954d2a0
# ╟─07631bb2-e345-4316-9be1-057e005f3ed3
# ╠═19a6460c-339d-41d9-8f1e-1b4854199b94
# ╠═d7aa862d-5473-4dfe-8735-73fea5d846e1
# ╠═b447ba1b-f175-4d0e-9055-8244f77660ee
# ╟─56cfc16f-90b3-4fa8-a024-636477ba9405
# ╠═69dd2de4-acc5-4667-88e4-dcbb9ba69c45
# ╠═33afe984-437a-4a59-9925-01ddac3448f0
# ╠═dba299f3-e61a-4443-9afd-a824454757a2
# ╟─50060b29-cbbc-4bbe-88c1-dac350110440
# ╟─77f85f11-f493-4e0b-8019-b8459ce188b2
# ╠═344e838e-1388-44ba-82e5-01849e51a471
# ╠═bccb5c3e-91e3-477f-8f3a-93431e5dd7c7
# ╠═a87c6615-6370-4970-af06-ce09e3929f1a
# ╠═99d01ec6-3bc5-43c6-a034-cb91be92614f
# ╠═624e983f-d097-4601-97c8-f492f125fd72
# ╠═e93e02a1-f3f9-4ede-9ee1-85a6aa940f1d
# ╟─9ce9ed34-efcc-4c35-9694-0bb15866686c
# ╠═6bd228d0-340f-4a34-9453-38c79d6a57f1
# ╠═a93e5422-b380-4883-93d6-f0ebc7fe21c7
# ╟─f48fca65-91dc-433a-9c83-4898f23bef5b
# ╟─ea285113-fff7-4251-be33-b0212ad7f3a0
# ╟─d1e347d6-bc69-4ec1-ad3e-dfa3ff793856
# ╠═35cdef23-db34-434d-90d2-72ed10ce6ec8
# ╠═303e6a38-c1c0-4edf-87d9-d0e8537712d8
# ╟─7b91862a-ccf8-4806-8d88-2d9cca77360e
# ╟─0d86fa27-93b4-4079-9022-e0e221ec1a69
# ╠═ec54d597-f344-48a4-841c-962264cd2aa9
# ╟─ac88b48e-aea0-4192-9743-22133bf74eb9
# ╟─9a6ed22a-eb0f-4eab-a044-e02e1a160893
