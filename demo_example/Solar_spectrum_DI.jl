### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ ce52e871-c888-46f7-9d50-91acd18cf0c6
import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")

# ╔═╡ 451964d2-8357-41f1-92e6-2ff98876bb3b
using DelimitedFiles, Plots

# ╔═╡ f4bef1dd-c1dd-4fbe-beef-f671af3727b3
using NCDatasets

# ╔═╡ ccd1d279-08ba-476f-845c-f798b09c5f60
using LegendrePolynomials, ForwardDiff, DiffResults, LinearAlgebra

# ╔═╡ 05c53115-d6a4-4540-b659-2a5b3442124c
using JLD2, Interpolations

# ╔═╡ ddf1db70-21ee-4184-bc49-efcbac560bad
using Einsum, Statistics

# ╔═╡ f22c2f08-e4ae-44f2-8f90-873a7b6ff069
using PACE_SIF, vSmartMOM,  vSmartMOM.Absorption

# ╔═╡ d016099e-d1ad-11f0-0d23-39cd260c5720
md"""
# Non-linearity in transmission
✍️ 2025-12-10
- include solar transmittance (DI) spectrum
- convolution before/after addition of τ
"""

# ╔═╡ 9800302f-6b8d-4586-871d-6720714c2374
md"""
#### High resolution solar spectra
"""

# ╔═╡ 61483014-4094-44c1-b0d7-49b422ed39d0
begin
	data = readdlm("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/solar_merged_20240731_600_33300_100.out", skipstart=2);
	wavenumber         = data[2:end, 1];
	trans_HighRes_wnum = data[2:end, 2];

	# wavenumber to wavelength
	wavelength   = reverse(1e7 ./ wavenumber);
	trans_HighRes_wlen = reverse(trans_HighRes_wnum);

	@info wavenumber
end

# ╔═╡ 6dc7c779-91e3-4ec0-8c3d-0a2fafabbc2f
md"""
#### High resolution xSection
"""

# ╔═╡ ae4ae48e-372d-4718-b14c-eee6ef1a46b3
function read_rescale(itp_filename::String)
	model = load_interpolation_model(itp_filename);
	ν_grid = model.ν_grid;
	p_grid = model.p_grid;
	t_grid = model.t_grid;
	itp    = model.itp;
	sitp = scale(itp, ν_grid, p_grid, t_grid);
	println("scaled! $itp_filename")
	return sitp
end

# ╔═╡ 0a9ce007-93e6-4516-92ca-439ce1fe6c28
begin
	o2_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_O2.jld2";
	o2_sitp = read_rescale(o2_jld2);
	h2o_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_H2O.jld2";
	h2o_sitp = read_rescale(h2o_jld2);

	metadata = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01.log"

	ν_grid, p_grid_hPa, t_grid = o2_sitp.ranges;
	@info ν_grid
end

# ╔═╡ 137df024-e28e-42c3-acb5-17c16435eac3
# select transmittance files
begin
	dir  = "/home/zhe2/data/MERRA2_reanalysis/";
	file = "MERRA2_400.inst6_3d_ana_Nv.20231230.nc4";

	# --- To extract upper atmospheric layers ---#
	n_layers_min = 50;
	n_layers_max = 72;
	
	# --- Here it get automated! ---
	
	# Readin ncDataset, select profiles representative of global atmosphere
	
	println("profile using: $file")
	ds   = Dataset(dir * file);
	nLon = ds.dim["lon"];
	nLat = ds.dim["lat"];
	lon  = ds["lon"][:];
	lat  = ds["lat"][:];
	T    = ds["T"][:];        # temperature
	q    = ds["QV"][:];    # specfic humidity
	psurf = ds["PS"][:];   # surface pressure
	ak   = ds.attrib["ak"][:];   # ak
	bk   = ds.attrib["bk"][:];   # bk
	close(ds)
	
	n_layers = 72;
	
	# swap
	T_swap = permutedims(T, (1, 2, 4, 3) );
	q_swap = permutedims(q, (1, 2, 4, 3) );
	p_swap = psurf;
	
	# reshape
	T_rshp = reshape(T_swap, :, n_layers);
	q_rshp = reshape(q_swap, :, n_layers);
	p_rshp = reshape(p_swap, :);
	
	# step
	MyStep = 50000;
	
	# select
	temp = T_rshp[1:MyStep:end, :];
	qw   = q_rshp[1:MyStep:end, :];
	ps   = p_rshp[1:MyStep:end];
	
	# pressure profile
	# half level pressure (Pa)
	@einsum p_half[i,j] := (ak[j] + bk[j] * ps[i]);
	p_full = (p_half[:, 2:end] + p_half[:, 1:end-1]) / 2;
	
	# rescale
	p_half /= 100;
	p_full /= 100;
	println("👍 pressure profiled rescaled to hPa!")
	
	# stdout
	println("$(size(p_full)[1]) profiles selected! 🤞")
end

# ╔═╡ 21151113-350d-4e12-9e2c-84f7d8dd18a2
begin
	# high resolution tau and trans
	# VCD => xSection
	vmr_o2  = .21
	len     = size(temp)[1]
	vmr_h2o = zeros(Float64, size(temp))
	vcd_dry = zeros(Float64, size(temp))
	vcd_h2o = zeros(Float64, size(temp))
	τ       = zeros(Float64, (len, length(ν_grid)))
	
	# generate idx
	idx = rand(n_layers_min:n_layers_max, len);
	println("number of layers to be included randomly selected from $n_layers_min to $n_layers_max, done!")	
	
	@info "Starting"
	for i in 1:len
		target_layer = idx[i];    # from the TOA down to the targeted layer
		p_half_slice = p_half[i, 1:(target_layer+1)];
		p_full_slice = p_full[i, 1:target_layer];
		T_slice = temp[i, 1:target_layer];
		q_slice = qw[i, 1:target_layer];
		
		# Apply the function
		vcd_dry_tmp, vcd_h2o_tmp, vmr_h2o_tmp = layer_VCD(
			p_half_slice, q_slice, n_layers=target_layer);
	
		if(length(p_full_slice)==length(T_slice))
			# get spectra & optical depth
			xSec_slice_o2 = [o2_sitp(ν_grid, j, k) for (j, k) in zip(p_full_slice, T_slice)];
			xSec_tmp_o2 = hcat(xSec_slice_o2...)
			τ_o2_tmp    = xSec_tmp_o2 * vcd_dry_tmp * vmr_o2;
	
			
			xSec_slice_h2o = [h2o_sitp(ν_grid, j, k) for (j, k) in zip(p_full_slice, T_slice)];
			xSec_tmp_h2o = hcat(xSec_slice_h2o...)
			τ_h2o_tmp    = xSec_tmp_h2o * vcd_h2o_tmp;
			
			# Store the result
			vcd_dry[i, 1:target_layer] = vcd_dry_tmp;
			vcd_h2o[i, 1:target_layer] = vcd_h2o_tmp;
			vmr_h2o[i, 1:target_layer] = vmr_h2o_tmp;
	
			τ[i, :]       = τ_o2_tmp .+ τ_h2o_tmp;
		else
			println("DimensionMismatch!")
		end
	
		if i % 200 == 0
			println("Processed $i / $len samples")
		end
	
	end
	@info "Completed!"
	
	
	# --- Taking into account the variations of SZA and VZA ---
	vza = [0, 15, 35, 45, 60, 75];
	AMF = 1 ./ cosd.(vza);
	# broadcast to # of profiles
	num_rep = floor(Int, len / length(AMF)) + 1;
	# repeat and truncate
	AMF_bc  = repeat(AMF, num_rep)[1:len];
	# multiply to get a slant optical depth => transmission
	trans_hres = exp.(- AMF_bc .* τ );
end

# ╔═╡ b49f4b49-f61c-4cdb-b8c8-d805b3048ae9
begin
	ν = collect(11111.0:0.01:16600.0);
	plot(
		ν, trans_hres[1,:], 
		xlims=(11111.0, 16600.0),
		size=(800, 300),
		label="High resolution",
		xlabel="wavenumber"
	)
	plot!(
		wavenumber, trans_HighRes_wnum, label="solar spec"
	)
end

# ╔═╡ 3646c5ef-75f8-405d-9fab-11ffd36fab5b
md"""
#### multiplication @ High res
"""

# ╔═╡ 76ea3abc-f907-48d7-86f7-d89cfa485334
begin
	# truncate into wavelength of interest
	λ_min = 625.0; λ_max = 860.0;
	ν_min = 1e7 ./ λ_max; ν_max = 1e7 ./ λ_min;
	# select
	ind_air = (ν .>= ν_min) .& (ν .<= ν_max);
	ind_sun = (wavenumber .>= ν_min) .& (wavenumber .<= ν_max);
	@info ν_sel   = ν[ind_air];
	@info wavenumber_sel = wavenumber[ind_sun];
	# \tau at selected wavenumbers
	τ_sun = - log.(trans_HighRes_wnum[ind_sun]);
	τ_air = (AMF_bc .* τ[:, ind_air]);
	τ_sun = reshape(τ_sun, (1, size(τ_sun)...));
	τ_tot = τ_sun .+ τ_air;
	# high res total trans
	trans_HighRes_wnum_tot = exp.( -τ_tot );
end

# ╔═╡ 442c8ada-de52-47b7-b1c5-f7ff56270381
md"""
#### Convolve @ both res :)
"""

# ╔═╡ 16200fae-68ac-46fe-b695-e3c029db0c76
begin
	# prepare data
	# pace SRF
	filename = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI/PACE_OCI_RSRs.nc";
	pace = Dataset(filename, "r");
	
	pace_wavlen = pace["wavelength"][:];
	RSR = pace["RSR"];
	band = pace["bands"];

	# convolve
	ind₁   = findall( λ_min .<= pace_wavlen .<= λ_max );
	ind₂   = findall( λ_min .<= band   .<= λ_max );
	ind₃   = findall( λ_min .<= wavelength   .<= λ_max );
	λ      = wavelength[ind₃];
	
	# get kernel instrument
	MyKernel = KernelInstrument(
		band=band[ind₂],          # band
		wvlen=pace_wavlen[ind₁],  # wavelength (to get SRF)
		RSR=RSR[ind₁, ind₂],      # [SRF dim, lowRes dim]
		wvlen_out=λ               # highRes
	);
end

# ╔═╡ c9076f82-b7fb-46cf-a93c-a78723aaa32e
begin
	# compare before and after convolution
	# multiply => convolve
	trans_HighRes_wlen_tot = reverse(trans_HighRes_wnum_tot, dims=2);
	trans_LowRes_wlen_tot_MulConv  = MyKernel.RSR_out * trans_HighRes_wlen_tot';

	# convolve => multiply
	trans_HighRes_wlen_sun  = reverse(trans_HighRes_wnum[ind_sun]);
	trans_HighRes_wvlen_air = reverse(trans_hres[:, ind_air]', dims=1);
	trans_LowRes_wlen_sun  = MyKernel.RSR_out * reverse(trans_HighRes_wnum[ind_sun]);
	trans_LowRes_wlen_air  = MyKernel.RSR_out * trans_HighRes_wvlen_air;
	trans_LowRes_wlen_tot_ConvMul = trans_LowRes_wlen_air .* trans_LowRes_wlen_sun;
end

# ╔═╡ 99e6142e-c59d-4d98-91c0-b8f1b625771c
begin
	k        = 6;
	oci_band = band[ind₂];
	abs_diff = trans_LowRes_wlen_tot_MulConv[:, k] .- trans_LowRes_wlen_tot_ConvMul[:, k];
	res_diff = abs_diff ./ trans_LowRes_wlen_tot_ConvMul[:, k] .* 100;
	
	p1 = plot()
	plot!(p1, oci_band, trans_LowRes_wlen_sun, label="solar", lw=.5, ls=:dash)
	plot!(p1, oci_band, trans_LowRes_wlen_air[:,k], label="atm", lw=.5, ls=:dash)
	plot!(p1, oci_band, trans_LowRes_wlen_tot_MulConv[:, k], label="(f·g)ₗ", lw=2)
	plot!(p1, oci_band, trans_LowRes_wlen_tot_ConvMul[:, k], label="fₗ·gₗ", lw=2)

	p2 = plot(oci_band, abs_diff, label="residual")

	p3 = plot(oci_band, res_diff, label="relative diff (%)")

	plot(p1, p2, p3, layout=(3, 1), size=(800, 800), xticks = (620:10:860, string.(620:10:860)),)
end

# ╔═╡ 77ee27de-87ec-45bf-ba94-0fbf970817d6
md"""
#### Fitting E spec using 2 forward models
"""

# ╔═╡ cddbe493-c6e3-4129-a15d-4e4570c67fa3
begin
	# load OCI
	granule_name = "sample_granule_20250808T204353_new_chl"
	path_oci = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/$(granule_name).nc"

	oci = Dataset(path_oci);
	red_band = oci["red_wavelength"][:];
	# ind_oci  = findall(λ_min .<= red_band .<= λ_max)
	E   = oci["red_solar_irradiance"]# [ind_oci];

	# bands matching up? - are they nominal band?
	@info red_band# [ind_oci]
	@info oci_band

	# interpolate E according to oci_band
	E = LinearInterpolation(red_band, E, extrapolation_bc=0).(oci_band);
end

# ╔═╡ 23d13dea-16c3-4428-b7f9-8cf6ee0d8266
function center_wavelength(λ)
    
    λ_max = ceil(maximum(λ))
    λ_min = floor(minimum(λ))
    range = (λ_max - λ_min) / 2
    λ_middle = (λ_max + λ_min) / 2
    λc = (λ .- λ_middle) ./ range
    
    return λc
end

# ╔═╡ c8fdeb6f-88d6-4fe8-8a95-528cb673a8bc
function Jacobian(
        x, 
        model, 
        len ::Int  # length of measured spectrum
    )
	res = DiffResults.JacobianResult(zeros(len), x);
	ForwardDiff.jacobian!(res, model, x);
	K   = DiffResults.jacobian(res);
	val = DiffResults.value(res);
	return K, val
end

# ╔═╡ d6f6442e-0fb0-4314-9538-7007e418c5b6
begin
	arr_len = length(E);
	Sₑ      = I(arr_len) .* 1.0;
	# edge is erroneous
	degraded_indices = [1, 2, 3, 4, arr_len-3, arr_len-2, arr_len-1, arr_len];
	for idx in degraded_indices
    	Sₑ[idx, idx] = 1e10
	end
end

# ╔═╡ d8591e6a-b25d-4f5f-929d-20251591bef0
begin
	# fitting polynoimial + solar irr
	nPoly = 5;
	# normalizing wavelength
	λc_LowRes  = center_wavelength(oci_band);
	# LegendrePolynomials vector
	v_LowRes   = collectPl.(λc_LowRes, lmax=nPoly);
	# invert
	inv_Sₑ = inv(Sₑ);
	K₀ = hcat(v_LowRes...)' .* trans_LowRes_wlen_sun;
	G₀ = inv( K₀'*inv_Sₑ*K₀ )K₀'inv_Sₑ;
	x₀ = G₀ * E;
	# reconstruct
	E_LowRes_rec = hcat(v_LowRes...)' * x₀ .* trans_LowRes_wlen_sun;
end

# ╔═╡ 1d044a95-4527-401a-b9ff-2a71f6e53f52
function forward_model_HighResConv_sun(x)
	# centering wavelength
	λc_HighRes = center_wavelength(λ);
	# HighRes reflectance
	v_HighRes  = collectPl.(λc_HighRes, lmax=nPoly);
	ρ_HighRes  = hcat(v_HighRes...)' * x[1 : nPoly+1];
	# HighRes trans
	E_HighRes  = ρ_HighRes .* trans_HighRes_wlen_sun;
	# convolve
	return MyKernel.RSR_out * E_HighRes
end
	

# ╔═╡ dc93e6a2-d8b2-4489-8477-98019c80009e
function MyIter(xₐ, Myforwar_model = forward_model_HighResConv_sun)
	# iteratively get estimation of highres
    # initial
    xₙ  = copy(xₐ);
    Kₙ, yₙ  = Jacobian(xₙ, Myforwar_model, arr_len);
    RMSE₀   = 1e20; 
    RMSE₁   = root_mean_square(E, yₙ);
    ΔRMSE   = RMSE₁ - RMSE₀;
	iter_label = 1;

	while ( abs(ΔRMSE) > 1e-6 ) & ( iter_label < 20 )
        # k += 1
        # get Gain matrix
        Gₙ     = inv( K₀'*inv_Sₑ*K₀ )K₀'*inv_Sₑ;
        # retrieval
        xₙ₊₁   = xₐ .+ Gₙ * (E .- yₙ .+ Kₙ * ( xₙ .- xₐ ) );
        # update x and y
		# @show "previous", xₙ
        Kₙ₊₁, yₙ₊₁ = Jacobian(xₙ₊₁, Myforwar_model, arr_len);
		xₙ     = xₙ₊₁;
        yₙ     = yₙ₊₁;
        Kₙ     = Kₙ₊₁;
		# @show "updated", xₙ
        # iter ++
        iter_label += 1;
        # test convergence
        RMSE₀  = RMSE₁;
        RMSE₁  = root_mean_square(E, yₙ);
        ΔRMSE  = RMSE₁ - RMSE₀;
    end

	print("number of iterations: $iter_label")

	return xₙ, yₙ
end

# ╔═╡ 3026b16a-b52c-4438-8bdd-22b7e987ae3f
x₀_HighRes, E_HighRes_rec = MyIter(Float64.(x₀));

# ╔═╡ f3b20e9a-9888-4171-b582-b1e590ce48f1
@info x₀, x₀_HighRes

# ╔═╡ da5a3f37-950a-48c2-bc46-bdc26014ebe3
begin
	plot(
		oci_band, E, 
		label="E - TSIS",
		size=(800, 300), 
		xticks = (620:10:860, string.(620:10:860))
	)
	plot!(
		oci_band, E_LowRes_rec, label="LowRes Retrieval, nPoly=$nPoly", color=:orange
	)
	plot!(
		oci_band, E_HighRes_rec, label="High Retrieval, nPoly=$nPoly", color=:purple
	)
	
end

# ╔═╡ 54650f5d-e8d7-406e-bd87-549fed04374b
# plot residual
begin
	p = plot(
		oci_band, E .- E_LowRes_rec, 
		label="LowRes Retrieval Residual",
		size=(800, 300), 
		xticks = (620:10:860, string.(620:10:860)),
		color=:orange,
		title="nPoly=$nPoly",
		ylims=(-4.0, 4.0)
	)
	plot!(p,
		oci_band, E .- E_HighRes_rec, label="HighRes Retrieval Residual", color=:purple
	)
	
	# Shade the degraded bands
	for idx in degraded_indices
		band_start = oci_band[idx] - 1.0
		band_end = oci_band[idx] + 1.0
		plot!(p, [band_start, band_end], [minimum(E .- E_LowRes_rec), maximum(E .- E_LowRes_rec)], color=:gray, alpha=0.5, label=false, lw=0)
	end
	p
end


# ╔═╡ 830ef1bb-ef76-4ef1-aeeb-02327d3e61e6
md"""
#### Same idea, test if high resolution optical depth is well captured by SVD?
"""

# ╔═╡ 0908a1a3-ccfc-4203-8a5b-909f0a40646f
begin
	# construct some two way transmittance in HighRes space
	k1 = 7; k2 = 7;
	τ₁ = τ_air[k1,:];
	τ₂ = τ_air[k2,:];
	τ_TwoWay = @. τ₁ + τ₂;
	# adding -> convolve
	trans_HighRes_wlen_TwoWay = @. exp( - τ_TwoWay );
	trans_LowRes_wlen_TwoWay_MulConv  = MyKernel.RSR_out * reverse(trans_HighRes_wlen_TwoWay);
	# convolve -> adding
	trans_LowRes_wlen_TwoWay_ConvMul  = @. trans_LowRes_wlen_air[:, k1] * trans_LowRes_wlen_air[:, k2];
	T₁ = MyKernel.RSR_out * reverse(exp.(-τ₁));
	T₂ = MyKernel.RSR_out * reverse(exp.(-τ₂));
	trans_LowRes_wlen_TwoWay_ConvMul1 = T₁ .* T₂;

	# plot
	p_twoway = plot(
		xticks = (620:10:860, string.(620:10:860)),
		size=(800, 300), title="Two-way transmittance"
	)
	plot!(p_twoway, oci_band, trans_LowRes_wlen_TwoWay_MulConv, label="(f·g)ₗ", lw=2)
	plot!(p_twoway, oci_band, trans_LowRes_wlen_TwoWay_ConvMul, label="fₗ·gₗ")
	plot!(p_twoway, oci_band, trans_LowRes_wlen_TwoWay_ConvMul1, label="fₗ·gₗ\n(same, just for double check)", color=:red, lw=2)	
	p_twoway
	
end

# ╔═╡ 51e0f2a4-1add-4ff0-8531-e9fc62c200df
begin
	# SVD to 
	path_transmittance_summer = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc"
	path_transmittance_winter = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc"
	
	println("Loading data...")
	
	# Load transmittance data
	trans_spec, bands = Dataset(path_transmittance_summer) do summer
	    Dataset(path_transmittance_winter) do winter
	        (cat(summer["transmittance"][:, :], winter["transmittance"][:, :], dims=1),
	         summer["band"][:])
	    end
	end
	println("Transmittance: $(size(trans_spec, 1)) spectra")
end

# ╔═╡ fe4108a2-9e32-4064-b015-2142908895f5
begin
	if_log = true;
	# res    = Spectral_NMF(
	# 			trans_spec, bands, 
	# 			Float64.(collect(skipmissing(oci_band))); 
 #               rank=15, if_log=if_log);
	res    = Spectral_SVD(
				Float64.(trans_spec'), 
				bands, 
				Float64.(collect(skipmissing(oci_band))), 
				if_log=if_log
			);
end

# ╔═╡ 9af96fee-497e-4bb6-9201-28dfabb573ee
function SVD_fit(PrinComp, y)
	K₀ = PrinComp;
	G₀ = inv( K₀'*K₀ )K₀';
	x₀ = G₀ * y;
	return x₀
end

# ╔═╡ d3642b1d-0206-4e28-88ec-f6190d1a2cdf
begin
	# y is log.(trans)
	nPC1 = 10; nPC2 = 15;
	PrinComp = res.PrinComp;
	
	x₁_MulConv = SVD_fit(PrinComp[:, 1:nPC1], -log.(trans_LowRes_wlen_TwoWay_MulConv));
	x₂_MulConv = SVD_fit(PrinComp[:, 1:nPC2], -log.(trans_LowRes_wlen_TwoWay_MulConv));

	x₁_ConvMul = SVD_fit(PrinComp[:, 1:nPC1], -log.(trans_LowRes_wlen_TwoWay_ConvMul1));
	x₂_ConvMul = SVD_fit(PrinComp[:, 1:nPC2], -log.(trans_LowRes_wlen_TwoWay_ConvMul1));

	# reconstruction
	spec_rec_MulConv1 = exp.(-PrinComp[:, 1:nPC1] * x₁_MulConv)
	spec_rec_MulConv2 = exp.(-PrinComp[:, 1:nPC2] * x₂_MulConv)

	spec_rec_ConvMul1 = exp.(-PrinComp[:, 1:nPC1] * x₁_ConvMul)
	spec_rec_ConvMul2 = exp.(-PrinComp[:, 1:nPC2] * x₂_ConvMul)

end

# ╔═╡ ae36b3bc-aa90-4f52-8041-7784be92086b
begin
	p_rec = plot(
		xticks = (620:10:860, string.(620:10:860)),
		title  = "log-SVD"
	)
	plot!(p_rec, oci_band, trans_LowRes_wlen_TwoWay_MulConv, label="(f·g)ₗ", ls=:dash, lw=2)
	plot!(p_rec, oci_band, trans_LowRes_wlen_TwoWay_ConvMul1, label="fₗ·gₗ", ls=:dash, lw=2)	
	# recons
	plot!(p_rec, oci_band, spec_rec_MulConv1, label="(f·g)ₗ - recon. w/ nPC=$nPC1")
	plot!(p_rec, oci_band, spec_rec_MulConv2, label="(f·g)ₗ - recon. w/ nPC=$nPC2")
	plot!(p_rec, oci_band, spec_rec_ConvMul1, label="fₗ·gₗ - recon. w/ nPC=$nPC1")
	plot!(p_rec, oci_band, spec_rec_ConvMul2, label="fₗ·gₗ - recon. w/ nPC=$nPC2")

	# residual
	p_resd = plot(xticks = (620:10:860, string.(620:10:860)),)
	plot!(
		p_resd, oci_band,
		trans_LowRes_wlen_TwoWay_MulConv .- spec_rec_MulConv1,
		label="(f·g)ₗ - recon. w/ nPC=$nPC1"
	)
	plot!(
		p_resd, oci_band,
		trans_LowRes_wlen_TwoWay_MulConv .- spec_rec_MulConv2,
		label="(f·g)ₗ - recon. w/ nPC=$nPC2"
	)
	plot!(
		p_resd, oci_band,
		trans_LowRes_wlen_TwoWay_ConvMul1 .- spec_rec_ConvMul1,
		label="fₗ·gₗ - recon. w/ nPC=$nPC1"
	)
	plot!(
		p_resd, oci_band,
		trans_LowRes_wlen_TwoWay_ConvMul1 .- spec_rec_ConvMul2,
		label="fₗ·gₗ - recon. w/ nPC=$nPC2"
	)

	# combine
	plot(p_rec, p_resd, layout=(2,1), size=(1000, 600))
end

# ╔═╡ cbc396ae-ec5c-4325-bb60-af2a572df1d1
md"""
#### Pseudo sim of (solar+atmospheric) transmittance 

- forward model: multiplication of (T_sun and T_air) at HighRes, then convolve  

- Fitting: T_sun @ LowRes + SVD
"""

# ╔═╡ 479c2a80-1ad7-4521-9caa-a02a8cf47979
function Solar_SVD_fit(PrinComp, y)
	K₀ = PrinComp;
	G₀ = inv( K₀'*K₀ )K₀';
	x₀ = G₀ * (y .+ log.(trans_LowRes_wlen_sun));
	return x₀
end

# ╔═╡ 7f9bfdda-2903-4157-ac0e-ba445a005bd7
begin
	# forward construction - done from: trans_LowRes_wlen_tot_MulConv
	# fitting: convolved solar ref spec + SVD
	# y should still be -log.(trans), notice the negative 
	k3 = 12;
	x₁_SolarSVD = Solar_SVD_fit(
		PrinComp[:, 1:nPC1], -log.(trans_LowRes_wlen_tot_MulConv[:, k3])
	);
	x₂_SolarSVD = Solar_SVD_fit(
		PrinComp[:, 1:nPC2], -log.(trans_LowRes_wlen_tot_MulConv[:, k3])
	);
	# recon
	spec_rec₁_SolarSVD = trans_LowRes_wlen_sun .* exp.(- PrinComp[:, 1:nPC1] * x₁_SolarSVD);
	spec_rec₂_SolarSVD = trans_LowRes_wlen_sun .* exp.(- PrinComp[:, 1:nPC2] * x₂_SolarSVD);
end

# ╔═╡ d0eb07f2-c60a-4fe4-ac6b-21478bb8bbea
begin
	p_rec1 = plot(
		xticks = (620:10:860, string.(620:10:860)),
		title  = "log-SVD, k=$k3"
	)
	plot!(p_rec1, oci_band, trans_LowRes_wlen_tot_MulConv[:, k3], label="(f·g)ₗ", ls=:dash, lw=2)
	# recons
	plot!(p_rec1, oci_band, spec_rec₁_SolarSVD, label="SVD recon. w/ nPC=$nPC1")
	plot!(p_rec1, oci_band, spec_rec₂_SolarSVD, label="SVD recon. w/ nPC=$nPC2")
	# ref spec
	plot!(p_rec1, oci_band, trans_LowRes_wlen_sun, label="solar", lw=.5, ls=:dash)
	plot!(p_rec1, oci_band, trans_LowRes_wlen_air[:, k3], label="atm", lw=.5, ls=:dash)

	# residual
	p_resd1 = plot(xticks = (620:10:860, string.(620:10:860)),)
	plot!(
		p_resd1, oci_band,
		trans_LowRes_wlen_tot_MulConv[:, k3] .- spec_rec₁_SolarSVD,
		label="(f·g)ₗ - recon. w/ nPC=$nPC1"
	)
	plot!(
		p_resd1, oci_band,
		trans_LowRes_wlen_tot_MulConv[:, k3] .- spec_rec₂_SolarSVD,
		label="(f·g)ₗ - recon. w/ nPC=$nPC2"
	)

	# combine
	plot(p_rec1, p_resd1, layout=(2,1), size=(1000, 600))
end

# ╔═╡ 0e7914dc-5fc0-4a99-82e1-70000d4127b5
begin
	# a scatter plot: x - solar, y - atm, c - residual
	residual = trans_LowRes_wlen_tot_MulConv[:, k3] .- spec_rec₂_SolarSVD;
	scatter(
		trans_LowRes_wlen_sun, trans_LowRes_wlen_air[:, k3],
	    zcolor=residual,
	    marker=:circle,
	    markersize=5,
	    xlabel="Solar trans.",
	    ylabel="Atmospheric trans.",
	    # title="Retrieval: Solar vs Atmospheric",
	    colorbar_title="Residual",
	    legend=false,
	    size=(300, 300),
	    color=:viridis,
		markerstrokewidth=0,
	)
end

# ╔═╡ 907d4c5f-4247-4fdc-8b7f-4fd11cb1268b


# ╔═╡ Cell order:
# ╟─d016099e-d1ad-11f0-0d23-39cd260c5720
# ╠═451964d2-8357-41f1-92e6-2ff98876bb3b
# ╠═f4bef1dd-c1dd-4fbe-beef-f671af3727b3
# ╠═ccd1d279-08ba-476f-845c-f798b09c5f60
# ╠═ce52e871-c888-46f7-9d50-91acd18cf0c6
# ╠═05c53115-d6a4-4540-b659-2a5b3442124c
# ╠═ddf1db70-21ee-4184-bc49-efcbac560bad
# ╠═f22c2f08-e4ae-44f2-8f90-873a7b6ff069
# ╟─9800302f-6b8d-4586-871d-6720714c2374
# ╠═61483014-4094-44c1-b0d7-49b422ed39d0
# ╟─6dc7c779-91e3-4ec0-8c3d-0a2fafabbc2f
# ╠═ae4ae48e-372d-4718-b14c-eee6ef1a46b3
# ╠═0a9ce007-93e6-4516-92ca-439ce1fe6c28
# ╠═137df024-e28e-42c3-acb5-17c16435eac3
# ╟─21151113-350d-4e12-9e2c-84f7d8dd18a2
# ╟─b49f4b49-f61c-4cdb-b8c8-d805b3048ae9
# ╟─3646c5ef-75f8-405d-9fab-11ffd36fab5b
# ╠═76ea3abc-f907-48d7-86f7-d89cfa485334
# ╟─442c8ada-de52-47b7-b1c5-f7ff56270381
# ╠═16200fae-68ac-46fe-b695-e3c029db0c76
# ╠═c9076f82-b7fb-46cf-a93c-a78723aaa32e
# ╠═99e6142e-c59d-4d98-91c0-b8f1b625771c
# ╟─77ee27de-87ec-45bf-ba94-0fbf970817d6
# ╠═cddbe493-c6e3-4129-a15d-4e4570c67fa3
# ╟─23d13dea-16c3-4428-b7f9-8cf6ee0d8266
# ╟─c8fdeb6f-88d6-4fe8-8a95-528cb673a8bc
# ╟─1d044a95-4527-401a-b9ff-2a71f6e53f52
# ╠═d6f6442e-0fb0-4314-9538-7007e418c5b6
# ╠═d8591e6a-b25d-4f5f-929d-20251591bef0
# ╠═dc93e6a2-d8b2-4489-8477-98019c80009e
# ╠═3026b16a-b52c-4438-8bdd-22b7e987ae3f
# ╠═f3b20e9a-9888-4171-b582-b1e590ce48f1
# ╟─da5a3f37-950a-48c2-bc46-bdc26014ebe3
# ╟─54650f5d-e8d7-406e-bd87-549fed04374b
# ╟─830ef1bb-ef76-4ef1-aeeb-02327d3e61e6
# ╠═0908a1a3-ccfc-4203-8a5b-909f0a40646f
# ╟─51e0f2a4-1add-4ff0-8531-e9fc62c200df
# ╠═fe4108a2-9e32-4064-b015-2142908895f5
# ╠═9af96fee-497e-4bb6-9201-28dfabb573ee
# ╠═d3642b1d-0206-4e28-88ec-f6190d1a2cdf
# ╟─ae36b3bc-aa90-4f52-8041-7784be92086b
# ╟─cbc396ae-ec5c-4325-bb60-af2a572df1d1
# ╠═479c2a80-1ad7-4521-9caa-a02a8cf47979
# ╠═7f9bfdda-2903-4157-ac0e-ba445a005bd7
# ╟─d0eb07f2-c60a-4fe4-ac6b-21478bb8bbea
# ╟─0e7914dc-5fc0-4a99-82e1-70000d4127b5
# ╠═907d4c5f-4247-4fdc-8b7f-4fd11cb1268b
