### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 52c5cc8a-620b-4132-9586-cf061ef3ee6e
begin
	import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE");
	using vSmartMOM,  vSmartMOM.Absorption
	using JLD2
	using NCDatasets, Einsum, Statistics, Random
	using Parameters, ProgressMeter, ProgressLogging, Base.Threads
	using PACE_SIF
end

# ╔═╡ 2d8f0c3f-1512-479f-99c3-dfeee855f70f
using Interpolations

# ╔═╡ 9dcf34aa-4b92-4000-936f-a0a5146c7c08
using LinearAlgebra

# ╔═╡ 67db2fd7-e281-41c3-89c5-ec396c166ee1
using Distributions

# ╔═╡ 9f25f8bc-4085-4cff-bbf9-1f6ed3e7197c
using Plots

# ╔═╡ 8b35e192-f134-11f0-3097-37bfabf9cf36
md"""
### Distributed two-way transmittance
----
Created: 2026-01-14


- Step 1: Chi-distribution
- Step 2: Load air profiles
- Step 3: Generate layer-resolved optical depth for O2 (same for H2O)
- Step 4: Weighted by distributions
- Step 5: Take the exponential, convolve => compare the exponentials (transmittance)
- Step 6: Does RMS of SVD fit change with distributions?

"""

# ╔═╡ 28da7053-0ead-49b9-9bec-902f38c0e933
md"""
### Explore distributions
"""

# ╔═╡ 47d83de9-43f4-4d75-8665-3d187fc5fddc
begin
	# degree of freedom: nu
	ν = 4;
	# noncentral param: \lambda
	# λ = 0;
	# generate distribution
	distribution = Chisq(ν);
	# pdf of a distribution
	n_layer = 72;
	x       = 1:n_layer;
	y_pdf   = pdf.(distribution, x);
	# to normalize the distribution
	@show sum_of_y = sum(y_pdf)
	y_pdf = y_pdf ./ sum_of_y;
	# plot
	p1 = plot(
		x, y_pdf,
		size=(800, 300),
		margin=8Plots.mm
		);
	xlabel!("Layer number")
	ylabel!("% of sunlight being \n reflected from the layer")
	p1
end

# ╔═╡ 5fa9abaa-250a-432e-b21f-1d51056ebd5b
begin
	# generate noncentral chi square
	ν₂ = 1; λ₂ = 5;
	distribution2 = NoncentralChisq(ν₂, λ₂);
	y_pdf_2 = pdf.(distribution2, x);
	@show sum_of_y2 = sum(y_pdf_2);
	y_pdf_2 = y_pdf_2 ./ sum_of_y2;
	p11 = plot(
		x, y_pdf_2,
		size=(800, 300),
		margin=8Plots.mm, label="non-central Chi",
		);
	plot!(p11, x, y_pdf,
		size=(800, 300),
		label="Chi distribution"
	)
	xlabel!("Layer number")
	ylabel!("% of sunlight being \n reflected from the layer")
	p11
end

# ╔═╡ 8c10e29a-f461-4551-bbff-6103d8a3a4e7
# set which distribution to be used
y_pdf_target = y_pdf_2;

# ╔═╡ 0e887dd1-a554-4c29-9035-8633b4482c5b
md"""
### Generate profiles
"""

# ╔═╡ 72d19300-14cb-41e7-8416-f632586c09f5
function read_rescale(itp_filename::String)
	model = load_interpolation_model(itp_filename);
	ν_grid = model.ν_grid;
	p_grid = model.p_grid;
	t_grid = model.t_grid;
	itp    = model.itp;
	sitp   = Interpolations.scale(itp, ν_grid, p_grid, t_grid);
	println("scaled! $itp_filename")
	return sitp
end

# ╔═╡ fb724321-762b-4bac-8e45-2d47d23e56c8
begin
	# --- reanalysis data used --- #
	dir  = "/home/zhe2/data/MERRA2_reanalysis/";
	file = "MERRA2_400.inst6_3d_ana_Nv.20231230.nc4";
	
	# --- specify the output NetCDF filename ---
	output_filename = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Nov23.nc";
	
	o2_jld2  = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_O2.jld2";
	o2_sitp  = read_rescale(o2_jld2);
	h2o_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_H2O.jld2";
	h2o_sitp = read_rescale(h2o_jld2);
	
	metadata = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01.log"
	
	ν_grid, p_grid_hPa, t_grid = o2_sitp.ranges;
end

# ╔═╡ 7677f2e9-6262-45ba-804c-5e8c19a79c87
begin
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
	MyStep = 60000;
	
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

# ╔═╡ 87b63aea-4852-405e-9a4a-47be24f0c784
begin
	# visualize temperature profile
	k = 1;
	# pressure goes from 0-1000 hPa
	p2 = plot(x, p_full[k,:],
			xlabel="Layer number",
		    ylabel="Pressure",
		    margin=8Plots.mm
			)
	p3 = plot(x, temp[k,:],
			xlabel="Layer number",
		    ylabel="Temperature",
		    size=(800, 500),
		    margin=8Plots.mm
	)
	plot(p2, p3, layout=(2,1))
end

# ╔═╡ c1102d6a-b729-444b-9fd7-e90e66d2331b
begin
	# pressure against temperature
	p4 = plot(
		temp[k,:], p_full[k,:],
		yflip=true,
		size=(300, 450),
		margin=8Plots.mm,
		xlabel="temperature",
		ylabel="pressure"
	)
	# reverse y_axis
	
end

# ╔═╡ cb7227d2-3e4e-42eb-b208-6573a4aff59e
md"""
### Optical depth (high resolutions)
"""

# ╔═╡ 53a3ea6e-cfa9-4d54-ab39-b82aa48eaa2c
begin
	# only using one profile for testing
	# 1. Set the fixed index
	ind = 13;
	
	# 2. Define constants
	vmr_o2 = 0.21;
	
	# 3. Extract slices (ensure these variables are defined in other cells)
	# These will now update together whenever p_half, temp, etc., change.
	p_half_slice = p_half[ind, :]
	p_full_slice = p_full[ind, :]
	T_slice      = temp[ind, :]
	q_slice      = qw[ind, :]

	# Apply the function: vertical column density
	vcd_dry_tmp, vcd_h2o_tmp, vmr_h2o_tmp = layer_VCD(
		p_half_slice, q_slice);

	# optical depth
	xSec_slice_o2  = [o2_sitp(ν_grid, j, k) for (j, k) in zip(p_full_slice, T_slice)];
	
	xSec_tmp_o2 = hcat(xSec_slice_o2...);
	
	xSec_slice_h2o = [h2o_sitp(ν_grid, j, k) for (j, k) in zip(p_full_slice, T_slice)];
	
	xSec_tmp_h2o = hcat(xSec_slice_h2o...);
end

# ╔═╡ f04edcf6-20a5-4f89-8a55-79ce2637313e
begin

	# show shapes
	@show size(xSec_tmp_o2)   # cross section
	@show size(vcd_dry_tmp * vmr_o2)   # layer-resolved vertical column density
	
	# get optical depth of each layer for o2
	vcd_o2_tmp = vcd_dry_tmp * vmr_o2;
	τ_o2       = xSec_tmp_o2 .* vcd_o2_tmp';

	# get optical depth of each layer for h2o
	τ_h2o      = xSec_tmp_h2o .* vcd_h2o_tmp';

	# order of magnitude
	idx = 1:10:72;
	plt_list = [
	    plot([τ_o2[:, k], τ_h2o[:, k]],
	         title = "Pressure: $(p_full_slice[k]) hPa", 
		     titlefontsize=8,
	         legend = false) 
	    for k in idx
	]

	plot(
		plt_list..., layout = (length(idx), 1), size = (600, 1000),
		plot_title = "optical depth per layer",
	    plot_titlevspan = 0.05,      # Allocates 5% of height to the title area
	    margin = .5Plots.mm          # Adds breathing room around panels
	)
end

# ╔═╡ 3ad6f119-2bba-45b7-a05b-0133a5522e37
begin
	# accumulated optical depth (cumsum) from TOA=0
	τ_o2_acc  = cumsum(τ_o2, dims=2);
	τ_h2o_acc = cumsum(τ_h2o, dims=2);
	# test whether the last term is same as doing a matrix product - yayy!
	@show mean(xSec_tmp_o2 * vcd_dry_tmp * vmr_o2)
	@show mean(τ_o2_acc[:, end])

	# hey
	plt_list1 = [
		plot([τ_o2_acc[:, k], τ_h2o_acc[:,k]],
			 title = "Pressure: $(p_full_slice[k]) hPa", 
			 titlefontsize=8,
			 legend = false) 
		for k in idx
		]

	plot(
		plt_list1..., layout = (length(idx), 1), size = (600, 1000),
		plot_title = "optical depth over the layers above",
	    plot_titlevspan = 0.05, # Allocates 5% of height to the title area
	    margin = .5Plots.mm          # Adds breathing room around panels
	)
	
end

# ╔═╡ 94a5b030-4a24-4529-9d14-44bb67c7b863
begin
	# weight by pdf
	weights  = reverse(y_pdf_target);
	τ_o2_tot_weight  = τ_o2_acc * weights
	τ_o2_tot  = τ_o2_acc[:, end];
	τ_h2o_tot_weight = τ_h2o_acc * weights;
	τ_h2o_tot = τ_h2o_acc[:, end];
	
	# compare weighted and unweighted
	plot(
		[
			τ_o2_tot_weight, τ_o2_tot,
			τ_h2o_tot_weight, τ_h2o_tot
		], 
		label=[
			"weighted - O₂" "unweighted - O₂" "weighted - H₂O" "unweighted - H₂O"
		],
		size=(800, 300)
	)
end

# ╔═╡ 62470600-a00d-4d73-9188-36ec778e436b
begin
	# Spectral response function
	# read data
	filename = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI/PACE_OCI_RSRs.nc";
	pace = Dataset(filename, "r");
	
	wavlen = pace["wavelength"][:];
	RSR = pace["RSR"];
	band = pace["bands"];
	
	if_wavenumber = true;
	
	λ = if_wavenumber ? reverse(1 ./ collect(ν_grid) .* 1E7) : (1 ./ collect(ν_grid) .* 1E7);
	
	ind₁   = findall( λ[1] .< wavlen .< λ[end]);
	ind₂   = findall( λ[1] .< band   .< λ[end]);
	λ_msr  = wavlen[ind₁];
	MyRSR  = RSR[ind₁, ind₂];
	
	MyKernel = KernelInstrument(
		band=band[ind₂],
		wvlen=λ_msr,
		RSR=RSR[ind₁, ind₂],
		wvlen_out=λ
	);
	println(size(MyKernel.RSR_out));

end

# ╔═╡ 483d44e8-89b6-43f1-8321-7906a059ac19
begin
	# compare two-way and 1.x-way transmittance: 
	# 1) optical depth should be doubled for both cases
	τ_o2_tot_weight2 = 2 * τ_o2_tot_weight;
	τ_o2_tot2        = 2 * τ_o2_tot;
	τ_h2o_tot_weight2 = 2 * τ_h2o_tot_weight;
	τ_h2o_tot2        = 2 * τ_h2o_tot;
	# 2) take the exponential term
	trans_hres_weighted = exp.( - τ_o2_tot_weight2);
	trans_hres          = exp.( - τ_o2_tot2);
	trans_hres_h2o_weighted = exp.( - τ_h2o_tot_weight2);
	trans_hres_h2o          = exp.( - τ_h2o_tot2);
		# 🌟 add them up (@ higher resolution)
	trans_hres_tot_weighted = exp.( - τ_o2_tot_weight2 .- τ_h2o_tot_weight2);
	trans_hres_tot          = exp.( - τ_o2_tot2 .- τ_h2o_tot2 );
		# also one-way
	trans_hres_tot_1        = exp.( - τ_o2_tot .- τ_h2o_tot );
	# 3) convolve to low resolution
	# 3.1) reserved order: wavenumber -> wavelength
	if if_wavenumber
		trans_revs = reverse(trans_hres);
		trans_revs_weighted = reverse(trans_hres_weighted);
		trans_revs_h2o = reverse(trans_hres_h2o);
		trans_revs_h2o_weighted = reverse(trans_hres_h2o_weighted);
		# 🌟 reverse hres_tot
		trans_revs_tot = reverse(trans_hres_tot);
		trans_revs_tot_weighted = reverse(trans_hres_tot_weighted);
		# one-way
		trans_revs_tot_1 = reverse(trans_hres_tot_1);
	end
	# 3.2) convolve
	trans_lres              = MyKernel.RSR_out * trans_revs;
	trans_lres_weighted     = MyKernel.RSR_out * trans_revs_weighted;
	trans_lres_h2o          = MyKernel.RSR_out * trans_revs_h2o;
	trans_lres_h2o_weighted = MyKernel.RSR_out * trans_revs_h2o_weighted;
		# 🌟
	trans_lres_tot          = MyKernel.RSR_out * trans_revs_tot;
	trans_lres_tot_weighted = MyKernel.RSR_out * trans_revs_tot_weighted;
		# one-way
	trans_lres_tot_1        = MyKernel.RSR_out * trans_revs_tot_1;
	# 4) visualize
	plot(
		[
			trans_lres_weighted, trans_lres,
			trans_lres_h2o_weighted, trans_lres_h2o,
		], 
		label=[
			"weighted - O₂" "unweighted - O₂" "weighted - H₂O" "unweighted - H₂O"
		],
		size=(800, 300),
		lw=1, linestyle=:dash,
		title="O₂ + H₂O, Chi-distribution (ν=$ν₂, λ=$λ₂)"
	)
	plot!(
		[trans_lres_tot_weighted, trans_lres_tot],
		label=["weighted - total" "unweighted - total"],
		lw=2,
	)
end

# ╔═╡ 6e6642c7-20e2-4709-8193-bd93cb74b4eb
begin
	# obtain sample transmittance
	DecompositionMethod = :SVD;  
	path_transmittance_summer = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc";
	path_transmittance_winter = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc";
	trans, bands = Dataset(path_transmittance_summer) do summer
	    Dataset(path_transmittance_winter) do winter
	        (cat(summer["transmittance"][:, :], winter["transmittance"][:, :], dims=1),
	         summer["band"][:])
	    end
	end
end

# ╔═╡ ca1e324e-c8c3-49b5-a4af-38740afe4f26
begin
	oci_band = band[ind₂];
	if_log   = true;
	nPC      = 15;
	
	loading_ave_trans, loading_var_trans, cov_matx, PrinComp = if DecompositionMethod == :NMF
	    res = Spectral_NMF(trans, bands, Float64.(collect(skipmissing(oci_band))); 
	                       rank=ranks, if_log=if_log)
	    ([mean(res.Loading[:, i]) for i in 1:ranks],
	     [var(res.Loading[:, i]) for i in 1:ranks],
	     cov(res.Loading, dims=1),
	     res.PrinComp')
	elseif DecompositionMethod == :SVD
	    res = Spectral_SVD(Float64.(trans'), bands, 
	                       Float64.(collect(skipmissing(oci_band))), if_log=if_log)
	    ([mean(res.Loading[i, :]) for i in 1:nPC],
	     [var(res.Loading[i, :]) for i in 1:nPC],
	     cov(res.Loading[1:nPC, :], dims=2),
	     res.PrinComp[:, 1:nPC])
	end
end

# ╔═╡ 5c4a20c9-f00d-4197-970d-8bec0935d5db
begin
	# visualize the shape of first nPC principle components
	plot(
		oci_band, res.PrinComp[:, 1:nPC],
		size=(800, 300),
		xlabel="wavelength",
		margin=8Plots.mm
	)
end

# ╔═╡ 7fedef1f-b0e1-480f-814b-b01dddfd62eb
begin
	# OLS fit
	params_weighted = (PrinComp' * PrinComp) * (PrinComp' * log.(trans_lres_tot_weighted) );
	trans_lres_tot_fit = exp.( PrinComp * params_weighted );
	# visualize, get residual
	residual_tot = trans_lres_tot_weighted .- trans_lres_tot_fit;
	@show dot(residual_tot, residual_tot)
	# plot
	p5 = plot(
		oci_band, trans_lres_tot_weighted, size=(800, 200), 
		label="original - weighted (up+down)", lw=2
	)
	# one-way
	plot!(p5, oci_band, trans_lres_tot_1, label="up")
	# fitted
	plot!(p5,
		oci_band, trans_lres_tot_fit, label="fit"
	)
end

# ╔═╡ 8154d4d2-7ba1-461b-b158-d24435b49c1d
md"""
what matters is the ratio of one-way and two-way transmittance.

How does the ratio change with different "light distribution"?
"""

# ╔═╡ 7ebf8bc8-8ef3-4c9a-b208-da5bb54f56e3
begin
	# take the ratio
	ratio = trans_lres_tot_weighted ./ trans_lres_tot_1
	# plot
	plot(
		oci_band, ratio, 
		label="two-way / one-way ratio",
		size=(600, 150)
	)
end

# ╔═╡ e132db67-7f4e-4fd6-a686-7920fb1c4a74
begin
	# take the ratio （in log space）
	ratio_log_weighted = log.(trans_lres_tot_weighted) ./ log.(trans_lres_tot_1)
	ratio_log          = log.(trans_lres_tot) ./ log.(trans_lres_tot_1)
	# plot
	plot(
		oci_band, [ratio_log_weighted, ratio_log], 
		title="two-way / one-way ratio in log space",
		label=["weighted" "unweighted"],
		titlefontsize=10,
		size=(600, 200)
	)
end

# ╔═╡ 4ee273ba-e908-458d-aa2f-2b0eb7d8ec99


# ╔═╡ 4a1460d5-d31b-4f84-84a4-ef643ae3b952
# ╠═╡ disabled = true
#=╠═╡
begin
	# VCD => xSection
	vmr_o2  = .21
	len     = size(temp)[1]
	vmr_h2o = zeros(Float64, size(temp))
	vcd_dry = zeros(Float64, size(temp))
	vcd_h2o = zeros(Float64, size(temp))
	τ       = zeros(Float64, (len, length(ν_grid)))
	
	@info "Starting"
	@threads for i in 1:len
		p_half_slice = p_half[i, :];
		p_full_slice = p_full[i, :];
		T_slice = temp[i, :];
		q_slice = qw[i, :];
		
		# Apply the function
		vcd_dry_tmp, vcd_h2o_tmp, vmr_h2o_tmp = layer_VCD(
			p_half_slice, q_slice);
	
		if(length(p_full_slice)==length(T_slice))
			# get spectra & optical depth
			xSec_slice_o2 = [o2_sitp(ν_grid, j, k) for (j, k) in zip(p_full_slice, T_slice)];
			xSec_tmp_o2 = hcat(xSec_slice_o2...)
			τ_o2_tmp    = xSec_tmp_o2 * vcd_dry_tmp * vmr_o2;
	
			
			xSec_slice_h2o = [h2o_sitp(ν_grid, j, k) for (j, k) in zip(p_full_slice, T_slice)];
			xSec_tmp_h2o = hcat(xSec_slice_h2o...)
			τ_h2o_tmp    = xSec_tmp_h2o * vcd_h2o_tmp;
			
			# Store the result
			vcd_dry[i, :] = vcd_dry_tmp;
			vcd_h2o[i, :] = vcd_h2o_tmp;
			vmr_h2o[i, :] = vmr_h2o_tmp;
	
			τ[i, :]       = τ_o2_tmp .+ τ_h2o_tmp;
		else
			println("DimensionMismatch!")
		end
	
		if i % 200 == 0
			println("Processed $i / $len samples")
		end
	
	end
	@info "Completed!"
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─8b35e192-f134-11f0-3097-37bfabf9cf36
# ╠═2d8f0c3f-1512-479f-99c3-dfeee855f70f
# ╠═9dcf34aa-4b92-4000-936f-a0a5146c7c08
# ╠═67db2fd7-e281-41c3-89c5-ec396c166ee1
# ╠═9f25f8bc-4085-4cff-bbf9-1f6ed3e7197c
# ╠═52c5cc8a-620b-4132-9586-cf061ef3ee6e
# ╟─28da7053-0ead-49b9-9bec-902f38c0e933
# ╟─47d83de9-43f4-4d75-8665-3d187fc5fddc
# ╟─5fa9abaa-250a-432e-b21f-1d51056ebd5b
# ╠═8c10e29a-f461-4551-bbff-6103d8a3a4e7
# ╟─0e887dd1-a554-4c29-9035-8633b4482c5b
# ╠═72d19300-14cb-41e7-8416-f632586c09f5
# ╠═fb724321-762b-4bac-8e45-2d47d23e56c8
# ╠═7677f2e9-6262-45ba-804c-5e8c19a79c87
# ╟─87b63aea-4852-405e-9a4a-47be24f0c784
# ╟─c1102d6a-b729-444b-9fd7-e90e66d2331b
# ╟─cb7227d2-3e4e-42eb-b208-6573a4aff59e
# ╠═53a3ea6e-cfa9-4d54-ab39-b82aa48eaa2c
# ╟─f04edcf6-20a5-4f89-8a55-79ce2637313e
# ╟─3ad6f119-2bba-45b7-a05b-0133a5522e37
# ╟─94a5b030-4a24-4529-9d14-44bb67c7b863
# ╠═62470600-a00d-4d73-9188-36ec778e436b
# ╟─483d44e8-89b6-43f1-8321-7906a059ac19
# ╠═6e6642c7-20e2-4709-8193-bd93cb74b4eb
# ╠═ca1e324e-c8c3-49b5-a4af-38740afe4f26
# ╟─5c4a20c9-f00d-4197-970d-8bec0935d5db
# ╟─7fedef1f-b0e1-480f-814b-b01dddfd62eb
# ╟─8154d4d2-7ba1-461b-b158-d24435b49c1d
# ╟─7ebf8bc8-8ef3-4c9a-b208-da5bb54f56e3
# ╠═e132db67-7f4e-4fd6-a686-7920fb1c4a74
# ╠═4ee273ba-e908-458d-aa2f-2b0eb7d8ec99
# ╠═4a1460d5-d31b-4f84-84a4-ef643ae3b952
