### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ e1202b2c-da00-4eee-8038-5458ed77bd58
import Pkg; Pkg.activate("..");

# ╔═╡ 2f2f6558-5c22-11f0-3713-9f03610b9eb3
begin
	using vSmartMOM,  vSmartMOM.Absorption
	using JLD2
	using Interpolations
end

# ╔═╡ 085c477d-4b36-4b83-8c95-743c3ab6746a
using NCDatasets, Einsum, Statistics

# ╔═╡ c380f4da-5d23-4c62-9fe1-2d5afbcfa8d3
using ProgressMeter, ProgressLogging

# ╔═╡ f0e56b38-2708-43c8-b70f-997fe35836e2
using Plots

# ╔═╡ 4bcd4623-81a3-46af-b20f-b07dd5da44f8
include("../PACE_SIF.jl")

# ╔═╡ 68bf5501-bd7b-4141-938c-35bf37c3cd9c
md"""
### load interpolation table and scale it
"""

# ╔═╡ d5d1f02d-2ca3-492a-a3c0-c39a8b135e41
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

# ╔═╡ 45c4303a-e814-4414-895c-eeb7e206e47a
md"""
### Change the saved file name!!
"""

# ╔═╡ 81836133-3ff1-4a4a-80ea-61c875e776e9
begin
	# --- reanalysis data used --- #
	dir  = "/home/zhe2/data/MERRA2_reanalysis/";
	file = "MERRA2_400.inst6_3d_ana_Nv.20231230.nc4";
	# --- specify the output NetCDF filename ---
	output_filename = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/transmittance_winter_DefaultModel.nc";
end

# ╔═╡ 8602f29f-bf67-41f1-b3e5-2076181ea5ce
begin
	o2_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Default_grid_Jul10/Default_grid_Jul10_O2.jld2";
	o2_sitp = read_rescale(o2_jld2);
	h2o_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Default_grid_Jul10/Default_grid_Jul10_H2O.jld2";
	h2o_sitp = read_rescale(h2o_jld2);

	metadata = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Default_grid_Jul10/Default_grid_Jul10.log"

	ν_grid, p_grid_hPa, t_grid = o2_sitp.ranges;
end

# ╔═╡ d467d0dc-7478-410e-9bc8-5c6365fde7d0
md"""
### Readin ncDataset, select profiles representative of global atmosphere
"""

# ╔═╡ 6347ff7c-f462-4998-9b25-1e7981354f87
begin
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
end

# ╔═╡ 4fb2a568-0c98-408b-9fe5-f0e126ba78c1
begin
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
	MyStep = 120;

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
	println("$(size(p_full)[1]) profiles selected! ✌️")
end

# ╔═╡ 2ae99ea5-46ce-4f45-8d36-f460cb4ad865
md"""
### Spectral response function
"""

# ╔═╡ 8eeb88f6-abc7-424c-8d14-80b1db89429a
begin
	# read data
	filename = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_RSRs.nc";
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
	
	MyKernel = PACE_SIF.KernelInstrument(
		band=band[ind₂],
		wvlen=λ_msr,
		RSR=RSR[ind₁, ind₂],
		wvlen_out=λ
	);
	println(size(MyKernel.RSR_out));
end

# ╔═╡ 1c293d5e-089f-4eb6-a0e4-27b7ac084cb7
md"""
### VCD => xSection
"""

# ╔═╡ dc3487d6-4ae3-47ea-9787-6b543b17135c
begin
	vmr_o2  = .21
	len     = size(temp)[1];
	vmr_h2o = zeros(Float64, size(temp))
    vcd_dry = zeros(Float64, size(temp))
    vcd_h2o = zeros(Float64, size(temp))
	τ       = zeros(Float64, (len, length(ν_grid)))

	@info "Starting"
	@progress "Processing data..." for i in 1:len
		p_half_slice = p_half[i, :];
		p_full_slice = p_full[i, :];
		T_slice = temp[i, :];
		q_slice = qw[i, :];
		
	    # Apply the function
	    vcd_dry_tmp, vcd_h2o_tmp, vmr_h2o_tmp = PACE_SIF.layer_VCD(
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
		    vcd_dry[i, :] = vcd_dry_tmp
			vcd_h2o[i, :] = vcd_h2o_tmp
			vmr_h2o[i, :] = vmr_h2o_tmp
	
			τ[i, :]       = τ_o2_tmp .+ τ_h2o_tmp;
		else
			println("DimensionMismatch!")
		end
	end
	@info "Completed!"
end

# ╔═╡ 0d52055d-e548-445c-ab15-915aa1eff14f
md"""
### Taking into account the variations of SZA and VZA
Here one-direction transmission (upwelling radiation) and thus VZA is considered.
"""

# ╔═╡ 40b55460-fa8f-4de8-80a1-3f34d2840494
begin
	vza = [0, 15, 35, 45, 75];
	AMF = 1 ./ cosd.(vza);
	# broadcast to# of profiles
	num_rep = floor(Int, len / length(AMF)) + 1;
	# repeat and truncate
	AMF_bc  = repeat(AMF, num_rep)[1:len];
	# multiply to get a slant optical depth => transmission
	trans_hres = exp.(- AMF_bc .* τ );
	# reverse if if_wavenumber is true
	if if_wavenumber
		trans_anly = reverse(trans_hres, dims=2);
		println("reversed!");
	else
		trans_anly = trans_hres;
	end
end

# ╔═╡ 80de2ed8-d0c5-49b0-abdf-17b3d75cf369
begin
	trans = zeros(Float64, (len, length(MyKernel.band)));
	for i in 1:len
		trans[i,:] = MyKernel.RSR_out * trans_anly[i,:];
	end
end

# ╔═╡ d4368739-7301-45e8-9135-7cc642df1cf2
begin
	temp_rd = round.(temp, digits=0);
	ps_rd   = round.(ps, digits=0);
	qw_rd   = round.(qw, digits=4);
	AMF_rd  = round.(AMF_bc, digits=2);
	
	p = plot(MyKernel.band, trans[1,:],
		label="T=$(temp_rd[1,end]) K, p=$(ps_rd[1]) Pa, q=$(qw_rd[1,end]) kg/kg",
		legend=:outerbottom,
		size=(700, 600),    
	);
	for k in 2:377:len
		plot!(p, MyKernel.band, trans[k,:], label="T=$(temp_rd[k,end]) K, p=$(ps_rd[k]) Pa, q=$(qw_rd[k,end]) kg/kg, AMF=$(AMF_rd[k])");
	end

	xlabel!("wv [nm]")
	ylabel!("transmission [-]")
	
	current(p)
	
end

# ╔═╡ a3716d42-e07d-4eab-aa8b-febadb475131
md"""
### Generate a dataset and save
"""

# ╔═╡ e3ab0e36-fc96-498a-8e09-7e3040758f43
begin

	# --- Create the NetCDF Dataset and define its structure ---
	# The `Dataset` function is used with a `do` block. This ensures the file is
	# automatically closed even if errors occur.
	# "c" mode means create (will overwrite if file exists).
	Dataset(output_filename, "c") do ds # `ds` is the Dataset object
	
	    # --- Global Attributes (metadata for the entire file) ---
	    # These provide general information about the dataset.
	    ds.attrib["profile_source"] = "profiles generated from file '$file', step=$MyStep";
		ds.attrib["RT_source"] = "metadata of the interpolation model: '$metadata',  o2 cross section table '$o2_jld2' and h2o cross section table '$h2o_jld2'"
	    ds.attrib["ak"]     = ak;
	    ds.attrib["bk"]     = bk;

		ds.attrib["ν_grid (cm^-1)"] = "Min=$(ν_grid[1]), Max=$(ν_grid[end]), res=$(ν_grid[2]-ν_grid[1])";
		ds.attrib["p_grid (hPa)"]   = p_grid_hPa;
		ds.attrib["T_grid (K)"]     = t_grid;
		ds.attrib["vza (˚)"]        = vza;
	
	    # --- Define Dimensions ---
	    # defDim(dataset_object, dimension_name_string, length)
	    # The length should match the length of your coordinate vectors.
	    defDim(ds, "profile", len)
	    defDim(ds, "layer", n_layers)
	    defDim(ds, "band", length(MyKernel.band))
	
	    # --- Define and Write Variables ---
	    # defVar(dataset_object, variable_name_string, data_type, dimension_names_tuple)
	    # The `dimension_names_tuple` must match the dimensions defined above.
	
	    # --- Coordinate Variables (usually 1D, named after their dimension) ---
	    band_var = defVar(ds, "band", Float64, ("band",))
	    band_var.attrib["long_name"] = "Wavelength"
	    band_var.attrib["units"] = "nm"
	    band_var[:] = MyKernel.band 
	
	    # --- Data Variables ---
	
	    # Variable
	
		trans_var = defVar(ds, "transmission", Float64, ("profile", "band", ))
		trans_var.attrib["long_name"] = "atmospheric transmission (one-way)"
		trans_var.attrib["unit"] = "unitless"
		trans_var.attrib["description"] = "only gas absorption - O2 and H2O"
		trans_var[:, :] = trans
		
	    temp_var = defVar(ds, "temperature", Float64, ("profile", "layer", ))
	    temp_var.attrib["long_name"] = "air Temperature"
	    temp_var.attrib["units"] = "K"
	    temp_var[:, :] = temp # Write the actual data to the variable
	
		pres_var = defVar(ds, "pressure", Float64, ("profile", ))
		pres_var.attrib["long_name"] = "surface pressure"
		pres_var.attrib["units"] = "hPa"
		pres_var[:] = ps / 100.
	
		q_var = defVar(ds, "q", Float64, ("profile", "layer", ))
		q_var.attrib["long_name"] = "specific humidity"
		q_var.attrib["unit"] = "kg/kg"
		q_var[:, :] = qw
	
		vcd_dry_var = defVar(ds, "vcd_dry", Float64, ("profile", "layer", ))
		vcd_dry_var.attrib["long_name"] = "vertical column density of dry air per layer"
		vcd_dry_var.attrib["unit"] = "molec/cm^2"
		vcd_dry_var[:, :] = vcd_dry
	
		vcd_h2o_var = defVar(ds, "vcd_h2o", Float64, ("profile", "layer", ))
		vcd_h2o_var.attrib["long_name"] = "vertical column density of water vapor per layer"
		vcd_h2o_var.attrib["unit"] = "molec/cm^2"
		vcd_h2o_var[:, :] = vcd_h2o
	
		vmr_h2o_var = defVar(ds, "vmr_h2o_var", Float64, ("profile", "layer", ))
		vmr_h2o_var.attrib["long_name"] = "vertical mixing ratio of H2O of water vapor per layer"
		vmr_h2o_var.attrib["unit"] = "molec/molec"
		vmr_h2o_var[:, :] = vmr_h2o
	
		AMF_var = defVar(ds, "AMF", Float64, ("profile", ))
		AMF_var.attrib["long_name"] = "air mass factor (one-way)"
		AMF_var.attrib["unit"] = "unitless"
		AMF_var[:] = AMF_bc
		
	    println("\nNetCDF Dataset created successfully and saved to '$output_filename'.")
	    println("File details (from NCDatasets.jl):")
	    show(ds) 
	end
end

# ╔═╡ Cell order:
# ╠═e1202b2c-da00-4eee-8038-5458ed77bd58
# ╟─68bf5501-bd7b-4141-938c-35bf37c3cd9c
# ╠═2f2f6558-5c22-11f0-3713-9f03610b9eb3
# ╟─d5d1f02d-2ca3-492a-a3c0-c39a8b135e41
# ╟─45c4303a-e814-4414-895c-eeb7e206e47a
# ╠═81836133-3ff1-4a4a-80ea-61c875e776e9
# ╠═8602f29f-bf67-41f1-b3e5-2076181ea5ce
# ╟─d467d0dc-7478-410e-9bc8-5c6365fde7d0
# ╠═085c477d-4b36-4b83-8c95-743c3ab6746a
# ╠═6347ff7c-f462-4998-9b25-1e7981354f87
# ╠═4fb2a568-0c98-408b-9fe5-f0e126ba78c1
# ╟─2ae99ea5-46ce-4f45-8d36-f460cb4ad865
# ╠═4bcd4623-81a3-46af-b20f-b07dd5da44f8
# ╠═8eeb88f6-abc7-424c-8d14-80b1db89429a
# ╟─1c293d5e-089f-4eb6-a0e4-27b7ac084cb7
# ╠═c380f4da-5d23-4c62-9fe1-2d5afbcfa8d3
# ╠═dc3487d6-4ae3-47ea-9787-6b543b17135c
# ╟─0d52055d-e548-445c-ab15-915aa1eff14f
# ╠═40b55460-fa8f-4de8-80a1-3f34d2840494
# ╠═80de2ed8-d0c5-49b0-abdf-17b3d75cf369
# ╠═f0e56b38-2708-43c8-b70f-997fe35836e2
# ╠═d4368739-7301-45e8-9135-7cc642df1cf2
# ╟─a3716d42-e07d-4eab-aa8b-febadb475131
# ╠═e3ab0e36-fc96-498a-8e09-7e3040758f43
