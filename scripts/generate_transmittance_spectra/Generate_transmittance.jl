import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE");
include("/home/zhe2/FraLab/PACE_redSIF_PACE/PACE_SIF.jl")

using vSmartMOM, vSmartMOM.Absorption
using JLD2
using Interpolations
using NCDatasets, Einsum, Statistics
using Parameters, ProgressMeter, ProgressLogging

function read_rescale(itp_filename::String)
	model = load_interpolation_model(itp_filename);
	spec_grid = if hasproperty(model, :ν_grid)
		model.ν_grid
	elseif hasproperty(model, :λ_grid)
		model.λ_grid
	else
		error("Interpolation model must contain either :ν_grid or :λ_grid")
	end
	p_grid = model.p_grid;
	t_grid = model.t_grid;
	itp    = model.itp;
	sitp = scale(itp, spec_grid, p_grid, t_grid);
	println("scaled! $itp_filename")
	return sitp
end

if abspath(PROGRAM_FILE) == @__FILE__

	# Files path

	# --- reanalysis data used --- #
	dir  = "/home/zhe2/data/MERRA2_reanalysis/";
	file = "MERRA2_400.inst6_3d_ana_Nv.20240705.nc4";
	# --- specify the output NetCDF filename ---
	output_filename = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/transmittance_summer_FineWvResModel_FullRange_Aug01.nc";

	xsec_subdir = "wavelength-run";
	xsec_dir = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/$xsec_subdir";
	o2_jld2 = joinpath(xsec_dir, "$(xsec_subdir)_O2.jld2");
	o2_sitp = read_rescale(o2_jld2);
	h2o_jld2 = joinpath(xsec_dir, "$(xsec_subdir)_H2O.jld2");
	h2o_sitp = read_rescale(h2o_jld2);

	metadata = joinpath(xsec_dir, "$(xsec_subdir).log")

	spec_grid_raw, p_grid_hPa, t_grid = o2_sitp.ranges;
	spec_grid_raw = collect(spec_grid_raw);

	# Detect spectral axis unit and reorder to increasing wavelength for high-res forward runs.
	# Note: For wavelength-run LUTs, the first axis is wavelength even if field name is ν_grid.
	is_wavenumber_axis = maximum(spec_grid_raw) > 3000.0;
	if is_wavenumber_axis
		λ_from_raw = PACE_SIF.ν_to_λ.(spec_grid_raw);
		axis_unit = "wavenumber_cm^-1";
	else
		λ_from_raw = spec_grid_raw;
		axis_unit = "wavelength_nm";
	end

	if λ_from_raw[1] <= λ_from_raw[end]
		spec_eval_grid = spec_grid_raw;
		λ_hres = λ_from_raw;
	else
		spec_eval_grid = reverse(spec_grid_raw);
		λ_hres = reverse(λ_from_raw);
	end
	println("high-res axis unit: $axis_unit, n=$(length(spec_eval_grid))")

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
	MyStep = 95;

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


	# Spectral response function
	# read data
	filename = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI_RSRs.nc";
	pace = Dataset(filename, "r");

	wavlen = pace["wavelength"][:];
	RSR = pace["RSR"];
	band = pace["bands"];

	ind₁   = findall( λ_hres[1] .< wavlen .< λ_hres[end]);
	ind₂   = findall( λ_hres[1] .< band   .< λ_hres[end]);
	λ_msr  = wavlen[ind₁];
	MyRSR  = RSR[ind₁, ind₂];

	MyKernel = PACE_SIF.KernelInstrument(
		band[ind₂],
		λ_msr,
		MyRSR,
		λ_hres,
		spec_eval_grid,
	);
	close(pace)
	println(size(MyKernel.RSR_out));


	# VCD => xSection

	vmr_o2  = .21
	len     = size(temp)[1];
	vmr_h2o = zeros(Float64, size(temp))
	vcd_dry = zeros(Float64, size(temp))
	vcd_h2o = zeros(Float64, size(temp))
	τ       = zeros(Float64, (len, length(λ_hres)))

	@info "Starting"
	@showprogress 5 "Processing data..." for i in 1:len
		p_half_slice = p_half[i, :];
		p_full_slice = p_full[i, :];
		T_slice = temp[i, :];
		q_slice = qw[i, :];
		
		# Apply the function
		vcd_dry_tmp, vcd_h2o_tmp, vmr_h2o_tmp = PACE_SIF.layer_VCD(
			p_half_slice, q_slice);

		if(length(p_full_slice)==length(T_slice))
			# get spectra & optical depth
			xSec_slice_o2 = [o2_sitp(spec_eval_grid, j, k) for (j, k) in zip(p_full_slice, T_slice)];
			xSec_tmp_o2 = hcat(xSec_slice_o2...)
			τ_o2_tmp    = xSec_tmp_o2 * vcd_dry_tmp * vmr_o2;

			
			xSec_slice_h2o = [h2o_sitp(spec_eval_grid, j, k) for (j, k) in zip(p_full_slice, T_slice)];
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


	# --- Taking into account the variations of SZA and VZA ---
	vza = [0, 15, 35, 45, 60, 75];
	AMF = 1 ./ cosd.(vza);
	# broadcast to# of profiles
	num_rep = floor(Int, len / length(AMF)) + 1;
	# repeat and truncate
	AMF_bc  = repeat(AMF, num_rep)[1:len];
	# multiply to get a slant optical depth => transmission
	trans_hres = exp.(- AMF_bc .* τ );


	trans = zeros(Float64, (len, length(MyKernel.band)));
	for i in 1:len
		trans[i,:] = MyKernel.RSR_out * trans_hres[i,:];
	end


	# --- Create the NetCDF Dataset and define its structure ---
	# The `Dataset` function is used with a `do` block. This ensures the file is
	# automatically closed even if errors occur.
	# "c" mode means create (will overwrite if file exists).
	Dataset(output_filename, "c") do ds # `ds` is the Dataset object

		# --- Global Attributes (metadata for the entire file) ---
		# These provide general information about the dataset.
		ds.attrib["profile_source"] = "profiles generated from file '$file', step=$MyStep";
		ds.attrib["RT_source"]      = "metadata of the interpolation model: '$metadata',  o2 cross section table '$o2_jld2' and h2o cross section table '$h2o_jld2'"
		ds.attrib["ak"]     = ak;
		ds.attrib["bk"]     = bk;

		ds.attrib["highres_axis_unit"] = axis_unit;
		ds.attrib["highres_axis"]      = "Min=$(spec_eval_grid[1]), Max=$(spec_eval_grid[end]), n=$(length(spec_eval_grid))";
		ds.attrib["λ_hres (nm)"]       = "Min=$(λ_hres[1]), Max=$(λ_hres[end]), n=$(length(λ_hres))";
		ds.attrib["p_grid (hPa)"]      = p_grid_hPa;
		ds.attrib["T_grid (K)"]        = t_grid;
		ds.attrib["vza (˚)"]           = vza;

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

		trans_var = defVar(ds, "transmittance", Float64, ("profile", "band", ))
		trans_var.attrib["long_name"] = "atmospheric transmittance (one-way)"
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
	end

end
