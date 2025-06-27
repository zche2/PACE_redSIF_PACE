### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 1252addf-8f35-4b86-8c0b-6c842aeee51b
import Pkg; Pkg.activate("..")

# ╔═╡ fd46ebf6-5122-11f0-20b4-0f8c9af052ee
using Markdown, InteractiveUtils, Plots, NCDatasets

# ╔═╡ 5d6daeae-98e9-4f9c-ac4a-bd6e4e32eb75
using vSmartMOM,  vSmartMOM.Absorption

# ╔═╡ e98384a5-11e8-4d03-91ac-7b386ff60cf1
using Einsum

# ╔═╡ 2cf3bad7-d64d-4ca7-8aef-5cd3249dba48
using PlutoUI

# ╔═╡ 9a1a87fe-6632-46ef-9de0-15d78b5fa82f
include("/home/zhe2/FraLab/PACE_redSIF_PACE/tools.jl")

# ╔═╡ 15209f18-eb42-4cbf-bd2e-e4bd8755decd
md"""
> ### O$_2$ transmission spectra
To activate packages in PACE project
"""

# ╔═╡ c49503cb-de14-478f-b5a2-6e8d823d6129
begin
	res   = 0.05;
	ν_min = 11627 
	ν_max = 16129
	ν = ν_min:res:ν_max;
	λ = 1 ./ ν .* 1E7;  # cm^-1 to nm
	if_wavenumber = true;
end

# ╔═╡ 4fae23a8-a0a8-4dd9-a610-a43211a4f6f9
md"""
> ##### Read atmospheric profile and get xsection & transmission for arbitrary number of absorbers 
"""

# ╔═╡ 68c5fbe1-45d2-41e2-a6ab-ff95c1a06804
function layer_VCD(
	p,    # pressure (Pa)
	q,    # specific humidity
	n_layers :: Int;       # number of layers
	g₀ = 9.8196,    # Avogradro's number
	Na = 6.0221415e23
)
	
    # Dry and wet mass
    dryMass = 28.9647e-3  / Na  # in kg/molec, weighted average for N2 and O2
    wetMass = 18.01528e-3 / Na  # just H2O
    ratio = dryMass / wetMass 

	vmr_h2o = zeros(Float64, n_layers)
    vcd_dry = zeros(Float64, n_layers)
    vcd_h2o = zeros(Float64, n_layers)
	
    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p[i + 1] - p[i]
        vmr_h2o[i] = q[i] * ratio
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dryMass + vmr_h2o[i] * wetMass
        vcd_dry[i] = vmr_dry * Δp / (M * g₀ * 100.0^2)   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * Δp / (M * g₀ * 100^2)
    end

	return vcd_dry, vcd_h2o, vmr_h2o
end

# ╔═╡ 90a38cea-dcd5-4892-a2d3-bfc18f208ac5
md"""
read hitran table and make hitran model
"""

# ╔═╡ 722d3879-7e43-4e34-8d65-04b888047afb
begin
	FT = Float64;
	
	# TimeIndex = 1;
    @einsum p_half[i,j,r,k] := (ak[r] + bk[r] * psurf[i,j,k])
    p_full = (p_half[:,:, 2:end ,:] + p_half[:,:,1:end - 1,:]) / 2
	
    # also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, size(T))
    vcd_dry = zeros(FT, size(T))
    vcd_h2o = zeros(FT, size(T))

	n_layers = Int(size(T, 3))  # 72
	
	lon_dim = length(lon);
	lat_dim = length(lat);
	time_dim = length(time)

	for (lon_idx, lat_idx, time_idx) in Iterators.product(
		1:lon_dim, 1:lat_dim, 1:time_dim)
	    
		p_slice = p_half[lon_idx, lat_idx, :, time_idx]
		q_slice = q[lon_idx, lat_idx, :, time_idx]
	    
	    # Apply the function
	    vcd_dry_tmp, vcd_h2o_tmp, vmr_h2o_tmp = layer_VCD(p_slice, q_slice, n_layers)
	    
	    # Store the result
	    vcd_dry[lon_idx, lat_idx, :, time_idx] = vcd_dry_tmp
		vcd_h2o[lon_idx, lat_idx, :, time_idx] = vcd_h2o_tmp
		vmr_h2o[lon_idx, lat_idx, :, time_idx] = vmr_h2o_tmp
	end

end

# ╔═╡ 027c01f8-54ee-400f-a54e-612bb8224ab8
# ╠═╡ disabled = true
#=╠═╡
begin
	# save！
	
	# --- Specify the output NetCDF filename ---
	output_filename = "../sample_data/atm_profile_vmrvcd.nc"
	
	# --- Create the NetCDF Dataset and define variables/dimensions ---
	# The `Dataset` function is used with a `do` block, which ensures the file is closed automatically.
	Dataset(output_filename, "c") do ds # "c" for create mode (overwrites if exists)
	
	    # Global attributes (metadata for the entire file)
		ds.attrib["description"] = "vertical mixing ratio (vmr), vertical column density (vcd) calculated from MERRA3 atmospheric profile, see ESE156 Lect04"
	
	    # Define dimensions
	    # defDim(dataset_object, dimension_name_string, length)
	    defDim(ds, "longitude", length(lon))
	    defDim(ds, "latitude", length(lat))
		defDim(ds, "z_layer", n_layers)
	    defDim(ds, "time", length(time))
	
	    # Define Variables (including coordinate variables)
	    # defVar(dataset_object, variable_name_string, data_type, dimension_names_tuple)
	
	    # Coordinate Variable
		lon_var = defVar(ds, "longitude", Float64, ("longitude",))
		lon_var[:] = lon

		lat_var = defVar(ds, "latitude", Float64, ("latitude",))
		lat_var[:] = lat

		time_var = defVar(ds, "time", Float64, ("time",))
		time_var[:] = collect(0:6:18);

		# var1: air temperature
		temp_var = defVar(ds, "temperature",
					Float64, ("longitude", "latitude", "z_layer", "time"))
		temp_var.attrib["long_name"] = "air temperature";
		temp_var.attrib["units"] = "K";
		temp_var[:,:,:,:] = T;
		
		# var2: specific humidity
		q_var = defVar(ds, "q",
					Float64, ("longitude", "latitude", "z_layer", "time"))
		q_var.attrib["long_name"] = "specific humidity";
		q_var.attrib["units"] = "kg/kg";
		q_var[:,:,:,:] = q;
		
		# var3: pressure
		pres_var = defVar(ds, "p",
					Float64, ("longitude", "latitude", "z_layer", "time"))
		pres_var.attrib["long_name"] = "mean pressure of the layer (not half level / boundary pressures!)";
		pres_var.attrib["units"] = "Pa";
		pres_var[:,:,:,:] = p_full;

		# var4: VCD
		dry_var = defVar(ds, "vcd_dry",
					Float64, ("longitude", "latitude", "z_layer", "time"))
		dry_var.attrib["long_name"] = "vertical column density of dry air"
		dry_var.attrib["units"] = "molec/cm^2";
		dry_var[:,:,:,:] = vcd_dry;

		h2o_var = defVar(ds, "vcd_h2o",
					Float64, ("longitude", "latitude", "z_layer", "time"))
		h2o_var.attrib["long_name"] = "vertical column density of water vapor";
		h2o_var.attrib["units"] = "molec/cm^2";
		h2o_var[:,:,:,:] = vcd_h2o;

		# var5: vmr of H2O
		vmr_var = defVar(ds, "vmr_h2o",
					Float64, ("longitude", "latitude", "z_layer", "time"))
		vmr_var.attrib["long_name"] = "vertical mixing ratio of H2O";
		vmr_var.attrib["units"] = "molec/molec";
		vmr_var[:,:,:,:] = vmr_h2o;
		
		
	    println("\nNetCDF Dataset created successfully and saved to '$output_filename'.")
	    println("File details (from NCDatasets.jl):")
	    show(ds) # Show the dataset summary in the REPL
	end

end
  ╠═╡ =#

# ╔═╡ dbcf9eff-f9cd-4425-abc6-692442502d97
q

# ╔═╡ f1cc3038-2473-4372-9ddf-fc6616f908af
# ╠═╡ disabled = true
#=╠═╡
begin
	a = sum(vcd_h2o4, dims=3);
	@show minimum(a)
	@show maximum(a)
end
  ╠═╡ =#

# ╔═╡ 3fe4b0e4-73ea-4ff9-b845-a6f805d0dda0
# ╠═╡ disabled = true
#=╠═╡
begin
	# Create matrix of cross sections for each atmospheric layer (takes quite some time!!):
	# This 
	cs_matrix_o2  = zeros((lon_dim, lat_dim, length(ν), n_layers, time_dim))
	cs_matrix_h2o = zeros((lon_dim, lat_dim, length(ν), n_layers, time_dim))

	for (lon_idx, lat_idx, time_idx) in Iterators.product(
		1:lon_dim, 1:lat_dim, 1:time_dim)
		
		p_slice = p_full[lon_idx, lat_idx, :, time_idx]
		T_slice = T[lon_idx, lat_idx, :, time_idx]
		
		# Loop over each layer 
		for i=1:n_layers
		    p_ = p_slice[i] / 100 # in hPa
		    T_ = T_slice[i]
		    cs_matrix_o2[lon_idx, lat_idx,:,i,time_idx] = 		absorption_cross_section(o2_voigt, ν, p_, T_);
		    cs_matrix_h2o[lon_idx, lat_idx,:,i,time_idx] = absorption_cross_section(h2o_voigt, ν, p_, T_);
		end
	end
end
  ╠═╡ =#

# ╔═╡ 3445b1e9-58f6-448b-9400-0c0c290a6af6
md"""
================ Place Holder ================
\
\
\
\
\
\
\
==============================================
"""

# ╔═╡ 7e7f5dd8-df67-44d6-a01f-9039e4d38e3a
md"""
> #### generate gridded T, p, and hitran model of O$_2$ / H$_2$O
"""

# ╔═╡ 0756e914-0cd1-449b-9681-327fc8b582e7
begin
	dT = 10.      # K
	dp = 20.      # hPa
	T_min = 100   # K
	T_max = 400   # K
	p_min = 100   # hPa
	p_max = 1200  # hPa
	T_grid = T_min:dT:T_max;
	p_grid = p_min:dp:p_max;
end

# ╔═╡ 579423f0-cde7-4acc-bd4f-02d78350a156
#=╠═╡
println("interpolation model generated: $(size(itp_model.itp))")
  ╠═╡ =#

# ╔═╡ e31af66b-cb28-46c7-8061-6370d15c7b09
md"""
> #### Scale to transmission and convolve to PACE resolution
"""

# ╔═╡ 957d5a2e-878f-4d6c-9d2d-24e501fc29e6
# ╠═╡ disabled = true
#=╠═╡
begin
	# get back to atm. transmission => a.u. as it does not matter for SVD
	τ_scale = 1.0e22;
		#  4.29e24;  # sum(vcd_dry' * .21); VCD from measurements of O2
	atm_trans = exp.( - itp_model.itp * τ_scale );
end
  ╠═╡ =#

# ╔═╡ 30fd851c-bb61-49e2-8f0b-ecf1ea8b0d92
# ╠═╡ disabled = true
#=╠═╡
begin
	# preallocate
	spec_conv_loop = 
		zeros(eltype(atm_trans_ref),
			(length(MyKernel.band), length(p_grid), length(T_grid))) * NaN;
	for j=1:length(p_grid)
		for k=1:length(T_grid)
			spec_conv_loop[:, j, k] = MyKernel.RSR_out * atm_trans_ref[:, j, k]
		end
	end
	println("successfully convolved!, size: $(size(spec_conv_loop))")
end
  ╠═╡ =#

# ╔═╡ a28ba46c-f9b8-4b71-83ab-11f9758adc55
@bind i Slider(1:56)

# ╔═╡ 88b21dad-d2b9-4e38-a4c8-ccc3d163236d
@bind j Slider(1:31)

# ╔═╡ bf1bff02-32f0-41ff-91a5-c3943c6ff285
# ╠═╡ disabled = true
#=╠═╡
begin
	plot(λ_ref, atm_trans_ref[:, i, j],
		label="original res=$res [cm-1]",
		alpha=.5)
	plot!(MyKernel.band, spec_conv_loop[:, i, j],
		label="conv",
		linewidth=2.)
	xlabel!("[nm]")
	ylabel!("transmission")
	title!("P=$(p_grid[i]) hPa, T=$(T_grid[j]) K")
end
  ╠═╡ =#

# ╔═╡ 55451be9-ac3d-4758-9958-7b54852d028e
md"""
> #### Create dataset and save
"""

# ╔═╡ 3bdf2f12-0408-46ca-bf59-aa97db3d7399
# ╠═╡ disabled = true
#=╠═╡
begin
	# --- Specify the output NetCDF filename ---
	# output_filename = "../sample_data/H2O_transmission.nc"
	
	# --- 4. Create the NetCDF Dataset and define variables/dimensions ---
	# The `Dataset` function is used with a `do` block, which ensures the file is closed automatically.
	Dataset(output_filename, "c") do ds # "c" for create mode (overwrites if exists)
	
	    # Global attributes (metadata for the entire file)
		ds.attrib["attr"] = "prescribed vmr, τ_scale = 1.0e22, single layer, see ESE156 Lect04"
	
	    # Define dimensions
	    # defDim(dataset_object, dimension_name_string, length)
	    defDim(ds, "temperature", length(T_grid))
	    defDim(ds, "pressure", length(p_grid))
	    defDim(ds, "band", length(MyKernel.band))
	
	    # Define Variables (including coordinate variables)
	    # defVar(dataset_object, variable_name_string, data_type, dimension_names_tuple)
	
	    # Coordinate Variable
	    temp_var = defVar(ds, "temperature", Float64, ("temperature",))
	    temp_var.attrib["units"] = "K" 
	    temp_var[:] = T_grid # Write the coordinates data
	
	    pres_var = defVar(ds, "pressure", Float64, ("pressure",))
	    pres_var.attrib["units"] = "hPa" 
	    pres_var[:] = p_grid # Write the coordinates data
	
		band_var = defVar(ds, "band", Float64, ("band",))
	    band_var.attrib["units"] = "nm" 
	    band_var[:] = MyKernel.band # Write the coordinates data
	
	    # Data Variable: Air Temperature
	    O2_trans_var = defVar(ds, "H2O_trans", Float64, ("band", "pressure", "temperature"))
	    O2_trans_var.attrib["long_name"] = "convolved transmission spectrum of O2"
	    O2_trans_var.attrib["units"] = "unitless"
	    O2_trans_var.attrib["_FillValue"] = NaN # Optional: specify a fill value
	    O2_trans_var[:, :, :] = spec_conv_loop # Write the data
	
	    println("\nNetCDF Dataset created successfully and saved to '$output_filename'.")
	    println("File details (from NCDatasets.jl):")
	    show(ds) # Show the dataset summary in the REPL
	end

end
  ╠═╡ =#

# ╔═╡ 5d186afe-cc00-42a9-aa88-bbaed86a5908
#=╠═╡
begin
	MF = "../sample_data/atm_profile.nc"
	
	ds = Dataset(MF)
	
	# See how easy it is to actually extract data? Note the [:] in the end reads in ALL the data in one step
	lat   = ds["YDim"][:]
	lon   = ds["XDim"][:]
	# Temperature profile
	T     = ds["T"][:]
	# specific humidity profile
	q     = ds["QV"][:]
	# mean pressure profile:
	# p     = ds["Height"][:]
	# Surafce pressure
	psurf = ds["PS"][:]
	# Time in UTC
	time  = ds["TIME"][:]
	
	# AK and BK global attributes (important to calculate pressure half-levels)
	ak = ds.attrib["HDF_GLOBAL.ak"][:]
	bk = ds.attrib["HDF_GLOBAL.bk"][:]
	
	close(ds)
end
  ╠═╡ =#

# ╔═╡ 11a03d19-0302-4dd6-8db5-aaaadf18a360
# ╠═╡ disabled = true
#=╠═╡
begin
	o2_par = 
		Absorption.read_hitran(artifact("H2O"), mol=1, iso=1, ν_min=ν_min, ν_max=ν_max);
		# Absorption.read_hitran(artifact("O2"), mol=7, iso=1, ν_min=ν_min, ν_max=ν_max);
	
	itp_model = 
		make_interpolation_model(
			o2_par,
			Voigt(),
			ν,
			p_grid,
			T_grid,
			wavelength_flag=false,
			wing_cutoff=10,
			architecture=CPU()
		)
end
  ╠═╡ =#

# ╔═╡ d1148803-1909-4706-a885-afe6916318e9
#=╠═╡
begin
	# read data
	filename = "/home/zhe2/FraLab/PACE_redSIF_PACE/sample_data/PACE_OCI_RSRs.nc";
	ds = Dataset(filename, "r");

	wavlen = ds["wavelength"][:];
	RSR = ds["RSR"];
	band = ds["bands"];

	λ_ref = if_wavenumber ? reverse(λ) : λ;
	atm_trans_ref = if_wavenumber ? reverse(atm_trans, dims=1) : atm_trans ;
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
  ╠═╡ =#

# ╔═╡ c03fbd2e-6798-4c6d-adbd-250861774c37
#=╠═╡
begin
	# (we have to know the HITRAN molecule numbers, given in http://hitran.org/docs/molec-meta/)
	# Read in HITRAN tables
	o2_par  = Absorption.read_hitran(artifact("O2"), mol=7, iso=1, ν_min=ν_min, ν_max=ν_max);
	h2o_par = Absorption.read_hitran(artifact("H2O"), mol=1, iso=1, ν_min=ν_min, ν_max=ν_max);

	o2_voigt   = make_hitran_model(o2_par, Voigt(), wing_cutoff=10, architecture=CPU());
	h2o_voigt   = make_hitran_model(h2o_par, Voigt(), wing_cutoff=10, architecture=CPU());
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─15209f18-eb42-4cbf-bd2e-e4bd8755decd
# ╠═1252addf-8f35-4b86-8c0b-6c842aeee51b
# ╠═fd46ebf6-5122-11f0-20b4-0f8c9af052ee
# ╠═5d6daeae-98e9-4f9c-ac4a-bd6e4e32eb75
# ╠═9a1a87fe-6632-46ef-9de0-15d78b5fa82f
# ╠═e98384a5-11e8-4d03-91ac-7b386ff60cf1
# ╠═c49503cb-de14-478f-b5a2-6e8d823d6129
# ╟─4fae23a8-a0a8-4dd9-a610-a43211a4f6f9
# ╠═5d186afe-cc00-42a9-aa88-bbaed86a5908
# ╟─68c5fbe1-45d2-41e2-a6ab-ff95c1a06804
# ╟─90a38cea-dcd5-4892-a2d3-bfc18f208ac5
# ╠═c03fbd2e-6798-4c6d-adbd-250861774c37
# ╠═722d3879-7e43-4e34-8d65-04b888047afb
# ╠═027c01f8-54ee-400f-a54e-612bb8224ab8
# ╠═dbcf9eff-f9cd-4425-abc6-692442502d97
# ╠═f1cc3038-2473-4372-9ddf-fc6616f908af
# ╠═3fe4b0e4-73ea-4ff9-b845-a6f805d0dda0
# ╟─3445b1e9-58f6-448b-9400-0c0c290a6af6
# ╟─7e7f5dd8-df67-44d6-a01f-9039e4d38e3a
# ╠═0756e914-0cd1-449b-9681-327fc8b582e7
# ╠═11a03d19-0302-4dd6-8db5-aaaadf18a360
# ╟─579423f0-cde7-4acc-bd4f-02d78350a156
# ╟─e31af66b-cb28-46c7-8061-6370d15c7b09
# ╠═957d5a2e-878f-4d6c-9d2d-24e501fc29e6
# ╠═d1148803-1909-4706-a885-afe6916318e9
# ╠═30fd851c-bb61-49e2-8f0b-ecf1ea8b0d92
# ╠═2cf3bad7-d64d-4ca7-8aef-5cd3249dba48
# ╠═a28ba46c-f9b8-4b71-83ab-11f9758adc55
# ╠═88b21dad-d2b9-4e38-a4c8-ccc3d163236d
# ╠═bf1bff02-32f0-41ff-91a5-c3943c6ff285
# ╟─55451be9-ac3d-4758-9958-7b54852d028e
# ╠═3bdf2f12-0408-46ca-bf59-aa97db3d7399
