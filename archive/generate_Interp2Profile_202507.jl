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

# ╔═╡ 87be1d01-00fc-41b6-9ce8-5d0a4518c807
using ProgressMeter

# ╔═╡ 9dd18d77-3d0f-4665-9dde-fe76d8923fa4
using ProgressLogging

# ╔═╡ 1e4371af-86fd-41b0-b0e2-432a7c39fc6e
using Interpolations

# ╔═╡ a2a8261e-0161-42ad-87f4-aa0c461cabf5
using Markdown, InteractiveUtils, Plots, NCDatasets

# ╔═╡ c7c5d86b-001a-4337-9e56-6821a6f7298b
using PlutoUI, Statistics, Einsum

# ╔═╡ 80a94e86-fd79-40c8-b058-ce21de36918f
using vSmartMOM,  vSmartMOM.Absorption

# ╔═╡ a35441c6-7551-419c-8934-013dd35a8277
using JLD2

# ╔═╡ bc2971ec-5768-11f0-3648-c1692dc3b9da
# import Pkg; Pkg.activate("..");

# ╔═╡ 7503ae08-6dcc-45c9-95a9-ae9ace863e66
md"""
> ### Read analysis profile and randomly select profiles that are represntative of global atmosphere
Season / Day & Night/ Land & Ocean?
"""

# ╔═╡ 9458cdfd-3a36-4b9d-aee4-f5794bf4e70d
begin
	dir  = "/home/zhe2/data/MERRA2_reanalysis/";
	file = "MERRA2_400.inst6_3d_ana_Nv.20240705.nc4";
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

# ╔═╡ fa75b745-dbc6-44c8-aa0f-4a8c5c21d9a6
begin
	# determine the range of interpolation?
	println("======== temperature (K) ========")
	println("minimum T: $(minimum(T))");
	println("5th T: $(quantile(T[:], 0.05))");
	println("95th T: $(quantile(T[:], 0.95))");
	println("maximum T: $(maximum(T))");

	@einsum p_half[i,j,r,k] := (ak[r] + bk[r] * psurf[i,j,k]);
	p_full = (p_half[:,:, 2:end ,:] + p_half[:,:,1:end - 1,:]) / 2
	println("======== pressure (Pa) ========")
	println("minimum p: $(minimum(p_full))");
	println("5th p: $(quantile(p_full[:], 0.05))");
	println("95th p: $(quantile(p_full[:], 0.95))");
	println("maximum p: $(maximum(p_full))");
end

# ╔═╡ ed7e418f-500d-4289-81a2-5468c5a963db
@bind TimeIndex Slider(1:4, default=1)

# ╔═╡ 5d4bb7fa-3c60-4a06-8757-d72848a01464
begin
	heatmap(lon, lat, T[:,:,end,TimeIndex]' , size = (800, 500))
end

# ╔═╡ cc4c4ddf-805e-4c78-ba5d-baf7743be139
md"""
> #### Load the `itp_model` and interpolate
now able to generate a matrix for each (lon, lat) point
"""

# ╔═╡ be766ea1-94b1-4ec1-9f6d-30a2da1b8faa
o2_model = load_interpolation_model("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection_O2_try1");

# ╔═╡ a844cb2f-f892-46c6-8825-f1fd35767e46
h2o_model = load_interpolation_model("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection_H2O_try1");

# ╔═╡ 1a73db90-e275-4c36-a793-697a5e8ed9ce
# ╠═╡ disabled = true
#=╠═╡
begin
	# single layer
	
	# given an arbitrary T and p
	T_local = T[300, 120, end, 3];
	p_local = psurf[300, 120, 3];
	# Pa -> hPa
	p_local /= 100;
	@show T_local, p_local
	# evaluate xSection at this T and p
	xSec = sitp(ν_grid, p_local, T_local);
end
  ╠═╡ =#

# ╔═╡ 8fe5df77-2ae2-4ba6-a239-3c927d55920e
# ╠═╡ disabled = true
#=╠═╡
begin
	# interpolate the whole column at the same time
	T_column = T[245, 168, :, 4];
	p_column = p_full[245, 168, :,4] / 100;  # hPa
	# list
	xSec_column = [sitp(ν_grid, i, j) for (i,j) in zip(p_column, T_column)];
	@show size(xSec_column)
	# unpack
	xSec_mat = hcat(xSec_column...)
	@show size(xSec_mat')

	# hcat: each vector becomes a column
	# vcat: each vector becomes a row
end
  ╠═╡ =#

# ╔═╡ 7c3d7f47-83bc-48bc-ad21-d49209094468
md"""
> #### vertical profile of VCD => optical depth
"""

# ╔═╡ 85b0a780-ee31-4491-a914-7d64aad289ea
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

# ╔═╡ 23892322-e14a-4155-adbf-6e7a00cf5a93
begin
	dgrid = 20;
	# select profile
	ind1 = 1:dgrid:nLon;
	ind2 = 1:dgrid:nLat;
	T_subset = T[ind1, 1:dgrid:nLat, :, :];
	q_subset = q[ind1, 1:dgrid:nLat, :, :];
	p_subset = p_half[ind1, 1:dgrid:nLat, :, :];
	p_ave_subset = p_full[ind1, 1:dgrid:nLat, :, :];
	n_layers = 72;
end

# ╔═╡ 18a3a4d8-d6e0-441c-9957-ee36c825b012
begin
	# swap
	swap_dim = (1, 2, 4, 3) 
	T_swap = permutedims(T_subset, swap_dim);
	q_swap = permutedims(q_subset, swap_dim);
	p_swap = permutedims(p_subset, swap_dim);
	p_ave_swap = permutedims(p_ave_subset, swap_dim);
	# reshape
	T_rshp = reshape(T_swap, :, n_layers);
	q_rshp = reshape(q_swap, :, n_layers);
	p_rshp = reshape(p_swap, :, n_layers+1);
	p_ave_rshp = reshape(p_ave_swap, :, n_layers);
	# finished!
	println("reshaped!")
end

# ╔═╡ aedd5952-3897-43c0-9b77-e4e827f5bb9b
# ╠═╡ disabled = true
#=╠═╡
begin
	len, _  = size(T_rshp);

	vmr_o2  = .21
	vmr_h2o = zeros(Float64, size(T_rshp))
    vcd_dry = zeros(Float64, size(T_rshp))
    vcd_h2o = zeros(Float64, size(T_rshp))
	τ       = zeros(Float64, (len, length(ν_grid)))

	@info "Starting"
	@progress "Processing data..." for i in 1:len
		p_slice = p_rshp[i, :];
		p_ave_slice = p_ave_rshp[i,:];
		T_slice = T_rshp[i, :];
		q_slice = q_rshp[i, :];
		
	    # Apply the function
	    vcd_dry_tmp, vcd_h2o_tmp, vmr_h2o_tmp = layer_VCD(p_slice, q_slice, n_layers)

		if(length(p_ave_slice)==length(T_slice))
			# get spectra & optical depth
			xSec_slice_o2 = [sitp_o2(ν_grid, ip/100, it) for (ip,it) in zip(p_ave_slice, T_slice)];
			xSec_tmp_o2   = hcat(xSec_slice_o2...)
			τ_o2_tmp  = xSec_tmp_o2 * vcd_dry_tmp * vmr_o2;
	
			
			xSec_slice_h2o = [sitp_h2o(ν_grid, ip/100,it) for (ip,it) in zip(p_ave_slice, T_slice)];
			xSec_tmp_h2o   = hcat(xSec_slice_h2o...)
			τ_h2o_tmp  = xSec_tmp_h2o * vcd_h2o_tmp;
			
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
  ╠═╡ =#

# ╔═╡ f3c92079-0058-487f-bb72-f2fa6824c743
md"""
> #### T, p profile => determine a proper grid
"""

# ╔═╡ 85ebaadc-84d9-4ff8-b2bf-b96672fc52cf
begin
	plot(p_ave_rshp[150:50:650, :]', -(1:n_layers))
end

# ╔═╡ c7fd0ad9-1f84-4729-a450-65c61a5ffbed
begin
	plot(T_rshp[210:50:680, :]', -(1:n_layers))
end

# ╔═╡ 1d896f69-82d8-482a-8754-81a6c00dae96
begin
	# upper layer: evaluate xSec at some specific P value
	seg = 42;
	@show mean(p_ave_rshp[:, 1:seg], dims=1);
	@show median(p_ave_rshp[:, 1:seg], dims=1);
	@show std(p_ave_rshp[:, 1:seg], dims=1);
end

# ╔═╡ 8d01ab2e-3adb-4a65-be72-d039d7fc3a4a
begin
	# lower level:
	@show mean(p_ave_rshp[:, seg:end], dims=1);
	@show median(p_ave_rshp[:, seg:end], dims=1);
	@show std(p_ave_rshp[:, seg:end], dims=1);
end

# ╔═╡ 05873d7b-b34c-4c7e-9264-17fa3055df03
md"""
> ### Now make a pressure JLD2 file to evaluate xSec at these specific values
"""

# ╔═╡ 3af49e85-d110-40da-9886-f1f5cbdb4925
a = mean(p_ave_rshp[:, 1:seg], dims=1)

# ╔═╡ a7412bf0-f885-4fd1-b73a-40bf38e8c559
b = mean(p_rshp[:, 1:seg+1], dims=1)

# ╔═╡ d694d5c3-bd85-4d10-8892-3d4a3e89c902
c = Float64.(hcat(a,b))

# ╔═╡ 7e18f981-1191-4869-b8a2-b29acff7f6ea
# ╠═╡ disabled = true
#=╠═╡
# save
@save "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Pressure_UpperAtm.jld2" p_grid
  ╠═╡ =#

# ╔═╡ daa5e0f7-ab72-4bdb-b769-2f927680299f
# ╠═╡ disabled = true
#=╠═╡
begin
	ν_grid = o2_model.ν_grid;
	p_grid = o2_model.p_grid;
	t_grid = o2_model.t_grid;
	itp    = o2_model.itp;
	@show ν_grid, p_grid, t_grid
	# scale
	sitp = scale(itp, ν_grid, p_grid, t_grid);
	println("scaled!")
end
  ╠═╡ =#

# ╔═╡ 291d54ef-cfa5-4eb1-9027-80678efee8c8
#=╠═╡
p_grid = sort(c[:])
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═bc2971ec-5768-11f0-3648-c1692dc3b9da
# ╠═87be1d01-00fc-41b6-9ce8-5d0a4518c807
# ╠═9dd18d77-3d0f-4665-9dde-fe76d8923fa4
# ╠═1e4371af-86fd-41b0-b0e2-432a7c39fc6e
# ╠═a2a8261e-0161-42ad-87f4-aa0c461cabf5
# ╠═c7c5d86b-001a-4337-9e56-6821a6f7298b
# ╠═80a94e86-fd79-40c8-b058-ce21de36918f
# ╟─7503ae08-6dcc-45c9-95a9-ae9ace863e66
# ╠═9458cdfd-3a36-4b9d-aee4-f5794bf4e70d
# ╠═fa75b745-dbc6-44c8-aa0f-4a8c5c21d9a6
# ╠═ed7e418f-500d-4289-81a2-5468c5a963db
# ╟─5d4bb7fa-3c60-4a06-8757-d72848a01464
# ╟─cc4c4ddf-805e-4c78-ba5d-baf7743be139
# ╠═be766ea1-94b1-4ec1-9f6d-30a2da1b8faa
# ╠═a844cb2f-f892-46c6-8825-f1fd35767e46
# ╠═daa5e0f7-ab72-4bdb-b769-2f927680299f
# ╠═1a73db90-e275-4c36-a793-697a5e8ed9ce
# ╠═8fe5df77-2ae2-4ba6-a239-3c927d55920e
# ╟─7c3d7f47-83bc-48bc-ad21-d49209094468
# ╠═85b0a780-ee31-4491-a914-7d64aad289ea
# ╠═23892322-e14a-4155-adbf-6e7a00cf5a93
# ╠═18a3a4d8-d6e0-441c-9957-ee36c825b012
# ╠═aedd5952-3897-43c0-9b77-e4e827f5bb9b
# ╟─f3c92079-0058-487f-bb72-f2fa6824c743
# ╠═85ebaadc-84d9-4ff8-b2bf-b96672fc52cf
# ╠═c7fd0ad9-1f84-4729-a450-65c61a5ffbed
# ╠═1d896f69-82d8-482a-8754-81a6c00dae96
# ╠═8d01ab2e-3adb-4a65-be72-d039d7fc3a4a
# ╟─05873d7b-b34c-4c7e-9264-17fa3055df03
# ╠═3af49e85-d110-40da-9886-f1f5cbdb4925
# ╠═a7412bf0-f885-4fd1-b73a-40bf38e8c559
# ╠═d694d5c3-bd85-4d10-8892-3d4a3e89c902
# ╠═291d54ef-cfa5-4eb1-9027-80678efee8c8
# ╠═a35441c6-7551-419c-8934-013dd35a8277
# ╠═7e18f981-1191-4869-b8a2-b29acff7f6ea
