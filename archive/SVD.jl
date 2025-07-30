### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ b024f7c9-59b1-487f-84fd-7fafcc0410da
import Pkg; Pkg.activate("../")

# ╔═╡ 8faca87b-24cf-4fda-bf98-404c8d5fa8c6
using Markdown, InteractiveUtils, Plots, PlutoUI

# ╔═╡ edb0179b-7e0f-480c-8460-0104c6f3342e
using NCDatasets

# ╔═╡ 6a92804a-4d15-402c-9bc4-20701c2b16a0
using LinearAlgebra

# ╔═╡ ed6f0dfc-56d6-4691-885d-769b9c808ecf
md"""
read in NetCDF data
"""

# ╔═╡ 33ccc904-a3e0-488e-b6c0-dacfeaa68e4a
begin
	file = "../sample_data/H2O_transmission.nc";
	ds   = Dataset(file, "r");
	println(ds)
	O2_trans = ds["H2O_trans"][:];
	temp = ds["temperature"][:];
	pres = ds["pressure"][:];
	band = ds["band"][:];

	λ_min = 620.;
	λ_max = 850.;
	ind = findall( λ_min .< band .< λ_max);
	
	mat_reshape = reshape(O2_trans, size(O2_trans, 1), :);
	
end

# ╔═╡ 1107134a-bed4-4721-bc46-90400eb50a1f
# log transform
a = log.(mat_reshape[ind, :]); # mat_reshape[ind, :]；

# ╔═╡ d13c371a-d91b-4f24-884c-c7f58cbb0235
begin
	F = svd(a)

	# Accessing the components:
	U = F.U       # Left singular vectors (m x m)
	S = F.S       # Singular values (1D vector)
	V = F.V       # Right singular vectors (n x n)
	Vt = F.Vt     # Conjugate transpose of V (V')
	
	println("\n--- SVD Components ---")
	println("U (Left Singular Vectors):\n", size(U))
	println("V (Right Singular Vectors):\n", size(V))
	println("Vt (Transpose of V):\n", size(Vt))

	# normalize S
	S_norm = S ./ sum(S)
	println("S (Normalized Singular Values):\n", S_norm[1:20])
end

# ╔═╡ 099ec820-188b-48b7-9d1f-27d132b36070
begin 
	m   = 4;
	PC1 = U[:, 1];
	PC2 = U[:, 2];
	PC3 = U[:, 3];
	loading1 = Vt[1, :];
	loading2 = Vt[2, :];
	loading3 = Vt[3, :];

	p = plot(
		band[ind], PC1, 
		linewidth=2.,
		label="PC1, $(S_norm[1]*100)%", legend=:outerbottom
	)
	for i=2:m
		plot!(p, band[ind], U[:, i],
			linewidth=2.,
			label="PC$i, $(S_norm[i]*100)%",
			alpha=.8)
	end
	title!("principle components in log space\n λ [$(λ_min), $(λ_max)] nm")
	# xlims!(650, 695)
	p
end

# ╔═╡ c3ab4b27-6f9c-4b95-b7d6-96fb6eb368a9
md"""
Is it possible to relate loading with T/p?
"""

# ╔═╡ 77596640-8708-42f3-9d21-06c472493529
begin
	# broadcast
	p_mesh = [p for p in pres, t in temp];
	T_mesh = [t for p in pres, t in temp];
	# reshape
	p_1d = reshape(p_mesh, :);
	T_1d = reshape(T_mesh, :);
end

# ╔═╡ ffb9408c-9e70-4e0b-93f8-4de3cf4b2f0d
begin
	k = 3;
	heatmap(
	    pres, temp, reshape(Vt[k, :], length(pres), :)',
	    colorbar=true,       # Display the colorbar
	    title="loading of PC$k ",
	    xlabel="Pressure (hPa)",
	    ylabel="Temperature (K)",
	    c=:viridis,          # Colormap (e.g., :viridis, :plasma, :jet, :coolwarm)
	)
end

# ╔═╡ Cell order:
# ╠═b024f7c9-59b1-487f-84fd-7fafcc0410da
# ╠═8faca87b-24cf-4fda-bf98-404c8d5fa8c6
# ╠═edb0179b-7e0f-480c-8460-0104c6f3342e
# ╠═6a92804a-4d15-402c-9bc4-20701c2b16a0
# ╟─ed6f0dfc-56d6-4691-885d-769b9c808ecf
# ╟─33ccc904-a3e0-488e-b6c0-dacfeaa68e4a
# ╠═1107134a-bed4-4721-bc46-90400eb50a1f
# ╠═d13c371a-d91b-4f24-884c-c7f58cbb0235
# ╠═099ec820-188b-48b7-9d1f-27d132b36070
# ╟─c3ab4b27-6f9c-4b95-b7d6-96fb6eb368a9
# ╠═77596640-8708-42f3-9d21-06c472493529
# ╠═ffb9408c-9e70-4e0b-93f8-4de3cf4b2f0d
