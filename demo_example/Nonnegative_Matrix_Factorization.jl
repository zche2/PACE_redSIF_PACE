### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ fa9af6ff-5676-4eed-87d5-6dba5bf13a45
import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE");

# ╔═╡ 0a50c2e7-e18f-4a52-8a65-718f2074e953
using Polynomials, ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics, Interpolations

# ╔═╡ 0984787b-f797-41c4-bd97-69f3bdfccc3c
using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools, NMF

# ╔═╡ 05277310-ee51-43ca-9ca8-6e5651054503
using JLD2

# ╔═╡ 9612db3e-ae0d-11f0-1799-f1ca8fae2819
md"""
## Compare SVD and NMF on transmittance spectra
2025-10-20

"""

# ╔═╡ fe96c2dd-799e-497f-9355-94c9376f7953
begin
	# MERRA2 generated
	summer = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc");
	winter = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc");
	println("Opened datasets.")
	
	trans = cat(summer["transmittance"][:,:], winter["transmittance"][:,:], dims=1);
	println("\nConcatenated!")

	bands  = summer["band"][:];

	λ_min = 620.;
	λ_max = 780.;
	ind   = findall( λ_min .< bands .< λ_max);

	# downscale to the specific wavelength range
	wvlen = bands[ind];
	spec  = trans[:, ind];
	
	close(summer);
	close(winter);
end

# ╔═╡ 817fd719-17d2-4aed-8128-4ed6652cbd2b
md"""
#### 1️⃣ NMF
---
"""

# ╔═╡ 01350b23-2c99-42d2-ad4b-8f6002c3e9b0
begin
	k = 20
	
	 # initialize
	W, H = NMF.spa(spec, k)
	
	 # optimize
	NMF.solve!(NMF.SPA{Float64}(obj=:mse), spec, W, H)
end

# ╔═╡ 25bbdfa6-054f-4ae4-8810-c2f3931de0df
begin
    mean_val  = [round(mean(W[:, i]), digits=2) for i in 1:size(W, 2)];
    max_val   = [round(maximum(W[:, i]), digits=2) for i in 1:size(W, 2)];
    min_val   = [round(minimum(W[:, i]), digits=2) for i in 1:size(W, 2)];
    mean_spec = mean(spec, dims=1);

    # Create a plot with k panels (one for each row)
    plot(
        [begin
            p = plot(wvlen, H[i, :], label="",title="$(mean_val[i]) ($(min_val[i]), $(max_val[i]))", lw=2.)
            plot!(p, wvlen, vec(mean_spec), color=:silver, label="", lw=2., alpha=.3)
            p
        end for i in 1:size(H, 1)]..., 
        layout=(k÷2, 2),     # k÷2 rows, 2 columns layout
        size=(1200, 800),     # Adjust the size of the plot
    )
end

# ╔═╡ bf99b6d3-18bf-426b-bddb-f017508c6c8e
begin
	# variance explained by each single component?
	# I don't think this variance make sense, as it is NOT comparable with SVD.
	total_var = norm(spec, 2)^2;
    S_norm = [norm(W[:, i:i] * H[i:i, :])^2 / total_var * 100 for i in 1:k]
end

# ╔═╡ 008f6eec-2647-41f6-a07a-3ff8e4ad5edc
md"""
visualize the change in residual with rank
"""

# ╔═╡ 8f6d3156-14b3-41c5-809c-34501fc72992
begin
	# Perform NMF for different ranks and compute residuals
	ranks = 5:2:20  # Range of ranks to test
	residuals = Float64[]
	
	for r in ranks
	    Wᵢ, Hᵢ   = NMF.spa(spec, r)  
	    residual = norm(spec - Wᵢ * Hᵢ, 2)  # Frobenius norm of the residual
	    push!(residuals, residual)
	end
	
	# Plot the residuals vs rank
	plot(ranks, residuals, 
	     xlabel="Rank", 
	     ylabel="Residual (Frobenius Norm)", 
	     title="Residual vs Rank",
	     marker=:circle,
	     lw=2,
	     legend=false,
		 # yaxis=:log,
	     size=(800, 400)
	)
end

# ╔═╡ Cell order:
# ╟─9612db3e-ae0d-11f0-1799-f1ca8fae2819
# ╠═fa9af6ff-5676-4eed-87d5-6dba5bf13a45
# ╠═0a50c2e7-e18f-4a52-8a65-718f2074e953
# ╠═0984787b-f797-41c4-bd97-69f3bdfccc3c
# ╠═05277310-ee51-43ca-9ca8-6e5651054503
# ╠═fe96c2dd-799e-497f-9355-94c9376f7953
# ╟─817fd719-17d2-4aed-8128-4ed6652cbd2b
# ╠═01350b23-2c99-42d2-ad4b-8f6002c3e9b0
# ╠═25bbdfa6-054f-4ae4-8810-c2f3931de0df
# ╠═bf99b6d3-18bf-426b-bddb-f017508c6c8e
# ╟─008f6eec-2647-41f6-a07a-3ff8e4ad5edc
# ╠═8f6d3156-14b3-41c5-809c-34501fc72992
