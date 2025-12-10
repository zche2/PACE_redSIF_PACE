### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ afe34bfe-d19f-11f0-0f09-25a407e3deed
begin
	# Activate project environment
	import Pkg
	Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")
	
	# Load packages
	using StatsBase
	using JLD2, Interpolations, Revise
	using Base.Threads, Dates
	using BlockDiagonals, LinearAlgebra
	using ForwardDiff, DiffResults, Plots, DelimitedFiles, NCDatasets, Statistics
	using Polynomials, Random
	using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools
	using PACE_SIF
end

# ╔═╡ f122b6a6-b5ec-4cb7-b7df-d9b1c23f95c7
begin
	# load retrieval data
	filename = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/retrieval_from_realData/retrieval_for_sample_granule_20250808T204353_new_chl_20251204_215503.jld2";
	@load filename results ref_nflh ref_Rtoa params configuration
end

# ╔═╡ 13dbdf0a-9a42-4795-bc8a-02c1228725cc
begin
	@info params
	@info configuration
	
	# obtain some parameters
	λ 		= params.λ;
	nSIF    = params.nSIF;
	ThisModel = params.forward_model;
	if_log  = true;
	DecompositionMethod = :SVD;
	SIF_PC1 = params.SIFComp[:, 1];
end

# ╔═╡ 692bfae1-d586-4100-8bae-461bd80185ab
begin
	# reconstruct each components
	n_sample = length(results);
	n_λ      = length(λ);
	# preallocate
	ρ_rec_all   = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	T₁_rec_all  = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	T₂_rec_all  = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	SIF_rec_all = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	residual_all = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	SIF_rec_loading_all = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, nSIF);
	# distribution of the uncertainty
	Ŝ_SIF_loaing_all    = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, nSIF);
	
	# loop
	for i in 1:n_sample
		MyPixel = results[i];
		if ismissing(MyPixel)
			residual_all[i,:] .= missing;
			ρ_rec_all[i,:]    .= missing;
			T₁_rec_all[i,:]   .= missing;
			T₂_rec_all[i,:]   .= missing;
			SIF_rec_all[i,:]  .= missing;
			SIF_rec_loading_all[i,:] .= missing;
			Ŝ_SIF_loaing_all[i,:]    .= missing;
			continue
		end
		_, ρ, T₁, T₂, SIF = forward_model(MyPixel.x, MyPixel, return_components=true, if_log=if_log)
		# residual normalized by standard deviation
		residual_all[i,:] = @. MyPixel.y - MyPixel.R_toa;
		# sd = sqrt.(diag(MyPixel.Sₑ))
		# residual_all[i,:] = residual_all[i,:]./sd;
		ρ_rec_all[i,:]    = ρ;
		T₁_rec_all[i,:]   = T₁;
		T₂_rec_all[i,:]   = T₂;
		SIF_rec_all[i,:]  = SIF;
		SIF_rec_loading_all[i,:] = MyPixel.x[end-nSIF+1:end];
		Ŝ_SIF_loaing_all[i,:] = diag(MyPixel.Ŝ[end-nSIF+1:end, end-nSIF+1:end]);
	end
end

# ╔═╡ f17fe912-fa1b-4f61-959a-e92424031db7
begin
	# compare SIF @ 678 nm
	# ind = argmin(abs.(λ .- 678.2));
	# compare peak SIF against nFLH
	ind = argmax(SIF_PC1);
	peak_SIF = SIF_rec_all[:, ind];
	
	# hist
	histogram2d(
		ref_nflh, peak_SIF,
	    bins=200,
	    xlabel="nfLH",
		ylabel="Retrieved SIF",
	    title="Retrieved SIF @ peak wavelength (~$(λ[ind]) nm)",
	    colorbar_title="Count",
	)
	# add 1:1 plot
	plot!([0.1, 0.7], [0.1, 0.7], color=:silver, ls=:dash)
end

# ╔═╡ d778b872-e444-412d-b0d6-609862757709
begin
	∆n = 5000
	colors = palette(:tab10)
	p = plot(
		size=(800, 250), title="T₁ & T₂", legend=false, dpi=300,
		xticks = (620:10:860, string.(620:10:860)),
	)
	
	for (i, idx) in enumerate(1:∆n:n_sample)
		color = colors[mod1(i, length(colors))]
		plot!(λ, T₁_rec_all[idx,:], color=color)
		plot!(λ, T₂_rec_all[idx,:], ls=:dash, color=color)
	end
	p	
end

# ╔═╡ c167d153-b1e2-4fc6-8582-e7c50d8ccb51
begin
	p1 = plot(
		size=(800, 250), title="reconstructed SIF", legend=false, dpi=300,
		xticks = (620:10:860, string.(620:10:860)),
	)
	
	for (i, idx) in enumerate(1:∆n:n_sample)
		color = colors[mod1(i, length(colors))]
		plot!(λ, SIF_rec_all[idx,:], color=color)
	end
	p1	
end

# ╔═╡ 9471228a-d40e-4ae9-a727-ff06501bf12b
begin
	p2 = plot(
		size=(800, 250), title="residual", legend=false, dpi=300,
		xticks = (620:10:860, string.(620:10:860)),)
	
	for (i, idx) in enumerate(1:∆n:n_sample)
		color = colors[mod1(i, length(colors))]
		plot!(λ, residual_all[idx,:], color=color)
	end
	p2
end

# ╔═╡ 897e39a9-f1f7-4a1b-8f75-0d6a482e2af6
begin
	# overlaid by priori
	var_SIF₁ = params.Sₐ[end-nSIF+1, end-nSIF+1];
	var_SIF₂ = params.Sₐ[end, end];
	
	# distribution of uncertainty
	h3 = histogram(
		log10.(Ŝ_SIF_loaing_all[:, 1]), 
	    bins = 100,
	    xlabel = "log Value",
	    ylabel = "Frequency",
	    title  = "Variance of SIF PC1 (log10(Sₐ)=$(log10(var_SIF₁)))", titlefontsize=10,
	    legend = false,
		margin=8Plots.mm
	)

	h4 = histogram(
		log10.(Ŝ_SIF_loaing_all[:, 2]), 
	    bins = 100,
	    xlabel = "log Value",
	    ylabel = "Frequency",
	    title  = "Variance of SIF PC2 (log10(Sₐ)=$(log10(var_SIF₂)))", titlefontsize=10,
	    legend = false,
		margin=8Plots.mm
	)

	# vline!(h3, [mean(log10.(var_SIF₁))], color=:red, linewidth=2)
	# vline!(h4, [mean(log10.(var_SIF₂))], color=:red, linewidth=2)
	
	plot(h3, h4, layout=(2,1), size=(800, 500))
end

# ╔═╡ Cell order:
# ╠═afe34bfe-d19f-11f0-0f09-25a407e3deed
# ╠═f122b6a6-b5ec-4cb7-b7df-d9b1c23f95c7
# ╠═13dbdf0a-9a42-4795-bc8a-02c1228725cc
# ╠═692bfae1-d586-4100-8bae-461bd80185ab
# ╟─f17fe912-fa1b-4f61-959a-e92424031db7
# ╟─d778b872-e444-412d-b0d6-609862757709
# ╟─c167d153-b1e2-4fc6-8582-e7c50d8ccb51
# ╟─9471228a-d40e-4ae9-a727-ff06501bf12b
# ╟─897e39a9-f1f7-4a1b-8f75-0d6a482e2af6
