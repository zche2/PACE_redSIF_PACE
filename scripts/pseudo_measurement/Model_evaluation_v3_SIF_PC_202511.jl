### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 5df61699-9460-4b55-90e8-0a1ba8fc01de
import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")

# ╔═╡ a7e9e7cf-9788-44a0-a81a-3c4ed8424e1d
using PACE_SIF

# ╔═╡ c5b7922a-3a58-44a1-9f4b-032c6bdeceeb
using StatsBase

# ╔═╡ f04030a6-9d76-4d8f-b05e-13812053c442
using JLD2, Plots

# ╔═╡ 53ffd4a4-1b68-4367-8fc2-981ac10fd6cf
using LegendrePolynomials, LinearAlgebra, BlockDiagonals

# ╔═╡ 4ba2eb08-51de-4cef-81bb-399e170038db
md"""
### Evaluate Retrieval from pseudo-measurements
"""

# ╔═╡ a0e0c732-d392-45a1-85a0-1120d673bddf
# load data
@load "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/retrieval_from_pseudoObs/retrieval_results_v3_2_1.jld2" Retrieval_all ground_truth params message

# ╔═╡ 80afb5d5-d9db-4232-8dda-54272516f3c5
println(message)

# ╔═╡ ae141818-d9ac-4af5-937b-db206eb66e07
if_log = true

# ╔═╡ 070867ac-5d08-4fee-b94d-c02c2915513b
begin
	ρ_all   = ground_truth.ρ_all;
	T₂_all  = ground_truth.T₂_all;
	T₁_all  = ground_truth.T₁_all;
	SIF_all = ground_truth.SIF_all;
	SIF_loading_all = ground_truth.SIF_loadings_all;
	println("components loaded")
end

# ╔═╡ a7321f67-d360-4011-912e-bee5d6df95ee
begin
	λ = params.λ;
	n_sample = length(Retrieval_all);
	n_λ      = length(λ);
	nSIF     = params.nSIF;
	# preallocate
	ρ_rec_all   = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	T₁_rec_all  = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	T₂_rec_all  = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	SIF_rec_all = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	residual_all = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	# loadings
	SIF_rec_loading_all = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, nSIF);
	# reconstruct
	for i in 1:n_sample
		MyPixel = Retrieval_all[i];
		if ismissing(MyPixel)
			residual_all[i,:] .= missing;
			ρ_rec_all[i,:]    .= missing;
			T₁_rec_all[i,:]   .= missing;
			T₂_rec_all[i,:]   .= missing;
			SIF_rec_all[i,:]  .= missing;
			SIF_rec_loading_all[i,:] .= missing;
			continue
		end
		_, ρ, T₁, T₂, SIF = forward_model(MyPixel.x, MyPixel, return_components=true, if_log=if_log)
		# residual normalized by standard deviation
		residual_all[i,:] = @. MyPixel.y - MyPixel.R_toa;
		sd = sqrt.(diag(MyPixel.Sₑ))
		# residual_all[i,:] = residual_all[i,:]./sd;
		ρ_rec_all[i,:]    = ρ;
		T₁_rec_all[i,:]   = T₁;
		T₂_rec_all[i,:]   = T₂;
		SIF_rec_all[i,:]  = SIF;
		SIF_rec_loading_all[i,:] = MyPixel.x[end-nSIF+1:end];
	end
end

# ╔═╡ 43661e88-d46e-4043-8a0d-ecaf0896d192
Δn = 100;

# ╔═╡ 48e326f5-f20c-4994-80bb-b7d5a8517876
begin
	# check pseudo measurement
	p = plot(size=(800, 300), legend=false, title="created ground truth measurement")
	for i in 1:Δn:n_sample
		MyPixel = Retrieval_all[i];
		if !ismissing(MyPixel)
			plot!(λ, MyPixel.R_toa)
		end
	end
	p
end

# ╔═╡ 37975e90-11c6-4aa1-84a1-a637b2a1e3b4
md"""
##### SIF
"""

# ╔═╡ 864e68b0-293c-4ce4-a454-374c65a34970
begin
    p_SIF = plot(
			size=(800, 300), legend=false, dpi=300,
			title="SIF (solid: truth, dash: reconstructed)", titlefontsize=8,
			xticks = (620:10:860, string.(620:10:860))
		)
    
    # Define color palette
    colors = palette(:tab10)  # or :default, :Set1, etc.
    
    # Plot true and retrieved with matching colors
    for (i, idx) in enumerate(1:Δn:n_sample)
        color = colors[mod1(i, length(colors))]
        plot!(p_SIF, λ, SIF_all[idx, :], color=color, lw=1)
        plot!(p_SIF, λ, SIF_rec_all[idx, :], color=color, ls=:dash, lw=2)
    end
    
    p_SIF
end

# ╔═╡ d9107ea5-f9c9-417f-ab8f-002c2c121909
plot(
	λ, 
	(SIF_all[1:Δn:n_sample, :] .- SIF_rec_all[1:Δn:n_sample, :])',
	size=(800, 300),
	xticks = (620:10:860, string.(620:10:860)),
	legend=false,
	title="SIF (truth-reconstructed)", titlefontsize=8
)

# ╔═╡ 562ab347-5913-41fd-8dd8-02577ce53157
plot(
	λ, 
	((SIF_all[1:Δn:n_sample, :] .- SIF_rec_all[1:Δn:n_sample, :]) ./ SIF_all[1:Δn:n_sample, :] .* 100)',
	size=(800, 300),
	xticks = (620:10:860, string.(620:10:860)),
	# xlims=(650, 850),
	ylims=(-200, 200),
	legend=false,
	# ylabel="%",
	title="SIF: (truth-reconstructed) / truth (% change)", titlefontsize=8,
	margin=8Plots.mm,
)

# ╔═╡ 2e6935e5-eb36-4231-83fe-7bafec557f7f
# SIF peak ind
ind = argmax(SIF_all[1,:])

# ╔═╡ 1e1db338-7d49-456b-8acc-e2e38d64f70c
begin
	a = (SIF_all .- SIF_rec_all) ./ SIF_all .* 100;
	a = a[:, ind];
	histogram(
		a, 
	    bins = 200,
	    xlabel = "Value",
	    ylabel = "Frequency",
	    title = "Distribution of relative error (%) in SIF retrieval @ $(λ[ind]) nm", titlefontsize=8,
		xlims = (-20, 40),
	    legend = false,
		size  = (800, 300),
		margin=8Plots.mm
	)
end

# ╔═╡ 2ea4e56d-c82b-45ed-aa14-efdd4f58a933
begin
	check_SIF_PC = 1;
	x = SIF_loading_all[:,check_SIF_PC];
	y = SIF_rec_loading_all[:, check_SIF_PC];

	# Create vectors of valid pairs only
	valid_pairs = collect(zip(skipmissing(x), skipmissing(y)))
	x_clean = [p[1] for p in valid_pairs]
	y_clean = [p[2] for p in valid_pairs]

	x1 = SIF_loading_all[:,2];
	y1 = SIF_rec_loading_all[:, 2];

	# Create vectors of valid pairs only
	valid_pairs1 = collect(zip(skipmissing(x1), skipmissing(y1)))
	x_clean1 = [p[1] for p in valid_pairs1]
	y_clean1 = [p[2] for p in valid_pairs1]

	h1 = histogram2d(
		x_clean, y_clean,
	    bins=100,
	    xlabel="SIF loading (truth)",
	    ylabel="SIF loading (recon.)",
	    title="PC$(check_SIF_PC)",
	    colorbar_title="Count",
		# ylims=(-0.15, 0.05)
	)

	h2 = histogram2d(
		x_clean1, y_clean1,
	    bins=100,
	    xlabel="SIF loading (truth)",
	    ylabel="SIF loading (recon.)",
	    title="PC2",
	    colorbar_title="Count",
	)

	plot(h1, h2, layout=(1,2), size=(800, 300))
	
end

# ╔═╡ 80c4cade-aa9a-4d82-a9c5-85226e576025
@info size(x_clean)

# ╔═╡ ec2b1bab-9512-4131-bd60-c9c994c3c220
begin
	# posterior covariance matrix
	k = 410;
	if ismissing(Retrieval_all[k])
		@info "missing"
	else
		Ŝ    = Retrieval_all[k].Ŝ;
		logŜ = log10.(abs.(Ŝ));
		cmin = -10.5; cmax=-1.0;
		heatmap(logŜ,
	    	xlabel="Parameter index",
		    ylabel="Parameter index",
		    title="Posterior Covariance (abs, k=$k)", 
		    colorbar_title="Covariance",
		    aspect_ratio=:equal,
		    color=:viridis,
		    clims=(cmin, cmax)  # symmetric color scale
		)
	end
end

# ╔═╡ 7f65e7fa-3273-4855-8ba7-7f3f010f1c5f
begin
	# prior covariance matrix
	local Sₐ, logSₐ  # Declare in outer scope

	try
	    blocks = [Matrix(b) for b in Retrieval_all[k].Sₐ.blocks]
	    Sₐ = Matrix(BlockDiagonal(blocks))
	    logSₐ = log10.(abs.(Sₐ))
	catch e
	    @warn e
	    Sₐ    = Retrieval_all[k].Sₐ
	    logSₐ = log10.(abs.(Sₐ))
	end
	heatmap(logSₐ,
  	    xlabel="state vector",
	    ylabel="state vector",
	    title="priori Covariance", 
	    colorbar_title="Covariance (log)",
	    aspect_ratio=:equal,
	    color=:viridis,
	    clims=(cmin, cmax),  # symmetric color scale
		size=(500,500)
	)
end

# ╔═╡ 1442ba93-305e-4f10-8795-4492c24fb10f
Retrieval_all[k].Sₑ

# ╔═╡ 62785b35-482d-43c2-a89f-346f24e09d6c
begin
	# distribution of the last term of S\hat
	# preallocate
	cov_SIF1   = Vector{Union{Missing, AbstractFloat}}(undef, n_sample);
	cov_SIF2   = Vector{Union{Missing, AbstractFloat}}(undef, n_sample);
	cov_γ      = Vector{Union{Missing, AbstractFloat}}(undef, n_sample);
	
	# reconstruct
	for i in 1:n_sample
		MyPixel = Retrieval_all[i];
		if ismissing(MyPixel)
			cov_SIF1[i] = missing;
			cov_SIF2[i] = missing;
			cov_γ[i]    = missing;
			continue
		end
		cov_SIF1[i] = MyPixel.Ŝ[(end-1),(end-1)];
		cov_SIF2[i] = MyPixel.Ŝ[end,end];
		cov_γ[i]    = MyPixel.Ŝ[(end-2),(end-2)];
	end
	
	h3 = histogram(
		log10.(cov_SIF1), 
	    bins = 100,
	    xlabel = "log Value",
	    ylabel = "Frequency",
	    title  = "Variance of SIF PC1", titlefontsize=10,
	    legend = false,
		margin=8Plots.mm
	)

	h4 = histogram(
		log10.(cov_SIF2), 
	    bins = 100,
	    xlabel = "log Value",
	    ylabel = "Frequency",
	    title  = "Variance of SIF PC2", titlefontsize=10,
	    legend = false,
		margin=8Plots.mm
	)

	h5 = histogram(
		log10.(cov_γ), 
	    bins = 100,
	    xlabel = "log Value",
	    ylabel = "Frequency",
	    title  = "Variance of the factor α", titlefontsize=10,
	    legend = false,
		margin=8Plots.mm
	)

	plot(h3,h4,h5, layout=(3,1), size=(800, 800))
	
end

# ╔═╡ 14ca76c0-1f90-4606-8ad4-e84b1fe578d3
md"""
##### Reconstructed PCs
"""

# ╔═╡ c09df677-540c-497c-9d2a-bf8b68b47d6a
begin
	path_sif_shapes = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/SIF_singular_vector.jld2";
	scale_factor_SIF = 20;
	# Load and interpolate SIF shapes
	println("Loading SIF shapes...")
	SIF_shape_dict = JLD2.load(path_sif_shapes)
	
	# SVD to SIF
	SIF_SVD = Spectral_SVD(
	    SIF_shape_dict["SIF_shapes"]*scale_factor_SIF,
	    SIF_shape_dict["SIF_wavelen"],
	    Float64.(collect(skipmissing(λ))),
	    if_log = false
	)
	
	SIF_new = SIF_SVD.PrinComp * SIF_SVD.Loading
	println("SIF shapes SVD completed!")
end

# ╔═╡ 09baea12-c86c-4b6e-bfac-047262d5eb6a
begin
	SIF_PC1_true  = SIF_loading_all[:,1] * SIF_SVD.PrinComp[:,1]';
	SIF_PC1_recon = SIF_rec_loading_all[:,1] * SIF_SVD.PrinComp[:,1]';

	p_SIF_PC1 = plot(
			size=(800, 300), legend=false, dpi=300,
			title="PC₁ (solid: truth, dash: reconstructed)", titlefontsize=8
		)

    # Plot true and retrieved with matching colors
    for (i, idx) in enumerate(1:Δn:n_sample)
        color = colors[mod1(i, length(colors))]
        plot!(p_SIF_PC1, λ, SIF_PC1_true[idx, :], color=color, lw=1)
        plot!(p_SIF_PC1, λ, SIF_PC1_recon[idx, :], color=color, ls=:dash, lw=2)
    end
    
    p_SIF_PC1
end

# ╔═╡ 87114d37-1035-48d5-9045-f054fad0462a
begin
	SIF_PC2_true  = SIF_loading_all[:,2] * SIF_SVD.PrinComp[:,2]';
	SIF_PC2_recon = SIF_rec_loading_all[:,2] * SIF_SVD.PrinComp[:,2]';

	p_SIF_PC2 = plot(
			size=(800, 300), legend=false, dpi=300,
			title="PC₂ (solid: truth, dash: reconstructed)", titlefontsize=8,
			xticks = (620:10:860, string.(620:10:860))
		)

    # Plot true and retrieved with matching colors
    for (i, idx) in enumerate(1:Δn:n_sample)
        color = colors[mod1(i, length(colors))]
        plot!(p_SIF_PC2, λ, SIF_PC2_true[idx, :], color=color, lw=1)
        plot!(p_SIF_PC2, λ, SIF_PC2_recon[idx, :], color=color, ls=:dash, lw=2)
    end
    
    p_SIF_PC2
end

# ╔═╡ 8c97ca1e-eabd-4b72-9457-0f6a0dacd877
md"""
##### T₂
"""

# ╔═╡ 58fc7afa-bdb2-438d-96da-53ce752e2f81
begin
    p_T₂ = plot(
			size=(800, 300), legend=false, dpi=300,
			title="T₂(solid: truth, dash: reconstructed)", titlefontsize=8
		)

    # Plot true and retrieved with matching colors
    for (i, idx) in enumerate(1:Δn:n_sample)
        color = colors[mod1(i, length(colors))]
        plot!(p_T₂, λ, T₂_all[idx, :], color=color, lw=1)
        plot!(p_T₂, λ, T₂_rec_all[idx, :], color=color, ls=:dash, lw=2)
    end
    
    p_T₂
end

# ╔═╡ 93d648e0-21de-4893-9d81-df5979be333f
plot(
	λ, 
	(T₂_all[1:Δn:n_sample, :] .- T₂_rec_all[1:Δn:n_sample, :])',
	size=(800, 300),
	legend=false,
	title="T₂ (truth-reconstructed)", titlefontsize=8
)

# ╔═╡ 59d7d7db-e7d3-49d4-acc4-58bea87e8099
md"""
##### T₁
"""

# ╔═╡ cb1bbaf5-698c-4157-95bd-c1a85fae0947
begin
    p_T₁ = plot(
			size=(800, 300), legend=false, dpi=300,
			title="T₁ (solid: truth, dash: reconstructed)", titlefontsize=8
		)

    # Plot true and retrieved with matching colors
    for (i, idx) in enumerate(1:Δn:n_sample)
        color = colors[mod1(i, length(colors))]
        plot!(p_T₁, λ, T₁_all[idx, :], color=color, lw=1)
        plot!(p_T₁, λ, T₁_rec_all[idx, :], color=color, ls=:dash, lw=1)
    end
    
    p_T₁
end

# ╔═╡ 09951a0e-c6f1-4f60-8c42-a47e0dc50196
plot(
	λ, 
	(T₁_all[1:Δn:n_sample, :] .- T₁_rec_all[1:Δn:n_sample, :])',
	size=(800, 300),
	legend=false,
	title="T₁ (truth-reconstructed)", titlefontsize=8,
	xticks = (620:10:860, string.(620:10:860))
)

# ╔═╡ b0146c13-766f-4e5c-a063-3ac5175de4b2
md"""
##### ρ
"""

# ╔═╡ 04200a20-e507-4375-a9ab-15ccf6df9031
begin
    p_ρ = plot(
			size=(800, 300), legend=false, dpi=300,
			title="ρ (solid: truth, dash: reconstructed)", titlefontsize=8
		)

    # Plot true and retrieved with matching colors
    for (i, idx) in enumerate(1:Δn:n_sample)
        color = colors[mod1(i, length(colors))]
        plot!(p_ρ, λ, ρ_all[idx, :], color=color, lw=1)
        plot!(p_ρ, λ, ρ_rec_all[idx, :], color=color, ls=:dash, lw=2)
    end
    
    p_ρ
end

# ╔═╡ 546dc52d-0484-4a70-9a6b-6e31bdd5965a
plot(
	λ, 
	(ρ_all[1:Δn:n_sample, :] ./ ρ_rec_all[1:Δn:n_sample, :])',
	size=(800, 300),
	legend=false,
	title="ρ (truth/reconstructed)", titlefontsize=8
)

# ╔═╡ 186394b5-ca04-4b35-bf0e-642ad4fe613d
begin
    p_RTOA = plot(
			size=(800, 300), legend=false, dpi=300,
			title="R_toa (solid: truth, dash: reconstructed)", titlefontsize=8,
			xticks = (620:10:860, string.(620:10:860))
		)

    
    # Plot true and retrieved with matching colors
    for (i, idx) in enumerate(1:Δn:n_sample)
        color = colors[mod1(i, length(colors))]
		MyPixel = Retrieval_all[i];
		if ismissing(MyPixel)
			continue
		end
        plot!(p_RTOA, λ, MyPixel.y, color=color, lw=1)
        plot!(p_RTOA, λ, MyPixel.R_toa, color=color, ls=:dash, lw=2)
    end
    
    p_RTOA
end

# ╔═╡ 2e528284-7418-4b22-b251-502cdc6d5e08
md"""
##### Residual
"""

# ╔═╡ a2486341-0ad3-475a-aae2-bfdf0d81e91e
plot(
	λ, 
	residual_all[1:Δn:n_sample, :]',
	size=(800, 300),
	legend=false,
	title="residual [$(λ[1]), $(λ[end])] nm", 
	titlefontsize=8,
	# title="residual [$(λ[1]), $(λ[end])] nm (normalized by Se)", 
	xticks = (620:10:860, string.(620:10:860))
)

# ╔═╡ 9e02dfc6-693f-4c42-867e-22426bab1d55
md"""
##### Other test
"""

# ╔═╡ a155418b-b63b-40bf-9a91-f2ed8a18c776
plot(
	λ, 
	(@. log(T₂_all[1:Δn:n_sample, :]) / log(T₁_all[1:Δn:n_sample, :]))',
	size=(800, 200),
	legend=false,
	title="log(T₂)/log(T₁): truth", titlefontsize=8,
	xticks = (620:10:860, string.(620:10:860)),
	ylims=(0.0, 5.0)
)

# ╔═╡ 0c0f9051-80ad-4387-8f42-e88bae3b62ef
plot(
	λ, 
	(@. log(T₂_rec_all[1:Δn:n_sample, :]) / log(T₁_rec_all[1:Δn:n_sample, :]))',
	size   = (800, 200),
	legend = false,
	title  = "log(T₂)/log(T₁): reconstructed", titlefontsize=8,
	xticks = (620:10:860, string.(620:10:860)),
	# ylims  = (1.9, 2.1)
)

# ╔═╡ 2fa0294e-ce86-4d8d-9341-cad12a6ec81a
md"""
##### E (solar irr.)
"""

# ╔═╡ 03414f6b-e693-460e-a310-ecfd43011bc5
plot(
	λ, Retrieval_all[k].E, 
	size=(800, 200),
	xticks = (620:10:860, string.(620:10:860))
)

# ╔═╡ Cell order:
# ╠═5df61699-9460-4b55-90e8-0a1ba8fc01de
# ╠═a7e9e7cf-9788-44a0-a81a-3c4ed8424e1d
# ╠═c5b7922a-3a58-44a1-9f4b-032c6bdeceeb
# ╠═f04030a6-9d76-4d8f-b05e-13812053c442
# ╠═53ffd4a4-1b68-4367-8fc2-981ac10fd6cf
# ╟─4ba2eb08-51de-4cef-81bb-399e170038db
# ╠═a0e0c732-d392-45a1-85a0-1120d673bddf
# ╠═80afb5d5-d9db-4232-8dda-54272516f3c5
# ╠═ae141818-d9ac-4af5-937b-db206eb66e07
# ╠═070867ac-5d08-4fee-b94d-c02c2915513b
# ╠═a7321f67-d360-4011-912e-bee5d6df95ee
# ╠═43661e88-d46e-4043-8a0d-ecaf0896d192
# ╟─48e326f5-f20c-4994-80bb-b7d5a8517876
# ╟─37975e90-11c6-4aa1-84a1-a637b2a1e3b4
# ╟─864e68b0-293c-4ce4-a454-374c65a34970
# ╟─d9107ea5-f9c9-417f-ab8f-002c2c121909
# ╟─562ab347-5913-41fd-8dd8-02577ce53157
# ╟─2e6935e5-eb36-4231-83fe-7bafec557f7f
# ╟─1e1db338-7d49-456b-8acc-e2e38d64f70c
# ╟─2ea4e56d-c82b-45ed-aa14-efdd4f58a933
# ╟─80c4cade-aa9a-4d82-a9c5-85226e576025
# ╟─ec2b1bab-9512-4131-bd60-c9c994c3c220
# ╟─7f65e7fa-3273-4855-8ba7-7f3f010f1c5f
# ╠═1442ba93-305e-4f10-8795-4492c24fb10f
# ╟─62785b35-482d-43c2-a89f-346f24e09d6c
# ╟─14ca76c0-1f90-4606-8ad4-e84b1fe578d3
# ╟─c09df677-540c-497c-9d2a-bf8b68b47d6a
# ╟─09baea12-c86c-4b6e-bfac-047262d5eb6a
# ╟─87114d37-1035-48d5-9045-f054fad0462a
# ╟─8c97ca1e-eabd-4b72-9457-0f6a0dacd877
# ╟─58fc7afa-bdb2-438d-96da-53ce752e2f81
# ╟─93d648e0-21de-4893-9d81-df5979be333f
# ╟─59d7d7db-e7d3-49d4-acc4-58bea87e8099
# ╟─cb1bbaf5-698c-4157-95bd-c1a85fae0947
# ╟─09951a0e-c6f1-4f60-8c42-a47e0dc50196
# ╟─b0146c13-766f-4e5c-a063-3ac5175de4b2
# ╟─04200a20-e507-4375-a9ab-15ccf6df9031
# ╟─546dc52d-0484-4a70-9a6b-6e31bdd5965a
# ╟─186394b5-ca04-4b35-bf0e-642ad4fe613d
# ╟─2e528284-7418-4b22-b251-502cdc6d5e08
# ╟─a2486341-0ad3-475a-aae2-bfdf0d81e91e
# ╟─9e02dfc6-693f-4c42-867e-22426bab1d55
# ╟─a155418b-b63b-40bf-9a91-f2ed8a18c776
# ╟─0c0f9051-80ad-4387-8f42-e88bae3b62ef
# ╟─2fa0294e-ce86-4d8d-9341-cad12a6ec81a
# ╠═03414f6b-e693-460e-a310-ecfd43011bc5
