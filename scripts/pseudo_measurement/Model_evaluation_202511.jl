### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ cf7c06c9-9b7f-467b-a72c-9fbef4fd00dd
import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")

# ╔═╡ 62ade827-b36b-43be-8368-b1219d158b02
using PACE_SIF

# ╔═╡ 1e7ee390-2fcb-475b-81f9-47d96589fc45
using JLD2, Plots

# ╔═╡ 5784065a-bb43-11f0-3465-4109cdf30bb7
md"""
### Evaluate Retrieval from pseudo-measurements
"""

# ╔═╡ 0043b5a8-7df7-4044-9fef-e1fd874ba9fb
# load data
@load "retrieval_results_v1.jld2" Retrieval_all pseudo_obs_all ρ_all T₁_all T₂_all SIF_all params

# ╔═╡ 0785ec75-0e55-4c7a-856c-47a11f6fef3e
begin
	λ = params.λ;
	n_sample = length(Retrieval_all);
	n_λ      = length(λ);

	# preallocate
	ρ_rec_all   = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	T₁_rec_all  = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	T₂_rec_all  = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	SIF_rec_all = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	residual_all = Matrix{Union{Missing, AbstractFloat}}(undef, n_sample, n_λ);
	# reconstruct
	for i in 1:n_sample
		MyPixel = Retrieval_all[i];
		if ismissing(MyPixel)
			continue
		end
		_, ρ, T₁, T₂, SIF = forward_model(MyPixel.x, MyPixel, return_components=true)
		residual_all[i,:] = @. MyPixel.y - MyPixel.R_toa;
		ρ_rec_all[i,:]    = ρ;
		T₁_rec_all[i,:]   = T₁;
		T₂_rec_all[i,:]   = T₂;
		SIF_rec_all[i,:]  = SIF;
	end
end

# ╔═╡ cbc20a9d-c08e-470b-9cb3-143e8100fc06
Δn = 150;

# ╔═╡ 3c70a6da-0dcf-46b0-9303-62dcc45db133
#show params
println("nPoly, nPC = ($(params.nPoly), $(params.nPC))")

# ╔═╡ d611690f-cf36-4b34-8031-fd4e2354cbc0
md"""
##### SIF
"""

# ╔═╡ adc231eb-8101-4766-93e2-dbcb514759d5
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

# ╔═╡ 814bc0e3-250d-427d-b4ee-885df2b0d041
plot(
	λ, 
	(SIF_all[1:Δn:n_sample, :] .- SIF_rec_all[1:Δn:n_sample, :])',
	size=(800, 300),
	legend=false,
	title="SIF (truth-reconstructed)", titlefontsize=8
)

# ╔═╡ a935fc93-17ac-4ace-a277-e8ca23ca00c1
plot(
	λ, 
	((SIF_all[1:Δn:n_sample, :] .- SIF_rec_all[1:Δn:n_sample, :]) ./ SIF_all[1:Δn:n_sample, :])',
	size=(800, 300),
	legend=false,
	title="SIF: (truth-reconstructed) / truth", titlefontsize=8
)

# ╔═╡ 27105614-d35b-414b-a8ee-ad843b15dceb
md"""
##### T₂
"""

# ╔═╡ da5690a2-ed42-465b-ae13-2b2a32218f24
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

# ╔═╡ 8deeced5-bb89-4284-845d-da0a28fbfd3e
plot(
	λ, 
	(T₂_all[1:Δn:n_sample, :] .- T₂_rec_all[1:Δn:n_sample, :])',
	size=(800, 300),
	legend=false,
	title="T₂ (truth-reconstructed)", titlefontsize=8
)

# ╔═╡ e8d766da-ee45-4c39-a15b-10dd2332263b
md"""
##### T₁
"""

# ╔═╡ d3d2552e-6fc4-4582-834b-ab08db4a9cf6
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

# ╔═╡ 17d4cac8-8208-4bd7-ac54-3e2f2d55255c
plot(
	λ, 
	(T₁_all[1:Δn:n_sample, :] .- T₁_rec_all[1:Δn:n_sample, :])',
	size=(800, 300),
	legend=false,
	title="T₁ (truth-reconstructed)", titlefontsize=8
)

# ╔═╡ d0167c75-0497-4ac8-82b7-31232284fc01
md"""
##### ρ
"""

# ╔═╡ 3549415e-bc3b-42e4-87b2-39ff8946c3ad
plot(
	λ, 
	(ρ_all[1:Δn:n_sample, :] ./ ρ_rec_all[1:Δn:n_sample, :])',
	size=(800, 300),
	legend=false,
	title="ρ (truth/reconstructed)", titlefontsize=8
)

# ╔═╡ b7db293c-8845-4e81-8442-992da119c390
md"""
##### Residual
"""

# ╔═╡ c0c40fb4-d6df-47e8-87e3-8103c7373dbd
plot(
	λ, 
	residual_all[1:Δn:n_sample, :]',
	size=(800, 300),
	legend=false,
	title="residual [$(λ[1]), $(λ[end])] nm", titlefontsize=8,
	xticks = (620:10:860, string.(620:10:860))
)

# ╔═╡ Cell order:
# ╟─5784065a-bb43-11f0-3465-4109cdf30bb7
# ╠═cf7c06c9-9b7f-467b-a72c-9fbef4fd00dd
# ╠═62ade827-b36b-43be-8368-b1219d158b02
# ╠═1e7ee390-2fcb-475b-81f9-47d96589fc45
# ╠═0043b5a8-7df7-4044-9fef-e1fd874ba9fb
# ╠═0785ec75-0e55-4c7a-856c-47a11f6fef3e
# ╠═cbc20a9d-c08e-470b-9cb3-143e8100fc06
# ╟─3c70a6da-0dcf-46b0-9303-62dcc45db133
# ╟─d611690f-cf36-4b34-8031-fd4e2354cbc0
# ╟─adc231eb-8101-4766-93e2-dbcb514759d5
# ╟─814bc0e3-250d-427d-b4ee-885df2b0d041
# ╟─a935fc93-17ac-4ace-a277-e8ca23ca00c1
# ╟─27105614-d35b-414b-a8ee-ad843b15dceb
# ╟─da5690a2-ed42-465b-ae13-2b2a32218f24
# ╟─8deeced5-bb89-4284-845d-da0a28fbfd3e
# ╟─e8d766da-ee45-4c39-a15b-10dd2332263b
# ╟─d3d2552e-6fc4-4582-834b-ab08db4a9cf6
# ╟─17d4cac8-8208-4bd7-ac54-3e2f2d55255c
# ╟─d0167c75-0497-4ac8-82b7-31232284fc01
# ╟─3549415e-bc3b-42e4-87b2-39ff8946c3ad
# ╟─b7db293c-8845-4e81-8442-992da119c390
# ╟─c0c40fb4-d6df-47e8-87e3-8103c7373dbd
