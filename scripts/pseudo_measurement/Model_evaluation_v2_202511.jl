### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ cf7c06c9-9b7f-467b-a72c-9fbef4fd00dd
import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")

# ╔═╡ 62ade827-b36b-43be-8368-b1219d158b02
using PACE_SIF

# ╔═╡ 1e7ee390-2fcb-475b-81f9-47d96589fc45
using JLD2, Plots, NCDatasets

# ╔═╡ e5d1f463-2598-4c9f-a8bb-0713138d8791
using LegendrePolynomials

# ╔═╡ 5784065a-bb43-11f0-3465-4109cdf30bb7
md"""
### Evaluate Retrieval from pseudo-measurements
"""

# ╔═╡ 0043b5a8-7df7-4044-9fef-e1fd874ba9fb
# load data
@load "retrieval_results_v2_4.jld2" Retrieval_all pseudo_obs_all ρ_all T₁_all T₂_all SIF_all params_to_save message MyIter

# ╔═╡ 71807df9-dbde-4b49-a956-cb96a8f8e675
println(message)

# ╔═╡ af8ab435-e6ac-46a6-8e1f-4409ce14a2c1
md"""
### How many pixels retrieved successfully?
"""

# ╔═╡ 8e50fcfa-6d8f-4c9e-9a6a-085b96c667d9
begin

	λ_min = 620.0  # Define if not already set
	λ_max = 860.0  # Define if not already set

	
	println("Load LUT for cross sections...")
	o2_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_O2.jld2";
	o2_sitp = read_rescale(o2_jld2);
	h2o_jld2 = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01_H2O.jld2";
	h2o_sitp = read_rescale(h2o_jld2);
	metadata = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/interp_xSection/Finer_Wavenumber_grid_FullRange_Aug01/Finer_Wavenumber_grid_FullRange_Aug01.log";
	ν_grid_o2, p_grid_hPa, t_grid = o2_sitp.ranges;
	println("LUT loaded.")
	
	kernel_dir = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/"
	@load kernel_dir*"KernelInstrument.jld2" MyKernel
	println("Kernel instrument loaded.")
end

# ╔═╡ c5c347e7-c12e-459d-b207-bb674495a5f7
params = RetrievalParams_xSecFit(
    λ  = params_to_save.λ,
    λc = params_to_save.λc,
    E  = params_to_save.E,
    c₁ = params_to_save.c₁,
    c₂ = params_to_save.c₂,
    o2_sitp  = o2_sitp,
    h2o_sitp = h2o_sitp,
    InstrumentKernel = MyKernel,
    forward_model    = forward_model,
    nPoly  = params_to_save.nPoly,
    nLayer = params_to_save.nLayer,
    nSIF   = params_to_save.nSIF,
    Sₐ = params_to_save.Sₐ,
    βₐ = params_to_save.βₐ,
    γₐ = params_to_save.γₐ,
    SIFComp  = params_to_save.SIFComp,      
    iteration_method = MyIter,
    thr_Converge = params_to_save.thr_Converge
);

# ╔═╡ 7a2d63f0-4784-4e15-96ae-dc14c8169d65
MyModel = (x, px) -> forward_model(
	x,
	px, 
	params;
	return_components=true
);

# ╔═╡ 0785ec75-0e55-4c7a-856c-47a11f6fef3e
begin
	λ = params_to_save.λ;
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
			residual_all[i,:] .= missing;
			ρ_rec_all[i,:]    .= missing;
			T₁_rec_all[i,:]   .= missing;
			T₂_rec_all[i,:]   .= missing;
			SIF_rec_all[i,:]  .= missing;
		else
			_, ρ, T₁, T₂, SIF = MyModel(MyPixel.x, MyPixel)
			residual_all[i,:] = @. MyPixel.y - MyPixel.R_toa;
			ρ_rec_all[i,:]    = ρ;
			T₁_rec_all[i,:]   = T₁;
			T₂_rec_all[i,:]   = T₂;
			SIF_rec_all[i,:]  = SIF;
		end
	end
end

# ╔═╡ 52f5f17b-5c42-4913-8fab-e02887d7af1b
k = 60; _, ρ, T₁_k, T₂_k, SIF = MyModel(Retrieval_all[k].x, Retrieval_all[k])

# ╔═╡ fc4e824a-91a3-4328-9263-8c9bcf43acc6
begin
	
	# manual reconstruction of T2
	this_vec  = Retrieval_all[k].x;
	thisPixel = Retrieval_all[k];
	
	# reflectance
    xᵨ    = this_vec[1 : thisPixel.nPoly+1]
    v     = collectPl.(thisPixel.λc, lmax=thisPixel.nPoly);
    thisρ = hcat(v...)' * xᵨ;
	
	# cos(sza) and cos(vza)
    μ₁    = 1 / cosd(thisPixel.sza);
    μ₂    = 1 / cosd(thisPixel.vza);

    # T↑ transmittance for SIF
    thisx₁    = this_vec[(thisPixel.nPoly+2):(thisPixel.nPoly+7)];
    thisT₁    = Retrieval.compute_transmittance(
                thisx₁,
                MyKernel,
                o2_sitp,
                h2o_sitp,
            );

    # T↓↑ transmittance for reflected radiance
    thisx₂    = this_vec[(thisPixel.nPoly+8):(thisPixel.nPoly+13)];
    thisT₂    = Retrieval.compute_transmittance(
                thisx₂,
                MyKernel,
                o2_sitp,
                h2o_sitp,
            );

    # take light path into account
    @. thisT₂ = thisT₂^( μ₁ + μ₂ );
    @. thisT₁ = thisT₁^( μ₂ );

	# SIF
	thisxₛ     = this_vec[end - thisPixel.nSIF + 1 : end];
    thisSIF    = thisPixel.SIF_shape * thisxₛ;

end

# ╔═╡ 1f9b016f-01b1-4dc0-986e-aa58aa48bb73
begin
	plot(
		λ, ρ_all[k,:], size=(800, 300),
		xticks = (620:10:860, string.(620:10:860)),
		title="ρ",
		label="truth", lw=1.,
		)
	plot!(
		λ, ρ, label="reconstructed",  lw=1.,
	);
	plot!(
		λ, thisρ, label="manually reconstructed",  lw=2., ls=:dash, color=:olive
	);
end

# ╔═╡ 44368e86-8b96-4ec3-89fd-8111ad61d780
begin
	plot(
		λ, T₂_all[k,:], size=(800, 300),
		xticks = (620:10:860, string.(620:10:860)),
		title="T₂",
		label="truth", lw=1.,
		)
	plot!(
		λ, T₂_k, label="reconstructed",  lw=1.,
	);
	plot!(
		λ, thisT₂, label="manually reconstructed",  lw=2., ls=:dash, color=:olive
	);
end

# ╔═╡ 3a09dc6e-cf67-43f4-9068-c5c9fb5e93db
begin
	plot(
		λ, T₁_all[k,:], size=(800, 300),
		xticks = (620:10:860, string.(620:10:860)),
		title="T₁",
		label="truth", lw=1.,
		)
	plot!(
		λ, T₁_k, label="reconstructed",  lw=1.,
	);
	plot!(
		λ, thisT₁, label="manually reconstructed",  lw=2., ls=:dash, color=:olive
	);
end

# ╔═╡ 657ac6d3-5d10-4269-86a1-4cbccdb04451
begin
	plot(
		λ, SIF_all[k,:], size=(800, 300),
		xticks = (620:10:860, string.(620:10:860)),
		title="SIF",
		label="truth", lw=1.,
		)
	plot!(
		λ, SIF, label="reconstructed",  lw=1.,
	);
	plot!(
		λ, thisSIF, label="manually reconstructed",  lw=2., ls=:dash, color=:olive
	);
end

# ╔═╡ cbc20a9d-c08e-470b-9cb3-143e8100fc06
@info n_sample ;Δn = 10;

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
	( (SIF_all[1:Δn:n_sample, :] .- SIF_rec_all[1:Δn:n_sample, :]) ./ SIF_all[1:Δn:n_sample, :] .* 100)',
	size=(800, 300),
	legend=false,
	ylabel="%",
	title="SIF: (truth-reconstructed) / truth", titlefontsize=8,
	margin=8Plots.mm
)

# ╔═╡ f2652a7b-b54f-42d8-8a33-08b5822101d1
# SIF peak ind
ind = argmax(SIF_all[1,:])

# ╔═╡ 06940505-ab59-4972-b3d3-d2a480c9361f
begin
	a = - (SIF_all .- SIF_rec_all) ./ SIF_all .* 100;
	a = a[:, ind];
	histogram(
		a, 
	    bins = 100,
	    xlabel = "Value",
	    ylabel = "Frequency",
	    title = "Distribution of relative error (%) in SIF retrieval @ $(λ[ind]) nm", titlefontsize=8,
		# xlims = (-30, 30),
	    legend = false,
		size  = (800, 300),
		margin=8Plots.mm
	)
end

# ╔═╡ 27105614-d35b-414b-a8ee-ad843b15dceb
md"""
##### T₂
"""

# ╔═╡ da5690a2-ed42-465b-ae13-2b2a32218f24
begin
    p_T₂ = plot(
			size=(800, 300), legend=false, dpi=300,
			title="T₂(solid: truth, dash: reconstructed)", titlefontsize=8,
			xticks = (620:10:860, string.(620:10:860))
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
	xticks = (620:10:860, string.(620:10:860)),
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

# ╔═╡ b82e5b1c-1e55-454a-9936-b187061a84b1
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

# ╔═╡ 4735a6bd-7853-47c2-a9b3-c253a2f836a0
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
        plot!(p_RTOA, λ, MyPixel.y, color=color, lw=.5)
        plot!(p_RTOA, λ, MyPixel.R_toa, color=color, ls=:dash, lw=1)
    end
    
    p_RTOA
end

# ╔═╡ 404988a2-14ca-4ec2-876d-8237f7165428
md"""
##### Other test
"""

# ╔═╡ 5fb64951-1f67-405d-ad32-94247cb8807d
plot(
	λ, 
	(@. log(T₂_all[1:Δn:n_sample, :]) / log(T₁_all[1:Δn:n_sample, :]))',
	size=(800, 200),
	legend=false,
	title="log(T₂)/log(T₁): truth", titlefontsize=8,
	xticks = (620:10:860, string.(620:10:860)),
	ylims=(0.0, 135.0)
)

# ╔═╡ 23b72be9-2999-47ec-8d8b-ea58ad120823
plot(
	λ, 
	(@. log(T₂_rec_all[1:Δn:n_sample, :]) / log(T₁_rec_all[1:Δn:n_sample, :]))',
	size   = (800, 200),
	legend = false,
	title  = "log(T₂)/log(T₁): reconstructed", titlefontsize=8,
	xticks = (620:10:860, string.(620:10:860)),
	ylims  = (0., 135.)
)

# ╔═╡ 33147397-cb05-4eb5-ac93-ace9e3cf434e
# save the Kernel (no need to compute everytime)
println("MyKernel: $(Base.summarysize(MyKernel) / 1024^2) MB")

# ╔═╡ a8dfb70f-4fc6-4659-a008-484bed3cd0da
# save_dir = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/"

# ╔═╡ 889d02a2-b1d0-46d1-bc1b-030064a322ef
# @save save_dir*"KernelInstrument.jld2" MyKernel

# ╔═╡ Cell order:
# ╟─5784065a-bb43-11f0-3465-4109cdf30bb7
# ╠═cf7c06c9-9b7f-467b-a72c-9fbef4fd00dd
# ╠═62ade827-b36b-43be-8368-b1219d158b02
# ╠═1e7ee390-2fcb-475b-81f9-47d96589fc45
# ╠═0043b5a8-7df7-4044-9fef-e1fd874ba9fb
# ╠═71807df9-dbde-4b49-a956-cb96a8f8e675
# ╟─af8ab435-e6ac-46a6-8e1f-4409ce14a2c1
# ╠═8e50fcfa-6d8f-4c9e-9a6a-085b96c667d9
# ╠═c5c347e7-c12e-459d-b207-bb674495a5f7
# ╠═7a2d63f0-4784-4e15-96ae-dc14c8169d65
# ╠═0785ec75-0e55-4c7a-856c-47a11f6fef3e
# ╠═52f5f17b-5c42-4913-8fab-e02887d7af1b
# ╠═e5d1f463-2598-4c9f-a8bb-0713138d8791
# ╠═fc4e824a-91a3-4328-9263-8c9bcf43acc6
# ╟─1f9b016f-01b1-4dc0-986e-aa58aa48bb73
# ╟─44368e86-8b96-4ec3-89fd-8111ad61d780
# ╟─3a09dc6e-cf67-43f4-9068-c5c9fb5e93db
# ╟─657ac6d3-5d10-4269-86a1-4cbccdb04451
# ╠═cbc20a9d-c08e-470b-9cb3-143e8100fc06
# ╟─d611690f-cf36-4b34-8031-fd4e2354cbc0
# ╟─adc231eb-8101-4766-93e2-dbcb514759d5
# ╟─814bc0e3-250d-427d-b4ee-885df2b0d041
# ╟─a935fc93-17ac-4ace-a277-e8ca23ca00c1
# ╠═f2652a7b-b54f-42d8-8a33-08b5822101d1
# ╠═06940505-ab59-4972-b3d3-d2a480c9361f
# ╟─27105614-d35b-414b-a8ee-ad843b15dceb
# ╟─da5690a2-ed42-465b-ae13-2b2a32218f24
# ╟─8deeced5-bb89-4284-845d-da0a28fbfd3e
# ╟─e8d766da-ee45-4c39-a15b-10dd2332263b
# ╟─d3d2552e-6fc4-4582-834b-ab08db4a9cf6
# ╟─17d4cac8-8208-4bd7-ac54-3e2f2d55255c
# ╟─d0167c75-0497-4ac8-82b7-31232284fc01
# ╟─b82e5b1c-1e55-454a-9936-b187061a84b1
# ╟─3549415e-bc3b-42e4-87b2-39ff8946c3ad
# ╟─b7db293c-8845-4e81-8442-992da119c390
# ╟─c0c40fb4-d6df-47e8-87e3-8103c7373dbd
# ╟─4735a6bd-7853-47c2-a9b3-c253a2f836a0
# ╟─404988a2-14ca-4ec2-876d-8237f7165428
# ╟─5fb64951-1f67-405d-ad32-94247cb8807d
# ╟─23b72be9-2999-47ec-8d8b-ea58ad120823
# ╠═33147397-cb05-4eb5-ac93-ace9e3cf434e
# ╠═a8dfb70f-4fc6-4659-a008-484bed3cd0da
# ╠═889d02a2-b1d0-46d1-bc1b-030064a322ef
