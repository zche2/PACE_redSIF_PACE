### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 9ea9d0f6-a9e5-11f0-2550-d9b7c5e15ac9
import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE");

# ╔═╡ 34913258-2c71-42aa-a6d6-c8670a8390f3
using Polynomials, ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics

# ╔═╡ 6d49da55-aee1-4b05-8f99-f5c706c30bb2
using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools

# ╔═╡ e3fe55f3-9b43-4b3d-bf6b-3951df33e5d0
using JLD2, Interpolations

# ╔═╡ 00efc4d3-01e4-494d-af87-02a828506e90
include("/home/zhe2/FraLab/PACE_redSIF_PACE/PACE_SIF.jl")

# ╔═╡ d4c3f543-5fbf-48e7-b832-e340bc4a4ebb
md"""
## Try different transmittance scaling methods / SVD method to avoid incorrect fit!
---
"""

# ╔═╡ 3307ec96-9bbe-42c8-bfaa-7aca72d3c6c1
md"""
##### Find baseline
---
"""

# ╔═╡ b69efca5-94b8-44b7-aeb8-089ab7f8f79a
λ_min = 610.; λ_max = 880.;

# ╔═╡ c017125c-c923-40dc-b4e1-11482a102adc
begin
	# MERRA2 generated
	summer = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc");
	winter = Dataset("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc");
	println("Opened datasets.")
	
	trans = cat(summer["transmittance"][:,:], winter["transmittance"][:,:], dims=1);
	println("\nConcatenated!")

	bands  = summer["band"][:];
	
	close(summer);
	close(winter);
	
	# SVD
	HighResSVD = PACE_SIF.Spectral_SVD(trans, bands, λ_min=λ_min, λ_max=λ_max);

	# PACE data
	oci = Dataset(
	"/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/sample_granule_20250501T183011_new_chl.nc");
	red_band = oci["red_wavelength"][:];
	nflh     = oci["nflh"][:, :];
	vza      = oci["sensor_zenith"][:, :];
	sza      = oci["solar_zenith"][:, :];
	println("\nRead in PACE Dataset")

	# select band (continuum spectrum)
	ind      = findall( λ_min .< red_band .< λ_max );
	E        = oci["red_solar_irradiance"][ind];
	R_toa    = oci["radiance_red"][:, :, ind];
	oci_band = red_band[ind];
	println("\nBand selected: $oci_band")
end

# ╔═╡ 2ec5f788-cba2-4195-9655-fb746c74cca1
md"""
##### Find wavelengths where there are NO absorption features
---
"""

# ╔═╡ 455963d2-2252-4a1b-becc-c982361efed7
begin
	T̄  = vec(mean(trans, dims=1)');
	# wavelengths where the mean transmittance is high enough!
	thr = 0.9993;
	λᵦ  = findall(T̄ .> thr);
	println(round.(bands[λᵦ], digits=2))
end

# ╔═╡ 247a768f-3746-4bfd-a423-ba78fee8e071
begin
	# plot
	# temp_ind = findall( λ_min .< bands .< λ_max );
	plot(bands, trans[1:400:end, :]',
		label="",
		size=(600, 250),
		lw=.3,
		xlabel="wavelength [nm]",
		ylabel="T"
	)
	plot!(bands, T̄, label="average", lw=2, color=:black)
	scatter!(bands[λᵦ], T̄[λᵦ],
		label="baseline pts > $thr",
		markersize=2,
	)
	# xlims!(λ_min, λ_max)
end

# ╔═╡ 33e59f87-9069-4cdd-a347-158196ef3483
begin
	# a sample transmittance spec to be scaled
	sample₁ = HighResSVD.PrinComp[:, 1:3] * [-4., 1., 0.25];     # bad sample
	sample₂ = HighResSVD.PrinComp[:, 1:3] * [-4., -.1, -0.25];   # good sample
end

# ╔═╡ 143852bd-48f3-4684-acfd-866111203d3b
plot(oci_band, [sample₁, sample₂],
	label=["T₁ - bad" "T₂ - good"],
	size=(600, 250),
	lw=1,
	xlabel="wavelength [nm]",
	ylabel="T"
)

# ╔═╡ 1fe0b727-5eb8-439d-a64a-7a27a7b3edb5
md"""
> 😶‍🌫️ Comment
I don't think by scaling T should work as either a good or bad fit has "flat" baseline, so scaling it only results in SIF magnitude being changed.
"""

# ╔═╡ dee4f0b7-35e1-4f78-94fb-74f53db187cd
md"""
##### Subtract by unity to produce 0 in absence of absorption (Joiner, 2016)
---
"""

# ╔═╡ 4c11cd42-2814-4102-be31-cf5e8b9ebdeb
trans₀ = trans .- 1.0;

# ╔═╡ 5aaf49bd-e2f3-4382-9526-9c7a70090b72
plot(bands, trans₀[1:400:end, :]',
	label="",
	size=(600, 250),
	lw=.3,
	title="T subtracted by 1",
	titlefontsize=10
)

# ╔═╡ 39d759b2-42f1-4e9d-9712-f09cc20dc29b
HighResSVD₀ = PACE_SIF.Spectral_SVD(trans₀, bands, λ_min=λ_min, λ_max=λ_max);

# ╔═╡ 91384f10-7017-4886-b417-e9983882b5e7
std(HighResSVD₀.VarExp .* HighResSVD₀.Loading, dims=2)

# ╔═╡ a44b84e3-48de-4b7f-9cae-fc311baafdde
minimum(HighResSVD₀.VarExp .* HighResSVD₀.Loading, dims=2)

# ╔═╡ a670149e-74d2-49bd-aee3-ee2f0d5d49b5
md"""
##### Or alternatively, log transformation
---
"""

# ╔═╡ a8560b65-104f-43be-a0b0-65f1b905bc11
log_trans = log.(trans);

# ╔═╡ ef2bb0b2-db33-4cbc-a117-ea2f51cf0846
HighResSVD_log = PACE_SIF.Spectral_SVD(log_trans, bands, λ_min=λ_min, λ_max=λ_max);

# ╔═╡ 18c03d04-1482-4d59-b96b-d0cd270e61a8
# what do they look like?
begin
	nPC = 10;
	gr()
	MyLayout = (div(nPC, 2), 2);
	
	edit_title2 =  "log SVD";

	p = plot(oci_band, HighResSVD_log.PrinComp[:,1],
			 title = "PC1 ($(round.(HighResSVD_log.VarExp[1], digits=3))%)",
			 subplot = 1,
			 legend=:topleft,
			 layout = MyLayout,
			 linewidth=2,
			 link = :x,
		     titlefontsize=18,
			 size=(1500, 1000),
			 left_margin=15Plots.mm,
			 xticks=:none,
			 dpi=600
	)
	
	for i in 2:nPC
		if i==5
			plot!(p, oci_band, HighResSVD_log.PrinComp[:,i],
			 		 title = "PC$i ($(round.(HighResSVD_log.VarExp[i], digits=3))%)",
					 subplot = i,
					 legend=:topleft,
					 linewidth=2,
					 titlefontsize=18,
					 layout = MyLayout,
					 ylabel="loading",
				     xticks=:none
			)
		elseif ( i==9 ) | ( i==10 )
			plot!(p, oci_band, HighResSVD_log.PrinComp[:,i],
			 		 title = "PC$i ($(round.(HighResSVD_log.VarExp[i], digits=3))%)",
					 subplot = i,
					 legend=:topleft,
					 linewidth=2,
					 titlefontsize=18,
					 layout = MyLayout,
					 xlabel="wavelength [nm]",
				     bottom_margin=10Plots.mm    
			)
		else
			plot!(p, oci_band, HighResSVD_log.PrinComp[:,i],
			 		 title = "PC$i ($(round.(HighResSVD_log.VarExp[i], digits=3))%)",
					 subplot = i,
					 legend=:topleft,
					 linewidth=2,
					 titlefontsize=18,
					 layout = MyLayout,
				     xticks=:none
			)
		end
	end
	plot!(plot_title=edit_title2, 
		  plot_titlefontsize=20, 
		  legend=false,
		  xtickfontsize=15,
		  ytickfontsize=12,
		  xlabelfontsize=15,
		  ylabelfontsize=15,
		  plot_xlabel="Wavelength [nm]",
		  plot_ylabel="Wavelength [nm]",
    )
	current(p)
end

# ╔═╡ 355ce487-1dd6-4abb-ada3-19b5c4b78bb6
# what do they look like?
begin
	edit_title1 =  "SVD (subtracted by 1)" # if_log ? "log SVD" : "SVD"
	# edit_title *= " sampled from $(size(a)[2]) transmittance spectra"

	p1 = plot(oci_band, HighResSVD₀.PrinComp[:,1],
			 title = "PC1 ($(round.(HighResSVD₀.VarExp[1], digits=3))%)",
			 subplot = 1,
			 legend=:topleft,
			 layout = MyLayout,
			 linewidth=2,
			 link = :x,
		     titlefontsize=18,
			 size=(1500, 1000),
			 left_margin=15Plots.mm,
			 xticks=:none,
			 dpi=600
	)
	
	for i in 2:nPC
		if i==5
			plot!(p1, oci_band, HighResSVD₀.PrinComp[:,i],
			 		 title = "PC$i ($(round.(HighResSVD₀.VarExp[i], digits=3))%)",
					 subplot = i,
					 legend=:topleft,
					 linewidth=2,
					 titlefontsize=18,
					 layout = MyLayout,
					 ylabel="loading",
				     xticks=:none
			)
		elseif ( i==9 ) | ( i==10 )
			plot!(p1, oci_band, HighResSVD₀.PrinComp[:,i],
			 		 title = "PC$i ($(round.(HighResSVD₀.VarExp[i], digits=3))%)",
					 subplot = i,
					 legend=:topleft,
					 linewidth=2,
					 titlefontsize=18,
					 layout = MyLayout,
					 xlabel="wavelength [nm]",
				     bottom_margin=10Plots.mm    
			)
		else
			plot!(p1, oci_band, HighResSVD₀.PrinComp[:,i],
			 		 title = "PC$i ($(round.(HighResSVD₀.VarExp[i], digits=3))%)",
					 subplot = i,
					 legend=:topleft,
					 linewidth=2,
					 titlefontsize=18,
					 layout = MyLayout,
				     xticks=:none
			)
		end
	end
	plot!(plot_title=edit_title1, 
		  plot_titlefontsize=20, 
		  legend=false,
		  xtickfontsize=15,
		  ytickfontsize=12,
		  xlabelfontsize=15,
		  ylabelfontsize=15,
		  plot_xlabel="Wavelength [nm]",
		  plot_ylabel="Wavelength [nm]",
    )
	current(p1)
end

# ╔═╡ 08d143c3-e4b1-46a8-9d94-da238269f93a
md"""
> 🤩 Comment
Cool！ There is no need to rescale as the baseline is 0！！
"""

# ╔═╡ c1522108-6f16-42cf-ac43-a0430a68e0c3
md"""
## Baseline Polynomial fit?
---
"""

# ╔═╡ Cell order:
# ╟─d4c3f543-5fbf-48e7-b832-e340bc4a4ebb
# ╠═9ea9d0f6-a9e5-11f0-2550-d9b7c5e15ac9
# ╠═34913258-2c71-42aa-a6d6-c8670a8390f3
# ╠═6d49da55-aee1-4b05-8f99-f5c706c30bb2
# ╠═e3fe55f3-9b43-4b3d-bf6b-3951df33e5d0
# ╠═00efc4d3-01e4-494d-af87-02a828506e90
# ╟─3307ec96-9bbe-42c8-bfaa-7aca72d3c6c1
# ╠═b69efca5-94b8-44b7-aeb8-089ab7f8f79a
# ╠═c017125c-c923-40dc-b4e1-11482a102adc
# ╟─2ec5f788-cba2-4195-9655-fb746c74cca1
# ╠═455963d2-2252-4a1b-becc-c982361efed7
# ╠═247a768f-3746-4bfd-a423-ba78fee8e071
# ╠═33e59f87-9069-4cdd-a347-158196ef3483
# ╟─143852bd-48f3-4684-acfd-866111203d3b
# ╟─1fe0b727-5eb8-439d-a64a-7a27a7b3edb5
# ╟─dee4f0b7-35e1-4f78-94fb-74f53db187cd
# ╠═4c11cd42-2814-4102-be31-cf5e8b9ebdeb
# ╟─5aaf49bd-e2f3-4382-9526-9c7a70090b72
# ╠═39d759b2-42f1-4e9d-9712-f09cc20dc29b
# ╟─355ce487-1dd6-4abb-ada3-19b5c4b78bb6
# ╠═91384f10-7017-4886-b417-e9983882b5e7
# ╠═a44b84e3-48de-4b7f-9cae-fc311baafdde
# ╟─a670149e-74d2-49bd-aee3-ee2f0d5d49b5
# ╠═a8560b65-104f-43be-a0b0-65f1b905bc11
# ╠═ef2bb0b2-db33-4cbc-a117-ea2f51cf0846
# ╠═18c03d04-1482-4d59-b96b-d0cd270e61a8
# ╟─08d143c3-e4b1-46a8-9d94-da238269f93a
# ╟─c1522108-6f16-42cf-ac43-a0430a68e0c3
