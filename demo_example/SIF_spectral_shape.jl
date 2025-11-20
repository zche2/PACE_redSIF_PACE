### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 7b951624-b362-4c44-b364-157ab5e7373d
begin
	import Pkg
	Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")
	
	# Load packages
	using JLD2, Interpolations, Revise
	using Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics
	using Polynomials, Random
	using PACE_SIF
end

# ╔═╡ 1cd2d4da-c63f-11f0-0848-59692eec694b
md"""
## SVD to SIF shapes, dynamic wavelength range
✍️ 2025-11-20
"""

# ╔═╡ c01835e7-9c80-4e86-891e-5133ed6ca733
begin
	λ_min = 620.0
	λ_max = 860.0
	
	# SIF scaling
	scale_factor_SIF = 20
	
	# File paths
	path_oci = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/sample_granule_20240830T131442_new_chl.nc"
	path_sif_shapes = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/SIF_singular_vector.jld2"

	println("Loading PACE OCI data...")
	oci = Dataset(path_oci)
	red_band = oci["red_wavelength"][:]
	# Select spectral band
	ind = findall(λ_min .< red_band .< λ_max)
	oci_band = red_band[ind]
	println("PACE data loaded: Band selected to [$λ_min, $λ_max] nm")
end

# ╔═╡ a0781372-987f-4052-8743-e8c9b3299169
begin
	# load SIF
	println("Loading SIF shapes...")
	SIF_shape_dict = JLD2.load(path_sif_shapes)

	# Create interpolator
	itp = interpolate(SIF_shape_dict["SIF_shapes"], (BSpline(Linear()), NoInterp()))
	range₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end]
	range₂ = 1:size(itp, 2)
	sitp = scale(itp, range₁, range₂)
	setp0 = extrapolate(sitp, 0)
	
	# Interpolate to OCI bands
	SIF_new = reduce(hcat, [setp0.(oci_band, i) for i in range₂])
	SIF_new *= scale_factor_SIF
	println("SIF shapes scaled by factor $scale_factor_SIF: $(size(SIF_new, 2)) spectra, interpolated to oci bands (n_bands = $(size(SIF_new, 1))).")

end

# ╔═╡ 2e987128-abf1-49a9-b3b6-5e3bd0618c42
begin
	# svd
	F = svd(SIF_new);
    PrinComp = F.U;       # Left singular vectors
    S        = F.S;       # Singular values (1D vector)
    Vt       = F.Vt;      # Right singular vectors 
    VarExp_  = S ./ sum(S) * 100;
	Loading  = diagm(S) * Vt;

	println("Shape of PrinComp: $(size(PrinComp))\n" *
	        "Shape of VarExp: $(size(S))\n"*
			"Shape of Loading: $(size(Vt))"
	)
	
end

# ╔═╡ 62725e91-1c83-445d-a293-7f760424d1bf
begin
	title1 = "SIF"
	p1     = plot(size=(800, 300), title=title1, legend=false)
	plot(p1, oci_band, SIF_new[:,1:10:end], lw=2)
end

# ╔═╡ 05c6cb1a-e3cf-4da8-97e0-91869d44417c
begin
	# wrap it up in a new function `Spectral_SVD`
	function Spectral_SVD1(
	    profile::Matrix{FT},     
		# [wavelength x sample]
	    bandᵢₙ::Vector{FT},
	    bandₒᵤₜ::Vector{FT};
		if_log::Bool = false,
	    ) where {FT <: AbstractFloat}

		n_sample    = size(profile, 2)
		profile_new = zeros(FT, length(bandₒᵤₜ), n_sample) 
				# [wavelength x sample]
    
	    for i in 1:n_sample
	        itpᵢ = LinearInterpolation(bandᵢₙ, profile[:, i], extrapolation_bc=0)
	        profile_new[:, i] = itpᵢ.(bandₒᵤₜ)
	    end

		print("Spectra interpolated to target bands: from $(length(bandᵢₙ)) to $(length(bandₒᵤₜ)).\n")

		# --- SVD ---
		F        = svd(profile_new);
		PrinComp = F.U;    
	    S        = F.S;  
	    Vt       = F.Vt; 
	    VarExp   = S ./ sum(S) * 100;
		Loading  = diagm(S) * Vt;
	
		println("Shape of PrinComp: $(size(PrinComp))\n" *
		        "Shape of VarExp: $(size(S))\n"*
				"Shape of Loading (S x V'): $(size(Vt))"
		)

		# return as a struct
	    return SpectraOfPC(
	        band     = bandₒᵤₜ,
	        PrinComp = PrinComp,
	        VarExp   = VarExp,
	        Loading  = Loading,
	        if_log   = if_log
	    )
	end
	
end

# ╔═╡ dabd0781-165b-4d8b-be7d-7eb909d21af4
SIF_PC = Spectral_SVD1(
	SIF_shape_dict["SIF_shapes"]*scale_factor_SIF, SIF_shape_dict["SIF_wavelen"], 
    Float64.(collect(skipmissing(oci_band)))
)

# ╔═╡ 6a59463e-1d7d-4e36-9ab7-2ee741d28782
begin
	# time to visualize it!
	Title = "PrinComp - percent of variance explained (%)"
	p = plot(size=(800, 300), title=Title)
	plot!(p, oci_band, PrinComp[:,1:3], lw=1, label=VarExp_[1:3]')
	plot!(p, oci_band, SIF_PC.PrinComp[:,1:3], lw=2, ls=:dash)
end

# ╔═╡ 1f2455ac-4cf5-4285-b69b-8bffd936a366
begin
	# reconstruct use the first two PCs
	nPC = 3;
	SIF_recon = SIF_PC.PrinComp[:,1:nPC] * SIF_PC.Loading[1:nPC,:];

	# compare
	∆n     = 10;
	title2 = "SIF vs. reconstructed SIF, nPC=$nPC"
	p2     = plot(size=(800, 300), title=title2, legend=false)
	plot!(p2, oci_band, SIF_new[:,1:∆n:end], lw=1)
	plot!(p2, oci_band, SIF_recon[:,1:∆n:end], ls=:dash, lw=2)
end

# ╔═╡ 86ff4e64-a7c4-4ad5-a7eb-46f5a08f8b37
var(SIF_PC.Loading[1:5,:], dims=2)

# ╔═╡ 65812732-f3b3-4fbd-9fa9-29320eddcf6d
mean(SIF_PC.Loading[1:5,:], dims=2)'

# ╔═╡ b73e00b9-1bd3-4ba2-bae2-71711fb1f493
[1.66685; 0.00876223;  -0.000753521;  0.00012267;  1.24771e-5]

# ╔═╡ Cell order:
# ╟─1cd2d4da-c63f-11f0-0848-59692eec694b
# ╠═7b951624-b362-4c44-b364-157ab5e7373d
# ╠═c01835e7-9c80-4e86-891e-5133ed6ca733
# ╠═a0781372-987f-4052-8743-e8c9b3299169
# ╠═2e987128-abf1-49a9-b3b6-5e3bd0618c42
# ╟─62725e91-1c83-445d-a293-7f760424d1bf
# ╠═05c6cb1a-e3cf-4da8-97e0-91869d44417c
# ╠═dabd0781-165b-4d8b-be7d-7eb909d21af4
# ╟─6a59463e-1d7d-4e36-9ab7-2ee741d28782
# ╟─1f2455ac-4cf5-4285-b69b-8bffd936a366
# ╠═86ff4e64-a7c4-4ad5-a7eb-46f5a08f8b37
# ╠═65812732-f3b3-4fbd-9fa9-29320eddcf6d
# ╠═b73e00b9-1bd3-4ba2-bae2-71711fb1f493
