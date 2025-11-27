### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 48759e8d-15c3-4276-aab1-ebd2b3a7967b
begin
	# Activate project environment
	import Pkg
	Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")
	
	# Load packages
	using JLD2, Interpolations, Revise
	using Base.Threads, Dates
	using BlockDiagonals, LinearAlgebra
	using ForwardDiff, DiffResults, Plots, DelimitedFiles, NCDatasets, Statistics
	using Polynomials, Random
	using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools
	using PACE_SIF
end

# ╔═╡ 768249a7-0333-4331-a1c6-349f172efef2
using StatsBase

# ╔═╡ 2afc3990-caf7-11f0-18ce-7977b6d46c74
md"""
### Single Pixel Retrieval
--- 
2025.11.26
"""

# ╔═╡ 5e2dcdbd-ef97-426d-8195-4ef363c8e14f
md"""
##### Model config.
"""

# ╔═╡ 36a8df7f-3459-440f-886a-24f07da5e34d
begin
	# ===========================================
	# CONFIGURATION PARAMETERS
	# ===========================================
	
	#====== Wavelength range ======#
	λ_min = 620.0
	λ_max = 860.0
	
	λ_remove_min = 750.0
	λ_remove_max = 749.0
	
	#====== Pseudo observation generation ======#
	n_sample = 5000
	random_seed = 512
	
	# Baseline wavelengths for fitting
	λ_bl_ref = [607.99, 610.36, 612.73, 615.14, 617.6, 620.06, 622.53, 669.52, 
	            670.76, 671.99, 673.24, 674.51, 675.73, 676.97, 678.21, 679.45, 
	            754.3, 779.33, 867.11, 869.61, 872.13]
	
	# SIF scaling
	scale_factor_SIF = 20
	
	# File paths
	path_transmittance_summer = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_summer_FineWvResModel_FullRange_Aug01.nc"
	path_transmittance_winter = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/transmittance_winter_FineWvResModel_FullRange_Aug01.nc"
	path_oci = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/sample/sample_granule_20240830T131442_new_chl.nc"
	path_sif_shapes = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/SIF_singular_vector.jld2"
	path_snr = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt"
	
	nflh_threshold = 0.05  # Threshold for valid NFLH
	
	#====== retrieval ======#
	DecompositionMethod = :NMF;    # "NMF" or "SVD"
	if_log = false;                 # whether to do log-SVD for transmittance
	n     = 10;
	ranks = 15;
	nPC   = ranks;
	nSIF  = 2;
	nIter = 25;
	thr_Converge = 1e-6;
	
	println("\n=== Configuration ===")
	println("Wavelength range: $λ_min - $λ_max nm")
	println("SNR degradation range: $λ_remove_min - $λ_remove_max nm")
	println("SIF scale factor: $scale_factor_SIF")
	println("Number of samples: $n_sample")
	println("Random seed: $random_seed")
	println("NFLH threshold: $nflh_threshold")
	println("Decomposition method: $DecompositionMethod with if_log=$if_log")
	println("NMF rank (effective only when method=NMF): $ranks")
	println("Order of polynomials to fit: $n, Number of retrieval PCs: $nPC, SIF PCs: $nSIF")
	println("Number of iterations: $nIter, Convergence threshold: $thr_Converge")
	println("====================\n")
	
	
	# ===========================================
	# Load Data
	# ===========================================
	
	println("Loading data...")
	
	# Load MERRA2 transmittance data
	println("Loading transmittance data...")
	summer = Dataset(path_transmittance_summer)
	winter = Dataset(path_transmittance_winter)
	
	trans = cat(summer["transmittance"][:, :], winter["transmittance"][:, :], dims=1)
	bands = summer["band"][:]
	
	close(summer)
	close(winter)
	println("Transmittance data loaded: $(size(trans, 1)) spectra")
	
	# Load PACE OCI data
	println("Loading PACE OCI data...")
	oci = Dataset(path_oci)
	
	red_band = oci["red_wavelength"][:]
	nflh = oci["nflh"][:, :]
	vza = oci["sensor_zenith"][:, :]
	sza = oci["solar_zenith"][:, :]
	chlor_a = oci["chlor_a"][:, :]
	
	# Select spectral band
	ind = findall(λ_min .< red_band .< λ_max)
	
	E = oci["red_solar_irradiance"][ind]
	R_toa = oci["radiance_red"][:, :, ind]
	oci_band = red_band[ind]

	println("PACE data loaded: Band selected to [$λ_min, $λ_max] nm")
	
	# Interpolate transmittance to OCI bands
	println("Interpolating transmittance to OCI bands...")
	trans_new = zeros(size(trans, 1), length(oci_band))
	
	@threads for i in 1:size(trans, 1)
	    itp_row = LinearInterpolation(bands, trans[i, :], extrapolation_bc=0)
	    trans_new[i, :] = itp_row.(oci_band)
	end
	println("Transmittance interpolated: $(size(trans_new, 1)) spectra")
	
	# Load and interpolate SIF shapes
	println("Loading SIF shapes...")
	SIF_shape_dict = JLD2.load(path_sif_shapes)
	
	# SVD to SIF
	SIF_SVD = Spectral_SVD(
	    SIF_shape_dict["SIF_shapes"]*scale_factor_SIF,
	    SIF_shape_dict["SIF_wavelen"],
	    Float64.(collect(skipmissing(oci_band))),
	    if_log = false
	)
	
	SIF_new = SIF_SVD.PrinComp * SIF_SVD.Loading
	println("SIF shapes SVD completed!")
	
	# Load SNR
	println("Loading SNR data...")
	lines = readlines(path_snr)
	end_header_index = findfirst(x -> x == "/end_header", lines)
	data = readdlm(path_snr, String, skipstart=end_header_index)
	
	FPA = data[:, 1]
	wvlen = parse.(Float64, data[:, 2])
	
	wv_val = (λ_min .< wvlen .< λ_max)
	snr_ind = findall((FPA .== "Red") .& wv_val)
	
	c1 = parse.(Float64, data[snr_ind, 4])
	c2 = parse.(Float64, data[snr_ind, 5])
	
	println("Data loading complete!\n")
	
	# Print summary
	println("Data Summary:")
	println("  Wavelength range: $λ_min - $λ_max nm ($(length(oci_band)) bands)")
	println("  TOA radiance: $(size(R_toa)) pixels")
	println("  Transmittance: $(size(trans_new, 1)) spectra")
	println("  SIF shapes: $(size(SIF_new, 2)) components")
	println("  SNR coefficients: $(length(c1)) bands")
	println()

end

# ╔═╡ c394c120-f74f-4266-b3ce-5e71b9bd1029
md"""
##### Retrieval set-up
"""

# ╔═╡ 75ddcbcd-6e3b-48c6-9783-69a2202d32cb
begin
	# Wavelength setup
	λ = oci_band
	λc = center_wavelength(oci_band)
	
	# Find transmittance baseline
	bl_ind = map(λ_bl_ref -> argmin(abs.(λ .- λ_bl_ref)), λ_bl_ref)

	# ===========================================
	# Retrieval setup
	# ===========================================
	
	#= Principal Components =#
	
	if DecompositionMethod == :NMF
	    # NMF
	    HighResNMF = Spectral_NMF(
	        trans, 
	        bands,
	        Float64.(collect(skipmissing(oci_band))); 
	        rank=ranks
	    );
	    # W and H
	    λ₀ = HighResNMF.band;
	    W₀ = HighResNMF.Loading;
	    H₀ = HighResNMF.PrinComp;
	
	    # matrics
	    mean_val  = [round(mean(W₀[:, i]), digits=2) for i in 1:ranks];
	    max_val   = [round(maximum(W₀[:, i]), digits=2) for i in 1:ranks];
	    min_val   = [round(minimum(W₀[:, i]), digits=2) for i in 1:ranks];
	
	    # s.d. for the loading term
	    loading_ave_trans = [mean(W₀[:, i]) for i in 1:ranks];
	    loading_var_trans = [var(W₀[:, i]) for i in 1:ranks];
	    cov_matx          = cov(W₀, dims=1);
	
	    # pass PrinComp
	    PrinComp = H₀';
	    println("\t NMF decomposition complete!")
	
	elseif DecompositionMethod == :SVD
	    # SVD
	    HighResSVD = Spectral_SVD(
	        Float64.(trans'),
	        bands,
	        Float64.(collect(skipmissing(oci_band))),
	        if_log = if_log
	    )
	
	    # matrics
	    loading_ave_trans = [mean(HighResSVD.Loading[i, :]) for i in 1:nPC];
	    loading_var_trans = [var(HighResSVD.Loading[i, :]) for i in 1:nPC];
	    cov_matx          = cov(HighResSVD.Loading[1:nPC, :], dims=2);
	
	    # pass PrinComp
	    PrinComp = HighResSVD.PrinComp[:, 1:nPC];
	    println("\t SVD decomposition complete!")
	end
	
	# SVD 
	loading_var_sif   = var(SIF_SVD.Loading[1:nSIF,:], dims=2) .* 2 ;  # 2 as a scale factor?
	
	# MakePriori
	Sₐ = BlockDiagonal([
	    diagm(fill(1e10, n+1)),           # ρ block
	    cov_matx,                         # transmittance block  
	    diagm([2.0]),                     # smooth_x block
	    diagm(loading_var_sif[1:nSIF])    # SIF covariance block
	])
	
	println("\t Sₐ updated, with diagonal terms: $(diag(Sₐ))")
	
	# remove bands from retrieval evaluation by manually degrading ther SNR
	c2_modified    = copy(c2);
	wv_degrade_ind = findall((λ .>= λ_remove_min) .& (λ .<= λ_remove_max));
	c2_modified[wv_degrade_ind] .= 1e6;  
	
	# forward model
	forward_model_here = (x, px) -> forward_model(
	    x,
	    px;
	    if_log = if_log,
	    return_components = false
	)
	
	# Create the retrieval parameters
	params = RetrievalParams(
	    # Measurement specific
	    λ  = oci_band,                   # Wavelength array
	    λc = λc, 						 # Centered wavelength
	    λ_bl_ind = bl_ind,               # Baseline band indices
	    E        = E,                    # Solar irradiance
		c₁       = c1, 					 # PACE SNR 
		c₂       = c2_modified, 			       	
	    
	    # Forward model settings
	    forward_model = forward_model_here,
	    nPoly = n,                       # Degree of Legendre polynomial
	    nPC   = nPC,                     # Number of transmittance PCs
	    nSIF  = nSIF,                    # Number of SIF PCs
	    Sₐ    = Sₐ,   					 # Prior covariance
	    βₐ    = loading_ave_trans,       # Prior state
	    PrinComp = PrinComp,             # Principal components
	    SIFComp  = SIF_SVD.PrinComp,     # SIF components
	    
	    # Iteration settings (optional, have defaults)
	    iteration_method = LM_Iteration!,
	    nIter = nIter,
	    thr_Converge = thr_Converge
	)
	println("=========\n Parameters setup complete!")
end

# ╔═╡ 63e4ca2d-9e0b-4ed2-80d1-1474bc846c74
md"""
##### Take a look of the real data

take a smaller step, try single pixel retrieval :)
"""

# ╔═╡ a0df2f33-5a3e-4d2c-9f72-02ea27e08dc8
findall(coalesce.((nflh .> 0.5) .& (nflh .< 0.6), false))

# ╔═╡ b7704aae-be09-4eb8-a2a7-b12c525a68e1
begin
	# choose one pixel
	pixel = 1033; scan = 21;
	@info "lat: ", oci["latitude"][pixel, scan];
	@info "lon: ", oci["longitude"][pixel, scan];

	# get R_toa, sza, vza, nflh (as priori)
	px_Rtoa = R_toa[pixel, scan, :];
	px_sza  = sza[pixel, scan];
	px_vza  = vza[pixel, scan];
	px_nflh = nflh[pixel, scan];
	@info "sza, vza, and nflh", px_sza, px_vza, px_nflh

	# struct pixel
	MyPixel = Retrieval_for_Pixel(
		px_Rtoa,
		px_sza,
		px_vza,
		px_nflh,
		1.0,
		1.0,
		params
	);
	println("==== Retrieval completed!====\n")

	# reconstruct
	_, ρ, T₁, T₂, SIF = forward_model(MyPixel.x, MyPixel, return_components=true, if_log=if_log)
	println("==== reconstruction completed!====\n")
end

# ╔═╡ cdffdefa-b643-40f5-b9ec-d349f28ac093
begin
	# \alpha
	smooth_x  = 10. / (1 + exp( -MyPixel.x[end-MyPixel.nSIF]) ) + 1.;
		# tanh(-MyPixel.x[end-MyPixel.nSIF]) + 2.0;
	
	# visualize
	ylabels = ["R_toa", "ρ", "T", "SIF"];
	components = [
		[MyPixel.R_toa, MyPixel.y], 
		ρ,
		[T₁, T₂],
		SIF
	];
	labels  = [
		["obs" "reconstructed"],
		"", 
		["T₁" "T₂, α=$(round(smooth_x, digits=2))"],
		"nFLH=$(MyPixel.nflh)"
	]
	
	
	p = plot(layout=(4,1), size=(600, 600), link=:x)
	
	for (i, (ylabel, y)) in enumerate(zip(ylabels, components))
	    plot!(p[i], λ, y, 
	        ylabel=ylabel, 
	        xlabel=(i==5 ? "Wavelength (nm)" : ""),
	        label=(ismissing(labels[i]) ? "" : labels[i])
		)
	end
	
	p	
end

# ╔═╡ 95481bb3-deea-418d-8a0c-c5fd83808fff
begin
	# S\hat
	Ŝ    = MyPixel.Ŝ;
	logŜ = log10.(abs.(Ŝ));
	# prior covariance matrix
	logSₐ = log10.(abs.(Sₐ));
	
	cmin = -10.5; cmax=-1.0;
	l = @layout [a b{0.55w}]
	
	hmap = plot(layout=l, size=(1200, 500))
	
	# heatmap!(hmap[1], logSₐ,
	#     title="Prior",
	#     aspect_ratio=:equal,
	#     color=:viridis,
	#     clims=(cmin, cmax),
	#     colorbar=false)
	
	# heatmap!(hmap[2], logŜ,
	#     title="Posterior",
	#     aspect_ratio=:equal,
	#     color=:viridis,
	#     clims=(cmin, cmax),
	#     colorbar=true,
	#     colorbar_title="log₁₀|Cov|",
	#     framestyle=:none)
	
	# hmap

	# Top plot: Prior
	heatmap!(hmap[1], logSₐ,
	    ylabel="Parameter index",
	    title="Prior Covariance", 
	    aspect_ratio=:equal,
	    color=:viridis,
	    clims=(cmin, cmax),
	    colorbar=true
	)
	
	# Bottom plot: Posterior
	heatmap!(hmap[2], logŜ,
	    xlabel="Parameter index",
	    ylabel="Parameter index",
	    title="Posterior Covariance", 
	    aspect_ratio=:equal,
	    color=:viridis,
	    clims=(cmin, cmax),
	    colorbar=true
	)
	
	# Shared colorbar (right side)
	# heatmap!(hmap[3], zeros(2,2),  # Dummy plot for colorbar
	#     color=:viridis,
	#     clims=(cmin, cmax),
	#     colorbar=true,
	#     colorbar_title="log₁₀|Cov|",
	#     framestyle=:none,
	#     grid=false,
	#     showaxis=false
	# )
end

# ╔═╡ f1f2a32d-a4d8-40b7-bd5c-5ab507f8e635
begin
	# shows the distribution of \alpha (SVF)
	SVF = @. (secd(sza) + secd(vza)) / secd(vza);
	SVF_vza0 = @. (secd(sza) + secd(0.)) / secd(0.);

	@info "min: $(minimum(SVF))", "max: $(maximum(SVF))"

	# histogram
	histogram(
		SVF[:],
	    xlabel="(secd(sza) + secd(vza)) / secd(vza)",
	    ylabel="Density",
	    normalize=:pdf,
		size=(600, 200),
	    bins=100,
	    legend=false,
		margin=8Plots.mm
	)
	histogram!(
		SVF_vza0[:],
		normalize=:pdf,
	)
end

# ╔═╡ a98f2bf5-6ccd-4cde-b71d-7c747443f3d2
md"""
##### Loop over all pixels
"""

# ╔═╡ 1732eb9b-c23b-4afa-9bcb-b80b46051af7
begin
	# Find valid indices
	valid_indices = findall(coalesce.((nflh .> 0.1) .& (nflh .< 0.6), false))
	
	# Preallocate results
	results       = Vector{Union{Missing, typeof(MyPixel)}}(
		missing, length(valid_indices));
	retrieved_SIF = Vector{Union{Missing, AbstractFloat}}(
		missing, length(valid_indices));
	truth_SIF = Vector{Union{Missing, AbstractFloat}}(
		missing, length(valid_indices))
	
	# Loop through valid indices only
	for (i, idx) in enumerate(valid_indices)
	    try
	        # Extract data for this pixel
	        R_toa_pixel = R_toa[idx, :]
	        sza_pixel = sza[idx]
	        vza_pixel = vza[idx]
	        nflh_pixel = nflh[idx]
	        chlor_a_pixel = chlor_a[idx]

			truth_SIF[i] = nflh_pixel;
	        
	        # Do retrieval
	        results[i] = Retrieval_for_Pixel(
	            R_toa_pixel, 
	            sza_pixel, 
	            vza_pixel, 
	            nflh_pixel, 
	            chlor_a_pixel,  # chlor_a
				1.0,            # flag
	            params
	        )
			
			# reconstruct & store the peak SIF
			if !ismissing(results[i])
			_, _, _, _, SIF_temp = forward_model(
												results[i].x, 
												results[i], 
												return_components=true, 
												if_log=if_log
												);
				retrieved_SIF[i] = maximum(SIF_temp);
			end
	        
	        if i % 100 == 0
	            println("Processed $i / $(length(valid_indices)) pixels")
	        end
			
	    catch e
	        @warn "Failed at index $idx: $e"
	        results[i] = missing
	    end
	end
	
	println("Completed $(count(!ismissing, results)) / $(length(valid_indices)) retrievals")
end

# ╔═╡ eb519f09-ebea-4fa9-a92a-c7e63eac73cc
begin

	histogram2d(
		truth_SIF, retrieved_SIF,
	    bins=100,
	    xlabel="nfLH",
		ylabel="Retrieved SIF",
	    title="Retrieved SIF @ peak wavelength (~681 nm)",
	    colorbar_title="Count",
	)
	
	# # Calculate 2D kernel density
	# k = kde((truth_SIF, retrieved_SIF))
	# density = pdf(k, truth_SIF, retrieved_SIF)
	
	# # Create density scatter plot
	# scatter(truth_SIF, retrieved_SIF,
	#     zcolor=density,
	#     marker=:circle,
	#     markersize=3,
	#     xlabel="True SIF",
	#     ylabel="Retrieved SIF",
	#     title="SIF Retrieval Density",
	#     colorbar_title="Density",
	#     legend=false,
	#     size=(600, 600)
	# )
	
	# # Add 1:1 line
	plot!([minimum(truth_SIF), maximum(truth_SIF)], 
	      [minimum(truth_SIF), maximum(truth_SIF)],
	      color=:red, linestyle=:dash, linewidth=2, label="1:1")
end

# ╔═╡ Cell order:
# ╟─2afc3990-caf7-11f0-18ce-7977b6d46c74
# ╠═48759e8d-15c3-4276-aab1-ebd2b3a7967b
# ╠═768249a7-0333-4331-a1c6-349f172efef2
# ╟─5e2dcdbd-ef97-426d-8195-4ef363c8e14f
# ╟─36a8df7f-3459-440f-886a-24f07da5e34d
# ╟─c394c120-f74f-4266-b3ce-5e71b9bd1029
# ╟─75ddcbcd-6e3b-48c6-9783-69a2202d32cb
# ╟─63e4ca2d-9e0b-4ed2-80d1-1474bc846c74
# ╠═a0df2f33-5a3e-4d2c-9f72-02ea27e08dc8
# ╠═b7704aae-be09-4eb8-a2a7-b12c525a68e1
# ╟─cdffdefa-b643-40f5-b9ec-d349f28ac093
# ╟─95481bb3-deea-418d-8a0c-c5fd83808fff
# ╟─f1f2a32d-a4d8-40b7-bd5c-5ab507f8e635
# ╟─a98f2bf5-6ccd-4cde-b71d-7c747443f3d2
# ╠═1732eb9b-c23b-4afa-9bcb-b80b46051af7
# ╟─eb519f09-ebea-4fa9-a92a-c7e63eac73cc
