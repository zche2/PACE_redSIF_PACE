### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 65a175df-fd44-4776-ab1b-056ee104a46c
import Pkg; Pkg.activate("/home/zhe2/FraLab/PACE_redSIF_PACE")

# ╔═╡ 36ecf93e-f903-4328-8483-e8807816b99d
using PACE_SIF

# ╔═╡ 61a60cb4-a633-420d-a14f-c0da6d8bf6b5
using StatsBase, NCDatasets

# ╔═╡ 20911d61-deaf-422a-99fa-688cc00320fe
using JLD2, Plots

# ╔═╡ 94d4f3c2-517a-4c39-9ffc-049c4fb32be5
using LegendrePolynomials, LinearAlgebra

# ╔═╡ ae62f81e-d074-11f0-3d34-ff4a9b0a1fc1
md"""
## Test NMF to transmittance spectra / optical thickness (log)
Created: 2025-12-03
"""

# ╔═╡ 164ab7ea-8d72-424f-ba2a-151741d39f2a
md"""
##### load spectra
"""

# ╔═╡ 285bf434-11c5-4165-b89e-281538981829
rank  = 20;

# ╔═╡ d011a97a-e571-4297-ae22-7e43151b738e
begin
	#====== Wavelength range ======#
	λ_min = 620.0
	λ_max = 860.0
	
	λ_remove_min = 750.0
	λ_remove_max = 749.0
	
	#====== Pseudo observation generation ======#
	n_sample = 5000
	random_seed = 512
	
	# Polynomial fitting
	order = 4
	
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
	nPC   = rank;
	nSIF  = 2;
	nIter = 25;
	thr_Converge = 1e-6;
	
	println("\n=== Configuration ===")
	println("Wavelength range: $λ_min - $λ_max nm")
	println("SNR degradation range: $λ_remove_min - $λ_remove_max nm")
	println("Polynomial order: $order")
	println("SIF scale factor: $scale_factor_SIF")
	println("Number of samples: $n_sample")
	println("Random seed: $random_seed")
	println("NFLH threshold: $nflh_threshold")
	println("Decomposition method: $DecompositionMethod with if_log=$if_log")
	println("NMF rank (effective only when method=NMF): $rank")
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
	
	close(oci)
	println("PACE data loaded: Band selected to [$λ_min, $λ_max] nm")
end

# ╔═╡ 8e54ee43-a51d-4741-b457-ee4358c7da71
md"""
##### NMF to transmittance
"""

# ╔═╡ a44bc4c4-127d-4c7f-abf6-15255d2cf034
begin
	HighResNMF = Spectral_NMF(
        trans, 
        bands,
        Float64.(collect(skipmissing(oci_band))); 
        rank=rank
    );
    # W and H
    λ₀ = HighResNMF.band;
    W₀ = HighResNMF.Loading;
    H₀ = HighResNMF.PrinComp;

    # matrics
    mean_val  = [round(mean(W₀[:, i]), digits=2) for i in 1:rank];
    max_val   = [round(maximum(W₀[:, i]), digits=2) for i in 1:rank];
    min_val   = [round(minimum(W₀[:, i]), digits=2) for i in 1:rank];

    # s.d. for the loading term
    loading_ave_trans = [mean(W₀[:, i]) for i in 1:rank];
    loading_var_trans = [var(W₀[:, i]) for i in 1:rank];
    cov_matx          = cov(W₀, dims=1);

    # pass PrinComp
    PrinComp = H₀';
    println("NMF decomposition complete!")
end

# ╔═╡ a450e94b-4b49-4cce-b149-05670ea3f0af
plot(bands, trans[1:500:end, :]', label=false, size=(800, 300))

# ╔═╡ 1d78b744-398a-4764-a337-1da80cc3f4d5
begin
	n_rows = ceil(Int, rank / 3)
	p = plot(layout=(n_rows, 3), size=(1200, 150*n_rows), dpi=300)
	
	for i in 1:rank
	    pc = PrinComp[:, i]
	    avg_loading = round(loading_ave_trans[i], sigdigits=3)
	    
	    plot!(p[i], λ₀, pc,
	        title="PC $i | μ = $avg_loading",
	        xlabel=(i > rank-2 ? "Wavelength (nm)" : ""),
	        ylabel="Amplitude",
	        legend=false,
	        linewidth=2,
	        color=:steelblue,
	        titlefontsize=9
	    )
	end
	
	if rank % 3 == 1
	    plot!(p[rank+1], axis=false, grid=false)
	end
	
	p
end

# ╔═╡ 18e8ce85-9a4d-43eb-b98f-ab3ddbeaad84
md"""
##### NMF to -log(trans)
"""

# ╔═╡ da26ddb5-6900-425a-a83c-bfa331e13969
begin
	HighResNMF_tau = Spectral_NMF(
        -log.(trans), 
        bands,
        Float64.(collect(skipmissing(oci_band))); 
        rank=rank
    );
    # W and H
    λ₀_tau = HighResNMF_tau.band;
    W₀_tau = HighResNMF_tau.Loading;
    H₀_tau = HighResNMF_tau.PrinComp;

    # s.d. for the loading term
    loading_ave_trans_tau = [mean(W₀_tau[:, i]) for i in 1:rank];
    loading_var_trans_tau = [var(W₀_tau[:, i]) for i in 1:rank];
    cov_matx_tau          = cov(W₀_tau, dims=1);

    # pass PrinComp
    PrinComp_tau = H₀_tau';
    println("NMF decomposition to optical thickness - complete!")
end

# ╔═╡ 7a5d9550-cb88-461e-b064-21531f12b5d2
plot(bands, -log.(trans[1:500:end, :]'), label=false, size=(800, 300))

# ╔═╡ c8e0451f-15ef-468b-a6b4-4cfd6f19c12a
begin
	p1 = plot(layout=(n_rows, 3), size=(1200, 150*n_rows), dpi=300)
	
	for i in 1:rank
	    pc = PrinComp_tau[:, i]
	    avg_loading = round(loading_ave_trans_tau[i], sigdigits=3)
	    
	    plot!(p1[i], λ₀, pc,
	        title="PC $i | μ = $avg_loading",
	        xlabel=(i > rank-2 ? "Wavelength (nm)" : ""),
	        ylabel="Amplitude",
	        legend=false,
	        linewidth=2,
	        color=:steelblue,
	        titlefontsize=9
	    )
	end
	
	if rank % 3 == 1
	    plot!(p1[rank+1], axis=false, grid=false)
	end
	
	p1
end

# ╔═╡ Cell order:
# ╟─ae62f81e-d074-11f0-3d34-ff4a9b0a1fc1
# ╠═65a175df-fd44-4776-ab1b-056ee104a46c
# ╠═36ecf93e-f903-4328-8483-e8807816b99d
# ╠═61a60cb4-a633-420d-a14f-c0da6d8bf6b5
# ╠═20911d61-deaf-422a-99fa-688cc00320fe
# ╠═94d4f3c2-517a-4c39-9ffc-049c4fb32be5
# ╟─164ab7ea-8d72-424f-ba2a-151741d39f2a
# ╠═285bf434-11c5-4165-b89e-281538981829
# ╟─d011a97a-e571-4297-ae22-7e43151b738e
# ╟─8e54ee43-a51d-4741-b457-ee4358c7da71
# ╟─a44bc4c4-127d-4c7f-abf6-15255d2cf034
# ╟─a450e94b-4b49-4cce-b149-05670ea3f0af
# ╟─1d78b744-398a-4764-a337-1da80cc3f4d5
# ╟─18e8ce85-9a4d-43eb-b98f-ab3ddbeaad84
# ╟─da26ddb5-6900-425a-a83c-bfa331e13969
# ╟─7a5d9550-cb88-461e-b064-21531f12b5d2
# ╟─c8e0451f-15ef-468b-a6b4-4cfd6f19c12a
