using JLD2
using PACE_SIF
using Plots
using Statistics

sif_path = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/SIF_singular_vector.jld2"
sif = JLD2.load(sif_path)

sif_u = convert.(Float64, sif["SIF_U"])
sif_shapes = convert.(Float64, sif["SIF_shapes"]')
λ_ref = collect(Float64.(sif["SIF_wavelen"]))

# all shapes
plot(λ_ref, sif_shapes[1:10, :]')
# principle components
plot(λ_ref, sif_u[:, 1:4])

# manual SVD
sif_svd = PACE_SIF.Spectral_SVD(sif_shapes, λ_ref, λ_min=350.0, λ_max=2500.0)
plot(sif_svd.PrinComp[:, 1:4])
# compute covariance matrix
cov(sif_svd.Loading[1:4, :], dims=2) .* 5

# normalize shapes to max=1
sif_shapes_norm = sif_shapes ./ maximum(sif_shapes, dims=2)
plot(λ_ref, sif_shapes_norm[1:end, :]')

println(sif_svd.Loading)