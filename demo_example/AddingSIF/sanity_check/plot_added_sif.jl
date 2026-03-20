using JLD2
using Plots

sif_path = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/reference_spectra/SIF_singular_vector.jld2"
sif = JLD2.load(sif_path)

sif_u = convert.(Float64, sif["SIF_U"])
sif_shapes = convert.(Float64, sif["SIF_shapes"])
λ_ref = collect(Float64.(sif["SIF_wavelen"]))

# plot the first 4 shapes
plot(λ_ref, sif_shapes[:, 1:4])

# scale the shapes to max=0.5
n_bands = size(sif_shapes, 1)
n_shapes = size(sif_shapes, 2)

sif_shapes_lres = Vector{Vector{Float64}}(undef, n_shapes)
for i in 1:n_shapes
    s_val = maximum(abs.(sif_shapes[:, i]))
    sif_shapes_lres[i] = sif_shapes[:, i] ./ s_val .* 0.5
end

# plot the first 4 shapes
plot(λ_ref, sif_shapes_lres[1:14:end])

# alternatively, scale by λ=678 nm
sif_shapes_lres_678 = Vector{Vector{Float64}}(undef, n_shapes)
idx_678 = argmin(abs.(λ_ref .- 678.2))
for i in 1:n_shapes
    s_lres = sif_shapes[:, i]
    s_val = abs(s_lres[idx_678] > 0 ? s_lres[idx_678] : maximum(abs.(s_lres)))
    sif_shapes_lres_678[i] = sif_shapes[:, i] ./ s_val .* 0.5
end

# plot the first 4 shapes
plot(λ_ref, sif_shapes_lres_678[1:14:end])

# add noise to the shapes
sif_shapes_lres_678_noisy = Vector{Vector{Float64}}(undef, n_shapes)
for i in 1:n_shapes
    sif_shapes_lres_678_noisy[i] = sif_shapes_lres_678[i] .+ randn(n_bands) .* 0.01
end

# plot the first 4 shapes
plot(λ_ref, sif_shapes_lres_678_noisy[1:14:end])
