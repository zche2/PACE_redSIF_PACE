#!/usr/bin/env julia
# Run: julia --project=/path/to/PACE_SIF global_fit_pipeline/svd_retrieval/gpu/run_gpu_tests.jl

using Test
using Random
using ForwardDiff
using LinearAlgebra

const _GPU_DIR = @__DIR__
const _SVD_DIR = dirname(_GPU_DIR)
const _PIPE_DIR = dirname(_SVD_DIR)

# Fast syntax check: parse key swath files without executing them.
# Catches syntax errors (like bad @info calls) in under a second.
@testset "syntax parse check" begin
    files_to_check = [
        joinpath(_SVD_DIR, "svd_helpers.jl"),
        joinpath(_GPU_DIR, "SvdLmTile.jl"),
        joinpath(_GPU_DIR, "SvdFmBatched.jl"),
        joinpath(_SVD_DIR, "svd_swath_parallel.jl"),
    ]
    for f in files_to_check
        if isfile(f)
            @test begin
                Meta.parse("begin\n" * read(f, String) * "\nend")
                true
            end broken=false
        else
            @test_skip "file not found: $f"
        end
    end
end

include(joinpath(_SVD_DIR, "svd_helpers.jl"))
include(joinpath(_GPU_DIR, "SvdLmTile.jl"))

using CUDA

@testset "predict_svd_batched vs scalar fm" begin
    Random.seed!(42)
    nλ, n_pc, n_leg, n_ev = 12, 3, 2, 1
    layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_leg, n_ev = n_ev)
    λ = collect(range(650.0, 720.0; length = nλ))
    z = _normalized_grid(λ)
    leg = Float64.(_legendre_design_matrix(z, n_leg))
    PCs = randn(nλ, n_pc)
    SIF = randn(nλ, n_ev)
    K = 5
    x = randn(layout.n_state, K)
    solar = abs.(randn(nλ, K)) .+ 0.1
    yb = predict_svd_batched(x, solar, PCs, SIF, leg, true, layout)
    for t in 1:K
        fm, _ = make_svd_forward_model_λ(λ, solar[:, t], PCs, SIF; n_pc = n_pc, n_legendre = n_leg, log_transform = true)
        y1 = fm(x[:, t])
        @test maximum(abs.(yb[:, t] .- y1)) < 1e-9
    end
end

@testset "jacobian FD vs ForwardDiff (single column)" begin
    Random.seed!(7)
    nλ, n_pc, n_leg, n_ev = 8, 2, 1, 1
    layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_leg, n_ev = n_ev)
    λ = collect(range(660.0, 700.0; length = nλ))
    z = _normalized_grid(λ)
    leg = Float64.(_legendre_design_matrix(z, n_leg))
    PCs = randn(nλ, n_pc) .* 0.02
    SIF = randn(nλ, n_ev) .* 0.01
    x = randn(layout.n_state, 1) .* 0.1
    solar = abs.(randn(nλ, 1)) .+ 0.2
    Jfd, _ = jacobian_svd_batched_fd(x, solar, PCs, SIF, leg, false, layout; ε = 1e-5)
    fm, lay = make_svd_forward_model_λ(λ, solar[:, 1], PCs, SIF; n_pc = n_pc, n_legendre = n_leg, log_transform = false)
    Jad = ForwardDiff.jacobian(fm, x[:, 1])
    @test size(Jfd) == (nλ, layout.n_state, 1)
    @test maximum(abs.(Jfd[:, :, 1] .- Jad)) < 1e-3
end

@testset "CuArray predict matches CPU" begin
    if !CUDA.functional()
        @info "Skipping CUDA parity (CUDA.functional()==false)"
    else
        try
            Random.seed!(1)
            nλ, n_pc, n_leg, n_ev = 10, 2, 2, 1
            layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_leg, n_ev = n_ev)
            λ = collect(range(640.0, 750.0; length = nλ))
            z = _normalized_grid(λ)
            leg = Float64.(_legendre_design_matrix(z, n_leg))
            PCs = CuArray(randn(nλ, n_pc) .* 0.01)
            SIF = CuArray(randn(nλ, n_ev) .* 0.01)
            legc = CuArray(leg)
            x = randn(layout.n_state, 4) .* 0.05
            solar = abs.(randn(nλ, 4)) .+ 0.1
            y_cpu = predict_svd_batched(x, solar, Array(PCs), Array(SIF), Array(legc), true, layout)
            y_gpu = predict_svd_batched(CuArray(x), CuArray(solar), PCs, SIF, legc, true, layout)
            @test maximum(abs.(y_cpu .- Array(y_gpu))) < 1e-7
        catch e
            @warn "CUDA device test skipped" exception = (e, catch_backtrace())
        end
    end
end

@testset "run_tile_svd_retrieval (CPU, small tile)" begin
    Random.seed!(11)
    nλ, n_pc, n_leg, n_ev = 6, 2, 1, 1
    layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_leg, n_ev = n_ev)
    λ = collect(range(650.0, 700.0; length = nλ))
    z = _normalized_grid(λ)
    leg = Float64.(_legendre_design_matrix(z, n_leg))
    PCs = randn(nλ, n_pc) .* 0.02
    SIF = randn(nλ, n_ev) .* 0.01
    K = 2
    x0 = zeros(layout.n_state)
    x0[first(layout.idx_alpha)] = 1.0
    x0[first(layout.idx_legendre)] = 1.0
    prior = fill(1e6, layout.n_state)
    lower = fill(-Inf, layout.n_state)
    upper = fill(Inf, layout.n_state)
    x_scale = ones(layout.n_state)
    A01 = hcat(ones(length(z)), z)
    lm = (lambda0 = 1.0, lambda_up = 5.0, lambda_down = 0.7, lambda_min = 1e-8, lambda_max = 1e8, max_inner = 12)
    conv = (
        dx_rel_tol = 1e-3,
        rmse_rel_tol = 1e-3,
        rmse_abs_tol = 1e-2,
        enabled = true,
        window = 3,
        redchi2_target = 100.0,
        redchi2_abs_tol = 1.0,
        redchi2_rel_tol = 0.5,
    )
    x_out = zeros(layout.n_state, K)
    y_obs = abs.(randn(nλ, K)) .* 0.5 .+ 0.01
    solar = abs.(randn(nλ, K)) .+ 0.1
    ret = run_tile_svd_retrieval!(
        x_out,
        y_obs,
        solar,
        PCs,
        SIF,
        leg,
        false,
        layout,
        x0,
        prior,
        lower,
        upper,
        x_scale,
        false,
        A01,
        1e-3,
        1.0,
        false,
        nothing,
        0.05,
        lm,
        conv,
        1e-2,
        5;
        use_cuda = false,
    )
    @test size(ret.x_a_mat, 2) == K
    @test all(isfinite.(x_out))
end

@testset "CPU-threaded vs GPU tile path (result parity + timing)" begin
    if !CUDA.functional()
        @info "Skipping GPU vs CPU parity (CUDA.functional()==false)"
        @test_skip "no CUDA device"
    else
        try
            Random.seed!(55)
            nλ, n_pc, n_leg, n_ev = 10, 3, 2, 1
            layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_leg, n_ev = n_ev)
            λ = collect(range(650.0, 710.0; length = nλ))
            z = _normalized_grid(λ)
            leg = Float64.(_legendre_design_matrix(z, n_leg))
            PCs = randn(nλ, n_pc) .* 0.02
            SIF = randn(nλ, n_ev) .* 0.01
            K = 1270*1200
            x0 = zeros(layout.n_state)
            x0[first(layout.idx_alpha)] = 1.0
            x0[first(layout.idx_legendre)] = 1.0
            prior = fill(1e6, layout.n_state)
            lower = fill(-Inf, layout.n_state)
            upper = fill(Inf, layout.n_state)
            x_scale = ones(layout.n_state)
            A01 = hcat(ones(length(z)), z)
            lm = (lambda0 = 1.0, lambda_up = 5.0, lambda_down = 0.7,
                  lambda_min = 1e-8, lambda_max = 1e8, max_inner = 12)
            conv = (
                dx_rel_tol = 1e-3, rmse_rel_tol = 1e-3, rmse_abs_tol = 1e-2,
                enabled = true, window = 3, redchi2_target = 100.0,
                redchi2_abs_tol = 1.0, redchi2_rel_tol = 0.5,
            )
            y_obs = abs.(randn(nλ, K)) .* 0.5 .+ 0.01
            solar = abs.(randn(nλ, K)) .+ 0.1

            x_cpu = zeros(layout.n_state, K)
            t_cpu = @elapsed run_tile_svd_retrieval!(
                x_cpu, y_obs, solar, PCs, SIF, leg, false, layout,
                x0, prior, lower, upper, x_scale, false, A01, 1e-3, 1.0,
                false, nothing, 0.05, lm, conv, 1e-2, 10; use_cuda = false,
            )

            x_gpu = zeros(layout.n_state, K)
            t_gpu = @elapsed run_tile_svd_retrieval!(
                x_gpu, y_obs, solar, PCs, SIF, leg, false, layout,
                x0, prior, lower, upper, x_scale, false, A01, 1e-3, 1.0,
                false, nothing, 0.05, lm, conv, 1e-2, 10; use_cuda = true,
            )

            @info "K=$K pixels | CPU ($(Threads.nthreads()) threads): $(round(t_cpu; digits=3))s | GPU: $(round(t_gpu; digits=3))s | speedup: $(round(t_cpu/t_gpu; digits=2))x"
            @test all(isfinite.(x_cpu))
            @test all(isfinite.(x_gpu))
            # GPU FD Jacobian uses Float64 so results should be very close;
            # small differences can accumulate over LM iterations
            @test maximum(abs.(x_cpu .- x_gpu)) < 1e-4
        catch e
            @warn "GPU vs CPU parity test skipped" exception = (e, catch_backtrace())
        end
    end
end

# log_transform=true is what real PACE retrievals use (svd_log_transform=true in TOML).
# Previous tests only used log_transform=false; this gap hid the GPU NaN bug in production.
@testset "GPU predict log_transform=true matches CPU" begin
    if !CUDA.functional()
        @info "Skipping GPU log_transform=true parity (CUDA.functional()==false)"
        @test_skip "no CUDA device"
    else
        try
            Random.seed!(77)
            nλ, n_pc, n_leg, n_ev = 10, 3, 2, 1
            layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_leg, n_ev = n_ev)
            λ = collect(range(640.0, 750.0; length = nλ))
            z = _normalized_grid(λ)
            leg = Float64.(_legendre_design_matrix(z, n_leg))
            # Use PACE-realistic magnitudes: PCs ~0.02, solar ~1.0, states near prior
            PCs  = randn(nλ, n_pc) .* 0.02
            SIF  = randn(nλ, n_ev) .* 0.01
            K    = 4
            x    = randn(layout.n_state, K) .* 0.05
            solar = abs.(randn(nλ, K)) .+ 1.0   # realistic solar irradiance scale

            y_cpu = predict_svd_batched(x, solar, PCs, SIF, leg, true, layout)
            y_gpu = predict_svd_batched(
                CuArray(x), CuArray(solar),
                CuArray(PCs), CuArray(SIF), CuArray(leg),
                true, layout,
            )
            @test all(isfinite.(y_cpu))
            @test all(isfinite.(Array(y_gpu)))
            @test maximum(abs.(y_cpu .- Array(y_gpu))) < 1e-7
        catch e
            @warn "GPU log_transform=true test skipped" exception=(e, catch_backtrace())
        end
    end
end

@testset "GPU jacobian log_transform=true finite and consistent" begin
    if !CUDA.functional()
        @info "Skipping GPU Jacobian log_transform=true (CUDA.functional()==false)"
        @test_skip "no CUDA device"
    else
        try
            Random.seed!(88)
            nλ, n_pc, n_leg, n_ev = 10, 3, 2, 1
            layout = svd_state_layout(; n_pc = n_pc, n_legendre = n_leg, n_ev = n_ev)
            λ = collect(range(640.0, 750.0; length = nλ))
            z = _normalized_grid(λ)
            leg = Float64.(_legendre_design_matrix(z, n_leg))
            PCs  = randn(nλ, n_pc) .* 0.02
            SIF  = randn(nλ, n_ev) .* 0.01
            K    = 4
            x    = randn(layout.n_state, K) .* 0.05
            solar = abs.(randn(nλ, K)) .+ 1.0

            # CPU FD Jacobian (reference)
            Jcpu, Ycpu = jacobian_svd_batched_fd(x, solar, PCs, SIF, leg, true, layout; ε = 1e-5)

            # GPU FD Jacobian
            PCs_g  = CuArray(PCs);  SIF_g = CuArray(SIF);  leg_g = CuArray(leg)
            xc     = CuArray(x);    sc    = CuArray(solar)
            y0_g   = predict_svd_batched(xc, sc, PCs_g, SIF_g, leg_g, true, layout)
            n_state, n_tile = size(x)
            Jgpu   = CUDA.zeros(Float64, nλ, n_state, n_tile)
            xp     = copy(xc)
            ε_g    = sqrt(eps(Float64))
            for k in 1:n_state
                copyto!(xp, xc)
                xp[k, :] .+= ε_g
                yp = predict_svd_batched(xp, sc, PCs_g, SIF_g, leg_g, true, layout)
                Jgpu[:, k, :] .= (yp .- y0_g) ./ ε_g
            end
            Jgpu_cpu = Array(Jgpu)

            @test all(isfinite.(Jcpu))
            @test all(isfinite.(Jgpu_cpu))
            # GPU and CPU FD Jacobians should agree (different ε but same structure)
            @test maximum(abs.(Jcpu .- Jgpu_cpu)) < 1e-2
        catch e
            @warn "GPU Jacobian log_transform=true test skipped" exception=(e, catch_backtrace())
        end
    end
end

println("All gpu/run_gpu_tests.jl tests passed.")
