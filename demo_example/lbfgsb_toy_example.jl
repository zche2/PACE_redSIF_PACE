#=
  Minimal bound-constrained quasi-Newton demo (L-BFGS-style + box constraints).

  Pure Julia (stdlib + LinearAlgebra + Plots): L-BFGS two-loop recursion, Armijo
  line search along a **feasible** step length (max step to stay inside bounds),
  then optional projection of the iterate.

  Run from repo root:
    julia --project=. demo_example/lbfgsb_toy_example.jl

  Writes PNGs under demo_example/:
    - lbfgsb_toy_spectra_manual.png
    - lbfgsb_toy_spectra_optim.png
    - lbfgsb_toy_rmse_compare.png
=#

using LinearAlgebra
using Printf
using Random
using Statistics
using Plots

const HAVE_OPTIM = try
    @eval using Optim
    true
catch
    false
end

Random.seed!(42)

const N_WL = 80
wavelengths() = range(650.0, 850.0; length=N_WL)

# --- Forward model: nonlinear spectral toy y = f(x), x ∈ ℝ³ ---
function forward_spectral(x::AbstractVector{<:Real})
    a, μ, w = Float64(x[1]), Float64(x[2]), Float64(x[3])
    λ = wavelengths()
    return @. a * exp(-0.5 * ((λ - μ) / w)^2)
end

spectrum_rmse(x::AbstractVector{<:Real}, y_obs::AbstractVector{<:Real}) =
    sqrt(mean((forward_spectral(x) .- y_obs) .^ 2))

function misfit(x::Vector{Float64}, y_obs::Vector{Float64})::Float64
    y = forward_spectral(x)
    r = y .- y_obs
    return dot(r, r)
end

function misfit_grad(x::Vector{Float64}, y_obs::Vector{Float64})::Vector{Float64}
    y = forward_spectral(x)
    r = y .- y_obs
    a, μ, w = x[1], x[2], x[3]
    λ = wavelengths()
    z = @. (λ - μ) / w
    ex = @. exp(-0.5 * z^2)
    J_a = ex
    J_μ = @. a * ex * z / w
    J_w = @. a * ex * z^2 / w
    J = hcat(J_a, J_μ, J_w)
    return 2 * vec(sum(J .* r, dims=1))
end

"""L-BFGS search direction -H∇f (two-loop), limited memory `m`."""
function lbfgs_direction!(
    g::Vector{Float64},
    s_hist::Vector{Vector{Float64}},
    y_hist::Vector{Vector{Float64}},
    rho::Vector{Float64},
    m::Int,
)::Vector{Float64}
    q = copy(g)
    α = zeros(length(s_hist))
    for i in length(s_hist):-1:1
        α[i] = rho[i] * dot(s_hist[i], q)
        q .-= α[i] .* y_hist[i]
    end
    if isempty(s_hist)
        return -q
    end
    γ = dot(s_hist[end], y_hist[end]) / dot(y_hist[end], y_hist[end])
    γ = max(γ, 1e-10)
    r = γ .* q
    for i in eachindex(s_hist)
        β = rho[i] * dot(y_hist[i], r)
        r .+= s_hist[i] .* (α[i] - β)
    end
    return -r
end

"""Largest α>0 such that lb ≤ x + α p ≤ ub (componentwise)."""
function max_step_to_bounds(x::Vector{Float64}, p::Vector{Float64}, lb::Vector{Float64}, ub::Vector{Float64})
    α_max = Inf
    @inbounds for i in eachindex(x)
        pi = p[i]
        if pi > 1e-30
            α_max = min(α_max, (ub[i] - x[i]) / pi)
        elseif pi < -1e-30
            α_max = min(α_max, (lb[i] - x[i]) / pi)
        end
    end
    return isfinite(α_max) ? max(α_max, 0.0) : 0.0
end

"""
  Minimize `f` with gradient `g` and box `[lb, ub]`.
  L-BFGS direction + Armijo backtracking on `φ(α) = f(P(x + α p))` with `P` = clamp.
"""
function minimize_box_lbfgs(
    f,
    g,
    x0::Vector{Float64},
    lb::Vector{Float64},
    ub::Vector{Float64};
    m::Int=10,
    maxiter::Int=200,
    tol_g::Float64=1e-8,
    c1::Float64=1e-4,
    max_backtrack::Int=40,
    obs_for_trace::Union{Nothing,Vector{Float64}}=nothing,
    trace_rmse::Bool=false,
    snapshot_every::Int=0,
)
    if (trace_rmse || snapshot_every > 0) && obs_for_trace === nothing
        error("obs_for_trace is required when trace_rmse=true or snapshot_every > 0")
    end
    do_trace = trace_rmse || (snapshot_every > 0 && obs_for_trace !== nothing)
    rmse_trace = do_trace ? Float64[] : nothing
    snap_iters = Int[]
    snap_y = Vector{Vector{Float64}}()

    x = clamp.(copy(x0), lb, ub)
    s_hist = Vector{Vector{Float64}}()
    y_hist = Vector{Vector{Float64}}()
    rho = Float64[]

    fx = f(x)
    gx = g(x)

    if do_trace && obs_for_trace !== nothing
        push!(rmse_trace, spectrum_rmse(x, obs_for_trace))
        if snapshot_every > 0
            push!(snap_iters, 0)
            push!(snap_y, copy(forward_spectral(x)))
        end
    end

    for k in 1:maxiter
        gnorm = norm(gx, Inf)
        @printf("  L-BFGS-B iter %3d  f = %.6e  ‖g‖∞ = %.3e  x = [%.4f, %.2f, %.2f]\n",
                k, fx, gnorm, x[1], x[2], x[3])
        if gnorm < tol_g
            return (
                x=x,
                fx=fx,
                iter=k,
                converged=true,
                rmse_trace=rmse_trace,
                snap_iters=snap_iters,
                snap_y=snap_y,
            )
        end

        p = lbfgs_direction!(gx, s_hist, y_hist, rho, m)
        # Descent unless badly scaled; flip if necessary
        if dot(gx, p) > 0
            p .= -gx
        end

        α_cap = max_step_to_bounds(x, p, lb, ub)
        α_cap <= 0 && error("zero feasible step; check bounds vs. x0")

        α = min(1.0, 0.99 * α_cap)
        suff = false
        for _ in 1:max_backtrack
            x_try = clamp.(x .+ α .* p, lb, ub)
            f_try = f(x_try)
            if f_try ≤ fx + c1 * dot(gx, x_try .- x)  # Armijo on projected update (common surrogate)
                suff = true
                break
            end
            α *= 0.5
            if α < 1e-16 * α_cap
                break
            end
        end
        suff || @warn "line search failed at iter $k; stopping"
        suff || return (
            x=x,
            fx=fx,
            iter=k,
            converged=false,
            rmse_trace=rmse_trace,
            snap_iters=snap_iters,
            snap_y=snap_y,
        )

        x_new = clamp.(x .+ α .* p, lb, ub)
        f_new = f(x_new)
        g_new = g(x_new)

        s = x_new .- x
        yv = g_new .- gx
        ys = dot(yv, s)
        if ys > 1e-12
            while length(s_hist) ≥ m
                popfirst!(s_hist)
                popfirst!(y_hist)
                popfirst!(rho)
            end
            push!(s_hist, s)
            push!(y_hist, yv)
            push!(rho, 1.0 / ys)
        end

        x .= x_new
        fx = f_new
        gx = g_new

        if do_trace && obs_for_trace !== nothing
            push!(rmse_trace, spectrum_rmse(x, obs_for_trace))
            if snapshot_every > 0 && k % snapshot_every == 0
                push!(snap_iters, k)
                push!(snap_y, copy(forward_spectral(x)))
            end
        end
    end
    return (
        x=x,
        fx=fx,
        iter=maxiter,
        converged=false,
        rmse_trace=rmse_trace,
        snap_iters=snap_iters,
        snap_y=snap_y,
    )
end

function minimize_optim_lbfgsb(
    f,
    g,
    x0::Vector{Float64},
    lb::Vector{Float64},
    ub::Vector{Float64};
    m::Int=10,
    maxiter::Int=200,
    tol_g::Float64=1e-8,
    f_tol::Float64=0.0,
    x_tol::Float64=0.0,
    obs_for_trace::Union{Nothing,Vector{Float64}}=nothing,
)
    HAVE_OPTIM || error(
        "Optim.jl is required for comparison. Install with:\n" *
        "  julia --project=. -e 'using Pkg; Pkg.add(\"Optim\")'",
    )

    x_trace = Vector{Vector{Float64}}()
    rmse_trace = Float64[]
    iter_count = Ref(0)

    x0c = clamp.(copy(x0), lb, ub)
    push!(x_trace, copy(x0c))
    if obs_for_trace !== nothing
        push!(rmse_trace, spectrum_rmse(x0c, obs_for_trace))
    end

    function g!(G, x)
        G .= g(collect(Float64.(x)))
        return nothing
    end
    f_obj(x) = f(collect(Float64.(x)))

    od = Optim.OnceDifferentiable(f_obj, g!, x0c)
    callback = function (state)
        iter_count[] += 1
        xv = nothing
        if hasproperty(state, :x)
            xv = collect(Float64.(getproperty(state, :x)))
        elseif hasproperty(state, :metadata)
            md = getproperty(state, :metadata)
            if md isa AbstractDict
                if haskey(md, "x")
                    xv = collect(Float64.(md["x"]))
                elseif haskey(md, :x)
                    xv = collect(Float64.(md[:x]))
                end
            end
        end
        xv === nothing && return false

        push!(x_trace, copy(xv))
        if obs_for_trace !== nothing
            push!(rmse_trace, spectrum_rmse(xv, obs_for_trace))
        end
        fval = hasproperty(state, :value) ? getproperty(state, :value) : f_obj(xv)
        @printf(
            "  Optim L-BFGS-B iter %3d  f = %.6e  x = [%.4f, %.2f, %.2f]\n",
            iter_count[],
            fval,
            xv[1],
            xv[2],
            xv[3],
        )
        return false
    end

    options = Optim.Options(
        iterations=maxiter,
        g_tol=tol_g,
        f_reltol=f_tol,
        x_abstol=x_tol,
        show_trace=false,
        callback=callback,
    )
    res = Optim.optimize(od, lb, ub, x0c, Optim.Fminbox(Optim.LBFGS(m=m)), options)
    xf = collect(Float64.(Optim.minimizer(res)))

    return (
        x=xf,
        fx=Optim.minimum(res),
        iter=iter_count[],
        converged=Optim.converged(res),
        rmse_trace=rmse_trace,
        x_trace=x_trace,
    )
end

function plot_spectrum_progress(
    λ::AbstractVector{Float64},
    y_truth::Vector{Float64},
    y_obs::Vector{Float64},
    snap_iters::Vector{Int},
    snap_y::Vector{Vector{Float64}};
    outpath::String,
)
    n = length(snap_iters)
    n == 0 && return nothing
    ps = []
    for i in 1:n
        lab = snap_iters[i] == 0 ? "initial" : "iter $(snap_iters[i])"
        p = plot(
            λ, y_truth;
            label="truth (noise-free)",
            color=:black,
            lw=2,
            title=lab,
            xlabel="λ (nm)",
            ylabel="radiance (arb.)",
            legend=:outerright,
            legendfontsize=7,
        )
        plot!(p, λ, y_obs; label="observed", color=:gray, ls=:dash, lw=1.5)
        plot!(p, λ, snap_y[i]; label="reconstructed", color=:royalblue, lw=1.5)
        push!(ps, p)
    end
    h = plot(ps...; layout=(n, 1), size=(720, min(220 * n, 2400)), link=:x)
    savefig(h, outpath)
    println("Saved spectrum figure: ", outpath)
    return nothing
end

function main()
    x_true = [1.2, 740.0, 25.0]
    y_obs = forward_spectral(x_true) .+ 0.02 .* randn(length(forward_spectral(x_true)))
    y_truth = forward_spectral(x_true)
    λ = collect(wavelengths())

    x0 = [0.5, 720.0, 18.0]
    lb = [0.1, 700.0, 8.0]
    ub = [3.0, 780.0, 40.0]

    println("Truth: x* = ", x_true)
    println("Start: x0 = ", x0)
    println("Bounds: lb = ", lb, ", ub = ", ub)
    println("Running manual box L-BFGS ...\n")

    f(x) = misfit(x, y_obs)
    g(x) = misfit_grad(x, y_obs)
    m = 10
    maxiter = 150
    tol_g = 1e-8

    res_manual = minimize_box_lbfgs(
        f, g, x0, lb, ub;
        m=m,
        maxiter=maxiter,
        tol_g=tol_g,
        obs_for_trace=y_obs,
        trace_rmse=true,
        snapshot_every=5,
    )

    println("\nRunning Optim.jl Fminbox(LBFGS) ...\n")
    res_optim = minimize_optim_lbfgsb(
        f, g, x0, lb, ub;
        m=m,
        maxiter=maxiter,
        tol_g=tol_g,
        obs_for_trace=y_obs,
    )

    println("\n--- Manual Result ---")
    println("converged: ", res_manual.converged)
    @printf("final f = %.6e\n", res_manual.fx)
    println("x_hat     = ", res_manual.x)
    println("x_true    = ", x_true)
    println("‖x_hat - x*‖ = ", norm(res_manual.x .- x_true))

    println("\n--- Optim.jl Result ---")
    println("converged: ", res_optim.converged)
    @printf("final f = %.6e\n", res_optim.fx)
    println("x_hat     = ", res_optim.x)
    println("x_true    = ", x_true)
    println("‖x_hat - x*‖ = ", norm(res_optim.x .- x_true))

    # Ensure final reconstructed spectrum appears in the multi-panel plot
    si_manual = copy(res_manual.snap_iters)
    sy_manual = copy(res_manual.snap_y)
    if res_manual.iter > 0 && !(res_manual.iter in si_manual)
        push!(si_manual, res_manual.iter)
        push!(sy_manual, forward_spectral(res_manual.x))
    end
    si_optim = Int[]
    sy_optim = Vector{Vector{Float64}}()
    for k in 0:5:res_optim.iter
        idx = min(k + 1, length(res_optim.x_trace))
        push!(si_optim, k)
        push!(sy_optim, forward_spectral(res_optim.x_trace[idx]))
    end
    if isempty(si_optim) || si_optim[end] != res_optim.iter
        push!(si_optim, res_optim.iter)
        push!(sy_optim, forward_spectral(res_optim.x))
    end

    out_dir = @__DIR__
    plot_spectrum_progress(
        λ, y_truth, y_obs, si_manual, sy_manual;
        outpath=joinpath(out_dir, "lbfgsb_toy_spectra_manual.png"),
    )
    plot_spectrum_progress(
        λ, y_truth, y_obs, si_optim, sy_optim;
        outpath=joinpath(out_dir, "lbfgsb_toy_spectra_optim.png"),
    )

    if res_manual.rmse_trace !== nothing && !isempty(res_manual.rmse_trace) &&
       res_optim.rmse_trace !== nothing && !isempty(res_optim.rmse_trace)
        iters_manual = 0:(length(res_manual.rmse_trace) - 1)
        iters_optim = 0:(length(res_optim.rmse_trace) - 1)
        p_rmse = plot(
            collect(iters_manual),
            res_manual.rmse_trace;
            marker=:circle,
            ms=3,
            lw=2,
            color=:darkgreen,
            xlabel="iteration",
            ylabel="RMSE vs observations",
            title="RMSE comparison: manual vs Optim.jl L-BFGS-B",
            label="manual",
        )
        plot!(
            p_rmse,
            collect(iters_optim),
            res_optim.rmse_trace;
            marker=:diamond,
            ms=3,
            lw=2,
            color=:royalblue,
            label="Optim.jl",
        )
        rmse_path = joinpath(out_dir, "lbfgsb_toy_rmse_compare.png")
        savefig(p_rmse, rmse_path)
        println("Saved RMSE figure: ", rmse_path)
    end

    return nothing
end

main()
