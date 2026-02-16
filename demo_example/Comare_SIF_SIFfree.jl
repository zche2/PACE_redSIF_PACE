#!/usr/bin/env julia
# Compare_SIF_SIFfree.jl

using TOML
using LinearAlgebra
using Statistics
using ForwardDiff
using SparseArrays
using Plots
using JLD2
using CSV
using DataFrames

include(joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE_Functions.jl"))
using .SimplePACEXSecFitMWEFunctions

# Import needed functions explicitly
import .SimplePACEXSecFitMWEFunctions: prepare_mwe_inputs, load_solar_spectrum_on_grid, load_pace_spectrum_on_grid

include(joinpath(@__DIR__, "toy_forward_model_flexibleSIF.jl"))

# Wrap Fit_toy_forward_model.jl in its own module
module FitHelpers
    include(joinpath(@__DIR__, "Fit_toy_forward_model.jl"))
end

# Now use only the specific function I need
const lm_one_step = FitHelpers.lm_one_step
const _spdiag_invvar = FitHelpers._spdiag_invvar

"""
Run LM retrieval and return results
"""
function run_lm_retrieval(
    fm,
    x0,
    y_obs,
    layout;
    x_a,
    S_a_inv,
    S_e_inv,
    lambda0=1.0,
    max_steps=20,
    jacobian_eval=nothing,
    x_scale=nothing,
    lower_bounds=nothing,
    upper_bounds=nothing,
    verbose=true,
)
    x_curr = copy(x_a)
    λ = lambda0
    
    history = []
    
    for istep in 1:max_steps
        step = lm_one_step(
            fm,
            x_curr,
            y_obs;
            x_a=x_a,
            S_e_inv=S_e_inv,
            S_a_inv=S_a_inv,
            lambda=λ,
            lambda_up=2.0,
            lambda_down=0.7,
            lambda_min=1e-8,
            lambda_max=1e8,
            max_inner=8,
            jacobian_eval=jacobian_eval,
            x_scale=x_scale,
            lower_bounds=lower_bounds,
            upper_bounds=upper_bounds,
        )
        
        λ = step.lambda_next
        
        if !step.accepted
            if verbose
                println("  Step $istep: No accepted update")
            end
            break
        end
        
        x_curr = step.x_next
        
        push!(history, (
            iter = istep,
            rmse = step.rmse_next,
            chi2 = 2.0 * step.cost_next,
            lambda = λ,
            dx_norm = norm(step.dx),
        ))
        
        if verbose
            println("  Step $istep: RMSE=$(step.rmse_next), χ²=$(2.0*step.cost_next), λ=$λ")
        end
        
        # Convergence check
        if step.rmse_next < 1e-6 || norm(step.dx) < 1e-6
            if verbose
                println("  Converged!")
            end
            break
        end
    end
    
    y_final = fm(x_curr)
    residual = y_obs .- y_final
    
    return (
        x_final = x_curr,
        y_final = y_final,
        residual = residual,
        rmse = sqrt(mean(residual.^2)),
        history = history,
    )
end

function main()
    println("Current directory: $(abspath(@__DIR__))")
    println("="^70)
    println("Comparing SIF vs SIF-free retrievals")
    println("="^70)
    
    # Load configuration
    config_path = joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE_zcheVer.toml")
    cfg = TOML.parsefile(config_path)
    
    fit_cfg = get(cfg, "fit", Dict{String, Any}())
    n_legendre = Int(get(fit_cfg, "n_legendre", 2))
    meas_sigma = Float64(get(fit_cfg, "meas_sigma", 0.01))
    
    # Get output directory from config or use default
    out_dir = String(get(fit_cfg, "out_dir", joinpath(@__DIR__, "toy_fit")))
    mkpath(out_dir)
    println("\nOutput directory: $out_dir")
    
    # Prepare inputs
    println("\nLoading data and models...")
    ctx = prepare_mwe_inputs(config_path)
    
    # Load solar spectrum
    data_cfg = get(cfg, "data", Dict{String, Any}())
    solar_file = get(data_cfg, "solar_file", "solar_merged_20240731_600_33300_100.out")
    solar_path = joinpath(ctx.paths.base_dir, solar_file)
    solar_hres, _ = load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
    
    # Load observation
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())
    pace_file = get(pace_cfg, "pace_file", "sample_granule_20240830T131442_new_chl.nc")
    pace_path = joinpath(ctx.paths.base_dir, pace_file)
    pixel_idx = Int(get(pace_cfg, "pixel_index", 600))
    scan_idx = Int(get(pace_cfg, "scan_index", 800))
    
    y_obs, obs_info = load_pace_spectrum_on_grid(
        pace_path,
        ctx.λ;
        pixel_idx=pixel_idx,
        scan_idx=scan_idx,
        wavelength_var=String(get(pace_cfg, "wavelength_var", "red_wavelength")),
        spectrum_var=String(get(pace_cfg, "spectrum_var", "radiance_red")),
    )
    
    println("  Observation: pixel=$pixel_idx, scan=$scan_idx")
    println("  n_wavelengths: $(length(y_obs))")
    
    # =========================================================================
    # Retrieval 1: WITH SIF (n_ev > 0)
    # =========================================================================
    
    println("\n" * "="^70)
    println("RETRIEVAL 1: WITH SIF (n_ev = $(ctx.sif_nev))")
    println("="^70)
    
    # Build forward model with SIF
    fm_sif = make_forward_model_simple(
        ctx,
        solar_hres;
        n_legendre=n_legendre,
        preallocate_float64=true,
    )
    
    layout_sif  = state_layout_simple(ctx; n_legendre=n_legendre)
    state_names = state_names_simple(ctx; n_legendre=n_legendre)
    x0_sif = initial_state_simple(ctx; n_legendre=n_legendre)
    
    println("  State vector size: $(layout_sif.n_state)")
    println("  Including: $(layout_sif.n_ev) SIF components")
    
    # Set up priors
    x_a_sif = copy(x0_sif)
    prior_sigma_sif = fill(1e30, length(x0_sif))
    
    # Priors on p/T
    prior_sigma_sif[layout_sif.idx_p_o2_hpa] = 200.0
    prior_sigma_sif[layout_sif.idx_p_h2o_hpa] = 200.0
    prior_sigma_sif[layout_sif.idx_t_o2_k] = 20.0
    prior_sigma_sif[layout_sif.idx_t_h2o_k] = 20.0
    
    # Priors on VCD
    prior_sigma_sif[layout_sif.idx_vcd_o2_intercept] = 1e23
    prior_sigma_sif[layout_sif.idx_vcd_h2o_intercept] = 3e22
    
    if layout_sif.n_ev > 0
        prior_sigma_sif[layout_sif.idx_vcd_o2_sif] = 1e23
        prior_sigma_sif[layout_sif.idx_vcd_h2o_sif] = 3e22
        prior_sigma_sif[layout_sif.idx_sif] .= 1e12
    end
    
    S_a_inv_sif = _spdiag_invvar(prior_sigma_sif)
    S_e_inv = spdiagm(0 => fill(1.0 / (meas_sigma^2), length(y_obs)))
    
    # Run retrieval with SIF
    println("\nRunning retrieval WITH SIF...")
    result_sif = run_lm_retrieval(
        fm_sif,
        x0_sif,
        y_obs,
        layout_sif;
        x_a=x_a_sif,
        S_a_inv=S_a_inv_sif,
        S_e_inv=S_e_inv,
        lambda0=1.0,
        max_steps=20,
        verbose=true,
    )
    
    println("\nRETRIEVAL 1 RESULTS:")
    println("  Final RMSE: $(result_sif.rmse)")
    println("  Iterations: $(length(result_sif.history))")
    
    # =========================================================================
    # Retrieval 2: WITHOUT SIF (n_ev = 0)
    # =========================================================================
    
    println("\n" * "="^70)
    println("RETRIEVAL 2: WITHOUT SIF (n_ev = 0)")
    println("="^70)
    ctx_nosif = merge(ctx, Dict(:sif_nev => 0));

    # Build forward model without SIF
    fm_nosif = make_forward_model_simple(
        ctx_nosif,
        solar_hres;
        n_legendre=n_legendre,
        preallocate_float64=true,
    )
    
    layout_nosif = state_layout_simple(ctx_nosif; n_legendre=n_legendre)
    state_names_nosif = state_names_simple(ctx_nosif; n_legendre=n_legendre)
    x0_nosif = initial_state_simple(ctx_nosif; n_legendre=n_legendre)
    
    println("  State vector size: $(layout_nosif.n_state)")
    println("  No SIF components")
    
    # Set up priors (no SIF)
    x_a_nosif = copy(x0_nosif)
    prior_sigma_nosif = fill(1e30, length(x0_nosif))
    
    prior_sigma_nosif[layout_nosif.idx_p_o2_hpa] = 200.0
    prior_sigma_nosif[layout_nosif.idx_p_h2o_hpa] = 200.0
    prior_sigma_nosif[layout_nosif.idx_t_o2_k] = 20.0
    prior_sigma_nosif[layout_nosif.idx_t_h2o_k] = 20.0
    prior_sigma_nosif[layout_nosif.idx_vcd_o2_intercept] = 1e23
    prior_sigma_nosif[layout_nosif.idx_vcd_h2o_intercept] = 3e22
    
    S_a_inv_nosif = _spdiag_invvar(prior_sigma_nosif)
    
    # Run retrieval without SIF
    println("\nRunning retrieval WITHOUT SIF...")
    result_nosif = run_lm_retrieval(
        fm_nosif,
        x0_nosif,
        y_obs,
        layout_nosif;
        x_a=x_a_nosif,
        S_a_inv=S_a_inv_nosif,
        S_e_inv=S_e_inv,
        lambda0=1.0,
        max_steps=20,
        verbose=true,
    )
    
    println("\nRETRIEVAL 2 RESULTS:")
    println("  Final RMSE: $(result_nosif.rmse)")
    println("  Iterations: $(length(result_nosif.history))")
    
    # =========================================================================
    # Compare results
    # =========================================================================
    
    println("\n" * "="^70)
    println("COMPARISON")
    println("="^70)
    println("  RMSE with SIF:    $(round(result_sif.rmse, digits=6))")
    println("  RMSE without SIF: $(round(result_nosif.rmse, digits=6))")
    println("  RMSE improvement: $(round((result_nosif.rmse - result_sif.rmse)/result_nosif.rmse * 100, digits=2))%")
    
    # =========================================================================
    # Save results to CSV
    # =========================================================================
    
    println("\n" * "="^70)
    println("SAVING RESULTS")
    println("="^70)
    
    # Save WITH SIF results
    df_sif = DataFrame(
        wavelength = ctx.λ,
        observation = y_obs,
        fitted = result_sif.y_final,
        residual = result_sif.residual,
    )
    csv_sif = joinpath(out_dir, "retrieval_with_SIF.csv")
    CSV.write(csv_sif, df_sif)
    println("  Saved WITH SIF: $csv_sif")
    
    # Save WITHOUT SIF results
    df_nosif = DataFrame(
        wavelength = ctx.λ,
        observation = y_obs,
        fitted = result_nosif.y_final,
        residual = result_nosif.residual,
    )
    csv_nosif = joinpath(out_dir, "retrieval_without_SIF.csv")
    CSV.write(csv_nosif, df_nosif)
    println("  Saved WITHOUT SIF: $csv_nosif")
    
    # Save comparison statistics
    df_comparison = DataFrame(
        retrieval = ["With SIF", "Without SIF"],
        n_state = [layout_sif.n_state, layout_nosif.n_state],
        n_sif_components = [layout_sif.n_ev, layout_nosif.n_ev],
        final_rmse = [result_sif.rmse, result_nosif.rmse],
        n_iterations = [length(result_sif.history), length(result_nosif.history)],
    )
    csv_comparison = joinpath(out_dir, "retrieval_comparison.csv")
    CSV.write(csv_comparison, df_comparison)
    println("  Saved comparison: $csv_comparison")
    
    # =========================================================================
    # Plot comparison
    # =========================================================================
    
    println("\n" * "="^70)
    println("PLOTTING COMPARISON")
    println("="^70)
    
    # Plot 1: Spectra
    p1 = plot(ctx.λ, y_obs;
              label="Observation",
              xlabel="Wavelength [nm]",
              ylabel="Radiance",
              title="Fitted Spectra",
              lw=2,
              color=:black,
              legend=:topright)
    plot!(p1, ctx.λ, result_sif.y_final;
          label="With SIF (RMSE=$(round(result_sif.rmse, digits=4)))",
          lw=2, ls=:dash, color=:blue)
    plot!(p1, ctx.λ, result_nosif.y_final;
          label="Without SIF (RMSE=$(round(result_nosif.rmse, digits=4)))",
          lw=2, ls=:dot, color=:red)
    
    # Plot 2: Residuals
    p2 = plot(ctx.λ, result_sif.residual;
              label="With SIF",
              xlabel="Wavelength [nm]",
              ylabel="Residual",
              title="Fit Residuals",
              lw=2,
              color=:blue,
              legend=:topright)
    plot!(p2, ctx.λ, result_nosif.residual;
          label="Without SIF",
          lw=2, color=:red)
    hline!(p2, [0.0]; color=:black, ls=:dash, alpha=0.5, label="")
    
    # Plot 3: Absolute residuals
    p3 = plot(ctx.λ, abs.(result_sif.residual);
              label="With SIF",
              xlabel="Wavelength [nm]",
              ylabel="|Residual|",
              title="Absolute Residuals",
              lw=2,
              color=:blue,
              yscale=:log10,
              legend=:topright)
    plot!(p3, ctx.λ, abs.(result_nosif.residual);
          label="Without SIF",
          lw=2, color=:red)
    
    # Plot 4: Histogram of residuals
    p4 = histogram(result_sif.residual;
                   label="With SIF",
                   xlabel="Residual",
                   ylabel="Frequency",
                   title="Residual Distribution",
                   alpha=0.6,
                   color=:blue,
                   bins=30)
    histogram!(p4, result_nosif.residual;
              label="Without SIF",
              alpha=0.6,
              color=:red,
              bins=30)
    
    # Combine plots
    p_all = plot(p1, p2, p3, p4;
                layout=(2, 2),
                size=(1400, 1000),
                dpi=150)
    
    plot_file = joinpath(out_dir, "SIF_vs_SIFfree_comparison.png")
    savefig(p_all, plot_file)
    println("  Saved plot: $plot_file")
    
    # =========================================================================
    # Summary
    # =========================================================================
    
    println("\n" * "="^70)
    println("SUMMARY")
    println("="^70)
    println("  Observation: pixel=$pixel_idx, scan=$scan_idx")
    println("  Wavelength range: $(minimum(ctx.λ)) - $(maximum(ctx.λ)) nm")
    println("\n  WITH SIF:")
    println("    State size: $(layout_sif.n_state)")
    println("    SIF components: $(layout_sif.n_ev)")
    println("    Final RMSE: $(result_sif.rmse)")
    println("    Iterations: $(length(result_sif.history))")
    println("    State vector and names:")
    for (i, name) in enumerate(state_names)
        println("      [$i] $name: $(result_sif.x_final[i])")
    end
    println("\n  WITHOUT SIF:")
    println("    State size: $(layout_nosif.n_state)")
    println("    Final RMSE: $(result_nosif.rmse)")
    println("    Iterations: $(length(result_nosif.history))")
    println("    State vector and names:")
    for (i, name) in enumerate(state_names_nosif)
        println("      [$i] $name: $(result_nosif.x_final[i])")
    end
    println("\n  IMPROVEMENT:")
    println("    ΔRMSE: $(result_nosif.rmse - result_sif.rmse)")
    println("    % improvement: $(round((result_nosif.rmse - result_sif.rmse)/result_nosif.rmse * 100, digits=2))%")
    println("="^70)
    
    return (sif=result_sif, nosif=result_nosif)
end

# Run comparison
if abspath(PROGRAM_FILE) == @__FILE__
    results = main()
end