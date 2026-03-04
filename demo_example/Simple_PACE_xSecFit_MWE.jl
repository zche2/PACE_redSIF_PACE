#!/usr/bin/env julia

using PACE_SIF

include(joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE_Functions.jl"))
using .SimplePACEXSecFitMWEFunctions

"""
    main()

Entry point for the MWE setup workflow.
Loads config, prepares LUT/kernel/SIF inputs, and prints diagnostics.
No retrieval or inversion is run here.
"""
function main()
    config_path = get(
        ENV,
        "PACE_MWE_CONFIG",
        joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE.toml"),
    )

    ctx = prepare_mwe_inputs(config_path)

    println("Simple PACE xSec MWE (no fit)")
    println("config: ", config_path)
    println("xsec run: ", ctx.paths.xsec_run)
    println("xsec dir: ", ctx.paths.xsec_dir)
    println("o2 lut: ", ctx.paths.o2_path)
    println("h2o lut: ", ctx.paths.h2o_path)
    println("sif file: ", ctx.paths.sif_path)
    println("pace rsr file: ", ctx.paths.pace_rsr_path)
    println("stored kernel file: ", ctx.paths.kernel_path)
    println("axis unit: ", ctx.axis_unit)
    println("high-res Δλ [nm]: ", ctx.delta_lambda_nm)
    println("high-res λ range [nm]: ", (minimum(ctx.λ_hres), maximum(ctx.λ_hres)))
    println("high-res λ grid size: ", length(ctx.λ_hres))
    println("fit window bands: ", length(ctx.λ))
    println("low-res λ source: ", ctx.lambda_source)
    println("sif basis size high-res (λ_hres × EV): ", size(ctx.sif_basis_hres))
    println("sif basis size low-res (λ × EV): ", size(ctx.sif_basis))
    println()

    println("Kernel generation")
    println("  sampled RSR wavelengths: ", ctx.generation_info.n_wavelength_samples)
    println("  selected bands: ", ctx.generation_info.n_bands)
    println("  selected λ range [nm]: ", ctx.generation_info.λ_range_nm)
    println("  row-sum min/max: ", (ctx.regenerated_kernel_quality.row_sum_min, ctx.regenerated_kernel_quality.row_sum_max))
    println("  negative weight fraction: ", ctx.regenerated_kernel_quality.negative_weight_fraction)
    println("  center offset mean [nm]: ", ctx.regenerated_kernel_quality.center_offset_mean_nm)
    println("  center offset abs max [nm]: ", ctx.regenerated_kernel_quality.center_offset_abs_max_nm)
    println("  peak offset abs max [nm]: ", ctx.regenerated_kernel_quality.peak_offset_abs_max_nm)

    if !isnothing(ctx.kernel_comparison)
        cmp = ctx.kernel_comparison
        println()
        println("Kernel comparison (stored vs regenerated)")
        println("  common bands: ", cmp.n_common)
        println("  max band mismatch [nm]: ", cmp.band_max_delta_nm)
        println("  L1 mean/max: ", (cmp.l1_mean, cmp.l1_max))
        println("  L2 mean/max: ", (cmp.l2_mean, cmp.l2_max))
        println("  max|Δ| mean/max: ", (cmp.max_abs_mean, cmp.max_abs_max))
    else
        println()
        println("Kernel comparison: skipped (stored kernel missing or disabled)")
    end

    # ------------------------------------------------------------------
    # Your own model/inversion starts here:
    #
    # λ        = ctx.λ
    # λc       = ctx.λc
    # K        = ctx.kernel                     # regenerated kernel
    # o2_sitp  = ctx.o2_sitp
    # h2o_sitp = ctx.h2o_sitp
    # SIF_U_hres = ctx.sif_basis_hres            # for high-res forward simulation
    # SIF_U      = ctx.sif_basis                 # low-res (convolved) convenience matrix
    #
    # Define forward model:
    #   y = f(x, λ, K, o2_sitp, h2o_sitp, SIF_U, ...)
    #
    # Then run your own optimizer/inversion.
    # ------------------------------------------------------------------
end

main()
