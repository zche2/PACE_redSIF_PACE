#!/usr/bin/env julia
# Synthetic data experiment: vertical profiles (72-layer MERRA2) + SIF shapes
# -> transmittance at high spectral resolution -> TOA radiance at high res
# -> convolve to OCI resolution -> run retrieval. Compare retrieved state to truth.

using TOML
using NCDatasets
using LinearAlgebra
using Statistics
using Interpolations

const PSEUDO_DIR = @__DIR__
const DEMO_DIR = joinpath(PSEUDO_DIR, "..")
include(joinpath(DEMO_DIR, "Fit_toy_forward_model.jl"))

const MWEF = SimplePACEXSecFitMWEFunctions

"""
    compute_trans_hres_from_profile(trans_nc_path, profile_index, ctx, amf; o2_sitp, h2o_sitp)

Load one profile from the MERRA2 transmittance NetCDF (72 layers) and compute
one-way transmittance at high spectral resolution using LUTs.
Returns trans_hres (vector on ctx.λ_hres).
"""
function compute_trans_hres_from_profile(
    trans_nc_path::AbstractString,
    profile_index::Int,
    ctx,
    amf::Float64;
    o2_sitp,
    h2o_sitp,
)
    ds = Dataset(trans_nc_path)
    n_profiles = ds.dim["profile"]
    n_layers = ds.dim["layer"]
    1 <= profile_index <= n_profiles || error("profile_index=$profile_index out of range 1:$n_profiles")

    vcd_dry = ds["vcd_dry"][profile_index, :]   # (layer,)
    vcd_h2o = ds["vcd_h2o"][profile_index, :]
    temp = ds["temperature"][profile_index, :]
    ps_hpa = Float64(ds["pressure"][profile_index])  # surface pressure in hPa
    ak = ds.attrib["ak"]
    bk = ds.attrib["bk"]
    close(ds)

    # Half-level pressures (Pa -> hPa)
    p_half = (ak .+ bk .* (ps_hpa * 100)) ./ 100
    p_full = (p_half[1:end-1] .+ p_half[2:end]) ./ 2
    length(p_full) == n_layers || error("p_full length mismatch")

    spectral_axis = ctx.spectral_axis
    n_λ = length(spectral_axis)
    vmr_o2 = 0.21

    τ = zeros(Float64, n_λ)
    for l in 1:n_layers
        xsec_o2 = vec(o2_sitp(spectral_axis, p_full[l], temp[l]))
        xsec_h2o = vec(h2o_sitp(spectral_axis, p_full[l], temp[l]))
        τ .+= xsec_o2 .* (vcd_dry[l] * vmr_o2) .+ xsec_h2o .* vcd_h2o[l]
    end
    trans_hres = exp.(-amf .* τ)
    return trans_hres
end

"""
    generate_synthetic_spectrum(ctx, solar_hres, trans_hres; sif_coeff, leg_coeff=nothing, n_legendre=3)

Build TOA radiance at high resolution and convolve to OCI bands.
TOA_hres = solar_hres .* trans_2way + trans_1way .* (sif_basis_hres * sif_coeff).
Returns y_obs on ctx.λ. If leg_coeff provided (length n_legendre+1), multiply convolved radiance by Legendre baseline.
"""
function generate_synthetic_spectrum(
    ctx,
    solar_hres::AbstractVector{<:Real},
    trans_hres::AbstractVector{<:Real};
    sif_coeff::AbstractVector{<:Real},
    leg_coeff::Union{Nothing, AbstractVector{<:Real}}=nothing,
    n_legendre::Int=3,
)
    trans_2way = trans_hres .* trans_hres
    trans_1way = trans_hres
    sif_hres = ctx.sif_basis_hres * sif_coeff
    toa_hres = solar_hres .* trans_2way .+ trans_1way .* sif_hres
    K = ctx.kernel_rsr_out
    y_lres = K * toa_hres
    if leg_coeff !== nothing
        z = _normalized_grid(ctx.λ)
        leg_basis = _legendre_design_matrix(z, n_legendre)
        y_lres .*= leg_basis * leg_coeff
    end
    return y_lres
end

function main_synthetic()
    pseudo_cfg_path = joinpath(PSEUDO_DIR, "pseudo_measurement_config.toml")
    pseudo_cfg = TOML.parsefile(pseudo_cfg_path)
    pm = get(pseudo_cfg, "pseudo_measurement", Dict{String, Any}())
    main_config = normpath(joinpath(PSEUDO_DIR, get(pm, "main_config_path", "../Simple_PACE_xSecFit_MWE_zcheVer.toml")))
    trans_nc = get(pm, "trans_nc", "")
    isempty(trans_nc) && error("Set pseudo_measurement.trans_nc in config")
    trans_nc = isabspath(trans_nc) ? trans_nc : normpath(joinpath(DEMO_DIR, trans_nc))
    profile_index = Int(get(pm, "profile_index", 1))
    amf = Float64(get(pm, "amf", 1.0))
    sif_coeff_true = Float64.(get(pm, "sif_coeff_true", [1.0]))
    noise_scale = Float64(get(pm, "noise_scale", 0.0))
    n_samples = Int(get(pm, "n_samples", 5))

    println("Pseudo-measurement: Create_synthetic_radiance")
    println("  main_config: ", main_config)
    println("  trans_nc: ", trans_nc)
    println("  profile_index: ", profile_index, "  amf: ", amf)
    println("  sif_coeff_true: ", sif_coeff_true, "  n_samples: ", n_samples)

    cfg = TOML.parsefile(main_config)
    ctx = MWEF.prepare_mwe_inputs(main_config)
    data_cfg = get(cfg, "data", Dict{String, Any}())
    solar_file = get(data_cfg, "solar_file", "solar_merged_20240731_600_33300_100.out")
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)
    solar_hres, _ = MWEF.load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)
    solar_hres = Float64.(solar_hres)

    # Transmittance at high res from profile
    trans_hres = compute_trans_hres_from_profile(
        trans_nc,
        profile_index,
        ctx,
        amf;
        o2_sitp=ctx.o2_sitp,
        h2o_sitp=ctx.h2o_sitp,
    )
    println("  trans_hres range: ", minimum(trans_hres), " .. ", maximum(trans_hres))

    fit_cfg = get(cfg, "fit", Dict{String, Any}())
    n_legendre = Int(get(fit_cfg, "n_legendre", 3))
    n_sif = length(sif_coeff_true)
    length(sif_coeff_true) == size(ctx.sif_basis_hres, 2) || error("sif_coeff_true length must match SIF basis columns")

    # Build forward model and layout for retrieval
    fm = make_forward_model_simple(
        ctx,
        solar_hres;
        n_legendre=n_legendre,
        preallocate_float64=true,
    )
    layout = state_layout_simple(ctx; n_legendre=n_legendre)
    x0 = initial_state_simple(ctx; n_legendre=n_legendre, T=Float64)
    x_true = copy(x0)
    x_true[layout.idx_sif] .= sif_coeff_true

    # Generate synthetic observations and run retrieval (Legendre P0=1, others 0)
    leg_coeff_true = zeros(n_legendre + 1)
    leg_coeff_true[1] = 1.0

    retrieved_sif = Float64[]
    for i in 1:n_samples
        y_synth = generate_synthetic_spectrum(
            ctx,
            solar_hres,
            trans_hres;
            sif_coeff=sif_coeff_true,
            leg_coeff=leg_coeff_true,
            n_legendre=n_legendre,
        )
        if noise_scale > 0
            σ = noise_scale * (0.01 .+ 0.001 .* y_synth)
            y_synth .+= randn(length(y_synth)) .* σ
        end

        # Run retrieval (single-pixel flow)
        jacobian_eval = make_jacobian_evaluator(fm, x0; use_preallocated=false)
        x_a = copy(x0)
        x_a[layout.idx_sif] .= 0.0
        prior_sigma = fill(1e30, length(x0))
        prior_sigma[layout.idx_sif] .= 1e2
        S_a_inv = _spdiag_invvar(prior_sigma)
        x_scale = ones(Float64, length(x0))
        x_scale[layout.idx_sif] .= 1.0
        meas_sigma = 0.01
        S_e_inv = spdiagm(0 => fill(1.0 / (meas_sigma^2), length(y_synth)))

        x_curr = copy(x_a)
        λ = 1.0
        for _ in 1:15
            step = lm_one_step(
                fm,
                x_curr,
                y_synth;
                x_a=x_a,
                S_a_inv=S_a_inv,
                lambda=λ,
                lambda_up=5.0,
                lambda_down=0.7,
                lambda_min=1e-8,
                lambda_max=1e8,
                max_inner=12,
                jacobian_eval=jacobian_eval,
                x_scale=x_scale,
                meas_sigma=meas_sigma,
            )
            λ = step.lambda_next
            x_curr = step.x_next
            step.accepted || break
        end
        sif_ret = x_curr[layout.idx_sif]
        append!(retrieved_sif, sif_ret)
    end

    n_ev = length(sif_coeff_true)
    retrieved_sif = reshape(retrieved_sif, n_ev, n_samples)
    println("\n--- Synthetic experiment results ---")
    println("  True SIF coeff: ", sif_coeff_true)
    for ev in 1:n_ev
        println("  Retrieved SIF ev$ev: mean = ", mean(retrieved_sif[ev, :]), "  std = ", std(retrieved_sif[ev, :]))
    end
    return (sif_true=sif_coeff_true, sif_retrieved=retrieved_sif, ctx=ctx, layout=layout)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_synthetic()
end
