#!/usr/bin/env julia

using TOML
using Statistics
using LinearAlgebra

include(joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE_Functions.jl"))
using .SimplePACEXSecFitMWEFunctions

"""
Small test snippet:
1) Load high-res grid + kernel + cross-section interpolators from MWE context.
2) Read high-res solar spectrum.
3) Compute atmospheric transmission:
   T = exp(-(VCD_H2O * xs_H2O + VCD_O2 * xs_O2))
4) Multiply high-res solar by T.
5) Convolve to PACE low-resolution bands.
6) Read one observed PACE spectrum and compare residual statistics.
"""
function main()
    config_path = get(
        ENV,
        "PACE_MWE_CONFIG",
        joinpath(@__DIR__, "Simple_PACE_xSecFit_MWE.toml"),
    )
    cfg = TOML.parsefile(config_path)
    ctx = prepare_mwe_inputs(config_path)

    data_cfg = get(cfg, "data", Dict{String, Any}())
    test_cfg = get(cfg, "test_convolution", Dict{String, Any}())
    pace_cfg = get(cfg, "pace_observation", Dict{String, Any}())

    solar_file = get(data_cfg, "solar_file", "solar_merged_20200720_600_33300_100.out")
    solar_path = isabspath(solar_file) ? solar_file : joinpath(ctx.paths.base_dir, solar_file)

    # 1) Interpolate solar spectrum to the high-resolution wavelength grid FIRST.
    E_hres, solar_info = load_solar_spectrum_on_grid(solar_path, ctx.λ_hres; header_lines=3)

    p_o2 = Float64(get(test_cfg, "p_o2_hpa", 800.0))
    t_o2 = Float64(get(test_cfg, "t_o2_k", 250.0))
    vcd_o2 = Float64(get(test_cfg, "vcd_o2", 3.981071705534985e24)) # 10^24.6

    p_h2o = Float64(get(test_cfg, "p_h2o_hpa", 780.0))
    t_h2o = Float64(get(test_cfg, "t_h2o_k", 285.0))
    vcd_h2o = Float64(get(test_cfg, "vcd_h2o", 3.981071705534969e22)) # 10^22.6

    xs_o2 = ctx.o2_sitp(ctx.spectral_axis, p_o2, t_o2)
    xs_h2o = ctx.h2o_sitp(ctx.spectral_axis, p_h2o, t_h2o)
    T_atm = @. exp(-(vcd_h2o * xs_h2o + vcd_o2 * xs_o2))

    # 2) Multiply on the same high-res wavelength grid.
    solar_toa_hres = E_hres .* T_atm
    # 3) Convolve once to PACE low-resolution bands.
    solar_toa_pace = ctx.kernel.RSR_out * solar_toa_hres

    # 4) Read one real PACE spectrum and compare on the same low-res grid.
    do_obs_compare = Bool(get(pace_cfg, "enabled", true))
    pace_compare = nothing
    if do_obs_compare
        pace_file = get(pace_cfg, "pace_file", "sample_granule_20240830T131442_new_chl.nc")
        pace_path = isabspath(pace_file) ? pace_file : joinpath(ctx.paths.base_dir, pace_file)
        wavelength_var = String(get(pace_cfg, "wavelength_var", "red_wavelength"))
        spectrum_var = String(get(pace_cfg, "spectrum_var", "radiance_red"))
        pixel_idx = Int(get(pace_cfg, "pixel_index", 600))
        scan_idx = Int(get(pace_cfg, "scan_index", 800))

        y_obs, obs_info = load_pace_spectrum_on_grid(
            pace_path,
            ctx.λ;
            pixel_idx=pixel_idx,
            scan_idx=scan_idx,
            wavelength_var=wavelength_var,
            spectrum_var=spectrum_var,
        )

        y_mod = solar_toa_pace
        d_raw = y_obs .- y_mod
        rmse_raw = sqrt(sum(d_raw .^ 2) / length(d_raw))

        # Scale-only comparison for shape agreement.
        α = dot(y_obs, y_mod) / dot(y_mod, y_mod)
        d_scale = y_obs .- α .* y_mod
        rmse_scale = sqrt(sum(d_scale .^ 2) / length(d_scale))

        # Scale + offset comparison for first-order bias removal.
        A = hcat(y_mod, ones(length(y_mod)))
        θ = A \ y_obs
        y_fit = A * θ
        d_affine = y_obs .- y_fit
        rmse_affine = sqrt(sum(d_affine .^ 2) / length(d_affine))
        corr = cor(y_obs, y_mod)

        pace_compare = (
            obs_info = obs_info,
            y_obs = y_obs,
            y_mod = y_mod,
            rmse_raw = rmse_raw,
            rmse_scale = rmse_scale,
            rmse_affine = rmse_affine,
            scale_only = α,
            scale_affine = θ[1],
            offset_affine = θ[2],
            correlation = corr,
        )
    end

    println("Convolution test complete")
    println("  solar file: ", solar_path)
    println("  solar source range [nm]: ", solar_info.λ_source_range_nm)
    println("  high-res grid size: ", length(ctx.λ_hres))
    println("  pace bands: ", length(ctx.λ))
    println("  transmission min/max: ", extrema(T_atm))
    println("  high-res signal min/max: ", extrema(solar_toa_hres))
    println("  pace signal min/max: ", extrema(solar_toa_pace))
    println("  first 5 pace bands [nm]: ", ctx.λ[1:5])
    println("  first 5 convolved values: ", solar_toa_pace[1:5])

    if isnothing(pace_compare)
        println("  PACE observation compare: disabled")
    else
        cmp = pace_compare
        println("PACE observation comparison")
        println("  spectrum var: ", cmp.obs_info.spectrum_var)
        println("  wavelength var: ", cmp.obs_info.wavelength_var)
        println("  selected pixel/scan: ", (cmp.obs_info.pixel_idx, cmp.obs_info.scan_idx))
        println("  observation min/max: ", extrema(cmp.y_obs))
        println("  model min/max: ", extrema(cmp.y_mod))
        println("  corr(obs, model): ", cmp.correlation)
        println("  rmse raw: ", cmp.rmse_raw)
        println("  rmse after scale-only: ", cmp.rmse_scale, " (scale=", cmp.scale_only, ")")
        println("  rmse after scale+offset: ", cmp.rmse_affine, " (scale=", cmp.scale_affine, ", offset=", cmp.offset_affine, ")")
    end
end

main()
