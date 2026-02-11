"""
Create pseudo measurement - Understand transmittance construction
Standalone version (converted from Pluto notebook)
Created on 2025-11-13
Modified: 2026-02-06 (standalone version with profiling)
"""



using JLD2, Interpolations
using ForwardDiff, DiffResults, Plots, LinearAlgebra, DelimitedFiles, NCDatasets, Statistics
using Polynomials, Random
using LegendrePolynomials, Parameters, NonlinearSolve, BenchmarkTools
using PACE_SIF

println("=" ^ 80)
println("Starting Create_pseudo_v2_visualize_xSecFit (Standalone)")
println("=" ^ 80)

# ============================================================================
# Configuration - Update these paths as needed
# ============================================================================
const DATA_DIR = "/Users/cfranken/data"
const PACE_DATA_DIR = DATA_DIR
const XSECTION_DIR = joinpath(DATA_DIR, "default_run")
const REFERENCE_DIR = DATA_DIR
const KERNEL_DIR = DATA_DIR

# ============================================================================
# Parameters
# ============================================================================
λ_min = 630.
λ_max = 750.

# ============================================================================
# Load Transmittance and PACE Data
# ============================================================================
println("\n[1/8] Loading transmittance and PACE data...")
@time begin
    # MERRA2 generated
    summer_path = joinpath(PACE_DATA_DIR, "transmittance_summer_FineWvResModel_FullRange_Aug01.nc")
    winter_path = joinpath(PACE_DATA_DIR, "transmittance_winter_FineWvResModel_FullRange_Aug01.nc")
    
    if !isfile(summer_path)
        error("Summer transmittance file not found: $summer_path")
    end
    
    summer = Dataset(summer_path)
    
    # For now, use summer data twice if winter doesn't exist
    if isfile(winter_path)
        winter = Dataset(winter_path)
        trans = cat(summer["transmittance"][:,:], winter["transmittance"][:,:], dims=1)
        close(winter)
    else
        @warn "Winter transmittance file not found, using summer data only"
        trans = summer["transmittance"][:,:]
    end
    
    println("  Opened datasets.")
    println("  Concatenated! Shape: $(size(trans))")

    bands  = summer["band"][:]
    
    close(summer)
    
    # PACE data
    oci_path = joinpath(PACE_DATA_DIR, "sample_granule_20240830T131442_new_chl.nc")
    if !isfile(oci_path)
        error("PACE OCI file not found: $oci_path")
    end
    
    oci = Dataset(oci_path)
    red_band = oci["red_wavelength"][:]
    
    # select band (continuum spectrum) first to know which wavelengths to read
    ind      = findall( λ_min .< red_band .< λ_max )
    E        = oci["red_solar_irradiance"][ind]
    oci_band = red_band[ind]
    println("  Band selected to [$(λ_min), $(λ_max)], n_bands = $(length(oci_band))")
    
    # Read only the metadata arrays (small)
    nflh     = oci["nflh"][:, :]
    vza      = oci["sensor_zenith"][:, :]
    sza      = oci["solar_zenith"][:, :]
    chlor_a  = oci["chlor_a"][:, :]
    println("  Read metadata from PACE Dataset (nflh, angles, chlor_a)")
    
    # Don't read R_toa yet - we'll read only needed pixels later
end

# ============================================================================
# Interpolate Transmittance to OCI Bands
# ============================================================================
println("\n[2/8] Interpolating transmittance to OCI bands...")
@time begin
    trans_new = zeros(size(trans, 1), length(oci_band))
    
    for i in 1:size(trans, 1)
        # Create interpolator for this row (try to make this a spline later)
        itp_row = LinearInterpolation(bands, trans[i, :], extrapolation_bc=0)
        
        # Evaluate at OCI bands
        trans_new[i, :] = itp_row.(oci_band)
    end
    
    println("  Transmittance interpolated: $(size(trans_new, 1)) spectra × $(size(trans_new, 2)) bands")
end

# Visualization
p1 = plot(
    oci_band,
    trans_new[1:400:end,:]',
    label="",
    color=:silver,
    size=(600, 250),
    title="Transmittance spectra: $(size(trans_new,1))",
    titlefontsize=8,
)
plot!(oci_band, mean(trans_new', dims=2), lw=2.5, label="mean")
display(p1)

# ============================================================================
# Load and Interpolate SIF Shapes
# ============================================================================
println("\n[3/8] Loading and interpolating SIF shapes...")
@time begin
    # load SIF
    sif_path = joinpath(REFERENCE_DIR, "SIF_singular_vector.jld2")
    if !isfile(sif_path)
        error("SIF file not found: $sif_path")
    end
    SIF_shape_dict = JLD2.load(sif_path)
    
    # Create interpolator
    itp = interpolate(SIF_shape_dict["SIF_shapes"], (BSpline(Linear()), NoInterp()))
    range₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end]
    range₂ = 1:size(itp, 2)
    sitp = scale(itp, range₁, range₂)
    setp0 = extrapolate(sitp, 0)
    
    # Interpolate to OCI bands
    SIF_new = reduce(hcat, [setp0.(oci_band, i) for i in range₂])
    # Scale SIF
    scale_factor_SIF = 20
    SIF_new *= scale_factor_SIF
    println("  SIF shapes scaled by factor $scale_factor_SIF: $(size(SIF_new, 2)) shapes")
    
    # interpolation in the first dimension and no interp. in the second
    itp₂    = interpolate(SIF_shape_dict["SIF_U"], (BSpline(Linear()), NoInterp()))
    # scale
    r₁ = SIF_shape_dict["SIF_wavelen"][1]:SIF_shape_dict["SIF_wavelen"][end]
    r₂ = 1:size(itp₂, 2)
    sitp₂   = scale(itp₂, r₁, r₂)
    # set extrapolation filling value = 0
    setp0₂  = extrapolate(sitp₂, 0)
    # interpolation
    SIF_PC      = reduce(hcat, [setp0₂.(oci_band, i) for i in range₂])
    SIF_PC1_max = maximum(SIF_PC[:,1])
    # scale SIF principle components
    SIF_PC      /= SIF_PC1_max
    println("  SIF shape interpolated, PC scaled")
end

# Visualization
p2 = plot(
     oci_band,
     SIF_new[:, 1:5:end], 
     label="",
     size=(600, 250),
     title="wl SIF radiance (scaled), n_shapes: $(size(itp, 2))",
     titlefontsize=8
)
display(p2)

# ============================================================================
# Fit Reflectance from Non-Absorption Windows
# ============================================================================
println("\n[4/8] Fitting reflectance from baseline windows...")
@time begin
    # wavelength
    λ = oci_band
    
    # get center wavelength
    λc = center_wavelength(oci_band)

    # find transmittance baseline
    λ_bl_ref = [607.99, 610.36, 612.73, 615.14, 617.6, 620.06, 622.53, 669.52, 670.76, 671.99, 673.24, 674.51, 675.73, 676.96, 678.21, 679.45, 754.3, 779.33, 867.11, 869.61, 872.13]
    bl_ind = map(λ_bl_ref -> argmin(abs.(λ .- λ_bl_ref)), λ_bl_ref)
    
    # Extract indices where nFLH is detectable but negligible
    valid_mask  = findall(coalesce.(nflh .< 0.005, false))
    println("  Found $(length(valid_mask)) pixels with nFLH < 0.005")
    
    # Now read only the needed pixels from radiance_red (much more efficient!)
    println("  Reading radiance data for valid pixels only...")
    # Convert linear indices to 2D subscripts
    dims = size(nflh)
    valid_rows = [CartesianIndices(dims)[i][1] for i in valid_mask]
    valid_cols = [CartesianIndices(dims)[i][2] for i in valid_mask]
    
    # Read only the needed data
    n_valid = length(valid_mask)
    R_noSIF = zeros(n_valid, length(ind))
    for (idx, (row, col)) in enumerate(zip(valid_rows, valid_cols))
        R_noSIF[idx, :] = oci["radiance_red"][row, col, ind]
        if idx % 1000 == 0
            println("    Read $idx / $n_valid pixels")
        end
    end
    println("  Completed reading radiance data")
    
    sza_noSIF   = sza[valid_mask]
    vza_noSIF   = vza[valid_mask]
    R_baseline  = R_noSIF[:,bl_ind]
    E_baseline  = E[bl_ind]
    λ_baseline  = λ[bl_ind]
    λc_baseline = λc[bl_ind]

    println("  Baseline wavelengths: $(length(λ_baseline)) windows")
end

# Polynomial fit
println("  Performing polynomial fit (order=4)...")
@time begin
    order = 4
    n_pixels = size(R_baseline, 1)
    K₀       = hcat(collectPl.(λc_baseline, lmax=order)...)'
    K₀_recon = hcat(collectPl.(λc, lmax=order)...)'
    
    # Preallocate for fitted values and coefficients
    R_fitted = zeros(n_pixels, length(λc))
    coeffs_record = zeros(n_pixels, order + 1)  # Store coefficients for each pixel
    
    # Fit and record coefficients
    for i in 1:n_pixels
        # Fit (NOTE: This could be optimized - see profiling notes)
        coeff₀ = inv( K₀'K₀ )K₀' * ( R_baseline[i, :] .* pi ./ ( E_baseline .* cosd(sza_noSIF[i]) ) )

        # Evaluate at all wavelengths
        R_fitted[i, :] = K₀_recon * coeff₀ .* ( E .* cosd(sza_noSIF[i]) ./ pi)
        
        # Store coefficients
        coeffs_record[i, :] = coeff₀
    end
    println("  Reflectance fitted for $n_pixels pixels")
end

# Visualization
p3 = plot(
    λ,
    R_noSIF[1:250:end,:]',
    label="",
    size=(600, 250),
    title="Background pixel TOA radiance & Polynomial fit (order=$order), n_px: $n_pixels",
    titlefontsize=8
)
plot!(
    λ,
    R_fitted[1:250:end,:]',
    ls=:dash,
    label="",
)
display(p3)

# ============================================================================
# Load SNR Data
# ============================================================================
println("\n[5/8] Loading SNR data...")
@time begin
    filename = joinpath(PACE_DATA_DIR, "PACE_OCI_L1BLUT_baseline_SNR_1.1.txt")
    if !isfile(filename)
        error("SNR file not found: $filename")
    end
    lines = readlines(filename)
    end_header_index = findfirst(x -> x == "/end_header", lines)
    data  = readdlm(filename, String, skipstart=end_header_index)

    FPA   = data[:, 1]                   # 1st column: band
    wvlen = parse.(Float64, data[:, 2])  # 2nd column: center wavelength

    wv_val  = (λ_min .< wvlen .< λ_max)
    snr_ind = findall((FPA .== "Red") .& wv_val)

    # get c1 and c2 at that range
    c1    = parse.(Float64, data[snr_ind, 4])  # 4th column: c1
    c2    = parse.(Float64, data[snr_ind, 5])  # 5th column: c2
    println("  SNR coefficients retrieved: $(length(c1)) bands")
end

# ============================================================================
# Generate Pseudo Measurements
# ============================================================================
println("\n[6/8] Generating pseudo measurements...")
nᵨ = n_pixels   # also select SZA from
nₛ = size(SIF_new, 2)
nₜ = size(trans_new, 1)

# generate this amount of samples
n_sample = 5000
Random.seed!(512)
indᵨ     = rand(1:nᵨ, n_sample)
indₛ     = rand(1:nₛ, n_sample)
indₜ₁    = rand(1:nₜ, n_sample)
indₜ₂    = rand(1:nₜ, n_sample)
ind_sza  = rand(1:nᵨ, n_sample)   # sza
ind_vza  = rand(1:nᵨ, n_sample)   # vza

println("  Generating $n_sample pseudo measurements...")
@time begin
    # Preallocate storage for each component
    len_λ   = length(λ)
    ρ_all   = zeros(n_sample, len_λ)
    μ₁_all  = zeros(n_sample)
    μ₂_all  = zeros(n_sample)
    T₁_all  = zeros(n_sample, len_λ)
    T₂_all  = zeros(n_sample, len_λ)
    SIF_all = zeros(n_sample, len_λ)
    pseudo_obs_all = zeros(n_sample, len_λ)
    
    # Loop over samples and store each component
    for i in 1:n_sample
        # ----- rho -----
        ρ_all[i, :] = K₀_recon * coeffs_record[indᵨ[i], :]
        
        # ----- cos(sza) and cos(vza) -----
        μ₁_all[i] = 1.   # cosd(sza_noSIF[ind_sza[i]])
        μ₂_all[i] = 1.   # cosd(vza_noSIF[ind_vza[i]])
        
        # ----- Transmittance -----
        σ₁ = @. - 1 / μ₁_all[i] * log( trans_new[indₜ₁[i], :] )
        σ₂ = @. - 1 / μ₂_all[i] * log( trans_new[indₜ₂[i], :] )
        T₁_all[i, :] .= exp.(- σ₁)
        T₂_all[i, :] .= exp.(- σ₁ - σ₂)
        
        # ----- water-leaving SIF -----
        SIF_all[i, :] = SIF_new[:, indₛ[i]]
        
        # ----- TOA -----
        pseudo_obs_all[i, :] = @. E / pi * μ₁_all[i] * ρ_all[i, :] * T₂_all[i, :] + SIF_all[i, :] * T₁_all[i, :]

        # ----- noise -----
        stds = sqrt.(c1 .+ c2 .* pseudo_obs_all[i, :])
        pseudo_obs_all[i, :] += randn(len_λ) .* stds
        
        if i % 1000 == 0
            println("  Processed $i / $n_sample samples")
        end
    end
    println("  Completed $n_sample pseudo measurements")
end

# Visualization
Δn = 200
p4 = plot(layout=(5, 1), size=(800, 1000), legend=false)
plot!(p4[1], λ, ρ_all[1:Δn:end,:]', title="Reflectance", ylabel="ρ", xlabel="")
plot!(p4[2], λ, T₁_all[1:Δn:end,:]', title="Transmittance (one-way)", ylabel="T", xlabel="")
plot!(p4[3], λ, T₂_all[1:Δn:end,:]', title="Transmittance (two-way)", ylabel="T", xlabel="")
plot!(p4[4], λ, SIF_all[1:Δn:end,:]', title="SIF", ylabel="SIF", xlabel="")
plot!(p4[5], λ, pseudo_obs_all[1:Δn:end,:]', title="TOA", ylabel="Radiance", xlabel="Wavelength (nm)")
display(p4)

# ============================================================================
# Setup Retrieval Parameters
# ============================================================================
println("\n[7/8] Setting up retrieval parameters...")
@time begin
    println("  Loading LUT for cross sections...")
    
    o2_jld2 = joinpath(XSECTION_DIR, "default_run_O2.jld2")
    if !isfile(o2_jld2)
        error("O2 cross section file not found: $o2_jld2")
    end
    o2_sitp = read_rescale(o2_jld2)
    h2o_jld2 = joinpath(XSECTION_DIR, "default_run_H2O.jld2")
    if !isfile(h2o_jld2)
        error("H2O cross section file not found: $h2o_jld2")
    end
    h2o_sitp = read_rescale(h2o_jld2)
    metadata = joinpath(XSECTION_DIR, "Finer_Wavenumber_grid_FullRange_Aug01.log")
    ν_grid_o2, p_grid_hPa, t_grid = o2_sitp.ranges
    println("  LUT loaded.")

    
end

filename = "/Users/cfranken/data/PACE_OCI_RSRs.nc"
pace = Dataset(filename, "r")

wavlen = convert.(Float64, pace["wavelength"][:])
RSR_    = convert.(Float64, pace["RSR"][:])
band    = convert.(Float64, pace["bands"][:])
indiB   = findall(λ_min .< band .< λ_max)
indiW   = findall(λ_min .- 5 .< wavlen .< λ_max .+ 5)

# Create the wavenumber grid
ν_grid = λ_to_ν.(band[indiB])  # Convert your wavelength grid to wavenumber

MyKernel = KernelInstrument(
    band[indiB],           # band central wavelengths
    wavlen[indiW],         # fine wavelength grid for RSR
    RSR_[indiW, indiB],    # RSR matrix subset
    λ,                     # output wavelength grid (your main spectrum grid)
    ν_grid                 # wavenumber grid
)

close(pace)

# Setup retrieval parameters
println("  Setting up retrieval parameters...")
nPoly  = 6
nLayer = 1
nSIF   = 1

Sₐ   = I(nPoly+nSIF+13) .* 0.;
# update diagonal term
for i=1:(nPoly+1)
    Sₐ[i,i] = 1e10;
    # large variance applies no constrain to these polynomial term
end
# SIF magnitude
Sₐ[end, end] = 1;
# pressure terms
index_pressure = [1, 4, 7, 10] .+ nPoly .+ 1;
Sₐ[CartesianIndex.(index_pressure, index_pressure)] .= 1e4;
# temperature terms
index_temperature = [2, 5, 8, 11] .+ nPoly .+ 1;
Sₐ[CartesianIndex.(index_temperature, index_temperature)] .= 1e3;
# vcd for O2 (log10 transformed)
index_vcd_o2 = [3, 9] .+ nPoly .+ 1;
Sₐ[CartesianIndex.(index_vcd_o2, index_vcd_o2)]   .= 0.01;   # variance for O2 statistically is ~0.0025002
# vcd for H2O (log10 transformed)
index_vcd_h2o = [6, 12] .+ nPoly .+ 1;
Sₐ[CartesianIndex.(index_vcd_h2o, index_vcd_h2o)] .= 0.5;    # variance for H2O statistically is ~0.381829

println("  Sₐ diagonal: $(diag(Sₐ))")

# βₐ and γₐ
βₐ = [800.0; 250.0; 24.6; 800.0; 290.0; 22.6]              # see demo for the stats
γₐ = βₐ
println("  βₐ and γₐ set to climatology values")

# iteration method
thr_Converge = 1e-4
MyIter = (px, MyModel) -> LM_Iteration!(
    px, 
    MyModel;
    thr_Converge = thr_Converge,
    nIter = 20,
    γ_init = 1000.0,
    γ⁺ = 10.0,
    γ⁻ = 2.0,
    max_runtime = 400.0
)

params = RetrievalParams_xSecFit(
    λ  = λ,
    λc = λc,
    E  = E,
    c₁ = c1,
    c₂ = c2,
    o2_sitp  = o2_sitp,
    h2o_sitp = h2o_sitp,
    InstrumentKernel = MyKernel,
    forward_model    = forward_model,
    nPoly  = nPoly,
    nLayer = nLayer,
    nSIF   = nSIF,
    Sₐ = Sₐ,
    βₐ = βₐ,
    γₐ = γₐ,
    SIFComp  = SIF_PC,      
    iteration_method = MyIter,
    thr_Converge = thr_Converge
)
println("  Retrieval parameters set.")

# ============================================================================
# Single Pixel Retrieval
# ============================================================================
println("\n[8/8] Running single pixel retrieval...")
i = 51
println("  Running retrieval for sample $i...")
@time begin
    SinglePixel = Retrieval_for_Pixel(
        pseudo_obs_all[i, :],
        0.,  # sza_noSIF[ind_sza[i]],
        0.,  # vza_noSIF[ind_sza[i]],
        maximum(SIF_all[i, :]),
        1.0,
        1.0,
        params
    )
end

println("\nRetrieval complete!")
println("SinglePixel result:")
println(SinglePixel)

# ============================================================================
# Reconstruct Components
# ============================================================================
println("\nReconstructing components...")
MyReconModel = (x, px) -> forward_model(
    x,
    px, 
    params;
    return_components=true
)

_, ρ, T₁, T₂, SIF = MyReconModel(SinglePixel.x, SinglePixel)

# Manual reconstruction check
xᵨ    = SinglePixel.x[1 : SinglePixel.nPoly+1];
v     = collectPl.(SinglePixel.λc, lmax=SinglePixel.nPoly);
thisρ = hcat(v...)' * xᵨ;

# ============================================================================
# Visualization of Results
# ============================================================================
println("\nGenerating result plots...")

# Reflectance
p_rho = plot(
    λ, ρ_all[i,:], size=(800, 300),
    xticks = (620:10:860, string.(620:10:860)),
    label="truth",
    title="Reflectance"
)
plot!(λ, ρ, label="reconst.")
plot!(λ, thisρ, label="manual reconstruct")
display(p_rho)

# T₂
p_t2 = plot(
    λ, [T₂_all[i,:] T₂], size=(800, 300),
    xticks = (620:10:860, string.(620:10:860)),
    title="T₂",
    label=["truth" "recon."]
)
display(p_t2)

# T₁
p_t1 = plot(
    λ, [T₁_all[i,:] T₁], size=(800, 300),
    xticks = (620:10:860, string.(620:10:860)),
    title="T₁",
    label=["truth" "recon."]
)
display(p_t1)

# SIF
p_sif = plot(
    λ, [SIF_all[i,:] SIF], size=(800, 300),
    xticks = (620:10:860, string.(620:10:860)),
    title="SIF",
    label=["truth" "recon."]
)
display(p_sif)

# TOA radiance
p_toa = plot(
    λ, [SinglePixel.R_toa SinglePixel.y], size=(800, 300),
    xticks = (620:10:860, string.(620:10:860)),
    title="TOA radiance",
    label=["truth" "recon."]
)
display(p_toa)

# Residual
p_res = plot(
    λ, (SinglePixel.R_toa .- SinglePixel.y), size=(800, 300),
    xticks = (620:10:860, string.(620:10:860)),
    title="residual"
)
display(p_res)

println("\n" * "=" ^ 80)
println("Script completed successfully!")
println("=" ^ 80)
