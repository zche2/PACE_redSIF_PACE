# Pseudo-measurement synthetic data experiments

This folder contains scripts for synthetic data experiments following the idea of
`scripts/pseudo_measurement`: generate TOA radiance from vertical profiles (72-layer MERRA2)
and SIF shapes, convolve to OCI resolution, then run the retrieval.

## Scripts

- **Create_synthetic_radiance.jl**  
  - Load MERRA2-derived transmittance NetCDF (72 layers; see `scripts/generate_transmittance_spectra/Generate_transmittance.jl`).
  - Pick profile(s), compute transmittance at high spectral resolution using LUTs.
  - Choose SIF shapes and coefficients.
  - Compute TOA radiance at high resolution: `solar * trans_2way + sif_radiance * trans_1way`.
  - Convolve to OCI resolution with the instrument kernel.
  - Optionally add noise.
  - Run the retrieval (same pipeline as `Fit_toy_forward_model.jl`) and compare retrieved state to truth.

- **SIF_addition.jl**  
  - Load real PACE spectrum from a .nc file.
  - Add a known SIF radiance (SIF basis × chosen coefficients, convolved to OCI).
  - Run the retrieval on the modified spectrum.
  - Compare retrieved SIF coefficients to the added values to check if pseudo SIF is recovered.

## Config

- **pseudo_measurement_config.toml** — Paths for transmittance NetCDF, main MWE config, and options.
- The main retrieval config (e.g. `../Simple_PACE_xSecFit_MWE_zcheVer.toml`) is used for LUTs, kernel, SIF basis, and fit settings.

## Data requirements

- MERRA2 transmittance NetCDF produced by `Generate_transmittance.jl` (variables: `transmittance`, `temperature`, `pressure`, `vcd_dry`, `vcd_h2o`, `band`; attributes: `ak`, `bk`).
- Same cross-section LUTs, RSR, solar, and SIF basis as the main demo (via the main config’s `[data]`).
