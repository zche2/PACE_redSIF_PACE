# Simple PACE xSecFit MWE

This MWE is intentionally limited to **data preparation and kernel QA**.
It does **not** run retrieval or optimization.

## Files

- `demo_example/Simple_PACE_xSecFit_MWE.jl`
  - Minimal driver script.
  - Reads TOML config and prints prepared inputs + kernel diagnostics.
- `demo_example/Simple_PACE_xSecFit_MWE_Functions.jl`
  - Utility functions for path resolution, SIF EV interpolation, kernel generation, and kernel comparison.
- `demo_example/Simple_PACE_xSecFit_MWE.toml`
  - All path and processing settings.
  - Includes high-res spectral grid controls (`lambda_min_nm`, `lambda_max_nm`, `delta_lambda_nm`).

## Run

```bash
julia --project=. demo_example/Simple_PACE_xSecFit_MWE.jl
```

Optional overrides:

```bash
PACE_MWE_CONFIG=/path/to/your_config.toml julia --project=. demo_example/Simple_PACE_xSecFit_MWE.jl
PACE_DATA_DIR=/path/to/data PACE_XSEC_RUN=wavelength-run julia --project=. demo_example/Simple_PACE_xSecFit_MWE.jl
```

## Kernel QA Metrics

`Kernel generation` block:

- `row-sum min/max`
  - Should be close to `1.0`; rows are normalized in `KernelInstrument`.
- `negative weight fraction`
  - Ideally near `0`; nonzero values suggest negative RSR entries propagated.
- `center offset mean [nm]`
  - Mean of `(weighted-center - nominal band center)`.
- `center offset abs max [nm]`
  - Worst absolute weighted-center offset.
- `peak offset abs max [nm]`
  - Worst absolute peak-location offset.
- `low-res λ source`
  - Either PACE nominal centers or kernel-weighted effective centers (configurable via `use_effective_band_centers`).

`Kernel comparison (stored vs regenerated)` block:

- `common bands`
  - Number of band centers matched within tolerance.
- `L1 / L2 / max|Δ|`
  - Row-wise shape difference after mapping stored kernel rows onto regenerated wavelength grid and renormalizing each row.

## For Your Custom Fit

Use prepared context from the driver:

- `ctx.λ`
- `ctx.λc`
- `ctx.kernel`
- `ctx.o2_sitp`
- `ctx.h2o_sitp`
- `ctx.λ_hres` (high-res forward-model wavelength grid)
- `ctx.sif_basis_hres` (SIF EVs on full high-res forward-model wavelength grid)
- `ctx.sif_basis` (low-res convolved SIF EVs aligned with `ctx.λ`)

Then implement your own forward model and inversion routine on top of this.
