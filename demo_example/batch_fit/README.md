# Full-swath batch retrieval

`Run_batch_full_nc.jl` runs the same retrieval logic as `../Run_batch_PACE_fit.jl` over the **entire** PACE .nc file, with a flexible pixel filter so only selected pixels are fitted. The script is **self-contained**: it does not include `Run_batch_PACE_fit.jl`; it only includes `../Fit_toy_forward_model.jl` and defines its own batch helpers and retrieval core (same behavior, no cross-file dependency).

## Usage

From the repo root or `demo_example`:

```bash
julia demo_example/batch_fit/Run_batch_full_nc.jl
```

Or set the config path:

```bash
PACE_MWE_CONFIG=/path/to/config.toml julia demo_example/batch_fit/Run_batch_full_nc.jl
```

## Config: `[batch_fit]`

Use the same config as the regular batch fit. Additional / overridden keys:

- **`pixel_filter_vars`** (default: `["nflh"]`) — List of variable names. A pixel `(i, j)` is eligible for retrieval only if **all** listed variables are non-missing at that pixel. Use `[]` to fit every pixel (no filter).
- Full swath is always used (no `pixel_start`, `pixel_end`, `scan_start`, `scan_end`).

Example in your config:

```toml
[batch_fit]
# Only fit pixels where nflh is present (not missing)
pixel_filter_vars = ["nflh"]
output_dir = "batch_output"
output_suffix = "_retrieval_full.nc"
# Optional: dark/ocean filters as in Run_batch_PACE_fit
# dark_filter_enabled = false
# ocean_filter_enabled = false
```

To add more criteria later, extend `pixel_filter_vars`, e.g. `["nflh", "quality_flag"]`, and ensure each variable is 2D with dimensions `(pixels, scans)` or `(scans, pixels)` and that “eligible” means non-missing (or adapt `build_pixel_eligible_mask` for other rules).

## Status codes

Same as `Run_batch_PACE_fit.jl`, plus:

- **7** = pixel_filter_skipped (pixel did not pass `pixel_filter_vars`).
