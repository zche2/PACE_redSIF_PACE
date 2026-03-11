# Global Fit Pipeline

Pipeline for batch retrieval over raw PACE L1B + L2AOP + L2BGC data. Converts TOA reflectance (`rhot`) to radiance (`Rtoa`), merges into interim NetCDF, and runs the batch_fit retrieval engine (Fit_toy_forward_model, LM).

## Modes

| Mode | Description |
|------|--------------|
| `date` | Process all granules for a given date (YYYYMMDD) |
| `granule` | Process a single granule by ID (e.g. `20250927T131442`) |

## Config

Set `PACE_GLOBAL_FIT_CONFIG` or pass the config path as the first argument.

### `[global_fit]` section

| Key | Description |
|-----|-------------|
| `mode` | `"date"` or `"granule"` |
| `date` | YYYYMMDD when `mode=date` |
| `granule_id` | e.g. `"20250927T131442"` when `mode=granule` |
| `L1B_dir` | Directory containing L1B granules |
| `L2AOP_dir` | Directory containing L2 AOP granules |
| `L2BGC_dir` | Directory containing L2 BGC granules |
| `output_dir` | Output directory for retrieval NetCDF files |
| `interim_dir` | Directory for interim subset/merge files |
| `parallel_granules` | If true, parallelize over granules (threads) |
| `use_in_memory_merge` | If true, read L1B/L2 directly and write single interim (no subset files). Default: true for granule mode, false for date+parallel |

### `[pace_observation]` overrides

For preprocessed interim files, the pipeline overrides:

- `spectrum_var = "Rtoa_red"`
- `wavelength_var = "red_wavelength"`

### `[batch_fit]`

- `pixel_filter_vars = ["nflh"]` — only pixels with valid `nflh` are fitted
- `use_threads` — pixel-level parallelism (disabled when `parallel_granules=true` to avoid oversubscription)

## Usage

```bash
# Date mode (all granules for 2025-09-27)
julia -t 8 demo_example/global_fit/Run_global_fit.jl demo_example/global_fit/global_fit_config.toml
```

```bash
# Granule mode
# Set mode=granule and granule_id="20250927T131442" in config
julia -t 8 demo_example/global_fit/Run_global_fit.jl demo_example/global_fit/global_fit_config.toml
```

## Data layout

- **L1B**: `PACE_OCI.{date}T{time}.L1B.V3.nc` → granule_id = `"{date}T{time}"`
- **L2AOP**: `PACE_OCI.{granule_id}.L2.OC_AOP.V3_1.nc`
- **L2BGC**: `PACE_OCI.{granule_id}.L2.OC_BGC.V3_1.nc`

## Output

Retrieval NetCDF files are written to `output_dir` with the same naming as the batch_fit parallel output:

- `interim_{granule_id}_retrieval_full_parallel.nc`

## Dependencies

- Uses `demo_example/batch_fit/Run_batch_full_nc_parallel.jl` (via include)
- NCDatasets, TOML, Glob, Base.Threads
