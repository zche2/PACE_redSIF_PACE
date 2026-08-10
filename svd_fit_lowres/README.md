# SVD fit low-res pipeline

Unified TOML for **PACE download** (Python, earthaccess) and **SVD swath retrieval** (Julia). Orchestration lives in `[general]`, `[download]`, `[svd_retrieval]`; physics and file paths in `[data]`, `[spectral]`, `[kernel]`, `[fit]`, `[batch_fit]`, etc. The Julia driver builds the effective retrieval settings **in memory** (no separate merged config file).

Local Julia helpers: [`julia/sif_basis.jl`](julia/sif_basis.jl) (SIF singular vectors), [`julia/pace_io.jl`](julia/pace_io.jl) (data path resolution), [`julia/merge_interim.jl`](julia/merge_interim.jl) (L1B+L2 → interim NetCDF).

## Configuration

Copy [`global_fit_pipeline.example.toml`](global_fit_pipeline.example.toml) to **`global_fit_pipeline.toml`** (gitignored) and edit:

- **`[general]`** — `L1B_dir`, `L2AOP_dir`, `L2BGC_dir` shared by download and SVD retrieval.
- **`[download.*]`** — earthaccess temporal/spatial/collections/paths/options.
- **`[svd_retrieval]`** — retrieval mode (below), `output_dir`, `interim_dir`, threading flags.
- **Retrieval tables** — `[data]`, `[spectral]`, `[fit]`, `[fit.svd]`, `[batch_fit]`, … (see example). Requires `[data].summer_nc`, `winter_nc`, `sif_file`, and (if band SNR is on) `pace_snr_file`.

## Retrieval modes

### (a) Granule-wise

```toml
[svd_retrieval]
mode = "granule"
granule_id = "20250701T000124"
# optional subset:
# pixel_start = 1
# pixel_end   = 1272
# scan_start  = 1
# scan_end    = 100
```

```bash
julia --project=. -t 8 svd_fit_lowres/svd_retrieval/run_svd_fit.jl svd_fit_lowres/global_fit_pipeline.toml
```

### (b) Global — one day or date range

```toml
[svd_retrieval]
mode = "date"
date = "20250701"
# or:
# mode = "dates"
# date_start = "20250703"
# date_end   = "20250707"
# date_interval_days = 1
# dates = ["20250701", "20250702"]
```

```bash
julia --project=. -t 8 svd_fit_lowres/svd_retrieval/run_svd_fit.jl svd_fit_lowres/global_fit_pipeline.toml
```

Granules run **sequentially** in one process (HDF5/NetCDF are not safe with concurrent granules). Use `parallel_pixels` + `julia -t N` for intra-granule parallelism, or several Julia jobs for multi-granule throughput.

### (c) Single pixel + spectra plot

```toml
[svd_retrieval]
granule_id = "20250701T000124"
pixel = 640
scan  = 100
output_dir  = "/path/to/output"
interim_dir = "/path/to/interim"
```

```bash
julia --project=. svd_fit_lowres/svd_retrieval/run_single_pixel.jl svd_fit_lowres/global_fit_pipeline.toml
```

Writes under `output_dir`:

- `svd_single_<granule>_p<_>_s<_>.jld2` — spectra + state
- `svd_single_<granule>_p<pixel>_s<scan>_spectra.png` — **left axis:** observed & model; **right axis (`twinx`):** SIF (TOA contribution) & residual

## Downloading inputs (Python)

Requires Python **3.11+** and `pip install earthaccess`.

```bash
python svd_fit_lowres/download/download_pace_products.py svd_fit_lowres/global_fit_pipeline.toml
```

## SVD retrieval details (Julia)

Preprocess/merge → per-pixel SVD transmittance LM on **sorted interim red bands** inside `(spectral.lambda_min_nm, spectral.lambda_max_nm)`. Transmittance PCs and SIF basis are interpolated onto those wavelengths; measurement noise uses `pace_snr_file` when band SNR is enabled.

Outputs NetCDF under `[svd_retrieval].output_dir` with suffix `[batch_fit].output_suffix_parallel` (`x_hat`, diagnostics, `state_names_csv`, …).

Use `parallel_pixels = true` and `julia -t N` (e.g. 8–16) for intra-granule multithreading:

```bash
julia --project=. -t 8 svd_fit_lowres/svd_retrieval/run_svd_fit.jl svd_fit_lowres/global_fit_pipeline.toml
```

## Rasterize SVD retrievals

Separate config: copy [`rasterize/rasterize.example.toml`](rasterize/rasterize.example.toml).

```bash
julia --project=. svd_fit_lowres/rasterize/rasterize_svd_retrievals.jl \
  svd_fit_lowres/rasterize/rasterize.example.toml
```

**Inputs:** granules matching `interim_<YYYYMMDDTHHmmss>_svd_retrieval_*.nc`. Optional `[rasterize].l2aop_dir` for `nflh` masking.

**Outputs:** `sif678_raster_<YYYYMMDD>_<YYYYMMDD>.nc` with `sif_radiance_678nm(lon, lat)` and `counts(lon, lat)`.

**`status_code` values:**

| Code | `converged` | Meaning |
|------|-------------|---------|
| 0 | 0 | Hit `max_outer_steps` without tolerance stop |
| 1 | 1 | Converged |
| 2 | 0 | LM step never accepted |
| 3 | 0 | Bad / missing spectrum |
| 4 | 0 | Exception in LM |
| 5 | 0 | Dark filter skip |
| 6 | 0 | Ocean filter skip |
| 7 | 0 | Pixel filter skip (`pixel_filter_vars`) |

Smoke tests: `julia --project=. svd_fit_lowres/rasterize/run_rasterize_tests.jl`
