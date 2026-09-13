# Global fit pipeline

Unified TOML for **PACE download** (Python, earthaccess) and **SVD swath retrieval** (Julia). Retrieval uses a **single** pipeline file: orchestration in `[general]`, `[download]`, `[svd_retrieval]`; physics and file paths in `[data]`, `[spectral]`, `[kernel]`, `[fit]`, `[batch_fit]`, etc. The Julia driver builds the effective retrieval settings **in memory** (no separate merged config file); output NetCDF `**config_file`** global attribute points at this pipeline TOML path.

Local Julia code under this folder includes `[julia/Simple_PACE_xSecFit_MWE_Functions.jl](julia/Simple_PACE_xSecFit_MWE_Functions.jl)` (same module as the repo’s historical `demo_example` copy, vendored here so the pipeline does not `include` `demo_example`) and `[julia/merge_interim.jl](julia/merge_interim.jl)` for L1B+L2 → interim NetCDF.

## Configuration

Copy `[global_fit_pipeline.example.toml](global_fit_pipeline.example.toml)` to `**global_fit_pipeline.toml**` and edit:

- `**[general]**` — `**L1B_dir**`, `**L2AOP_dir**`, `**L2BGC_dir**` shared by download and SVD retrieval (required for retrieval; download uses these when `[download.paths]` leaves a slot empty). Download writes granules under `<dir>/YYYY/MM/DD/` (kiwi-style); retrieval accepts that layout or a flat directory.
- `**[download.*]**` — earthaccess temporal/spatial/collections/paths/options.
- `**[svd_retrieval]**` — `**mode**`: `granule` (`**granule_id**`), `date` (`**date**` = `YYYYMMDD`), or `dates` (explicit `**dates**` list *or* `**date_start*`* + `**date_end**` with optional `**date_interval_days**`); optional `**pixel_*` / `scan_***`, `**output_dir**` (retrievals land under `output_dir/YYYY/MM/DD/`), `**interim_dir**`, `**skip_existing**` (default `false`; when `true`, skip if any `interim_<granule_id>_*.nc` exists under `output_dir/YYYY/MM/DD/` or flat `output_dir`, regardless of `output_suffix_parallel`), `**show_progress**` (default on when stdout is a TTY: granule bar per calendar day + scan bar per granule; safe with `**parallel_pixels**` / `**julia -t N**` because only the serial scan loop updates the bar), `**parallel_pixels**` (+ `**julia -t N**`); optional GPU tile retrieval `**use_gpu**`, `**gpu_tile_pixels**` (CUDA; CPU fallback if unavailable). Granules run **sequentially** in one process (HDF5/NetCDF are not safe with multi-threaded concurrent granules); use several Julia jobs for multi-granule throughput if needed.
- **Retrieval tables** — `[data]`, `[spectral]`, `[fit]`, `[fit.svd]`, `[batch_fit]`, … (see example). SVD retrieval uses **only** interim red bands and `**[spectral]`** λ bounds (no high-res LUT). Requires `**[data].summer_nc**`, `**winter_nc**`, `**sif_file**`, and (if band SNR is on) `**pace_snr_file**`. SNR toggle: `**[fit].use_band_snr**` or fallback `**[kernel].use_band_snr**`.

## Downloading inputs (Python)

Requires Python **3.11+** and `pip install earthaccess`. Auth: `[download.options] login_strategy` default `**all`**.

```bash
python global_fit_pipeline/download/download_pace_products.py global_fit_pipeline/global_fit_pipeline.toml
```

## SVD retrieval (Julia)

Preprocess/merge → per-pixel SVD transmittance LM on **sorted interim red bands** inside `**(spectral.lambda_min_nm, spectral.lambda_max_nm)`** (same slice for `Rtoa_red` and L1B `red_solar_irradiance`). Transmittance PCs and SIF basis are interpolated onto those wavelengths; measurement noise uses `**pace_snr_file**` when band SNR is enabled. Per-thread scratch buffers are sized for `Threads.maxthreadid()`.

```bash
julia --project=. -t 8 global_svd_fit_pipeline/svd_retrieval/run_svd_fit.jl path/to/pipeline.toml
```

### Multitask many days/granules (same TOML)

One Julia process still runs granules **sequentially** (NetCDF/HDF5-safe). To use many CPUs across days, launch **separate OS processes** that share the same pipeline file:

```bash
# one process per calendar day in [svd_retrieval] date range / dates list
MAX_JOBS=4 JULIA_THREADS=8 \
  ./global_svd_fit_pipeline/scripts/run_svd_fit_multitask.sh path/to/pipeline.toml

# or one process per L1B granule (more processes; lower JULIA_THREADS)
SPLIT=granule MAX_JOBS=8 JULIA_THREADS=4 \
  ./global_svd_fit_pipeline/scripts/run_svd_fit_multitask.sh path/to/pipeline.toml

# preview commands only
DRY_RUN=1 ./global_svd_fit_pipeline/scripts/run_svd_fit_multitask.sh path/to/pipeline.toml
```

Each worker calls `run_svd_fit.jl … --date YYYYMMDD` or `--granule ID` (CLI overrides `[svd_retrieval].mode`). Set `[svd_retrieval].skip_existing = true` so restarts skip finished granules. Logs go under `<output_dir>/multitask_logs/` (or `LOG_DIR=`).

With GPU, prefer fewer processes and `JULIA_THREADS=1`–`2`.

Outputs NetCDF under `[svd_retrieval].output_dir/YYYY/MM/DD/` with suffix `[batch_fit].output_suffix_parallel`, with variables aligned to the batch-fit swath style (`x_hat`, diagnostics, `state_names_csv`). By default `[batch_fit].save_posterior = true` writes `S_posterior_diag` (posterior variance per state element; ~380 MB/granule vs ~3 GB for the full covariance matrix). Set `save_posterior = false` to omit it.

The tropomi helper [`fetch_pace_earthaccess.py`](../toolbox/pace_tropomi_coincidence/fetch_pace_earthaccess.py) still uses **environment-only** login.

### Optional GPU tile path (SVD retrieval)

When `**[svd_retrieval].use_gpu = true`** and `**CUDA.functional()**` is true, the swath driver processes eligible pixels in **tiles** of at most `**gpu_tile_pixels`** (default 256). Each tile uses a **batched forward model** and **finite-difference Jacobian** on the GPU where possible; the Levenberg–Marquardt inner loop still matches the CPU `**lm_one_step`** acceptance and lambda logic (per column). If `**use_gpu**` is true but no usable CUDA device is available, the driver **warns once** and continues on the **CPU** per-pixel path.

#### Running GPU retrieval

Enable GPU in your TOML:

```toml
[svd_retrieval]
use_gpu       = true
gpu_tile_pixels = 256   # reduce if VRAM is tight
```

Then launch with a **low thread count** (1–2 recommended with GPU; see threading note below):

```bash
# GPU path — keep threads low to avoid CPU–GPU contention
julia --project=. -t 2 global_fit_pipeline/svd_retrieval/run_svd_fit.jl global_fit_pipeline/global_fit_pipeline.toml
```

CPU-only multithreaded (no GPU):

```bash
# CPU path — scale threads to available cores
julia --project=. -t 8 global_fit_pipeline/svd_retrieval/run_svd_fit.jl global_fit_pipeline/global_fit_pipeline.toml
```

#### Threading and GPU interaction

The LM inner loop over active columns within each tile always uses `Base.Threads.@threads` regardless of the `use_cuda` flag. The `use_gpu` flag only controls whether the **batched forward model and FD Jacobian** are computed on the GPU; the per-column LM step decisions remain on CPU in both cases.


| Mode                         | Batched predict / Jacobian | LM inner loop           | Recommended `-t` |
| ---------------------------- | -------------------------- | ----------------------- | ---------------- |
| CPU only (`use_gpu = false`) | CPU                        | `@threads` over columns | 8–16             |
| GPU (`use_gpu = true`)       | GPU (CUDA)                 | `@threads` over columns | 1–2              |


- **Configuration:** set `**use_gpu`** and `**gpu_tile_pixels**` under `**[svd_retrieval]**` in your pipeline TOML. The Julia driver merges these into the effective `**[batch_fit]**` dict before the swath code runs. If the same keys appear in `**[batch_fit]**`, `**[svd_retrieval]**` wins when that key is present in the `**[svd_retrieval]**` table.
- **VRAM:** work scales with tile size and spectral/state dimensions (Jacobian-shaped work per outer iteration). Reduce `**gpu_tile_pixels`** if you hit out-of-memory errors.
- **Threading:** the GPU scan path accumulates tiles **serially** in the pixel loop; avoid pairing it with very large `**julia -t N`** on the same socket if you see CPU–GPU contention. `**[batch_fit].use_threads**` still applies on the **CPU** path; the GPU path does not use the threaded per-pixel loop over the same region.
- **Tests:** `julia --project=. global_fit_pipeline/svd_retrieval/gpu/run_gpu_tests.jl` (CPU parity always; CUDA device checks are skipped when the driver returns errors such as code 999). Includes a parity + timing test comparing the CPU-threaded path against the GPU tile path over a K=32 pixel tile.

## Rasterize SVD retrievals (Julia)

Post-process per-granule SVD NetCDFs into **regular lat/lon maps** of mean `sif_radiance_678nm` and per-cell observation counts. Uses a **separate** config file (not the main pipeline TOML): copy [`rasterize/rasterize.example.toml`](rasterize/rasterize.example.toml) and edit `[rasterize]`, `[rasterize.grid]`, and `[rasterize.filters]`.

```bash
julia --project=. global_fit_pipeline/rasterize/rasterize_svd_retrievals.jl \
  global_fit_pipeline/rasterize/rasterize.example.toml
```

**Inputs:** granules in `input_dir` matching `interim_<YYYYMMDDTHHmmss>_svd_retrieval_*.nc` (same files as SVD retrieval output), either flat under `input_dir/` or under `input_dir/YYYY/MM/DD/`. Optional `[rasterize].l2aop_dir` points at PACE L2 OC_AOP files for `nflh` masking (`PACE_OCI.<granule_id>.L2.OC_AOP.*.nc`). If the exact L2 filename is missing, rasterize falls back to the same granule time **ignoring the last two digits** (seconds), e.g. retrieval `…211157` can use L2 `…211148`.

**Time windows:** window centers run from `start_date` to `end_date` every `chunk_frequency_days` days. Each window spans `[center − half_chunk_days, center + half_chunk_days]` and includes all granules whose sensing date falls in that range.

**Outputs:** one file per window under `output_dir` (created with `mkpath` if missing):

`sif678_raster_<YYYYMMDD>_<YYYYMMDD>.nc`

Variables: `sif_radiance_678nm(lon, lat)` (unweighted mean of valid soundings), `counts(lon, lat)`. Global attribute `status_codes` records the filter used.

**Grid:** `[rasterize.grid].resolution` sets cell size as `180 / resolution` degrees (e.g. `resolution = 180` → 1° cells, 360×180 grid).

**Filters** (`[rasterize.filters]`):

| Key | Purpose |
|-----|---------|
| `status_codes` | List of allowed `status_code` values (see table below). Default via legacy `valid_status_only = true` is `[1]`. |
| `exclude_dark` | Skip `is_dark == 1` |
| `ocean_only` | Keep only `is_ocean == 1` (open ocean) |
| `exclude_ocean` | Skip `is_ocean == 1` (land/coast); cannot be `true` with `ocean_only` |
| `max_sif` | Drop if \|sif_radiance_678nm\| exceeds limit |
| `max_chi2` | Drop if `reduced_chi2` exceeds limit (omit key to disable) |
| `exclude_missing_nflh` | Drop pixels where L2 OC_AOP `nflh` is missing (requires `[rasterize].l2aop_dir`) |
| `nflh_var` | L2 AOP variable name (default `nflh`, often under `geophysical_data/`) |

**`status_code` values** (variable `status_code` in each granule; `converged` is 1 only for code 1):

| Code | `converged` | Meaning |
|------|-------------|---------|
| 0 | 0 | Fit ran; hit `max_outer_steps` without RMSE/dx tolerance stop |
| 1 | 1 | Converged (tolerance stop or stalled LM) |
| 2 | 0 | LM step never accepted |
| 3 | 0 | Bad/missing spectrum or SZA; default for unprocessed pixels |
| 4 | 0 | Exception in `lm_one_step` |
| 5 | 0 | Dark filter skip |
| 6 | 0 | Ocean filter skip |
| 7 | 0 | Pixel filter skip (`pixel_filter_vars`) |

Typical map QC: `status_codes = [1]`. For converged plus iteration-cap fits: `status_codes = [0, 1]`.

**Console:** each window prints how many granule files are included, how many were read successfully, total valid soundings, and the output path.

**Smoke tests:** `julia --project=. global_fit_pipeline/rasterize/run_rasterize_tests.jl` (optional env `RASTER_TEST_RET`, `RASTER_TEST_L2AOP_DIR`).
