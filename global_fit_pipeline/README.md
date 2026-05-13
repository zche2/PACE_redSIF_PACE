# Global fit pipeline

Unified TOML for **PACE download** (Python, earthaccess) and **SVD swath retrieval** (Julia). Retrieval uses a **single** pipeline file: orchestration in `[general]`, `[download]`, `[svd_retrieval]`; physics and file paths in `[data]`, `[spectral]`, `[kernel]`, `[fit]`, `[batch_fit]`, etc. The Julia driver builds the effective retrieval settings **in memory** (no separate merged config file); output NetCDF **`config_file`** global attribute points at this pipeline TOML path.

Local Julia code under this folder includes [`julia/Simple_PACE_xSecFit_MWE_Functions.jl`](julia/Simple_PACE_xSecFit_MWE_Functions.jl) (same module as the repo’s historical `demo_example` copy, vendored here so the pipeline does not `include` `demo_example`) and [`julia/merge_interim.jl`](julia/merge_interim.jl) for L1B+L2 → interim NetCDF.

## Configuration

Copy [`global_fit_pipeline.example.toml`](global_fit_pipeline.example.toml) to **`global_fit_pipeline.toml`** and edit:

- **`[general]`** — **`L1B_dir`**, **`L2AOP_dir`**, **`L2BGC_dir`** shared by download and SVD retrieval (required for retrieval; download uses these when `[download.paths]` leaves a slot empty).
- **`[download.*]`** — earthaccess temporal/spatial/collections/paths/options.
- **`[svd_retrieval]`** — **`mode`**: `granule` (**`granule_id`**), `date` (**`date`** = `YYYYMMDD`), or `dates` (explicit **`dates`** list *or* **`date_start`** + **`date_end`** with optional **`date_interval_days`**); optional **`pixel_*` / `scan_*`**, **`output_dir`**, **`interim_dir`**, **`parallel_pixels`** (+ **`julia -t N`**); optional GPU tile retrieval **`use_gpu`**, **`gpu_tile_pixels`** (CUDA; CPU fallback if unavailable). Granules run **sequentially** in one process (HDF5/NetCDF are not safe with multi-threaded concurrent granules); use several Julia jobs for multi-granule throughput if needed.
- **Retrieval tables** — `[data]`, `[spectral]`, `[fit]`, `[fit.svd]`, `[batch_fit]`, … (see example). SVD retrieval uses **only** interim red bands and **`[spectral]`** λ bounds (no high-res LUT). Requires **`[data].summer_nc`**, **`winter_nc`**, **`sif_file`**, and (if band SNR is on) **`pace_snr_file`**. SNR toggle: **`[fit].use_band_snr`** or fallback **`[kernel].use_band_snr`**.

## Downloading inputs (Python)

Requires Python **3.11+** and `pip install earthaccess`. Auth: `[download.options] login_strategy` default **`all`**.

```bash
python global_fit_pipeline/download/download_pace_products.py global_fit_pipeline/global_fit_pipeline.toml
```

## SVD retrieval (Julia)

Preprocess/merge → per-pixel SVD transmittance LM on **sorted interim red bands** inside **`(spectral.lambda_min_nm, spectral.lambda_max_nm)`** (same slice for `Rtoa_red` and L1B `red_solar_irradiance`). Transmittance PCs and SIF basis are interpolated onto those wavelengths; measurement noise uses **`pace_snr_file`** when band SNR is enabled. Per-thread scratch buffers are sized for `Threads.maxthreadid()`.

```bash
julia --project=. -t 8 global_fit_pipeline/svd_retrieval/run_svd_fit.jl global_fit_pipeline/global_fit_pipeline.toml
```

Outputs NetCDF under `[svd_retrieval].output_dir` with suffix `[batch_fit].output_suffix_parallel`, with variables aligned to the batch-fit swath style (`x_hat`, diagnostics, `state_names_csv`).

The tropomi helper [`fetch_pace_earthaccess.py`](../toolbox/pace_tropomi_coincidence/fetch_pace_earthaccess.py) still uses **environment-only** login.

### Optional GPU tile path (SVD retrieval)

When **`[svd_retrieval].use_gpu = true`** and **`CUDA.functional()`** is true, the swath driver processes eligible pixels in **tiles** of at most **`gpu_tile_pixels`** (default 256). Each tile uses a **batched forward model** and **finite-difference Jacobian** on the GPU where possible; the Levenberg–Marquardt inner loop still matches the CPU **`lm_one_step`** acceptance and lambda logic (per column). If **`use_gpu`** is true but no usable CUDA device is available, the driver **warns once** and continues on the **CPU** per-pixel path.

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

| Mode | Batched predict / Jacobian | LM inner loop | Recommended `-t` |
|---|---|---|---|
| CPU only (`use_gpu = false`) | CPU | `@threads` over columns | 8–16 |
| GPU (`use_gpu = true`) | GPU (CUDA) | `@threads` over columns | 1–2 |

- **Configuration:** set **`use_gpu`** and **`gpu_tile_pixels`** under **`[svd_retrieval]`** in your pipeline TOML. The Julia driver merges these into the effective **`[batch_fit]`** dict before the swath code runs. If the same keys appear in **`[batch_fit]`**, **`[svd_retrieval]`** wins when that key is present in the **`[svd_retrieval]`** table.
- **VRAM:** work scales with tile size and spectral/state dimensions (Jacobian-shaped work per outer iteration). Reduce **`gpu_tile_pixels`** if you hit out-of-memory errors.
- **Threading:** the GPU scan path accumulates tiles **serially** in the pixel loop; avoid pairing it with very large **`julia -t N`** on the same socket if you see CPU–GPU contention. **`[batch_fit].use_threads`** still applies on the **CPU** path; the GPU path does not use the threaded per-pixel loop over the same region.
- **Tests:** `julia --project=. global_fit_pipeline/svd_retrieval/gpu/run_gpu_tests.jl` (CPU parity always; CUDA device checks are skipped when the driver returns errors such as code 999). Includes a parity + timing test comparing the CPU-threaded path against the GPU tile path over a K=32 pixel tile.
