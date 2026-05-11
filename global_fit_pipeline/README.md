# Global fit pipeline

Unified TOML for **PACE download** (Python, earthaccess) and **SVD swath retrieval** (Julia). Retrieval uses a **single** pipeline file: orchestration in `[general]`, `[download]`, `[svd_retrieval]`; physics and file paths in `[data]`, `[spectral]`, `[kernel]`, `[fit]`, `[batch_fit]`, etc. The Julia driver builds the effective retrieval settings **in memory** (no separate merged config file); output NetCDF **`config_file`** global attribute points at this pipeline TOML path.

Local Julia code under this folder includes [`julia/Simple_PACE_xSecFit_MWE_Functions.jl`](julia/Simple_PACE_xSecFit_MWE_Functions.jl) (same module as the repo’s historical `demo_example` copy, vendored here so the pipeline does not `include` `demo_example`) and [`julia/merge_interim.jl`](julia/merge_interim.jl) for L1B+L2 → interim NetCDF.

## Configuration

Copy [`global_fit_pipeline.example.toml`](global_fit_pipeline.example.toml) to **`global_fit_pipeline.toml`** and edit:

- **`[general]`** — **`L1B_dir`**, **`L2AOP_dir`**, **`L2BGC_dir`** shared by download and SVD retrieval (required for retrieval; download uses these when `[download.paths]` leaves a slot empty).
- **`[download.*]`** — earthaccess temporal/spatial/collections/paths/options.
- **`[svd_retrieval]`** — **`mode`**: `granule` (**`granule_id`**), `date` (**`date`** = `YYYYMMDD`), or `dates` (explicit **`dates`** list *or* **`date_start`** + **`date_end`** with optional **`date_interval_days`**); optional **`pixel_*` / `scan_*`**, **`output_dir`**, **`interim_dir`**, **`parallel_pixels`** (+ **`julia -t N`**). Granules run **sequentially** in one process (HDF5/NetCDF are not safe with multi-threaded concurrent granules); use several Julia jobs for multi-granule throughput if needed.
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
