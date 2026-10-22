# retrieval_comparison

Tools to compare O2A SVD SIF retrievals of one PACE granule across retrieval
configurations, and to compute averaging kernels (AK) for selected pixels.
Default granule: `20250724T170203` (Gulf of Maine / Georges Bank).

Run every command from the repository root.

## Setup

- **Python**: `numpy`, `netCDF4`, `matplotlib`, `scipy`, `xarray`, `cartopy`
  (the last two only for `region_coeff_maps.py`).
- **Julia** (AK only): use the main checkout's environment, which has the
  pipeline packages installed. `Manifest.toml` is gitignored, so a fresh
  checkout or worktree has no instantiated environment of its own.

  ```bash
  julia --project=/home/zhe2/FraLab/PACE_redSIF_PACE -t 4 retrieval_comparison/ak_selected_pixels.jl
  ```

### Inputs

| What | Where |
|---|---|
| Retrieval output | `/home/zhe2/data/PACE/<run folder>/YYYY/MM/DD/interim_<granule>_svd_retrieval.nc` |
| L1B, L2 AOP (nFLH), L2 BGC (chl, flags) | `/kiwi-data/Data/satellite/PACE_OCI/{L1B_V3, L2_AOP_V3.2, L2_BGC_V3.2}/YYYY/MM/DD/` |
| Retrieval configs (AK only) | `svd_configs/global_fit_new_Sep26/`; if a config is missing there, `/home/zhe2/FraLab/PACE_redSIF_PACE/svd_configs/global_fit_new_Sep26/` (or pass `--config-dir`) |

### Run keys

All scripts refer to retrieval configurations by these keys. `compare_all_runs.py`
and `region_coeff_maps.py` include a run only when its output file is complete
(> 300 MB and readable), and otherwise skip it with a `[skip]` message.

The configs for `npoly2`, `newSVD_npoly2` and `newSVD_npoly3` are **not tracked in
git**. They exist only in the main checkout, which `ak_selected_pixels.jl` falls
back to.

| Key | Output folder | Config | SVD basis | Legendre terms |
|---|---|---|---|---|
| `npoly2` | `new_svd_retrieval_output_npoly2` | `new.singleGranule.npoly2nPC15_O2A.toml` | Aug01 | p0–p2 |
| `npoly3` | `new_svd_retrieval_output_O2A` | `new.singleGranule.npoly3nPC15_O2A.toml` | Aug01 | p0–p3 |
| `npoly5` | `new_svd_retrieval_output_npoly5` | `new.singleGranule.npoly5nPC15_O2A.toml` | Aug01 | p0–p5 |
| `newSVD_npoly2` | `new_svd_retrieval_output_O2A_newSVD_npoly2` | `new.newSVD.singleGranule.npoly2nPC15_O2A.toml` | scattering | p0–p2 |
| `newSVD_npoly3` | `new_svd_retrieval_output_O2A_newSVD` | `new.newSVD.singleGranule.npoly3nPC15_O2A.toml` | scattering | p0–p3 |

PC coefficients are only comparable between runs on the **same** SVD basis.

---

## 1. Averaging kernel for selected pixels

The swath output stores only `dof` and `diag(S_post)`. To get the full AK and
posterior covariance, the selected pixels are re-solved with the pipeline's own
`_run_one_svd_retrieval!` (`global_svd_fit_pipeline/svd_retrieval/svd_helpers.jl`),
using the same setup as `_process_one_pixel_svd!`. The script checks that the
re-run reproduces the stored retrieval.

### Step 1: select pixels

```bash
python retrieval_comparison/select_ak_pixels.py --run npoly3 --n 30 \
    --converged-in npoly2,npoly3,npoly5,newSVD_npoly2,newSVD_npoly3
```

| Option | Default | Meaning |
|---|---|---|
| `--run` | `npoly3` | run whose flags define the categories |
| `--converged-in` | same as `--run` | pixels must converge in **all** listed runs, so the same pixels can be compared across configurations |
| `--n` | 30 | pixels drawn at random per category (4 categories, so 120 pixels by default); fewer if a category has fewer candidates |
| `--granule` | `20250724T170203` | granule stamp |
| `--seed` | 0 | random seed for the pick |

Every pixel is clear (ocean, L2 flags ok, chl valid) open water (`is_dark == 1`).
The categories are `high_chl` (chl ≥ 1), `mid_chl` (0.1–0.3) and `low_chl` (< 0.05),
plus `nonconverged_lowchl` as a contrast case (not converged in `--run`).

Output: `out/ak/selected_pixels.csv`. You can also write this file by hand. The
Julia step needs the columns `category, scan, pixel, chl, nflh`, with **1-based**
`scan` and `pixel`.

### Step 2: compute AK and posterior covariance (Julia)

```bash
julia --project=/home/zhe2/FraLab/PACE_redSIF_PACE -t 4 retrieval_comparison/ak_selected_pixels.jl \
    --runs npoly3,newSVD_npoly3 [--granule 20250724T170203] \
    [--pixels retrieval_comparison/out/ak/selected_pixels.csv] [--config-dir DIR]
```

- The first run takes about 9 minutes. Most of that is Julia compilation and
  regenerating the pipeline's interim radiance file, because the pipeline
  deletes its own copy after each run. The file is cached in `out/_interim/`
  (about 1 GB). Later runs take 1–2 minutes per configuration.
- The log prints, for each pixel, the convergence flag, `DOF_sig`,
  `A[sif,sif]` and the re-run vs stored `sif_ev1`.
- Output: `out/ak/ak_<run>_<granule>.nc`.

| Variable | Meaning |
|---|---|
| `A` | averaging kernel; in Python `A[sel, i, j] = ∂x̂_i / ∂x_true,j` (`A = S_post·H_obs` at the final state) |
| `S_post` | posterior **covariance**, same layout; the diagonal holds **variances** |
| `x_hat`, `x_prior` | retrieved state and prior mean `[sel, state]` |
| `sigma_prior` | prior **standard deviation**, not variance |
| `dof_signal` | `tr(A)` |
| `dof_residual` | `n_band − tr(A)`; this is what the pipeline stores as `dof` |
| `reduced_chi2`, `converged`, `status` | from the re-run |
| `converged_stored`, `sif_ev1_stored`, `sif_ev1_rerun` | consistency check against the swath file |
| `wavelength`, `y_obs`, `y_fit` | measured and fitted spectra |
| `scan`, `pixel`, `category_id`, `chl`, `nflh` | pixel metadata |

The global attributes are `state_names_csv`, `category_names_csv`, `run`,
`config_file` and `granule`. In Julia (NCDatasets) the dimension order is reversed.

### Step 3: plots and summary

```bash
python retrieval_comparison/plot_ak.py             # all out/ak/ak_*.nc
python retrieval_comparison/plot_ak.py npoly3      # selected runs only
```

| Output | Content |
|---|---|
| `AK_heatmaps_<run>.png` | full AK, one pixel per category |
| `AK_sif_row_<run>.png` | (a) SIF row of the AK; (b) `A[sif,j]·σ_post,j`; (c) posterior `corr(sif, j)`; (d) `ρ_sif,j·σ_sif` in SIF₆₇₈ units, with ±1σ noise lines |
| `AK_run_comparison_<group>.png` | only with ≥ 2 runs that share categories: R² of SIF on the Legendre terms, σ(SIF₆₇₈), signal DOF, and (GCHPIO) retrieved − true SIF |
| `summary_<group>.txt` | per category: DOF, σ(SIF₆₇₈), posterior R² of SIF on the Legendre terms, corr(sif, p_k), corr(sif, α), max \|corr\| with the PCs (plus truth-based errors for GCHPIO) |

### Interpreting the results

- **`A[sif,sif]` is 1 by construction**, because SIF has a flat prior (σ = 1e12).
  The useful information is in the off-diagonals and in `S_post`.
- **Units.** σ_post,j = `sqrt(S_post[j,j])`, and SIF₆₇₈ = 0.9004 × `sif_ev1`.
- **Posterior correlation.** `ρ_sif,j = S_post[sif,j] / sqrt(S_post[sif,sif]·S_post[j,j])`.
  This is error coupling between SIF and j that comes from the measurement (spectral degeneracy).
- **Panel (b), prior leakage.** `A[sif,j]·σ_post,j = −ρ·σ_sif·(σ_post,j/σ_prior,j)²`.
  It is large only when the prior on j matters.
- **Panel (d), the SIF shift that goes with a 1σ_post error in j.** More generally,
  if j is forced or truncated to be off by δ, the best-fit SIF moves by
  `ΔSIF = (S_post[sif,j] / S_post[j,j])·δ`. For example, dropping p3 (npoly2 vs npoly3)
  shifts SIF by `−slope·p̂₃`, which reproduces the npoly2 − npoly3 SIF difference
  pixel by pixel (r = 1.000).
- **σ(SIF) with j known exactly.** This is `σ_sif·sqrt(1 − ρ²)` for one element.
  For several elements, delete their rows and columns from `inv(S_post)` and invert again.

### AK on the synthetic GCHPIO ensemble

The same AK analysis runs on the full-RT GCHPIO ensemble: a Lambertian ocean with GCHPIO
TOMAS aerosol optics on a 0.05 cm⁻¹ grid, 1000 samples. The true SIF is known there, so the
actual retrieval error can be compared with the posterior σ.

`ak_gchpio_ensemble.jl` re-solves the selected samples with `run_pseudo_svd_retrieval`
(`surrogate_meas/build_single_meas.jl`), the same call `retrieve_rt_ensemble.jl` uses, and
keeps `stats.averaging_kernel` and `stats.S_posterior`.

#### Runs

| Key | Stored ensemble retrieval (`rt_retrieval_1000_…`) | Config | SVD basis | O2 A band | Legendre |
|---|---|---|---|---|---|
| `gchpio_npoly2` | `O2A_SIF678BandMismatch_nPC15_npoly2.nc` | `svd_nPC15_npoly2_O2A.toml` | Aug01 | yes | p0–p2 |
| `gchpio_npoly3` | `O2A_SIF678BandMismatch_nPC15_npoly3.nc` | `svd_nPC15_npoly3_O2A.toml` | Aug01 | yes | p0–p3 |
| `gchpio_newSVD_npoly3` | `O2A_ScatteringIncludedSVD_nPC15_npoly3.nc` | `svd_nPC15_npoly3_O2A_newSVD.toml` | scattering | yes | p0–p3 |
| `gchpio_noO2A_npoly3` | `SIF678BandMismatch_nPC15_npoly3.nc` | `svd_nPC15_npoly3.toml` | Aug01 | no | p0–p3 |
| `gchpio_npoly5` | none (config-only run) | `svd_nPC15_npoly5_O2A.toml` | Aug01 | yes | p0–p5 |

- **Stored runs** take their config from the stored file's `svd_config` attribute, so it is
  exactly the config that produced them. The script checks that the re-run reproduces the
  stored SIF and convergence flag.
- **Config-only runs** (listed in `CONFIG_ONLY`) have no stored ensemble retrieval. The AK
  re-run is the retrieval, so there is no reproduction check, and the summary reports the
  run's convergence rate instead. They cannot drive the sample selection, so they need
  `--samples-from`.
- **Location.** Ensemble, stored retrievals and configs live in the main checkout (not in
  git): `surrogate_meas/full_RT_construction/output_new_Lambertian_gchpio_finerwavnum0p05/`
  and `surrogate_meas/configs/`.

#### Step 1: select samples and compute the AK for the stored runs

```bash
julia --project=/home/zhe2/FraLab/PACE_redSIF_PACE -t 4 retrieval_comparison/ak_gchpio_ensemble.jl \
    --runs gchpio_npoly2,gchpio_npoly3,gchpio_newSVD_npoly3,gchpio_noO2A_npoly3
```

| Option | Default | Meaning |
|---|---|---|
| `--runs` | `gchpio_npoly2,gchpio_npoly3,gchpio_newSVD_npoly3` | runs to compute; samples must converge in **all** of them (stored runs only) |
| `--n` | 30 | samples drawn at random per category |
| `--min-total` | 120 | if categories run short, top up round-robin to this many samples |
| `--seed` | 0 | random seed |
| `--samples-from` | (none) | reuse the samples and categories of an existing `ak_*_gchpio.nc`, skipping selection |

Categories come from the true ensemble state:

| Category | Condition |
|---|---|
| `clean_nadir` | AOD678 ≤ 25th percentile and VZA < 25° |
| `high_aod` | AOD678 ≥ 90th percentile |
| `high_vza` | VZA ≥ 50° |
| `low_sif` | true SIF678 < 0.05 |
| `nonconverged` | not converged in the **first** run of `--runs` |

Few samples converge in all the O2A runs, so some categories are small. With the four
stored runs: `clean_nadir` 15, `high_aod` 4, `high_vza` 34, `low_sif` 15. The top-up fills
the rest, mostly from `nonconverged`, giving 120 samples.

#### Step 2: add more runs on the same samples

```bash
julia --project=/home/zhe2/FraLab/PACE_redSIF_PACE -t 4 retrieval_comparison/ak_gchpio_ensemble.jl \
    --runs gchpio_npoly5 --samples-from retrieval_comparison/out/ak/ak_gchpio_npoly3_gchpio.nc
```

Any run can be added this way, stored or config-only, and it stays sample-for-sample
comparable with the runs from step 1.

#### Step 3: plots and summary

```bash
python retrieval_comparison/plot_ak.py gchpio_npoly2 gchpio_npoly3 gchpio_npoly5 \
    gchpio_newSVD_npoly3 gchpio_noO2A_npoly3
```

- **Output files.** `out/ak/ak_<run>_gchpio.nc` has the same layout as the granule files. There,
  `scan` is the 1-based sample index, `chl` holds AOD678 (attribute `meta_label`), and the
  file adds `sif678_true`, `sif678_ret`, `sigma678`, `vza`, `sza` and `surface_albedo`. For
  config-only runs, `sif_ev1_stored` is NaN and `converged_stored` is −1.
- **Grouping.** `plot_ak.py` keeps GCHPIO runs and granule runs in separate groups
  (`summary_gchpio.txt` vs `summary_granule.txt`, `AK_run_comparison_<group>.png`).
  Runs are compared only when they share the same categories.
- **Truth metrics.** The summary adds the mean error (retrieved − true SIF678) and
  `rms err/σ`. The ratio is about 1 when the posterior σ captures the actual error, and
  above 1 when forward-model (systematic) error dominates. The comparison figure has a
  fourth panel with retrieved − true SIF (mean ± sd) per category.
- **Run time.** About 20–25 minutes for four runs × 120 samples, including Julia
  compilation and one shared spectral setup per run. A single added run takes about 10 minutes.

#### Adding a GCHPIO run

- **With a stored ensemble retrieval:** add it to `STORED` in `ak_gchpio_ensemble.jl`
  (key → file name in the ensemble folder).
- **Config only:** add it to `CONFIG_ONLY` (key → config path) and run it with `--samples-from`.
- **Naming:** keep the `gchpio_` prefix, since `plot_ak.py` uses it to group runs.

---

## 2. Retrieval comparison

| Script | Command | Output | Time |
|---|---|---|---|
| `compare_all_runs.py` | `python retrieval_comparison/compare_all_runs.py [granule]` | `out/all_runs/summary.txt`, `F_all_runs_summary.png` | about 1 min |
| `region_coeff_maps.py` | `python retrieval_comparison/region_coeff_maps.py` | `out/region/A–E*.png`, `summary.txt` | about 2 min |
| `compare_O2A_oldSVD_vs_newSVD.py` | `python retrieval_comparison/compare_O2A_oldSVD_vs_newSVD.py [granule]` | `out/summary.txt`, `out/1–4*.png` | about 1 min |
| `map_O2A_oldSVD_vs_newSVD.py` | `python retrieval_comparison/map_O2A_oldSVD_vs_newSVD.py [granule]` | `out/5–7*.png` (whole-swath maps) | about 1 min |

**`compare_all_runs.py`** gives granule-wide statistics for every available run:
1. convergence status, reduced χ², residual and signal DOF, on four pixel sets:
   clear ocean, clear dark, clear bright, and all ocean;
2. posterior σ(SIF) and |SIF|/σ;
3. SIF@678 distribution and pixel-wise differences vs npoly3;
4. SIF by chl bin, with % negative and median nFLH;
5. SIF vs nFLH correlation, slope and intercept, overall and by chl.

**`region_coeff_maps.py`** makes Figure-2-style maps of the Georges Bank box for every run:
- `A`: nFLH, SIF, SNR, reduced χ², and convergence;
- `B1–B3`: PC 1–15 coefficients;
- `C`: α (physical: `10/(1+e^−raw)+1`) and Legendre p0–p3;
- `D`: Legendre p4–p5 and SIF;
- `E`: an A→B transect.

To change the region or the granule, edit the constants at the top of the
script (`GRANULE`, `NAME`, `CENTER`, `LAT_HALF`, `LON_HALF`, `TRANSECT`). The
first run reads the L1B RGB bands and caches the basemap in
`out/region/basemap_<granule>.npy`.

**`compare_O2A_oldSVD_vs_newSVD.py` and `map_O2A_oldSVD_vs_newSVD.py`** are the
earlier two-run versions (npoly3 vs newSVD_npoly3). `compare_all_runs.py`
supersedes them for statistics, but they still produce the whole-swath maps
(clear ocean, convergence, all pixels).

### Conventions

- **Clear** means ocean and L2 flags ok (none of ATMFAIL, LAND, HIGLINT, HILT,
  HISATZEN, STRAYLIGHT, CLDICE, HISOLZEN, CHLFAIL) and chl valid.
- **`is_dark == 1` means clear open water**, with max radiance ≤ 20. It does not
  mean night. `is_dark == 0` is bright scene: cloud, glint or swath edge.
- **SIF comparisons use pixels converged in all runs**, so adding a run shrinks the
  common sample.
- **The stored `dof` is the residual DOF** `n_band − tr(A)`, despite its NetCDF
  label. The signal DOF is `100 − dof` (100 bands in 640–780 nm).
- **Units.** `S_posterior_diag` holds variances, so σ = sqrt.

### Adding a run or a granule

- **New GCHPIO AK run:** see *Adding a GCHPIO run* above.
- **New run:** add its key in `compare_all_runs.py` (`RUNS_ALL`, `COLORS`),
  `region_coeff_maps.py` (`RUNS_ALL` with basis tag `"aug01"` or `"new"`, and the
  transect colours `C`), `select_ak_pixels.py` (`RUN_DIRS`) and
  `ak_selected_pixels.jl` (`RUN_CFG`, `RUN_OUTDIR`).
- **New granule:** pass the stamp as the argument (or `--granule`). The region
  script needs its constants edited. L1B and the V3.2 L2 files must exist on kiwi.

### Housekeeping

`*.png` and `*.nc` are gitignored. Only the scripts and the `summary.txt` files
are tracked. The caches `out/_interim/` (about 1 GB) and
`out/region/basemap_*.npy` (about 9 MB) are untracked and can be deleted at any
time; they are rebuilt when needed.
