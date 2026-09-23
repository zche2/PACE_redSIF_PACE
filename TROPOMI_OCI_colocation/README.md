# TROPOMI ↔ PACE OCI co-location

Build co-located TROPOMI BD5 and PACE OCI L1B matchups, then optionally retrieve SIF on those pairs.

**Preferred entry points (current):**

| Goal | Script |
|------|--------|
| Full pipeline → `products/matchup.nc` | `run_colocation.py` |
| Granule pairing JSON only | `co_locate.py` |
| SVD SIF on matchups (Julia) | `run_matchup_sif.py` |
| Join independent TROPOMI/PACE SIF onto matchups | `join_independent_sif.py` |

Notebooks under `notebooks/` are exploratory / legacy. Prefer the CLIs above.

---

## Layout

```
TROPOMI_OCI_colocation/
  run_colocation.py          # full staged pipeline
  co_locate.py               # step-1 pairing → JSON
  run_matchup_sif.py         # Julia SIF wrapper
  join_independent_sif.py    # join external SIF products
  configs/
    config.yaml              # active colocation config
    config.example.yaml      # template
    config.sif_retrieval.toml
  colocation/                # pipeline package (discover, pair, match, …)
  sif_retrieval/             # Julia matchup SIF runner
  notebooks/                 # older notebooks + pairing JSON dumps
```

Run commands from the **repo root** (`PACE_redSIF_PACE/`) or from this folder; both work because the CLIs add their own directory to `sys.path`.

---

## Data roots (dated layout)

| Sensor | Root | On-disk layout |
|--------|------|----------------|
| PACE OCI L1B V3 | `/kiwi-data/Data/satellite/PACE_OCI/L1B_V3` | `YYYY/MM/DD/PACE_OCI.*.L1B.V3.nc` |
| TROPOMI BD5 L1 | `/net/squid/data1/projects/TROPOMI/ESA/L1` | `YYYY/MM/S5P_*_L1B_RA_BD5_*.nc` |

`colocation/discover.py` prefers these dated trees and falls back to a flat root / `rglob` if the dated folders are missing.

Also required for the full pipeline:

- OCI RSR NetCDF (`rsr_path`)
- OCI SNR LUT text file (`snr_path`)

---

## 1. Full pipeline → `matchup.nc`

Stages (cached under `output_dir`):

```
pairs → matches → svd → rsr → simulate → products/matchup.nc
```

### Config

Copy and edit:

```bash
cp TROPOMI_OCI_colocation/configs/config.example.yaml \
   TROPOMI_OCI_colocation/configs/config.yaml
```

Important fields:

- `pace_dir`, `tropomi_dir`, `rsr_path`, `snr_path`
- `output_dir` — all caches + final product
- `year`, `month`, `day_start`, `day_end`
- pairing / pixel / SVD / noise settings

### Run from scratch

```bash
python TROPOMI_OCI_colocation/run_colocation.py \
  --config TROPOMI_OCI_colocation/configs/config.yaml
```

Creates `output_dir` if needed and writes:

```
output_dir/
  pairs/          # swath list + meta
  matches/        # pixel matches
  svd/            # TROPOMI SVD basis
  rsr/            # RSR weights
  products/
    matchup.nc    # final co-location product
  run_meta.json
```

### Resume / rebuild stages

Stages are fingerprint-cached. Rebuild only what you need:

```bash
# Start mid-pipeline (earlier stages must already exist)
python TROPOMI_OCI_colocation/run_colocation.py \
  --config TROPOMI_OCI_colocation/configs/config.yaml \
  --from simulate

# Force rebuild of named stages (and dependents as coded)
python TROPOMI_OCI_colocation/run_colocation.py \
  --config TROPOMI_OCI_colocation/configs/config.yaml \
  --force svd simulate

# Force everything
python TROPOMI_OCI_colocation/run_colocation.py \
  --config TROPOMI_OCI_colocation/configs/config.yaml \
  --force all
```

Valid stage names: `pairs`, `matches`, `svd`, `rsr`, `simulate`.

A week of data can take a long time; redirect logs if running in the background:

```bash
nohup python -u TROPOMI_OCI_colocation/run_colocation.py \
  --config TROPOMI_OCI_colocation/configs/config.yaml \
  > /path/to/output_dir/run.log 2>&1 &
tail -f /path/to/output_dir/run.log
```

---

## 2. Pairing only → JSON (`co_locate.py`)

Use this when you only need TROPOMI swaths with time+bbox-overlapping PACE granules (same idea as the cleanup notebook). **Does not** write `matchup.nc`.

```bash
# Full month JSON
python TROPOMI_OCI_colocation/co_locate.py \
  --year 2025 --month 1 --day-end 31 \
  --output TROPOMI_OCI_colocation/notebooks/co-location_results_202501.json

# Smoke test
python TROPOMI_OCI_colocation/co_locate.py \
  --year 2025 --month 1 --day-end 7 \
  --output /tmp/coloc_smoke.json
```

Defaults for `--pace-dir` / `--tropomi-dir` match the dated roots above. Override with flags if needed (`--help`).

---

## 3. SIF retrieval on matchups

After `products/matchup.nc` exists, edit `configs/config.sif_retrieval.toml` (`matchup_nc`, `output_root`, spectra, fit window, …), then:

```bash
python TROPOMI_OCI_colocation/run_matchup_sif.py \
  --config TROPOMI_OCI_colocation/configs/config.sif_retrieval.toml \
  -t 16

# Quick test on N matches
python TROPOMI_OCI_colocation/run_matchup_sif.py \
  --config TROPOMI_OCI_colocation/configs/config.sif_retrieval.toml \
  --max-matches 50
```

This wraps `sif_retrieval/run_matchup_sif.jl` with the repo Julia project.

---

## 4. Join independent SIF products

`join_independent_sif.py` attaches **already-run** TROPOMI and PACE SIF retrievals onto each row of a colocation `matchup.nc`. This is separate from `run_matchup_sif.py` (which retrieves SIF on the simulated / measured OCI spectra *inside* the matchup file).

### What you need

| Input | Role | Expected files |
|-------|------|----------------|
| `matchup.nc` | Co-location product from `run_colocation.py` | Must contain `trop_pix`, `trop_scan`, `pace_pix`, `pace_scan`, `tropomi_path`, `pace_path`, plus lat/lon |
| `--trop-dir` | Independent TROPOMI BD5 SIF | `fs_ret_<ORBIT>_b5_v10.h5` (orbit from TROPOMI L1 filename) |
| `--pace-dir` | Independent PACE SVD SIF | `interim_<YYYYMMDDTHHMMSS>_svd_retrieval_full_parallel.nc` (stamp from PACE L1B name) |

**TROPOMI HDF5** variables used (`RETRIEVAL_RESULT/`):

- `detector_pixel`, `scanline` — join keys
- `sif`, `chi2` — values written out
- (`latitude`, `longitude` — only for `--dry-run-validate`)

**PACE NetCDF** variables used:

- `source_pixel_index`, `source_scan_index` — join keys (1-based L1 indices)
- `sif_radiance_678nm`, `reduced_chi2`, `converged` — values written out
- (`latitude`, `longitude` — only for `--dry-run-validate`)

### Index convention

Matchup indices are stored 0-based relative to the L1 grids; the independent products use 1-based L1 indices. The join always adds `+1`:

| Matchup field | Product field |
|---------------|---------------|
| `trop_pix + 1` | `RETRIEVAL_RESULT/detector_pixel` |
| `trop_scan + 1` | `RETRIEVAL_RESULT/scanline` |
| `pace_pix + 1` | `source_pixel_index` |
| `pace_scan + 1` | `source_scan_index` |

Orbit / granule lookup:

- TROPOMI orbit is parsed from `tropomi_path` (`_NNNNN_` in the L1 name) → `fs_ret_NNNNN_b5_v10.h5`
- PACE time stamp is parsed from `pace_path` (`PACE_OCI.YYYYMMDDTHHMMSS`) → `interim_<stamp>_svd_retrieval_full_parallel.nc`

Hits that lack a file or a matching (pixel, scan) stay as NaN / `found=0`.

### Validate indices first (recommended)

Checks lat/lon agreement for a handful of hits under the `+1` convention, then exits (no NetCDF written):

```bash
python TROPOMI_OCI_colocation/join_independent_sif.py \
  --matchup-nc /path/to/products/matchup.nc \
  --trop-dir /home/zhe2/data/MyProjects/TROPOMI_SIF/results/b5 \
  --pace-dir /home/zhe2/data/PACE/svd_retrieval_output \
  --dry-run-validate
```

Expect small `dlat` / `dlon` (≪ 0.01°) when the convention is correct.

### Run the join

Defaults currently point at the July 2025 run; override paths for your output:

```bash
python TROPOMI_OCI_colocation/join_independent_sif.py \
  --matchup-nc /path/to/products/matchup.nc \
  --trop-dir /home/zhe2/data/MyProjects/TROPOMI_SIF/results/b5 \
  --pace-dir /home/zhe2/data/PACE/svd_retrieval_output \
  --output-dir /path/to/output_dir/matches \
  --output-name matchup_independent_sif.nc
```

Useful flags:

| Flag | Meaning |
|------|---------|
| `--max-matches N` | Join only the first `N` matches (`0` = all) |
| `--output-dir` | Directory for the result NetCDF + sidecar JSON |
| `--output-name` | Filename under `--output-dir` (default `matchup_independent_sif.nc`) |
| `--dry-run-validate` | Lat/lon check only; no write |

Smoke test:

```bash
python TROPOMI_OCI_colocation/join_independent_sif.py \
  --matchup-nc /path/to/products/matchup.nc \
  --max-matches 500
```

### Outputs

Written to `--output-dir/--output-name` (e.g. `…/matches/matchup_independent_sif.nc`):

| Variable | Contents |
|----------|----------|
| `sif_trop` | TROPOMI independent SIF |
| `chi2_trop` | TROPOMI retrieval χ² |
| `trop_found` | `1` if a TROPOMI hit was found |
| `sif_pace_678nm` | PACE SVD SIF at 678 nm (W m⁻² sr⁻¹ µm⁻¹) |
| `chi2_pace` | PACE reduced χ² |
| `converged_pace` | PACE convergence flag |
| `pace_found` | `1` if a finite PACE SIF was found |
| `lat_trop`, `lon_trop`, `lat_pace`, `lon_pace`, `trop_pix`, `trop_scan`, `pace_pix`, `pace_scan` | Copied from matchup for convenience |

A sidecar `matchup_independent_sif.json` records hit / miss counts, directories, and the index convention used.

If the output file is open elsewhere (e.g. a notebook), the script writes `*_new.nc` beside it instead of overwriting.

### Default paths (built into the CLI)

These are only conveniences for the July 2025 workspace layout — always pass flags for a new run:

```
--matchup-nc  …/TROPOMI_OCI_colocation_data/runs/july2025_d01-31/products/matchup.nc
--trop-dir    /home/zhe2/data/MyProjects/TROPOMI_SIF/results/b5
--pace-dir    /home/zhe2/data/PACE/svd_retrieval_output
--output-dir  …/TROPOMI_OCI_colocation_data/runs/july2025_d01-31/matches
```

---

## Dependencies

Python: `numpy`, `netCDF4`, `PyYAML` (and for join: `h5py`, `xarray`).

Julia (repo `Project.toml`) for `run_matchup_sif.py`.

---

## Notebooks

| Notebook / artifact | Role |
|---------------------|------|
| `notebooks/co-locate_cleanup.ipynb` | Source for `co_locate.py` logic |
| `notebooks/co-locate.ipynb` | Older pairing notebook |
| `notebooks/co-location_results_*.json` | Pairing dumps (from notebook or `co_locate.py`) |
| Other notebooks | Algorithm / scatter / comparison plots |

For production rebuilds of `matchup.nc`, use **`run_colocation.py`**, not the notebooks.
