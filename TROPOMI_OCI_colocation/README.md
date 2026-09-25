# TROPOMI ↔ PACE OCI co-location (SIF matchup)

Two-stage pipeline: granule pairing → pixel match + **lookup** of existing TROPOMI / PACE SIF.


| Stage | Goal                                       | Script               | Product                      |
| ----- | ------------------------------------------ | -------------------- | ---------------------------- |
| 1     | TROPOMI swath ↔ OCI L1B granules           | `co_locate.py`       | `co-location_results_*.json` |
| 2     | Pixel pairs + SIF from existing retrievals | `run_sif_matchup.py` | `products/matchup_data.nc`   |


Plots: `notebooks/compare_algorithm_on_oci_tropomi.py --matchup-nc …/matchup_data.nc`

Radiance simulation (SVD / RSR / simulate) and related CLIs live under `[archive/](archive/)` and are **not** part of this path.

---



## Layout

```
TROPOMI_OCI_colocation/
  co_locate.py                 # stage 1 → JSON
  run_sif_matchup.py           # stage 2 → matchup_data.nc
  configs/
    sif_matchup.smoke.toml     # smoke defaults (local TROPOMI b5 fs_ret)
    config.yaml                # archived radiance pipeline only
  colocation/                  # pair, match, sif_join, sif_matchup
  archive/                     # svd/rsr/simulate + legacy CLIs
  notebooks/                   # stage-1 JSON dumps + thin plot driver
```

Run from the **repo root** or this folder (CLIs add themselves to `sys.path`).

---



## Data roots


| Role           | Default (smoke)                                     | Layout                                                 |
| -------------- | --------------------------------------------------- | ------------------------------------------------------ |
| PACE OCI L1B   | `/kiwi-data/Data/satellite/PACE_OCI/L1B_V3`         | `YYYY/MM/DD/*.nc`                                      |
| TROPOMI BD5 L1 | `/net/squid/data1/projects/TROPOMI/ESA/L1`          | `YYYY/MM/*.nc`                                         |
| TROPOMI SIF    | `/home/zhe2/data/MyProjects/TROPOMI_SIF/results/b5` | `fs_ret_{orbit}_b5_v10.h5`                             |
| PACE SVD SIF   | `/home/zhe2/data/PACE/new_svd_retrieval_output`     | dated `YYYY/MM/DD/interim_*_svd_retrieval*.nc` or flat |


Index convention (validated vs lat/lon):

- TROPOMI: `detector_pixel == trop_pix+1`, `scanline == trop_scan+1` → `RETRIEVAL_RESULT/sif`
- PACE: `source_pixel_index == pace_pix+1`, `source_scan_index == pace_scan+1` → `sif_radiance_678nm`

---



## Stage 1 — granule JSON

```bash
python TROPOMI_OCI_colocation/co_locate.py \
  --year 2025 --month 1 --day-end 7 \
  --output TROPOMI_OCI_colocation/notebooks/co-location_results_202501_smoke.json
```

Optional: also write pipeline-style `swath_list.json` with `--pairs-dir DIR`.

---



## Stage 2 — matchup_data.nc

Uses existing retrievals only (no retrieval run, no SVD/RSR/simulate).

```bash
# Discover pairs from L1 (or set output.colocation_json in the TOML)
python TROPOMI_OCI_colocation/run_sif_matchup.py \
  --config TROPOMI_OCI_colocation/configs/sif_matchup.example.toml

# Or reuse a stage-1 JSON
python TROPOMI_OCI_colocation/run_sif_matchup.py \
  --config TROPOMI_OCI_colocation/configs/sif_matchup.example.toml \
  --colocation-json TROPOMI_OCI_colocation/notebooks/co-location_results_202501.json \
  --max-swaths 2
```

Writes under `output_dir`:

```
pairs/swath_list.json
matches/matches.npz (+ meta / paths)
products/matchup_data.nc
products/matchup_data_meta.json
```

`matchup_data.nc` includes geo, indices, paths, `sif_tropomi`, `sif_pace_678nm`, QC flags (`trop_found`, `pace_found`), Δt / distance. Matches are cached by fingerprint; pass `--force-matches` to rebuild.

### Plots

```bash
python TROPOMI_OCI_colocation/notebooks/compare_algorithm_on_oci_tropomi.py \
  --matchup-nc /path/to/products/matchup_data.nc
```

---



## Archive (legacy / radiance)

See `[archive/README.md](archive/README.md)`. Examples:

```bash
python TROPOMI_OCI_colocation/archive/run_colocation.py \
  --config TROPOMI_OCI_colocation/configs/config.yaml
python TROPOMI_OCI_colocation/archive/join_independent_sif.py ...
python TROPOMI_OCI_colocation/archive/run_matchup_sif.py \
  --config TROPOMI_OCI_colocation/configs/config.sif_retrieval.toml
```

Row-aware TROPOMI PC→678 compare (pre-cleanup) is frozen as `archive/compare_algorithm_row_aware.py`.