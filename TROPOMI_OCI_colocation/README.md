# TROPOMI ↔ PACE OCI co-location (SIF matchup)

Find TROPOMI BD5 and PACE OCI pixels that look at the same place at nearly the same time, then attach the SIF each instrument's **existing** retrieval produced there.

| Stage | What it does                                                     | Script               | Output                               |
| ----- | ---------------------------------------------------------------- | -------------------- | ------------------------------------ |
| 1     | Pair TROPOMI orbit files with OCI L1B granules (time + bbox)      | `co_locate.py`       | `co-location_results_*.json`         |
| 2     | Match pixels (distance + Δt + dark ocean), look up both SIFs      | `run_sif_matchup.py` | `<output_dir>/products/matchup_data.nc` |
| plot  | Scatter / comparison figures                                     | `notebooks/compare_algorithm_on_oci_tropomi.py` | figures |

No retrieval is run here. Radiance simulation (SVD / RSR / simulate) is in `archive/` (local only, not tracked) and is not part of this path.

---

## Quick start (one month)

Run everything from the **repo root** (`PACE_redSIF_PACE/`).

**1. Write a config.** `configs/` is gitignored, so create a new file such as `TROPOMI_OCI_colocation/configs/sif_matchup.jan2025.toml` in any editor, paste the [template below](#config-template), and edit `[window]` and `[output]`.

**2. Stage 1 — granule pairing** (~6 min per week of data):

```bash
python TROPOMI_OCI_colocation/co_locate.py \
  --year 2025 --month 1 --day-end 31 --dt-max-min 60 \
  --output TROPOMI_OCI_colocation/notebooks/co-location_results_202501_dt60.json
```

Use the same `--dt-max-min` as `[match] max_dt_min` in the config.

**3. Stage 2 — pixel match + SIF lookup** (hours; run in the background):

```bash
nohup python -u TROPOMI_OCI_colocation/run_sif_matchup.py \
  --config TROPOMI_OCI_colocation/configs/sif_matchup.jan2025.toml \
  > /path/to/runs/jan2025.log 2>&1 &
tail -f /path/to/runs/jan2025.log
```

`output.colocation_json` in the config points stage 2 at the stage-1 JSON. Always pass `--config` — the built-in default (`configs/sif_matchup.smoke.toml`) usually does not exist.

**4. Plot:**

```bash
python TROPOMI_OCI_colocation/notebooks/compare_algorithm_on_oci_tropomi.py \
  --matchup-nc /path/to/runs/<run>/products/matchup_data.nc \
  --fig-dir /path/to/figures
```

**Smoke test** before a full month: `--day-end 3` in stage 1, and `--max-swaths 2` in stage 2.

---

## Config template

Stage 2 reads a TOML file. Every key is optional (defaults shown); in practice set `[window]`, `[retrieval]` and `[output]`.

```toml
[window]                 # days of PACE granules to use
year = 2025
month = 1
day_start = 1
day_end = 31

[l1]                     # L1 inputs + granule pairing (only used when no colocation_json is given)
pace_dir = "/kiwi-data/Data/satellite/PACE_OCI/L1B_V3"
tropomi_dir = "/net/squid/data1/projects/TROPOMI/ESA/L1"
tropomi_product = "OFFL"          # OFFL | RPRO | BOTH
pace_stride = 1                   # use every Nth PACE granule
pace_duration_min = 5             # nominal OCI granule length
dt_max_min_granule = 60.0         # granule time gap allowed (never tighter than max_dt_min)
bbox_margin_deg = 0.5             # pad on granule / TROPOMI-chunk bounding boxes
tropo_bbox_scan_stride = 40       # subsampling used only to build TROPOMI bboxes
tropo_bbox_pix_stride = 20
tropo_chunk_scans = 400           # TROPOMI scans per bbox chunk

[match]                  # pixel-level science cuts
max_dist_km = 2.5                 # PACE pixel ↔ nearest TROPOMI pixel centre
max_dt_min = 60.0                 # |t_TROPOMI − t_PACE| (minutes)
lt_max_oci = 30.0                 # "dark" OCI pixel: max L_TOA (W m-2 sr-1 µm-1) over bands ≥ lt_mask_wl_min_oci
lt_mask_wl_min_oci = 600.0
pace_stride_pix = 150             # PACE subsampling, applied to BOTH scans and pixels
tropo_scan_stride = 1             # TROPOMI subsampling for the KD-tree (1 = full grid)
tropo_pix_stride = 1

[retrieval]              # existing SIF products to look up
tropomi_ret_dir = "/home/zhe2/data/MyProjects/TROPOMI_SIF/results/b5"   # fs_ret_<orbit>_b5_v10.h5
pace_ret_dir = "/home/zhe2/data/PACE/new_svd_retrieval_output"          # [YYYY/MM/DD/]interim_<stamp>_svd_retrieval*.nc

[output]
output_dir = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/TROPOMI_OCI_colocation_data/runs/jan2025"
colocation_json = "/home/zhe2/FraLab/PACE_redSIF_PACE/TROPOMI_OCI_colocation/notebooks/co-location_results_202501_dt60.json"
# max_swaths = 2                  # smoke test: first N TROPOMI swaths only
```

The code defaults (used when a key is missing) are the older, sparser settings: `max_dt_min = 30`, `tropo_*_stride = 2`. The template values above are the recommended ones — see [Choosing settings](#choosing-settings).

---

## Stage 1 — `co_locate.py`

For each PACE L1B granule in the window, keep every TROPOMI BD5 orbit file that

1. is **close in time**: the gap between the two files' time spans (0 if they overlap) is ≤ `--dt-max-min`, and
2. **overlaps in space**: the padded OCI granule bbox intersects at least one padded bbox of a 400-scan TROPOMI chunk.

The output groups PACE granules by TROPOMI swath:

```json
{"meta": {"year": 2025, "month": 1, "dt_max_min": 60.0, "pace_dir": "...", "n_swaths": 409, ...},
 "swaths": [{"orbit": 37407, "tropomi_path": "...", "n_pace_matches": 12, "pace_granules": ["PACE_OCI....nc", ...]}, ...]}
```

| Flag | Default | Meaning |
| ---- | ------- | ------- |
| `--year`, `--month`, `--day-start`, `--day-end` | 2025, 7, 1, 31 | PACE window |
| `--dt-max-min` | 60 | granule time-gap limit (min) |
| `--bbox-margin-deg` | 0.5 | bbox padding |
| `--pace-dir`, `--tropomi-dir` | kiwi / squid roots | L1 roots |
| `--tropomi-product` | `OFFL` | `OFFL`, `RPRO` or `BOTH` |
| `--pace-stride` | 1 | every Nth PACE granule |
| `--output`, `-o` | `notebooks/co-location_results_YYYYMM.json` | JSON path |
| `--pairs-dir` | – | also write a raw `swath_list.json` there |
| `--print-swaths` | off | list swaths to stdout |

Stage 1 is only a prefilter. It must be at least as loose as the pixel Δt cut in stage 2; otherwise granules are dropped before pixel matching can see them.

---

## Stage 2 — `run_sif_matchup.py`

For every swath in the stage-1 JSON:

1. Read the TROPOMI geolocation on the `tropo_*_stride` grid and build a KD-tree.
2. For each paired PACE granule, take pixels on the `pace_stride_pix` grid that are **dark** (max L_TOA over bands ≥ 600 nm below `lt_max_oci`: clear ocean, no clouds or land).
3. Find each PACE pixel's nearest TROPOMI pixel; keep it if the distance ≤ `max_dist_km` and |Δt| ≤ `max_dt_min` (scan-line times).
4. Look up SIF at the matched indices in the TROPOMI `fs_ret` file of that orbit and the PACE SVD file of that granule.

| Flag | Meaning |
| ---- | ------- |
| `--config PATH` | TOML config (always pass it) |
| `--colocation-json PATH` | stage-1 JSON; overrides `output.colocation_json` |
| `--max-swaths N` | only the first N swaths (smoke test) |
| `--force-matches` | recompute pixel matches even if a cached result looks valid |
| `--force-pairs` | rebuild pairs from L1 (only without a JSON) and recompute matches |

Without any stage-1 JSON, stage 2 runs the pairing itself from `[l1]`.

### Outputs (under `output_dir`)

```
pairs/swath_list.json, pairs_meta.json
matches/matches.npz, match_paths.json, matches_meta.json     # cached pixel matches
products/matchup_data.nc                                      # the product
products/matchup_data_meta.json                               # hit/miss counts, timings
run_meta.json                                                 # full config of the last run
```

`matchup_data.nc` has one row per matched PACE pixel (dimension `match`):

| Variables | Meaning |
| --------- | ------- |
| `pace_lat/lon`, `trop_lat/lon` | pixel centres |
| `dist_km`, `dt_min` | separation; `dt_min = t_TROPOMI − t_PACE` (positive = TROPOMI later) |
| `pace_sza/vza`, `trop_sza/vza` | geometry |
| `pace_scan/pix`, `trop_scan/pix`, `swath_id` | 0-based L1 indices |
| `sif_tropomi`, `chi2_tropomi`, `trop_found` | TROPOMI `fs_ret` lookup |
| `sif_pace_678nm`, `chi2_pace`, `converged_pace`, `pace_found` | PACE SVD lookup (W m-2 sr-1 µm-1) |
| `tropomi_path`, `pace_path` | L1 source files |

Filter on `trop_found == 1 & pace_found == 1` before comparing. `fs_ret` files only hold the pixels the TROPOMI retrieval processed, so roughly half of the matches have no TROPOMI SIF (Jan 2025: 2,968 of 6,318). That is expected, not a matching error.

### Caching

Pixel matching is the slow step and is cached in `matches/`. It is reused (log: `[matches] CACHE HIT`) when the fingerprint of the stage-1 pairs (each swath's TROPOMI file + PACE granule list) and the `[match]` settings is unchanged. Changing a `[match]` value or regenerating the JSON with different pairs triggers a rebuild automatically. Use `--force-matches` when inputs changed in a way the fingerprint cannot see (e.g. L1 files replaced in place). A new `output_dir` always starts clean.

The SIF lookup (step 4) always reruns, so a newer PACE or TROPOMI retrieval only needs a stage-2 rerun, not new matching.

### Index convention

Matchup indices are 0-based; both retrieval products are 1-based (validated against lat/lon):

- TROPOMI: `RETRIEVAL_RESULT/detector_pixel == trop_pix + 1`, `scanline == trop_scan + 1`
- PACE: `source_pixel_index == pace_pix + 1`, `source_scan_index == pace_scan + 1`

---

## Choosing settings

PACE crosses the equator at ~13:00 and TROPOMI at ~13:30 local time, so most coincidences have TROPOMI **later** by +20 to +80 min. The distribution peaks at +30 to +50 min, so a ±30 min window cuts it at the peak.

Pixel matches for Jan 4–10 2025 (`pace_stride_pix = 150`, 2.5 km, dark ocean):

| TROPOMI stride | ±30 min | ±45 min | ±60 min | ±90 min |
| -------------- | ------: | ------: | ------: | ------: |
| 2 (code default) | 1,419 | 2,183 | 2,760 | 3,611 |
| 1 (recommended)  | 5,327 | 8,203 | 10,417 | 13,743 |

- **`max_dt_min`**: ±60 min roughly doubles the matches of ±30. Wider windows add more but pair scenes further apart in time.
- **`tropo_*_stride = 1`**: 3.75× more matches at the same window. With stride 2, TROPOMI pixel centres are ~7–11 km apart, so most PACE pixels have no kept centre within 2.5 km even when they sit inside a TROPOMI footprint.
- **`pace_stride_pix`**: 150 keeps ~90 OCI pixels per granule (every 150th scan and pixel). Smaller strides add many more rows, but several OCI pixels then fall in the same TROPOMI pixel (TROPOMI ≈ 5.5 × 3.5 km, OCI ≈ 1.2 km); de-duplicate or average per `(swath_id, trop_scan, trop_pix)` before statistics.

**Runtime.** Stage 1: ~6 min per week. Stage 2 with stride 1: ~20–40 s per swath (~2.5–4.5 h for a month of ~400 swaths), dominated by strided reads of the OCI red bands.

---

## Data roots

| Role | Default | Layout |
| ---- | ------- | ------ |
| PACE OCI L1B | `/kiwi-data/Data/satellite/PACE_OCI/L1B_V3` | `YYYY/MM/DD/PACE_OCI.*.L1B.V3.nc` |
| TROPOMI BD5 L1 | `/net/squid/data1/projects/TROPOMI/ESA/L1` | `YYYY/MM/S5P_*_L1B_RA_BD5_*.nc` (folder = start date) |
| TROPOMI SIF | `/home/zhe2/data/MyProjects/TROPOMI_SIF/results/b5` | `fs_ret_<orbit>_b5_v10.h5` |
| PACE SVD SIF | `/home/zhe2/data/PACE/new_svd_retrieval_output` | `YYYY/MM/DD/interim_<stamp>_svd_retrieval*.nc` or flat |

When several PACE retrieval files exist for one granule, names containing `full_parallel` win; otherwise the first match in sorted order is used.

---

## Troubleshooting

- **`FileNotFoundError` on the config** — the path in `--config` does not exist; `configs/` is untracked, so check the file on this machine.
- **Wider window but no new matches** — confirm stage 1 was rerun with the wider `--dt-max-min` and that stage 2 points at the new JSON. If the log says `CACHE HIT`, the pairs and settings match an earlier run.
- **`n_file_miss` in the log / meta** — no `fs_ret` file for that orbit, or no PACE retrieval file for that granule, under the `[retrieval]` dirs.
- **Reusing a log file** — `>` truncates it, but text from a process that still has it open can remain; start a new log per run.

---

## Layout

```
TROPOMI_OCI_colocation/
  co_locate.py            # stage 1 → JSON
  run_sif_matchup.py      # stage 2 → matchup_data.nc
  colocation/
    discover.py           # list L1 files by name/time
    geo.py                # bboxes, time windows, lon/lat → xyz
    pair.py               # stage-1 pairing
    match.py              # pixel matching (cached)
    sif_join.py           # fs_ret / PACE SVD lookups
    sif_matchup.py        # stage-2 driver + NetCDF writer
    cache.py, config.py
  notebooks/              # stage-1 JSONs + plot driver
  sif_retrieval/          # Julia matchup retrieval (archived radiance path)
  configs/                # local only (gitignored)
  archive/                # local only (gitignored): radiance simulation + legacy CLIs
```

## Archive (legacy / radiance path)

Local only; see `archive/README.md` on the machine that has it:

```bash
python TROPOMI_OCI_colocation/archive/run_colocation.py --config TROPOMI_OCI_colocation/configs/config.yaml
python TROPOMI_OCI_colocation/archive/run_matchup_sif.py --config TROPOMI_OCI_colocation/configs/config.sif_retrieval.toml
```
