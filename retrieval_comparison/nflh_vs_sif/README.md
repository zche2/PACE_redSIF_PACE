# nFLH vs SVD SIF, per sounding, for single days

Compares the PACE L2 nFLH (V3.2 AOP) with the SVD SIF retrieval (default run
`new_svd_retrieval_output_O2A`, nPC15 npoly3, 640–780 nm) for every sounding of a day, with a filter that can
be changed without re-reading the granules.

## 1. Extract (once per date)

```bash
python retrieval_comparison/nflh_vs_sif/extract_soundings.py --dates 20250115 20250415 --workers 24
```

Reads each retrieval granule (`<run-dir>/YYYY/MM/DD/interim_<granule>_svd_retrieval.nc`) and the matching
L2 AOP / BGC V3.2 and L1B files on the same grid. Keeps every sounding that has an L2 nflh (no other cut) and
writes `/home/zhe2/data/PACE/nflh_vs_sif_soundings/<run>/soundings_<date>.parquet` (~2 GB per day, outside
the repository). Granules without AOP files (polar, no ocean-colour products) are skipped.

Columns: `granule scan pixel lat lon sza vza sif_678 sif_678_sigma converged status_code n_steps
reduced_chi2 dof is_dark is_ocean alpha l2_flags nflh chlor_a aot_865 rtoa_678`.

- `sif_678`: water-leaving SIF at the band nearest 678.2 nm (`sif_radiance_678nm`); `sif_678_sigma` its
  posterior 1σ.
- `nflh`: stored L2 value (680 nm height, 664.6–667.0 | 710.6–713.1 nm baseline, clipped at 0).
- `rtoa_678`: L1B TOA radiance at 678.18 nm, rhot_red·E0·cos(SZA)/π/esd, as the retrieval sees it.

## 2. Plot with a filter

```bash
python retrieval_comparison/nflh_vs_sif/plot_nflh_vs_sif.py                      # l2_flags == 0, converged
python retrieval_comparison/nflh_vs_sif/plot_nflh_vs_sif.py \
    --query "converged == 1 and is_dark == 1 and reduced_chi2 < 5" --tag l2zero_dark_chi2lt5
python retrieval_comparison/nflh_vs_sif/plot_nflh_vs_sif.py --l2 LAND,COASTZ,HISOLZEN \
    --query "is_dark == 1 and status_code in [1, 5] and abs(sif_678) <= 5.0 and reduced_chi2 <= 1.5" \
    --tag 3flags_dark_status15_chi2le1.5
```

| option | meaning |
|---|---|
| `--l2 zero` | `l2_flags == 0`, no flag set (default) |
| `--l2 any` | no L2 flag cut |
| `--l2 NAME,...` | only these flags must be off |
| `--query "..."` | pandas query on the columns (default `converged == 1`; `""` for none) |
| `--rel-min` | relative ΔSIF panel only where nFLH > this (default 0.05) |
| `--dates`, `--tag` | dates, output folder name |

Outputs in `out/<tag>/`:
- `maps.png`: every sounding (no gridding) of SIF, nFLH and ΔSIF = SIF − nFLH, plus relative ΔSIF =
  (SIF − nFLH)/nFLH against nFLH (where nFLH > `--rel-min`), per date.
- `scatter.png`: density plot, histogram of the difference, and the difference against nFLH and latitude.
- `summary.txt`: how many soundings survive each filter step, plus statistics.

## 3. ΔSIF vs air mass by chlorophyll

```bash
python retrieval_comparison/nflh_vs_sif/plot_dsif_airmass.py --l2 LAND,COASTZ,HISOLZEN \
    --query "is_dark == 1 and status_code in [1, 5] and abs(sif_678) <= 5.0 and reduced_chi2 <= 1.5" \
    --tag 3flags_dark_status15_chi2le1.5
```

Single-day, per-sounding version of `notebooks/figures/Figure2/FigureS2_low_chl_vza_bin_month.py`. It uses
the same filter options (shared in `common.py`), Chl-a classes (`--chl-classes`, default 0.02–0.4 mg m⁻³)
and geometric air-mass bins 1/cos SZA + 1/cos VZA (`--amf-edges`, default 2–6 in steps of 0.1).
For each date it plots the median nFLH, the median SIF, and the median and IQR of ΔSIF, against air
mass, plus a Chl-a × air-mass heatmap of ΔSIF. Bins with fewer than `--min-per-cell` soundings
(default 500) are left blank.
Output: `out/<tag>/dsif_airmass_chl.png` and `.csv`.

## 4. Actual vs normalized SIF

```bash
python retrieval_comparison/nflh_vs_sif/plot_sif_normalized.py --l2 LAND,COASTZ,HISOLZEN \
    --query "is_dark == 1 and status_code in [1, 5] and abs(sif_678) <= 5.0 and reduced_chi2 <= 1.5" \
    --tag 3flags_dark_status15_chi2le1.5
```

Normalizes SIF per sounding as SIF · esd / cos(SZA), with esd = (d/AU)² taken from each granule's L1B
`earth_sun_distance_correction`. That puts SIF on the same footing as nFLH, which is normalized to
an overhead Sun at 1 AU (see `document/nflh/`). The downward transmittance t↓ is ignored. For each
date it plots densities of actual and normalized SIF against nFLH, and median ΔSIF against air mass by
Chl-a class, for both. Output: `out/<tag>/sif_normalized.png`. No data files are written.

## Notes

- `l2_flags == 0` also removes `PRODWARN`. That flag is set on about 99 % of the soundings whose nFLH
  would be ≤ 0 (11–14 % of the soundings with L2 nflh). The default filter therefore keeps only positive
  nFLH. To include those soundings, list the flags explicitly without `PRODWARN` (they enter at the
  clipped value, ≈ 0).
- nFLH is normalized (normalized water-leaving radiance); the retrieved SIF is not. SIF/nFLH therefore
  depends on the solar zenith angle (column `sza`).
