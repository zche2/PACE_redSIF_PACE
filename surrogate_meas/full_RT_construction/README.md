# full_RT_construction

Build surrogate PACE OCI TOA spectra over open ocean with full vSmartMOM RT
(Cox-Munk + whitecaps, O2/H2O absorption, GCHP aerosols, water-leaving SIF),
then test the SVD retrieval on them.

## Pipeline at a glance

| Step | Script | Output |
|------|--------|--------|
| 1. Aerosol columns (GCHP/TOMAS) | `build_gchp_ocean_columns.py` (or `GCHP_aerosol.py`) | `output_aerosol_profiles/gchp_ocean_columns_n500.nc` |
| 2. Add GCHP T/q/p_half to step 1 (optional) | `augment_gchp_columns_met.py` | same file, updated in place |
| 3. Atmosphere columns (MERRA-2 sea T/p/q) | `build_sea_surface_Tpq_columns.py` | `merra2_sea_Tpq_columns_n1500.nc` |
| 4. Single test spectrum | `test_run.jl` | `output/<tagged>.nc` |
| 5. RT ensemble | `generate_test_spec.jl` (via `generate_test_spec.sh`) | `output/rt_toa_ensemble.nc` |
| 6. SVD retrieval on the ensemble | `retrieve_rt_ensemble.jl` | `output/rt_retrieval_<tag>.{nc,png}` |

Steps 1 and 3 only need to be run once; both files already exist here.

## Environments

- **Python** (steps 1–3): needs `numpy`, `scipy`, `netCDF4`, `matplotlib`,
  `cartopy`, `shapely`. Run from this directory so the local imports resolve.
- **Julia RT** (steps 4–5): `test_run.jl` and `generate_test_spec.jl` activate
  `/home/zhe2/FraLab/vSmartMOM.jl` themselves.
- **Julia retrieval** (step 6): uses this repo's project; run from the repo
  root with `--project=.`.

Shared inputs (hard-coded paths):

- RT config: `../configs/ocean_coxmunk_0912.yaml`
- OCI RSRs: `/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/PACE_OCI_RSRs.nc`
- SIF shapes: `.../Files_in_use/SIF_singular_vector.jld2`
- OCI SNR model: `.../Files_in_use/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt`

## 1. Aerosol columns

```bash
cd surrogate_meas/full_RT_construction
python build_gchp_ocean_columns.py          # 500 open-ocean columns
# or, with sea-salt maps/profiles as well:
python GCHP_aerosol.py
```

Reads GCHP output `/kiwi-data/Data/model/GeosChem/GEOSChem.Custom.20190702_0000z.nc4`.
Computes wet radii (GEOS-Chem TOMAS GETDP) and layer AOD at 550 nm (Mie),
split into the YAML species `so4, sala, salc, ocpi, bcpi, dust1..7`.
Levels are stored **BOA → TOA**.

- `N_SAMPLES` (default 500). Don't smoke-test with a small `N_SAMPLES`, because it
  overwrites the same `gchp_ocean_columns_n500.nc`.

Optional: `python augment_gchp_columns_met.py` adds GCHP `Met_T`, `Met_SPHU`
and `p_half` to that file. The ensemble does not use GCHP T/q (it uses MERRA);
`p_half` is used for the vertical AOD remap when present.

## 2. Atmosphere columns

```bash
python build_sea_surface_Tpq_columns.py
```

Random MERRA-2 sea-surface profiles with |lat| ≤ 60°, reduced from 72 to about 24
layers, stored **TOA → BOA**.

| Env var | Default |
|---------|---------|
| `MERRA_PATH` | `/home/zhe2/data/MERRA2_reanalysis/MERRA2_400.inst6_3d_ana_Nv.20240705.nc4` |
| `N_SAMPLES` | 1500 |
| `SEED` | 20260919 |
| `PROFILE_STRIDE` | 3 (keep in sync with `generate_test_spec.jl`) |
| `OUT_NC` | `merra2_sea_Tpq_columns_n<N_SAMPLES>.nc` |

## 3. Single test spectrum

```bash
julia surrogate_meas/full_RT_construction/test_run.jl
```

Edit the constants at the top of the script: `SIF_PEAK`, `SIF_LIBRARY_INDEX`,
`ENABLE_AEROSOLS`, `COLUMN_INDEX`, `INCLUDE_WHITECAPS`, `WHITECAP_ALBEDO`,
`WIND_SPEED`. Writes one NetCDF to `output/`, with the settings in its filename.

## 4. RT ensemble

Recommended: the GPU wrapper, which waits for a GPU with ≥ 12 GiB free, then
runs 1000 samples in the background:

```bash
bash surrogate_meas/full_RT_construction/generate_test_spec.sh
# prints PID, log path and OUT_NC
tail -f surrogate_meas/full_RT_construction/output/generate_test_spec_seaTpq_n1000_<stamp>.log
```

Direct runs:

```bash
cd surrogate_meas/full_RT_construction
julia --project=/home/zhe2/FraLab/vSmartMOM.jl generate_test_spec.jl
N_SAMPLES=2 julia --project=/home/zhe2/FraLab/vSmartMOM.jl generate_test_spec.jl       # smoke test
ENABLE_AEROSOLS=false N_SAMPLES=2 julia --project=/home/zhe2/FraLab/vSmartMOM.jl generate_test_spec.jl
ARCH=CPU julia --project=/home/zhe2/FraLab/vSmartMOM.jl generate_test_spec.jl          # force CPU
```

For each sample it draws:
- a MERRA T/p/q profile
- a GCHP aerosol column, remapped onto the RT pressure grid while keeping the column AOD
- a SIF shape, scaled to 0–0.5 W m⁻² sr⁻¹ µm⁻¹ at 678 nm
- geometry: SZA 5–70°, VZA 0–60°, relative azimuth 0–180°
- surface: wind speed 0–10 m/s, whitecap albedo 0.1–0.5, whitecaps on or off

It then runs the RT, convolves to OCI bands and adds SNR-model noise.
The ranges are set by the `*_RANGE` constants in the script.

| Env var | Default | Meaning |
|---------|---------|---------|
| `N_SAMPLES` | 500 | ensemble size (the `.sh` wrapper uses 1000) |
| `OUT_NC` | `output/rt_toa_ensemble.nc` | output file |
| `ENABLE_AEROSOLS` | `true` | GCHP aerosols on or off |
| `ENSEMBLE_SEED` | 20260913 | RNG seed |
| `MERRA_SEA_TPQ_NC` | `merra2_sea_Tpq_columns_n1500.nc` | atmosphere columns |
| `PROFILE_STRIDE` | 3 | layer reduction |
| `ARCH` | auto | `GPU` or `CPU` |

Output variables in `rt_toa_ensemble.nc`, per band × sample:
- `radiance_clean`, `radiance_noisy`, `sigma_noise`, `sif_waterleaving`

Per sample:
- `sif_678`, `sza`, `vza`, `vaz`, `wind_speed`, `whitecap_albedo`,
  `include_whitecaps`, `column_index`, `profile_index`
- the layer profiles `p`, `T`, `q`

Notes:
- **Resume.** Results are saved after every sample. Re-running with the same
  `OUT_NC`, seed and `N_SAMPLES` fills in only the unfinished samples.
- **Runtime.** About 2 min per sample on an A100 with `polarization_type: Stokes_I()`
  in the YAML. Stokes_IQUV is about 2× slower and can run out of memory on a
  40 GB GPU. CPU runs take weeks for 1000 samples.
- **Reading while it runs.** `ncdump` on the live file can fail with an HDF
  error while the generator has it open. Use `retrieve_rt_ensemble.jl`, which
  reads a snapshot copy, or wait until the run finishes.
- **Stopping it.** `kill <PID>` (the PID is printed by the `.sh` wrapper, or find
  it with `pgrep -af generate_test_spec`).

## 5. SVD retrieval on the ensemble

```bash
cd /home/zhe2/FraLab/PACE_redSIF_PACE
julia --project=. surrogate_meas/full_RT_construction/retrieve_rt_ensemble.jl
julia --project=. surrogate_meas/full_RT_construction/retrieve_rt_ensemble.jl \
    surrogate_meas/configs/svd_nPC15_npoly3.toml
```

- Config: first CLI argument, else `SVD_TOML`, else `../configs/svd_nPC15_npoly5.toml`.
  `n_pc` / `n_legendre` come from `[fit.svd]`.
- `RT_NC` (default `output/rt_toa_ensemble.nc`). The script copies it to
  `output/rt_toa_ensemble_snap.nc` first, so it is safe to run while the
  generator is still writing.
- Fits each finished sample at its own SZA. Outputs are tagged
  `nPC<k>_npoly<n>`:
  - `output/rt_retrieval_<tag>.nc`
  - `output/rt_retrieval_toa_sif_<tag>.png`
  - `output/rt_retrieval_sif678_<tag>.png`

## Other files

- `ocean_column_aerosols.jl`: library used by `test_run.jl` and
  `generate_test_spec.jl`. It loads a GCHP column into vSmartMOM aerosols and
  remaps its layer AOD onto the RT pressure grid. Pressure and AOD are always
  reversed together (GCHP BOA→TOA vs RT TOA→BOA), and the remap keeps the column AOD.
- `aerosol_tutorial_copy.jl`: vSmartMOM Mie tutorial (reference only).
- Notebooks (exploration): `GCHP_aerosol.ipynb` (source of `GCHP_aerosol.py`),
  `gchp_vs_geoschem.ipynb`, `seasalt_profiles.ipynb`, `TOA_spectra_configs.ipynb`.
- Older output folders: `output_no_CoxMunk_sad/`, `output_test_realRT_noAerosol/`.
