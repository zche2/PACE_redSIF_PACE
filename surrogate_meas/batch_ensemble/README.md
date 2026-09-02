# Surrogate ensemble (5,000 measurements + SVD retrieval)

## Layout


| Path                               | Description                                            |
| ---------------------------------- | ------------------------------------------------------ |
| `../configs/svd_nPC15_npoly5.toml` | Independent SVD config (**n_pc=15**, **n_legendre=5**) |
| `../run_ensemble.jl`               | End-to-end generator + retrieval                       |
| `truth_ensemble.nc`                | Ground-truth TOA, SIF, T1, T2, decay, …                |
| `retrieval_ensemble.nc`            | SVD state, fitted SIF/TOA, QA metrics                  |
| `plots/`                           | Radiance / SIF / retrieval comparison figures          |
| `run.log`                          | Full run log                                           |


## Run

```bash
# full ensemble (default N_SAMPLES=5000); decay ∼ U[4,12] per sample; OUT_SUFFIX=_new
julia --project=. -t 8 surrogate_meas/run_ensemble.jl

# npoly3 / npoly5 (regenerate truth once, then reuse)
# Truth and retrieval both use ENSEMBLE_CONFIG (alias SVD_CONFIG); [spectral] drives λ window.
REUSE_TRUTH=0 SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly3.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl
REUSE_TRUTH=1 SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly5.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl

# O₂-A window (640–770 nm): edit lambda_max_nm in svd_nPC15_npoly3_O2A.toml, then:
REUSE_TRUTH=0 OUT_SUFFIX=_O2A770 SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly3_O2A.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl

# Degraded measurement σ in 683–695 nm (σ×100 on those bands):
REUSE_TRUTH=0 OUT_SUFFIX=_degrade683 SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly3_degrade683.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl
# Or ENV: DEGRADE_NOISE=1 DEGRADE_SIGMA_FACTOR=100 ...

# zero-SIF floor (separate truth file)
ZERO_SIF=1 REUSE_TRUTH=0 SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly3.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl

# fixed SIF shape (EV1 only; sif_678 = water-leaving, not SIF×T₁)
FIXED_SIF_EV1=1 REUSE_TRUTH=0 SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly3.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl
REUSE_TRUTH=1 FIXED_SIF_EV1=1 SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly5.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl
# or: bash surrogate_meas/batch_ensemble/run_fixed_sif_ev1.sh

# smoke test
N_SAMPLES=20 N_ATM_POOL=5 N_CONT_POOL=15 julia --project=. -t 4 surrogate_meas/run_ensemble.jl

# legacy output names (no _new suffix)
OUT_SUFFIX= julia --project=. -t 8 surrogate_meas/run_ensemble.jl
```

Outputs: `truth_ensemble_new.nc`, `retrieval_ensemble_npoly3_new.nc`, `plots_npoly3_new/`, …

Fixed-SIF EV1: `truth_ensemble_fixedSIF_new.nc`, `retrieval_ensemble_npoly3_fixedSIF_new.nc`, `plots_npoly3_fixedSIF_new/`. Retrieval NetCDF includes `sif_wl_ret`, `SIF_true`, and `sif_678_*` in water-leaving units (`sif_678_metric=water_leaving`).

Retrieval NetCDF also stores `averaging_kernel(state_ret, state_true, sample)` (= `S_post · H_obs`; row = retrieved, col = true), `ak_trace` (= `tr(A)`), and `state_names_csv` attribute. In Python/xarray, transpose the 2D slice (`.T`) to match Julia row/column semantics.

## Truth variables (`truth_ensemble.nc`)

- `R_toa_clean`, `R_toa_noisy`, `R_cont`, `SIF`, `R_sif_toa`
- `T1` (atm one-way), `T2` (atm×solar two-way, decay-weighted)
- `decay`, `sif_library_index`, `sif_strength`, `profile_index`



## Retrieval config check

Confirm in `svd_nPC15_npoly5.toml`:

```toml
[fit.svd]
n_pc       = 15
n_legendre = 3
```

