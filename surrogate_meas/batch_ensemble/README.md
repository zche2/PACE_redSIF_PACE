# Surrogate ensemble (5,000 measurements + SVD retrieval)

## Layout

| Path | Description |
|------|-------------|
| `../configs/svd_nPC15_npoly5.toml` | Independent SVD config (**n_pc=15**, **n_legendre=5**) |
| `../run_ensemble.jl` | End-to-end generator + retrieval |
| `truth_ensemble.nc` | Ground-truth TOA, SIF, T1, T2, decay, … |
| `retrieval_ensemble.nc` | SVD state, fitted SIF/TOA, QA metrics |
| `plots/` | Radiance / SIF / retrieval comparison figures |
| `run.log` | Full run log |

## Run

```bash
# full ensemble (default N_SAMPLES=5000)
julia --project=. -t 8 surrogate_meas/run_ensemble.jl

# smoke test
N_SAMPLES=20 N_ATM_POOL=5 N_CONT_POOL=15 julia --project=. -t 4 surrogate_meas/run_ensemble.jl
```

## Truth variables (`truth_ensemble.nc`)

- `R_toa_clean`, `R_toa_noisy`, `R_cont`, `SIF`, `R_sif_toa`
- `T1` (atm one-way), `T2` (atm×solar two-way, decay-weighted)
- `decay`, `sif_library_index`, `sif_strength`, `profile_index`

## Retrieval config check

Confirm in `svd_nPC15_npoly5.toml`:

```toml
[fit.svd]
n_pc       = 15
n_legendre = 5
```
