# PACE–TROPOMI L1B coincidence toolbox

Match **PACE OCI** L1B granules (single UTC time from the filename) with **TROPOMI L1B RA** granules when  
`T_start ≤ t_pace ≤ T_end` (see `MatchCoincidence.jl`). Inputs are local paths (directories and/or manifests); optional remote steps download into cache directories first.

Run everything from the **repository root** with the project environment:

```bash
cd /path/to/PACE_redSIF_PACE
julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl [arguments...]
```

### How to see available arguments

```bash
julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl --help
```

Same information lives in this README, in the comment header at the top of `run_match.jl`, and in `config.example.toml` for TOML keys.

## CLI modes

You can combine flags with a config path as shown below. **`--self-test` is special**: if it appears anywhere in the arguments, only tests run and the process exits (no config required).

| Mode | Command | Notes |
|------|---------|--------|
| **Self-test** | `julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl --self-test` | Runs `MatchCoincidence` + `RemoteFetch` unit checks. Exit `0` if both print `OK`. No CSV/output files. |
| **Match (default config path)** | `julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl` | Uses `toolbox/pace_tropomi_coincidence/config.toml` if present, else errors. |
| **Match (explicit config)** | `julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl /path/to/config.toml` | Reads `[paths]`, `[match]`, `[output]`, optional `[remote.*]`. Writes `csv_path` and optional `jsonl_path`. |
| **Match via env** | `MATCH_CONFIG=/path/to/config.toml julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl` | Same as passing the config path; env wins over positional config when set (non-empty). |
| **Remote fetch + match** | `julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl --fetch-remote /path/to/config.toml` | If `[remote.tropomi]` / `[remote.pace]` have `enabled = true`, downloads into their `cache_dir` values, appends those dirs to the match search paths, then runs the matcher. Requires Earthdata auth / `earthaccess` for PACE as documented in `config.example.toml`. |

Examples:

```bash
# Quick sanity check (no config file needed)
julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl --self-test

# Normal run
julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl my_coincidence.toml

# Download (per config) then match
julia --project=. toolbox/pace_tropomi_coincidence/run_match.jl --fetch-remote my_coincidence.toml
```

## Configuration

Copy [`config.example.toml`](config.example.toml) and edit `[paths]`, `[match]` (including optional `utc_window_start` / `utc_window_end` / `utc_day_stride`), `[output]`, and optionally `[remote.tropomi]` / `[remote.pace]`. **`[match] utc_window_*` is the shared default for remote fetch**: if `[remote.tropomi]` omits `date_start`/`date_end` or `[remote.pace]` omits `temporal_start`/`temporal_end`, those ranges are inferred from `utc_window_start`/`utc_window_end` (calendar days for GES DISC; full timestamps for Earthdata). Set dates explicitly under `[remote.*]` only when the fetch interval should differ from the coincidence filter window. For PACE, `use_bounding_box = false` gives a temporal-only `earthaccess` search; use `bounding_box = [...]` when `use_bounding_box = true`.

## Files

| File | Role |
|------|------|
| `run_match.jl` | CLI entrypoint |
| `MatchCoincidence.jl` | Parsing, UTC window filter, pairing, CSV/JSONL export |
| `RemoteFetch.jl` | GES DISC TROPOMI crawl/download; invokes `fetch_pace_earthaccess.py` for PACE |
| `fetch_pace_earthaccess.py` | `earthaccess` search + download (install with `pip install earthaccess`) |
