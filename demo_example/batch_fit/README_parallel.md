# Parallel Full-Swath Batch Retrieval

`Run_batch_full_nc_parallel.jl` is a threaded version of `Run_batch_full_nc.jl` that processes pixels within each scan in parallel. It uses the same config and produces the same output format, but writes to a different file by default to avoid overwriting sequential results.

## Usage

Run with multiple threads (recommended):

```bash
julia -t 8 demo_example/batch_fit/Run_batch_full_nc_parallel.jl
```

Or let Julia use all available cores:

```bash
julia -t auto demo_example/batch_fit/Run_batch_full_nc_parallel.jl
```

Config path can be set via environment variable:

```bash
PACE_MWE_CONFIG=/path/to/config.toml julia -t 8 demo_example/batch_fit/Run_batch_full_nc_parallel.jl
```

## Config: `[batch_fit]`

Same keys as `Run_batch_full_nc.jl`, plus:

- **`use_threads`** (default `true`) — Enable thread-based pixel parallelization. Set to `false` to run sequentially (useful for debugging).
- **`output_suffix_parallel`** (default `"_retrieval_full_parallel.nc"`) — Output file suffix for parallel runs. Ensures parallel and sequential outputs do not overwrite each other.

## How It Works

- Each thread gets its own retrieval core (forward model + Jacobian evaluator) with independent scratch buffers, satisfying the thread-safety requirements of the preallocated forward model.
- Pixels within each scan are processed in parallel via `Threads.@threads`.
- Scans are processed sequentially (read slab → parallel pixel retrieval → write slab).
- Output format and status codes are identical to `Run_batch_full_nc.jl`.

## Dependencies

Requires `Run_batch_full_nc.jl` and `Fit_toy_forward_model.jl` (included automatically). No additional packages beyond those used by the sequential script.
