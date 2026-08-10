# Full-Spectrum Red SIF Retrieval for PACE Ocean Color Instrument

## Prerequisites

- Julia **1.9+** (see `Project.toml`)
- Optional: Python **3.11+** with `earthaccess` for downloading PACE products
- External data (not in git): L1B / L2 AOP (+ optional L2 BGC), transmittance PC NetCDFs, `SIF_singular_vector.jld2`, PACE SNR LUT

## Quickstart

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
cp svd_fit_lowres/global_fit_pipeline.example.toml svd_fit_lowres/global_fit_pipeline.toml
# edit global_fit_pipeline.toml: [general] dirs, [data] files, [svd_retrieval] mode/output
```


| Mode                           | Command                                                                                                                                                          |
| ------------------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **(a) Granule**                | Set `mode = "granule"` + `granule_id`, then `julia --project=. -t 8 svd_fit_lowres/svd_retrieval/run_svd_fit.jl svd_fit_lowres/global_fit_pipeline.toml`         |
| **(b) Day / dates**            | Set `mode = "date"` + `date`, or `mode = "dates"` + range/list, same `run_svd_fit.jl`                                                                            |
| **(c) Single pixel + spectra** | Set `granule_id`, `pixel`, `scan` (and paths), then `julia --project=. svd_fit_lowres/svd_retrieval/run_single_pixel.jl svd_fit_lowres/global_fit_pipeline.toml` |


Download (optional):

```bash
python svd_fit_lowres/download/download_pace_products.py svd_fit_lowres/global_fit_pipeline.toml
```

Rasterize (post-process):

```bash
julia --project=. svd_fit_lowres/rasterize/rasterize_svd_retrievals.jl \
  svd_fit_lowres/rasterize/rasterize.example.toml
```

Copy local TOMLs are gitignored (`svd_fit_lowres/global_fit_pipeline.toml`). Commit only `*.example.toml`.

See [svd_fit_lowres/README.md](svd_fit_lowres/README.md) for threading, filters, and config details.