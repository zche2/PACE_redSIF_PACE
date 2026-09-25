#!/usr/bin/env python3
"""Thin plots driver over ``matchup_data.nc`` from the two-stage SIF pipeline.

Build the NetCDF with::

  python TROPOMI_OCI_colocation/run_sif_matchup.py \\
      --config TROPOMI_OCI_colocation/configs/sif_matchup.smoke.toml

Then plot::

  python TROPOMI_OCI_colocation/notebooks/compare_algorithm_on_oci_tropomi.py \\
      --matchup-nc .../products/matchup_data.nc

The full row-aware PC→678 compare stack is archived at
``archive/compare_algorithm_row_aware.py``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4 as nc  # noqa: E402
import numpy as np  # noqa: E402

COLOC_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(COLOC_ROOT))

DEFAULT_MATCHUP = Path(
    "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/TROPOMI_OCI_colocation_data/"
    "runs/sif_matchup_smoke/products/matchup_data.nc"
)
DEFAULT_FIG_DIR = Path(__file__).resolve().parent / "outputs" / "sif_matchup_plots"


def load_matchup(path: Path) -> dict:
    with nc.Dataset(path) as ds:
        data = {name: np.asarray(ds[name][:]) for name in ds.variables}
        data["_attrs"] = {
            "title": getattr(ds, "title", ""),
            "tropomi_retrieval_dir": getattr(ds, "tropomi_retrieval_dir", ""),
            "pace_retrieval_dir": getattr(ds, "pace_retrieval_dir", ""),
            "n_match": int(getattr(ds, "n_match", data.get("dist_km", np.array([])).size)),
        }
    return data


def scatter_stats(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    finite = np.isfinite(x) & np.isfinite(y)
    x, y = x[finite], y[finite]
    if x.size == 0:
        return {"n": 0, "bias": np.nan, "rmse": np.nan, "r": np.nan}
    return {
        "n": int(x.size),
        "bias": float(np.mean(y - x)),
        "rmse": float(np.sqrt(np.mean((y - x) ** 2))),
        "r": float(np.corrcoef(x, y)[0, 1]) if x.size > 1 else np.nan,
    }


def plot_all(data: dict, fig_dir: Path) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)

    trop = np.asarray(data["sif_tropomi"], dtype=float)
    pace = np.asarray(data["sif_pace_678nm"], dtype=float)
    both = (
        (np.asarray(data.get("trop_found", np.ones_like(trop))) > 0)
        & (np.asarray(data.get("pace_found", np.ones_like(pace))) > 0)
        & np.isfinite(trop)
        & np.isfinite(pace)
    )
    x, y = trop[both], pace[both]
    stats = scatter_stats(x, y)

    fig, ax = plt.subplots(figsize=(5.5, 5))
    if x.size:
        ax.scatter(x, y, s=8, alpha=0.5, c="#4C72B0", edgecolors="none")
        lim_lo = float(np.nanmin([x.min(), y.min()]))
        lim_hi = float(np.nanmax([x.max(), y.max()]))
        pad = 0.05 * (lim_hi - lim_lo + 1e-12)
        lims = (lim_lo - pad, lim_hi + pad)
        ax.plot(lims, lims, "k--", lw=1, alpha=0.6)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
    ax.set_xlabel("TROPOMI fs_ret sif")
    ax.set_ylabel("PACE sif_radiance_678nm")
    ax.set_title(
        f"SIF matchup  n={stats['n']:,}  "
        f"r={stats['r']:.3f}  bias={stats['bias']:.4g}  rmse={stats['rmse']:.4g}"
    )
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(fig_dir / "sif_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    if "dist_km" in data and "dt_min" in data:
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
        axes[0].hist(np.asarray(data["dist_km"], dtype=float), bins=40, color="#4C72B0", edgecolor="white")
        axes[0].set_xlabel("dist_km")
        axes[0].set_title("Distance")
        axes[1].hist(np.asarray(data["dt_min"], dtype=float), bins=40, color="#55A868", edgecolor="white")
        axes[1].set_xlabel("dt_min (TROPOMI − PACE)")
        axes[1].set_title("Δt")
        fig.suptitle(f"match QA (n={data['_attrs']['n_match']:,})", y=1.02)
        fig.tight_layout()
        fig.savefig(fig_dir / "match_qa_hist.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    if "pace_lon" in data and "pace_lat" in data:
        fig, ax = plt.subplots(figsize=(6, 5))
        sc = ax.scatter(
            np.asarray(data["pace_lon"], dtype=float),
            np.asarray(data["pace_lat"], dtype=float),
            c=np.asarray(data.get("dist_km", np.zeros(pace.shape)), dtype=float),
            s=8,
            cmap="viridis",
            alpha=0.8,
        )
        fig.colorbar(sc, ax=ax, label="dist_km")
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
        ax.set_title("co-located pixels")
        ax.set_aspect("equal", adjustable="datalim")
        fig.tight_layout()
        fig.savefig(fig_dir / "match_map.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    print(f"Wrote figures under {fig_dir}")
    print(f"stats: {stats}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Plot SIF matchup_data.nc")
    p.add_argument("--matchup-nc", type=Path, default=DEFAULT_MATCHUP)
    p.add_argument("--fig-dir", type=Path, default=DEFAULT_FIG_DIR)
    args = p.parse_args(argv)

    if not args.matchup_nc.is_file():
        raise FileNotFoundError(
            f"No matchup at {args.matchup_nc}; run run_sif_matchup.py first."
        )
    data = load_matchup(args.matchup_nc)
    plot_all(data, args.fig_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
