#!/usr/bin/env python3
"""Stage 2: pixel match + lookup existing TROPOMI/PACE SIF → matchup_data.nc.

Examples
--------
  # Smoke (uses configs/sif_matchup.smoke.toml defaults)
  python TROPOMI_OCI_colocation/run_sif_matchup.py \\
      --config TROPOMI_OCI_colocation/configs/sif_matchup.smoke.toml

  # From an existing stage-1 JSON
  python TROPOMI_OCI_colocation/run_sif_matchup.py \\
      --config TROPOMI_OCI_colocation/configs/sif_matchup.smoke.toml \\
      --colocation-json TROPOMI_OCI_colocation/notebooks/co-location_results_202501.json \\
      --max-swaths 2
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore

from colocation.config import Config  # noqa: E402
from colocation.sif_matchup import run_sif_matchup  # noqa: E402

_PLACEHOLDER = Path("/dev/null")


def load_sif_toml(path: Path) -> dict:
    raw = tomllib.loads(Path(path).read_text())
    if not isinstance(raw, dict):
        raise ValueError(f"TOML root must be a table: {path}")
    return raw


def config_from_toml(raw: dict) -> tuple[Config, Path, Path, Path | None]:
    """Build Config + retrieval roots from sif_matchup TOML."""
    l1 = raw.get("l1") or {}
    ret = raw.get("retrieval") or {}
    match = raw.get("match") or {}
    window = raw.get("window") or {}
    out = raw.get("output") or {}

    pace_dir = Path(l1.get("pace_dir", "/kiwi-data/Data/satellite/PACE_OCI/L1B_V3"))
    tropomi_dir = Path(l1.get("tropomi_dir", "/net/squid/data1/projects/TROPOMI/ESA/L1"))
    tropomi_ret = Path(
        ret.get(
            "tropomi_ret_dir",
            "/home/zhe2/data/MyProjects/TROPOMI_SIF/results/b5",
        )
    )
    pace_ret = Path(
        ret.get(
            "pace_ret_dir",
            "/home/zhe2/data/PACE/new_svd_retrieval_output",
        )
    )
    output_dir = Path(
        out.get(
            "output_dir",
            "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/TROPOMI_OCI_colocation_data/"
            "runs/sif_matchup_smoke",
        )
    )
    coloc_json = out.get("colocation_json")
    colocation_json = Path(coloc_json) if coloc_json else None

    cfg = Config(
        pace_dir=pace_dir,
        tropomi_dir=tropomi_dir,
        rsr_path=_PLACEHOLDER,
        snr_path=_PLACEHOLDER,
        output_dir=output_dir,
        year=int(window.get("year", 2025)),
        month=int(window.get("month", 1)),
        day_start=int(window.get("day_start", 1)),
        day_end=int(window.get("day_end", 7)),
        tropomi_product=str(l1.get("tropomi_product", "OFFL")),
        pace_stride=int(l1.get("pace_stride", 1)),
        pace_duration_min=int(l1.get("pace_duration_min", 5)),
        dt_max_min_granule=float(l1.get("dt_max_min_granule", 60.0)),
        bbox_margin_deg=float(l1.get("bbox_margin_deg", 0.5)),
        tropo_bbox_scan_stride=int(l1.get("tropo_bbox_scan_stride", 40)),
        tropo_bbox_pix_stride=int(l1.get("tropo_bbox_pix_stride", 20)),
        tropo_chunk_scans=int(l1.get("tropo_chunk_scans", 400)),
        max_dist_km=float(match.get("max_dist_km", 2.5)),
        max_dt_min=float(match.get("max_dt_min", 30.0)),
        lt_max_oci=float(match.get("lt_max_oci", 30.0)),
        lt_mask_wl_min_oci=float(match.get("lt_mask_wl_min_oci", 600.0)),
        pace_stride_pix=int(match.get("pace_stride_pix", 150)),
        tropo_scan_stride=int(match.get("tropo_scan_stride", 2)),
        tropo_pix_stride=int(match.get("tropo_pix_stride", 2)),
    )
    return cfg, tropomi_ret, pace_ret, colocation_json


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Stage 2: pixel match + SIF lookup → matchup_data.nc"
    )
    p.add_argument(
        "--config",
        type=Path,
        default=_HERE / "configs" / "sif_matchup.smoke.toml",
        help="TOML config (default: configs/sif_matchup.smoke.toml)",
    )
    p.add_argument(
        "--colocation-json",
        type=Path,
        default=None,
        help="Stage-1 JSON (overrides config output.colocation_json)",
    )
    p.add_argument("--max-swaths", type=int, default=None)
    p.add_argument("--force-pairs", action="store_true")
    p.add_argument("--force-matches", action="store_true")
    args = p.parse_args(argv)

    raw = load_sif_toml(args.config)
    cfg, trop_ret, pace_ret, coloc_json = config_from_toml(raw)
    if args.colocation_json is not None:
        coloc_json = args.colocation_json
    max_swaths = args.max_swaths
    if max_swaths is None and (raw.get("output") or {}).get("max_swaths") is not None:
        max_swaths = int(raw["output"]["max_swaths"])

    out_nc, meta = run_sif_matchup(
        cfg,
        tropomi_ret_dir=trop_ret,
        pace_ret_dir=pace_ret,
        colocation_json=coloc_json,
        force_pairs=args.force_pairs,
        force_matches=args.force_matches,
        max_swaths=max_swaths,
    )
    print(json.dumps(meta, indent=2))
    print(f"Done: {out_nc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
