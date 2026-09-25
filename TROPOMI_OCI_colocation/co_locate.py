#!/usr/bin/env python3
"""Stage 1: TROPOMI BD5 ↔ PACE OCI L1B time+bbox co-location → JSON.

Uses ``colocation.pair`` / ``colocation.discover``. Writes a notebook-schema
JSON (``meta`` + ``swaths``) and optionally a pipeline ``swath_list.json``.

Examples
--------
  python TROPOMI_OCI_colocation/co_locate.py --year 2025 --month 1 --day-end 31 \\
      --output TROPOMI_OCI_colocation/notebooks/co-location_results_202501.json

  # Quick smoke (first 7 days)
  python TROPOMI_OCI_colocation/co_locate.py --month 1 --day-end 7 --output /tmp/coloc_smoke.json
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from colocation.config import Config  # noqa: E402
from colocation.pair import build_swath_list  # noqa: E402

_DEFAULT_PACE_DIR = Path("/kiwi-data/Data/satellite/PACE_OCI/L1B_V3")
_DEFAULT_TROPOMI_DIR = Path("/net/squid/data1/projects/TROPOMI/ESA/L1")
# pair stage does not use RSR/SNR; placeholders satisfy Config
_PLACEHOLDER = Path("/dev/null")


def build_payload(
    swath_list: list[dict],
    *,
    year: int,
    month: int,
    day_start: int,
    day_end: int,
    dt_max_min: float,
    bbox_margin_deg: float,
    pace_dir: Path,
    tropomi_dir: Path,
    tropomi_product: str,
) -> dict:
    # Notebook schema: pace_granules names only (drop pace_paths for compact JSON)
    swaths_out = []
    for s in swath_list:
        swaths_out.append(
            {
                "orbit": s["orbit"],
                "tropomi": s["tropomi"],
                "tropomi_start": s["tropomi_start"],
                "tropomi_end": s["tropomi_end"],
                "tropomi_path": s["tropomi_path"],
                "n_pace_matches": s["n_pace_matches"],
                "pace_granules": list(s["pace_granules"]),
            }
        )
    return {
        "meta": {
            "created": datetime.now(timezone.utc).isoformat(),
            "description": "TROPOMI BD5 ↔ PACE OCI L1B time+bbox co-location (stage 1)",
            "year": year,
            "month": month,
            "day_start": day_start,
            "day_end": day_end,
            "dt_max_min": dt_max_min,
            "bbox_margin_deg": bbox_margin_deg,
            "pace_dir": str(pace_dir),
            "tropomi_dir": str(tropomi_dir),
            "tropomi_product": tropomi_product,
            "n_swaths": len(swaths_out),
        },
        "swaths": swaths_out,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Stage 1: TROPOMI BD5 ↔ PACE OCI L1B granule co-location → JSON"
    )
    p.add_argument("--year", type=int, default=2025)
    p.add_argument("--month", type=int, default=7)
    p.add_argument("--day-start", type=int, default=1)
    p.add_argument("--day-end", type=int, default=31)
    p.add_argument("--pace-dir", type=Path, default=_DEFAULT_PACE_DIR)
    p.add_argument("--tropomi-dir", type=Path, default=_DEFAULT_TROPOMI_DIR)
    p.add_argument(
        "--tropomi-product",
        choices=("OFFL", "RPRO", "BOTH"),
        default="OFFL",
    )
    p.add_argument("--pace-stride", type=int, default=1)
    p.add_argument("--pace-duration-min", type=float, default=5.0)
    p.add_argument("--dt-max-min", type=float, default=60.0)
    p.add_argument("--bbox-margin-deg", type=float, default=0.5)
    p.add_argument("--tropo-scan-stride", type=int, default=40)
    p.add_argument("--tropo-pix-stride", type=int, default=20)
    p.add_argument("--tropo-chunk-scans", type=int, default=400)
    p.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help="Output JSON path (default: notebooks/co-location_results_YYYYMM.json)",
    )
    p.add_argument(
        "--pairs-dir",
        type=Path,
        default=None,
        help="Also write pipeline-style swath_list.json under this directory",
    )
    p.add_argument(
        "--print-swaths",
        action="store_true",
        help="Print each co-located TROPOMI swath to stdout",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    out = args.output
    if out is None:
        out = _HERE / "notebooks" / f"co-location_results_{args.year}{args.month:02d}.json"

    cfg = Config(
        pace_dir=args.pace_dir,
        tropomi_dir=args.tropomi_dir,
        rsr_path=_PLACEHOLDER,
        snr_path=_PLACEHOLDER,
        output_dir=out.parent,
        year=args.year,
        month=args.month,
        day_start=args.day_start,
        day_end=args.day_end,
        tropomi_product=args.tropomi_product,
        pace_stride=args.pace_stride,
        pace_duration_min=int(args.pace_duration_min),
        dt_max_min_granule=args.dt_max_min,
        bbox_margin_deg=args.bbox_margin_deg,
        tropo_bbox_scan_stride=args.tropo_scan_stride,
        tropo_bbox_pix_stride=args.tropo_pix_stride,
        tropo_chunk_scans=args.tropo_chunk_scans,
    )

    swath_list = build_swath_list(cfg)

    if args.print_swaths:
        print("=" * 88)
        for i, s in enumerate(swath_list, 1):
            print(
                f"{i:3d}. orbit={s['orbit']}  "
                f"{s['tropomi_start'][:19]}Z → {s['tropomi_end'][11:19]}Z  "
                f"PACE matches={s['n_pace_matches']}"
            )
            print(f"     {s['tropomi']}")

    payload = build_payload(
        swath_list,
        year=args.year,
        month=args.month,
        day_start=args.day_start,
        day_end=args.day_end,
        dt_max_min=args.dt_max_min,
        bbox_margin_deg=args.bbox_margin_deg,
        pace_dir=args.pace_dir,
        tropomi_dir=args.tropomi_dir,
        tropomi_product=args.tropomi_product,
    )
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")
    print(f"Wrote {out}  (n_swaths={len(swath_list)})")

    if args.pairs_dir is not None:
        pairs_dir = Path(args.pairs_dir)
        pairs_dir.mkdir(parents=True, exist_ok=True)
        list_path = pairs_dir / "swath_list.json"
        with open(list_path, "w") as f:
            json.dump(swath_list, f, indent=2)
            f.write("\n")
        print(f"Wrote {list_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
