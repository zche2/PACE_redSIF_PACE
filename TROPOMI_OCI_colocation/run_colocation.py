#!/usr/bin/env python3
"""CLI for TROPOMI–OCI co-location pipeline.

Examples
--------
  python run_colocation.py --config config.yaml
  python run_colocation.py --config config.yaml --force svd simulate
  python run_colocation.py --config config.yaml --from simulate
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running from this folder without installing the package
_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from colocation.cache import STAGES, force_set  # noqa: E402
from colocation.config import load_config  # noqa: E402
from colocation.pipeline import run_pipeline  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="TROPOMI BD5 ↔ PACE OCI co-location pipeline")
    p.add_argument("--config", required=True, type=Path, help="YAML config path")
    p.add_argument(
        "--force",
        nargs="+",
        default=[],
        metavar="STAGE",
        help=f"Rebuild stage(s): {'|'.join(STAGES)}|all",
    )
    p.add_argument(
        "--from",
        dest="from_stage",
        choices=list(STAGES),
        default=None,
        help="Start from this stage (earlier stages must be cached)",
    )
    args = p.parse_args(argv)

    cfg = load_config(args.config)
    forced = force_set(args.force)
    out = run_pipeline(cfg, force=forced, from_stage=args.from_stage)
    print(f"Done: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
