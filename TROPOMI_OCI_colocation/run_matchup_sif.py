#!/usr/bin/env python3
"""Thin CLI for match-up SIF retrieval (calls Julia).

Examples
--------
  python run_matchup_sif.py --config config.sif_retrieval.toml
  python run_matchup_sif.py --config config.sif_retrieval.toml -t 16
  python run_matchup_sif.py --config config.sif_retrieval.toml --max-matches 50
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parent  # PACE_redSIF_PACE (Julia Project.toml)
_JL = _HERE / "sif_retrieval" / "run_matchup_sif.jl"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Match-up SVD SIF retrieval (Julia)")
    p.add_argument("--config", required=True, type=Path, help="config.sif_retrieval.toml")
    p.add_argument("-t", "--threads", type=int, default=max(1, (os.cpu_count() or 8) // 2))
    p.add_argument(
        "--max-matches",
        type=int,
        default=None,
        help="Override [matchup].max_matches for a quick test (rewrites a temp TOML)",
    )
    p.add_argument("--julia", default=shutil.which("julia") or "julia")
    args = p.parse_args(argv)

    cfg = args.config.resolve()
    if not cfg.is_file():
        print(f"Config not found: {cfg}", file=sys.stderr)
        return 1
    if not _JL.is_file():
        print(f"Missing Julia entry: {_JL}", file=sys.stderr)
        return 1

    cfg_run = cfg
    tmp: Path | None = None
    if args.max_matches is not None:
        text = cfg.read_text()
        # crude override: append after [matchup] block key
        if "max_matches" in text:
            import re

            text = re.sub(
                r"(?m)^max_matches\s*=\s*.*$",
                f"max_matches = {int(args.max_matches)}",
                text,
                count=1,
            )
        else:
            text = text.replace("[matchup]", f"[matchup]\nmax_matches = {int(args.max_matches)}\n", 1)
        fd, name = tempfile.mkstemp(prefix="sif_retrieval_", suffix=".toml")
        os.close(fd)
        tmp = Path(name)
        tmp.write_text(text)
        cfg_run = tmp

    cmd = [
        args.julia,
        f"--project={_REPO}",
        f"-t{args.threads}",
        str(_JL),
        str(cfg_run),
    ]
    print(" ".join(cmd))
    try:
        return subprocess.call(cmd)
    finally:
        if tmp is not None:
            tmp.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
