#!/usr/bin/env python3
"""Download PACE OCI L1B (or other collections) via earthaccess; invoked by RemoteFetch.jl."""
from __future__ import annotations

import json
import os
import sys


def main() -> None:
    if len(sys.argv) < 2:
        print(json.dumps({"ok": False, "error": "usage: fetch_pace_earthaccess.py <config.json>"}))
        sys.exit(2)
    path = sys.argv[1]
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    try:
        import earthaccess
    except ImportError:
        print(
            json.dumps(
                {
                    "ok": False,
                    "error": "earthaccess not installed for this Python (pip install earthaccess)",
                }
            )
        )
        sys.exit(1)
    earthaccess.login(strategy="environment")
    temporal = (cfg["temporal_start"], cfg["temporal_end"])
    kw: dict = {
        "short_name": cfg["short_name"],
        "version": str(cfg["version"]),
        "temporal": temporal,
    }
    if cfg.get("use_bounding_box", True):
        kw["bounding_box"] = tuple(cfg["bounding_box"])
    granules = earthaccess.search_data(**kw)
    out_dir = cfg["cache_dir"]
    os.makedirs(out_dir, exist_ok=True)
    if not granules:
        print(json.dumps({"ok": True, "paths": []}))
        return
    files = earthaccess.download(granules, local_path=out_dir)
    paths: list[str] = []
    if files is None:
        files = []
    if files:
        for p in files:
            paths.append(os.path.abspath(str(p)))
    print(json.dumps({"ok": True, "paths": paths}))


if __name__ == "__main__":
    main()
