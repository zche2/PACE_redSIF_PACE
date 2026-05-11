#!/usr/bin/env python3
"""
Download PACE OCI L1B (V3), L2 BGC (V3.1), L2 AOP (V3.1) via NASA earthaccess.

Usage:
  python download_pace_products.py /path/to/global_fit_pipeline.toml

Requires Python 3.11+ (stdlib tomllib) and: pip install earthaccess

Config layout:
  • [general] — optional shared L1B_dir, L2AOP_dir, L2BGC_dir used when [download.paths]
    omits l1b_dir / l2_aop_dir / l2_bgc_dir.
  • [download.*] — temporal, spatial, collections, paths, options.

Config [download.options] login_strategy (default all): environment → netrc → interactive.
"""
from __future__ import annotations

import inspect
import os
import sys
import tomllib
from pathlib import Path


def _die(msg: str, code: int = 1) -> None:
    print(msg, file=sys.stderr)
    raise SystemExit(code)


def _get_str(d: dict, key: str, default: str = "") -> str:
    v = d.get(key, default)
    return default if v is None else str(v).strip()


def _optional_version_kw(version_raw: str) -> dict:
    v = str(version_raw).strip()
    return {} if v == "" else {"version": v}


def _build_search_kwargs(
    *,
    short_name: str,
    version_raw: str,
    temporal: tuple[str, str],
    use_bbox: bool,
    bbox: tuple[float, float, float, float] | None,
    cloud_cover: tuple[float, float] | None,
    apply_cloud: bool,
) -> dict:
    kw: dict = {
        "short_name": short_name,
        "temporal": temporal,
    }
    kw.update(_optional_version_kw(version_raw))
    if use_bbox and bbox is not None:
        kw["bounding_box"] = bbox
    if apply_cloud and cloud_cover is not None:
        kw["cloud_cover"] = cloud_cover
    return kw


def _nonempty_str_dir(paths: dict, path_key: str, general: dict, gen_key: str) -> str:
    raw = paths.get(path_key)
    if raw is not None and str(raw).strip():
        return str(raw).strip()
    return str(general.get(gen_key, "") or "").strip()


def _download_force_kw(download_fn, force: bool) -> dict:
    sig = inspect.signature(download_fn)
    if "force" not in sig.parameters:
        return {}
    return {"force": force}


def main() -> None:
    if len(sys.argv) != 2:
        _die("usage: download_pace_products.py <global_fit_pipeline.toml>")
    cfg_path = Path(sys.argv[1]).expanduser().resolve()
    if not cfg_path.is_file():
        _die(f"config not found: {cfg_path}")

    with open(cfg_path, "rb") as f:
        cfg = tomllib.load(f)

    general = cfg.get("general")
    if general is not None and not isinstance(general, dict):
        _die("[general] must be a table if present")

    dl = cfg.get("download")
    if not isinstance(dl, dict):
        _die(
            "Missing [download] table. Use nested keys: "
            "[download.temporal], [download.spatial], [download.collections], "
            "[download.paths], [download.options] — see global_fit_pipeline.example.toml"
        )

    temporal = dl.get("temporal") or {}
    t0 = _get_str(temporal, "start")
    t1 = _get_str(temporal, "end")
    if not t0 or not t1:
        _die("[download.temporal] start and end are required")

    spatial = dl.get("spatial") or {}
    use_bbox = bool(spatial.get("use_bounding_box", False))
    bb_raw = spatial.get("bounding_box")
    bbox: tuple[float, float, float, float] | None = None
    if use_bbox:
        if not isinstance(bb_raw, list) or len(bb_raw) != 4:
            _die(
                "[download.spatial] bounding_box must be [west, south, east, north] "
                "when use_bounding_box=true"
            )
        bbox = tuple(float(x) for x in bb_raw)

    col = dl.get("collections") or {}
    l1b_sn = _get_str(col, "l1b_short_name", "PACE_OCI_L1B_SCI")
    l1b_ver = _get_str(col, "l1b_version", "3")
    bgc_sn = _get_str(col, "l2_bgc_short_name", "PACE_OCI_L2_BGC")
    bgc_ver = _get_str(col, "l2_bgc_version", "3.1")
    aop_sn = _get_str(col, "l2_aop_short_name", "PACE_OCI_L2_AOP")
    aop_ver = _get_str(col, "l2_aop_version", "3.1")

    paths = dl.get("paths") or {}
    gen = general if isinstance(general, dict) else {}
    l1b_s = _nonempty_str_dir(paths, "l1b_dir", gen, "L1B_dir")
    l2_bgc_s = _nonempty_str_dir(paths, "l2_bgc_dir", gen, "L2BGC_dir")
    l2_aop_s = _nonempty_str_dir(paths, "l2_aop_dir", gen, "L2AOP_dir")
    if not l1b_s or not l2_bgc_s or not l2_aop_s:
        _die(
            "Download destinations required: set [general] L1B_dir, L2AOP_dir, L2BGC_dir "
            "and/or [download.paths] l1b_dir, l2_aop_dir, l2_bgc_dir"
        )
    l1b_dir = Path(os.path.expanduser(l1b_s))
    l2_bgc_dir = Path(os.path.expanduser(l2_bgc_s))
    l2_aop_dir = Path(os.path.expanduser(l2_aop_s))

    opts = dl.get("options") or {}
    cc_raw = opts.get("cloud_cover")
    cloud_cover: tuple[float, float] | None = None
    if cc_raw is not None:
        if not isinstance(cc_raw, list) or len(cc_raw) != 2:
            _die("[download.options] cloud_cover must be [min_pct, max_pct]")
        cloud_cover = (float(cc_raw[0]), float(cc_raw[1]))
    skip_existing = opts.get("skip_existing", True)
    if not isinstance(skip_existing, bool):
        _die("[download.options] skip_existing must be a boolean")
    force_redownload = not skip_existing

    login_strategy = _get_str(opts, "login_strategy", "all").lower()
    allowed_login = {"all", "environment", "netrc", "interactive"}
    if login_strategy not in allowed_login:
        _die(f"[download.options] login_strategy must be one of: {', '.join(sorted(allowed_login))}")

    try:
        import earthaccess
    except ImportError:
        _die("earthaccess not installed (pip install earthaccess)")

    try:
        auth = earthaccess.login(strategy=login_strategy)
    except Exception as e:
        _die(
            "earthaccess.login failed.\n"
            f"  login_strategy={login_strategy!r}\n"
            "Credentials:\n"
            "  • environment: EARTHDATA_USERNAME + EARTHDATA_PASSWORD, or EARTHDATA_TOKEN\n"
            "  • netrc: ~/.netrc (or NETRC) with machine urs.earthdata.nasa.gov\n"
            "  • all: environment → netrc → interactive\n"
            f"Details: {type(e).__name__}: {e}"
        )
    if not getattr(auth, "authenticated", False):
        _die(
            "earthaccess login did not authenticate.\n"
            "Check credentials or set [download.options] login_strategy."
        )
    temporal_tuple = (t0, t1)

    jobs = [
        ("L1B", l1b_sn, l1b_ver, l1b_dir, False),
        ("L2_BGC", bgc_sn, bgc_ver, l2_bgc_dir, True),
        ("L2_AOP", aop_sn, aop_ver, l2_aop_dir, True),
    ]

    for label, sn, ver, out_dir, apply_cloud in jobs:
        out_dir.mkdir(parents=True, exist_ok=True)
        kw = _build_search_kwargs(
            short_name=sn,
            version_raw=ver,
            temporal=temporal_tuple,
            use_bbox=use_bbox,
            bbox=bbox,
            cloud_cover=cloud_cover,
            apply_cloud=apply_cloud,
        )
        try:
            granules = earthaccess.search_data(**kw)
        except Exception as e:
            _die(f"[{label}] search_data failed: {e}")
        n = len(granules) if granules else 0
        print(f"[{label}] {sn} — found {n} granule(s) → {out_dir}")
        if not granules:
            continue
        dl_kw = _download_force_kw(earthaccess.download, force_redownload)
        if dl_kw:
            print(f"[{label}] skip_existing={skip_existing} → download(..., force={dl_kw['force']})")
        else:
            print(
                f"[{label}] skip_existing={skip_existing} "
                "(earthaccess.download has no force=; relying on library defaults)"
            )
        try:
            files = earthaccess.download(granules, local_path=str(out_dir), **dl_kw)
        except Exception as e:
            _die(f"[{label}] download failed: {e}")
        if files is None:
            files = []
        print(
            f"[{label}] {len(files)} local file path(s) "
            "(present files are reused when skip_existing=true)"
        )

    print("Done.")


if __name__ == "__main__":
    main()
