"""Stage 2: pixel match + lookup existing TROPOMI/PACE SIF → matchup_data.nc."""

from __future__ import annotations

import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import netCDF4 as nc
import numpy as np

from . import cache
from .config import Config
from .match import run_matches_stage
from .pair import run_pairs_stage
from .sif_join import (
    fill_pace_svd,
    fill_tropomi_fs_ret,
    index_fs_ret_dir,
    orbit_from_tropomi_path,
    path_str,
    stamp_from_pace_path,
)


_PACE_STAMP_RE = re.compile(r"PACE_OCI\.(\d{8}T\d{6})")


def resolve_pace_l1b(pace_dir: Path, name: str) -> Path:
    """Resolve L1B granule name under dated YYYY/MM/DD or flat root."""
    pace_dir = Path(pace_dir)
    flat = pace_dir / name
    if flat.is_file():
        return flat
    m = _PACE_STAMP_RE.search(name)
    if m:
        stamp = m.group(1)
        dated = pace_dir / stamp[:4] / stamp[4:6] / stamp[6:8] / name
        if dated.is_file():
            return dated
    return flat


def load_swath_list_from_json(path: Path) -> tuple[list[dict], dict[str, Any]]:
    """Load notebook-schema or raw swath_list.json; ensure ``pace_paths`` present."""
    path = Path(path)
    raw = cache.read_json(path)
    meta: dict[str, Any] = {}
    if isinstance(raw, dict) and "swaths" in raw:
        meta = dict(raw.get("meta") or {})
        swaths = list(raw["swaths"])
        pace_dir = Path(meta["pace_dir"]) if meta.get("pace_dir") else None
    elif isinstance(raw, list):
        swaths = raw
        pace_dir = None
    else:
        raise ValueError(f"Unrecognized colocation JSON schema: {path}")

    out: list[dict] = []
    for s in swaths:
        s = dict(s)
        if "pace_paths" not in s or not s["pace_paths"]:
            if pace_dir is None:
                raise ValueError(
                    f"Swath orbit={s.get('orbit')} missing pace_paths and meta.pace_dir"
                )
            names = s.get("pace_granules") or []
            s["pace_paths"] = [str(resolve_pace_l1b(pace_dir, n)) for n in names]
            s["pace_granules"] = list(names)
            s["n_pace_matches"] = len(s["pace_paths"])
        else:
            # Re-resolve if stored as flat names under dated root
            if pace_dir is not None:
                fixed = []
                for p in s["pace_paths"]:
                    pp = Path(p)
                    if pp.is_file():
                        fixed.append(str(pp))
                    else:
                        fixed.append(str(resolve_pace_l1b(pace_dir, pp.name)))
                s["pace_paths"] = fixed
            if "pace_granules" not in s:
                s["pace_granules"] = [Path(p).name for p in s["pace_paths"]]
                s["n_pace_matches"] = len(s["pace_paths"])
        out.append(s)
    return out, meta


def write_matchup_data_nc(
    path: Path,
    matches: dict[str, np.ndarray],
    *,
    trop_dir: Path,
    pace_ret_dir: Path,
    sif_trop: np.ndarray,
    chi2_trop: np.ndarray,
    found_trop: np.ndarray,
    sif_pace: np.ndarray,
    chi2_pace: np.ndarray,
    converged_pace: np.ndarray,
    found_pace: np.ndarray,
    fingerprints: dict[str, str] | None = None,
    extra_attrs: dict[str, Any] | None = None,
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n = int(matches["dist_km"].size)
    tmp = path.with_suffix(path.suffix + ".tmp")
    if tmp.exists():
        tmp.unlink()

    with nc.Dataset(tmp, "w") as ds:
        ds.createDimension("match", n)

        def _v(name, data, dtype, **attrs):
            var = ds.createVariable(name, dtype, ("match",), zlib=True, complevel=4)
            var[:] = data
            for k, v in attrs.items():
                setattr(var, k, v)
            return var

        for name, dtype in (
            ("pace_lat", "f4"),
            ("pace_lon", "f4"),
            ("trop_lat", "f4"),
            ("trop_lon", "f4"),
            ("dist_km", "f4"),
            ("dt_min", "f4"),
            ("pace_sza", "f4"),
            ("pace_vza", "f4"),
            ("trop_sza", "f4"),
            ("trop_vza", "f4"),
        ):
            if name in matches:
                _v(name, np.asarray(matches[name][:n], dtype=np.float32), dtype)

        for name in ("trop_pix", "trop_scan", "pace_pix", "pace_scan", "swath_id"):
            if name in matches:
                _v(name, np.asarray(matches[name][:n], dtype=np.int32), "i4")

        _v(
            "sif_tropomi",
            sif_trop,
            "f4",
            long_name="TROPOMI BD5 fs_ret RETRIEVAL_RESULT/sif",
            source_dir=str(trop_dir),
            index_note="detector_pixel=trop_pix+1, scanline=trop_scan+1",
        )
        _v("chi2_tropomi", chi2_trop, "f4", long_name="TROPOMI fs_ret chi2")
        _v("trop_found", found_trop, "u1", long_name="1 if TROPOMI SIF found")
        _v(
            "sif_pace_678nm",
            sif_pace,
            "f4",
            long_name="PACE SVD sif_radiance_678nm",
            units="W m-2 sr-1 um-1",
            source_dir=str(pace_ret_dir),
            index_note="source_pixel_index=pace_pix+1, source_scan_index=pace_scan+1",
        )
        _v("chi2_pace", chi2_pace, "f4", long_name="PACE SVD reduced_chi2")
        _v("converged_pace", converged_pace, "u1", long_name="PACE SVD converged")
        _v("pace_found", found_pace, "u1", long_name="1 if PACE SIF found")

        # paths as fixed-length strings
        trop_paths = [path_str(v) for v in matches["tropomi_path"][:n]]
        pace_paths = [path_str(v) for v in matches["pace_path"][:n]]
        maxlen_t = max((len(p) for p in trop_paths), default=1)
        maxlen_p = max((len(p) for p in pace_paths), default=1)
        vt = ds.createVariable("tropomi_path", f"S{maxlen_t}", ("match",))
        vp = ds.createVariable("pace_path", f"S{maxlen_p}", ("match",))
        vt[:] = np.asarray(trop_paths, dtype=f"S{maxlen_t}")
        vp[:] = np.asarray(pace_paths, dtype=f"S{maxlen_p}")

        ds.title = "TROPOMI–PACE co-located SIF matchup (lookup from existing retrievals)"
        ds.history = datetime.now(timezone.utc).isoformat() + " run_sif_matchup"
        ds.tropomi_retrieval_dir = str(trop_dir)
        ds.pace_retrieval_dir = str(pace_ret_dir)
        ds.n_match = np.int32(n)
        if fingerprints:
            ds.fingerprint_pairs = fingerprints.get("pairs", "")
            ds.fingerprint_matches = fingerprints.get("matches", "")
            ds.fingerprint_sif = fingerprints.get("sif", "")
        if extra_attrs:
            for k, v in extra_attrs.items():
                setattr(ds, k, v)

    tmp.replace(path)
    return path


def run_sif_matchup(
    cfg: Config,
    *,
    tropomi_ret_dir: Path,
    pace_ret_dir: Path,
    colocation_json: Path | None = None,
    force_pairs: bool = False,
    force_matches: bool = False,
    max_swaths: int | None = None,
    log=print,
) -> tuple[Path, dict]:
    """Run pairs (optional) → matches → SIF lookup → products/matchup_data.nc."""
    dirs = cache.ensure_output_dirs(cfg.output_dir)
    tropomi_ret_dir = Path(tropomi_ret_dir)
    pace_ret_dir = Path(pace_ret_dir)

    if colocation_json is not None:
        swath_list, coloc_meta = load_swath_list_from_json(colocation_json)
        pairs_fp = cache.fingerprint(
            {
                "colocation_json": str(Path(colocation_json).resolve()),
                "n_swaths": len(swath_list),
                "meta": {k: coloc_meta.get(k) for k in ("year", "month", "day_start", "day_end", "pace_dir")},
            }
        )
        cache.write_json(dirs["pairs"] / "swath_list.json", swath_list)
        cache.write_json(
            dirs["pairs"] / "pairs_meta.json",
            {"fingerprint": pairs_fp, "n_swaths": len(swath_list), "source": str(colocation_json)},
        )
        log(f"[pairs] loaded {len(swath_list)} swaths from {colocation_json}")
    else:
        swath_list, pairs_fp = run_pairs_stage(
            cfg, dirs, force=force_pairs, log=log
        )

    if max_swaths is not None:
        swath_list = swath_list[: int(max_swaths)]
        log(f"[pairs] truncated to max_swaths={max_swaths}")

    matches, matches_fp = run_matches_stage(
        cfg,
        dirs,
        swath_list,
        pairs_fp,
        force=force_matches or force_pairs,
        log=log,
    )

    n = int(matches["dist_km"].size)
    trop_paths = [path_str(v) for v in matches["tropomi_path"][:n]]
    pace_paths = [path_str(v) for v in matches["pace_path"][:n]]
    orbits = np.asarray([orbit_from_tropomi_path(p) or "" for p in trop_paths])
    stamps = np.asarray([stamp_from_pace_path(p) or "" for p in pace_paths])

    sif_trop = np.full(n, np.nan, dtype=np.float32)
    chi2_trop = np.full(n, np.nan, dtype=np.float32)
    found_trop = np.zeros(n, dtype=np.uint8)
    sif_pace = np.full(n, np.nan, dtype=np.float32)
    chi2_pace = np.full(n, np.nan, dtype=np.float32)
    converged_pace = np.zeros(n, dtype=np.uint8)
    found_pace = np.zeros(n, dtype=np.uint8)

    orbit_to_file = index_fs_ret_dir(tropomi_ret_dir)
    t0 = time.time()
    stats_t = fill_tropomi_fs_ret(
        orbit_to_file,
        orbits,
        np.asarray(matches["trop_pix"][:n], dtype=np.int32),
        np.asarray(matches["trop_scan"][:n], dtype=np.int32),
        sif_trop,
        chi2_trop,
        found_trop,
    )
    t1 = time.time()
    stats_p = fill_pace_svd(
        pace_ret_dir,
        stamps,
        np.asarray(matches["pace_pix"][:n], dtype=np.int32),
        np.asarray(matches["pace_scan"][:n], dtype=np.int32),
        sif_pace,
        chi2_pace,
        converged_pace,
        found_pace,
    )
    t2 = time.time()
    log(
        f"[sif] TROPOMI hit={stats_t['n_hit']} file_miss={stats_t['n_file_miss']} "
        f"pix_miss={stats_t['n_pix_miss']} ({t1 - t0:.1f}s)"
    )
    log(
        f"[sif] PACE hit={stats_p['n_hit']} file_miss={stats_p['n_file_miss']} "
        f"pix_miss={stats_p['n_pix_miss']} ({t2 - t1:.1f}s)"
    )

    sif_fp = cache.fingerprint(
        {
            "matches_fp": matches_fp,
            "tropomi_ret_dir": str(tropomi_ret_dir.resolve()),
            "pace_ret_dir": str(pace_ret_dir.resolve()),
            "n_fs_ret": len(orbit_to_file),
        }
    )
    fingerprints = {"pairs": pairs_fp, "matches": matches_fp, "sif": sif_fp}

    out_nc = dirs["products"] / "matchup_data.nc"
    write_matchup_data_nc(
        out_nc,
        matches,
        trop_dir=tropomi_ret_dir,
        pace_ret_dir=pace_ret_dir,
        sif_trop=sif_trop,
        chi2_trop=chi2_trop,
        found_trop=found_trop,
        sif_pace=sif_pace,
        chi2_pace=chi2_pace,
        converged_pace=converged_pace,
        found_pace=found_pace,
        fingerprints=fingerprints,
        extra_attrs={
            "colocation_json": str(colocation_json) if colocation_json else "",
        },
    )
    log(f"[products] wrote {out_nc}")

    meta = {
        "output_nc": str(out_nc),
        "n_match": n,
        "tropomi": {**stats_t, "dir": str(tropomi_ret_dir), "seconds": round(t1 - t0, 2)},
        "pace": {**stats_p, "dir": str(pace_ret_dir), "seconds": round(t2 - t1, 2)},
        "fingerprints": fingerprints,
    }
    cache.write_json(dirs["products"] / "matchup_data_meta.json", meta)
    cache.write_run_meta(
        cfg.output_dir,
        {
            **cfg.to_plain_dict(),
            "tropomi_ret_dir": str(tropomi_ret_dir),
            "pace_ret_dir": str(pace_ret_dir),
        },
        extra={"stages": ["pairs", "matches", "sif_join"], "matchup_data": meta},
    )
    return out_nc, meta
