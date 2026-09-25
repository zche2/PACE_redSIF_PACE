"""Lookup TROPOMI fs_ret and PACE SVD SIF onto pixel match indices."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import h5py
import netCDF4 as nc
import numpy as np

_ORBIT_RE = re.compile(r"_(\d{5})_\d{2}_")
_PACE_STAMP_RE = re.compile(r"PACE_OCI\.(\d{8}T\d{6})")
_FS_RET_RE = re.compile(r"fs_ret_(\d{5})_b5_v10\.h5$")


def path_str(v) -> str:
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode("utf-8", errors="replace").strip("\x00")
    s = str(v)
    if len(s) >= 3 and s[0] == "b" and s[1] in "'\"" and s[-1] == s[1]:
        s = s[2:-1]
    return s


def orbit_from_tropomi_path(path: str) -> str | None:
    m = _ORBIT_RE.search(Path(path).name)
    return m.group(1) if m else None


def stamp_from_pace_path(path: str) -> str | None:
    m = _PACE_STAMP_RE.search(Path(path).name)
    return m.group(1) if m else None


def index_fs_ret_dir(trop_dir: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for p in Path(trop_dir).glob("fs_ret_*_b5_v10.h5"):
        m = _FS_RET_RE.search(p.name)
        if m:
            out[m.group(1)] = p
    return out


def pace_ret_path(pace_dir: Path, stamp: str) -> Path | None:
    """Resolve PACE SVD retrieval NC for a granule time stamp.

    Tries dated ``YYYY/MM/DD/`` then flat; prefers ``full_parallel`` names.
    """
    if len(stamp) < 8:
        return None
    pace_dir = Path(pace_dir)
    y, m, d = stamp[:4], stamp[4:6], stamp[6:8]
    candidates: list[Path] = []
    day_dir = pace_dir / y / m / d
    if day_dir.is_dir():
        candidates.extend(sorted(day_dir.glob(f"interim_{stamp}_svd_retrieval*.nc")))
    candidates.extend(sorted(pace_dir.glob(f"interim_{stamp}_svd_retrieval*.nc")))
    if not candidates:
        return None
    preferred = [p for p in candidates if "full_parallel" in p.name]
    return preferred[0] if preferred else candidates[0]


def _nc_float_array(var) -> np.ndarray:
    """Read a netCDF variable as float32 with masked/fill values → NaN.

    ``np.asarray(masked)`` densifies CF fills (e.g. 9.96921e36) into ordinary
    floats that still pass ``np.isfinite``, so unfilled retrieval pixels look
    like huge valid SIF. Always unmask into NaN and scrub classic fills.
    """
    arr = np.ma.array(var[:], copy=False)
    out = np.ma.filled(arr, np.nan).astype(np.float32, copy=False)
    out = np.where(np.isfinite(out), out, np.nan).astype(np.float32, copy=False)
    # Classic netCDF4 default float fill, in case auto-mask was off
    out[np.abs(out) > 1.0e35] = np.nan
    return out


def _nc_uint8_array(var) -> np.ndarray:
    """Read uint8 flag; CF fill 255 → 0 (not converged / unknown)."""
    arr = np.ma.array(var[:], copy=False)
    out = np.ma.filled(arr, 0).astype(np.uint8, copy=False)
    out[out == 255] = 0
    return out


def fill_tropomi_fs_ret(
    orbit_to_file: dict[str, Path],
    orbits: np.ndarray,
    trop_pix: np.ndarray,
    trop_scan: np.ndarray,
    sif_out: np.ndarray,
    chi2_out: np.ndarray,
    found_out: np.ndarray,
) -> dict:
    by_orbit: dict[str, list[int]] = defaultdict(list)
    for i, orb in enumerate(orbits):
        if orb:
            by_orbit[str(orb)].append(i)

    n_file_miss = 0
    n_pix_miss = 0
    n_hit = 0
    for orb, idxs in by_orbit.items():
        path = orbit_to_file.get(orb)
        if path is None or not path.is_file():
            n_file_miss += len(idxs)
            continue
        with h5py.File(path, "r") as f:
            det = np.asarray(f["RETRIEVAL_RESULT/detector_pixel"][:], dtype=np.int32)
            scn = np.asarray(f["RETRIEVAL_RESULT/scanline"][:], dtype=np.int32)
            sif = np.asarray(f["RETRIEVAL_RESULT/sif"][:], dtype=np.float32)
            chi2 = np.asarray(f["RETRIEVAL_RESULT/chi2"][:], dtype=np.float32)
        key = det.astype(np.int64) * 100_000 + scn.astype(np.int64)
        lut = {int(k): j for j, k in enumerate(key)}
        for i in idxs:
            k = int(trop_pix[i] + 1) * 100_000 + int(trop_scan[i] + 1)
            j = lut.get(k)
            if j is None:
                n_pix_miss += 1
                continue
            sif_out[i] = sif[j]
            chi2_out[i] = chi2[j]
            found_out[i] = 1
            n_hit += 1
    return {
        "n_hit": n_hit,
        "n_file_miss": n_file_miss,
        "n_pix_miss": n_pix_miss,
        "n_orbits_used": len(by_orbit),
        "n_orbit_files": sum(1 for o in by_orbit if o in orbit_to_file),
    }


def fill_pace_svd(
    pace_dir: Path,
    stamps: np.ndarray,
    pace_pix: np.ndarray,
    pace_scan: np.ndarray,
    sif_out: np.ndarray,
    chi2_out: np.ndarray,
    converged_out: np.ndarray,
    found_out: np.ndarray,
) -> dict:
    by_stamp: dict[str, list[int]] = defaultdict(list)
    for i, st in enumerate(stamps):
        if st:
            by_stamp[str(st)].append(i)

    n_file_miss = 0
    n_pix_miss = 0
    n_hit = 0
    for st, idxs in by_stamp.items():
        path = pace_ret_path(pace_dir, st)
        if path is None or not path.is_file():
            n_file_miss += len(idxs)
            continue
        with nc.Dataset(path) as ds:
            sp = np.asarray(ds["source_pixel_index"][:], dtype=np.int32)
            ss = np.asarray(ds["source_scan_index"][:], dtype=np.int32)
            sif = _nc_float_array(ds["sif_radiance_678nm"])
            chi2 = _nc_float_array(ds["reduced_chi2"])
            conv = _nc_uint8_array(ds["converged"])
        pix_lut = {int(v): j for j, v in enumerate(sp)}
        scan_lut = {int(v): j for j, v in enumerate(ss)}
        for i in idxs:
            ia = pix_lut.get(int(pace_pix[i] + 1))
            ja = scan_lut.get(int(pace_scan[i] + 1))
            if ia is None or ja is None:
                n_pix_miss += 1
                continue
            val = float(sif[ja, ia])
            if not np.isfinite(val):
                n_pix_miss += 1
                continue
            sif_out[i] = val
            chi2_out[i] = chi2[ja, ia]
            converged_out[i] = conv[ja, ia]
            found_out[i] = 1
            n_hit += 1
    return {
        "n_hit": n_hit,
        "n_file_miss": n_file_miss,
        "n_pix_miss": n_pix_miss,
        "n_stamps_used": len(by_stamp),
        "n_stamp_files": sum(
            1 for s in by_stamp if pace_ret_path(pace_dir, s) is not None
        ),
    }
