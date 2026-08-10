#!/usr/bin/env python3
"""Join independent TROPOMI / PACE SIF retrievals onto colocation matchups.

Index convention (validated against lat/lon):
  TROPOMI fs_ret_*.h5 : detector_pixel == trop_pix+1, scanline == trop_scan+1
  PACE svd NC         : source_pixel_index == pace_pix+1, source_scan_index == pace_scan+1

Examples
--------
  python join_independent_sif.py
  python join_independent_sif.py --max-matches 500 --dry-run-validate
"""

from __future__ import annotations

import argparse
import json
import os
import re
import time
from collections import defaultdict
from pathlib import Path

import h5py
import netCDF4 as nc
import numpy as np
import xarray as xr

_ORBIT_RE = re.compile(r"_(\d{5})_\d{2}_")
_PACE_STAMP_RE = re.compile(r"PACE_OCI\.(\d{8}T\d{6})")
_FS_RET_RE = re.compile(r"fs_ret_(\d{5})_b5_v10\.h5$")


def _path_str(v) -> str:
    if isinstance(v, (bytes, np.bytes_)):
        return v.decode("utf-8", errors="replace").strip("\x00")
    s = str(v)
    # defend against numpy stringification of bytes: "b'/path'"
    if len(s) >= 3 and s[0] == "b" and s[1] in "'\"" and s[-1] == s[1]:
        s = s[2:-1]
    return s


def _orbit_from_tropomi_path(path: str) -> str | None:
    m = _ORBIT_RE.search(Path(path).name)
    return m.group(1) if m else None


def _stamp_from_pace_path(path: str) -> str | None:
    m = _PACE_STAMP_RE.search(Path(path).name)
    return m.group(1) if m else None


def _index_fs_ret_dir(trop_dir: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for p in trop_dir.glob("fs_ret_*_b5_v10.h5"):
        m = _FS_RET_RE.search(p.name)
        if m:
            out[m.group(1)] = p
    return out


def _pace_ret_path(pace_dir: Path, stamp: str) -> Path:
    return pace_dir / f"interim_{stamp}_svd_retrieval_full_parallel.nc"


def _fill_tropomi(
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
            by_orbit[orb].append(i)

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
        # (detector_pixel, scanline) -> flat index
        key = det.astype(np.int64) * 100_000 + scn.astype(np.int64)
        # last wins if duplicates
        lut = {int(k): j for j, k in enumerate(key)}
        for i in idxs:
            # matchup stores 0-based-ish indices; product uses +1 (validated vs lat/lon)
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


def _fill_pace(
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
            by_stamp[st].append(i)

    n_file_miss = 0
    n_pix_miss = 0
    n_hit = 0
    for st, idxs in by_stamp.items():
        path = _pace_ret_path(pace_dir, st)
        if not path.is_file():
            n_file_miss += len(idxs)
            continue
        with nc.Dataset(path) as ds:
            sp = np.asarray(ds["source_pixel_index"][:], dtype=np.int32)
            ss = np.asarray(ds["source_scan_index"][:], dtype=np.int32)
            sif = np.asarray(ds["sif_radiance_678nm"][:], dtype=np.float32)  # (scans, pixels)
            chi2 = np.asarray(ds["reduced_chi2"][:], dtype=np.float32)
            conv = np.asarray(ds["converged"][:], dtype=np.uint8)
        pix_lut = {int(v): j for j, v in enumerate(sp)}
        scan_lut = {int(v): j for j, v in enumerate(ss)}
        for i in idxs:
            # matchup pace_pix/scan + 1 == source_* (1-based L1)
            ia = pix_lut.get(int(pace_pix[i] + 1))
            ja = scan_lut.get(int(pace_scan[i] + 1))
            if ia is None or ja is None:
                n_pix_miss += 1
                continue
            val = sif[ja, ia]
            sif_out[i] = val
            chi2_out[i] = chi2[ja, ia]
            converged_out[i] = conv[ja, ia]
            found_out[i] = 1 if np.isfinite(val) else 0
            if found_out[i]:
                n_hit += 1
            else:
                n_pix_miss += 1
    return {
        "n_hit": n_hit,
        "n_file_miss": n_file_miss,
        "n_pix_miss": n_pix_miss,
        "n_stamps_used": len(by_stamp),
        "n_stamp_files": sum(1 for s in by_stamp if _pace_ret_path(pace_dir, s).is_file()),
    }


def _validate_offsets(
    matchup_nc: Path,
    trop_dir: Path,
    pace_dir: Path,
    n_check: int = 20,
) -> None:
    """Print lat/lon agreement for +1 index convention on a few hits."""
    orbit_to_file = _index_fs_ret_dir(trop_dir)
    ds = xr.open_dataset(matchup_nc)
    n = min(int(ds.sizes["match"]), 5000)
    print(f"[validate] checking up to {n_check} hits in first {n} matches")
    n_ok_t = n_ok_p = 0
    for i in range(n):
        if n_ok_t >= n_check and n_ok_p >= n_check:
            break
        if n_ok_t < n_check:
            orb = _orbit_from_tropomi_path(_path_str(ds["tropomi_path"].values[i]))
            h5 = orbit_to_file.get(orb or "")
            if h5 and h5.is_file():
                pix = int(ds["trop_pix"].values[i]) + 1
                scan = int(ds["trop_scan"].values[i]) + 1
                with h5py.File(h5, "r") as f:
                    det = f["RETRIEVAL_RESULT/detector_pixel"][:]
                    scn = f["RETRIEVAL_RESULT/scanline"][:]
                    la = f["RETRIEVAL_RESULT/latitude"][:]
                    lo = f["RETRIEVAL_RESULT/longitude"][:]
                    m = (det == pix) & (scn == scan)
                    if m.any():
                        j = int(np.flatnonzero(m)[0])
                        dlat = abs(float(la[j]) - float(ds["lat_trop"].values[i]))
                        dlon = abs(float(lo[j]) - float(ds["lon_trop"].values[i]))
                        print(f"  trop i={i}: dlat={dlat:.6f} dlon={dlon:.6f}")
                        n_ok_t += 1
        if n_ok_p < n_check:
            st = _stamp_from_pace_path(_path_str(ds["pace_path"].values[i]))
            if st:
                p = _pace_ret_path(pace_dir, st)
                if p.is_file():
                    pp = int(ds["pace_pix"].values[i]) + 1
                    ps = int(ds["pace_scan"].values[i]) + 1
                    with nc.Dataset(p) as d:
                        sp = d["source_pixel_index"][:]
                        ss = d["source_scan_index"][:]
                        la = d["latitude"][:]
                        lo = d["longitude"][:]
                        ii = np.where(sp == pp)[0]
                        jj = np.where(ss == ps)[0]
                        if ii.size and jj.size:
                            ia, ja = int(ii[0]), int(jj[0])
                            dlat = abs(float(la[ja, ia]) - float(ds["lat_pace"].values[i]))
                            dlon = abs(float(lo[ja, ia]) - float(ds["lon_pace"].values[i]))
                            print(f"  pace i={i}: dlat={dlat:.6f} dlon={dlon:.6f}")
                            n_ok_p += 1
    ds.close()
    print(f"[validate] trop hits shown={n_ok_t} pace hits shown={n_ok_p}")


def join(
    matchup_nc: Path,
    trop_dir: Path,
    pace_dir: Path,
    output_nc: Path,
    max_matches: int = 0,
) -> dict:
    ds = xr.open_dataset(matchup_nc)
    n_all = int(ds.sizes["match"])
    n = n_all if max_matches <= 0 else min(n_all, max_matches)

    trop_pix = np.asarray(ds["trop_pix"].values[:n], dtype=np.int32)
    trop_scan = np.asarray(ds["trop_scan"].values[:n], dtype=np.int32)
    pace_pix = np.asarray(ds["pace_pix"].values[:n], dtype=np.int32)
    pace_scan = np.asarray(ds["pace_scan"].values[:n], dtype=np.int32)
    trop_paths = [_path_str(v) for v in ds["tropomi_path"].values[:n]]
    pace_paths = [_path_str(v) for v in ds["pace_path"].values[:n]]
    orbits = np.asarray([_orbit_from_tropomi_path(p) or "" for p in trop_paths])
    stamps = np.asarray([_stamp_from_pace_path(p) or "" for p in pace_paths])

    sif_trop = np.full(n, np.nan, dtype=np.float32)
    chi2_trop = np.full(n, np.nan, dtype=np.float32)
    found_trop = np.zeros(n, dtype=np.uint8)
    sif_pace = np.full(n, np.nan, dtype=np.float32)
    chi2_pace = np.full(n, np.nan, dtype=np.float32)
    converged_pace = np.zeros(n, dtype=np.uint8)
    found_pace = np.zeros(n, dtype=np.uint8)

    orbit_to_file = _index_fs_ret_dir(trop_dir)
    t0 = time.time()
    stats_t = _fill_tropomi(
        orbit_to_file, orbits, trop_pix, trop_scan, sif_trop, chi2_trop, found_trop
    )
    t1 = time.time()
    stats_p = _fill_pace(
        pace_dir,
        stamps,
        pace_pix,
        pace_scan,
        sif_pace,
        chi2_pace,
        converged_pace,
        found_pace,
    )
    t2 = time.time()

    output_nc.parent.mkdir(parents=True, exist_ok=True)
    tmp_nc = output_nc.with_suffix(output_nc.suffix + ".tmp")
    if tmp_nc.exists():
        tmp_nc.unlink()
    with nc.Dataset(tmp_nc, "w") as out:
        out.createDimension("match", n)
        def _v(name, data, dtype, **attrs):
            var = out.createVariable(name, dtype, ("match",), zlib=True, complevel=4)
            var[:] = data
            for k, v in attrs.items():
                setattr(var, k, v)
            return var

        _v(
            "sif_trop",
            sif_trop,
            "f4",
            long_name="Independent TROPOMI BD5 SIF (RETRIEVAL_RESULT/sif)",
            source_dir=str(trop_dir),
            index_note="joined with detector_pixel=trop_pix+1, scanline=trop_scan+1",
        )
        _v(
            "chi2_trop",
            chi2_trop,
            "f4",
            long_name="Independent TROPOMI BD5 retrieval chi2",
            source_dir=str(trop_dir),
        )
        _v("trop_found", found_trop, "u1", long_name="1 if TROPOMI independent SIF found")
        _v(
            "sif_pace_678nm",
            sif_pace,
            "f4",
            long_name="Independent PACE SVD sif_radiance_678nm",
            units="W m-2 sr-1 um-1",
            source_dir=str(pace_dir),
            index_note="joined with source_pixel_index=pace_pix+1, source_scan_index=pace_scan+1",
        )
        _v(
            "chi2_pace",
            chi2_pace,
            "f4",
            long_name="Independent PACE SVD reduced_chi2",
            source_dir=str(pace_dir),
        )
        _v(
            "converged_pace",
            converged_pace,
            "u1",
            long_name="Independent PACE SVD converged (1=yes)",
            source_dir=str(pace_dir),
        )
        _v("pace_found", found_pace, "u1", long_name="1 if PACE independent SIF found")

        # helpful copies for downstream without reopening matchup
        for name in ("lat_trop", "lon_trop", "lat_pace", "lon_pace", "trop_pix", "trop_scan", "pace_pix", "pace_scan"):
            if name in ds:
                _v(name, np.asarray(ds[name].values[:n]), "f4" if "lat" in name or "lon" in name else "i4")

        out.title = "Independent TROPOMI & PACE SIF joined to OCI colocation matchups"
        out.matchup_nc = str(matchup_nc)
        out.tropomi_retrieval_dir = str(trop_dir)
        out.pace_retrieval_dir = str(pace_dir)
        out.n_match = np.int32(n)
        out.history = time.strftime("%Y-%m-%dT%H:%M:%S") + " join_independent_sif.py"

    ds.close()
    try:
        os.replace(tmp_nc, output_nc)
    except OSError as exc:
        # File may be open in a notebook kernel — write beside it instead
        alt = output_nc.with_name(output_nc.stem + "_new" + output_nc.suffix)
        os.replace(tmp_nc, alt)
        print(f"[warn] could not replace {output_nc} ({exc}); wrote {alt}")
        output_nc = alt

    meta = {
        "matchup_nc": str(matchup_nc),
        "output_nc": str(output_nc),
        "n_match": n,
        "n_match_all": n_all,
        "tropomi": {**stats_t, "dir": str(trop_dir), "seconds": round(t1 - t0, 2)},
        "pace": {**stats_p, "dir": str(pace_dir), "seconds": round(t2 - t1, 2)},
        "index_convention": {
            "tropomi": "detector_pixel=trop_pix+1, scanline=trop_scan+1",
            "pace": "source_pixel_index=pace_pix+1, source_scan_index=pace_scan+1",
        },
    }
    meta_path = output_nc.with_suffix(".json")
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--matchup-nc",
        type=Path,
        default=Path(
            "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/TROPOMI_OCI_colocation_data/"
            "runs/july2025_d01-31/products/matchup.nc"
        ),
    )
    p.add_argument(
        "--trop-dir",
        type=Path,
        default=Path("/home/zhe2/data/MyProjects/TROPOMI_SIF/results/b5"),
    )
    p.add_argument(
        "--pace-dir",
        type=Path,
        default=Path("/home/zhe2/data/PACE/svd_retrieval_output"),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/TROPOMI_OCI_colocation_data/"
            "runs/july2025_d01-31/matches"
        ),
    )
    p.add_argument(
        "--output-name",
        default="matchup_independent_sif.nc",
        help="NetCDF filename under --output-dir",
    )
    p.add_argument("--max-matches", type=int, default=0, help="0 = all matches")
    p.add_argument(
        "--dry-run-validate",
        action="store_true",
        help="Only check lat/lon agreement for +1 index convention, then exit",
    )
    args = p.parse_args(argv)

    if args.dry_run_validate:
        _validate_offsets(args.matchup_nc, args.trop_dir, args.pace_dir)
        return 0

    out = args.output_dir / args.output_name
    meta = join(args.matchup_nc, args.trop_dir, args.pace_dir, out, max_matches=args.max_matches)
    print(json.dumps(meta, indent=2))
    print(f"wrote {out}")
    print(f"wrote {out.with_suffix('.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
