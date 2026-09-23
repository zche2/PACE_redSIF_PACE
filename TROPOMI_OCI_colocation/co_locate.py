#!/usr/bin/env python3
"""TROPOMI BD5 ↔ PACE OCI L1B time+bbox co-location (from co-locate_cleanup.ipynb).

Step 1 only: granule pairing by filename time + bounding-box overlap, then write
JSON of unique TROPOMI swaths with matched PACE granules.

Examples
--------
  # Jan 2025 (full month) → notebooks/co-location_results_202501.json
  python TROPOMI_OCI_colocation/co_locate.py --year 2025 --month 1 --day-end 31 \
      --output TROPOMI_OCI_colocation/notebooks/co-location_results_202501.json

  # Quick smoke (first 7 days)
  python TROPOMI_OCI_colocation/co_locate.py --month 1 --day-end 7 --output /tmp/coloc_smoke.json

Roots (dated layout):
  PACE:    /kiwi-data/Data/satellite/PACE_OCI/L1B_V3/YYYY/MM/DD/
  TROPOMI: /net/squid/data1/projects/TROPOMI/ESA/L1/YYYY/MM/
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np

_HERE = Path(__file__).resolve().parent
_DEFAULT_OUT = _HERE / "notebooks" / "co-location_results.json"


# ── geometry / time helpers ───────────────────────────────────────────────────

@dataclass(frozen=True)
class TimeWindow:
    path: Path
    start: datetime
    end: datetime
    orbit: int | None = None

    @property
    def mid(self) -> datetime:
        return self.start + (self.end - self.start) / 2


@dataclass(frozen=True)
class BBox:
    lat_min: float
    lat_max: float
    lon_min: float
    lon_max: float  # if lon_min > lon_max, box crosses antimeridian


def expand_bbox(b: BBox, margin: float) -> BBox:
    return BBox(
        lat_min=max(-90.0, b.lat_min - margin),
        lat_max=min(90.0, b.lat_max + margin),
        lon_min=((b.lon_min - margin + 180) % 360) - 180,
        lon_max=((b.lon_max + margin + 180) % 360) - 180,
    )


def lon_intervals_overlap(a0: float, a1: float, b0: float, b1: float) -> bool:
    def parts(lo: float, hi: float) -> list[tuple[float, float]]:
        if lo <= hi:
            return [(lo, hi)]
        return [(lo, 180.0), (-180.0, hi)]

    for x0, x1 in parts(a0, a1):
        for y0, y1 in parts(b0, b1):
            if x0 <= y1 and y0 <= x1:
                return True
    return False


def bboxes_overlap(a: BBox, b: BBox) -> bool:
    if a.lat_max < b.lat_min or b.lat_max < a.lat_min:
        return False
    return lon_intervals_overlap(a.lon_min, a.lon_max, b.lon_min, b.lon_max)


def bbox_from_latlon(lat: np.ndarray, lon: np.ndarray) -> BBox | None:
    m = np.isfinite(lat) & np.isfinite(lon)
    if not m.any():
        return None
    lat_v, lon_v = lat[m], lon[m]
    lon_rad = np.deg2rad(lon_v)
    mean = np.angle(np.mean(np.exp(1j * lon_rad)))
    centered = (lon_rad - mean + np.pi) % (2 * np.pi) - np.pi
    lo = ((np.rad2deg(mean + centered.min()) + 180) % 360) - 180
    hi = ((np.rad2deg(mean + centered.max()) + 180) % 360) - 180
    return BBox(float(lat_v.min()), float(lat_v.max()), float(lo), float(hi))


def time_close(a: TimeWindow, b: TimeWindow, dt_max_min: float) -> bool:
    if a.start <= b.end and b.start <= a.end:
        return True
    return abs((a.mid - b.mid).total_seconds()) / 60.0 <= dt_max_min


# ── discover ──────────────────────────────────────────────────────────────────
# PACE L1B root: .../L1B_V3/YYYY/MM/DD/PACE_OCI.*.nc
# TROPOMI BD5:   .../ESA/L1/YYYY/MM/S5P_*_L1B_RA_BD5_*.nc

PACE_NAME_RE = re.compile(r"PACE_OCI\.(\d{8}T\d{6})\.L1B\.V3\.nc$")
TROPO_NAME_RE = re.compile(
    r"S5P_(OFFL|RPRO)_L1B_RA_BD5_(\d{8}T\d{6})_(\d{8}T\d{6})_(\d+)_.+\.nc$"
)

_DEFAULT_PACE_DIR = Path("/kiwi-data/Data/satellite/PACE_OCI/L1B_V3")
_DEFAULT_TROPOMI_DIR = Path("/net/squid/data1/projects/TROPOMI/ESA/L1")


def list_pace(
    pace_dir: Path,
    *,
    year: int,
    month: int,
    day_start: int,
    day_end: int,
    pace_duration_min: float,
) -> list[TimeWindow]:
    """Discover PACE L1B under ``pace_dir/YYYY/MM/DD/`` (falls back to flat root)."""
    window_start = datetime(year, month, day_start, tzinfo=timezone.utc)
    window_end = datetime(year, month, day_end, 23, 59, 59, tzinfo=timezone.utc)
    out: list[TimeWindow] = []

    def _add(p: Path) -> None:
        m = PACE_NAME_RE.match(p.name)
        if not m:
            return
        start = datetime.strptime(m.group(1), "%Y%m%dT%H%M%S").replace(
            tzinfo=timezone.utc
        )
        end = start + timedelta(minutes=pace_duration_min)
        if start <= window_end and end >= window_start:
            out.append(TimeWindow(p, start, end))

    found_day_layout = False
    for day in range(day_start, day_end + 1):
        day_dir = pace_dir / f"{year:04d}" / f"{month:02d}" / f"{day:02d}"
        if not day_dir.is_dir():
            continue
        found_day_layout = True
        for p in sorted(day_dir.glob("PACE_OCI.*.L1B.V3.nc")):
            _add(p)

    if not found_day_layout:
        # Legacy flat: all granules directly under pace_dir
        for p in sorted(pace_dir.glob(f"PACE_OCI.{year}{month:02d}*.L1B.V3.nc")):
            _add(p)

    return out


def list_tropomi(
    tropomi_dir: Path,
    *,
    year: int,
    month: int,
    day_start: int,
    day_end: int,
    tropomi_product: str,
) -> list[TimeWindow]:
    """Discover TROPOMI BD5 under ``tropomi_dir/YYYY/MM/`` (falls back to rglob)."""
    window_start = datetime(year, month, day_start, tzinfo=timezone.utc)
    window_end = datetime(year, month, day_end, 23, 59, 59, tzinfo=timezone.utc)
    out: list[TimeWindow] = []

    def _consume(paths: list[Path]) -> None:
        for p in paths:
            m = TROPO_NAME_RE.match(p.name)
            if not m:
                continue
            prod, t0, t1, orbit = m.groups()
            if tropomi_product != "BOTH" and prod != tropomi_product:
                continue
            start = datetime.strptime(t0, "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
            end = datetime.strptime(t1, "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
            if start <= window_end and end >= window_start:
                out.append(TimeWindow(p, start, end, orbit=int(orbit)))

    month_dir = tropomi_dir / f"{year:04d}" / f"{month:02d}"
    search_dirs = [month_dir]
    # Orbits that start late on the previous month can still overlap day_start.
    if day_start <= 2:
        if month == 1:
            prev = tropomi_dir / f"{year - 1:04d}" / "12"
        else:
            prev = tropomi_dir / f"{year:04d}" / f"{month - 1:02d}"
        search_dirs.insert(0, prev)

    found_month_layout = False
    for d in search_dirs:
        if not d.is_dir():
            continue
        found_month_layout = True
        _consume(sorted(d.glob("S5P_*_L1B_RA_BD5_*.nc")))

    if not found_month_layout:
        _consume(sorted(tropomi_dir.rglob("S5P_*_L1B_RA_BD5_*.nc")))

    # Deduplicate by path (prev-month + current can overlap listing)
    seen: set[Path] = set()
    uniq: list[TimeWindow] = []
    for tw in out:
        if tw.path in seen:
            continue
        seen.add(tw.path)
        uniq.append(tw)
    uniq.sort(key=lambda t: (t.start, t.orbit or 0))
    return uniq

def read_pace_bbox(path: Path, margin_deg: float) -> BBox:
    with nc.Dataset(path) as ds:
        lat = np.asarray(ds["geolocation_data"]["latitude"][:], dtype=float)
        lon = np.asarray(ds["geolocation_data"]["longitude"][:], dtype=float)
    bb = bbox_from_latlon(lat, lon)
    if bb is None:
        raise RuntimeError(f"No valid OCI geolocation in {path.name}")
    return expand_bbox(bb, margin_deg)


def read_tropomi_chunk_bboxes(
    path: Path,
    *,
    scan_stride: int,
    pix_stride: int,
    chunk_scans: int,
    margin_deg: float,
) -> list[BBox]:
    with nc.Dataset(path) as ds:
        g = ds["BAND5_RADIANCE"]["STANDARD_MODE"]["GEODATA"]
        lat = np.asarray(
            g["latitude"][0, ::scan_stride, ::pix_stride], dtype=float
        )
        lon = np.asarray(
            g["longitude"][0, ::scan_stride, ::pix_stride], dtype=float
        )
    chunk = max(1, chunk_scans // scan_stride)
    boxes: list[BBox] = []
    for i0 in range(0, lat.shape[0], chunk):
        bb = bbox_from_latlon(lat[i0 : i0 + chunk], lon[i0 : i0 + chunk])
        if bb is not None:
            boxes.append(expand_bbox(bb, margin_deg))
    return boxes


# ── pairing ───────────────────────────────────────────────────────────────────

def build_pairs(
    *,
    pace_dir: Path,
    tropomi_dir: Path,
    year: int,
    month: int,
    day_start: int,
    day_end: int,
    tropomi_product: str,
    pace_stride: int,
    pace_duration_min: float,
    dt_max_min: float,
    bbox_margin_deg: float,
    tropo_scan_stride: int,
    tropo_pix_stride: int,
    tropo_chunk_scans: int,
    log=print,
) -> tuple[list[dict], list[dict]]:
    """Return (pairs, swath_list)."""
    pace_all = list_pace(
        pace_dir,
        year=year,
        month=month,
        day_start=day_start,
        day_end=day_end,
        pace_duration_min=pace_duration_min,
    )
    tropo_times = list_tropomi(
        tropomi_dir,
        year=year,
        month=month,
        day_start=day_start,
        day_end=day_end,
        tropomi_product=tropomi_product,
    )
    pace_sub = pace_all[::pace_stride]
    log(
        f"{year}-{month:02d}-{day_start:02d} → {day_end:02d}: "
        f"PACE={len(pace_all):,}  (stride={pace_stride} → {len(pace_sub)})  "
        f"TROPOMI BD5 ({tropomi_product})={len(tropo_times):,}"
    )
    if pace_sub:
        log(f"  first PACE: {pace_sub[0].path.name}")
        log(f"  last  PACE: {pace_sub[-1].path.name}")
    if tropo_times:
        log(f"  first TROPOMI: {tropo_times[0].path.name}")
        log(f"  last  TROPOMI: {tropo_times[-1].path.name}")

    log(f"Reading PACE bboxes for {len(pace_sub)} granules …")
    pace_meta: list[tuple[TimeWindow, BBox]] = []
    for i, tw in enumerate(pace_sub, 1):
        bb = read_pace_bbox(tw.path, bbox_margin_deg)
        pace_meta.append((tw, bb))
        if i == 1 or i == len(pace_sub) or i % 50 == 0:
            log(
                f"  [{i}/{len(pace_sub)}] {tw.path.name}: "
                f"{tw.start:%m-%d %H:%M}Z  "
                f"lat[{bb.lat_min:.1f},{bb.lat_max:.1f}] "
                f"lon[{bb.lon_min:.1f},{bb.lon_max:.1f}]"
            )

    tropo_chunk_cache: dict[Path, list[BBox]] = {}
    pairs: list[dict] = []
    for pace_tw, pace_bb in pace_meta:
        time_hits = [t for t in tropo_times if time_close(pace_tw, t, dt_max_min)]
        for tropo_tw in time_hits:
            if tropo_tw.path not in tropo_chunk_cache:
                log(f"  TROPOMI geo  orbit={tropo_tw.orbit}  {tropo_tw.path.name}")
                tropo_chunk_cache[tropo_tw.path] = read_tropomi_chunk_bboxes(
                    tropo_tw.path,
                    scan_stride=tropo_scan_stride,
                    pix_stride=tropo_pix_stride,
                    chunk_scans=tropo_chunk_scans,
                    margin_deg=bbox_margin_deg,
                )
            chunks = tropo_chunk_cache[tropo_tw.path]
            if any(bboxes_overlap(pace_bb, c) for c in chunks):
                dt_min = abs((pace_tw.mid - tropo_tw.mid).total_seconds()) / 60.0
                pairs.append(
                    {
                        "pace": pace_tw.path.name,
                        "tropomi": tropo_tw.path.name,
                        "orbit": tropo_tw.orbit,
                        "pace_start": pace_tw.start.isoformat(),
                        "tropomi_start": tropo_tw.start.isoformat(),
                        "tropomi_end": tropo_tw.end.isoformat(),
                        "dt_mid_min": round(dt_min, 2),
                        "tropomi_path": str(tropo_tw.path),
                    }
                )

    log(f"PACE granules tested : {len(pace_meta)}")
    log(f"TROPOMI orbits opened: {len(tropo_chunk_cache)}")
    log(f"PACE↔TROPOMI pairs   : {len(pairs)}")

    swaths: dict[str, dict] = {}
    for pr in pairs:
        key = pr["tropomi"]
        if key not in swaths:
            swaths[key] = {
                "orbit": pr["orbit"],
                "tropomi": pr["tropomi"],
                "tropomi_start": pr["tropomi_start"],
                "tropomi_end": pr["tropomi_end"],
                "tropomi_path": pr["tropomi_path"],
                "n_pace_matches": 0,
                "pace_granules": [],
            }
        swaths[key]["n_pace_matches"] += 1
        swaths[key]["pace_granules"].append(pr["pace"])

    swath_list = sorted(swaths.values(), key=lambda s: (s["tropomi_start"], s["orbit"]))
    log(f"Co-located TROPOMI swaths: {len(swath_list)}")
    return pairs, swath_list


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
    return {
        "meta": {
            "created": datetime.now(timezone.utc).isoformat(),
            "description": "TROPOMI BD5 ↔ PACE OCI L1B time+bbox co-location (step 1)",
            "year": year,
            "month": month,
            "day_start": day_start,
            "day_end": day_end,
            "dt_max_min": dt_max_min,
            "bbox_margin_deg": bbox_margin_deg,
            "pace_dir": str(pace_dir),
            "tropomi_dir": str(tropomi_dir),
            "tropomi_product": tropomi_product,
            "n_swaths": len(swath_list),
        },
        "swaths": swath_list,
    }


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="TROPOMI BD5 ↔ PACE OCI L1B time+bbox co-location (cleanup notebook → CLI)"
    )
    p.add_argument("--year", type=int, default=2025)
    p.add_argument("--month", type=int, default=7)
    p.add_argument("--day-start", type=int, default=1)
    p.add_argument("--day-end", type=int, default=31)
    p.add_argument(
        "--pace-dir",
        type=Path,
        default=_DEFAULT_PACE_DIR,
        help="PACE L1B root (granules under YYYY/MM/DD/)",
    )
    p.add_argument(
        "--tropomi-dir",
        type=Path,
        default=_DEFAULT_TROPOMI_DIR,
        help="TROPOMI L1 root (BD5 under YYYY/MM/)",
    )
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

    _, swath_list = build_pairs(
        pace_dir=args.pace_dir,
        tropomi_dir=args.tropomi_dir,
        year=args.year,
        month=args.month,
        day_start=args.day_start,
        day_end=args.day_end,
        tropomi_product=args.tropomi_product,
        pace_stride=args.pace_stride,
        pace_duration_min=args.pace_duration_min,
        dt_max_min=args.dt_max_min,
        bbox_margin_deg=args.bbox_margin_deg,
        tropo_scan_stride=args.tropo_scan_stride,
        tropo_pix_stride=args.tropo_pix_stride,
        tropo_chunk_scans=args.tropo_chunk_scans,
    )

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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
