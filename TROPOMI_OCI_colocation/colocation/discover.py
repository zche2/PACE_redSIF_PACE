"""Discover PACE OCI L1B and TROPOMI BD5 granules from filenames.

Layouts
-------
PACE:    ``pace_dir/YYYY/MM/DD/PACE_OCI.*.L1B.V3.nc`` (falls back to flat root)
TROPOMI: ``tropomi_dir/YYYY/MM/S5P_*_L1B_RA_BD5_*.nc`` (falls back to rglob)
"""

from __future__ import annotations

import re
from datetime import datetime, timedelta, timezone
from pathlib import Path

from .config import Config
from .geo import TimeWindow

PACE_NAME_RE = re.compile(r"PACE_OCI\.(\d{8}T\d{6})\.L1B\.V3\.nc$")
TROPO_NAME_RE = re.compile(
    r"S5P_(OFFL|RPRO)_L1B_RA_BD5_(\d{8}T\d{6})_(\d{8}T\d{6})_(\d+)_.+\.nc$"
)


def list_pace(cfg: Config) -> list[TimeWindow]:
    window_start = datetime(cfg.year, cfg.month, cfg.day_start, tzinfo=timezone.utc)
    window_end = datetime(cfg.year, cfg.month, cfg.day_end, 23, 59, 59, tzinfo=timezone.utc)
    out: list[TimeWindow] = []

    def _add(p: Path) -> None:
        m = PACE_NAME_RE.match(p.name)
        if not m:
            return
        start = datetime.strptime(m.group(1), "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
        end = start + timedelta(minutes=cfg.pace_duration_min)
        if start <= window_end and end >= window_start:
            out.append(TimeWindow(p, start, end))

    found_day_layout = False
    for day in range(cfg.day_start, cfg.day_end + 1):
        day_dir = cfg.pace_dir / f"{cfg.year:04d}" / f"{cfg.month:02d}" / f"{day:02d}"
        if not day_dir.is_dir():
            continue
        found_day_layout = True
        for p in sorted(day_dir.glob("PACE_OCI.*.L1B.V3.nc")):
            _add(p)

    if not found_day_layout:
        pattern = f"PACE_OCI.{cfg.year}{cfg.month:02d}*.L1B.V3.nc"
        for p in sorted(cfg.pace_dir.glob(pattern)):
            _add(p)

    return out


def list_tropomi(cfg: Config) -> list[TimeWindow]:
    window_start = datetime(cfg.year, cfg.month, cfg.day_start, tzinfo=timezone.utc)
    window_end = datetime(cfg.year, cfg.month, cfg.day_end, 23, 59, 59, tzinfo=timezone.utc)
    out: list[TimeWindow] = []

    def _consume(paths: list[Path]) -> None:
        for p in paths:
            m = TROPO_NAME_RE.match(p.name)
            if not m:
                continue
            prod, t0, t1, orbit = m.groups()
            if cfg.tropomi_product != "BOTH" and prod != cfg.tropomi_product:
                continue
            start = datetime.strptime(t0, "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
            end = datetime.strptime(t1, "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
            if start <= window_end and end >= window_start:
                out.append(TimeWindow(p, start, end, orbit=int(orbit)))

    month_dir = cfg.tropomi_dir / f"{cfg.year:04d}" / f"{cfg.month:02d}"
    search_dirs = [month_dir]
    if cfg.day_start <= 2:
        if cfg.month == 1:
            search_dirs.insert(0, cfg.tropomi_dir / f"{cfg.year - 1:04d}" / "12")
        else:
            search_dirs.insert(
                0, cfg.tropomi_dir / f"{cfg.year:04d}" / f"{cfg.month - 1:02d}"
            )

    found_month_layout = False
    for d in search_dirs:
        if not d.is_dir():
            continue
        found_month_layout = True
        _consume(sorted(d.glob("S5P_*_L1B_RA_BD5_*.nc")))

    if not found_month_layout:
        _consume(sorted(cfg.tropomi_dir.rglob("S5P_*_L1B_RA_BD5_*.nc")))

    seen: set[Path] = set()
    uniq: list[TimeWindow] = []
    for tw in out:
        if tw.path in seen:
            continue
        seen.add(tw.path)
        uniq.append(tw)
    uniq.sort(key=lambda t: (t.start, t.orbit or 0))
    return uniq
