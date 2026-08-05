"""Discover PACE OCI L1B and TROPOMI BD5 granules from filenames."""

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
    pattern = f"PACE_OCI.{cfg.year}{cfg.month:02d}*.L1B.V3.nc"
    for p in sorted(cfg.pace_dir.glob(pattern)):
        m = PACE_NAME_RE.match(p.name)
        if not m:
            continue
        start = datetime.strptime(m.group(1), "%Y%m%dT%H%M%S").replace(tzinfo=timezone.utc)
        end = start + timedelta(minutes=cfg.pace_duration_min)
        if start <= window_end and end >= window_start:
            out.append(TimeWindow(p, start, end))
    return out


def list_tropomi(cfg: Config) -> list[TimeWindow]:
    window_start = datetime(cfg.year, cfg.month, cfg.day_start, tzinfo=timezone.utc)
    window_end = datetime(cfg.year, cfg.month, cfg.day_end, 23, 59, 59, tzinfo=timezone.utc)
    out: list[TimeWindow] = []
    for p in sorted(cfg.tropomi_dir.rglob("S5P_*_L1B_RA_BD5_*.nc")):
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
    return out
