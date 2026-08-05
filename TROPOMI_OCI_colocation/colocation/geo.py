"""Geometry helpers shared by discovery / pairing / matching."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np


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


def lonlat_to_xyz(lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    lon_r = np.deg2rad(lon)
    lat_r = np.deg2rad(lat)
    cl = np.cos(lat_r)
    return np.column_stack([cl * np.cos(lon_r), cl * np.sin(lon_r), np.sin(lat_r)])
