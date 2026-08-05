"""Load and validate pipeline YAML config."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

import yaml


@dataclass
class SVDConfig:
    sza_max_deg: float = 85.0
    vza_max_deg: float = 85.0
    lt_max: float = 80.0
    lt_mask_wl_min: float = 660.0
    wl_min: float = 660.0
    wl_max: float = 720.0
    n_wl: int = 300
    max_n: int = 50_000
    n_pc: int = 10
    mean_center: bool = True
    # subsample stride when collecting training candidates (before max_n)
    sample_stride: int = 1


@dataclass
class CacheConfig:
    enabled: bool = True


@dataclass
class Config:
    pace_dir: Path
    tropomi_dir: Path
    rsr_path: Path
    snr_path: Path
    output_dir: Path

    year: int = 2025
    month: int = 7
    day_start: int = 1
    day_end: int = 7
    tropomi_product: str = "OFFL"  # OFFL | RPRO | BOTH
    pace_stride: int = 1
    pace_duration_min: int = 5

    dt_max_min_granule: float = 60.0
    bbox_margin_deg: float = 0.5
    tropo_bbox_scan_stride: int = 40
    tropo_bbox_pix_stride: int = 20
    tropo_chunk_scans: int = 400

    max_dist_km: float = 2.5
    max_dt_min: float = 30.0
    lt_max_oci: float = 30.0
    lt_mask_wl_min_oci: float = 600.0
    pace_stride_pix: int = 150
    tropo_scan_stride: int = 2
    tropo_pix_stride: int = 2
    r_earth_km: float = 6371.0

    oci_band_min: float = 660.0
    oci_band_max: float = 720.0
    clip_negative_rsr: bool = True
    snr_fpa: str = "Red"
    snr_wl_min: float = 600.0
    snr_wl_max: float = 900.0
    add_noise: bool = True
    noise_seed: int = 0

    svd: SVDConfig = field(default_factory=SVDConfig)
    cache: CacheConfig = field(default_factory=CacheConfig)

    def to_plain_dict(self) -> dict[str, Any]:
        d = asdict(self)
        for k, v in list(d.items()):
            if isinstance(v, Path):
                d[k] = str(v)
        return d


_PATH_KEYS = ("pace_dir", "tropomi_dir", "rsr_path", "snr_path", "output_dir")


def _merge_svd(raw: dict[str, Any] | None) -> SVDConfig:
    if not raw:
        return SVDConfig()
    known = {f.name for f in fields(SVDConfig)}
    return SVDConfig(**{k: v for k, v in raw.items() if k in known})


def load_config(path: Path | str) -> Config:
    path = Path(path)
    raw = yaml.safe_load(path.read_text()) or {}
    if not isinstance(raw, dict):
        raise ValueError(f"Config root must be a mapping: {path}")

    missing = [k for k in _PATH_KEYS if k not in raw]
    if missing:
        raise ValueError(f"Missing required config keys: {missing}")

    svd = _merge_svd(raw.pop("svd", None))
    cache_raw = raw.pop("cache", None) or {}
    cache = CacheConfig(**{k: v for k, v in cache_raw.items() if k in {"enabled"}})

    known = {f.name for f in fields(Config)} - {"svd", "cache"}
    kwargs: dict[str, Any] = {}
    for k, v in raw.items():
        if k not in known:
            continue
        if k in _PATH_KEYS:
            kwargs[k] = Path(v)
        else:
            kwargs[k] = v

    cfg = Config(**kwargs, svd=svd, cache=cache)
    if cfg.tropomi_product not in {"OFFL", "RPRO", "BOTH"}:
        raise ValueError(f"tropomi_product must be OFFL|RPRO|BOTH, got {cfg.tropomi_product}")
    return cfg


def pairing_fingerprint_payload(cfg: Config) -> dict[str, Any]:
    return {
        "pace_dir": str(cfg.pace_dir.resolve()),
        "tropomi_dir": str(cfg.tropomi_dir.resolve()),
        "year": cfg.year,
        "month": cfg.month,
        "day_start": cfg.day_start,
        "day_end": cfg.day_end,
        "tropomi_product": cfg.tropomi_product,
        "pace_stride": cfg.pace_stride,
        "pace_duration_min": cfg.pace_duration_min,
        "dt_max_min_granule": cfg.dt_max_min_granule,
        "bbox_margin_deg": cfg.bbox_margin_deg,
        "tropo_bbox_scan_stride": cfg.tropo_bbox_scan_stride,
        "tropo_bbox_pix_stride": cfg.tropo_bbox_pix_stride,
        "tropo_chunk_scans": cfg.tropo_chunk_scans,
    }


def matches_fingerprint_payload(cfg: Config) -> dict[str, Any]:
    return {
        "max_dist_km": cfg.max_dist_km,
        "max_dt_min": cfg.max_dt_min,
        "lt_max_oci": cfg.lt_max_oci,
        "lt_mask_wl_min_oci": cfg.lt_mask_wl_min_oci,
        "pace_stride_pix": cfg.pace_stride_pix,
        "tropo_scan_stride": cfg.tropo_scan_stride,
        "tropo_pix_stride": cfg.tropo_pix_stride,
        "r_earth_km": cfg.r_earth_km,
    }


def svd_fingerprint_payload(cfg: Config) -> dict[str, Any]:
    s = cfg.svd
    return {
        "sza_max_deg": s.sza_max_deg,
        "vza_max_deg": s.vza_max_deg,
        "lt_max": s.lt_max,
        "lt_mask_wl_min": s.lt_mask_wl_min,
        "wl_min": s.wl_min,
        "wl_max": s.wl_max,
        "n_wl": s.n_wl,
        "max_n": s.max_n,
        "mean_center": s.mean_center,
        "sample_stride": s.sample_stride,
        "tropo_scan_stride": cfg.tropo_scan_stride,
        "tropo_pix_stride": cfg.tropo_pix_stride,
    }


def rsr_fingerprint_payload(cfg: Config) -> dict[str, Any]:
    return {
        "rsr_path": str(cfg.rsr_path.resolve()),
        "oci_band_min": cfg.oci_band_min,
        "oci_band_max": cfg.oci_band_max,
        "clip_negative_rsr": cfg.clip_negative_rsr,
    }


def simulate_fingerprint_payload(cfg: Config) -> dict[str, Any]:
    return {
        "n_pc": cfg.svd.n_pc,
        "add_noise": cfg.add_noise,
        "noise_seed": cfg.noise_seed,
        "snr_path": str(cfg.snr_path.resolve()),
        "snr_fpa": cfg.snr_fpa,
        "snr_wl_min": cfg.snr_wl_min,
        "snr_wl_max": cfg.snr_wl_max,
    }
