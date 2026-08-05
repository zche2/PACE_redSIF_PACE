"""Granule-level time + bbox pairing → unique TROPOMI swaths."""

from __future__ import annotations

from pathlib import Path

import netCDF4 as nc
import numpy as np

from .config import Config, pairing_fingerprint_payload
from . import cache
from .discover import list_pace, list_tropomi
from .geo import BBox, TimeWindow, bbox_from_latlon, bboxes_overlap, expand_bbox, time_close


def read_pace_bbox(path: Path, margin_deg: float) -> BBox:
    with nc.Dataset(path) as ds:
        lat = np.asarray(ds["geolocation_data"]["latitude"][:], dtype=float)
        lon = np.asarray(ds["geolocation_data"]["longitude"][:], dtype=float)
    bb = bbox_from_latlon(lat, lon)
    if bb is None:
        raise RuntimeError(f"No valid OCI geolocation in {path.name}")
    return expand_bbox(bb, margin_deg)


def read_tropomi_chunk_bboxes(path: Path, cfg: Config) -> list[BBox]:
    scan_s = cfg.tropo_bbox_scan_stride
    pix_s = cfg.tropo_bbox_pix_stride
    with nc.Dataset(path) as ds:
        g = ds["BAND5_RADIANCE"]["STANDARD_MODE"]["GEODATA"]
        lat = np.asarray(g["latitude"][0, ::scan_s, ::pix_s], dtype=float)
        lon = np.asarray(g["longitude"][0, ::scan_s, ::pix_s], dtype=float)
    chunk = max(1, cfg.tropo_chunk_scans // scan_s)
    boxes: list[BBox] = []
    for i0 in range(0, lat.shape[0], chunk):
        bb = bbox_from_latlon(lat[i0 : i0 + chunk], lon[i0 : i0 + chunk])
        if bb is not None:
            boxes.append(expand_bbox(bb, cfg.bbox_margin_deg))
    return boxes


def build_swath_list(cfg: Config, log=print) -> list[dict]:
    pace_all = list_pace(cfg)
    tropo_times = list_tropomi(cfg)
    pace_sub = pace_all[:: cfg.pace_stride]
    log(
        f"Discover: PACE={len(pace_all):,} (stride={cfg.pace_stride} → {len(pace_sub)})  "
        f"TROPOMI={len(tropo_times):,}"
    )

    pace_meta: list[tuple[TimeWindow, BBox]] = []
    for i, tw in enumerate(pace_sub, 1):
        bb = read_pace_bbox(tw.path, cfg.bbox_margin_deg)
        pace_meta.append((tw, bb))
        if i == 1 or i == len(pace_sub) or i % 50 == 0:
            log(f"  PACE bbox [{i}/{len(pace_sub)}] {tw.path.name}")

    tropo_chunk_cache: dict[Path, list[BBox]] = {}
    pairs: list[dict] = []
    for pace_tw, pace_bb in pace_meta:
        time_hits = [t for t in tropo_times if time_close(pace_tw, t, cfg.dt_max_min_granule)]
        for tropo_tw in time_hits:
            if tropo_tw.path not in tropo_chunk_cache:
                log(f"  TROPOMI geo orbit={tropo_tw.orbit} {tropo_tw.path.name}")
                tropo_chunk_cache[tropo_tw.path] = read_tropomi_chunk_bboxes(tropo_tw.path, cfg)
            chunks = tropo_chunk_cache[tropo_tw.path]
            if any(bboxes_overlap(pace_bb, c) for c in chunks):
                dt_min = abs((pace_tw.mid - tropo_tw.mid).total_seconds()) / 60.0
                pairs.append(
                    {
                        "pace": pace_tw.path.name,
                        "pace_path": str(pace_tw.path),
                        "tropomi": tropo_tw.path.name,
                        "tropomi_path": str(tropo_tw.path),
                        "orbit": tropo_tw.orbit,
                        "pace_start": pace_tw.start.isoformat(),
                        "tropomi_start": tropo_tw.start.isoformat(),
                        "tropomi_end": tropo_tw.end.isoformat(),
                        "dt_mid_min": round(dt_min, 2),
                    }
                )

    log(f"PACE↔TROPOMI pairs: {len(pairs)}")
    swaths: dict[str, dict] = {}
    for pr in pairs:
        key = pr["tropomi_path"]
        if key not in swaths:
            swaths[key] = {
                "orbit": pr["orbit"],
                "tropomi": pr["tropomi"],
                "tropomi_path": pr["tropomi_path"],
                "tropomi_start": pr["tropomi_start"],
                "tropomi_end": pr["tropomi_end"],
                "n_pace_matches": 0,
                "pace_granules": [],
                "pace_paths": [],
            }
        swaths[key]["n_pace_matches"] += 1
        swaths[key]["pace_granules"].append(pr["pace"])
        swaths[key]["pace_paths"].append(pr["pace_path"])

    swath_list = sorted(swaths.values(), key=lambda s: (s["tropomi_start"], s["orbit"] or 0))
    log(f"Co-located TROPOMI swaths: {len(swath_list)}")
    return swath_list


def run_pairs_stage(
    cfg: Config,
    dirs: dict[str, Path],
    *,
    force: bool,
    log=print,
) -> tuple[list[dict], str]:
    fp_payload = pairing_fingerprint_payload(cfg)
    fp = cache.fingerprint(fp_payload)
    meta_path = dirs["pairs"] / "pairs_meta.json"
    list_path = dirs["pairs"] / "swath_list.json"

    if (not force) and cache.is_fresh(
        meta_path, fp, required_files=[list_path], enabled=cfg.cache.enabled
    ):
        log(f"[pairs] CACHE HIT  fp={fp[:12]}…")
        swath_list = cache.read_json(list_path)
        return swath_list, fp

    log(f"[pairs] CACHE MISS — building swath list")
    swath_list = build_swath_list(cfg, log=log)
    cache.write_json(list_path, swath_list)
    cache.write_json(
        meta_path,
        {
            "fingerprint": fp,
            "n_swaths": len(swath_list),
            "payload": fp_payload,
        },
    )
    return swath_list, fp
