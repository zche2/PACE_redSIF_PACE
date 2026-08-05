"""Pixel-wise TROPOMI↔OCI matching (OCI dark + KD-tree + dist/Δt)."""

from __future__ import annotations

from pathlib import Path

import netCDF4 as nc
import numpy as np
from netCDF4 import num2date
from scipy.spatial import cKDTree

from . import cache
from .config import Config, matches_fingerprint_payload
from .geo import lonlat_to_xyz


def _datetime64_to_unix(t: np.ndarray) -> np.ndarray:
    """datetime64[ms] → float64 Unix seconds."""
    # np.datetime64 epoch is 1970-01-01
    return t.astype("datetime64[ms]").astype(np.int64).astype(np.float64) / 1000.0


def match_one_swath(
    cfg: Config,
    tropomi_path: Path,
    pace_paths: list[Path],
    *,
    swath_id: int,
    log=print,
) -> dict | None:
    tropomi_path = Path(tropomi_path)
    pace_paths = [Path(p) for p in pace_paths]
    scan_s = cfg.tropo_scan_stride
    pix_s = cfg.tropo_pix_stride
    pace_s = cfg.pace_stride_pix

    with nc.Dataset(tropomi_path) as ds:
        g = ds["BAND5_RADIANCE"]["STANDARD_MODE"]["GEODATA"]
        obs = ds["BAND5_RADIANCE"]["STANDARD_MODE"]["OBSERVATIONS"]
        trop_lat = np.asarray(g["latitude"][0, ::scan_s, ::pix_s], dtype=float)
        trop_lon = np.asarray(g["longitude"][0, ::scan_s, ::pix_s], dtype=float)
        trop_sza = np.asarray(g["solar_zenith_angle"][0, ::scan_s, ::pix_s], dtype=float)
        trop_vza = np.asarray(g["viewing_zenith_angle"][0, ::scan_s, ::pix_s], dtype=float)
        t0 = float(np.asarray(obs["time"][:]).ravel()[0])
        dt_ms = np.asarray(obs["delta_time"][0, :], dtype=float)
        trop_t_units = obs["time"].units

    trop_shape = trop_lat.shape
    trop_scan_time = np.array(
        num2date(t0 + dt_ms / 1000.0, units=trop_t_units, only_use_cftime_datetimes=False),
        dtype="datetime64[ms]",
    )

    trop_ok = np.isfinite(trop_lat) & np.isfinite(trop_lon)
    trop_flat_i = np.flatnonzero(trop_ok)
    trop_xyz = lonlat_to_xyz(trop_lon.ravel()[trop_flat_i], trop_lat.ravel()[trop_flat_i])
    tree = cKDTree(trop_xyz)

    chunks: list[dict] = []
    n_keep_total = 0

    for i_pace, pace_path in enumerate(pace_paths):
        with nc.Dataset(pace_path) as ds:
            geo = ds["geolocation_data"]
            pace_lat = np.asarray(geo["latitude"][::pace_s, ::pace_s], dtype=float)
            pace_lon = np.asarray(geo["longitude"][::pace_s, ::pace_s], dtype=float)
            pace_sza = np.asarray(geo["solar_zenith"][::pace_s, ::pace_s], dtype=float)
            pace_vza = np.asarray(geo["sensor_zenith"][::pace_s, ::pace_s], dtype=float)
            pace_t_var = ds["scan_line_attributes"]["time"]
            pace_t_all = np.asarray(pace_t_var[:], dtype=float)
            pace_t_units = pace_t_var.units

            esd = float(getattr(ds, "earth_sun_distance_correction", 1.0))
            wl = np.asarray(ds["sensor_band_parameters"]["red_wavelength"][:], dtype=float)
            f0 = np.asarray(ds["sensor_band_parameters"]["red_solar_irradiance"][:], dtype=float)
            ib = np.flatnonzero(wl >= cfg.lt_mask_wl_min_oci)
            if ib.size == 0:
                raise RuntimeError(
                    f"No OCI red bands ≥ {cfg.lt_mask_wl_min_oci:g} nm in {pace_path.name}"
                )
            i0, i1 = int(ib[0]), int(ib[-1]) + 1
            rhot = np.asarray(
                ds["observation_data"]["rhot_red"][i0:i1, ::pace_s, ::pace_s],
                dtype=float,
            )
            rhot = np.where(rhot > 0, rhot, np.nan)
            f0_b = f0[i0:i1][:, None, None]
            mu0 = np.cos(np.deg2rad(pace_sza))[None, :, :]
            mu0 = np.where(mu0 > 1e-6, mu0, np.nan)
            lt = rhot * f0_b * mu0 / (np.pi * esd)
            lt_max = np.nanmax(lt, axis=0)
            dark = np.isfinite(lt_max) & (lt_max < cfg.lt_max_oci)

        pace_shape = pace_lat.shape
        pace_scan_time = np.array(
            num2date(pace_t_all, units=pace_t_units, only_use_cftime_datetimes=False),
            dtype="datetime64[ms]",
        )

        pace_ok = np.isfinite(pace_lat) & np.isfinite(pace_lon) & dark
        pace_flat_i = np.flatnonzero(pace_ok)
        if pace_flat_i.size == 0:
            continue

        pace_xyz = lonlat_to_xyz(
            pace_lon.ravel()[pace_flat_i], pace_lat.ravel()[pace_flat_i]
        )
        chord, nn = tree.query(pace_xyz, k=1, workers=-1)
        dist_km = 2.0 * cfg.r_earth_km * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0))
        close_dist = dist_km <= cfg.max_dist_km

        pace_ji_all = np.unravel_index(pace_flat_i, pace_shape)
        trop_ji_nn = np.unravel_index(trop_flat_i[nn], trop_shape)
        pace_scan_orig = np.clip(
            pace_ji_all[0].astype(np.int64) * pace_s, 0, len(pace_scan_time) - 1
        )
        trop_scan_orig = np.clip(
            trop_ji_nn[0].astype(np.int64) * scan_s, 0, len(trop_scan_time) - 1
        )

        dt_min = (
            (trop_scan_time[trop_scan_orig] - pace_scan_time[pace_scan_orig])
            .astype("timedelta64[ms]")
            .astype(np.float64)
            / 60_000.0
        )
        keep = close_dist & (np.abs(dt_min) <= cfg.max_dt_min)
        n_keep = int(keep.sum())
        if n_keep == 0:
            continue
        n_keep_total += n_keep

        pace_ji = np.unravel_index(pace_flat_i[keep], pace_shape)
        trop_ji = np.unravel_index(trop_flat_i[nn[keep]], trop_shape)
        pf = pace_flat_i[keep]
        tf = trop_flat_i[nn[keep]]
        t_pace = _datetime64_to_unix(pace_scan_time[pace_scan_orig[keep]])
        t_trop = _datetime64_to_unix(trop_scan_time[trop_scan_orig[keep]])

        chunks.append(
            {
                "swath_id": np.full(n_keep, swath_id, dtype=np.int32),
                "pace_file_idx": np.full(n_keep, i_pace, dtype=np.int16),
                "pace_scan": (pace_ji[0].astype(np.int32) * pace_s),
                "pace_pix": (pace_ji[1].astype(np.int32) * pace_s),
                "trop_scan": (trop_ji[0].astype(np.int32) * scan_s),
                "trop_pix": (trop_ji[1].astype(np.int32) * pix_s),
                "pace_lat": pace_lat.ravel()[pf].astype(np.float32),
                "pace_lon": pace_lon.ravel()[pf].astype(np.float32),
                "trop_lat": trop_lat.ravel()[tf].astype(np.float32),
                "trop_lon": trop_lon.ravel()[tf].astype(np.float32),
                "dist_km": dist_km[keep].astype(np.float32),
                "dt_min": dt_min[keep].astype(np.float32),
                "time_pace": t_pace.astype(np.float64),
                "time_trop": t_trop.astype(np.float64),
                "pace_sza": pace_sza.ravel()[pf].astype(np.float32),
                "pace_vza": pace_vza.ravel()[pf].astype(np.float32),
                "trop_sza": trop_sza.ravel()[tf].astype(np.float32),
                "trop_vza": trop_vza.ravel()[tf].astype(np.float32),
            }
        )
        # string path columns stored separately
        chunks[-1]["_pace_path"] = np.full(n_keep, str(pace_path))
        chunks[-1]["_tropomi_path"] = np.full(n_keep, str(tropomi_path))
        chunks[-1]["_pace_name"] = np.full(n_keep, pace_path.name)

    if not chunks:
        log(f"  swath {swath_id}: no matches ({tropomi_path.name})")
        return None

    keys = [k for k in chunks[0] if not k.startswith("_")]
    out = {k: np.concatenate([c[k] for c in chunks]) for k in keys}
    out["pace_path"] = np.concatenate([c["_pace_path"] for c in chunks])
    out["tropomi_path"] = np.concatenate([c["_tropomi_path"] for c in chunks])
    out["pace_name"] = np.concatenate([c["_pace_name"] for c in chunks])
    log(f"  swath {swath_id}: kept {n_keep_total:,}  ({tropomi_path.name})")
    return out


def run_matches_stage(
    cfg: Config,
    dirs: dict[str, Path],
    swath_list: list[dict],
    pairs_fp: str,
    *,
    force: bool,
    log=print,
) -> tuple[dict, str]:
    fp_payload = {
        "pairs_fp": pairs_fp,
        "match": matches_fingerprint_payload(cfg),
        "tropomi_paths": sorted(s["tropomi_path"] for s in swath_list),
    }
    fp = cache.fingerprint(fp_payload)
    npz_path = dirs["matches"] / "matches.npz"
    meta_path = dirs["matches"] / "matches_meta.json"
    paths_json = dirs["matches"] / "match_paths.json"

    if (not force) and cache.is_fresh(
        meta_path, fp, required_files=[npz_path, paths_json], enabled=cfg.cache.enabled
    ):
        log(f"[matches] CACHE HIT  fp={fp[:12]}…")
        arrays = cache.load_npz(npz_path)
        path_info = cache.read_json(paths_json)
        arrays["pace_path"] = np.asarray(path_info["pace_path"])
        arrays["tropomi_path"] = np.asarray(path_info["tropomi_path"])
        arrays["pace_name"] = np.asarray(path_info["pace_name"])
        return arrays, fp

    log(f"[matches] CACHE MISS — matching {len(swath_list)} swaths")
    parts: list[dict] = []
    for sid, swath in enumerate(swath_list):
        part = match_one_swath(
            cfg,
            Path(swath["tropomi_path"]),
            [Path(p) for p in swath["pace_paths"]],
            swath_id=sid,
            log=log,
        )
        if part is not None:
            parts.append(part)

    if not parts:
        raise RuntimeError("No pixel matches found for any swath.")

    keys = [k for k in parts[0] if k not in ("pace_path", "tropomi_path", "pace_name")]
    arrays = {k: np.concatenate([p[k] for p in parts]) for k in keys}
    pace_path = np.concatenate([p["pace_path"] for p in parts])
    tropomi_path = np.concatenate([p["tropomi_path"] for p in parts])
    pace_name = np.concatenate([p["pace_name"] for p in parts])

    # NPZ cannot store unicode safely without pickle — keep paths in JSON
    cache.save_npz_json(
        npz_path,
        arrays,
        meta_path,
        {
            "fingerprint": fp,
            "n_match": int(arrays["dist_km"].size),
            "n_swaths_with_matches": len(parts),
            "payload": fp_payload,
        },
    )
    cache.write_json(
        paths_json,
        {
            "pace_path": pace_path.tolist(),
            "tropomi_path": tropomi_path.tolist(),
            "pace_name": pace_name.tolist(),
        },
    )
    arrays["pace_path"] = pace_path
    arrays["tropomi_path"] = tropomi_path
    arrays["pace_name"] = pace_name
    log(f"[matches] stored n={arrays['dist_km'].size:,}")
    return arrays, fp
