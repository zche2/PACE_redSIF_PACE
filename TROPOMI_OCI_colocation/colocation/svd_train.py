"""Global loose-filtered SVD on TROPOMI BD5 spectra."""

from __future__ import annotations

from pathlib import Path

import netCDF4 as nc
import numpy as np

from . import cache
from .config import Config, svd_fingerprint_payload
from .radiance import tropomi_mol_to_mW


def _collect_training_indices_for_file(
    tropomi_path: Path,
    cfg: Config,
    log=print,
) -> tuple[np.ndarray, tuple[int, int]]:
    """Return flat indices on the *subsampled* geo grid that pass loose SVD filters."""
    scan_s = cfg.tropo_scan_stride
    pix_s = cfg.tropo_pix_stride
    s = cfg.svd

    with nc.Dataset(tropomi_path) as ds:
        g = ds["BAND5_RADIANCE"]["STANDARD_MODE"]["GEODATA"]
        trop_lat = np.asarray(g["latitude"][0, ::scan_s, ::pix_s], dtype=float)
        trop_lon = np.asarray(g["longitude"][0, ::scan_s, ::pix_s], dtype=float)
        trop_sza = np.asarray(g["solar_zenith_angle"][0, ::scan_s, ::pix_s], dtype=float)
        trop_vza = np.asarray(g["viewing_zenith_angle"][0, ::scan_s, ::pix_s], dtype=float)
        inst = ds["BAND5_RADIANCE"]["STANDARD_MODE"]["INSTRUMENT"]
        obs = ds["BAND5_RADIANCE"]["STANDARD_MODE"]["OBSERVATIONS"]
        trop_wl_all = np.asarray(inst["nominal_wavelength"][0, :, :], dtype=float)
        rad = obs["radiance"]

        trop_ok = np.isfinite(trop_lat) & np.isfinite(trop_lon)
        trop_shape = trop_lat.shape
        trop_valid_i = np.flatnonzero(trop_ok)
        sza_v = trop_sza.ravel()[trop_valid_i]
        vza_v = trop_vza.ravel()[trop_valid_i]
        ang_ok = (
            np.isfinite(sza_v)
            & np.isfinite(vza_v)
            & (sza_v < s.sza_max_deg)
            & (vza_v < s.vza_max_deg)
        )
        trop_ang_i = trop_valid_i[ang_ok]
        if trop_ang_i.size == 0:
            log(f"  SVD train skip (angles): {tropomi_path.name}")
            return np.array([], dtype=np.int64), trop_shape

        if s.sample_stride > 1:
            trop_ang_i = trop_ang_i[:: s.sample_stride]

        scan_sub_a, pix_sub_a = np.unravel_index(trop_ang_i, trop_shape)
        scan_orig_a = scan_sub_a * scan_s
        pix_orig_a = pix_sub_a * pix_s

        dark_mask = np.zeros(trop_ang_i.size, dtype=bool)
        wl_mid = trop_wl_all[trop_wl_all.shape[0] // 2]
        ib = np.flatnonzero((wl_mid >= s.lt_mask_wl_min) & np.isfinite(wl_mid))
        if ib.size == 0:
            raise RuntimeError(f"No TROPOMI channels ≥ {s.lt_mask_wl_min:g} nm")
        i0, i1 = int(ib[0]), int(ib[-1]) + 1

        for iscan in np.unique(scan_orig_a):
            sel = scan_orig_a == iscan
            gp_list = pix_orig_a[sel]
            mol_scan = np.asarray(rad[0, iscan, :, i0:i1], dtype=float)
            mol_scan = np.where(np.isfinite(mol_scan) & (mol_scan < 1e30), mol_scan, np.nan)
            for local_k, gp in zip(np.flatnonzero(sel), gp_list):
                wl = trop_wl_all[gp, i0:i1]
                lt = tropomi_mol_to_mW(mol_scan[gp], wl)
                m = (wl >= s.lt_mask_wl_min) & np.isfinite(lt)
                if np.any(m) and float(np.nanmax(lt[m])) < s.lt_max:
                    dark_mask[local_k] = True

        trop_dark_i = trop_ang_i[dark_mask]
        log(
            f"  SVD train {tropomi_path.name}: "
            f"ang={trop_ang_i.size:,}  dark={trop_dark_i.size:,}"
        )
        return trop_dark_i.astype(np.int64), trop_shape


def _spectra_on_grid(
    tropomi_path: Path,
    flat_idx: np.ndarray,
    trop_shape: tuple[int, int],
    cfg: Config,
) -> np.ndarray:
    scan_s = cfg.tropo_scan_stride
    pix_s = cfg.tropo_pix_stride
    s = cfg.svd
    wl_grid = np.linspace(s.wl_min, s.wl_max, s.n_wl)
    X = np.full((flat_idx.size, s.n_wl), np.nan, dtype=np.float64)
    if flat_idx.size == 0:
        return X

    scan_sub, pix_sub = np.unravel_index(flat_idx, trop_shape)
    scan_orig = scan_sub * scan_s
    pix_orig = pix_sub * pix_s

    with nc.Dataset(tropomi_path) as ds:
        inst = ds["BAND5_RADIANCE"]["STANDARD_MODE"]["INSTRUMENT"]
        rad = ds["BAND5_RADIANCE"]["STANDARD_MODE"]["OBSERVATIONS"]["radiance"]
        trop_wl_all = np.asarray(inst["nominal_wavelength"][0, :, :], dtype=float)

        for iscan in np.unique(scan_orig):
            sel = np.flatnonzero(scan_orig == iscan)
            mol_scan = np.asarray(rad[0, iscan, :, :], dtype=float)
            mol_scan = np.where(np.isfinite(mol_scan) & (mol_scan < 1e30), mol_scan, np.nan)
            for k in sel:
                gp = int(pix_orig[k])
                wl = trop_wl_all[gp]
                lt = tropomi_mol_to_mW(mol_scan[gp], wl)
                m = np.isfinite(wl) & np.isfinite(lt)
                if m.sum() < 10:
                    continue
                order = np.argsort(wl[m])
                X[k] = np.interp(
                    wl_grid, wl[m][order], lt[m][order], left=np.nan, right=np.nan
                )
    return X


def build_global_svd(
    cfg: Config,
    tropomi_paths: list[str],
    log=print,
) -> dict:
    s = cfg.svd
    wl_grid = np.linspace(s.wl_min, s.wl_max, s.n_wl)
    blocks: list[np.ndarray] = []
    n_per_file: dict[str, int] = {}

    for path_s in tropomi_paths:
        path = Path(path_s)
        idx, trop_shape = _collect_training_indices_for_file(path, cfg, log=log)
        if idx.size == 0:
            n_per_file[path_s] = 0
            continue
        X = _spectra_on_grid(path, idx, trop_shape, cfg)
        row_ok = np.isfinite(X).all(axis=1)
        X = X[row_ok]
        n_per_file[path_s] = int(X.shape[0])
        if X.shape[0]:
            blocks.append(X)

    if not blocks:
        raise RuntimeError("No TROPOMI spectra available for SVD training.")

    X_all = np.vstack(blocks)
    if X_all.shape[0] > s.max_n:
        take = np.linspace(0, X_all.shape[0] - 1, s.max_n, dtype=int)
        X_all = X_all[take]
        log(f"SVD subsampled to max_n={s.max_n:,}")

    log(f"SVD matrix: {X_all.shape[0]:,} × {X_all.shape[1]}")
    if X_all.shape[0] < 10:
        raise RuntimeError("Too few complete spectra for SVD.")

    mean_spec = X_all.mean(axis=0) if s.mean_center else np.zeros(s.n_wl)
    Xc = X_all - mean_spec
    # Do not keep U (large); only spectral basis
    _, singular, Vt = np.linalg.svd(Xc, full_matrices=False)
    var_frac = (singular ** 2) / (singular ** 2).sum()

    return {
        "wl": wl_grid.astype(np.float64),
        "mean": mean_spec.astype(np.float64),
        "s": singular.astype(np.float64),
        "Vt": Vt.astype(np.float64),
        "var_frac": var_frac.astype(np.float64),
        "n_train": int(X_all.shape[0]),
        "mean_center": bool(s.mean_center),
        "wl_min": float(s.wl_min),
        "wl_max": float(s.wl_max),
        "n_per_file": n_per_file,
    }


def run_svd_stage(
    cfg: Config,
    dirs: dict[str, Path],
    tropomi_paths: list[str],
    matches_fp: str,
    *,
    force: bool,
    log=print,
) -> tuple[dict, str]:
    paths_sorted = sorted(set(tropomi_paths))
    fp_payload = {
        "matches_fp": matches_fp,
        "svd": svd_fingerprint_payload(cfg),
        "tropomi_paths": paths_sorted,
    }
    fp = cache.fingerprint(fp_payload)
    npz_path = dirs["svd"] / "trop_svd.npz"
    meta_path = dirs["svd"] / "svd_meta.json"

    if (not force) and cache.is_fresh(
        meta_path, fp, required_files=[npz_path], enabled=cfg.cache.enabled
    ):
        log(f"[svd] CACHE HIT  fp={fp[:12]}…")
        arrays = cache.load_npz(npz_path)
        meta = cache.read_json(meta_path)
        trop_svd = {
            "wl": arrays["wl"],
            "mean": arrays["mean"],
            "s": arrays["s"],
            "Vt": arrays["Vt"],
            "var_frac": arrays["var_frac"],
            "n_train": int(meta["n_train"]),
            "mean_center": bool(meta["mean_center"]),
            "wl_min": float(meta["wl_min"]),
            "wl_max": float(meta["wl_max"]),
            "n_per_file": meta.get("n_per_file", {}),
        }
        return trop_svd, fp

    log(f"[svd] CACHE MISS — training on {len(paths_sorted)} TROPOMI files (loose filters)")
    trop_svd = build_global_svd(cfg, paths_sorted, log=log)
    arrays = {
        "wl": trop_svd["wl"],
        "mean": trop_svd["mean"],
        "s": trop_svd["s"],
        "Vt": trop_svd["Vt"],
        "var_frac": trop_svd["var_frac"],
    }
    cache.save_npz_json(
        npz_path,
        arrays,
        meta_path,
        {
            "fingerprint": fp,
            "n_train": trop_svd["n_train"],
            "mean_center": trop_svd["mean_center"],
            "wl_min": trop_svd["wl_min"],
            "wl_max": trop_svd["wl_max"],
            "n_per_file": trop_svd["n_per_file"],
            "var_pc1_10_pct": float(100 * trop_svd["var_frac"][:10].sum()),
            "payload": fp_payload,
        },
    )
    log(
        f"[svd] done n_train={trop_svd['n_train']:,}  "
        f"var PC1–10={100 * trop_svd['var_frac'][:10].sum():.2f}%"
    )
    return trop_svd, fp
