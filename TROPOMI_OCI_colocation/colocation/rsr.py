"""PACE OCI RSR kernels on TROPOMI ground-pixel wavelength grids."""

from __future__ import annotations

from pathlib import Path

import netCDF4 as nc
import numpy as np
from scipy.interpolate import interp1d

from . import cache
from .config import Config, rsr_fingerprint_payload


def build_kernel_for_wl_out(
    wl_out: np.ndarray,
    rsr_b: np.ndarray,
    wavlen: np.ndarray,
    *,
    clip_negative: bool = True,
) -> np.ndarray:
    m_out = np.isfinite(wl_out)
    n_out = wl_out.size
    n_b = rsr_b.shape[0]
    K = np.zeros((n_b, n_out), dtype=np.float64)
    m_src = np.isfinite(wavlen) & np.all(np.isfinite(rsr_b), axis=0)
    x = wavlen[m_src]
    for i in range(n_b):
        y = rsr_b[i, m_src]
        if clip_negative:
            y = np.maximum(y, 0.0)
        f = interp1d(x, y, kind="linear", bounds_error=False, fill_value=0.0)
        row = np.zeros(n_out, dtype=np.float64)
        row[m_out] = f(wl_out[m_out])
        row = np.maximum(row, 0.0)
        s = row.sum()
        if s > 0:
            row /= s
        K[i] = row
    return K


def load_rsr_source(cfg: Config) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with nc.Dataset(cfg.rsr_path) as ds:
        rsr_wavlen = np.asarray(ds["wavelength"][:], dtype=float)
        rsr_bands_all = np.asarray(ds["bands"][:], dtype=float)
        rsr_all = np.asarray(ds["RSR"][:, :], dtype=float)
    return rsr_wavlen, rsr_bands_all, rsr_all


def build_rsr_bundle_for_tropomi(
    tropomi_path: Path,
    cfg: Config,
    rsr_wavlen: np.ndarray,
    rsr_bands_all: np.ndarray,
    rsr_all: np.ndarray,
) -> dict:
    with nc.Dataset(tropomi_path) as ds:
        trop_wl = np.asarray(
            ds["BAND5_RADIANCE"]["STANDARD_MODE"]["INSTRUMENT"]["nominal_wavelength"][
                0, :, :
            ],
            dtype=float,
        )
    n_gp, n_ch = trop_wl.shape
    # Use config band window only (same OCI bands for every TROPOMI file in a run)
    band_ok = (rsr_bands_all >= cfg.oci_band_min) & (rsr_bands_all <= cfg.oci_band_max)
    ib = np.flatnonzero(band_ok)
    if ib.size == 0:
        raise RuntimeError(
            f"No OCI RSR bands in [{cfg.oci_band_min}, {cfg.oci_band_max}] nm"
        )
    oci_bands = rsr_bands_all[ib]
    rsr_sel = rsr_all[ib, :]
    if cfg.clip_negative_rsr:
        rsr_sel = np.maximum(rsr_sel, 0.0)

    n_oci = oci_bands.size
    K = np.zeros((n_gp, n_oci, n_ch), dtype=np.float64)
    for gp in range(n_gp):
        K[gp] = build_kernel_for_wl_out(
            trop_wl[gp], rsr_sel, rsr_wavlen, clip_negative=cfg.clip_negative_rsr
        )

    return {
        "tropomi_path": str(tropomi_path),
        "oci_bands": oci_bands.astype(np.float64),
        "trop_wl": trop_wl.astype(np.float64),
        "K": K,
        "band_min": cfg.oci_band_min,
        "band_max": cfg.oci_band_max,
    }


def _stem_key(path: str | Path) -> str:
    return Path(path).stem


def run_rsr_stage(
    cfg: Config,
    dirs: dict[str, Path],
    tropomi_paths: list[str],
    matches_fp: str,
    *,
    force: bool,
    log=print,
) -> tuple[dict[str, dict], str]:
    paths_sorted = sorted(set(tropomi_paths))
    fp_payload = {
        "matches_fp": matches_fp,
        "rsr": rsr_fingerprint_payload(cfg),
        "tropomi_paths": paths_sorted,
    }
    fp = cache.fingerprint(fp_payload)
    index_path = dirs["rsr"] / "rsr_index.json"

    if (not force) and cfg.cache.enabled and index_path.is_file():
        try:
            index = cache.read_json(index_path)
        except Exception:
            index = None
        if index and index.get("fingerprint") == fp:
            missing = False
            bundles: dict[str, dict] = {}
            for tp, entry in index["files"].items():
                npz_p = dirs["rsr"] / entry["npz"]
                if not npz_p.is_file():
                    missing = True
                    break
                arr = cache.load_npz(npz_p)
                bundles[tp] = {
                    "tropomi_path": tp,
                    "oci_bands": arr["oci_bands"],
                    "trop_wl": arr["trop_wl"],
                    "K": arr["K"],
                    "band_min": cfg.oci_band_min,
                    "band_max": cfg.oci_band_max,
                }
            if not missing and set(bundles) == set(paths_sorted):
                log(f"[rsr] CACHE HIT  fp={fp[:12]}…  n_files={len(bundles)}")
                return bundles, fp

    log(f"[rsr] CACHE MISS — building kernels for {len(paths_sorted)} files")
    rsr_wavlen, rsr_bands_all, rsr_all = load_rsr_source(cfg)
    bundles = {}
    files_meta = {}
    for tp in paths_sorted:
        log(f"  RSR {Path(tp).name}")
        bundle = build_rsr_bundle_for_tropomi(
            Path(tp), cfg, rsr_wavlen, rsr_bands_all, rsr_all
        )
        npz_name = f"rsr_{_stem_key(tp)}.npz"
        npz_path = dirs["rsr"] / npz_name
        np.savez_compressed(
            npz_path,
            oci_bands=bundle["oci_bands"],
            trop_wl=bundle["trop_wl"],
            K=bundle["K"],
        )
        bundles[tp] = bundle
        files_meta[tp] = {"npz": npz_name, "n_gp": int(bundle["K"].shape[0]), "n_oci": int(bundle["K"].shape[1])}

    cache.write_json(
        index_path,
        {"fingerprint": fp, "files": files_meta, "payload": fp_payload},
    )
    log(f"[rsr] stored {len(bundles)} kernel sets")
    return bundles, fp
