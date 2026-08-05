"""PACE OCI baseline SNR LUT loader."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import Config


def load_pace_band_snr_table(
    path: Path,
    fpa: str = "Red",
    λ_min: float = 600.0,
    λ_max: float = 900.0,
) -> dict:
    lines = Path(path).read_text().splitlines()
    header_end = next(i for i, ln in enumerate(lines) if "/end_header" in ln)
    rows = []
    for ln in lines[header_end + 1 :]:
        ln = ln.strip()
        if not ln or ln.startswith("!"):
            continue
        parts = ln.split()
        if len(parts) < 5 or parts[0] != fpa:
            continue
        wl = float(parts[1])
        if not (λ_min <= wl <= λ_max):
            continue
        rows.append((wl, int(float(parts[2])), float(parts[3]), float(parts[4])))
    if not rows:
        raise RuntimeError(f"No {fpa} SNR rows in ({λ_min}, {λ_max}) nm in {path}")
    rows.sort(key=lambda r: r[0])
    return {
        "wavelength": np.array([r[0] for r in rows], dtype=float),
        "band_index": np.array([r[1] for r in rows], dtype=int),
        "c1": np.array([r[2] for r in rows], dtype=float),
        "c2": np.array([r[3] for r in rows], dtype=float),
        "fpa": fpa,
    }


def interp_snr_coeffs(snr_tab: dict, λ_target: np.ndarray) -> dict:
    wl = snr_tab["wavelength"]
    c1_t = np.interp(λ_target, wl, snr_tab["c1"], left=snr_tab["c1"][0], right=snr_tab["c1"][-1])
    c2_t = np.interp(λ_target, wl, snr_tab["c2"], left=snr_tab["c2"][0], right=snr_tab["c2"][-1])
    return {"c1": c1_t, "c2": c2_t, "wavelength": np.asarray(λ_target, dtype=float)}


def load_pace_snr(cfg: Config, oci_bands: np.ndarray | None = None) -> dict:
    table = load_pace_band_snr_table(
        cfg.snr_path, fpa=cfg.snr_fpa, λ_min=cfg.snr_wl_min, λ_max=cfg.snr_wl_max
    )
    coeffs = None
    if oci_bands is not None:
        coeffs = interp_snr_coeffs(table, np.asarray(oci_bands, dtype=float))
    return {"path": str(cfg.snr_path), "table": table, "coeffs": coeffs}


def snr_sigma_from_L(L: np.ndarray, c1: np.ndarray, c2: np.ndarray) -> np.ndarray:
    L = np.maximum(np.asarray(L, dtype=float), 0.0)
    return np.sqrt(np.maximum(c1 + c2 * L, 0.0))
