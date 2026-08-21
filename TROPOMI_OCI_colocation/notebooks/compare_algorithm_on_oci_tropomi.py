#!/usr/bin/env python3
"""OCI vs TROPOMI row-aware red SIF intercomparison (one or both July 2025 cycles).

Matches co-located L1 pixels from ``co-location_results_202507.json``, joins
independent retrievals (per-orbit cycle1/cycle2 dirs), converts TROPOMI
``sif_red_pc1`` + ``sif_red_pc2`` to SIF at 678 nm, and writes figures plus a
NetCDF matchup table for later scatter replots.

Run from this directory:
  python compare_algorithm_on_oci_tropomi.py              # both cycles
  python compare_algorithm_on_oci_tropomi.py --cycle cycle1
  python compare_algorithm_on_oci_tropomi.py --plots-only
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4 as nc  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import constants  # noqa: E402

COLOC_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(COLOC_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from colocation.config import Config  # noqa: E402
from colocation.match import match_one_swath  # noqa: E402
from tropomi_sif_at_wavelength import load_sif_basis_weights, sif_from_coefficients  # noqa: E402

DEFAULT_THRESHOLDS = {
    "max_dist_km": 2.5,
    "max_dt_min": 30.0,
    "lt_max_oci": 30.0,
    "lt_mask_wl_min_oci": 600.0,
    "pace_stride_pix": 50,
    "tropo_scan_stride": 2,
    "tropo_pix_stride": 2,
    "max_abs_dvza_deg": 20.0,
    "max_sza_deg": 70.0,
    "require_both_retrieved": True,
    "trop_qc": "paper",
}

RSR_PATH = Path("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI/PACE_OCI_RSRs.nc")
SNR_PATH = Path(
    "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/PACE_OCI/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt"
)
TROP_RET_ROOT = Path("/net/squid/data1/projects/TROPOMI/red_sif_july2025_row_aware_v2")
TROP_BASIS_ROOT = Path(
    "/net/squid/data1/projects/TROPOMI/red_sif_july2025_map_v1/diagnostics/row_aware_production"
)
PACE_RET_DIR = Path("/home/zhe2/data/PACE/svd_retrieval_output")
COLOCATION_JSON = COLOC_ROOT / "notebooks" / "co-location_results_202507.json"
OUTPUT_ROOT = Path(__file__).resolve().parent / "outputs" / "compare_oci_tropomi_202507"

N_EXAMPLE_SPECTRA = 6
WL_PLOT_MIN, WL_PLOT_MAX = 665.0, 715.0
REPORTING_WL_NM = 678.0
_PACE_STAMP_RE = re.compile(r"PACE_OCI\.(\d{8}T\d{6})")

# July 2025 TROPOMI red-SIF production cycle orbit spans (inclusive).
CYCLE_ORBIT_RANGES: dict[str, tuple[int, int]] = {
    "cycle1": (39989, 40201),
    "cycle2": (40202, 40414),
}


def parse_cycles(cycle_arg: str) -> list[str]:
    """Expand --cycle into one or more cycle names."""
    token = str(cycle_arg).strip().lower()
    if token in {"all", "both"}:
        return ["cycle1", "cycle2"]
    if token in CYCLE_ORBIT_RANGES:
        return [token]
    raise ValueError(
        f"Unknown cycle {cycle_arg!r}; use cycle1, cycle2, or all"
    )


def orbit_to_cycle(orbit: int) -> str | None:
    orb = int(orbit)
    for name, (lo, hi) in CYCLE_ORBIT_RANGES.items():
        if lo <= orb <= hi:
            return name
    return None


def default_basis_path(cycle: str) -> Path:
    return TROP_BASIS_ROOT / cycle / "band5_row_aware_procrustes4_basis.h5"


def trop_ret_path_for_orbit(orbit: int) -> Path | None:
    cycle = orbit_to_cycle(orbit)
    if cycle is None:
        return None
    path = TROP_RET_ROOT / cycle / "retrieval" / f"TROPOMI_RED_SIF_orbit_{int(orbit)}.h5"
    return path if path.is_file() else None


def stamp_from_pace_path(path: str | Path) -> str:
    match = _PACE_STAMP_RE.search(Path(path).name)
    if not match:
        raise ValueError(f"Cannot parse PACE stamp from {path}")
    return match.group(1)


def resolve_pace_ret(stamp: str) -> Path | None:
    gpu = PACE_RET_DIR / f"interim_{stamp}_svd_retrieval_full_gpu.nc"
    if gpu.is_file():
        return gpu
    parallel = PACE_RET_DIR / f"interim_{stamp}_svd_retrieval_full_parallel.nc"
    if parallel.is_file():
        return parallel
    candidates = sorted(PACE_RET_DIR.glob(f"interim_{stamp}_svd_retrieval_*.nc"))
    return candidates[0] if candidates else None


def build_jobs(
    cycles: list[str],
    orbit: int | None,
    max_swaths: int | None,
    colocation_json: Path,
) -> tuple[list[dict], Path]:
    with open(colocation_json) as handle:
        coloc = json.load(handle)

    meta = coloc["meta"]
    pace_dir = Path(meta["pace_dir"])
    allowed = set(cycles)

    # Orbit → retrieval file for selected cycles only.
    orbit_ret: dict[int, tuple[str, Path]] = {}
    for cycle in cycles:
        ret_dir = TROP_RET_ROOT / cycle / "retrieval"
        if not ret_dir.is_dir():
            raise FileNotFoundError(ret_dir)
        for path in ret_dir.glob("TROPOMI_RED_SIF_orbit_*.h5"):
            orb = int(path.stem.split("_")[-1])
            if orbit_to_cycle(orb) != cycle:
                continue
            orbit_ret[orb] = (cycle, path)
    if not orbit_ret:
        raise FileNotFoundError(
            f"no TROPOMI_RED_SIF_orbit_*.h5 for cycles={cycles} under {TROP_RET_ROOT}"
        )

    swaths = coloc["swaths"]
    if orbit is not None:
        swaths = [swath for swath in swaths if int(swath["orbit"]) == int(orbit)]
        if not swaths:
            available = sorted({int(swath["orbit"]) for swath in coloc["swaths"]})
            raise ValueError(
                f"Orbit {orbit} not in {colocation_json.name}; "
                f"available e.g. {available[:10]}..."
            )

    jobs: list[dict] = []
    skipped = 0
    for swath in swaths:
        orb = int(swath["orbit"])
        mapped = orbit_ret.get(orb)
        if mapped is None:
            skipped += 1
            continue
        cycle, trop_ret = mapped
        if cycle not in allowed:
            skipped += 1
            continue
        trop_l1 = Path(swath["tropomi_path"])
        pace_l1 = [pace_dir / name for name in swath["pace_granules"]]
        missing_pace = [path for path in pace_l1 if not path.is_file()]
        if not trop_l1.is_file() or missing_pace or not trop_ret.is_file():
            skipped += 1
            continue
        jobs.append(
            {
                "orbit": orb,
                "cycle": cycle,
                "swath": swath,
                "tropomi_path": trop_l1,
                "pace_paths": pace_l1,
                "trop_ret_path": trop_ret,
            }
        )

    if max_swaths is not None:
        # Per selected cycle so --max-swaths works for --cycle all smoke tests.
        limit = int(max_swaths)
        kept: list[dict] = []
        counts: dict[str, int] = defaultdict(int)
        for job in jobs:
            if counts[job["cycle"]] >= limit:
                continue
            kept.append(job)
            counts[job["cycle"]] += 1
        jobs = kept
    if not jobs:
        raise RuntimeError(
            f"No usable swaths for cycles={cycles} orbit={orbit!r} "
            f"(skipped {skipped}: missing L1 or retrieval)"
        )

    by_cycle = defaultdict(int)
    for job in jobs:
        by_cycle[job["cycle"]] += 1
    print(
        f"cycles={cycles}  orbit={orbit!r}  swaths={len(jobs)}  skipped={skipped}\n"
        f"  per-cycle: {dict(by_cycle)}\n"
        f"  TROPOMI ret root: {TROP_RET_ROOT}\n"
        f"  colocation: {colocation_json}"
    )
    return jobs, pace_dir


def match_swaths(jobs: list[dict], pace_dir: Path, thresholds: dict) -> dict:
    cfg = Config(
        pace_dir=pace_dir,
        tropomi_dir=Path("."),
        rsr_path=RSR_PATH,
        snr_path=SNR_PATH,
        output_dir=Path("/tmp/oci_tropomi_compare"),
        max_dist_km=float(thresholds["max_dist_km"]),
        max_dt_min=float(thresholds["max_dt_min"]),
        lt_max_oci=float(thresholds["lt_max_oci"]),
        lt_mask_wl_min_oci=float(thresholds["lt_mask_wl_min_oci"]),
        pace_stride_pix=int(thresholds["pace_stride_pix"]),
        tropo_scan_stride=int(thresholds["tropo_scan_stride"]),
        tropo_pix_stride=int(thresholds["tropo_pix_stride"]),
    )

    parts: list[dict] = []
    for sid, job in enumerate(jobs):
        print(
            f"  matching orbit {job['orbit']} [{job['cycle']}] "
            f"({sid + 1}/{len(jobs)})...",
            flush=True,
        )
        part = match_one_swath(
            cfg, job["tropomi_path"], job["pace_paths"], swath_id=sid
        )
        if part is None:
            print(f"    no matches")
            continue
        n = int(part["dist_km"].size)
        part["orbit"] = np.full(n, job["orbit"], dtype=np.int32)
        part["cycle"] = np.full(n, job["cycle"], dtype=object)
        part["trop_ret_path"] = np.full(n, str(job["trop_ret_path"]), dtype=object)
        parts.append(part)

    if not parts:
        raise RuntimeError("No raw matches with current thresholds")

    keys = list(parts[0].keys())
    merged = {key: np.concatenate([part[key] for part in parts]) for key in keys}
    print(f"Raw matches: {merged['dist_km'].size:,} from {len(parts)} swaths")
    return merged


def apply_geo_filters(match: dict, thresholds: dict) -> tuple[np.ndarray, np.ndarray]:
    n_raw = int(match["dist_km"].size)
    dvza = np.abs(
        np.asarray(match["pace_vza"], dtype=float) - np.asarray(match["trop_vza"], dtype=float)
    )
    sel = np.ones(n_raw, dtype=bool)

    max_dvza = thresholds.get("max_abs_dvza_deg")
    if max_dvza is not None:
        sel &= dvza <= float(max_dvza)
        print(f"  after max_abs_dvza_deg={max_dvza}: {int(sel.sum()):,}")

    max_sza = thresholds.get("max_sza_deg")
    if max_sza is not None:
        sel &= (np.asarray(match["pace_sza"], dtype=float) <= float(max_sza)) & (
            np.asarray(match["trop_sza"], dtype=float) <= float(max_sza)
        )
        print(f"  after max_sza_deg={max_sza}: {int(sel.sum()):,}")

    print(f"Pairs after geo post-filters (before retrieval gate): {int(sel.sum()):,}")
    return sel, dvza


def join_retrievals(
    match: dict,
    bases: dict[str, object],
    thresholds: dict,
) -> dict[str, np.ndarray]:
    n_raw = int(match["dist_km"].size)
    trop_scan = np.asarray(match["trop_scan"], dtype=np.int32)
    trop_pix = np.asarray(match["trop_pix"], dtype=np.int32)
    pace_scan = np.asarray(match["pace_scan"], dtype=np.int32)
    pace_pix = np.asarray(match["pace_pix"], dtype=np.int32)
    pace_path_arr = np.asarray(match["pace_path"])
    orbits = np.asarray(match["orbit"], dtype=np.int32)
    cycles = np.asarray(match["cycle"], dtype=object)
    trop_ret_paths = np.asarray(match["trop_ret_path"], dtype=object)

    sif_trop_pc1 = np.full(n_raw, np.nan, dtype=np.float64)
    sif_trop_pc2 = np.full(n_raw, np.nan, dtype=np.float64)
    sif_trop_678 = np.full(n_raw, np.nan, dtype=np.float64)
    chi2_trop = np.full(n_raw, np.nan, dtype=np.float64)
    retrieved_trop = np.zeros(n_raw, dtype=bool)

    trop_qc = str(thresholds.get("trop_qc", "paper")).lower()
    n_trop_miss = 0

    by_orbit: dict[int, list[int]] = defaultdict(list)
    for i, orb in enumerate(orbits):
        by_orbit[int(orb)].append(i)

    for orb, idxs in by_orbit.items():
        trop_ret_path = Path(str(trop_ret_paths[idxs[0]]))
        cycle = str(cycles[idxs[0]])
        basis = bases[cycle]
        with h5py.File(trop_ret_path, "r") as file:
            names = [
                item.decode() if isinstance(item, (bytes, np.bytes_)) else str(item)
                for item in file["state_name"][:]
            ]
            i_pc1 = names.index("sif_red_pc1")
            i_pc2 = names.index("sif_red_pc2")
            coeff = np.asarray(file["coefficient"][:], dtype=np.float64)
            pc1_flat = coeff[:, i_pc1]
            pc2_flat = coeff[:, i_pc2]
            chi2_flat = np.asarray(file["reduced_chi_square"][:], dtype=np.float64)
            det_row = np.asarray(file["detector_row"][:], dtype=np.int32)
            det_scan = np.asarray(file["scanline"][:], dtype=np.int32)
            if trop_qc == "paper" and "accepted_paper" in file:
                qc_flat = np.asarray(file["accepted_paper"][:], dtype=np.uint8) > 0
            elif trop_qc == "exploratory" and "accepted_exploratory" in file:
                qc_flat = np.asarray(file["accepted_exploratory"][:], dtype=np.uint8) > 0
            else:
                qc_flat = np.isfinite(pc1_flat)
            lut = {
                int(row) * 100_000 + int(scan): j
                for j, (row, scan) in enumerate(zip(det_row, det_scan))
            }

        sif_orbit = sif_from_coefficients(det_row, pc1_flat, pc2_flat, basis)

        for i in idxs:
            j = lut.get(int(trop_pix[i]) * 100_000 + int(trop_scan[i]))
            if j is None:
                n_trop_miss += 1
                continue
            sif_trop_pc1[i] = pc1_flat[j]
            sif_trop_pc2[i] = pc2_flat[j]
            sif_trop_678[i] = sif_orbit[j]
            chi2_trop[i] = chi2_flat[j]
            retrieved_trop[i] = bool(qc_flat[j]) and np.isfinite(sif_orbit[j])

    print(
        f"TROPOMI_RED_SIF: trop_qc={trop_qc!r}  orbits={len(by_orbit)}  "
        f"index misses={n_trop_miss:,}"
    )

    sif_pace = np.full(n_raw, np.nan, dtype=np.float64)
    chi2_pace = np.full(n_raw, np.nan, dtype=np.float64)
    converged_pace = np.zeros(n_raw, dtype=bool)
    pace_ret_path = np.array([""] * n_raw, dtype=object)

    by_stamp: dict[str, list[int]] = defaultdict(list)
    for i in range(n_raw):
        by_stamp[stamp_from_pace_path(pace_path_arr[i])].append(i)

    n_pace_file_miss = 0
    for stamp, idxs in by_stamp.items():
        ret_path = resolve_pace_ret(stamp)
        idxs_arr = np.asarray(idxs, dtype=np.int64)
        if ret_path is None:
            n_pace_file_miss += int(idxs_arr.size)
            continue
        with nc.Dataset(ret_path) as dataset:
            sif = np.asarray(dataset["sif_radiance_678nm"][:], dtype=np.float64)
            chi2 = np.asarray(dataset["reduced_chi2"][:], dtype=np.float64)
            if "converged" in dataset.variables:
                conv = np.asarray(dataset["converged"][:], dtype=np.uint8) > 0
            else:
                conv = np.isfinite(sif)
            scans = pace_scan[idxs_arr]
            pixels = pace_pix[idxs_arr]
            ok = (
                (scans >= 0)
                & (pixels >= 0)
                & (scans < sif.shape[0])
                & (pixels < sif.shape[1])
            )
            jj = idxs_arr[ok]
            sif_pace[jj] = sif[scans[ok], pixels[ok]]
            chi2_pace[jj] = chi2[scans[ok], pixels[ok]]
            converged_pace[jj] = conv[scans[ok], pixels[ok]]
            pace_ret_path[jj] = str(ret_path)

    both_ok = (
        retrieved_trop
        & converged_pace
        & np.isfinite(sif_trop_678)
        & np.isfinite(sif_pace)
    )
    print(
        f"Retrieval hits: trop={int(retrieved_trop.sum()):,}  "
        f"pace={int(converged_pace.sum()):,}  both={int(both_ok.sum()):,}"
        + (f"  pace file miss={n_pace_file_miss:,}" if n_pace_file_miss else "")
    )

    return {
        "sif_trop_pc1": sif_trop_pc1,
        "sif_trop_pc2": sif_trop_pc2,
        "sif_trop_678": sif_trop_678,
        "chi2_trop": chi2_trop,
        "retrieved_trop": retrieved_trop,
        "sif_pace": sif_pace,
        "chi2_pace": chi2_pace,
        "converged_pace": converged_pace,
        "both_ok": both_ok,
        "pace_ret_path": pace_ret_path,
        "trop_scan": trop_scan,
        "trop_pix": trop_pix,
        "pace_scan": pace_scan,
        "pace_pix": pace_pix,
        "pace_path_arr": pace_path_arr,
        "orbits": orbits,
        "cycles": cycles,
    }


def tropomi_mol_to_mw(radiance_mol: np.ndarray, wavelength_nm: np.ndarray) -> np.ndarray:
    h, c_light, na = constants.h, constants.c, constants.N_A
    return radiance_mol * (h * c_light * na / wavelength_nm) * 1e9 * 1e3


def load_tropomi_spectrum(path: Path, scan: int, pix: int) -> tuple[np.ndarray, np.ndarray]:
    with nc.Dataset(path) as dataset:
        obs = dataset["BAND5_RADIANCE"]["STANDARD_MODE"]["OBSERVATIONS"]
        inst = dataset["BAND5_RADIANCE"]["STANDARD_MODE"]["INSTRUMENT"]
        wl = np.asarray(inst["nominal_wavelength"][0, pix, :], dtype=float)
        rad = np.asarray(obs["radiance"][0, scan, pix, :], dtype=float)
    return wl, tropomi_mol_to_mw(rad, wl)


def load_oci_spectrum(path: Path, scan: int, pix: int) -> tuple[np.ndarray, np.ndarray, float]:
    with nc.Dataset(path) as dataset:
        wl = np.asarray(dataset["sensor_band_parameters"]["red_wavelength"][:], dtype=float)
        f0 = np.asarray(dataset["sensor_band_parameters"]["red_solar_irradiance"][:], dtype=float)
        esd = float(getattr(dataset, "earth_sun_distance_correction", 1.0))
        sza = float(dataset["geolocation_data"]["solar_zenith"][scan, pix])
        rhot = np.asarray(dataset["observation_data"]["rhot_red"][:, scan, pix], dtype=float)
    rhot = np.where(rhot > 0, rhot, np.nan)
    mu0 = np.cos(np.deg2rad(sza))
    mu0 = mu0 if mu0 > 1e-6 else np.nan
    lt = rhot * f0 * mu0 / (np.pi * esd)
    return wl, lt, sza


def save_matchup_netcdf(
    path: Path,
    *,
    cycle_label: str,
    orbit: int | None,
    thresholds: dict,
    basis_paths: dict[str, Path],
    colocation_json: Path,
    match: dict,
    joined: dict,
    sel: np.ndarray,
    dvza: np.ndarray,
) -> np.ndarray:
    idx = np.flatnonzero(sel)
    path.parent.mkdir(parents=True, exist_ok=True)

    if path.is_file():
        path.unlink()

    with nc.Dataset(path, "w") as dataset:
        dataset.createDimension("pair", idx.size)
        dataset.title = "OCI–TROPOMI co-located SIF matchup"
        dataset.history = f"Created {datetime.now(timezone.utc).isoformat()}"
        dataset.cycle = cycle_label
        dataset.orbit = "all" if orbit is None else str(int(orbit))
        dataset.reporting_wavelength_nm = float(REPORTING_WL_NM)
        dataset.tropomi_sif_definition = (
            "sif_red_pc1*basis_pc1(wavelength)+sif_red_pc2*basis_pc2(wavelength)"
        )
        dataset.source_basis = json.dumps({k: str(v) for k, v in basis_paths.items()})
        dataset.colocation_json = str(colocation_json)
        dataset.thresholds_json = json.dumps(thresholds)
        dataset.tropomi_root = str(TROP_RET_ROOT)
        dataset.pace_root = str(PACE_RET_DIR)

        def var(name: str, values: np.ndarray, dtype="f8", **attrs):
            variable = dataset.createVariable(name, dtype, ("pair",))
            variable[:] = values
            for key, value in attrs.items():
                setattr(variable, key, value)
            return variable

        var("orbit", joined["orbits"][idx].astype(np.int32), dtype="i4")
        cycle_var = dataset.createVariable("cycle", str, ("pair",))
        cycle_var[:] = np.asarray(joined["cycles"], dtype=object)[idx]
        var("trop_scan", joined["trop_scan"][idx].astype(np.int32), dtype="i4")
        var("trop_pix", joined["trop_pix"][idx].astype(np.int32), dtype="i4")
        var("pace_scan", joined["pace_scan"][idx].astype(np.int32), dtype="i4")
        var("pace_pix", joined["pace_pix"][idx].astype(np.int32), dtype="i4")
        var("dist_km", match["dist_km"][idx].astype(np.float64), long_name="great-circle distance")
        var("dt_min", match["dt_min"][idx].astype(np.float64), long_name="TROPOMI minus PACE [min]")
        var("dvza_deg", dvza[idx].astype(np.float64), long_name="abs(pace_vza - trop_vza)")
        var("pace_lat", match["pace_lat"][idx].astype(np.float64), units="degrees_north")
        var("pace_lon", match["pace_lon"][idx].astype(np.float64), units="degrees_east")
        var("pace_sza", match["pace_sza"][idx].astype(np.float64), units="degree")
        var("trop_sza", match["trop_sza"][idx].astype(np.float64), units="degree")
        var("sif_trop_pc1", joined["sif_trop_pc1"][idx], long_name="TROPOMI sif_red_pc1 coefficient")
        var("sif_trop_pc2", joined["sif_trop_pc2"][idx], long_name="TROPOMI sif_red_pc2 coefficient")
        var(
            "sif_trop_678nm",
            joined["sif_trop_678"][idx],
            long_name="TROPOMI SIF radiance at 678 nm",
            units="W m-2 sr-1 um-1",
        )
        var(
            "sif_pace_678nm",
            joined["sif_pace"][idx],
            long_name="PACE SIF radiance at 678 nm",
            units="W m-2 sr-1 um-1",
        )
        var("chi2_trop", joined["chi2_trop"][idx], long_name="TROPOMI reduced chi-square")
        var("chi2_pace", joined["chi2_pace"][idx], long_name="PACE reduced chi-square")

        pace_path = dataset.createVariable("pace_l1_path", str, ("pair",))
        pace_path[:] = np.asarray(match["pace_path"], dtype=object)[idx]
        trop_path = dataset.createVariable("trop_l1_path", str, ("pair",))
        trop_path[:] = np.asarray(match["tropomi_path"], dtype=object)[idx]
        pace_ret = dataset.createVariable("pace_ret_path", str, ("pair",))
        pace_ret[:] = joined["pace_ret_path"][idx]

    print(f"Saved matchup table: {path}  (n={idx.size:,})")
    return idx


def load_matchup_netcdf(path: Path) -> dict:
    with nc.Dataset(path) as dataset:
        data = {name: np.asarray(dataset[name][:]) for name in dataset.variables}
        meta = {
            "cycle": getattr(dataset, "cycle", ""),
            "orbit": getattr(dataset, "orbit", ""),
            "reporting_wavelength_nm": float(getattr(dataset, "reporting_wavelength_nm", REPORTING_WL_NM)),
        }
    data.update(meta)
    return data


def scatter_stats(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    if x.size == 0:
        return {"n": 0, "bias": np.nan, "rmse": np.nan, "r": np.nan}
    bias = float(np.mean(y - x))
    rmse = float(np.sqrt(np.mean((y - x) ** 2)))
    r = float(np.corrcoef(x, y)[0, 1]) if x.size > 1 else np.nan
    return {"n": int(x.size), "bias": bias, "rmse": rmse, "r": r}


def plot_match_qa(
    fig_dir: Path,
    *,
    cycle: str,
    orbit: int | None,
    match: dict,
    idx: np.ndarray,
    dvza: np.ndarray,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    axes[0].hist(match["dist_km"][idx], bins=40, color="#4C72B0", edgecolor="white")
    axes[0].set_xlabel("dist_km")
    axes[0].set_ylabel("count")
    axes[0].set_title("Distance")

    axes[1].hist(match["dt_min"][idx], bins=40, color="#55A868", edgecolor="white")
    axes[1].set_xlabel("dt_min (TROPOMI − PACE)")
    axes[1].set_title("Δt")

    axes[2].hist(dvza[idx], bins=40, color="#C44E52", edgecolor="white")
    axes[2].set_xlabel(r"$|\Delta$VZA$|$ (deg)")
    axes[2].set_title("|ΔVZA|")

    orbit_label = f"orbit {orbit}" if orbit is not None else "all co-located orbits"
    fig.suptitle(f"cycle={cycle} ({orbit_label}): match QA (n={idx.size:,})", y=1.02)
    fig.tight_layout()
    fig.savefig(fig_dir / "match_qa_hist.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(
        match["pace_lon"][idx],
        match["pace_lat"][idx],
        c=match["dist_km"][idx],
        s=8,
        cmap="viridis",
        alpha=0.8,
    )
    fig.colorbar(sc, ax=ax, label="dist_km")
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_title(f"cycle={cycle}: co-located pixels")
    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    fig.savefig(fig_dir / "match_map.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_example_spectra(
    fig_dir: Path,
    *,
    cycle: str,
    orbit: int | None,
    match: dict,
    joined: dict,
    idx: np.ndarray,
    rng: np.random.Generator,
) -> None:
    trop_path_arr = np.asarray(match["tropomi_path"])
    pace_path_arr = joined["pace_path_arr"]
    n_ex = min(N_EXAMPLE_SPECTRA, idx.size)
    if n_ex == 0:
        return

    example_i = np.sort(rng.choice(idx, size=n_ex, replace=False))
    examples = []
    for i in example_i:
        wl_t, lt_t = load_tropomi_spectrum(
            Path(str(trop_path_arr[i])), int(joined["trop_scan"][i]), int(joined["trop_pix"][i])
        )
        wl_o, lt_o, _ = load_oci_spectrum(
            Path(str(pace_path_arr[i])), int(joined["pace_scan"][i]), int(joined["pace_pix"][i])
        )
        examples.append(
            {
                "dist_km": float(match["dist_km"][i]),
                "dt_min": float(match["dt_min"][i]),
                "sif_trop": float(joined["sif_trop_678"][i]),
                "sif_pace": float(joined["sif_pace"][i]),
                "chi2_trop": float(joined["chi2_trop"][i]),
                "chi2_pace": float(joined["chi2_pace"][i]),
                "wl_trop": wl_t,
                "lt_trop": lt_t,
                "wl_oci": wl_o,
                "lt_oci": lt_o,
            }
        )

    ncols = 2
    nrows = int(np.ceil(len(examples) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 3.2 * nrows), squeeze=False)
    for ax, ex in zip(axes.ravel(), examples):
        mt = (
            (ex["wl_trop"] >= WL_PLOT_MIN)
            & (ex["wl_trop"] <= WL_PLOT_MAX)
            & np.isfinite(ex["lt_trop"])
        )
        mo = (
            (ex["wl_oci"] >= WL_PLOT_MIN)
            & (ex["wl_oci"] <= WL_PLOT_MAX)
            & np.isfinite(ex["lt_oci"])
        )
        ax.plot(ex["wl_trop"][mt], ex["lt_trop"][mt], color="#1f77b4", lw=1.0, label="TROPOMI BD5")
        ax.plot(ex["wl_oci"][mo], ex["lt_oci"][mo], color="#EC4B26", lw=1.2, label="OCI red")
        ax.axvline(REPORTING_WL_NM, color="gray", ls=":", lw=0.8)
        ax.set_xlabel("λ (nm)")
        ax.set_ylabel(r"$L_t$ (mW m$^{-2}$ nm$^{-1}$ sr$^{-1}$)")
        ax.set_title(
            f"d={ex['dist_km']:.2f} km  Δt={ex['dt_min']:.1f} min\n"
            f"SIF$_T$={ex['sif_trop']:.3g}  SIF$_P$={ex['sif_pace']:.3g}  "
            f"χ²$_T$={ex['chi2_trop']:.2f}  χ²$_P$={ex['chi2_pace']:.2f}",
            fontsize=9,
        )
        ax.legend(fontsize=8, loc="upper right")

    for ax in axes.ravel()[len(examples) :]:
        ax.set_visible(False)

    orbit_label = f"orbit {orbit}" if orbit is not None else f"cycle {cycle}"
    fig.suptitle(
        f"{orbit_label}: co-located L1 spectra ({WL_PLOT_MIN:.0f}–{WL_PLOT_MAX:.0f} nm)",
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(fig_dir / "example_spectra.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_sif_comparison(
    fig_dir: Path,
    *,
    cycle: str,
    x: np.ndarray,
    y: np.ndarray,
    chi2_trop: np.ndarray,
    chi2_pace: np.ndarray,
    stats: dict[str, float],
) -> None:
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    chi2_trop = chi2_trop[finite]
    chi2_pace = chi2_pace[finite]

    lo = float(np.nanpercentile(np.concatenate([x, y]), 1))
    hi = float(np.nanpercentile(np.concatenate([x, y]), 99))
    pad = 0.05 * (hi - lo + 1e-9)
    lo, hi = lo - pad, hi + pad

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    hb = ax.hexbin(x, y, gridsize=50, cmap="viridis", mincnt=1)
    ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="1:1")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("PACE sif_radiance_678nm")
    ax.set_ylabel(f"TROPOMI SIF @ {REPORTING_WL_NM:g} nm")
    ax.set_title(
        f"Independent SIF  n={stats['n']:,}\n"
        f"bias={stats['bias']:.3g}  RMSE={stats['rmse']:.3g}  R={stats['r']:.3f}"
    )
    ax.legend(loc="upper left", fontsize=8)
    fig.colorbar(hb, ax=ax, label="count")

    ax = axes[1]
    bins = np.linspace(0, 10, 80)
    ax.hist(chi2_pace[np.isfinite(chi2_pace)], bins=bins, alpha=0.55, color="#2ca02c", label="PACE χ²", density=True)
    ax.hist(chi2_trop[np.isfinite(chi2_trop)], bins=bins, alpha=0.55, color="#ff7f0e", label="TROPOMI χ²", density=True)
    ax.set_xlabel("reduced χ²")
    ax.set_ylabel("density")
    ax.set_title("χ² on co-located pairs")
    ax.legend()

    fig.suptitle(f"cycle={cycle}: retrieval comparison", y=1.02)
    fig.tight_layout()
    fig.savefig(fig_dir / "sif_scatter_chi2.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.hist(y - x, bins=50, color="#4C72B0", edgecolor="white")
    ax.axvline(0, color="k", ls="--", lw=1)
    ax.set_xlabel("SIF$_T$ − SIF$_P$")
    ax.set_ylabel("count")
    ax.set_title(f"ΔSIF  mean={stats['bias']:.3g}")
    fig.tight_layout()
    fig.savefig(fig_dir / "sif_difference_hist.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_from_matchup_file(matchup_path: Path, fig_dir: Path) -> dict[str, float]:
    data = load_matchup_netcdf(matchup_path)
    cycle = str(data.get("cycle", ""))
    x = data["sif_pace_678nm"]
    y = data["sif_trop_678nm"]
    stats = scatter_stats(x, y)
    fig_dir.mkdir(parents=True, exist_ok=True)
    plot_sif_comparison(
        fig_dir,
        cycle=cycle,
        x=x,
        y=y,
        chi2_trop=data["chi2_trop"],
        chi2_pace=data["chi2_pace"],
        stats=stats,
    )
    print(
        f"Replotted from {matchup_path.name}: n={stats['n']:,}  "
        f"bias={stats['bias']:.4g}  RMSE={stats['rmse']:.4g}  R={stats['r']:.4f}"
    )
    return stats


def run_pipeline(args: argparse.Namespace) -> None:
    thresholds = DEFAULT_THRESHOLDS.copy()
    cycles = parse_cycles(args.cycle)
    cycle_label = "all" if len(cycles) > 1 else cycles[0]
    output_dir = OUTPUT_ROOT / cycle_label
    fig_dir = output_dir / "figs"
    matchup_path = output_dir / "matchup_data.nc"
    fig_dir.mkdir(parents=True, exist_ok=True)
    colocation_json = Path(args.colocation_json)

    basis_paths: dict[str, Path] = {}
    bases: dict[str, object] = {}
    for cycle in cycles:
        if args.basis and len(cycles) == 1:
            basis_path = Path(args.basis)
        elif args.basis and len(cycles) > 1:
            raise ValueError(
                "--basis override is only allowed with a single --cycle "
                "(cycle1 or cycle2), not --cycle all"
            )
        else:
            basis_path = default_basis_path(cycle)
        if not basis_path.is_file():
            raise FileNotFoundError(basis_path)
        print(f"Loading SIF basis weights [{cycle}] from {basis_path} @ {REPORTING_WL_NM:g} nm...")
        basis_paths[cycle] = basis_path
        bases[cycle] = load_sif_basis_weights(basis_path, REPORTING_WL_NM)

    jobs, pace_dir = build_jobs(
        cycles, args.orbit, args.max_swaths, colocation_json
    )
    match = match_swaths(jobs, pace_dir, thresholds)
    sel, dvza = apply_geo_filters(match, thresholds)
    joined = join_retrievals(match, bases, thresholds)

    if thresholds.get("require_both_retrieved", True):
        sel &= joined["both_ok"]
        print(f"After require_both_retrieved: {int(sel.sum()):,}")

    if int(sel.sum()) == 0:
        raise RuntimeError("No pairs left after filters — loosen THRESHOLDS")

    idx = save_matchup_netcdf(
        matchup_path,
        cycle_label=cycle_label,
        orbit=args.orbit,
        thresholds=thresholds,
        basis_paths=basis_paths,
        colocation_json=colocation_json,
        match=match,
        joined=joined,
        sel=sel,
        dvza=dvza,
    )

    rng = np.random.default_rng(args.seed)
    plot_match_qa(fig_dir, cycle=cycle_label, orbit=args.orbit, match=match, idx=idx, dvza=dvza)
    plot_example_spectra(
        fig_dir,
        cycle=cycle_label,
        orbit=args.orbit,
        match=match,
        joined=joined,
        idx=idx,
        rng=rng,
    )

    stats = scatter_stats(joined["sif_pace"][idx], joined["sif_trop_678"][idx])
    plot_sif_comparison(
        fig_dir,
        cycle=cycle_label,
        x=joined["sif_pace"][idx],
        y=joined["sif_trop_678"][idx],
        chi2_trop=joined["chi2_trop"][idx],
        chi2_pace=joined["chi2_pace"][idx],
        stats=stats,
    )

    rows = []
    example_i = np.sort(rng.choice(idx, size=min(N_EXAMPLE_SPECTRA, idx.size), replace=False))
    for i in example_i:
        rows.append(
            {
                "cycle": str(joined["cycles"][i]),
                "orbit": int(joined["orbits"][i]),
                "trop_scan": int(joined["trop_scan"][i]),
                "trop_pix": int(joined["trop_pix"][i]),
                "pace_scan": int(joined["pace_scan"][i]),
                "pace_pix": int(joined["pace_pix"][i]),
                "dist_km": round(float(match["dist_km"][i]), 3),
                "dt_min": round(float(match["dt_min"][i]), 2),
                "lat": round(float(match["pace_lat"][i]), 3),
                "lon": round(float(match["pace_lon"][i]), 3),
                "sif_trop_678nm": float(joined["sif_trop_678"][i]),
                "sif_pace_678nm": float(joined["sif_pace"][i]),
                "chi2_trop": float(joined["chi2_trop"][i]),
                "chi2_pace": float(joined["chi2_pace"][i]),
            }
        )
    pd.DataFrame(rows).to_csv(output_dir / "example_pairs.csv", index=False)

    print(
        f"\nDone. Final n={idx.size:,}  "
        f"SIF bias(T−P)={stats['bias']:.4g}  RMSE={stats['rmse']:.4g}  R={stats['r']:.4f}\n"
        f"  matchup: {matchup_path}\n"
        f"  figures: {fig_dir}/"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cycle",
        default="all",
        choices=("cycle1", "cycle2", "all"),
        help="Which TROPOMI retrieval cycle(s) to include (default: all)",
    )
    parser.add_argument("--orbit", type=int, default=None, help="Single orbit or all co-located")
    parser.add_argument(
        "--max-swaths",
        type=int,
        default=None,
        help="Smoke-test limit: max swaths per selected cycle",
    )
    parser.add_argument("--colocation-json", type=Path, default=COLOCATION_JSON)
    parser.add_argument(
        "--basis",
        default=None,
        help="Override row-aware basis h5 (single-cycle runs only)",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Regenerate scatter/χ² figures from saved matchup_data.nc",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cycles = parse_cycles(args.cycle)
    cycle_label = "all" if len(cycles) > 1 else cycles[0]
    output_dir = OUTPUT_ROOT / cycle_label
    matchup_path = output_dir / "matchup_data.nc"
    fig_dir = output_dir / "figs"

    if args.plots_only:
        if not matchup_path.is_file():
            raise FileNotFoundError(f"No saved matchup at {matchup_path}; run full pipeline first.")
        plot_from_matchup_file(matchup_path, fig_dir)
        return

    run_pipeline(args)


if __name__ == "__main__":
    main()
