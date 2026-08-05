"""Denoise TROPOMI → convolve to OCI bands; load measured OCI Lt."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import netCDF4 as nc
import numpy as np

from .config import Config
from .radiance import fill_nan_1d, tropomi_mol_to_mW
from .snr import load_pace_snr, snr_sigma_from_L


def denoise_with_svd(spec_on_svd_wl: np.ndarray, trop_svd: dict, n_pc: int) -> np.ndarray:
    mean = trop_svd["mean"]
    Vt = trop_svd["Vt"]
    k = min(int(n_pc), Vt.shape[0])
    y = np.asarray(spec_on_svd_wl, dtype=float)
    xc = y - mean
    coeffs = xc @ Vt[:k].T
    return mean + coeffs @ Vt[:k]


def tropomi_to_oci(
    lt_native: np.ndarray,
    gp: int,
    n_pc: int,
    *,
    trop_svd: dict,
    rsr_bundle: dict,
    pace_snr: dict,
    add_noise: bool = True,
    rng: np.random.Generator | None = None,
) -> dict:
    wl_svd = trop_svd["wl"]
    trop_wl_gp = rsr_bundle["trop_wl"][gp]
    K = rsr_bundle["K"][gp]
    oci_bands = rsr_bundle["oci_bands"]

    lt_native = fill_nan_1d(lt_native)
    m = np.isfinite(trop_wl_gp)
    order = np.argsort(trop_wl_gp[m])
    spec_svd = np.interp(
        wl_svd, trop_wl_gp[m][order], lt_native[m][order], left=np.nan, right=np.nan
    )
    spec_svd = fill_nan_1d(spec_svd)
    lt_den = denoise_with_svd(spec_svd, trop_svd, n_pc)
    spec_native_den = fill_nan_1d(
        np.interp(trop_wl_gp, wl_svd, lt_den, left=np.nan, right=np.nan)
    )
    lt_oci = K @ spec_native_den

    if pace_snr.get("coeffs") is not None and np.array_equal(
        np.asarray(pace_snr["coeffs"]["wavelength"]), oci_bands
    ):
        c1 = np.asarray(pace_snr["coeffs"]["c1"], dtype=float)
        c2 = np.asarray(pace_snr["coeffs"]["c2"], dtype=float)
    else:
        tab = pace_snr["table"]
        c1 = np.interp(
            oci_bands, tab["wavelength"], tab["c1"], left=tab["c1"][0], right=tab["c1"][-1]
        )
        c2 = np.interp(
            oci_bands, tab["wavelength"], tab["c2"], left=tab["c2"][0], right=tab["c2"][-1]
        )

    sigma = snr_sigma_from_L(lt_oci, c1, c2)
    snr = np.divide(
        np.maximum(lt_oci, 0.0),
        sigma,
        out=np.full_like(lt_oci, np.nan),
        where=sigma > 0,
    )

    if add_noise:
        if rng is None:
            rng = np.random.default_rng()
        noise = rng.normal(0.0, sigma)
        lt_oci_noisy = lt_oci + noise
    else:
        noise = np.zeros_like(lt_oci)
        lt_oci_noisy = lt_oci.copy()

    return {
        "oci_bands": oci_bands,
        "lt_den_svd": lt_den,
        "lt_oci": lt_oci,
        "lt_oci_noisy": lt_oci_noisy,
        "sigma": sigma,
        "snr": snr,
        "noise": noise,
        "gp": int(gp),
        "n_pc": int(n_pc),
    }


def _simulate_unique_tropomi(
    matches: dict,
    trop_svd: dict,
    rsr_bundles: dict[str, dict],
    cfg: Config,
    log=print,
) -> dict:
    n_match = int(matches["dist_km"].size)
    trop_scan = np.asarray(matches["trop_scan"], dtype=np.int64)
    trop_pix = np.asarray(matches["trop_pix"], dtype=np.int64)
    trop_paths = np.asarray(matches["tropomi_path"])

    # group by TROPOMI file
    by_file: dict[str, list[int]] = defaultdict(list)
    for i in range(n_match):
        by_file[str(trop_paths[i])].append(i)

    # determine oci bands from first bundle (assumes same band window)
    first_bundle = rsr_bundles[next(iter(rsr_bundles))]
    oci_bands = np.asarray(first_bundle["oci_bands"], dtype=np.float64)
    n_oci = oci_bands.size
    pace_snr = load_pace_snr(cfg, oci_bands)

    lt_oci = np.full((n_match, n_oci), np.nan, dtype=np.float32)
    lt_oci_noisy = np.full((n_match, n_oci), np.nan, dtype=np.float32)
    sigma = np.full((n_match, n_oci), np.nan, dtype=np.float32)
    trop_gp = np.zeros(n_match, dtype=np.int32)

    rng = np.random.default_rng(cfg.noise_seed)
    n_pc = cfg.svd.n_pc

    for tp, idxs in by_file.items():
        bundle = rsr_bundles[tp]
        idxs_arr = np.asarray(idxs, dtype=np.int64)
        pairs = np.stack([trop_scan[idxs_arr], trop_pix[idxs_arr]], axis=1)
        uniq, inv = np.unique(pairs, axis=0, return_inverse=True)
        log(f"  simulate {Path(tp).name}: matches={len(idxs):,} unique={uniq.shape[0]:,}")

        oci_u = np.full((uniq.shape[0], n_oci), np.nan, dtype=np.float32)
        noisy_u = np.full((uniq.shape[0], n_oci), np.nan, dtype=np.float32)
        sigma_u = np.full((uniq.shape[0], n_oci), np.nan, dtype=np.float32)
        gp_u = np.zeros(uniq.shape[0], dtype=np.int32)
        key_to_u = {(int(s), int(p)): i for i, (s, p) in enumerate(uniq)}

        with nc.Dataset(tp) as ds:
            trop_wl_all = np.asarray(
                ds["BAND5_RADIANCE"]["STANDARD_MODE"]["INSTRUMENT"]["nominal_wavelength"][
                    0, :, :
                ],
                dtype=float,
            )
            rad = ds["BAND5_RADIANCE"]["STANDARD_MODE"]["OBSERVATIONS"]["radiance"]
            for iscan in np.unique(uniq[:, 0]):
                mol_scan = np.asarray(rad[0, int(iscan), :, :], dtype=float)
                mol_scan = np.where(
                    np.isfinite(mol_scan) & (mol_scan < 1e30), mol_scan, np.nan
                )
                pix_this = uniq[uniq[:, 0] == iscan, 1]
                for ipix in pix_this:
                    ipix = int(ipix)
                    u = key_to_u[(int(iscan), ipix)]
                    gp = int(np.clip(ipix, 0, bundle["K"].shape[0] - 1))
                    wl_gp = trop_wl_all[gp]
                    lt = tropomi_mol_to_mW(mol_scan[gp], wl_gp)
                    out = tropomi_to_oci(
                        lt,
                        gp,
                        n_pc,
                        trop_svd=trop_svd,
                        rsr_bundle=bundle,
                        pace_snr=pace_snr,
                        add_noise=cfg.add_noise,
                        rng=rng,
                    )
                    oci_u[u] = out["lt_oci"]
                    noisy_u[u] = out["lt_oci_noisy"]
                    sigma_u[u] = out["sigma"]
                    gp_u[u] = gp

        lt_oci[idxs_arr] = oci_u[inv]
        lt_oci_noisy[idxs_arr] = noisy_u[inv]
        sigma[idxs_arr] = sigma_u[inv]
        trop_gp[idxs_arr] = gp_u[inv]

    return {
        "oci_bands": oci_bands.astype(np.float32),
        "lt_oci": lt_oci,
        "lt_oci_noisy": lt_oci_noisy,
        "sigma": sigma,
        "trop_gp": trop_gp,
    }


def _load_measured_oci(matches: dict, oci_bands: np.ndarray, log=print) -> tuple[np.ndarray, np.ndarray]:
    n_match = int(matches["dist_km"].size)
    n_oci = oci_bands.size
    lt_oci_meas = np.full((n_match, n_oci), np.nan, dtype=np.float32)
    mu_scale = np.full(n_match, np.nan, dtype=np.float32)

    # group by pace file for open-once
    by_pace: dict[str, list[int]] = defaultdict(list)
    for i in range(n_match):
        by_pace[str(matches["pace_path"][i])].append(i)

    for ppath, idxs in by_pace.items():
        with nc.Dataset(ppath) as ds_p:
            esd = float(getattr(ds_p, "earth_sun_distance_correction", 1.0))
            wl_oci = np.asarray(ds_p["sensor_band_parameters"]["red_wavelength"][:], dtype=float)
            f0 = np.asarray(ds_p["sensor_band_parameters"]["red_solar_irradiance"][:], dtype=float)
            rhot_var = ds_p["observation_data"]["rhot_red"]
            for j in idxs:
                pscan = int(matches["pace_scan"][j])
                ppix = int(matches["pace_pix"][j])
                sza_o = float(matches["pace_sza"][j])
                sza_t = float(matches["trop_sza"][j])
                mu_o = max(float(np.cos(np.deg2rad(sza_o))), 1e-6)
                mu_t = max(float(np.cos(np.deg2rad(sza_t))), 1e-6)
                mu_scale[j] = mu_t / mu_o
                rhot = np.asarray(rhot_var[:, pscan, ppix], dtype=float)
                rhot = np.where(rhot > 0, rhot, np.nan)
                lt_full = rhot * f0 * mu_o / (np.pi * esd)
                m = np.isfinite(wl_oci) & np.isfinite(lt_full)
                if m.sum() < 3:
                    continue
                order = np.argsort(wl_oci[m])
                lt_oci_meas[j] = np.interp(
                    oci_bands, wl_oci[m][order], lt_full[m][order], left=np.nan, right=np.nan
                )
        log(f"  measured OCI from {Path(ppath).name}: n={len(idxs):,}")

    return lt_oci_meas, mu_scale


def run_simulate(
    cfg: Config,
    matches: dict,
    trop_svd: dict,
    rsr_bundles: dict[str, dict],
    log=print,
) -> dict:
    log("[simulate] TROPOMI → OCI")
    sim = _simulate_unique_tropomi(matches, trop_svd, rsr_bundles, cfg, log=log)
    log("[simulate] load measured OCI Lt")
    lt_meas, mu_scale = _load_measured_oci(matches, sim["oci_bands"], log=log)
    sim["lt_oci_meas"] = lt_meas
    sim["mu_trop_over_mu_oci"] = mu_scale
    sim["n_pc"] = cfg.svd.n_pc
    sim["add_noise"] = cfg.add_noise
    sim["noise_seed"] = cfg.noise_seed
    return sim
