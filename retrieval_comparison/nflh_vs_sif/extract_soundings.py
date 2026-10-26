"""Collect per-sounding SVD SIF and L2 nFLH for whole days, for nflh_vs_sif plots with any later filter.

    python retrieval_comparison/nflh_vs_sif/extract_soundings.py [--dates 20250115 20250415] [--workers 16]
        [--run-dir /home/zhe2/data/PACE/new_svd_retrieval_output_O2A]

For every retrieval granule of the date (run-dir/YYYY/MM/DD/interim_<granule>_svd_retrieval.nc) the matching
L2 AOP / BGC V3.2 and L1B files are read on the same scan/pixel grid (source_scan_index/source_pixel_index).
Kept: every sounding with a finite L2 nflh (i.e. AOP produced ocean-colour products there); no other cut,
so filtering (l2_flags, convergence, is_dark, χ², ...) is left to plot_nflh_vs_sif.py.

Columns: granule, scan, pixel, lat, lon, sza, vza, sif_678 (water-leaving SIF at the band nearest 678.2 nm,
sif_radiance_678nm), sif_678_sigma (posterior 1σ), converged, status_code, n_steps, reduced_chi2, dof,
is_dark, is_ocean, alpha, l2_flags, nflh (stored L2: 680 nm definition, clipped at 0), chlor_a, aot_865,
rtoa_678 (L1B TOA radiance at 678.18 nm = rhot_red · E0 · cos(SZA) / π / earth_sun_distance_correction).
Output: <out-dir>/<run name>/soundings_<date>.parquet (large; kept outside the repository).
"""
import argparse
import glob
import os
import re
from multiprocessing import Pool
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

K = "/kiwi-data/Data/satellite/PACE_OCI"
SIF_NM = 678.2                                   # rtoa_678: L1B band nearest this (678.18 nm)

p = argparse.ArgumentParser()
p.add_argument("--dates", nargs="+", default=["20250115", "20250415"])
p.add_argument("--run-dir", default="/home/zhe2/data/PACE/new_svd_retrieval_output_O2A")
p.add_argument("--aop-dir", default=f"{K}/L2_AOP_V3.2")
p.add_argument("--bgc-dir", default=f"{K}/L2_BGC_V3.2")
p.add_argument("--l1b-dir", default=f"{K}/L1B_V3")
p.add_argument("--out-dir", default="/home/zhe2/data/PACE/nflh_vs_sif_soundings")
p.add_argument("--workers", type=int, default=16)
A = p.parse_args()


def dated(root, g, name):
    return os.path.join(root, g[:4], g[4:6], g[6:8], name)


def filled(v, sl=slice(None)):
    return np.ma.filled(v[sl], np.nan).astype(np.float32)


def one(path):
    g = re.search(r"interim_(\d{8}T\d{6})_", os.path.basename(path)).group(1)
    try:
        with netCDF4.Dataset(path) as r:
            si = np.asarray(r["source_scan_index"][:]) - 1
            pi = np.asarray(r["source_pixel_index"][:]) - 1
            names = r.state_names_csv.split(",")
            k_sif, k_alpha = names.index("sif_ev1"), names.index("alpha")
            sif = filled(r["sif_radiance_678nm"])
            coeff = filled(r["sif_coeffs"])[0]
            var_sif = filled(r["S_posterior_diag"], (k_sif,))
            raw_alpha = filled(r["x_hat"], (k_alpha,))
            ret = dict(sif_678=sif, converged=np.ma.filled(r["converged"][:], 0).astype(np.uint8),
                       status_code=np.ma.filled(r["status_code"][:], -1).astype(np.int16),
                       n_steps=np.ma.filled(r["n_steps"][:], -1).astype(np.int16),
                       reduced_chi2=filled(r["reduced_chi2"]), dof=filled(r["dof"]),
                       is_dark=np.ma.filled(r["is_dark"][:], 0).astype(np.uint8),
                       is_ocean=np.ma.filled(r["is_ocean"][:], 0).astype(np.uint8),
                       lat=filled(r["latitude"]), lon=filled(r["longitude"]))
        # SIF at 678 = b678 · coeff, so its 1σ = |b678| · sqrt(var(coeff))
        ok = np.isfinite(sif) & np.isfinite(coeff) & (np.abs(coeff) > 1e-6)
        b678 = np.nanmedian(sif[ok] / coeff[ok]) if ok.any() else np.nan
        ret["sif_678_sigma"] = (np.abs(b678) * np.sqrt(np.clip(var_sif, 0, None))).astype(np.float32)
        ret["alpha"] = (10.0 / (1.0 + np.exp(-raw_alpha)) + 1.0).astype(np.float32)
        S, P = np.ix_(si, pi)

        with netCDF4.Dataset(dated(A.aop_dir, g, f"PACE_OCI.{g}.L2.OC_AOP.V3_2.nc")) as a:
            gd = a["geophysical_data"]
            nflh = filled(gd["nflh"])[S, P]
            keep = np.isfinite(nflh)
            ret["nflh"] = nflh
            ret["l2_flags"] = (np.ma.filled(gd["l2_flags"][:], 0).astype(np.int64) & 0xFFFFFFFF)[S, P].astype(np.uint32)
            ret["aot_865"] = filled(gd["aot_865"])[S, P]
        bgc = dated(A.bgc_dir, g, f"PACE_OCI.{g}.L2.OC_BGC.V3_2.nc")
        if os.path.exists(bgc):
            with netCDF4.Dataset(bgc) as b:
                ret["chlor_a"] = filled(b["geophysical_data"]["chlor_a"])[S, P]
        l1b = dated(A.l1b_dir, g, f"PACE_OCI.{g}.L1B.V3.nc")
        if os.path.exists(l1b):
            with netCDF4.Dataset(l1b) as l1:
                geo = l1["geolocation_data"]
                ret["sza"] = filled(geo["solar_zenith"])[S, P]
                ret["vza"] = filled(geo["sensor_zenith"])[S, P]
                # TOA radiance at the SIF band, as the retrieval builds it from rhot_red
                sbp = l1["sensor_band_parameters"]
                wl = np.ma.filled(sbp["red_wavelength"][:], np.nan)
                k = int(np.nanargmin(np.abs(wl - SIF_NM)))
                E0 = float(sbp["red_solar_irradiance"][k])
                es = float(l1.getncattr("earth_sun_distance_correction"))
                rhot = filled(l1["observation_data"]["rhot_red"], (k,))[S, P]
                ret["rtoa_678"] = (rhot * E0 * np.cos(np.radians(ret["sza"])) / np.pi / es).astype(np.float32)
        ss, pp = np.meshgrid(si, pi, indexing="ij")
        df = pd.DataFrame({"granule": g, "scan": ss[keep].astype(np.int16), "pixel": pp[keep].astype(np.int16),
                           **{k: v[keep] for k, v in ret.items()}})
        return g, df, int(keep.sum()), int(sif.size), None
    except Exception as e:                                             # noqa: BLE001
        return g, None, 0, 0, repr(e)


run = Path(A.run_dir).name
out_dir = Path(A.out_dir) / run
out_dir.mkdir(parents=True, exist_ok=True)
for date in A.dates:
    files = sorted(glob.glob(os.path.join(A.run_dir, date[:4], date[4:6], date[6:8], "interim_*_svd_retrieval*.nc")))
    if not files:
        print(f"{date}: no retrieval files under {A.run_dir}")
        continue
    with Pool(A.workers) as pool:
        res = pool.map(one, files)
    bad = [(g, e) for g, _, _, _, e in res if e]
    dfs = [d for _, d, _, _, e in res if d is not None and len(d)]
    df = pd.concat(dfs, ignore_index=True)
    df["granule"] = df.granule.astype("category")
    out = out_dir / f"soundings_{date}.parquet"
    df.to_parquet(out, index=False)
    n_all = sum(n for _, _, _, n, _ in res)
    print(f"{date}: {len(files)} granules ({len(bad)} failed), {n_all:,} soundings, {len(df):,} with L2 nflh "
          f"({(df.l2_flags == 0).sum():,} with l2_flags == 0) -> {out} ({out.stat().st_size / 1e6:.0f} MB)")
    for g, e in bad[:5]:
        print(f"   failed {g}: {e}")
