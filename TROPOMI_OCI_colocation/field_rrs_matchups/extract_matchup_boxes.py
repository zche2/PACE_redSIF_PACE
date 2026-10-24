"""Locate every field/PACE matchup on the L1B swath and extract a 5x5 box with all retrieval inputs.

    python TROPOMI_OCI_colocation/field_rrs_matchups/extract_matchup_boxes.py [fullmatchups_Rrs_GUI.csv]

Inputs: fullmatchups_Rrs_GUI.csv (default: main checkout) and, per granule, the L1B V3, L2 AOP V3.2
and L2 BGC V3.2 files in /kiwi-data/Data/satellite/PACE_OCI/<product>/YYYY/MM/DD/.

Centre pixel: the matchup's own PACE centre (nearest L1B pixel to PACE_lat/PACE_lon; matches to <60 m).
The CSV's PACE_center_line/pixel are subscene indices and are not used. The pixel nearest to the
field location is stored as well (it differs from the matchup centre in ~40 % of cases).

Pixel screening inside the 5x5 box (Bailey & Werdell 2006 style):
  - L2 flags (AOP l2_flags): ATMFAIL LAND HIGLINT MODGLINT HILT HISATZEN STRAYLIGHT CLDICE COCCOLITH
    HISOLZEN LOWLW NAVWARN NAVFAIL. This is the set that best reproduces the matchup tool's
    PACE_valid_pixels (72 % exact; the tool's full rule is not documented). The stricter Bailey &
    Werdell (2006) set, which also excludes CHLFAIL ABSAER MAXAERITER ATMWARN, is kept as flag_ok_bw2006.
  - VZA < 60 deg and SZA < 75 deg (L1B geolocation)
  - finite Rrs at all AOP bands
  valid = all of the above. The matchup-level protocol also needs n_valid >= 13 (half the box),
  |pixel time - field time| <= 3 h, and CV <= 0.15 (filtered Rrs CV over 400-600 nm).

Retrieval inputs (as global_svd_fit_pipeline builds them):
  Rtoa_<λ> = rhot_red * E0 * cos(SZA) / pi / earth_sun_distance_correction(L1B)   [W m-2 sr-1 um-1]
  for the red bands with 640 < λ < 780 nm; plus SZA, watermask (is_ocean), is_dark (max Rtoa <= 20),
  rtoa_complete (no missing red band, which the pipeline requires), nflh, chlor_a, aot_865, angstrom.
  nflh is the stored L2 value (V3.2: height 680 nm, base 664.6-667.0 and 710.5-713.0 nm, clipped at 0);
  nflh_678 is recomputed from Rrs with the V3.0 definition (height 678.2 nm, base 649.6-659.6 and
  710.5-719.2 nm, flh_offset 0), unclipped; nflh_680 recomputes the stored one, unclipped, as a check.
  The FLH attributes of each AOP file are kept in the summary (aop_flh_*).
  E0 per band is in red_bands_640_780.csv, the Earth-Sun factor per granule in the summary.

Rrs double check: per pixel, OCI Rrs is smoothed to the CSV field wavelengths with a 10 nm FWHM
Gaussian, then the Bailey & Werdell filtered mean (pixels within 1.5 sigma of the median) is compared
with the CSV's PACE_Rrs. Only field wavelengths bracketed by OCI bands at most 5 nm apart are kept: the
CSV values below 346 nm (extrapolated) and at 590-610 nm (interpolated across the 588-613 nm gap of the
L2 Rrs bands) have no OCI measurement behind them.

Outputs (this folder; *.csv is gitignored):
  matchup_boxes_5x5.csv   one row per box pixel: indices, geometry, time, flags, validity, all inputs
  matchup_summary.csv     one row per matchup: centres, protocol checks, valid-pixel means, Rrs check
  red_bands_640_780.csv   red band index, wavelength and solar irradiance E0 of the retrieval window
"""
import re
import sys
import warnings
from datetime import datetime, timedelta, timezone
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore", category=RuntimeWarning)
CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path("/home/zhe2/FraLab/PACE_redSIF_PACE/TROPOMI_OCI_colocation/fullmatchups_Rrs_GUI.csv")
OUT = Path(__file__).resolve().parent
K = "/kiwi-data/Data/satellite/PACE_OCI"
HALF = 2                                   # 5x5 box
VZA_MAX, SZA_MAX, DT_MAX_H = 60.0, 75.0, 3.0
LAM_MIN, LAM_MAX = 640.0, 780.0            # retrieval window (strict inequality, as the pipeline)
DARK_MAX = 20.0
FLAGS = ("ATMFAIL LAND HIGLINT MODGLINT HILT HISATZEN STRAYLIGHT CLDICE COCCOLITH HISOLZEN LOWLW "
         "NAVWARN NAVFAIL").split()
BW_FLAGS = ("ATMFAIL LAND HIGLINT HILT HISATZEN STRAYLIGHT CLDICE COCCOLITH HISOLZEN LOWLW CHLFAIL "
            "NAVWARN ABSAER MAXAERITER ATMWARN NAVFAIL").split()
# nFLH definitions (l2gen flh_base_wavelengths / flh_height_wavelength). V3.1 and V3.2 use FLH_680, and
# their stored nflh is clipped at 0 (valid_min); V3.0 used FLH_678, the one the SIF comparison assumes.
FLH_678 = ([649.599976, 650.900024, 652.099976, 653.299988, 654.599976, 655.799988, 657.099976, 658.299988,
            659.599976, 710.5, 711.799988, 713.0, 714.299988, 716.799988, 719.200012], 678.20001)
FLH_680 = ([664.564026, 665.794983, 667.0, 710.5, 711.799988, 713.0], 680.0)
CHECK_FWHM = 10.0
MAX_GAP = 5.0                              # nm; OCI L2 Rrs spacing is 2-3 nm, with no bands at 588-613 nm


def path(sub, tag, g):
    return f"{K}/{sub}/{g[:4]}/{g[4:6]}/{g[6:8]}/PACE_OCI.{g}.{tag}.nc"


def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def filled(v, sl=slice(None)):
    return np.ma.filled(v[sl], np.nan).astype(float)


def gauss_matrix(lam, centres, fwhm):
    s = fwhm / 2.355
    W = np.exp(-0.5 * ((lam[None, :] - centres[:, None]) / s) ** 2)
    return W / W.sum(1, keepdims=True)


def measured(lam, centres, max_gap=MAX_GAP):
    """Field wavelengths bracketed by two OCI bands at most max_gap apart (no extrapolation, no gap filling)."""
    j = np.searchsorted(lam, centres)
    lo, hi = lam[np.clip(j - 1, 0, len(lam) - 1)], lam[np.clip(j, 0, len(lam) - 1)]
    exact = np.isin(centres, lam)
    return exact | ((j > 0) & (j < len(lam)) & (hi - lo <= max_gap))


def nflh_from_rrs(rrs, lam, F0, base, height):
    """l2gen FLH, unclipped: nLw = Rrs * F0 [W m-2 um-1 sr-1]; least-squares line through the bands nearest
    the base wavelengths, evaluated at the height wavelength, subtracted from nLw at the band nearest it.
    Reproduces the stored V3.2 nflh (FLH_680) to the packing precision wherever that is not clipped at 0."""
    X = rrs * F0
    bi = [int(np.argmin(np.abs(lam - w))) for w in base]
    hi = int(np.argmin(np.abs(lam - height)))
    A = np.column_stack([np.ones(len(bi)), lam[bi]])
    w = np.array([1.0, height]) @ np.linalg.pinv(A)      # LSQ line at the height = w . nLw(base)
    return X[..., hi] - X[..., bi] @ w


def bw_filtered(x):
    """Bailey & Werdell filtered mean / std over pixels (rows) for each column."""
    med = np.nanmedian(x, 0)
    sd = np.nanstd(x, 0, ddof=1)
    keep = np.abs(x - med) <= 1.5 * sd
    xm = np.where(keep, x, np.nan)
    return np.nanmean(xm, 0), np.nanstd(xm, 0, ddof=1)


df = pd.read_csv(CSV, low_memory=False)
df["granule"] = [re.search(r"(\d{8}T\d{6})", str(s)).group(1) for s in df["PACE_parent_file"]]
df["field_time"] = pd.to_datetime(df["FIELD_datetime"], utc=True)
lam_csv = np.array(sorted(float(m.group(1)) for c in df.columns if (m := re.fullmatch(r"PACE_Rrs_(\d+)", c))))
csv_rrs_cols = [f"PACE_Rrs_{int(w)}" for w in lam_csv]
chk = None                                 # CSV wavelengths with OCI Rrs support, set from the first AOP file

pix_rows, sum_rows = [], []
band_info = None
for n, (g, grp) in enumerate(df.groupby("granule"), 1):
    print(f"[{n:2d}/{df.granule.nunique()}] {g}  ({len(grp)} matchups)", flush=True)
    l1 = netCDF4.Dataset(path("L1B_V3", "L1B.V3", g))
    aop = netCDF4.Dataset(path("L2_AOP_V3.2", "L2.OC_AOP.V3_2", g))
    bgc = netCDF4.Dataset(path("L2_BGC_V3.2", "L2.OC_BGC.V3_2", g))
    for d in (l1, aop, bgc):
        d.set_auto_mask(True)
    geo = l1["geolocation_data"]
    lat, lon = filled(geo["latitude"]), filled(geo["longitude"])
    nS, nP = lat.shape
    ok = np.isfinite(lat) & np.isfinite(lon)
    idx = np.argwhere(ok)
    tree = cKDTree(xyz(lat[ok], lon[ok]))
    # scan times: "seconds since YYYY-MM-DD"
    tv = l1["scan_line_attributes"]["time"]
    t0 = datetime.strptime(tv.units.split("since")[1].strip()[:10], "%Y-%m-%d").replace(tzinfo=timezone.utc)
    scan_time = filled(tv)
    # red bands of the retrieval window, sorted by wavelength as the pipeline does
    wl = filled(l1["sensor_band_parameters"]["red_wavelength"])
    E0 = filled(l1["sensor_band_parameters"]["red_solar_irradiance"])
    order = np.argsort(wl)
    win = order[(wl[order] > LAM_MIN) & (wl[order] < LAM_MAX)]
    es = float(l1.getncattr("earth_sun_distance_correction"))
    if band_info is None:
        band_info = pd.DataFrame({"red_band_index": win, "wavelength_nm": wl[win], "E0_W_m2_um": E0[win]})
    sbp = aop["sensor_band_parameters"]
    lam_aop = filled(sbp["wavelength_3d"])
    # float wavelength and F0 of each Rrs band (wavelength_3d is the rounded integer)
    ref = np.array([np.flatnonzero(filled(sbp["wavelength"]) == w)[0] for w in lam_aop])
    fwave_aop, F0_aop = filled(sbp["fwave"])[ref], filled(sbp["F0"])[ref]
    ip = aop["processing_control"]["input_parameters"]
    flh_attrs = {k: ip.getncattr(k) for k in ("flh_offset", "flh_base_wavelengths", "flh_height_wavelength")}
    fl = aop["geophysical_data"]["l2_flags"]
    names = fl.flag_meanings.split()
    masks = np.array(fl.flag_masks).astype(np.int64) & 0xFFFFFFFF
    bad = 0
    for f in FLAGS:
        bad |= int(masks[names.index(f)])
    bad_bw = 0
    for f in BW_FLAGS:
        bad_bw |= int(masks[names.index(f)])
    Gm = gauss_matrix(lam_aop, lam_csv, CHECK_FWHM)
    if chk is None:
        chk = measured(lam_aop, lam_csv)

    for i, r in grp.iterrows():
        dc, kc = tree.query(xyz(np.array([r.PACE_lat]), np.array([r.PACE_lon])))
        df_, kf = tree.query(xyz(np.array([r.FIELD_lat]), np.array([r.FIELD_lon])))
        s, p = idx[kc[0]]
        sf, pf = idx[kf[0]]
        s0, s1, p0, p1 = max(s - HALF, 0), min(s + HALF + 1, nS), max(p - HALF, 0), min(p + HALF + 1, nP)
        box = (slice(s0, s1), slice(p0, p1))
        sza, vza = filled(geo["solar_zenith"], box), filled(geo["sensor_zenith"], box)
        saa, vaa = filled(geo["solar_azimuth"], box), filled(geo["sensor_azimuth"], box)
        wm = np.ma.filled(geo["watermask"][box], -1).astype(int)
        flags = np.ma.filled(fl[box], 0).astype(np.int64) & 0xFFFFFFFF
        gd = aop["geophysical_data"]
        rrs = filled(gd["Rrs"], (slice(s0, s1), slice(p0, p1), slice(None)))
        nflh, aot, ang = filled(gd["nflh"], box), filled(gd["aot_865"], box), filled(gd["angstrom"], box)
        nflh678 = nflh_from_rrs(rrs, fwave_aop, F0_aop, *FLH_678)
        nflh680 = nflh_from_rrs(rrs, fwave_aop, F0_aop, *FLH_680)
        chl = filled(bgc["geophysical_data"]["chlor_a"], box)
        rhot_all = filled(l1["observation_data"]["rhot_red"], (slice(None), slice(s0, s1), slice(p0, p1)))
        rhot = rhot_all[win]
        rtoa = rhot * E0[win][:, None, None] * np.cos(np.radians(sza))[None] / np.pi / es
        ptime = [t0 + timedelta(seconds=float(scan_time[ss])) for ss in range(s0, s1)]

        flag_ok = (flags & bad) == 0
        flag_ok_bw = (flags & bad_bw) == 0
        geom_ok = (vza < VZA_MAX) & (sza < SZA_MAX)
        rrs_ok = np.isfinite(rrs).all(axis=2)
        valid = flag_ok & geom_ok & rrs_ok
        rtoa_complete = np.isfinite(rhot_all).all(axis=0)
        is_dark = np.nanmax(rtoa, axis=0) <= DARK_MAX

        for a in range(s1 - s0):
            dt_min = (ptime[a] - r.field_time).total_seconds() / 60.0
            for b in range(p1 - p0):
                row = {"matchup_id": i, "granule": g, "scan": s0 + a, "pixel": p0 + b,
                       "d_scan": s0 + a - s, "d_pixel": p0 + b - p, "is_centre": (s0 + a == s) and (p0 + b == p),
                       "lat": lat[s0 + a, p0 + b], "lon": lon[s0 + a, p0 + b],
                       "pixel_time_utc": ptime[a].isoformat(), "dt_minutes": dt_min,
                       "sza": sza[a, b], "vza": vza[a, b], "saa": saa[a, b], "vaa": vaa[a, b],
                       "watermask": wm[a, b], "l2_flags": int(flags[a, b]),
                       "flag_ok": bool(flag_ok[a, b]), "flag_ok_bw2006": bool(flag_ok_bw[a, b]),
                       "geom_ok": bool(geom_ok[a, b]),
                       "rrs_finite": bool(rrs_ok[a, b]), "valid": bool(valid[a, b]),
                       "rtoa_complete": bool(rtoa_complete[a, b]), "is_dark": bool(is_dark[a, b]),
                       "nflh": nflh[a, b], "nflh_678": nflh678[a, b], "nflh_680": nflh680[a, b],
                       "chlor_a": chl[a, b], "aot_865": aot[a, b], "angstrom": ang[a, b]}
                row.update({f"Rtoa_{w:.2f}": rtoa[q, a, b] for q, w in enumerate(wl[win])})
                row.update({f"Rrs_{int(w)}": rrs[a, b, q] for q, w in enumerate(lam_aop)})
                pix_rows.append(row)

        # matchup summary
        v = valid.ravel()
        nv = int(v.sum())
        dt_centre = (ptime[s - s0] - r.field_time).total_seconds() / 60.0
        summ = {"matchup_id": i, "granule": g, "FIELD_datetime": r.FIELD_datetime, "FIELD_lat": r.FIELD_lat,
                "FIELD_lon": r.FIELD_lon, "FIELD_experiment": r.FIELD_experiment, "FIELD_cruise": r.FIELD_cruise,
                "data_type": r.data_type,
                "centre_scan": s, "centre_pixel": p, "centre_dist_to_PACE_latlon_m": dc[0] * 6.371e6,
                "fieldnearest_scan": sf, "fieldnearest_pixel": pf, "field_to_nearest_pixel_km": df_[0] * 6.371e3,
                "fieldnearest_in_box": abs(sf - s) <= HALF and abs(pf - p) <= HALF,
                "centre_time_utc": ptime[s - s0].isoformat(), "dt_minutes": dt_centre,
                "within_3h": abs(dt_centre) <= DT_MAX_H * 60,
                "centre_sza": sza[s - s0, p - p0], "centre_vza": vza[s - s0, p - p0],
                "earth_sun_distance_correction_L1B": es,
                "n_box": int(valid.size), "n_flag_ok": int(flag_ok.sum()), "n_geom_ok": int(geom_ok.sum()),
                "n_valid": nv, "n_valid_bw2006": int((flag_ok_bw & geom_ok & rrs_ok).sum()),
                "PACE_valid_pixels_csv": r.PACE_valid_pixels, **{f"aop_{k}": v for k, v in flh_attrs.items()}}
        if nv > 0:
            C = rrs.reshape(-1, rrs.shape[2])[v] @ Gm.T          # per pixel, smoothed to the CSV grid
            fm, fsd = bw_filtered(C)
            csv = r[csv_rrs_cols].to_numpy(float)
            m = chk & np.isfinite(csv) & np.isfinite(fm)
            rel = np.abs(fm[m] - csv[m]) / np.maximum(np.abs(csv[m]), 1e-6)
            bwin = (lam_csv >= 400) & (lam_csv <= 600)
            cv = np.nanmedian(fsd[bwin] / np.abs(fm[bwin])) if bwin.any() else np.nan
            summ.update({"filtered_cv_400_600": cv,
                         "rrs_check_median_reldiff": np.median(rel), "rrs_check_p90_reldiff": np.percentile(rel, 90),
                         "rrs_check_median_absdiff": np.median(np.abs(fm[m] - csv[m]))})
            for name, arr in (("nflh", nflh), ("nflh_678", nflh678), ("nflh_680", nflh680), ("chlor_a", chl),
                              ("aot_865", aot), ("angstrom", ang)):
                summ[f"{name}_mean_valid"] = np.nanmean(arr.ravel()[v])
                summ[f"{name}_centre"] = arr[s - s0, p - p0]
            rt = rtoa.reshape(rtoa.shape[0], -1)[:, v]
            summ.update({f"Rtoa_mean_{w:.2f}": np.nanmean(rt[q]) for q, w in enumerate(wl[win])})
            summ.update({f"PACE_Rrs_rep_{int(w)}": fm[q] for q, w in enumerate(lam_csv) if chk[q]})
        summ["passes_protocol"] = bool(summ["within_3h"] and nv >= 13 and
                                       summ.get("filtered_cv_400_600", np.inf) <= 0.15)
        sum_rows.append(summ)
    for d in (l1, aop, bgc):
        d.close()

pix = pd.DataFrame(pix_rows)
summ = pd.DataFrame(sum_rows).sort_values("matchup_id")
pix.to_csv(OUT / "matchup_boxes_5x5.csv", index=False)
summ.to_csv(OUT / "matchup_summary.csv", index=False)
band_info.to_csv(OUT / "red_bands_640_780.csv", index=False)

print(f"\nwrote {OUT}/matchup_boxes_5x5.csv ({len(pix)} pixels x {pix.shape[1]} columns), "
      f"matchup_summary.csv ({len(summ)} matchups), red_bands_640_780.csv ({len(band_info)} bands)")
print(f"centre found from PACE_lat/lon: max distance {summ.centre_dist_to_PACE_latlon_m.max():.1f} m")
print(f"field point to nearest pixel: median {summ.field_to_nearest_pixel_km.median():.2f} km; "
      f"nearest pixel inside the 5x5 box for {summ.fieldnearest_in_box.mean() * 100:.0f} % of matchups")
print(f"|dt| <= 3 h: {summ.within_3h.sum()} / {len(summ)};  max |dt| {summ.dt_minutes.abs().max():.0f} min")
print(f"n_valid: median {summ.n_valid.median():.0f};  matchups with n_valid >= 13: {(summ.n_valid >= 13).sum()};  "
      f"n_valid == CSV PACE_valid_pixels: {(summ.n_valid == summ.PACE_valid_pixels_csv).mean() * 100:.0f} %")
print(f"passes protocol (|dt|<=3h, n_valid>=13, CV<=0.15): {summ.passes_protocol.sum()} / {len(summ)}")
ok = summ.rrs_check_median_reldiff.notna()
for k in ("flh_base_wavelengths", "flh_height_wavelength"):
    print(f"AOP {k}: " + "; ".join(f"'{v}' x{n}" for v, n in summ[f"aop_{k}"].value_counts().items()))
vp = pix[pix.valid & pix.nflh.notna()]
pos = vp.nflh_680 > 0.001
print(f"nflh_680 (recomputed) vs stored nflh, valid pixels with nflh_680 > 0: median |diff| "
      f"{(vp.nflh_680 - vp.nflh)[pos].abs().median():.2e}, max {(vp.nflh_680 - vp.nflh)[pos].abs().max():.2e}; "
      f"stored nflh clipped at 0 for {(vp.nflh_680 < 0).sum()} of {len(vp)} valid pixels")
print(f"nflh_678 vs nflh_680 (valid pixels): median ratio {(vp.nflh_678 / vp.nflh_680)[pos].median():.3f}, "
      f"median diff {(vp.nflh_678 - vp.nflh_680).median():+.4f} W m-2 um-1 sr-1")
print(f"CSV wavelengths kept (OCI Rrs support): {chk.sum()} / {len(lam_csv)}; dropped {lam_csv[~chk].astype(int).tolist()}")
print(f"Rrs check vs CSV PACE_Rrs (supported wavelengths): median of per-matchup median rel. diff "
      f"{summ.rrs_check_median_reldiff[ok].median() * 100:.2f} %, 90th pct {summ.rrs_check_median_reldiff[ok].quantile(0.9) * 100:.2f} %, "
      f"worst {summ.rrs_check_median_reldiff[ok].max() * 100:.1f} %")
