"""Matchups where PACE SIF is high while PACE L2 nflh is low: where they are, and their L1B TOA radiance and
L2 Rrs over 600-800 nm.

    python TROPOMI_OCI_colocation/field_rrs_matchups/plot_sif_outliers.py

Selection (from insitu_vs_pace.csv, the 160 matchups of compare_insitu_nflh_sif.py):
    box-mean SIF − box-mean nflh > EXCESS_MIN and box-mean nflh < NFLH_MAX
A reference matchup (dark water, nflh near the outliers' median, SIF ≈ nflh) is plotted for contrast.
TOA radiance = rhot_red · E0 · cos(SZA) / π / earth_sun_distance_correction, read from the L1B for the box.
Rrs = L2 AOP V3.2 (no bands at 588-613 nm), in-situ FIELD_Rrs in green.

Outputs: out/11_sif_outlier_map.png, out/12_sif_outlier_spectra.png
"""
import re
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import netCDF4
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
CSV = Path("/home/zhe2/FraLab/PACE_redSIF_PACE/TROPOMI_OCI_colocation/fullmatchups_Rrs_GUI.csv")
K = "/kiwi-data/Data/satellite/PACE_OCI"
EXCESS_MIN, NFLH_MAX = 0.15, 0.3
LAM = (600.0, 800.0)

cmp_ = pd.read_csv(HERE / "insitu_vs_pace.csv")
pix = pd.read_csv(HERE / "matchup_sif_nPC15_npoly3_O2A.csv")
box = pd.read_csv(HERE / "matchup_boxes_5x5.csv")
orig = pd.read_csv(CSV, low_memory=False)
vc = pix[pix.valid & pix.converged]
agg = vc.groupby("matchup_id").agg(chi2=("reduced_chi2", "median"), alpha=("alpha", "median"),
                                   dark=("is_dark", "mean"), chl=("chlor_a", "mean"), aot=("aot_865", "mean"),
                                   sza=("sza", "mean"), vza=("vza", "mean"), lat=("lat", "mean"),
                                   lon=("lon", "mean"), granule=("granule", "first"))
u = cmp_[cmp_.used].merge(agg, left_on="matchup_id", right_index=True)
u["excess"] = u.sif_mean - u.nflh_mean
out = u[(u.excess > EXCESS_MIN) & (u.nflh_mean < NFLH_MAX)].sort_values("excess", ascending=False)
ref_pool = u[(u.dark == 1) & (u.excess.abs() < 0.02)]
ref = ref_pool.iloc[(ref_pool.nflh_mean - out.nflh_mean.median()).abs().argsort().iloc[0]]
print(f"{len(out)} outliers (SIF − nflh > {EXCESS_MIN}, nflh < {NFLH_MAX}) of {len(u)} matchups:")
print(out[["matchup_id", "data_type", "FIELD_experiment", "lat", "lon", "nflh680_insitu", "nflh_mean", "sif_mean",
           "chi2", "alpha", "dark", "chl", "aot", "sza", "vza"]].round(3).to_string(index=False))
print(f"reference: matchup {ref.matchup_id} ({ref.data_type}, {ref.FIELD_experiment}), nflh {ref.nflh_mean:.3f}, "
      f"SIF {ref.sif_mean:.3f}, χ² {ref.chi2:.2f}")

# ------------------------------------------------------------------ map
clusters = out.assign(key=list(zip(out.lat.round(0), out.lon.round(0)))).groupby("key")
fig = plt.figure(figsize=(6 + 5 * len(clusters), 5.2))
gs = fig.add_gridspec(1, 1 + len(clusters), width_ratios=[1.3] + [1] * len(clusters))
ax = fig.add_subplot(gs[0], projection=ccrs.PlateCarree())
ax.set_extent([-128, -60, 20, 52])
ax.add_feature(cfeature.LAND, color="0.93")
ax.add_feature(cfeature.LAKES, facecolor="white", edgecolor="0.6", lw=0.3)
ax.coastlines(lw=0.4)
sc = ax.scatter(u.lon, u.lat, c=u.excess, cmap="RdBu_r", vmin=-0.3, vmax=0.3, s=14, edgecolor="0.3", lw=0.3,
                transform=ccrs.PlateCarree(), zorder=5)
ax.scatter(out.lon, out.lat, s=60, facecolor="none", edgecolor="k", lw=1.0, transform=ccrs.PlateCarree(), zorder=6)
ax.set_title(f"{len(u)} matchups, colour = PACE SIF − nflh; circled: {len(out)} outliers\n"
             f"(matchups outside the US are off this map)", fontsize=9)
fig.colorbar(sc, ax=ax, shrink=0.7, label="box-mean SIF − nflh [W m$^{-2}$ µm$^{-1}$ sr$^{-1}$]")
for j, (key, g) in enumerate(clusters, 1):
    a = fig.add_subplot(gs[j], projection=ccrs.PlateCarree())
    allp = u[(u.lat - key[0]).abs().lt(1.0) & (u.lon - key[1]).abs().lt(1.0)]
    pad = 0.12
    a.set_extent([allp.lon.min() - pad, allp.lon.max() + pad, allp.lat.min() - pad, allp.lat.max() + pad])
    a.add_feature(cfeature.NaturalEarthFeature("physical", "land", "10m"), color="0.93")
    a.add_feature(cfeature.NaturalEarthFeature("physical", "lakes", "10m"), facecolor="white", edgecolor="0.5",
                  lw=0.4)
    a.coastlines("10m", lw=0.4)
    a.scatter(allp.lon, allp.lat, c=allp.excess, cmap="RdBu_r", vmin=-0.3, vmax=0.3, s=30, edgecolor="0.3",
              lw=0.3, transform=ccrs.PlateCarree(), zorder=5)
    for _, r in g.iterrows():
        a.scatter(r.lon, r.lat, s=110, facecolor="none", edgecolor="k", lw=1.0, transform=ccrs.PlateCarree(), zorder=6)
        a.text(r.lon + 0.01, r.lat + 0.01, f"{r.matchup_id}", fontsize=8, transform=ccrs.PlateCarree(), zorder=7)
    a.set_title(f"{', '.join(sorted(g.FIELD_experiment.unique()))}\n{len(g)} outliers, "
                f"{len(allp)} matchups in this area", fontsize=9)
    gl = a.gridlines(draw_labels=True, lw=0.2)
    gl.top_labels = gl.right_labels = False
fig.tight_layout()
fig.savefig(OUT / "11_sif_outlier_map.png", dpi=150)
plt.close(fig)

# ------------------------------------------------------------------ spectra
rrs_cols = [c for c in box.columns if c.startswith("Rrs_")]
lam_rrs = np.array([float(c[4:]) for c in rrs_cols])
mr = (lam_rrs >= LAM[0]) & (lam_rrs <= LAM[1])
lam_f = np.array(sorted(float(m.group(1)) for c in orig.columns if (m := re.fullmatch(r"FIELD_Rrs_(\d+)", c))))
mf = (lam_f >= LAM[0]) & (lam_f <= LAM[1])
summ = pd.read_csv(HERE / "matchup_summary.csv").set_index("matchup_id")


def l1b_box(g, b):
    p = f"{K}/L1B_V3/{g[:4]}/{g[4:6]}/{g[6:8]}/PACE_OCI.{g}.L1B.V3.nc"
    s0, s1, p0, p1 = b.scan.min(), b.scan.max() + 1, b.pixel.min(), b.pixel.max() + 1
    with netCDF4.Dataset(p) as d:
        d.set_auto_mask(True)
        sbp = d["sensor_band_parameters"]
        wl = np.ma.filled(sbp["red_wavelength"][:], np.nan).astype(float)
        E0 = np.ma.filled(sbp["red_solar_irradiance"][:], np.nan).astype(float)
        es = float(d.getncattr("earth_sun_distance_correction"))
        k = np.flatnonzero((wl >= LAM[0]) & (wl <= LAM[1]))
        k = k[np.argsort(wl[k])]
        rh = np.ma.filled(d["observation_data"]["rhot_red"][k, s0:s1, p0:p1], np.nan).astype(float)
        sza = np.ma.filled(d["geolocation_data"]["solar_zenith"][s0:s1, p0:p1], np.nan).astype(float)
    L = rh * E0[k][:, None, None] * np.cos(np.radians(sza))[None] / np.pi / es
    return wl[k], {(s0 + a, p0 + c): L[:, a, c] for a in range(s1 - s0) for c in range(p1 - p0)}


rows = list(out.itertuples()) + [ref]
fig, axs = plt.subplots(len(rows), 2, figsize=(14, 2.7 * len(rows)), squeeze=False)
for i, r in enumerate(rows):
    mid = int(r.matchup_id)
    b = box[box.matchup_id == mid]
    pp = pix[pix.matchup_id == mid].set_index(["scan", "pixel"])
    wl, L = l1b_box(str(r.granule), b)
    ax = axs[i, 0]
    for _, q in b[b.valid].iterrows():
        sif_q = pp.loc[(q.scan, q.pixel)]
        c = "k" if q.is_centre else ("0.55" if sif_q.converged else "tab:orange")
        ax.plot(wl, L[(q.scan, q.pixel)], color=c, lw=1.3 if q.is_centre else 0.5, zorder=3 if q.is_centre else 2)
    for w0, w1 in ((686, 697), (757, 771)):
        ax.axvspan(w0, w1, color="0.93", zorder=0)
    ax.axvspan(640, 780, ymin=0, ymax=0.03, color="tab:blue", alpha=0.5)
    ax.axvline(678.2, color="tab:green", lw=0.6, ls=":")
    ax.set_xlim(*LAM)
    ax.set_ylabel("TOA radiance\n[W m$^{-2}$ sr$^{-1}$ µm$^{-1}$]", fontsize=8)
    tag = "REFERENCE: " if i == len(rows) - 1 else ""
    ax.set_title(f"{tag}matchup {mid}  {r.data_type} / {r.FIELD_experiment}  ({r.lat:.3f}, {r.lon:.3f})\n"
                 f"SIF {r.sif_mean:.3f}, L2 nflh {r.nflh_mean:.3f}, in-situ nFLH {r.nflh680_insitu:.3f}; "
                 f"χ²ᵣ {r.chi2:.1f}, α {r.alpha:.2f}, dark {r.dark:.0%}; chl {r.chl:.1f}, AOT865 {r.aot:.3f}, "
                 f"SZA {r.sza:.0f}, VZA {r.vza:.0f}", fontsize=8, loc="left")
    ax = axs[i, 1]
    for _, q in b[b.valid].iterrows():
        ax.plot(lam_rrs[mr], q[rrs_cols].to_numpy(float)[mr], color="k" if q.is_centre else "0.6",
                lw=1.3 if q.is_centre else 0.5, zorder=3 if q.is_centre else 2)
    fr = orig.loc[mid, [f"FIELD_Rrs_{int(w)}" for w in lam_f[mf]]].to_numpy(float)
    ok = np.isfinite(fr)
    ax.plot(lam_f[mf][ok], fr[ok], "o-", color="tab:green", ms=2.5, lw=1.1, zorder=4)
    ax.axvspan(588, 613, color="0.93", zorder=0)
    ax.axvline(678.2, color="tab:green", lw=0.6, ls=":")
    ax.axhline(0, color="0.8", lw=0.5)
    ax.set_xlim(LAM[0], 725)                      # L2 Rrs ends at 719 nm
    ax.set_ylabel("Rrs [sr$^{-1}$]", fontsize=8)
    for a in axs[i]:
        a.tick_params(labelsize=7)
        a.spines[["top", "right"]].set_visible(False)
axs[-1, 0].set_xlabel("Wavelength (nm)")
axs[-1, 1].set_xlabel("Wavelength (nm)")
axs[0, 0].plot([], [], "k", lw=1.3, label="centre pixel")
axs[0, 0].plot([], [], color="0.55", lw=0.5, label="valid pixel, converged")
axs[0, 0].plot([], [], color="tab:orange", lw=0.5, label="valid pixel, not converged")
axs[0, 0].legend(fontsize=7, loc="upper right")
axs[0, 1].plot([], [], "k", lw=1.3, label="PACE L2 Rrs, centre pixel")
axs[0, 1].plot([], [], color="0.6", lw=0.5, label="PACE L2 Rrs, valid pixels")
axs[0, 1].plot([], [], "o-", color="tab:green", ms=2.5, label="in-situ FIELD_Rrs")
axs[0, 1].legend(fontsize=7, loc="upper right")
fig.suptitle("High PACE SIF at low PACE nflh: L1B TOA radiance (left; grey bands = O₂-B and O₂-A, blue bar = "
             "640-780 nm fit window) and L2 Rrs (right; grey = 588-613 nm, no L2 bands); dotted = 678.2 nm",
             fontsize=10, y=1.0)
fig.tight_layout()
fig.savefig(OUT / "12_sif_outlier_spectra.png", dpi=130, bbox_inches="tight")
print("wrote out/11_sif_outlier_map.png, out/12_sif_outlier_spectra.png")
