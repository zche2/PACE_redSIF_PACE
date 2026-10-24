"""nFLH at 680 nm from the in-situ Rrs of the matchup CSV, compared with PACE.

    python TROPOMI_OCI_colocation/field_rrs_matchups/insitu_nflh680.py [fullmatchups_Rrs_GUI.csv] [V3.2 AOP file]

In situ (FIELD_Rrs, on the CSV's 5 nm field grid):
    nLw = Rrs * F0, F0 = OCI band-mean solar flux (AOP sensor_band_parameters) interpolated to the field λ
    baseline = least-squares line through nLw at all finite bands in SHORT (664-667.5 nm) and LONG
               (710-713.5 nm), evaluated at the height band
    nFLH = nLw(height band) - baseline;  height band = 680 nm, else the nearest finite band within 5 nm
On the 5 nm grid these windows hold one band each (665 and 710 nm). Multispectral in-situ spectra (< 40
finite bands: aircraft and some casts) have no 680 nm band; they use 685 nm and are marked in the plots.
The uncertainty propagates FIELD_Rrs_unc through the same linear combination (bands independent).
In-situ nFLH is not clipped.

PACE, for comparison:
    nflh_L2         stored V3.2 L2 nflh (680 nm definition, clipped at 0), mean of the valid box pixels
    nflh_pace_grid  the in-situ formula applied to PACE_Rrs_rep (box Rrs smoothed to the field grid), same
                    bands as the in-situ spectrum, unclipped

Outputs: insitu_nflh680.csv, out/9_insitu_nflh680.png
"""
import re
import sys
from pathlib import Path

import matplotlib
import netCDF4
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path("/home/zhe2/FraLab/PACE_redSIF_PACE/TROPOMI_OCI_colocation/fullmatchups_Rrs_GUI.csv")
AOP = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(
    "/kiwi-data/Data/satellite/PACE_OCI/L2_AOP_V3.2/2025/07/08/PACE_OCI.20250708T204612.L2.OC_AOP.V3_2.nc")
SHORT, LONG = (664.0, 667.5), (710.0, 713.5)
HEIGHT, HEIGHT_TOL = 680.0, 5.0
HYPER_MIN = 40                                   # finite FIELD_Rrs bands for "hyperspectral"
UNC_MAX_PLOT = 1.0                               # error bars above this are not drawn (aircraft: median ~13)

with netCDF4.Dataset(AOP) as d:
    sb = d["sensor_band_parameters"]
    fw, F0 = np.ma.filled(sb["fwave"][:]).astype(float), np.ma.filled(sb["F0"][:]).astype(float)
o = np.argsort(fw)
fw, F0 = fw[o], F0[o]

orig = pd.read_csv(CSV, low_memory=False)
summ = pd.read_csv(HERE / "matchup_summary.csv").set_index("matchup_id")
lam = np.array(sorted(float(m.group(1)) for c in orig.columns if (m := re.fullmatch(r"FIELD_Rrs_(\d+)", c))))
F0f = np.interp(lam, fw, F0)
R = orig.loc[summ.index, [f"FIELD_Rrs_{int(w)}" for w in lam]].to_numpy(float)
U = orig.loc[summ.index, [f"FIELD_Rrs_unc_{int(w)}" for w in lam]].to_numpy(float)
rep_lam = np.array([w for w in lam if f"PACE_Rrs_rep_{int(w)}" in summ.columns])


def flh(r, u=None, lam=lam, F=F0f):
    """nFLH, its uncertainty, and the bands used (short, long, height) for one spectrum on grid lam."""
    ok = np.isfinite(r)
    sb_ = np.flatnonzero(ok & (lam >= SHORT[0]) & (lam <= SHORT[1]))
    lb_ = np.flatnonzero(ok & (lam >= LONG[0]) & (lam <= LONG[1]))
    hc = np.flatnonzero(ok & (np.abs(lam - HEIGHT) <= HEIGHT_TOL))
    if not (len(sb_) and len(lb_) and len(hc)):
        return np.nan, np.nan, (np.nan, np.nan, np.nan)
    hi = hc[np.lexsort((-lam[hc], np.abs(lam[hc] - HEIGHT)))[0]]      # nearest; a tie goes to the longer band
    bi = np.r_[sb_, lb_]
    A = np.column_stack([np.ones(len(bi)), lam[bi]])
    w = np.array([1.0, lam[hi]]) @ np.linalg.pinv(A)
    nLw = r * F
    val = nLw[hi] - nLw[bi] @ w
    unc = np.nan if u is None else np.sqrt((F[hi] * u[hi]) ** 2 + np.sum((w * F[bi] * u[bi]) ** 2))
    return val, unc, (lam[sb_].mean(), lam[lb_].mean(), lam[hi])


rows = []
for k, mid in enumerate(summ.index):
    s = summ.loc[mid]
    nb = int(np.isfinite(R[k]).sum())
    v, u, (bs, bl, bh) = flh(R[k], U[k])
    # PACE with the same bands, from the box Rrs smoothed to the field grid
    pr = np.array([s.get(f"PACE_Rrs_rep_{int(w)}", np.nan) for w in lam], float)
    pr[~np.isin(lam, [bs, bl, bh])] = np.nan
    vp = flh(pr)[0] if np.isfinite(bh) else np.nan
    rows.append({"matchup_id": mid, "data_type": s.data_type, "FIELD_experiment": s.FIELD_experiment,
                 "n_field_bands": nb, "hyperspectral": nb >= HYPER_MIN,
                 "short_band": bs, "long_band": bl, "height_band": bh,
                 "nflh680_insitu": v, "nflh680_insitu_unc": u,
                 "nflh_L2": s.nflh_mean_valid, "nflh_pace_grid": vp, "passes_protocol": s.passes_protocol})
df = pd.DataFrame(rows)
df.to_csv(HERE / "insitu_nflh680.csv", index=False)

use = df.passes_protocol & df.nflh680_insitu.notna()
print(f"{len(df)} matchups; in-situ nFLH computed for {df.nflh680_insitu.notna().sum()}; "
      f"used (passes protocol): {use.sum()}")
print("bands used (short, long, height) x count:")
print(df[df.nflh680_insitu.notna()].groupby(["hyperspectral", "short_band", "long_band", "height_band"])
      .size().to_string())
print("no in-situ nFLH (missing band):", df.loc[df.nflh680_insitu.isna(), "data_type"].value_counts().to_dict())
for h, g in df[use].groupby("hyperspectral"):
    for ref in ("nflh_L2", "nflh_pace_grid"):
        x, y = g[ref], g.nflh680_insitu
        print(f"  {'hyper' if h else 'multi'}spectral (n={len(g)}) vs {ref:14s}: r = {np.corrcoef(x, y)[0, 1]:.3f}, "
              f"median (in situ − PACE) {np.median(y - x):+.4f}, median ratio {np.median(y / x):.3f}")

# ------------------------------------------------------------------ figure
types = sorted(df.data_type.unique())
col = dict(zip(types, plt.cm.tab10.colors))
fig, axs = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw=dict(width_ratios=[1, 1, 1.15]))
for ax, ref, lab in ((axs[0], "nflh_L2", "PACE L2 nflh (V3.2, clipped at 0), mean of valid pixels"),
                     (axs[1], "nflh_pace_grid", "PACE, in-situ formula on the field grid (unclipped)")):
    g = df[use]
    for t in types:
        for h, mk in ((True, "o"), (False, "^")):
            m = (g.data_type == t) & (g.hyperspectral == h)
            if m.any():
                ue = g.nflh680_insitu_unc[m].where(g.nflh680_insitu_unc[m] <= UNC_MAX_PLOT)
                ax.errorbar(g[ref][m], g.nflh680_insitu[m], yerr=ue, fmt=mk, ms=5,
                            mfc=col[t] if h else "none", mec=col[t], ecolor=col[t], elinewidth=0.5, alpha=0.8,
                            label=f"{t}{'' if h else ' (multispectral, height 685)'}")
    vals = np.r_[g[ref], g.nflh680_insitu]
    lim = (min(np.nanmin(vals), 0) - 0.05, np.nanmax(vals) * 1.08)
    ax.plot(lim, lim, "k--", lw=0.8)
    ax.axhline(0, color="0.8", lw=0.5)
    ax.axvline(0, color="0.8", lw=0.5)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(lab + "\n[W m$^{-2}$ µm$^{-1}$ sr$^{-1}$]", fontsize=9)
    ax.set_ylabel("in-situ nFLH$_{680}$ (± FIELD_Rrs_unc)")
    ax.spines[["top", "right"]].set_visible(False)
axs[0].legend(fontsize=7, loc="upper left")
nbig = (df[use].nflh680_insitu_unc > UNC_MAX_PLOT).sum()
axs[0].set_title(f"{use.sum()} matchups passing the protocol;  filled o: hyperspectral (665 | 710, height 680),\n"
                 f"open △: multispectral (665 | 710, height 685);  {nbig} error bars > {UNC_MAX_PLOT:g} not drawn",
                 fontsize=9)
axs[1].set_title("same bands for PACE and in situ\n(separates the band definition from the measurement)", fontsize=9)

# example spectra: nLw in the red with the bands and baseline
ax = axs[2]
ex = [df[use & df.hyperspectral].sort_values("nflh680_insitu").iloc[len(df[use & df.hyperspectral]) * 3 // 4],
      df[use & ~df.hyperspectral].sort_values("nflh680_insitu").iloc[len(df[use & ~df.hyperspectral]) * 3 // 4]]
for e, c in zip(ex, ("tab:blue", "tab:orange")):
    k = summ.index.get_loc(e.matchup_id)
    f_ = np.isfinite(R[k]) & (lam >= 640) & (lam <= 720)
    nl = R[k] * F0f
    ax.plot(lam[f_], nl[f_], "-o", color=c, ms=3, lw=1,
            label=f"matchup {e.matchup_id} ({e.data_type}, {'hyper' if e.hyperspectral else 'multi'}spectral, "
                  f"{e.n_field_bands} bands)")
    bands = [e.short_band, e.long_band]
    yb = [nl[lam == b][0] for b in bands]
    sl = (yb[1] - yb[0]) / (bands[1] - bands[0])
    xx = np.array([660, 715])
    ax.plot(xx, yb[0] + sl * (xx - bands[0]), "--", color=c, lw=0.9)
    ax.plot(bands, yb, "s", color=c, ms=7, mfc="none")
    yh = nl[lam == e.height_band][0]
    ax.annotate("", xy=(e.height_band, yh), xytext=(e.height_band, yb[0] + sl * (e.height_band - bands[0])),
                arrowprops=dict(arrowstyle="<->", color=c, lw=1.2))
    ax.text(e.height_band + 1.5, (yh + yb[0] + sl * (e.height_band - bands[0])) / 2,
            f"nFLH {e.nflh680_insitu:.3f}\n(height {e.height_band:.0f})", color=c, fontsize=8, va="center")
for w0, w1 in (SHORT, LONG):
    ax.axvspan(w0, w1, color="0.9", zorder=0)
ax.axvline(HEIGHT, color="0.7", lw=0.6, ls=":")
ax.set_xlim(640, 720)
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("in-situ nLw = Rrs·F0  [W m$^{-2}$ µm$^{-1}$ sr$^{-1}$]")
ax.set_title("examples: grey = baseline windows, squares = bands used,\ndashed = baseline, arrow = nFLH",
             fontsize=9)
ax.legend(fontsize=7, loc="upper right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "9_insitu_nflh680.png", dpi=150)
print("wrote insitu_nflh680.csv, out/9_insitu_nflh680.png")
