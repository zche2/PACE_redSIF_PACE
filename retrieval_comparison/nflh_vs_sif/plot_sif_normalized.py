"""Actual vs normalized SIF against L2 nFLH, for single days.

nFLH is normalized (Rrs·F0: Sun at zenith, 1 AU); the retrieved SIF is the actual water-leaving radiance.
Normalized SIF here: SIF_norm = SIF · esd / cos(SZA), esd = (d / 1 AU)² from the L1B attribute
earth_sun_distance_correction of each granule (read on the fly; the downwelling transmittance is ignored).

    python retrieval_comparison/nflh_vs_sif/plot_sif_normalized.py --l2 LAND,COASTZ,HISOLZEN \\
        --query "is_dark == 1 and status_code in [1, 5] and abs(sif_678) <= 5.0 and reduced_chi2 <= 1.5" \\
        --tag 3flags_dark_status15_chi2le1.5

One row per date: density of actual SIF vs nFLH, density of SIF_norm vs nFLH, and median ΔSIF against
air mass by Chl-a class for actual and normalized SIF (same y axis).
Output: out/<tag>/sif_normalized.png (no data files)
"""
import argparse
import sys
from pathlib import Path

import matplotlib
import netCDF4
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import UNIT, add_filter_args, filter_text, load_filtered, out_dir  # noqa: E402

L1B = "/kiwi-data/Data/satellite/PACE_OCI/L1B_V3"
p = argparse.ArgumentParser()
add_filter_args(p)
p.add_argument("--chl-classes", type=float, nargs="+", default=[0.02, 0.04, 0.07, 0.12, 0.2, 0.4])
p.add_argument("--min-per-cell", type=int, default=500)
A = p.parse_args()
OUT = out_dir(HERE, A)
CHL = np.array(A.chl_classes)
NC = len(CHL) - 1
AMF = np.round(np.arange(2.0, 6.0 + 1e-9, 0.1), 2)
XM = 0.5 * (AMF[:-1] + AMF[1:])
COLORS = plt.get_cmap("viridis")(np.linspace(0.1, 0.9, NC))
LABELS = [f"{CHL[i]:g}–{CHL[i + 1]:g}" for i in range(NC)]
COLS = ["granule", "l2_flags", "sif_678", "nflh", "chlor_a", "sza", "vza", "status_code", "reduced_chi2",
        "is_dark", "converged", "lat"]


def esd_of(g):
    with netCDF4.Dataset(f"{L1B}/{g[:4]}/{g[4:6]}/{g[6:8]}/PACE_OCI.{g}.L1B.V3.nc") as d:
        return float(d.getncattr("earth_sun_distance_correction"))


def stats(x, y):
    s, b = np.polyfit(x, y, 1)
    big = x > 0.05
    return dict(r=np.corrcoef(x, y)[0, 1], slope=s, icpt=b, med=np.median(y - x), ratio=np.median(y[big] / x[big]))


def amf_medians(df, col):
    amf = 1 / np.cos(np.radians(df.sza)) + 1 / np.cos(np.radians(df.vza))
    kc = np.digitize(df.chlor_a, CHL) - 1
    kx = np.digitize(amf, AMF) - 1
    ok = (kc >= 0) & (kc < NC) & (kx >= 0) & (kx < len(XM))
    d = pd.DataFrame({"kc": kc[ok], "kx": kx[ok], "v": (df[col] - df.nflh).to_numpy()[ok]})
    g = d.groupby(["kc", "kx"]).v
    med, n = g.median(), g.size()
    out = np.full((NC, len(XM)), np.nan)
    for (ic, ix), m in med.items():
        if n[(ic, ix)] >= A.min_per_cell:
            out[ic, ix] = m
    return out


fig, axs = plt.subplots(len(A.dates), 4, figsize=(24, 5.3 * len(A.dates)), squeeze=False)
for i, date in enumerate(A.dates):
    df, _ = load_filtered(A, date, columns=COLS)
    esd = {g: esd_of(g) for g in df.granule.astype(str).unique()}
    e = df.granule.astype(str).map(esd).to_numpy(float)
    df["sif_norm"] = df.sif_678.to_numpy() * e / np.cos(np.radians(df.sza.to_numpy()))
    x = df.nflh.to_numpy(float)
    lim = (np.percentile(x, 0.1) - 0.02, np.percentile(np.r_[x, df.sif_norm], 99.9))
    for j, (col, name) in enumerate((("sif_678", "actual SIF"), ("sif_norm", "normalized SIF · esd / cos(SZA)"))):
        y = df[col].to_numpy(float)
        s = stats(x, y)
        ax = axs[i, j]
        h = ax.hist2d(x, y, bins=200, range=[lim, lim], norm=mcolors.LogNorm(), cmap="viridis", rasterized=True)
        fig.colorbar(h[3], ax=ax, label="soundings")
        ax.plot(lim, lim, "k--", lw=0.8)
        xx = np.array(lim)
        ax.plot(xx, s["slope"] * xx + s["icpt"], "r", lw=1, label=f"OLS: y = {s['slope']:.2f} nFLH {s['icpt']:+.3f}")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel(f"L2 nFLH [{UNIT}]")
        ax.set_ylabel(f"{name} [{UNIT}]")
        ax.set_title(f"{date}  {name}\nr = {s['r']:.3f}, median ΔSIF {s['med']:+.4f}, median SIF/nFLH "
                     f"(nFLH > 0.05) {s['ratio']:.2f}", fontsize=9)
        ax.legend(fontsize=7, loc="upper left")
        print(f"{date} {col:9s}: r={s['r']:.3f} slope={s['slope']:.3f} icpt={s['icpt']:+.4f} "
              f"median ΔSIF={s['med']:+.4f} median ratio={s['ratio']:.3f}")
    meds = {c: amf_medians(df, c) for c in ("sif_678", "sif_norm")}
    allv = np.concatenate([m[np.isfinite(m)] for m in meds.values()])
    yl = (np.nanmin(allv) - 0.01, np.nanmax(allv) + 0.01) if allv.size else (-0.1, 0.1)
    for j, (c, name) in enumerate((("sif_678", "actual SIF"), ("sif_norm", "normalized SIF")), start=2):
        ax = axs[i, j]
        for ic in range(NC):
            ax.plot(XM, meds[c][ic], "o-", ms=3, lw=1.4, color=COLORS[ic], label=LABELS[ic])
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_ylim(yl)
        ax.set_xlabel(r"geometric air mass $1/\cos\theta_s + 1/\cos\theta_v$")
        ax.set_ylabel(f"median ΔSIF = {'SIF' if c == 'sif_678' else 'SIF_norm'} − nFLH [{UNIT}]")
        ax.set_title(f"{date}  {name} − nFLH vs air mass, by Chl-a class", fontsize=9)
    axs[i, 2].legend(title=r"Chl-a (mg m$^{-3}$)", fontsize=7, title_fontsize=7, ncol=2)
    for ax in axs[i]:
        ax.spines[["top", "right"]].set_visible(False)
    print(f"   esd: {min(esd.values()):.4f}-{max(esd.values()):.4f}; {len(df):,} soundings")
fig.suptitle(filter_text(A) + "   |   SIF_norm = SIF · esd / cos(SZA), esd = (d/AU)² from L1B", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "sif_normalized.png", dpi=140)
print(f"wrote {OUT}/sif_normalized.png")
