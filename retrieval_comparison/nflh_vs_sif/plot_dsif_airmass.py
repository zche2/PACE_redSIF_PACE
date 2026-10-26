"""ΔSIF = SIF − nFLH against geometric air mass 1/cos(SZA) + 1/cos(VZA), grouped by chlorophyll, for single days.

Day / per-sounding version of notebooks/figures/Figure2/FigureS2_low_chl_vza_bin_month.py (air-mass figure),
on the soundings from extract_soundings.py with the shared --l2 / --query filter (common.py).

    python retrieval_comparison/nflh_vs_sif/plot_dsif_airmass.py --l2 LAND,COASTZ,HISOLZEN \\
        --query "is_dark == 1 and status_code in [1, 5] and abs(sif_678) <= 5.0 and reduced_chi2 <= 1.5" \\
        --tag 3flags_dark_status15_chi2le1.5

Soundings also need finite, positive chlor_a inside the class edges (--chl-classes). Per date (one row):
(a) median nFLH, (b) median SIF, (c) median and IQR of ΔSIF against air mass, one line per Chl-a class;
(d) median ΔSIF by Chl-a class × air mass. Cells with fewer than --min-per-cell soundings are left blank.

Output: out/<tag>/dsif_airmass_chl.png, out/<tag>/dsif_airmass_chl.csv (per class × air-mass bin quantiles)
"""
import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import UNIT, add_filter_args, filter_text, load_filtered, out_dir  # noqa: E402

p = argparse.ArgumentParser()
add_filter_args(p)
p.add_argument("--chl-classes", type=float, nargs="+", default=[0.02, 0.04, 0.07, 0.12, 0.2, 0.4],
               help="Chl-a class edges (mg m-3), as in FigureS2")
p.add_argument("--amf-edges", type=float, nargs=3, default=[2.0, 6.0, 0.1], metavar=("MIN", "MAX", "STEP"))
p.add_argument("--min-per-cell", type=int, default=500)
A = p.parse_args()
OUT = out_dir(HERE, A)

CHL = np.array(A.chl_classes)
NC = len(CHL) - 1
AMF = np.round(np.arange(A.amf_edges[0], A.amf_edges[1] + 1e-9, A.amf_edges[2]), 4)
XM = 0.5 * (AMF[:-1] + AMF[1:])
LABELS = [f"{CHL[i]:g}–{CHL[i + 1]:g}" for i in range(NC)]
COLORS = plt.get_cmap("viridis")(np.linspace(0.1, 0.9, NC))
COLS = ["granule", "l2_flags", "sif_678", "nflh", "chlor_a", "sza", "vza", "converged", "status_code",
        "reduced_chi2", "is_dark", "is_ocean", "sif_678_sigma", "alpha", "aot_865", "lat", "lon", "rtoa_678", "dof",
        "n_steps"]


def binned(df):
    """Quantiles (25/50/75) of nFLH, SIF and ΔSIF per Chl-a class × air-mass bin, as arrays [class, bin, q]."""
    amf = 1.0 / np.cos(np.radians(df.sza)) + 1.0 / np.cos(np.radians(df.vza))
    kc = np.digitize(df.chlor_a, CHL) - 1
    kx = np.digitize(amf, AMF) - 1
    ok = (kc >= 0) & (kc < NC) & (kx >= 0) & (kx < len(XM)) & np.isfinite(amf)
    d = pd.DataFrame({"kc": kc[ok], "kx": kx[ok], "nflh": df.nflh.to_numpy()[ok], "sif": df.sif_678.to_numpy()[ok]})
    d["diff"] = d.sif - d.nflh
    g = d.groupby(["kc", "kx"])
    n = g.size()
    q = g[["nflh", "sif", "diff"]].quantile([0.25, 0.5, 0.75]).unstack()
    res = {k: np.full((NC, len(XM), 3), np.nan) for k in ("nflh", "sif", "diff")}
    cnt = np.zeros((NC, len(XM)), int)
    for (ic, ix), m in n.items():
        cnt[ic, ix] = m
        if m >= A.min_per_cell:
            for k in res:
                res[k][ic, ix] = q.loc[(ic, ix), k].to_numpy()
    return res, cnt, int(ok.sum())


fig, axs = plt.subplots(len(A.dates), 4, figsize=(24, 5.2 * len(A.dates)), squeeze=False,
                        gridspec_kw=dict(width_ratios=[1, 1, 1, 1.15]))
rows = []
for i, date in enumerate(A.dates):
    df, steps = load_filtered(A, date, columns=COLS)
    df = df[np.isfinite(df.chlor_a) & (df.chlor_a > 0)]
    res, cnt, n_used = binned(df)
    for ic in range(NC):
        for ix in range(len(XM)):
            rows.append(dict(date=date, chl_lo=CHL[ic], chl_hi=CHL[ic + 1], amf_lo=AMF[ix], amf_hi=AMF[ix + 1],
                             n=cnt[ic, ix], **{f"{k}_q{qq}": res[k][ic, ix, j] for k in res
                                               for j, qq in enumerate((25, 50, 75))}))
    ax_f, ax_s, ax_d, ax_h = axs[i]
    for ic in range(NC):
        ax_f.plot(XM, res["nflh"][ic, :, 1], "o-", ms=3, lw=1.4, color=COLORS[ic], label=LABELS[ic])
        ax_s.plot(XM, res["sif"][ic, :, 1], "o-", ms=3, lw=1.4, color=COLORS[ic])
        ax_d.fill_between(XM, res["diff"][ic, :, 0], res["diff"][ic, :, 2], color=COLORS[ic], alpha=0.12, lw=0)
        ax_d.plot(XM, res["diff"][ic, :, 1], "o-", ms=3, lw=1.4, color=COLORS[ic])
    ax_s.sharey(ax_f)
    ax_f.set_title(f"{date}  (a) L2 nFLH (median) vs air mass", loc="left", fontsize=10)
    ax_s.set_title("(b) SIF 678 nm (median) vs air mass", loc="left", fontsize=10)
    ax_d.set_title("(c) ΔSIF = SIF − nFLH (median and IQR)", loc="left", fontsize=10)
    ax_f.set_ylabel(f"radiance [{UNIT}]")
    ax_d.set_ylabel(f"ΔSIF [{UNIT}]")
    ax_d.axhline(0, color="0.6", lw=0.8)
    ax_f.legend(title=r"Chl-a (mg m$^{-3}$)", frameon=False, fontsize=8, title_fontsize=8, loc="upper center", ncol=3)
    occ = np.flatnonzero(np.isfinite(res["diff"][..., 1]).any(axis=0))
    med = np.ma.masked_invalid(res["diff"][..., 1])
    vmax = np.nanpercentile(np.abs(res["diff"][..., 1]), 98) if occ.size else 1.0
    pc = ax_h.pcolormesh(AMF, np.arange(NC + 1), med, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax_h.set_yticks(np.arange(NC) + 0.5, LABELS)
    ax_h.set_ylabel(r"Chl-a class (mg m$^{-3}$)")
    ax_h.set_title(f"(d) median ΔSIF by Chl-a × air mass (blank: < {A.min_per_cell} soundings)", loc="left",
                   fontsize=10)
    fig.colorbar(pc, ax=ax_h, label=f"ΔSIF [{UNIT}]", extend="both")
    for ax in (ax_f, ax_s, ax_d, ax_h):
        ax.set_xlabel(r"geometric air mass $1/\cos\theta_s + 1/\cos\theta_v$")
        ax.spines[["top", "right"]].set_visible(False)
        if occ.size:
            ax.set_xlim(AMF[occ[0]], AMF[occ[-1] + 1])
    ax_s.text(0.98, 0.02, f"{n_used:,} soundings in the Chl-a classes", transform=ax_s.transAxes, ha="right",
              fontsize=8, color="0.3")
    print(f"{date}: {steps[-1][1]:,} after filter, {n_used:,} with Chl-a in {CHL[0]:g}–{CHL[-1]:g} mg m-3")
fig.suptitle(filter_text(A) + "   (soundings with Chl-a inside the class edges)", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "dsif_airmass_chl.png", dpi=150)
pd.DataFrame(rows).to_csv(OUT / "dsif_airmass_chl.csv", index=False)
print(f"wrote {OUT}/dsif_airmass_chl.png, dsif_airmass_chl.csv")
