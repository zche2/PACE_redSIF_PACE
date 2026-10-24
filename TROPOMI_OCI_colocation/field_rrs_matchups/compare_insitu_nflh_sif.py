"""In-situ nFLH (680 nm) against PACE L2 nflh and PACE SVD SIF at the matchups.

    python TROPOMI_OCI_colocation/field_rrs_matchups/compare_insitu_nflh_sif.py [matchup_sif_<config>.csv]

Inputs (this folder): insitu_nflh680.csv (insitu_nflh680.py), matchup_sif_nPC15_npoly3_O2A.csv
(retrieve_matchup_sif.jl; per box pixel, carries the L2 nflh too), matchup_summary.csv.
PACE per matchup = mean ± std over the box pixels:
    nflh  L2 V3.2 (680 nm, clipped at 0) over valid pixels (L2 flags + geometry)
    SIF   water-leaving SIF at 678.2 nm over valid pixels whose retrieval converged
Matchups used: pass the protocol (|dt| <= 3 h, >= 13 valid pixels, CV <= 0.15), in-situ nFLH available,
>= MIN_CONV converged valid pixels; EXCLUDE lists in-situ spectra rejected by eye.
Markers: filled o = hyperspectral in situ (665 | 710, height 680); open △ = multispectral (height 685).

Output: out/10_insitu_vs_pace_nflh_sif.png, insitu_vs_pace.csv
"""
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
DARK = "--dark" in sys.argv                      # PACE over dark pixels only (max Rtoa 641-779 nm <= 20)
SIF_CSV = Path(ARGS[0]) if ARGS else HERE / "matchup_sif_nPC15_npoly3_O2A.csv"
TAG = SIF_CSV.stem.replace("matchup_sif_", "")
SFX = ("" if TAG == "nPC15_npoly3_O2A" else f"_{TAG}") + ("_dark" if DARK else "")   # default keeps old names
MIN_CONV = 5
EXCLUDE = {10: "in-situ Rrs(710) 7x Rrs(705): broken NIR in the cast"}
UNC_MAX_PLOT = 1.0

ins = pd.read_csv(HERE / "insitu_nflh680.csv")
pix = pd.read_csv(SIF_CSV)


def box_stats(v):
    vc = v[v.converged]
    gv, gc = v.groupby("matchup_id"), vc.groupby("matchup_id")
    return pd.DataFrame({"nflh_mean": gv.nflh.mean(), "nflh_std": gv.nflh.std(),
                         "sif_mean": gc.sif_678.mean(), "sif_std": gc.sif_678.std(),
                         "sif_sigma_mean": gc.sif_678_sigma.mean(), "n_valid": gv.size(), "n_conv": gc.size()})


pace_all = box_stats(pix[pix.valid])
pace = box_stats(pix[pix.valid & pix.is_dark]) if DARK else pace_all
df = ins.merge(pace.reset_index(), on="matchup_id", how="left")
df = df.merge(pace_all.add_suffix("_allpix").reset_index(), on="matchup_id", how="left")
base_ok = df.passes_protocol & df.nflh680_insitu.notna() & ~df.matchup_id.isin(EXCLUDE)
df["used"] = base_ok & (df.n_conv.fillna(0) >= MIN_CONV)
# --dark: matchups without enough dark pixels, plotted at their all-pixel values and not fitted
df["marked_not_dark"] = (base_ok & ~df.used & (df.n_conv_allpix >= MIN_CONV)) if DARK else False
df.to_csv(HERE / f"insitu_vs_pace{SFX}.csv", index=False)
u = df[df.used]
nd = df[df.marked_not_dark]
print(f"{len(df)} matchups; used {len(u)} ({u.hyperspectral.sum()} hyperspectral, {(~u.hyperspectral).sum()} "
      f"multispectral){'; PACE over dark pixels only, ' + str(len(nd)) + ' matchups with < ' + str(MIN_CONV) + ' dark converged pixels marked' if DARK else ''}; "
      f"excluded by hand: {EXCLUDE}")
if DARK:
    print("  not-dark matchups by type:", nd.data_type.value_counts().to_dict())


def stats(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    s, b = np.polyfit(x, y, 1)
    return dict(n=m.sum(), r=np.corrcoef(x, y)[0, 1], slope=s, icpt=b, bias=np.median(y - x),
                ratio=np.median(y[x > 0.05] / x[x > 0.05]))


groups = {"hyperspectral": u[u.hyperspectral], "multispectral": u[~u.hyperspectral], "all": u}
for yv, nm in (("nflh_mean", "PACE L2 nflh"), ("sif_mean", "PACE SIF678")):
    for g, d in groups.items():
        s = stats(d.nflh680_insitu, d[yv])
        print(f"  {nm:13s} vs in-situ nFLH, {g:13s}: n={s['n']:3d}  r={s['r']:.3f}  y = {s['slope']:.3f}x "
              f"{s['icpt']:+.3f}  median(y-x) {s['bias']:+.3f}  median y/x {s['ratio']:.2f}")
s = stats(u.nflh_mean, u.sif_mean)
print(f"  PACE SIF678 vs PACE nflh (all): r={s['r']:.3f}  SIF = {s['slope']:.3f} nflh {s['icpt']:+.3f}")

# ------------------------------------------------------------------ figure
types = sorted(df.data_type.unique())
col = dict(zip(types, plt.cm.tab10.colors))
fig, axs = plt.subplots(1, 3, figsize=(18, 5.8))


def scatter(ax, x, y, xe=None, ye=None, d=u):
    for t in types:
        for h, mk in ((True, "o"), (False, "^")):
            m = (d.data_type == t) & (d.hyperspectral == h)
            if not m.any():
                continue
            xerr = None if xe is None else d[xe][m].where(d[xe][m] <= UNC_MAX_PLOT)
            yerr = None if ye is None else d[ye][m]
            ax.errorbar(d[x][m], d[y][m], xerr=xerr, yerr=yerr, fmt=mk, ms=5, mfc=col[t] if h else "none",
                        mec=col[t], ecolor=col[t], elinewidth=0.5, alpha=0.8,
                        label=f"{t}{'' if h else ' (multispectral)'}")


def mark_not_dark(ax, x, y):
    """--dark: matchups without enough dark pixels, at their all-pixel PACE values (not in the fits)."""
    if not len(nd):
        return
    xs = nd[x + "_allpix"] if x + "_allpix" in nd else nd[x]
    ys = nd[y + "_allpix"] if y + "_allpix" in nd else nd[y]
    ax.scatter(xs, ys, marker="x", s=22, color="0.55", lw=0.9, zorder=1,
               label=f"not dark (< {MIN_CONV} dark pixels; all-pixel mean, not fitted): {len(nd)}")


def fitline(ax, x, y, d, color, label):
    s = stats(d[x], d[y])
    xx = np.array(ax.get_xlim())
    ax.plot(xx, s["slope"] * xx + s["icpt"], color=color, lw=1.2,
            label=f"{label}: r={s['r']:.2f}, y={s['slope']:.2f}x{s['icpt']:+.3f} (n={s['n']})")


def lims(ax, *cols_, d=u):
    vals = np.concatenate([d[c].to_numpy(float) for c in cols_] +
                          [nd[c + "_allpix" if c + "_allpix" in nd else c].to_numpy(float) for c in cols_])
    lo, hi = min(np.nanmin(vals), 0) - 0.03, np.nanmax(vals) * 1.08
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
    ax.axhline(0, color="0.85", lw=0.5)
    ax.axvline(0, color="0.85", lw=0.5)


unit = "[W m$^{-2}$ µm$^{-1}$ sr$^{-1}$]"
ax = axs[0]
mark_not_dark(ax, "nflh680_insitu", "nflh_mean")
scatter(ax, "nflh680_insitu", "nflh_mean", "nflh680_insitu_unc", "nflh_std")
lims(ax, "nflh680_insitu", "nflh_mean")
fitline(ax, "nflh680_insitu", "nflh_mean", u[u.hyperspectral], "k", "hyperspectral")
fitline(ax, "nflh680_insitu", "nflh_mean", u[~u.hyperspectral], "0.5", "multispectral")
ax.set_ylabel(f"PACE L2 nflh (V3.2, 680 nm), box mean ± std {unit}")
ax.set_title("PACE L2 nflh vs in situ", fontsize=10)

ax = axs[1]
mark_not_dark(ax, "nflh680_insitu", "sif_mean")
scatter(ax, "nflh680_insitu", "sif_mean", "nflh680_insitu_unc", "sif_std")
lims(ax, "nflh680_insitu", "sif_mean")
fitline(ax, "nflh680_insitu", "sif_mean", u[u.hyperspectral], "k", "hyperspectral")
fitline(ax, "nflh680_insitu", "sif_mean", u[~u.hyperspectral], "0.5", "multispectral")
ax.set_ylabel(f"PACE SIF 678 nm (SVD, {SIF_CSV.stem.replace('matchup_sif_', '')}), box mean ± std {unit}",
              fontsize=9)
ax.set_title("PACE SIF vs in-situ nFLH", fontsize=10)

for ax in axs[:2]:
    ax.set_xlabel(f"in-situ nFLH$_{{680}}$ (665 | 710 nm baseline) {unit}")
    h, lab = ax.get_legend_handles_labels()
    fits = [i for i, t in enumerate(lab) if t.startswith(("hyper", "multi", "not dark"))]
    ax.legend([h[i] for i in fits], [lab[i] for i in fits], fontsize=7.5, loc="lower right")

ax = axs[2]
mark_not_dark(ax, "nflh_mean", "sif_mean")
scatter(ax, "nflh_mean", "sif_mean", "nflh_std", "sif_std")
lims(ax, "nflh_mean", "sif_mean")
fitline(ax, "nflh_mean", "sif_mean", u, "k", "all")
ax.set_xlabel(f"PACE L2 nflh, box mean {unit}")
ax.set_ylabel(f"PACE SIF 678 nm, box mean {unit}")
ax.set_title("PACE SIF vs PACE nflh (same matchups)", fontsize=10)
h, lab = ax.get_legend_handles_labels()
ax.legend(h, lab, fontsize=7, loc="upper left")
for ax in axs:
    ax.spines[["top", "right"]].set_visible(False)
pix_txt = "DARK valid pixels only (max TOA radiance 641-779 nm ≤ 20)" if DARK else "valid pixels"
fig.suptitle(f"PACE box means over {pix_txt}.  {len(u)} matchups fitted (protocol passed, ≥{MIN_CONV} converged "
             f"such pixels; matchup {', '.join(map(str, EXCLUDE))} excluded)"
             f"{f'; grey ×: {len(nd)} not-dark matchups' if DARK else ''}\n"
             f"filled o: hyperspectral in situ (height 680 nm), open △: multispectral (height 685 nm); "
             f"in-situ error bars > {UNC_MAX_PLOT:g} not drawn", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / f"10_insitu_vs_pace_nflh_sif{SFX}.png", dpi=150)
print(f"wrote out/10_insitu_vs_pace_nflh_sif{SFX}.png, insitu_vs_pace{SFX}.csv")
