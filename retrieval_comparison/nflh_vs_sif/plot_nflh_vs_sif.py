"""Compare L2 nFLH with SVD SIF for every sounding of a day, with a configurable filter.

    python retrieval_comparison/nflh_vs_sif/plot_nflh_vs_sif.py                       # defaults below
    python ... --query "converged == 1 and is_dark == 1 and reduced_chi2 < 5" --tag dark_chi2lt5
    python ... --l2 LAND,COASTZ,HISOLZEN

Filter = L2 flag cut AND pandas query, applied to the soundings from extract_soundings.py:
    --l2 zero          l2_flags == 0 (no L2 flag set at all; default)
    --l2 any           no L2 flag cut
    --l2 NAME,NAME     these flags must be off (names as in l2_flags:flag_meanings)
    --query "..."      any pandas query on the columns (default "converged == 1"); "" for none
    columns: lat lon sza vza sif_678 sif_678_sigma converged status_code n_steps reduced_chi2 dof is_dark
             is_ocean alpha l2_flags nflh chlor_a aot_865 rtoa_678 granule scan pixel
nFLH is the stored L2 nflh (680 nm definition, clipped at 0).

Outputs: out/<tag>/maps.png (every sounding: SIF, nFLH, ΔSIF = SIF − nFLH; and relative ΔSIF vs nFLH
for nFLH > --rel-min),
out/<tag>/scatter.png (density, difference histogram, difference vs nFLH and vs latitude), summary.txt
"""
import argparse
import sys
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common import UNIT, add_filter_args, load_filtered, out_dir  # noqa: E402

p = argparse.ArgumentParser()
add_filter_args(p)
p.add_argument("--rel-min", type=float, default=0.05, help="relative ΔSIF only where nFLH > this")
A = p.parse_args()
OUT = out_dir(HERE, A)


def stats(x, y):
    s, b = np.polyfit(x, y, 1)
    d = y - x
    big = x > 0.05
    return dict(n=len(x), r=np.corrcoef(x, y)[0, 1], slope=s, icpt=b, mean_d=d.mean(), median_d=np.median(d),
                std_d=d.std(), p10_d=np.percentile(d, 10), p90_d=np.percentile(d, 90),
                ratio=np.median(y[big] / x[big]) if big.any() else np.nan, neg_sif=(y < 0).mean(),
                neg_nflh=(x <= 0).mean())


def binned(x, y, edges):
    c, q = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi)
        if m.sum() >= 50:
            c.append(0.5 * (lo + hi))
            q.append(np.percentile(y[m], [10, 25, 50, 75, 90]))
    return np.array(c), np.array(q)


data, lines = {}, [f"filter: --l2 {A.l2}  --query \"{A.query}\"", ""]
for date in A.dates:
    df, steps = load_filtered(A, date)
    data[date] = df
    s = stats(df.nflh.to_numpy(float), df.sif_678.to_numpy(float))
    lines += [f"{date}:"] + [f"   {k:<45s} {v:>12,d}" for k, v in steps]
    lines += [f"   SIF vs nflh: r = {s['r']:.3f}, SIF = {s['slope']:.3f} nFLH {s['icpt']:+.4f}",
              f"   SIF − nFLH: mean {s['mean_d']:+.4f}, median {s['median_d']:+.4f}, std {s['std_d']:.4f}, "
              f"10-90 % [{s['p10_d']:+.4f}, {s['p90_d']:+.4f}]  [{UNIT.replace('$', '')}]",
              f"   median SIF / nFLH (nFLH > 0.05) {s['ratio']:.3f};  SIF < 0: {s['neg_sif']:.1%};  "
              f"nFLH <= 0: {s['neg_nflh']:.1%}", ""]
(OUT / "summary.txt").write_text("\n".join(lines))
print("\n".join(lines))

# ------------------------------------------------------------------ maps (raw soundings) + relative ΔSIF vs nFLH
rng = np.random.default_rng(0)
fig = plt.figure(figsize=(23, 4.4 * len(A.dates)))
gs = fig.add_gridspec(len(A.dates), 4, width_ratios=[1, 1, 1, 0.8])
for i, (date, df) in enumerate(data.items()):
    o = rng.permutation(len(df))                 # random drawing order, so no granule hides another
    lon, lat = df.lon.to_numpy(float)[o], df.lat.to_numpy(float)[o]
    panels = [("SIF 678 nm", df.sif_678.to_numpy(float)[o], "viridis", (0, 0.5)),
              ("L2 nFLH", df.nflh.to_numpy(float)[o], "viridis", (0, 0.5)),
              ("ΔSIF = SIF − nFLH", (df.sif_678 - df.nflh).to_numpy(float)[o], "RdBu_r", (-0.2, 0.2))]
    for j, (title, val, cmap, (lo, hi)) in enumerate(panels):
        ax = fig.add_subplot(gs[i, j], projection=ccrs.PlateCarree())
        ax.set_global()
        ax.coastlines(lw=0.3)
        sc = ax.scatter(lon, lat, c=val, s=0.05, marker=",", lw=0, cmap=cmap, vmin=lo, vmax=hi,
                        transform=ccrs.PlateCarree(), rasterized=True)
        fig.colorbar(sc, ax=ax, shrink=0.75, pad=0.02, label=UNIT)
        ax.set_title(f"{date}  {title}  ({len(df):,} soundings)", fontsize=9)
    ax = fig.add_subplot(gs[i, 3])
    big = df.nflh.to_numpy(float) > A.rel_min
    x = df.nflh.to_numpy(float)[big]
    rel = 100 * ((df.sif_678 - df.nflh) / df.nflh).to_numpy(float)[big]
    ok = np.isfinite(x) & np.isfinite(rel)
    x, rel = x[ok], rel[ok]
    xr = (np.percentile(x, 0.5), np.percentile(x, 99.5))
    yr = (max(np.percentile(rel, 1), -200), min(np.percentile(rel, 99), 200))
    h = ax.hist2d(x, rel, bins=150, range=[xr, yr], norm=mcolors.LogNorm(), cmap="viridis", rasterized=True)
    fig.colorbar(h[3], ax=ax, label="soundings")
    c, q = binned(x, rel, np.linspace(*xr, 30))
    if len(c):
        ax.plot(c, q[:, 2], "w-", lw=1.6, label="median")
        ax.plot(c, q[:, 1], "w--", lw=0.8, label="25-75 %")
        ax.plot(c, q[:, 3], "w--", lw=0.8)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xlabel(f"L2 nFLH [{UNIT}]")
    ax.set_ylabel("relative ΔSIF = (SIF − nFLH) / nFLH  [%]")
    ax.set_title(f"{date}  relative ΔSIF vs nFLH (nFLH > {A.rel_min:g}: {big.sum():,})", fontsize=9)
    ax.legend(fontsize=7, loc="upper right")
fig.suptitle(f"filter: --l2 {A.l2}   --query \"{A.query}\"", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "maps.png", dpi=150)
plt.close(fig)

# ------------------------------------------------------------------ scatter / differences
fig, axs = plt.subplots(len(A.dates), 4, figsize=(21, 4.6 * len(A.dates)), squeeze=False)
for i, (date, df) in enumerate(data.items()):
    x, y = df.nflh.to_numpy(float), df.sif_678.to_numpy(float)
    d = y - x
    s = stats(x, y)
    ax = axs[i, 0]
    lim = (min(np.percentile(x, 0.1), np.percentile(y, 0.1), 0) - 0.02, max(np.percentile(x, 99.9),
                                                                            np.percentile(y, 99.9)))
    h = ax.hist2d(x, y, bins=200, range=[lim, lim], norm=mcolors.LogNorm(), cmap="viridis", rasterized=True)
    fig.colorbar(h[3], ax=ax, label="soundings")
    ax.plot(lim, lim, "k--", lw=0.8)
    xx = np.array(lim)
    ax.plot(xx, s["slope"] * xx + s["icpt"], "r", lw=1.0, label=f"OLS: SIF = {s['slope']:.2f} nFLH {s['icpt']:+.3f}")
    c, q = binned(x, y, np.linspace(lim[0], lim[1], 40))
    if len(c):
        ax.plot(c, q[:, 2], "w-", lw=1.5, label="median SIF per nFLH bin")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel(f"L2 nFLH [{UNIT}]")
    ax.set_ylabel(f"SIF 678 nm [{UNIT}]")
    ax.set_title(f"{date}: n = {s['n']:,}, r = {s['r']:.3f}", fontsize=10)
    ax.legend(fontsize=7, loc="upper left")

    ax = axs[i, 1]
    rng = (np.percentile(d, 0.5), np.percentile(d, 99.5))
    ax.hist(d, bins=200, range=rng, color="tab:blue")
    ax.axvline(0, color="k", lw=0.7)
    ax.axvline(s["median_d"], color="r", lw=1)
    ax.set_xlabel(f"SIF − nFLH [{UNIT}]")
    ax.set_ylabel("soundings")
    ax.set_title(f"median {s['median_d']:+.3f}, mean {s['mean_d']:+.3f}, std {s['std_d']:.3f}", fontsize=10)

    for ax, xv, xl in ((axs[i, 2], x, "L2 nFLH"), (axs[i, 3], df.lat.to_numpy(float), "latitude (°)")):
        edges = np.linspace(np.percentile(xv, 0.5), np.percentile(xv, 99.5), 40)
        c, q = binned(xv, d, edges)
        if len(c):
            ax.fill_between(c, q[:, 0], q[:, 4], color="tab:blue", alpha=0.15, label="10-90 %")
            ax.fill_between(c, q[:, 1], q[:, 3], color="tab:blue", alpha=0.3, label="25-75 %")
            ax.plot(c, q[:, 2], color="tab:blue", lw=1.5, label="median")
        ax.axhline(0, color="k", lw=0.7)
        ax.set_xlabel(xl)
        ax.set_ylabel(f"SIF − nFLH [{UNIT}]")
        ax.legend(fontsize=7)
    axs[i, 2].set_title("difference vs nFLH", fontsize=10)
    axs[i, 3].set_title("difference vs latitude", fontsize=10)
    for ax in axs[i]:
        ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"filter: --l2 {A.l2}   --query \"{A.query}\"", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "scatter.png", dpi=130)
print(f"wrote {OUT}/maps.png, scatter.png, summary.txt")
