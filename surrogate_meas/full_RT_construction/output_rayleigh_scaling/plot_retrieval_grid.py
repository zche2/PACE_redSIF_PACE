"""Retrieved vs true SIF on the Rayleigh × aerosol × SIF × view grid (output of retrieve_rt_grid.jl).

    python plot_retrieval_grid.py [--case noisy|clean] [--aer-column 100]

rt_retrieval_grid_AerColumn_<k>_scatter_colored_<case>.png
    retrieved vs true water-leaving SIF at 678.2 nm (truth: point value at the band centre, not the
    band average); one row per SVD config, columns coloured by
    VZA, aerosol scale and Rayleigh scale (filled: converged, hollow: not converged)
rt_retrieval_grid_AerColumn_<k>_dsif_vs_od_<case>.png
    ΔSIF against AOD + Rayleigh OD, one column per SIF intensity; top row TOA − water-leaving SIF,
    then retrieved − true SIF per config; colour = Rayleigh OD, marker = VZA, each line spans the aerosol scales
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.lines import Line2D

HERE = Path(__file__).absolute().parent
CONFIGS = {
    "nPC15_npoly3": "640–756 nm, nPC15, npoly3",
    "nPC15_npoly3_O2A": "640–780 nm (O$_2$-A), nPC15, npoly3",
    "nPC15_npoly5": "640–756 nm, nPC15, npoly5",
    "nPC15_npoly5_O2A": "640–780 nm (O$_2$-A), nPC15, npoly5",
}
UNITS = r"W m$^{-2}$ sr$^{-1}$ µm$^{-1}$"
ORDER = ("sif", "view", "rayleigh", "aerosol")
MARKERS = ["o", "s", "^", "D", "v"]


def load(aer_column):
    out = {}
    for tag, label in CONFIGS.items():
        path = HERE / f"rt_retrieval_grid_AerColumn_{aer_column}_{tag}.nc"
        if path.exists():
            out[label] = xr.open_dataset(path).transpose(*ORDER)
        else:
            print(f"missing {path.name}, skipped")
    if not out:
        raise SystemExit("no retrieval files found")
    return out


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)


def scene_text(d):
    return (f"SZA {d.attrs['sza']:g}°, albedo 0.02, GCHPIO column {d.attrs['aer_column']}, "
            f"true SIF = point value at {d.attrs['band_nm']:.1f} nm")


def plot_scatter(R, case, out):
    d0 = next(iter(R.values()))
    truth = np.broadcast_to(d0.sif_true.values[:, None, None, None], d0.sif_toa.shape)
    idx = np.indices(truth.shape)
    groups = [("VZA", 1, [f"{v:g}°" for v in d0.vza.values], "viridis"),
              ("aerosol scale (AOD)", 3, [f"×{a:g} ({t:.3f})" for a, t in zip(d0.aerosol_scale.values, d0.aod.values)], "plasma"),
              ("Rayleigh scale (τ$_R$)", 2, [f"×{r:g} ({t:.3g})" for r, t in zip(d0.rayleigh_scale.values, d0.rayleigh_od.values)], "cividis")]
    lim = (-0.12, 1.15)
    fig, axs = plt.subplots(len(R), 3, figsize=(17, 5.2 * len(R)), layout="constrained",
                            sharex=True, sharey=True, squeeze=False)
    for row, (name, d) in enumerate(R.items()):
        ret = d[f"sif_ret_{case}"].values
        ok = d[f"status_{case}"].values == 1
        for col, (gname, axis, labels, cmap) in enumerate(groups):
            ax = axs[row, col]
            colors = plt.get_cmap(cmap)(np.linspace(0.05, 0.85, len(labels)))
            for k, c in enumerate(colors):
                m = idx[axis] == k
                x, y, e = truth[m], ret[m], ret[m] - truth[m]
                ax.scatter(x[ok[m]], y[ok[m]], s=22, color=c, alpha=0.8, lw=0)
                ax.scatter(x[~ok[m]], y[~ok[m]], s=22, facecolor="none", edgecolor=c, alpha=0.8, lw=0.9)
                ax.scatter([], [], color=c, label=f"{labels[k]}: bias {e.mean():+.3f}, RMSE {np.sqrt((e ** 2).mean()):.3f}")
            ax.plot(lim, lim, "k--", lw=0.8)
            ax.scatter([], [], color="0.4", label="filled: converged")
            ax.scatter([], [], facecolor="none", edgecolor="0.4", label="hollow: not converged")
            ax.legend(title=gname, fontsize=7.5, title_fontsize=8, frameon=False, loc="upper left")
            ax.set_title(f"({chr(97 + 3 * row + col)}) {name}\ncoloured by {gname.split(' (')[0]}", loc="left", fontsize=10)
            ax.set_xlim(*lim)
            ax.set_ylim(*lim)
            style(ax)
            if row == len(R) - 1:
                ax.set_xlabel(f"True water-leaving SIF ({UNITS})")
            if col == 0:
                ax.set_ylabel(f"Retrieved SIF ({UNITS})")
    fig.suptitle(f"SVD retrieval vs truth on the RT grid ({case} spectra; {scene_text(d0)})", fontsize=11)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def plot_dsif(R, case, out):
    d0 = next(iter(R.values()))
    truth = d0.sif_true.values
    tau_r, aod = d0.rayleigh_od.values, d0.aod.values
    tau_tot = tau_r[:, None] + aod[None, :]                       # (rayleigh, aerosol)
    truth_band = d0.sif_true_bandavg.values                       # sif_toa is band-averaged too
    rows = [("TOA SIF − water-leaving SIF\n(band-averaged)", d0.sif_toa.values - truth_band[:, None, None, None])] + \
           [(f"Retrieved − true SIF\n{name}", d[f"sif_ret_{case}"].values - truth[:, None, None, None])
            for name, d in R.items()]
    n_s, n_r, n_v = len(truth), len(tau_r), d0.sizes["view"]
    colors = plt.get_cmap("cividis")(np.linspace(0.05, 0.85, n_r))
    fig, axs = plt.subplots(len(rows), n_s, figsize=(3.1 * n_s, 3.0 * len(rows)), layout="constrained",
                            sharex=True, sharey="row", squeeze=False)
    for r_i, (rname, dS) in enumerate(rows):
        for s_i in range(n_s):
            ax = axs[r_i, s_i]
            for ir in range(n_r):
                for iv in range(n_v):
                    ax.plot(tau_tot[ir], dS[s_i, iv, ir], color=colors[ir], marker=MARKERS[iv % len(MARKERS)],
                            ms=5, lw=1, alpha=0.9)
            ax.axhline(0, color="0.5", lw=0.6)
            ax.set_xscale("symlog", linthresh=0.01)
            ax.set_xlim(-0.002, 0.4)
            style(ax)
            if r_i == 0:
                ax.set_title(f"true SIF = {truth[s_i]:.3g}", fontsize=10)
            if r_i == len(rows) - 1:
                ax.set_xlabel(f"AOD + Rayleigh OD at {d0.attrs['band_nm']:.0f} nm")
            if s_i == 0:
                ax.set_ylabel(f"{rname}\n({UNITS})", fontsize=8.5)
    handles = [Line2D([], [], color=colors[i], lw=2, label=f"τ$_R$ = {t:.3g} (×{r:g})")
               for i, (r, t) in enumerate(zip(d0.rayleigh_scale.values, tau_r))] + \
              [Line2D([], [], color="0.3", marker=MARKERS[i % len(MARKERS)], ls="", label=f"VZA {v:g}°")
               for i, v in enumerate(d0.vza.values)]
    fig.legend(handles=handles, loc="outside lower center", ncol=len(handles), frameon=False, fontsize=9,
               title="colour: Rayleigh OD;  marker: VZA;  points along each line: AOD " + ", ".join(f"{a:.3f}" for a in aod))
    fig.suptitle(f"ΔSIF vs AOD + Rayleigh OD; retrieval rows use {case} spectra ({scene_text(d0)}; {UNITS})", fontsize=11)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def main(case, aer_column):
    R = load(aer_column)
    stem = f"rt_retrieval_grid_AerColumn_{aer_column}"
    plot_scatter(R, case, HERE / f"{stem}_scatter_colored_{case}.png")
    plot_dsif(R, case, HERE / f"{stem}_dsif_vs_od_{case}.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default="noisy", choices=["noisy", "clean"])
    ap.add_argument("--aer-column", type=int, default=100)
    args = ap.parse_args()
    main(args.case, args.aer_column)
