"""Plots and a text summary for the AK files written by ak_selected_pixels.jl.

    python retrieval_comparison/plot_ak.py [run ...]      # default: every out/ak/ak_*.nc

For each run: AK heatmaps (one pixel per category), the SIF row of the AK, and the
posterior correlation of SIF with every other state element. With several runs: a
side-by-side summary per category.

AK convention: A[i, j] = d xhat_i / d x_true_j (A = S_post * H_obs at the final state).
A[sif, sif] = 1 means the retrieved SIF responds fully to true SIF; A[sif, j] != 0
(j != sif) is leakage of other elements into SIF. The "SIF shift per 1 sigma_post of
element j" is A[sif, j] * sqrt(S_post[j, j]), in SIF-coefficient units (x 0.9004 for the
678 nm radiance).

Panel (d) of the SIF-row figure: rho_sif,j * sigma_sif = S_post[sif,j] / sigma_post,j, the SIF
shift that goes with a 1 sigma_post error in element j. More generally, if element j is forced
(or truncated) to be off by delta_j, the best-fit SIF shifts by (S_post[sif,j]/S_post[j,j]) * delta_j
(Gaussian conditional mean). Check: dropping p3 (npoly2) shifts SIF by -slope*p3_hat, which
reproduces npoly2 - npoly3 at the selected pixels (r = 1.000).
"""
W678 = 0.9004  # sif_radiance_678nm = W678 * sif_ev1 (single SIF basis vector)
import sys
from pathlib import Path

import matplotlib
import netCDF4 as nc
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path(__file__).resolve().parent / "out" / "ak"
files = sorted(OUT.glob("ak_*.nc"))
if len(sys.argv) > 1:
    files = [f for f in files if any(f.name.startswith(f"ak_{r}_") for r in sys.argv[1:])]
if not files:
    sys.exit("no ak_*.nc found; run ak_selected_pixels.jl or ak_gchpio_ensemble.jl first")
OPTIONAL = ["sif678_true", "sif678_ret", "sigma678", "vza"]


def load(f):
    ds = nc.Dataset(f)
    ds.set_auto_mask(False)
    d = {k: ds[k][:] for k in ["A", "S_post", "x_hat", "dof_signal", "reduced_chi2", "converged", "converged_stored",
                               "sif_ev1_rerun", "sif_ev1_stored", "scan", "pixel", "category_id", "chl", "nflh",
                               "sigma_prior"]}
    for k in OPTIONAL:
        if k in ds.variables:
            d[k] = ds[k][:]
    d["names"] = ds.getncattr("state_names_csv").split(",")
    d["cats"] = ds.getncattr("category_names_csv").split(",")
    d["run"] = ds.getncattr("run")
    d["meta"] = ds.getncattr("meta_label") if "meta_label" in ds.ncattrs() else "chl"
    ds.close()
    return d


def corr_from_cov(S):
    sd = np.sqrt(np.diag(S))
    return S / np.outer(sd, sd)


def r2_on(C, i, grp):
    """Posterior R^2 of element i on the group `grp` (multiple correlation from the correlation matrix)."""
    c = C[i, grp]
    return float(c @ np.linalg.solve(C[np.ix_(grp, grp)], c))


R = {}
for f in files:
    d = load(f)
    R[d["run"]] = d
groups = {}
for run, d in R.items():
    groups.setdefault("gchpio" if run.startswith("gchpio") else "granule", []).append(run)
colors = plt.get_cmap("tab10")

for tag, runs in groups.items():
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    # ----------------------------------------------------------- text summary
    for run in runs:
        d = R[run]
        names, nst = d["names"], len(d["names"])
        i_sif = names.index("sif_ev1")
        leg = [n for n in names if n.startswith("legendre_p")]
        lg = [names.index(n) for n in leg]
        has_truth = "sif678_true" in d
        say(f"\n===== {run} =====  ({len(d['scan'])} samples, state = {nst} elements)")
        if np.all(np.isnan(d["sif_ev1_stored"])):
            say(f"no stored retrieval for this run: the AK re-run is the retrieval;  converged {np.mean(d['converged'] == 1) * 100:.0f}% of samples")
        else:
            ok = np.allclose(d["sif_ev1_rerun"], d["sif_ev1_stored"], rtol=1e-3, atol=1e-4)
            say(f"re-run reproduces stored retrieval: sif max |diff| = {np.max(np.abs(d['sif_ev1_rerun'] - d['sif_ev1_stored'])):.2e}"
                f" (coefficient units; all close: {ok});  converged flag agrees for {np.mean(d['converged'] == d['converged_stored']) * 100:.0f}%")
        hdr = f"{'category':16s}{'n':>3s}{'DOF_sig':>8s}{'σ_SIF678':>10s}{'R2(Leg)':>9s}" + \
              "".join(f"{'ρ(sif,' + n.replace('legendre_', '') + ')':>11s}" for n in leg) + f"{'ρ(sif,α)':>10s}{'max|ρ|PC':>10s}"
        if has_truth:
            hdr += f"{'mean err':>10s}{'rms err/σ':>11s}"
        say(hdr)
        for ic, cname in enumerate(d["cats"]):
            m = np.where(d["category_id"] == ic)[0]
            if len(m) == 0:
                continue
            C = np.array([corr_from_cov(d["S_post"][k]) for k in m])
            sig = np.sqrt(d["S_post"][m, i_sif, i_sif]) * W678
            pcs = [names.index(n) for n in names if n.startswith("pc_")]
            row = (f"{cname:16s}{len(m):3d}{d['dof_signal'][m].mean():8.1f}{sig.mean():10.4f}"
                   f"{np.mean([r2_on(c, i_sif, lg) for c in C]):9.3f}")
            row += "".join(f"{C[:, i_sif, j].mean():11.3f}" for j in lg)
            row += f"{C[:, i_sif, names.index('alpha')].mean():10.3f}{np.abs(C[:, i_sif, pcs]).max(axis=1).mean():10.3f}"
            if has_truth:
                err = d["sif678_ret"][m] - d["sif678_true"][m]
                row += f"{err.mean():10.4f}{np.sqrt(np.mean((err / d['sigma678'][m]) ** 2)):11.2f}"
            say(row)
        if has_truth:
            say("  (mean err = retrieved - true SIF678; rms err/σ ≈ 1 if the posterior σ captures the actual error)")

    # ---------------------------------------------------------------- figures
    for run in runs:
        d = R[run]
        names, nst = d["names"], len(d["names"])
        i_sif = names.index("sif_ev1")
        ncat = len(d["cats"])
        fig, axs = plt.subplots(1, ncat, figsize=(4.6 * ncat, 4.4), squeeze=False)
        for ic, cname in enumerate(d["cats"]):
            ax = axs[0, ic]
            idx = np.where(d["category_id"] == ic)[0]
            if len(idx) == 0:
                ax.set_axis_off()
                continue
            k = idx[0]
            im = ax.imshow(d["A"][k], cmap="RdBu_r", vmin=-1, vmax=1)
            ax.set_xticks(range(nst), names, rotation=90, fontsize=6)
            ax.set_yticks(range(nst), names, fontsize=6)
            ax.set_title(f"{cname}\n{d['meta']}={d['chl'][k]:.3f}  DOF={d['dof_signal'][k]:.1f}  "
                         f"A[sif,sif]={d['A'][k, i_sif, i_sif]:.2f}", fontsize=8)
        fig.colorbar(im, ax=axs.ravel().tolist(), shrink=0.8, label="A[i,j] = d xhat_i / d x_true_j")
        fig.suptitle(f"{run}: averaging kernel (rows = retrieved, columns = true)", fontsize=10)
        fig.savefig(OUT / f"AK_heatmaps_{run}.png", dpi=130, bbox_inches="tight")
        plt.close(fig)

        fig, axs = plt.subplots(4, 1, figsize=(11, 11.5), sharex=True)
        x = np.arange(nst)
        la = 0.7 if len(d["scan"]) <= 20 else max(0.12, 8.0 / len(d["scan"]))   # fainter lines for many samples
        for ic, cname in enumerate(d["cats"]):
            for k in np.where(d["category_id"] == ic)[0]:
                Sk = d["S_post"][k]
                axs[0].plot(x, d["A"][k, i_sif, :], "o-", ms=3, lw=0.8, color=colors(ic), alpha=la,
                            label=cname if k == np.where(d["category_id"] == ic)[0][0] else None)
                axs[1].plot(x, d["A"][k, i_sif, :] * np.sqrt(np.diag(Sk)), "o-", ms=3, lw=0.8, color=colors(ic), alpha=la)
                axs[2].plot(x, corr_from_cov(Sk)[i_sif, :], "o-", ms=3, lw=0.8, color=colors(ic), alpha=la)
                rho_sig = corr_from_cov(Sk)[i_sif, :] * np.sqrt(Sk[i_sif, i_sif]) * W678
                rho_sig[i_sif] = np.nan   # the SIF-SIF entry is just sigma_sif
                axs[3].plot(x, rho_sig, "o-", ms=3, lw=0.8, color=colors(ic), alpha=la)
        sig_sif = np.median(np.sqrt(d["S_post"][:, i_sif, i_sif])) * W678
        axs[0].set_ylabel("A[sif, j]  (SIF row of AK)")
        axs[1].set_ylabel("SIF shift per 1σ_post of j\n(A[sif,j]·σ_post,j)")
        axs[2].set_ylabel("posterior corr(sif, j)")
        axs[3].set_ylabel("ρ_sif,j · σ_sif  [W m$^{-2}$ µm$^{-1}$ sr$^{-1}$]\n= SIF$_{678}$ shift per 1σ_post of j")
        axs[3].axhline(sig_sif, color="0.4", ls="--", lw=0.8)
        axs[3].axhline(-sig_sif, color="0.4", ls="--", lw=0.8)
        axs[3].text(0.01, 0.93, "dashed: ±median σ_sif (1σ noise on SIF$_{678}$)", transform=axs[3].transAxes, fontsize=8,
                    color="0.3", va="top")
        for a in axs:
            a.axhline(0, color="k", lw=0.5)
        axs[3].set_xticks(x, names, rotation=90, fontsize=7)
        leg_ = axs[0].legend(fontsize=8, ncol=ncat)
        for h in leg_.legend_handles:
            h.set_alpha(1)
        for a, lab in zip(axs, "abcd"):
            a.text(0.995, 0.93, f"({lab})", transform=a.transAxes, ha="right", va="top", fontsize=10, fontweight="bold")
        fig.suptitle(f"{run}: how the retrieved SIF depends on, and correlates with, the other state elements", fontsize=10)
        fig.tight_layout()
        fig.savefig(OUT / f"AK_sif_row_{run}.png", dpi=130)
        plt.close(fig)

    # ------------------------------------- run comparison (same categories only)
    same = len(runs) > 1 and all(R[r]["cats"] == R[runs[0]]["cats"] for r in runs)
    if same:
        cats = R[runs[0]]["cats"]
        truth = all("sif678_true" in R[r] for r in runs)
        npan = 4 if truth else 3
        fig, axs = plt.subplots(1, npan, figsize=(5 * npan, 4.2))
        w = 0.8 / len(runs)
        for ir, run in enumerate(runs):
            d = R[run]
            names = d["names"]
            i_sif = names.index("sif_ev1")
            lg = [names.index(n) for n in names if n.startswith("legendre_p")]
            for ic in range(len(cats)):
                m = np.where(d["category_id"] == ic)[0]
                if len(m) == 0:
                    continue
                C = [corr_from_cov(d["S_post"][k]) for k in m]
                xpos = ic + (ir - (len(runs) - 1) / 2) * w
                axs[0].bar(xpos, np.mean([r2_on(c, i_sif, lg) for c in C]), w, color=colors(ir),
                           label=run if ic == 0 else None)
                axs[1].bar(xpos, np.mean(np.sqrt(d["S_post"][m, i_sif, i_sif]) * W678), w, color=colors(ir))
                axs[2].bar(xpos, d["dof_signal"][m].mean(), w, color=colors(ir))
                if truth:
                    err = d["sif678_ret"][m] - d["sif678_true"][m]
                    axs[3].bar(xpos, err.mean(), w, color=colors(ir), yerr=err.std() if len(m) > 1 else None, capsize=2)
        titles = ["posterior R² of SIF on Legendre terms", "σ_post(SIF$_{678}$)", "signal DOF tr(A)"]
        if truth:
            titles.append("retrieved − true SIF$_{678}$ (mean ± sd)")
            axs[3].axhline(0, color="k", lw=0.5)
        for a, t in zip(axs, titles):
            a.set_xticks(range(len(cats)), cats, rotation=15, fontsize=8)
            a.set_title(t, fontsize=10)
        axs[0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT / f"AK_run_comparison_{tag}.png", dpi=130)
        plt.close(fig)

    (OUT / f"summary_{tag}.txt").write_text("\n".join(lines) + "\n")
    print("wrote figures and", OUT / f"summary_{tag}.txt")
