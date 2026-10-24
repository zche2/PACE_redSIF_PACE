"""Compare the PACE Rrs recomputed from the L2 AOP files with the PACE_Rrs carried in the original CSV.

    python TROPOMI_OCI_colocation/field_rrs_matchups/plot_rrs_check.py [fullmatchups_Rrs_GUI.csv]

Needs matchup_summary.csv and matchup_boxes_5x5.csv from extract_matchup_boxes.py.
Recomputed = Bailey & Werdell filtered mean over the valid 5x5 pixels of the OCI Rrs, each pixel first
smoothed to the CSV wavelengths with a 10 nm Gaussian (PACE_Rrs_rep_<λ> in matchup_summary.csv).
Native = the same filtered mean on the 172 OCI bands, without smoothing.
Only field wavelengths with OCI Rrs behind them are compared (350-715 nm without 590-610 nm).

Outputs (out/): 4_rrs_check_scatter.png, 5_rrs_check_vs_wavelength.png, 6_rrs_check_examples.png,
7_rrs_check_drivers.png
"""
import re
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)
CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path("/home/zhe2/FraLab/PACE_redSIF_PACE/TROPOMI_OCI_colocation/fullmatchups_Rrs_GUI.csv")
GAP = (588.0, 613.0)           # L2 AOP Rrs has no bands here (blue/red CCD overlap)

orig = pd.read_csv(CSV, low_memory=False)
summ = pd.read_csv(HERE / "matchup_summary.csv").set_index("matchup_id")
# only the field wavelengths with OCI Rrs behind them (extract_matchup_boxes.py writes no others):
# the CSV's values below 346 nm and at 590-610 nm are extrapolated / interpolated and are left out
lam = np.array(sorted(float(m.group(1)) for c in summ.columns if (m := re.fullmatch(r"PACE_Rrs_rep_(\d+)", c))))
LO, HI = lam.min(), lam.max()
O = orig.loc[summ.index, [f"PACE_Rrs_{int(w)}" for w in lam]].to_numpy(float)           # original CSV
Rp = summ[[f"PACE_Rrs_rep_{int(w)}" for w in lam]].to_numpy(float)                       # recomputed
inr = np.ones(lam.size, bool)
ok = np.isfinite(O) & np.isfinite(Rp)
D = Rp - O
rel = D / np.where(np.abs(O) > 1e-5, np.abs(O), np.nan)

# ------------------------------------------------------------ 1: scatter
fig, axs = plt.subplots(1, 2, figsize=(12, 5.4))
x, y, c = O[ok], Rp[ok], np.broadcast_to(lam, O.shape)[ok]
for ax, scale in zip(axs, ("linear", "log")):
    sc = ax.scatter(x, y, c=c, s=4, cmap="turbo", alpha=0.6, rasterized=True)
    lim = (max(np.nanmin(np.abs(x[x > 0])), 1e-5), np.nanmax([x.max(), y.max()]) * 1.05) if scale == "log" else \
        (min(x.min(), y.min()), max(x.max(), y.max()) * 1.05)
    ax.plot(lim, lim, "k--", lw=0.8)
    ax.set_xscale(scale)
    ax.set_yscale(scale)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("PACE Rrs in original CSV (sr$^{-1}$)")
    ax.set_ylabel("PACE Rrs recomputed from L2 AOP (sr$^{-1}$)")
    ax.set_title(f"{scale} axes", fontsize=10)
fig.colorbar(sc, ax=axs, label="wavelength (nm)", shrink=0.85)
r = np.corrcoef(x, y)[0, 1]
fig.suptitle(f"{ok.any(axis=1).sum()} matchups, {ok.sum()} points ({LO:.0f}-{HI:.0f} nm):  r = {r:.5f},  "
             f"median |rel. diff| = {np.nanmedian(np.abs(rel[ok])) * 100:.2f} %,  "
             f"median diff = {np.nanmedian(D[ok]):+.1e} sr$^{{-1}}$", fontsize=10)
fig.savefig(OUT / "4_rrs_check_scatter.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------- 2: difference vs wavelength
fig, axs = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
for ax, arr, lab, sc_ in ((axs[0], D, "recomputed − original (sr$^{-1}$)", 1.0),
                          (axs[1], rel * 100, "relative difference (%)", 1.0)):
    a = np.where(ok, arr, np.nan)
    q10, q50, q90 = (np.nanpercentile(a, q, axis=0) for q in (10, 50, 90))
    seg = np.where(lam > GAP[0], 1, 0)           # break the lines across the gap
    for s_ in (0, 1):
        m_ = seg == s_
        ax.fill_between(lam[m_], q10[m_], q90[m_], color="tab:blue", alpha=0.25,
                        label="10-90 % of matchups" if s_ == 0 else None)
        ax.plot(lam[m_], q50[m_], color="tab:blue", lw=1.5, label="median" if s_ == 0 else None)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylabel(lab)
    ax.spines[["top", "right"]].set_visible(False)
axs[1].set_ylim(-8, 8)
axs[1].set_xlabel("Wavelength (nm)")
axs[0].legend(fontsize=8)
for a_ in axs:
    a_.axvspan(*GAP, color="0.9", zorder=0)
axs[0].text(np.mean(GAP), axs[0].get_ylim()[1], "no OCI L2\nRrs bands\n588-613 nm", ha="center", va="top",
            fontsize=7, color="0.4")
fig.suptitle("Recomputed vs original PACE Rrs as a function of wavelength", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "5_rrs_check_vs_wavelength.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------- 3: example spectra
box = pd.read_csv(HERE / "matchup_boxes_5x5.csv")
rrs_cols = [c for c in box.columns if c.startswith("Rrs_")]
lam_oci = np.array([float(c.split("_")[1]) for c in rrs_cols])
lam_f = np.array(sorted(float(m.group(1)) for c in orig.columns if (m := re.fullmatch(r"FIELD_Rrs_(\d+)", c))))
med_rel = summ["rrs_check_median_reldiff"]
picks = [("best", med_rel.idxmin()), ("typical", (med_rel - med_rel.median()).abs().idxmin()), ("worst", med_rel.idxmax())]
fig, axs = plt.subplots(2, 3, figsize=(15, 7.5), gridspec_kw=dict(height_ratios=[2.2, 1]), sharex=True)
for j, (label, mid) in enumerate(picks):
    b = box[(box.matchup_id == mid) & box.valid]
    X = b[rrs_cols].to_numpy(float)
    med = np.median(X, 0)
    sd = X.std(0, ddof=1)
    nat = np.array([X[np.abs(X[:, q] - med[q]) <= 1.5 * sd[q], q].mean() for q in range(X.shape[1])])
    k = summ.index.get_loc(mid)
    ax = axs[0, j]
    gap_nan = lambda x_, y_: np.where(np.r_[False, np.diff(x_) > GAP[1] - GAP[0] - 1], np.nan, y_)  # noqa: E731
    ax.plot(lam_oci, gap_nan(lam_oci, nat), color="0.6", lw=0.8, label="native OCI bands (filtered mean)")
    ax.plot(lam, gap_nan(lam, O[k]), "o-", color="k", ms=2.5, lw=1.0, label="original CSV PACE_Rrs")
    ax.plot(lam, gap_nan(lam, Rp[k]), "-", color="tab:red", lw=1.2, label="recomputed (10 nm smoothing)")
    fr, fu = (orig.loc[mid, [f"{p_}_{int(w)}" for w in lam_f]].to_numpy(float) for p_ in ("FIELD_Rrs", "FIELD_Rrs_unc"))
    f_ = np.isfinite(fr)
    ax.fill_between(lam_f[f_], (fr - np.nan_to_num(fu))[f_], (fr + np.nan_to_num(fu))[f_], color="tab:green", alpha=0.25,
                    lw=0)
    ax.plot(lam_f[f_], fr[f_], "-", color="tab:green", lw=1.2, marker="o", ms=2.5,
            label="in situ FIELD_Rrs (dots = measured λ, ± unc)")
    ys = np.concatenate([O[k], Rp[k], fr[f_], nat])           # axis from the spectra; a wide unc band is clipped
    y0, y1 = np.nanmin(ys), np.nanmax(ys)
    ax.set_ylim(y0 - 0.08 * (y1 - y0), y1 + 0.08 * (y1 - y0))
    s_ = summ.loc[mid]
    ax.set_title(f"{label}: matchup {mid} ({s_.data_type}), median |rel. diff| {s_.rrs_check_median_reldiff * 100:.1f} %\n"
                 f"valid pixels {s_.n_valid} vs CSV {s_.PACE_valid_pixels_csv},  CV {s_.filtered_cv_400_600:.3f}", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    axs[1, j].plot(lam, gap_nan(lam, (Rp[k] - O[k]) / np.abs(O[k]) * 100), color="tab:red", lw=1)
    axs[1, j].axhline(0, color="k", lw=0.5)
    axs[1, j].set_ylim(-20, 20)
    for a_ in axs[:, j]:
        a_.axvspan(*GAP, color="0.92", zorder=0)
    axs[1, j].set_xlabel("Wavelength (nm)")
axs[0, 0].set_ylabel("Rrs (sr$^{-1}$)")
axs[1, 0].set_ylabel("rel. diff (%)")
axs[0, 0].legend(fontsize=7)
fig.suptitle(f"Example spectra: PACE on the {lam.size} field wavelengths with OCI Rrs ({LO:.0f}-{HI:.0f} nm; "
             f"grey = 588-613 nm, no OCI L2 Rrs bands); in situ on all its wavelengths.  "
             f"Bottom: recomputed vs CSV PACE_Rrs", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "6_rrs_check_examples.png", dpi=150)
plt.close(fig)

# -------------------------------------------------------------- 4: drivers
fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
same = summ.n_valid == summ.PACE_valid_pixels_csv
mr = summ.rrs_check_median_reldiff * 100
axs[0].boxplot([mr[same].dropna(), mr[~same].dropna()], tick_labels=["counts agree", "counts differ"], showfliers=True)
axs[0].set_ylabel("median |rel. diff| per matchup (%)")
axs[0].set_title(f"valid-pixel count vs CSV ({same.sum()} agree, {(~same).sum()} differ)", fontsize=9)
axs[1].scatter(summ.filtered_cv_400_600, mr, s=12, c=np.where(same, "tab:blue", "tab:orange"))
axs[1].set_xlabel("box CV (filtered, 400-600 nm)")
axs[1].set_title("blue = counts agree, orange = differ", fontsize=9)
for t, g in summ.groupby("data_type"):
    axs[2].scatter([t] * len(g), g.rrs_check_median_reldiff * 100, s=12, alpha=0.6)
axs[2].tick_params(axis="x", rotation=20)
axs[2].set_title("by platform type", fontsize=9)
for a_ in axs[1:]:
    a_.set_ylabel("median |rel. diff| (%)")
for a_ in axs:
    a_.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "7_rrs_check_drivers.png", dpi=150)
plt.close(fig)

a = np.where(ok, rel, np.nan)
print(f"points {ok.sum()} in {LO:.0f}-{HI:.0f} nm; r = {r:.5f}; median rel. diff {np.nanmedian(a) * 100:+.2f} % "
      f"(|.|: {np.nanmedian(np.abs(a)) * 100:.2f} %)")
for w0, w1 in ((350, 400), (400, 500), (500, 585), (615, 660), (660, 715)):
    m = (lam >= w0) & (lam <= w1)
    print(f"  {w0}-{w1} nm: median rel. diff {np.nanmedian(a[:, m]) * 100:+.2f} %, median |rel. diff| "
          f"{np.nanmedian(np.abs(a[:, m])) * 100:.2f} %")
print(f"median |rel. diff| per matchup: counts agree {mr[same].median():.2f} %, counts differ {mr[~same].median():.2f} %")
print("wrote", ", ".join(f.name for f in sorted(OUT.glob("[4-7]_rrs_check*.png"))))
