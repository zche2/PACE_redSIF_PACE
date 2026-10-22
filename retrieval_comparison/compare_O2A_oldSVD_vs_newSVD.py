"""Compare two O2A SVD retrievals of one PACE granule: old (Aug01 transmittance
SVD) vs newSVD (scattering-included transmittance SVD).

Metrics: convergence, posterior sigma / SNR of SIF@678, SIF magnitude,
low-chl SIF (expected ~0), and SIF vs nFLH.

    python retrieval_comparison/compare_O2A_oldSVD_vs_newSVD.py [YYYYMMDDTHHMMSS]

chl / nFLH come from the L2 products (same 1710x1272 grid) under
L2_ROOT/L2_{BGC,AOP}_<L2_VERSION>/YYYY/MM/DD/.
"""
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

STAMP = sys.argv[1] if len(sys.argv) > 1 else "20250724T170203"
DATE = f"{STAMP[:4]}/{STAMP[4:6]}/{STAMP[6:8]}"
ROOT = Path("/home/zhe2/data/PACE")
L2_VERSION = "V3.2"
L2_ROOT = Path("/kiwi-data/Data/satellite/PACE_OCI")
OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

RUNS = {
    "old": ROOT / "new_svd_retrieval_output_O2A" / DATE / f"interim_{STAMP}_svd_retrieval.nc",
    "newSVD": ROOT / "new_svd_retrieval_output_O2A_newSVD" / DATE / f"interim_{STAMP}_svd_retrieval.nc",
}
BGC = L2_ROOT / f"L2_BGC_{L2_VERSION}" / DATE / f"PACE_OCI.{STAMP}.L2.OC_BGC.{L2_VERSION.replace('.', '_')}.nc"
AOP = L2_ROOT / f"L2_AOP_{L2_VERSION}" / DATE / f"PACE_OCI.{STAMP}.L2.OC_AOP.{L2_VERSION.replace('.', '_')}.nc"

BAD_FLAGS = ["ATMFAIL", "LAND", "HIGLINT", "HILT", "HISATZEN", "STRAYLIGHT",
             "CLDICE", "HISOLZEN", "CHLFAIL"]
STATUS_NAMES = {0: "max-iter (not conv.)", 1: "converged", 2: "LM stalled", 4: "failed"}
SIF_STATE_IDX = 20  # sif_ev1 is the last of the 21 state elements

lines = []


def say(s=""):
    print(s)
    lines.append(s)


def load_run(path):
    ds = nc.Dataset(path)
    ds.set_auto_mask(False)
    d = {k: ds[k][:] for k in ["latitude", "longitude", "converged", "status_code", "n_steps",
                               "rmse", "reduced_chi2", "dof", "objective", "sif_radiance_678nm",
                               "sif_ev1", "is_dark", "is_ocean"]}
    d["var_sif"] = ds["S_posterior_diag"][SIF_STATE_IDX]
    d["x_sif"] = ds["x_hat"][SIF_STATE_IDX]
    ds.close()
    d["sigma_sif"] = np.abs(np.where(d["sif_ev1"] != 0, d["sif_radiance_678nm"] / d["sif_ev1"], np.nan)) \
        * np.sqrt(np.where(d["var_sif"] > 0, d["var_sif"], np.nan))
    return d


def load_l2():
    with nc.Dataset(BGC) as ds:
        g = ds["geophysical_data"]
        g.set_auto_mask(True)
        chl = np.ma.filled(g["chlor_a"][:].astype(float), np.nan)
        fl = g["l2_flags"]
        flags = np.ma.filled(fl[:], 0).astype(np.int64)
        names = fl.getncattr("flag_meanings").split()
        masks = np.atleast_1d(fl.getncattr("flag_masks"))
        lat = np.ma.filled(ds["navigation_data"]["latitude"][:], np.nan)
    bad = 0
    for n in BAD_FLAGS:
        if n in names:
            bad |= int(masks[names.index(n)])
    with nc.Dataset(AOP) as ds:
        g = ds["geophysical_data"]
        g.set_auto_mask(True)
        nflh = np.ma.filled(g["nflh"][:].astype(float), np.nan)
    chl[chl <= 0] = np.nan
    return chl, nflh, (flags & bad) == 0, lat


def pct(x, q=(5, 25, 50, 75, 95)):
    return " ".join(f"{v:9.4f}" for v in np.nanpercentile(x, q))


def corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return np.corrcoef(a[ok], b[ok])[0, 1], ok.sum()


def slope(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    A = np.vstack([x[ok], np.ones(ok.sum())]).T
    return np.linalg.lstsq(A, y[ok], rcond=None)[0]


R = {k: load_run(v) for k, v in RUNS.items()}
chl, nflh, l2_good, lat_l2 = load_l2()
O, N = R["old"], R["newSVD"]

say(f"Granule {STAMP}   L2 version {L2_VERSION}   grid {O['sif_ev1'].shape}")
say(f"geolocation check |lat_ret - lat_L2| median = {np.nanmedian(np.abs(O['latitude'] - lat_l2)):.2e} deg")
say(f"sanity: sif_ev1 == x_hat[sif]  old {np.allclose(O['sif_ev1'], O['x_sif'])}  new {np.allclose(N['sif_ev1'], N['x_sif'])}")
say(f"sanity: ocean/dark masks identical across runs: "
    f"{np.array_equal(O['is_ocean'], N['is_ocean']) and np.array_equal(O['is_dark'], N['is_dark'])}")

ocean = O["is_ocean"] == 1
# is_dark = max(radiance) <= dark_max_radiance: dark = clear open water, bright = cloud/glint/edge.
dark = O["is_dark"] == 1
base = ocean & l2_good & np.isfinite(chl)             # clear ocean, dark or bright
say(f"\npixels: total {ocean.size}  ocean {ocean.sum()} (dark {(ocean & dark).sum()}, bright {(ocean & ~dark).sum()})  "
    f"clear = L2 flags ok & chl valid {base.sum()} (dark {(base & dark).sum()}, bright {(base & ~dark).sum()})")

# ---------------------------------------------------------------- convergence
say("\n=== 1. CONVERGENCE (clear ocean: L2 flags ok & chl valid) ===")
say(f"{'':22s}{'old':>12s}{'newSVD':>12s}")
for code, name in STATUS_NAMES.items():
    say(f"{name:22s}" + "".join(f"{100 * np.mean(R[k]['status_code'][base] == code):11.1f}%" for k in R))
say(f"{'converged flag':22s}" + "".join(f"{100 * np.mean(R[k]['converged'][base] == 1):11.1f}%" for k in R))
for key, lab in [("n_steps", "n_steps (median)"), ("reduced_chi2", "reduced chi2 (median)"),
                 ("rmse", "rmse (median)"), ("dof", "residual DOF n-trA (median)")]:
    say(f"{lab:22s}" + "".join(f"{np.nanmedian(R[k][key][base]):12.4g}" for k in R))
say(f"{'reduced chi2 p90':22s}" + "".join(f"{np.nanpercentile(R[k]['reduced_chi2'][base], 90):12.4g}" for k in R))
say("\nconvergence by subset (bright pixels are mostly cloud / swath edge; dark = open water):")
subsets = {
    "ocean, all": ocean,
    "ocean & dark": ocean & dark,
    "ocean & bright (mostly cloud)": ocean & ~dark,
    "clear [= main set]": base,
    "  clear & dark": base & dark,
    "  clear & bright": base & ~dark,
}
say(f"{'subset':36s}{'n':>9s}  conv% old/new   redchi2 med old/new")
for lab, sel in subsets.items():
    say(f"{lab:36s}{sel.sum():9d}  {100 * np.mean(O['converged'][sel] == 1):5.1f}/{100 * np.mean(N['converged'][sel] == 1):5.1f}"
        f"      {np.nanmedian(O['reduced_chi2'][sel]):6.2f}/{np.nanmedian(N['reduced_chi2'][sel]):6.2f}")
both_conv = base & (O["converged"] == 1) & (N["converged"] == 1)
only_old = base & (O["converged"] == 1) & (N["converged"] != 1)
only_new = base & (O["converged"] != 1) & (N["converged"] == 1)
say(f"converged in both {both_conv.sum()}   only old {only_old.sum()}   only newSVD {only_new.sum()}")
say("reduced chi2 on converged-in-both: " + "  ".join(
    f"{k}={np.nanmedian(R[k]['reduced_chi2'][both_conv]):.4g}" for k in R))

# ----------------------------------------------------------- sigma / SNR
say("\n=== 2. POSTERIOR sigma(SIF@678) and SIF/sigma (clear ocean; W m-2 um-1 sr-1) ===")
say(f"{'':22s}{'p5':>9s} {'p25':>9s} {'p50':>9s} {'p75':>9s} {'p95':>9s}")
for k in R:
    say(f"{'sigma ' + k:22s}" + pct(R[k]["sigma_sif"][base]))
for k in R:
    snr = np.abs(R[k]["sif_radiance_678nm"]) / R[k]["sigma_sif"]
    say(f"{'|SIF|/sigma ' + k:22s}" + pct(snr[base]))
ratio = N["sigma_sif"] / O["sigma_sif"]
say(f"sigma ratio new/old per pixel: median {np.nanmedian(ratio[base]):.3f}  p5/p95 "
    f"{np.nanpercentile(ratio[base], 5):.3f}/{np.nanpercentile(ratio[base], 95):.3f}")

# -------------------------------------------------------- SIF magnitude
say("\n=== 3. SIF@678 MAGNITUDE (converged in both) ===")
so, sn = O["sif_radiance_678nm"], N["sif_radiance_678nm"]
say(f"{'':22s}{'p5':>9s} {'p25':>9s} {'p50':>9s} {'p75':>9s} {'p95':>9s}{'mean':>10s}{'std':>9s}")
for k, s in (("old", so), ("newSVD", sn)):
    v = s[both_conv]
    say(f"{k:22s}" + pct(v) + f"{np.nanmean(v):10.4f}{np.nanstd(v):9.4f}")
d = (sn - so)[both_conv]
r, n = corr(so[both_conv], sn[both_conv])
sl, ic = slope(so[both_conv], sn[both_conv])
say(f"new-old: mean {d.mean():.4f}  median {np.median(d):.4f}  rms {np.sqrt(np.mean(d**2)):.4f}  "
    f"corr {r:.3f}  fit new = {sl:.3f}*old + {ic:.4f}  (n={n})")
say(f"fraction SIF<0:  old {100 * np.mean(so[both_conv] < 0):.1f}%  newSVD {100 * np.mean(sn[both_conv] < 0):.1f}%")

# --------------------------------------------------------------- low chl
say("\n=== 4. LOW-CHL SIF (true SIF should be ~0 to small positive) ===")
say("chl bin [mg m-3]      n(conv-both)  nFLH med   median old/new      mean old/new        std old/new     %<0 old/new")
bins = [(0, 0.05), (0.05, 0.1), (0.1, 0.3), (0.3, 1), (1, 100)]
for lo, hi in bins:
    m = both_conv & (chl >= lo) & (chl < hi)
    if m.sum() < 50:
        continue
    a, b = so[m], sn[m]
    say(f"{lo:5.2f}-{hi:<6.2f}      {m.sum():9d}   {np.nanmedian(nflh[m]):8.4f}   {np.median(a):7.4f}/{np.median(b):7.4f}   "
        f"{a.mean():7.4f}/{b.mean():7.4f}   {a.std():7.4f}/{b.std():7.4f}   "
        f"{100 * np.mean(a < 0):4.1f}/{100 * np.mean(b < 0):4.1f}")
lc = both_conv & (chl < 0.1)
say(f"all chl<0.1: n={lc.sum()}  mean SIF old {so[lc].mean():.4f} new {sn[lc].mean():.4f};  "
    f"sigma median old {np.nanmedian(O['sigma_sif'][lc]):.4f} new {np.nanmedian(N['sigma_sif'][lc]):.4f}")

# ------------------------------------------------------------------ nFLH
say("\n=== 5. SIF vs nFLH (converged in both, nFLH valid) ===")
ok = both_conv & np.isfinite(nflh)
say(f"n={ok.sum()}  nFLH p5/50/95 = {np.nanpercentile(nflh[ok], [5, 50, 95])}")
for k, s in (("old", so), ("newSVD", sn)):
    r, _ = corr(s[ok], nflh[ok])
    sl, ic = slope(nflh[ok], s[ok])
    say(f"{k:8s} pearson r = {r:.3f}   SIF = {sl:.3f}*nFLH + {ic:.4f}")
for lo, hi in [(0, 0.1), (0.1, 0.3), (0.3, 100)]:
    m = ok & (chl >= lo) & (chl < hi)
    if m.sum() < 200:
        continue
    say(f"  chl {lo}-{hi}: " + "  ".join(
        f"{k} r={corr(s[m], nflh[m])[0]:.3f} slope={slope(nflh[m], s[m])[0]:.3f}" for k, s in (("old", so), ("newSVD", sn))))

# ---------------------------------------------------------------- figures
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
for k, c in (("old", "C0"), ("newSVD", "C3")):
    cnt = [np.mean(R[k]["status_code"][base] == code) * 100 for code in STATUS_NAMES]
    ax[0].bar(np.arange(len(cnt)) + (0.2 if k == "newSVD" else -0.2), cnt, 0.4, color=c, label=k)
    ax[1].hist(np.clip(R[k]["reduced_chi2"][base], 0, 10), 80, histtype="step", color=c, label=k)
    ax[2].hist(R[k]["n_steps"][base], np.arange(0, R[k]["n_steps"][base].max() + 2), histtype="step", color=c, label=k)
ax[0].set_xticks(range(len(STATUS_NAMES)), [STATUS_NAMES[c] for c in STATUS_NAMES], rotation=20, fontsize=8)
ax[0].set_ylabel("% of pixels")
ax[1].set_xlabel("reduced chi2 (clipped at 10)")
ax[2].set_xlabel("n_steps")
ax[0].legend()
fig.tight_layout()
fig.savefig(OUT / "1_convergence.png", dpi=130)

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for k, c in (("old", "C0"), ("newSVD", "C3")):
    ax[0].hist(np.log10(R[k]["sigma_sif"][base]), 100, histtype="step", color=c, label=k)
    ax[1].hist(np.clip(np.abs(R[k]["sif_radiance_678nm"][base]) / R[k]["sigma_sif"][base], 0, 30), 100,
               histtype="step", color=c, label=k)
ax[0].set_xlabel("log10 sigma(SIF@678)")
ax[1].set_xlabel("|SIF| / sigma (clipped at 30)")
ax[0].legend()
fig.tight_layout()
fig.savefig(OUT / "2_sigma_snr.png", dpi=130)

fig, ax = plt.subplots(1, 2, figsize=(10, 4.5))
h = ax[0].hist2d(so[both_conv], sn[both_conv], 150, range=[[-0.3, 0.6]] * 2, cmin=1, norm=matplotlib.colors.LogNorm())
ax[0].plot([-0.3, 0.6], [-0.3, 0.6], "k--", lw=0.8)
ax[0].set_xlabel("SIF@678 old")
ax[0].set_ylabel("SIF@678 newSVD")
for k, c in (("old", "C0"), ("newSVD", "C3")):
    ax[1].hist(R[k]["sif_radiance_678nm"][both_conv], 150, range=(-0.3, 0.6), histtype="step", color=c, label=k)
    ax[1].hist(R[k]["sif_radiance_678nm"][lc], 150, range=(-0.3, 0.6), histtype="step", color=c, ls=":",
               label=f"{k}, chl<0.1")
ax[1].axvline(0, color="k", lw=0.5)
ax[1].set_xlabel("SIF@678")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "3_sif_magnitude_lowchl.png", dpi=130)

fig, ax = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
edges = np.linspace(-0.05, 0.3, 36)
mid = 0.5 * (edges[1:] + edges[:-1])
for a, (k, s) in zip(ax, (("old", so), ("newSVD", sn))):
    idx = np.digitize(nflh[ok], edges) - 1
    med = np.array([np.median(s[ok][idx == i]) if np.sum(idx == i) > 50 else np.nan for i in range(len(mid))])
    q1 = np.array([np.percentile(s[ok][idx == i], 25) if np.sum(idx == i) > 50 else np.nan for i in range(len(mid))])
    q3 = np.array([np.percentile(s[ok][idx == i], 75) if np.sum(idx == i) > 50 else np.nan for i in range(len(mid))])
    a.plot(mid, med, "o-")
    a.fill_between(mid, q1, q3, alpha=0.25)
    a.set_title(k)
    a.set_xlabel("nFLH [W m-2 um-1 sr-1]")
ax[0].set_ylabel("SIF@678 median / IQR")
fig.tight_layout()
fig.savefig(OUT / "4_sif_vs_nflh.png", dpi=130)

(OUT / "summary.txt").write_text("\n".join(lines) + "\n")
print(f"\nwrote {OUT}/summary.txt and 4 figures")
