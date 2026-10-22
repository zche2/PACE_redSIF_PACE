"""Whole-granule comparison of every available O2A SVD retrieval of one PACE granule.

    python retrieval_comparison/compare_all_runs.py [YYYYMMDDTHHMMSS]

Runs (those whose output is complete are included automatically):
  npoly2 / npoly3 / npoly5   Aug01 transmittance SVD, 2 / 3 / 5 Legendre terms
  newSVD_npoly2 / _npoly3    scattering-included SVD, 2 / 3 Legendre terms (different PC basis)

Pixel selection: ocean & L2 flags ok (no cloud/glint/straylight/...) & chl valid.
`is_dark` = max radiance <= dark_max_radiance, i.e. dark = clear open water,
bright = cloud / glint / swath edge; both are reported.
Outputs: out/all_runs/summary.txt and out/all_runs/F_all_runs_summary.png
"""
import sys
from pathlib import Path

import matplotlib
import netCDF4 as nc
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

STAMP = sys.argv[1] if len(sys.argv) > 1 else "20250724T170203"
DATE = f"{STAMP[:4]}/{STAMP[4:6]}/{STAMP[6:8]}"
PACE = Path("/home/zhe2/data/PACE")
KIWI = Path("/kiwi-data/Data/satellite/PACE_OCI")
RUNS_ALL = [  # key, folder
    ("npoly2", "new_svd_retrieval_output_npoly2"),
    ("npoly3", "new_svd_retrieval_output_O2A"),
    ("npoly5", "new_svd_retrieval_output_npoly5"),
    ("newSVD_npoly2", "new_svd_retrieval_output_O2A_newSVD_npoly2"),
    ("newSVD_npoly3", "new_svd_retrieval_output_O2A_newSVD"),
]
COLORS = {"npoly2": "#e7298a", "npoly3": "#0B7A75", "npoly5": "#d95f02", "newSVD_npoly2": "#7570b3", "newSVD_npoly3": "#a6761d"}
OUT = Path(__file__).resolve().parent / "out" / "all_runs"
OUT.mkdir(parents=True, exist_ok=True)
BAD_FLAGS = ["ATMFAIL", "LAND", "HIGLINT", "HILT", "HISATZEN", "STRAYLIGHT", "CLDICE", "HISOLZEN", "CHLFAIL"]
STATUS = {1: "converged", 2: "LM stalled", 0: "max-iter", 4: "failed"}
lines = []


def say(s=""):
    print(s)
    lines.append(s)


def path_of(folder):
    return PACE / folder / DATE / f"interim_{STAMP}_svd_retrieval.nc"


def ready(folder):
    f = path_of(folder)
    try:
        with nc.Dataset(f) as d:
            return f.stat().st_size > 3e8 and d["x_hat"].shape[0] > 0
    except Exception:
        return False


def load(folder):
    ds = nc.Dataset(path_of(folder))
    ds.set_auto_mask(False)
    r = {k: ds[k][:] for k in ["converged", "status_code", "n_steps", "reduced_chi2", "dof", "rmse",
                               "sif_radiance_678nm", "sif_ev1", "is_dark", "is_ocean"]}
    var = ds["S_posterior_diag"][-1]
    ds.close()
    ev1, sif = r["sif_ev1"].astype(float), r["sif_radiance_678nm"].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        w = np.nanmedian((sif / ev1)[np.abs(ev1) > 1e-6])
    r["sif"], r["sigma"] = sif, np.abs(w) * np.sqrt(np.maximum(var, 0.0))
    return r


RUNS = [(k, f) for k, f in RUNS_ALL if ready(f)]
for k, f in RUNS_ALL:
    if (k, f) not in RUNS:
        print(f"[skip] {k}: output not complete / not readable yet")
keys = [k for k, _ in RUNS]
R = {k: load(f) for k, f in RUNS}

v = "V3_2"
with nc.Dataset(KIWI / "L2_BGC_V3.2" / DATE / f"PACE_OCI.{STAMP}.L2.OC_BGC.{v}.nc") as ds:
    g = ds["geophysical_data"]
    g.set_auto_mask(True)
    chl = np.ma.filled(g["chlor_a"][:].astype(float), np.nan)
    fl = g["l2_flags"]
    flags = np.ma.filled(fl[:], 0).astype(np.int64)
    names, masks = fl.getncattr("flag_meanings").split(), np.atleast_1d(fl.getncattr("flag_masks"))
with nc.Dataset(KIWI / "L2_AOP_V3.2" / DATE / f"PACE_OCI.{STAMP}.L2.OC_AOP.{v}.nc") as ds:
    g = ds["geophysical_data"]
    g.set_auto_mask(True)
    nflh = np.ma.filled(g["nflh"][:].astype(float), np.nan)
bad = sum(int(masks[names.index(n)]) for n in BAD_FLAGS if n in names)
chl[chl <= 0] = np.nan

ref = R[keys[0]]
ocean, dark = ref["is_ocean"] == 1, ref["is_dark"] == 1
for k in keys:
    assert np.array_equal(R[k]["is_ocean"], ref["is_ocean"]) and np.array_equal(R[k]["is_dark"], ref["is_dark"])
clear = ocean & ((flags & bad) == 0) & np.isfinite(chl)
say(f"Granule {STAMP}   runs: {keys}   L2 V3.2")
say(f"ocean {ocean.sum()} (dark {(ocean & dark).sum()}, bright {(ocean & ~dark).sum()});  "
    f"clear {clear.sum()} (dark {(clear & dark).sum()}, bright {(clear & ~dark).sum()})")
H = f"{'':34s}" + "".join(f"{k:>14s}" for k in keys)


def row(label, fn, fmt="{:14.3f}"):
    say(f"{label:34s}" + "".join(fmt.format(fn(R[k])) for k in keys))


# ----------------------------------------------------------- 1. convergence
say("\n=== 1. CONVERGENCE ===")
say(H)
for lab, sel in [("clear ocean", clear), ("clear & dark (open water)", clear & dark),
                 ("clear & bright", clear & ~dark), ("all ocean (incl. cloud)", ocean)]:
    say(f"-- {lab}  (n={sel.sum()})")
    row("  converged %", lambda r: 100 * np.mean(r["converged"][sel] == 1), "{:13.1f}%")
    row("  LM stalled %", lambda r: 100 * np.mean(r["status_code"][sel] == 2), "{:13.1f}%")
    row("  max-iter (not converged) %", lambda r: 100 * np.mean(r["status_code"][sel] == 0), "{:13.1f}%")
    row("  reduced chi2 median", lambda r: np.nanmedian(r["reduced_chi2"][sel]))
    row("  reduced chi2 p90", lambda r: np.nanpercentile(r["reduced_chi2"][sel], 90))
    row("  n_steps median", lambda r: np.nanmedian(r["n_steps"][sel]), "{:14.0f}")
    row("  DOF signal median", lambda r: np.nanmedian(r["dof"][sel]), "{:14.1f}")
allc = clear.copy()
for k in keys:
    allc &= R[k]["converged"] == 1
say(f"converged in ALL runs: {allc.sum()} of {clear.sum()} clear pixels")
row("  reduced chi2 median there", lambda r: np.nanmedian(r["reduced_chi2"][allc]))

# ------------------------------------------------------------- 2. sigma / SNR
say("\n=== 2. POSTERIOR sigma(SIF@678) and |SIF|/sigma  (clear ocean) ===")
say(H)
row("sigma median", lambda r: np.nanmedian(r["sigma"][clear]), "{:14.4f}")
row("sigma p95", lambda r: np.nanpercentile(r["sigma"][clear], 95), "{:14.4f}")
row("|SIF|/sigma median", lambda r: np.nanmedian(np.abs(r["sif"][clear]) / r["sigma"][clear]))
r0 = R["npoly3"] if "npoly3" in R else R[keys[0]]
for k in keys:
    if R[k] is not r0:
        ra = R[k]["sigma"][clear] / r0["sigma"][clear]
        say(f"  sigma {k}/{[kk for kk in keys if R[kk] is r0][0]}: median {np.nanmedian(ra):.3f}")

# ------------------------------------------------------------- 3. SIF magnitude
say(f"\n=== 3. SIF@678 MAGNITUDE (converged in all runs, n={allc.sum()}) ===")
say(H)
for q in (5, 25, 50, 75, 95):
    row(f"  p{q}", lambda r: np.nanpercentile(r["sif"][allc], q), "{:14.4f}")
row("  mean", lambda r: np.nanmean(r["sif"][allc]), "{:14.4f}")
row("  std", lambda r: np.nanstd(r["sif"][allc]), "{:14.4f}")
row("  fraction SIF<0 %", lambda r: 100 * np.mean(r["sif"][allc] < 0), "{:13.1f}%")
refk = "npoly3" if "npoly3" in R else keys[0]
say(f"vs {refk} (pixel-wise):")
for k in keys:
    if k == refk:
        continue
    d = R[k]["sif"][allc] - R[refk]["sif"][allc]
    sl, ic = np.polyfit(R[refk]["sif"][allc], R[k]["sif"][allc], 1)
    say(f"  {k:7s} mean diff {d.mean():+.4f}  median {np.median(d):+.4f}  rms {np.sqrt(np.mean(d**2)):.4f}  "
        f"r {np.corrcoef(R[refk]['sif'][allc], R[k]['sif'][allc])[0, 1]:.3f}  fit = {sl:.3f}*{refk} {ic:+.4f}")

# ------------------------------------------------------------- 4. low chl
say("\n=== 4. SIF BY CHL BIN (converged in all runs; median SIF, % negative) ===")
say(f"{'chl bin (mg m-3)':20s}{'n':>8s}{'nFLH med':>10s}" + "".join(f"{k:>18s}" for k in keys))
for lo, hi in [(0, 0.05), (0.05, 0.1), (0.1, 0.3), (0.3, 1), (1, 100)]:
    m = allc & (chl >= lo) & (chl < hi)
    if m.sum() < 50:
        continue
    say(f"{lo:5.2f}-{hi:<14.2f}{m.sum():8d}{np.nanmedian(nflh[m]):10.4f}" + "".join(
        f"{np.median(R[k]['sif'][m]):11.4f} /{100 * np.mean(R[k]['sif'][m] < 0):4.1f}%" for k in keys))

# ------------------------------------------------------------- 5. nFLH
say("\n=== 5. SIF vs nFLH (converged in all runs, nFLH valid) ===")
ok = allc & np.isfinite(nflh)
say(f"{'':26s}" + "".join(f"{k:>14s}" for k in keys))
for lab, sel in [("all chl", ok), ("chl<0.1", ok & (chl < 0.1)), ("0.1<=chl<0.3", ok & (chl >= 0.1) & (chl < 0.3)),
                 ("chl>=0.3", ok & (chl >= 0.3))]:
    say(f"r     {lab:20s}" + "".join(f"{np.corrcoef(nflh[sel], R[k]['sif'][sel])[0, 1]:14.3f}" for k in keys))
    say(f"slope {lab:20s}" + "".join(f"{np.polyfit(nflh[sel], R[k]['sif'][sel], 1)[0]:14.3f}" for k in keys))
say(f"{'intercept (all chl)':26s}" + "".join(f"{np.polyfit(nflh[ok], R[k]['sif'][ok], 1)[1]:14.4f}" for k in keys))

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(2, 2, figsize=(12, 9))
sel = clear & dark
w = 0.8 / len(keys)
for i, k in enumerate(keys):
    for j, (code, name) in enumerate(STATUS.items()):
        pass
    fr = [100 * np.mean(R[k]["status_code"][sel] == c) for c in STATUS]
    ax[0, 0].bar(np.arange(len(fr)) + (i - (len(keys) - 1) / 2) * w, fr, w, color=COLORS[k], label=k)
    x = np.sort(np.clip(R[k]["reduced_chi2"][sel], 0, 10))
    ax[0, 1].plot(x, np.arange(1, x.size + 1) / x.size, color=COLORS[k], label=k)
ax[0, 0].set_xticks(range(len(STATUS)), list(STATUS.values()))
ax[0, 0].set_ylabel("% of clear dark pixels")
ax[0, 0].legend()
ax[0, 1].set_xlabel("reduced chi2 (clipped at 10)")
ax[0, 1].set_ylabel("CDF (clear dark pixels)")
edges = np.linspace(-0.02, 0.3, 33)
mid = 0.5 * (edges[1:] + edges[:-1])
idx = np.digitize(nflh[ok], edges) - 1
for k in keys:
    med = [np.median(R[k]["sif"][ok][idx == b]) if np.sum(idx == b) > 100 else np.nan for b in range(len(mid))]
    ax[1, 0].plot(mid, med, "o-", ms=3, color=COLORS[k], label=k)
ax[1, 0].plot(mid, mid, "k:", lw=0.8, label="1:1")
ax[1, 0].set_xlabel("L2 nFLH")
ax[1, 0].set_ylabel("median SIF@678")
ax[1, 0].legend(fontsize=8)
lc = allc & (chl < 0.1)
for k in keys:
    ax[1, 1].hist(R[k]["sif"][lc], 120, range=(-0.15, 0.25), histtype="step", color=COLORS[k], label=f"{k}")
ax[1, 1].axvline(0, color="k", lw=0.5)
ax[1, 1].set_xlabel("SIF@678, chl<0.1 (converged in all runs)")
ax[1, 1].legend(fontsize=8)
fig.suptitle(f"PACE {STAMP} O2A retrievals: {', '.join(keys)}")
fig.tight_layout()
fig.savefig(OUT / "F_all_runs_summary.png", dpi=130)
(OUT / "summary.txt").write_text("\n".join(lines) + "\n")
print("wrote", OUT)
