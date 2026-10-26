"""Figures and numbers for document/nflh/nFLH_definition.md.

    python document/nflh/make_nflh_figures.py

Reads PACE OCI V3 L1B, L2 AOP V3.2 and the SVD SIF retrieval output (new_svd_retrieval_output_O2A) on kiwi /
~/data, and the SIF basis used by the retrieval (SIF_singular_vector.jld2). Writes four PNGs next to this
script and prints the numbers quoted in the document.

  fig1_nflh_construction.png   Rrs and nLw = Rrs·F0 of one pixel, the V3.2 baseline and line height, and the
                               same spectrum at actual illumination (× cosθs/esd)
  fig2_pure_sif_line_height.png  nFLH of a pure fluorescence spectrum (retrieval SIF basis): V3.2 vs V3.0
  fig3_granule_test.png        manual ΔnLw vs L2 nflh for a whole high-SZA granule; ratio against SZA
  fig4_sif_vs_nflh_boxes.png   retrieved SIF (raw / ×esd/cosθs) vs L2 nflh in four high-latitude 5×5 boxes
"""
from pathlib import Path

import h5py
import matplotlib
import netCDF4
import numpy as np
from scipy.interpolate import CubicSpline

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
K = "/kiwi-data/Data/satellite/PACE_OCI"
RET = "/home/zhe2/data/PACE/new_svd_retrieval_output_O2A"
SIF_LIB = "/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use/SIF_singular_vector.jld2"
U = "W m$^{-2}$ µm$^{-1}$ sr$^{-1}$"

# nFLH definitions (processing_control/input_parameters of the L2 AOP files)
V32 = dict(name="V3.1/V3.2", base=[664.564026, 665.794983, 667.0, 710.5, 711.799988, 713.0], height=680.0)
V30 = dict(name="V3.0", base=[649.599976, 650.900024, 652.099976, 653.299988, 654.599976, 655.799988, 657.099976,
                              658.299988, 659.599976, 710.5, 711.799988, 713.0, 714.299988, 716.799988, 719.200012],
           height=678.20001)
PIXEL = ("20250415T021101", 1580, 427)                 # Southern Ocean, high chl (fig 1)
GRANULE = "20250115T215515"                            # NE Pacific, SZA 58-77° (fig 3)
BOXES = [("North Atlantic, 2025-01-15", "20250115T152200", 796, 688),
         ("Southern Ocean, 2025-01-15", "20250115T181337", 269, 681),
         ("North Atlantic, 2025-04-15", "20250415T123110", 1146, 773),
         ("Southern Ocean, 2025-04-15", "20250415T021101", 1573, 411)]


def dated(sub, tag, g):
    return f"{K}/{sub}/{g[:4]}/{g[4:6]}/{g[6:8]}/PACE_OCI.{g}.{tag}.nc"


def aop_bands(a):
    """Float centre wavelength and F0 of each Rrs band (wavelength_3d is the rounded integer)."""
    sb = a["sensor_band_parameters"]
    w3 = np.ma.filled(sb["wavelength_3d"][:]).astype(int)
    wint = np.ma.filled(sb["wavelength"][:]).astype(int)
    ref = np.array([np.flatnonzero(wint == w)[0] for w in w3])
    return np.ma.filled(sb["fwave"][:]).astype(float)[ref], np.ma.filled(sb["F0"][:]).astype(float)[ref]


def flh_weights(lam, d):
    """Base-band indices, height-band index and LSQ weights so that baseline(height) = w · y[base]."""
    bi = [int(np.argmin(np.abs(lam - w))) for w in d["base"]]
    hi = int(np.argmin(np.abs(lam - d["height"])))
    w = np.array([1.0, d["height"]]) @ np.linalg.pinv(np.column_stack([np.ones(len(bi)), lam[bi]]))
    return bi, hi, w


def l1b_geometry(g):
    with netCDF4.Dataset(dated("L1B_V3", "L1B.V3", g)) as l1:
        return (np.ma.filled(l1["geolocation_data"]["solar_zenith"][:], np.nan).astype(float),
                float(l1.getncattr("earth_sun_distance_correction")))


def sif_basis(lam):
    """First SIF eigenvector (as the retrieval uses it), cubic-interpolated to `lam`, = 1 at 678.2 nm."""
    with h5py.File(SIF_LIB, "r") as f:
        wl, ev = np.array(f["SIF_wavelen"]), np.array(f["SIF_U"])[0, :]   # JLD2 column-major: row = EV
    ev = ev * np.sign(ev[np.argmax(np.abs(ev))])
    s = CubicSpline(wl, ev)
    return s(lam) / s(678.2)


# ---------------------------------------------------------------- fig 1: construction on one pixel
g, sc, px = PIXEL
with netCDF4.Dataset(dated("L2_AOP_V3.2", "L2.OC_AOP.V3_2", g)) as a:
    lam, F0 = aop_bands(a)
    rrs = np.ma.filled(a["geophysical_data"]["Rrs"][sc, px, :], np.nan).astype(float)
    nflh_l2 = float(a["geophysical_data"]["nflh"][sc, px])
    ip = a["processing_control"]["input_parameters"]
    print("AOP attributes:", ip.flh_base_wavelengths, "| height", ip.flh_height_wavelength, "| offset", ip.flh_offset,
          "| brdf_opt", getattr(ip, "brdf_opt", "?"), "| raman_opt", getattr(ip, "raman_opt", "?"))
sza_g, esd = l1b_geometry(g)
sza = float(sza_g[sc, px]); mu = np.cos(np.radians(sza))
nlw = rrs * F0
bi, hi, w = flh_weights(lam, V32)
base = nlw[bi] @ w
print(f"fig1 pixel {g} ({sc},{px}): SZA {sza:.1f}°, esd {esd:.4f}; ΔnLw = {nlw[hi] - base:.4f}, L2 nflh = {nflh_l2:.4f}; "
      f"height band {lam[hi]:.3f} nm; actual-illumination line height {(nlw[hi] - base) * mu / esd:.4f}")

fig, axs = plt.subplots(1, 2, figsize=(14, 4.8))
m = (lam >= 640) & (lam <= 722)
ax = axs[0]
ax.plot(lam[m], rrs[m], "k-", lw=1.4)
ax.plot(lam[bi], rrs[bi], "s", color="tab:blue", mfc="none", label="base bands")
ax.axvline(lam[hi], color="tab:green", ls=":", lw=1, label=f"height band {lam[hi]:.2f} nm")
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Rrs [sr$^{-1}$]"); ax.legend(fontsize=8)
ax.set_title("(a) L2 Rrs (remote-sensing reflectance) of one pixel", fontsize=10, loc="left")
ax = axs[1]
for y, c, lab in ((nlw, "tab:blue", "nLw = Rrs·F0  (normalized: Sun at zenith, 1 AU, no atmosphere)"),
                  (nlw * mu / esd, "tab:orange", f"Rrs·F0·cosθs/esd  (actual illumination, SZA {sza:.0f}°; no t↓)")):
    b = y[bi] @ w
    ax.plot(lam[m], y[m], color=c, lw=1.5, label=lab)
    ax.plot(lam[bi], y[bi], "s", color=c, mfc="none", ms=5)
    xx = np.array([660.0, 718.0]); ax.plot(xx, np.polyval(np.polyfit(lam[bi], y[bi], 1), xx), "--", color=c, lw=1)
    ax.annotate("", xy=(lam[hi], y[hi]), xytext=(lam[hi], b), arrowprops=dict(arrowstyle="<->", color=c, lw=1.6))
    ax.text(lam[hi] + 1.5, 0.5 * (y[hi] + b), f"{y[hi] - b:.3f}", color=c, fontsize=9, va="center")
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel(f"water-leaving radiance [{U}]"); ax.legend(fontsize=7.5, loc="upper right")
ax.set_title(f"(b) nFLH = nLw(679.4) − LSQ baseline(680);  L2 nflh = {nflh_l2:.3f}", fontsize=10, loc="left")
for a_ in axs:
    a_.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"How the L2 nflh is built — granule {g}, scan {sc}, pixel {px}", fontsize=10)
fig.tight_layout(); fig.savefig(HERE / "fig1_nflh_construction.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- fig 2: pure fluorescence line height
with netCDF4.Dataset(dated("L1B_V3", "L1B.V3", g)) as l1:
    red = np.sort(np.ma.filled(l1["sensor_band_parameters"]["red_wavelength"][:], np.nan).astype(float))
lamb = red[(red > 640) & (red < 780)]
sif = sif_basis(lamb)
i678 = int(np.argmin(np.abs(lamb - 678.2)))
fig, axs = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw=dict(width_ratios=[1.5, 1]))
ax = axs[0]
ax.plot(lamb, sif, "tab:green", lw=2, label="fluorescence emission (retrieval SIF basis), = 1 at 678.2 nm")
res = {}
for d, c in ((V32, "tab:blue"), (V30, "tab:purple")):
    bi_, hi_, w_ = flh_weights(lamb, d)
    b_ = sif[bi_] @ w_
    lh = sif[hi_] - b_
    res[d["name"]] = (lh, b_, sif[hi_], bi_, w_)
    ax.plot(lamb[bi_], sif[bi_], "s", color=c, mfc="none", ms=6)
    xx = np.array([min(lamb[bi_]) - 3, max(lamb[bi_]) + 3])
    ax.plot(xx, np.polyval(np.polyfit(lamb[bi_], sif[bi_], 1), xx), "--", color=c, lw=1.2,
            label=f"{d['name']} baseline (base {min(x for x in d['base'] if x < 690):.1f}–{max(x for x in d['base'] if x < 690):.1f} | "
                  f"{min(x for x in d['base'] if x > 690):.1f}–{max(x for x in d['base'] if x > 690):.1f} nm)")
    xo = lamb[hi_] + (0 if d is V32 else -1.5)
    ax.annotate("", xy=(xo, sif[hi_]), xytext=(xo, b_), arrowprops=dict(arrowstyle="<->", color=c, lw=1.6))
    ax.text(xo + (1.5 if d is V32 else -14), 0.5 * (sif[hi_] + b_), f"line height\n{lh:.2f}", color=c, fontsize=9, va="center")
ax.axhline(0, color="0.7", lw=0.6)
ax.set_xlim(640, 730); ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("fluorescence radiance / value at 678.2 nm")
ax.set_title("(a) a pure fluorescence spectrum seen by the nFLH algorithm", fontsize=10, loc="left"); ax.legend(fontsize=7.5)
ax = axs[1]
lh, b_, sh, bi_, w_ = res["V3.1/V3.2"]
contrib = w_ * sif[bi_]
ax.bar(range(6), contrib, color=["tab:blue"] * 3 + ["tab:cyan"] * 3)
ax.set_xticks(range(6), [f"{lamb[i]:.1f}" for i in bi_], rotation=30)
for k, (cv, i) in enumerate(zip(contrib, bi_)):
    ax.text(k, cv + 0.003, f"w={w_[k]:.3f}\nSIF={sif[i]:.2f}", ha="center", fontsize=7)
ax.set_ylim(0, contrib.max() * 1.3)
ax.set_ylabel("w · SIF(base band)  [fraction of SIF(678.2)]"); ax.set_xlabel("V3.2 base band (nm)")
ax.set_title(f"(b) V3.2 baseline = Σ w·SIF = {b_:.3f};  height band = {sh:.3f}\n→ line height = {sh:.3f} − {b_:.3f} = {lh:.3f} × SIF(678.2)",
             fontsize=9, loc="left")
for a_ in axs:
    a_.spines[["top", "right"]].set_visible(False)
fig.tight_layout(); fig.savefig(HERE / "fig2_pure_sif_line_height.png", dpi=150); plt.close(fig)
for k_, (lh, b_, sh, _, _) in res.items():
    print(f"fig2 {k_}: height band {sh:.3f}, baseline {b_:.3f}, line height {lh:.3f} × SIF(678.2); "
          f"peak {lamb[np.argmax(sif)]:.1f} nm = {sif.max():.3f} × SIF(678.2)")

# ---------------------------------------------------------------- fig 3: whole-granule test
with netCDF4.Dataset(dated("L2_AOP_V3.2", "L2.OC_AOP.V3_2", GRANULE)) as a:
    lam, F0 = aop_bands(a)
    bi, hi, w = flh_weights(lam, V32)
    lo, up = min(bi + [hi]), max(bi + [hi]) + 1
    R = np.ma.filled(a["geophysical_data"]["Rrs"][:, :, lo:up], np.nan).astype(np.float32)
    nflh = np.ma.filled(a["geophysical_data"]["nflh"][:], np.nan).astype(float)
    flags = np.ma.filled(a["geophysical_data"]["l2_flags"][:], 0).astype(np.int64) & 0xFFFFFFFF
sza_g, esd_g = l1b_geometry(GRANULE)
b2, h2 = [b - lo for b in bi], hi - lo
ok = np.isfinite(nflh) & np.isfinite(R[..., b2 + [h2]]).all(-1) & np.isfinite(sza_g)
Rv, L2, Z = R[ok], nflh[ok], sza_g[ok]
dnlw = Rv[:, h2] * F0[hi] - (Rv[:, b2] * F0[bi]) @ w
act = dnlw * np.cos(np.radians(Z)) / esd_g
pos = L2 > 1e-3
clip = ~pos
prodwarn = (flags[ok] & 4) != 0
print(f"fig3 {GRANULE}: {ok.sum():,} pixels with L2 nflh, SZA 5/50/95 % {np.percentile(Z, [5, 50, 95]).round(1)}, esd {esd_g:.4f}")
print(f"  ΔnLw / L2 (L2 > 0): median {np.median(dnlw[pos] / L2[pos]):.4f}; median |diff| {np.median(abs(dnlw[pos] - L2[pos])):.4f}, "
      f"max {np.max(abs(dnlw[pos] - L2[pos])):.4f}")
print(f"  actual / L2: median {np.median(act[pos] / L2[pos]):.3f}")
print(f"  L2 at the clip (< 1e-3): {clip.mean():.1%}; of these ΔnLw < 0: {np.mean(dnlw[clip] < 0):.1%}; PRODWARN set: {np.mean(prodwarn[clip]):.1%}")
fig, axs = plt.subplots(1, 3, figsize=(19, 5))
lim = (-0.1, 0.45)
for ax, y, ttl in ((axs[0], dnlw, "(a) manual ΔnLw (unclipped) vs L2 nflh"),
                   (axs[1], act, "(b) actual-illumination line height ΔnLw·cosθs/esd vs L2 nflh")):
    h = ax.hist2d(L2, y, bins=200, range=[lim, lim], norm=mcolors.LogNorm(), cmap="viridis", rasterized=True)
    fig.colorbar(h[3], ax=ax, label="pixels"); ax.plot(lim, lim, "k--", lw=0.8)
    ax.set_xlabel(f"L2 nflh [{U}]"); ax.set_ylabel(f"manual [{U}]"); ax.set_title(ttl, fontsize=10, loc="left")
for zz, c in ((60, "tab:orange"), (70, "tab:red")):
    axs[1].plot(lim, np.array(lim) * np.cos(np.radians(zz)) / esd_g, ":", color=c, label=f"y = x·cos({zz}°)/esd")
axs[1].legend(fontsize=8)
ax = axs[2]
edges = np.arange(58, 78, 1.0); cc = 0.5 * (edges[:-1] + edges[1:])
for y, col, lab in ((dnlw, "tab:blue", "ΔnLw / L2 nflh"), (act, "tab:orange", "ΔnLw·cosθs/esd / L2 nflh")):
    med = [np.median((y / L2)[pos & (Z >= e0) & (Z < e1)]) if np.sum(pos & (Z >= e0) & (Z < e1)) > 200 else np.nan
           for e0, e1 in zip(edges[:-1], edges[1:])]
    ax.plot(cc, med, "o-", ms=4, color=col, label=lab)
ax.plot(cc, np.cos(np.radians(cc)) / esd_g, "k:", label="cos(SZA)/esd")
ax.set_ylim(0, 1.15); ax.set_xlabel("SZA (°)"); ax.set_ylabel("median ratio (L2 > 0)"); ax.legend(fontsize=8)
ax.set_title("(c) the L2 nflh carries no cosθs or esd", fontsize=10, loc="left")
for a_ in axs:
    a_.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"Granule {GRANULE} (NE Pacific, 36–57°N), every pixel with an L2 nflh: V3.2 definition recomputed from Rrs·F0", fontsize=10)
fig.tight_layout(); fig.savefig(HERE / "fig3_granule_test.png", dpi=150); plt.close(fig)

# ---------------------------------------------------------------- fig 4: SIF vs nFLH in four boxes
LH = res["V3.1/V3.2"][0]
fig, axs = plt.subplots(1, 4, figsize=(21, 5))
rows = []
for ax, (lab, g, sc, px) in zip(axs, BOXES):
    s0, s1, p0, p1 = sc - 2, sc + 3, px - 2, px + 3
    with netCDF4.Dataset(f"{RET}/{g[:4]}/{g[4:6]}/{g[6:8]}/interim_{g}_svd_retrieval.nc") as r:
        sif = np.ma.filled(r["sif_radiance_678nm"][s0:s1, p0:p1], np.nan).astype(float)
        good = ((np.ma.filled(r["converged"][s0:s1, p0:p1], 0) == 1) & (np.ma.filled(r["is_dark"][s0:s1, p0:p1], 0) == 1)
                & (np.ma.filled(r["reduced_chi2"][s0:s1, p0:p1], np.inf) < 1.5))
    with netCDF4.Dataset(dated("L2_AOP_V3.2", "L2.OC_AOP.V3_2", g)) as a:
        nf = np.ma.filled(a["geophysical_data"]["nflh"][s0:s1, p0:p1], np.nan).astype(float)
        fl = np.ma.filled(a["geophysical_data"]["l2_flags"][s0:s1, p0:p1], 1).astype(np.int64) & 0xFFFFFFFF
    sz, es = l1b_geometry(g)
    mu = np.cos(np.radians(sz[s0:s1, p0:p1]))
    k = good & (fl == 0) & np.isfinite(nf) & np.isfinite(sif)
    x, y, yn, xa = nf[k], sif[k], sif[k] * es / mu[k], nf[k] * mu[k] / es
    ax.plot(x, y, "o", color="tab:green", ms=5, label=f"raw SIF (median {np.median(y):.3f})")
    ax.plot(x, yn, "^", color="tab:green", mfc="none", ms=6, label=f"SIF·esd/cosθs (median {np.median(yn):.3f})")
    ax.plot(x, xa, "o", color="tab:orange", ms=4, label=f"nFLH·cosθs/esd (median {np.median(xa):.3f})")
    top = max(yn.max(), x.max()) * 1.1; xx = np.array([0, top])
    ax.plot(xx, xx, "k--", lw=0.8, label="1:1"); ax.plot(xx, xx / LH, "k:", lw=1, label=f"nFLH / {LH:.2f}")
    ax.set_xlim(0, top); ax.set_ylim(0, top); ax.set_xlabel(f"L2 nflh [{U}]"); ax.set_ylabel(f"[{U}]")
    ax.set_title(f"{lab}\n{k.sum()} pixels, SZA {np.degrees(np.arccos(np.median(mu[k]))):.0f}°", fontsize=9, loc="left")
    ax.legend(fontsize=7, loc="upper left"); ax.spines[["top", "right"]].set_visible(False)
    rows.append((lab, k.sum(), np.median(y), np.median(x), np.median(yn), np.median(x) / LH, np.median(xa)))
fig.suptitle("Retrieved SIF vs L2 nflh in four randomly drawn high-latitude, high-chl 5×5 boxes "
             "(pixels: l2_flags == 0, dark, converged, χ² < 1.5)", fontsize=10)
fig.tight_layout(); fig.savefig(HERE / "fig4_sif_vs_nflh_boxes.png", dpi=150); plt.close(fig)
print("fig4 (box, n, median raw SIF, nflh, SIF·esd/cosθs, nflh/LH, nflh·cosθs/esd):")
for r_ in rows:
    print("  " + r_[0] + ": " + ", ".join(f"{v:.3f}" if isinstance(v, float) else str(v) for v in r_[1:]))
print("wrote fig1-4 in", HERE)
