"""Refit the V3.2 L2 nFLH (680 nm definition) from the 5x5-box Rrs and F0, and compare it with the stored
L2 nflh in matchup_boxes_5x5.csv.

    python TROPOMI_OCI_colocation/field_rrs_matchups/fit_nflh680.py [any V3.2 L2 OC_AOP file]

V3.2 processing_control/input_parameters:
    flh_offset = 0, flh_base_wavelengths = 664.564026, 665.794983, 667.0, 710.5, 711.799988, 713.0,
    flh_height_wavelength = 680.0
l2gen convention (reproduces the stored values):
    nLw = Rrs * F0                                     [W m-2 um-1 sr-1]
    baseline = least-squares line through nLw at the bands nearest the base wavelengths, evaluated at 680.0
    nFLH = max(nLw(band nearest 680.0, i.e. 679.448) - baseline - flh_offset, 0)
The clip at 0 matches the stored nflh (packed valid_min); the unclipped value is kept as well.
F0 and the float band centres come from the AOP file (identical across granules).

Outputs: out/8_nflh680_fit_vs_L2.png, nflh680_fit.csv (per valid pixel: ids, L2 nflh, fit, unclipped fit)
"""
import sys
from pathlib import Path

import matplotlib
import netCDF4
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
AOP = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/kiwi-data/Data/satellite/PACE_OCI/L2_AOP_V3.2/2025/07/08/PACE_OCI.20250708T204612.L2.OC_AOP.V3_2.nc")
BASE = [664.564026, 665.794983, 667.0, 710.5, 711.799988, 713.0]
HEIGHT = 680.0
OFFSET = 0.0

# ------------------------------------------------------------ bands and F0
with netCDF4.Dataset(AOP) as d:
    sb = d["sensor_band_parameters"]
    w3 = np.ma.filled(sb["wavelength_3d"][:]).astype(int)            # Rrs bands (rounded nm)
    wint = np.ma.filled(sb["wavelength"][:]).astype(int)
    ref = np.array([np.flatnonzero(wint == w)[0] for w in w3])
    fwave = np.ma.filled(sb["fwave"][:]).astype(float)[ref]
    F0 = np.ma.filled(sb["F0"][:]).astype(float)[ref]
    ip = d["processing_control"]["input_parameters"]
    print("AOP file:", AOP.name)
    print("  flh_base_wavelengths  =", ip.flh_base_wavelengths)
    print("  flh_height_wavelength =", ip.flh_height_wavelength, "  flh_offset =", ip.flh_offset)

bi = [int(np.argmin(np.abs(fwave - w))) for w in BASE]
hi = int(np.argmin(np.abs(fwave - HEIGHT)))
A = np.column_stack([np.ones(len(bi)), fwave[bi]])
wts = np.array([1.0, HEIGHT]) @ np.linalg.pinv(A)                    # baseline(680) = wts . nLw(base)
print(f"\nheight {HEIGHT} -> band {fwave[hi]:.3f} nm (F0 {F0[hi]:.2f})")
for w, i, k in zip(BASE, bi, wts):
    print(f"  base {w:8.3f} -> band {fwave[i]:8.3f} nm  F0 {F0[i]:8.2f}  weight {k:+.4f}")

# ------------------------------------------------------------ fit on the box pixels
box = pd.read_csv(HERE / "matchup_boxes_5x5.csv")
cols = [f"Rrs_{w}" for w in w3]
px = box[box.valid & box.nflh.notna()].copy()
nLw = px[cols].to_numpy(float) * F0
px["nflh680_fit_unclipped"] = nLw[:, hi] - nLw[:, bi] @ wts - OFFSET
px["nflh680_fit"] = px.nflh680_fit_unclipped.clip(lower=0)          # as l2gen stores it
d_ = px.nflh680_fit - px.nflh
clip = px.nflh680_fit_unclipped <= 0
px[["matchup_id", "granule", "scan", "pixel", "is_centre", "nflh", "nflh680_fit", "nflh680_fit_unclipped"]].to_csv(
    HERE / "nflh680_fit.csv", index=False)

print(f"\n{len(px)} valid pixels in {px.matchup_id.nunique()} matchups; fit clipped at 0 for {clip.sum()}")
print(f"fit vs L2: median |diff| {d_.abs().median():.2e}, 99th pct {d_.abs().quantile(.99):.2e}, max {d_.abs().max():.2e}; "
      f"|diff| < 0.005 for {(d_.abs() < 0.005).mean() * 100:.2f} %; "
      f"median ratio (unclipped) {(px.nflh680_fit / px.nflh)[~clip].median():.5f}; "
      f"r = {np.corrcoef(px.nflh680_fit, px.nflh)[0, 1]:.6f}")
print(f"clipped pixels: L2 nflh range [{px.nflh[clip].min():.2e}, {px.nflh[clip].max():.2e}], unclipped fit range "
      f"[{px.nflh680_fit_unclipped[clip].min():.3f}, {px.nflh680_fit_unclipped[clip].max():.3f}]")

# per-matchup means of the valid pixels, as in matchup_summary.csv
mm = px.groupby("matchup_id")[["nflh", "nflh680_fit"]].mean()

# ------------------------------------------------------------ figure
fig, axs = plt.subplots(1, 3, figsize=(16, 5))
ax = axs[0]
lim = (-0.02, max(px.nflh.max(), px.nflh680_fit.max()) * 1.05)
ax.scatter(px.nflh[~clip], px.nflh680_fit[~clip], s=5, alpha=0.5, label=f"fit > 0 ({(~clip).sum()})")
ax.scatter(px.nflh[clip], px.nflh680_fit[clip], s=12, alpha=0.7, color="tab:red",
           label=f"fit clipped at 0 ({clip.sum()})")
ax.plot(lim, lim, "k--", lw=0.8)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel("L2 nflh (V3.2, stored)  [W m$^{-2}$ µm$^{-1}$ sr$^{-1}$]")
ax.set_ylabel("refit nFLH$_{680}$ = max(nLw(679.4) − LSQ baseline(680), 0)")
ax.set_title("valid box pixels", fontsize=10)
ax.legend(fontsize=8)

ax = axs[1]
ax.hist(d_ * 1e3, bins=60, color="tab:blue")
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("refit − L2 (10$^{-3}$ W m$^{-2}$ µm$^{-1}$ sr$^{-1}$)")
ax.set_ylabel("pixels")
ax.set_title(f"median |diff| {d_.abs().median() * 1e3:.2f}e-3, max {d_.abs().max() * 1e3:.2f}e-3\n"
             f"(nLw packing ≈ Rrs step × F0)", fontsize=10)

ax = axs[2]
ax.scatter(mm.nflh, mm.nflh680_fit, s=14)
l2 = (-0.02, max(mm.nflh.max(), mm.nflh680_fit.max()) * 1.05)
ax.plot(l2, l2, "k--", lw=0.8)
ax.set_xlim(l2)
ax.set_ylim(l2)
ax.set_xlabel("L2 nflh, mean of valid pixels")
ax.set_ylabel("refit nFLH$_{680}$, mean of valid pixels")
ax.set_title(f"per matchup ({len(mm)})", fontsize=10)
for a in axs:
    a.spines[["top", "right"]].set_visible(False)
fig.suptitle("nFLH 680 refit from Rrs·F0 (base 664.6, 665.8, 667.0 | 710.6, 711.8, 713.1 nm; height band 679.4 nm, "
             "baseline at 680.0 nm) vs stored V3.2 L2 nflh", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "8_nflh680_fit_vs_L2.png", dpi=150)
print("wrote out/8_nflh680_fit_vs_L2.png, nflh680_fit.csv")
