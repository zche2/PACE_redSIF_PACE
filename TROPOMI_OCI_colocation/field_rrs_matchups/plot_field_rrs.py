"""Locations and Rrs spectra of the field / PACE OCI matchups in fullmatchups_Rrs_GUI.csv.

    python TROPOMI_OCI_colocation/field_rrs_matchups/plot_field_rrs.py [path/to/fullmatchups_Rrs_GUI.csv]

The CSV (gitignored) defaults to the main checkout's TROPOMI_OCI_colocation/ folder.
Field and PACE Rrs [sr^-1] are given at 315-715 nm (5 nm steps plus a few extra bands), so the
red window is 600-715 nm: the data stop at 715 nm, below the requested 800 nm.

Outputs (out/):
  1_field_locations.png   map of field locations, coloured by platform type
  2_rrs_spectra.png       field and PACE Rrs, full range and 600-715 nm, one line per matchup
  3_red_rrs_by_type.png   600-715 nm field Rrs per platform type (median and 10-90 % band)
  summary.txt             counts per platform / experiment and spectral coverage
"""
import re
import sys
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    Path("/home/zhe2/FraLab/PACE_redSIF_PACE/TROPOMI_OCI_colocation/fullmatchups_Rrs_GUI.csv")
OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)
RED = (600.0, 800.0)
lines = []


def say(s=""):
    print(s)
    lines.append(s)


def spectrum(df, prefix):
    """Wavelengths and the (n_matchup, n_lambda) Rrs matrix for columns like '<prefix>_675'."""
    cols = [(float(m.group(1)), c) for c in df.columns if (m := re.fullmatch(rf"{prefix}_(\d+(?:\.\d+)?)", c))]
    cols.sort()
    lam = np.array([w for w, _ in cols])
    return lam, df[[c for _, c in cols]].to_numpy(float)


def draw(ax, lam, r, sel=None, **kw):
    """Plot one spectrum through its valid points only (NaN gaps would otherwise hide sparse,
    e.g. multispectral aircraft, spectra); sparse spectra also get markers."""
    f = np.isfinite(r) & (True if sel is None else sel)
    if f.sum() < 1:
        return
    sparse = np.isfinite(r).sum() < 30
    ax.plot(lam[f], r[f], marker="o" if sparse else None, ms=2.2 if sparse else 0, **kw)


df = pd.read_csv(CSV, low_memory=False)
lam_f, R_f = spectrum(df, "FIELD_Rrs")
lam_p, R_p = spectrum(df, "PACE_Rrs")
types = df["data_type"].fillna("unknown")
type_names = types.value_counts().index.tolist()
cmap = plt.get_cmap("tab10")
col = {t: cmap(i) for i, t in enumerate(type_names)}
red_f = (lam_f >= RED[0]) & (lam_f <= RED[1])
red_p = (lam_p >= RED[0]) & (lam_p <= RED[1])

# ------------------------------------------------------------------ summary
say(f"{CSV.name}: {len(df)} matchups, {df['FIELD_experiment'].nunique()} experiments, "
    f"{df['FIELD_cruise'].nunique()} cruises/deployments")
say(f"dates {df['FIELD_datetime'].min()} -> {df['FIELD_datetime'].max()}")
say(f"field Rrs: {lam_f.min():.0f}-{lam_f.max():.0f} nm ({lam_f.size} bands); "
    f"PACE Rrs: {lam_p.min():.0f}-{lam_p.max():.0f} nm ({lam_p.size} bands)")
say(f"red window {RED[0]:.0f}-{RED[1]:.0f} nm: data available only up to {lam_f[red_f].max():.0f} nm "
    f"({red_f.sum()} field bands)")
has_red = np.isfinite(R_f[:, red_f]).any(axis=1)
reach = np.array([lam_f[np.isfinite(r)].max() if np.isfinite(r).any() else np.nan for r in R_f])
say(f"matchups with any field Rrs in 600-715 nm: {has_red.sum()} / {len(df)};  "
    f"reaching >= 700 nm: {(reach >= 700).sum()};  reaching 715 nm: {(reach >= 715).sum()}")
say("\nby platform type (n, n with red data, longest field wavelength):")
for t in type_names:
    m = (types == t).to_numpy()
    say(f"  {t:16s}{m.sum():4d}{has_red[m].sum():5d}   max λ {np.nanmax(reach[m]):.0f} nm")
say("\nby experiment:")
for e, n in df["FIELD_experiment"].value_counts().items():
    m = (df["FIELD_experiment"] == e).to_numpy()
    say(f"  {e:24s}{n:4d}   platform: {', '.join(sorted(types[m].unique()))}")

# --------------------------------------------------------------------- map
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson(central_longitude=-60))
ax.set_global()
ax.add_feature(cfeature.LAND, facecolor="0.88")
ax.coastlines(linewidth=0.5)
ax.gridlines(linewidth=0.3, alpha=0.5, linestyle=":")
for t in type_names:
    m = (types == t).to_numpy()
    ax.scatter(df["FIELD_lon"][m], df["FIELD_lat"][m], s=30, color=col[t], edgecolor="k", linewidth=0.4,
               transform=ccrs.PlateCarree(), label=f"{t} (n={m.sum()})", zorder=3)
ax.legend(loc="lower left", fontsize=8, frameon=True)
ax.set_title(f"Field Rrs matchups with PACE OCI ({len(df)} matchups, {df['FIELD_experiment'].nunique()} experiments)",
             fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "1_field_locations.png", dpi=150)
plt.close(fig)

# ----------------------------------------------------------------- spectra
fig, axs = plt.subplots(2, 2, figsize=(13, 8.5), sharey="row")
for j, (lam, R, src, red) in enumerate(((lam_f, R_f, "Field", red_f), (lam_p, R_p, "PACE OCI", red_p))):
    for t in type_names:
        for i in np.where((types == t).to_numpy())[0]:
            draw(axs[0, j], lam, R[i], color=col[t], lw=0.6, alpha=0.5)
            draw(axs[1, j], lam, R[i], sel=red, color=col[t], lw=0.7, alpha=0.6)
    axs[0, j].set_title(f"{src} Rrs, all wavelengths ({lam.min():.0f}-{lam.max():.0f} nm)", fontsize=10)
    axs[1, j].set_title(f"{src} Rrs, {RED[0]:.0f}-{lam[red].max():.0f} nm (no data beyond {lam.max():.0f} nm)",
                        fontsize=10)
    axs[1, j].set_xlabel("Wavelength (nm)")
    axs[1, j].axvspan(670, 700, color="0.9", zorder=0)
for a in axs[:, 0]:
    a.set_ylabel("Rrs (sr$^{-1}$)")
for a in axs.ravel():
    a.axhline(0, color="k", lw=0.4)
    a.spines[["top", "right"]].set_visible(False)
axs[1, 1].text(0.02, 0.95, "grey band: 670-700 nm (chl fluorescence)", transform=axs[1, 1].transAxes,
               fontsize=8, va="top", color="0.35")
handles = [plt.Line2D([], [], color=col[t], lw=2, label=f"{t} (n={(types == t).sum()})") for t in type_names]
fig.legend(handles=handles, loc="upper center", ncol=len(type_names), fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(OUT / "2_rrs_spectra.png", dpi=150)
plt.close(fig)

# --------------------------------------------- red window per platform type
nt = len(type_names)
ncol = 3
nrow = int(np.ceil(nt / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 3.2 * nrow), sharex=True, squeeze=False)
for k, t in enumerate(type_names):
    a = axs.ravel()[k]
    m = (types == t).to_numpy()
    Rr = R_f[m][:, red_f]
    for i in np.where(m)[0]:
        draw(a, lam_f, R_f[i], sel=red_f, color=col[t], lw=0.5, alpha=0.35)
    n_ok = np.isfinite(Rr).sum(axis=0)
    b = n_ok > 0                      # bands with any data for this platform type
    if b.any():
        med = np.array([np.median(Rr[np.isfinite(Rr[:, q]), q]) if b[q] else np.nan for q in range(Rr.shape[1])])
        lo = np.array([np.percentile(Rr[np.isfinite(Rr[:, q]), q], 10) if b[q] else np.nan for q in range(Rr.shape[1])])
        hi = np.array([np.percentile(Rr[np.isfinite(Rr[:, q]), q], 90) if b[q] else np.nan for q in range(Rr.shape[1])])
        a.fill_between(lam_f[red_f][b], lo[b], hi[b], color=col[t], alpha=0.2)
        a.plot(lam_f[red_f][b], med[b], color="k", lw=1.4, marker="o" if b.sum() < 10 else None, ms=3,
               label="median")
    a.axvspan(670, 700, color="0.92", zorder=0)
    a.axhline(0, color="k", lw=0.4)
    a.set_title(f"{t}: {m.sum()} matchups, {has_red[m].sum()} with red data", fontsize=9)
    a.spines[["top", "right"]].set_visible(False)
for a in axs.ravel()[nt:]:
    a.set_axis_off()
for a in axs[-1]:
    a.set_xlabel("Wavelength (nm)")
for a in axs[:, 0]:
    a.set_ylabel("Field Rrs (sr$^{-1}$)")
fig.suptitle("Field Rrs in the red (600-715 nm): lines = matchups, black = median, band = 10-90 %", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "3_red_rrs_by_type.png", dpi=150)
plt.close(fig)

(OUT / "summary.txt").write_text("\n".join(lines) + "\n")
print(f"\nwrote {OUT}/1_field_locations.png, 2_rrs_spectra.png, 3_red_rrs_by_type.png, summary.txt")
