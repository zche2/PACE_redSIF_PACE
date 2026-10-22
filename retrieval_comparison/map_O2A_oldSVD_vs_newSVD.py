"""Swath maps for the old vs newSVD O2A retrievals (same granule as the compare script).

    python retrieval_comparison/map_O2A_oldSVD_vs_newSVD.py [YYYYMMDDTHHMMSS]

Writes out/5_maps_clear_sky.png (SIF, difference, nFLH, chl, convergence on
clear-sky ocean), out/6_maps_convergence_all.png (convergence on all lit
ocean incl. cloud) and out/7_maps_all_pixels.png (no quality filtering: every
ocean pixel, incl. cloud, swath edge and dark).
"""
import sys
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm

STAMP = sys.argv[1] if len(sys.argv) > 1 else "20250724T170203"
DATE = f"{STAMP[:4]}/{STAMP[4:6]}/{STAMP[6:8]}"
ROOT = Path("/home/zhe2/data/PACE")
L2_ROOT = Path("/kiwi-data/Data/satellite/PACE_OCI")
L2_VERSION = "V3.2"
OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)
BAD_FLAGS = ["ATMFAIL", "LAND", "HIGLINT", "HILT", "HISATZEN", "STRAYLIGHT",
             "CLDICE", "HISOLZEN", "CHLFAIL"]


def run(d):
    ds = nc.Dataset(ROOT / d / DATE / f"interim_{STAMP}_svd_retrieval.nc")
    ds.set_auto_mask(False)
    r = {k: ds[k][:] for k in ["latitude", "longitude", "converged", "sif_radiance_678nm", "is_dark", "is_ocean"]}
    ds.close()
    return r


O, N = run("new_svd_retrieval_output_O2A"), run("new_svd_retrieval_output_O2A_newSVD")
v = L2_VERSION.replace(".", "_")
with nc.Dataset(L2_ROOT / f"L2_BGC_{L2_VERSION}" / DATE / f"PACE_OCI.{STAMP}.L2.OC_BGC.{v}.nc") as ds:
    g = ds["geophysical_data"]
    g.set_auto_mask(True)
    chl = np.ma.filled(g["chlor_a"][:].astype(float), np.nan)
    fl = g["l2_flags"]
    flags = np.ma.filled(fl[:], 0).astype(np.int64)
    names = fl.getncattr("flag_meanings").split()
    masks = np.atleast_1d(fl.getncattr("flag_masks"))
with nc.Dataset(L2_ROOT / f"L2_AOP_{L2_VERSION}" / DATE / f"PACE_OCI.{STAMP}.L2.OC_AOP.{v}.nc") as ds:
    g = ds["geophysical_data"]
    g.set_auto_mask(True)
    nflh = np.ma.filled(g["nflh"][:].astype(float), np.nan)
bad = sum(int(masks[names.index(n)]) for n in BAD_FLAGS if n in names)
chl[chl <= 0] = np.nan

lon, lat = O["longitude"], O["latitude"]
all_ocean = O["is_ocean"] == 1
clear = all_ocean & ((flags & bad) == 0) & np.isfinite(chl)   # dark (open water) and bright
co, cn = O["converged"] == 1, N["converged"] == 1
sif_o, sif_n = O["sif_radiance_678nm"].astype(float), N["sif_radiance_678nm"].astype(float)


def masked(a, m):
    return np.where(m, a, np.nan)


def draw(ax, field, title, cmap, norm=None, vmin=None, vmax=None, label=None, cbar=True):
    kw = dict(norm=norm) if norm is not None else dict(vmin=vmin, vmax=vmax)
    pm = ax.pcolormesh(lon, lat, field, cmap=cmap, shading="nearest", rasterized=True,
                       transform=ccrs.PlateCarree(), **kw)
    ax.add_feature(cfeature.LAND, facecolor="0.85", zorder=2)
    ax.coastlines(linewidth=0.5, zorder=3)
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)
    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = gl.ylabel_style = {"size": 7}
    ax.set_extent([np.nanmin(lon) - 1, np.nanmax(lon) + 1, np.nanmin(lat) - 1, np.nanmax(lat) + 1], ccrs.PlateCarree())
    ax.set_title(title, fontsize=10)
    if cbar:
        cb = plt.colorbar(pm, ax=ax, shrink=0.8, pad=0.03)
        cb.set_label(label or "", fontsize=8)
        cb.ax.tick_params(labelsize=7)
    return pm


def conv_class(m):
    """0 neither, 1 only old, 2 only newSVD, 3 both; NaN outside mask."""
    c = co.astype(int) + 2 * cn.astype(int)
    return np.where(m, c, np.nan)


conv_cmap = ListedColormap(["#bbbbbb", "#1f77b4", "#d62728", "#2ca02c"])
conv_norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], 4)


def conv_colorbar(fig, pm, ax):
    cb = fig.colorbar(pm, ax=ax, shrink=0.8, pad=0.03, ticks=[0, 1, 2, 3])
    cb.ax.set_yticklabels(["neither", "only old", "only newSVD", "both"], fontsize=7)


proj = ccrs.PlateCarree()
fig, axs = plt.subplots(3, 2, figsize=(12, 13), subplot_kw=dict(projection=proj))
sif_kw = dict(cmap="viridis", vmin=-0.05, vmax=0.3, label="SIF@678 [W m$^{-2}$ µm$^{-1}$ sr$^{-1}$]")
draw(axs[0, 0], masked(sif_o, clear), "SIF@678 — old", **sif_kw)
draw(axs[0, 1], masked(sif_n, clear), "SIF@678 — newSVD", **sif_kw)
draw(axs[1, 0], masked(sif_n - sif_o, clear & co & cn), "newSVD − old (converged in both)",
     "RdBu_r", vmin=-0.05, vmax=0.05, label="ΔSIF@678")
draw(axs[1, 1], masked(nflh, clear), f"nFLH ({L2_VERSION})", "viridis", vmin=-0.02, vmax=0.3,
     label="nFLH [W m$^{-2}$ µm$^{-1}$ sr$^{-1}$]")
draw(axs[2, 0], masked(chl, clear), f"chlor_a ({L2_VERSION})", "viridis", norm=LogNorm(0.02, 5), label="chl [mg m$^{-3}$]")
pm = draw(axs[2, 1], conv_class(clear), "Convergence (clear ocean)", conv_cmap, norm=conv_norm, cbar=False)
conv_colorbar(fig, pm, axs[2, 1])
fig.suptitle(f"PACE {STAMP}  O2A SVD retrievals, clear ocean (L2 flags ok, dark + bright)", fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "5_maps_clear_sky.png", dpi=130)

fig, axs = plt.subplots(1, 2, figsize=(12, 5.5), subplot_kw=dict(projection=proj))
draw(axs[0], conv_class(all_ocean), "Convergence — all ocean pixels (incl. cloud)", conv_cmap, norm=conv_norm, cbar=False)
conv_colorbar(fig, axs[0].collections[0], axs[0])
draw(axs[1], masked(np.ones_like(chl), all_ocean & ((flags & bad) != 0)), "Pixels removed by L2 flags (cloud etc.)",
     ListedColormap(["#444444"]), vmin=0, vmax=1, cbar=False)
fig.tight_layout()
fig.savefig(OUT / "6_maps_convergence_all.png", dpi=130)

fig, axs = plt.subplots(3, 2, figsize=(12, 13), subplot_kw=dict(projection=proj))
draw(axs[0, 0], masked(sif_o, all_ocean), "SIF@678 — old (all ocean pixels)", **sif_kw)
draw(axs[0, 1], masked(sif_n, all_ocean), "SIF@678 — newSVD (all ocean pixels)", **sif_kw)
draw(axs[1, 0], masked(sif_n - sif_o, all_ocean), "newSVD − old (all ocean pixels)",
     "RdBu_r", vmin=-0.1, vmax=0.1, label="ΔSIF@678")
pm = draw(axs[1, 1], conv_class(all_ocean), "Convergence (all ocean pixels)", conv_cmap, norm=conv_norm, cbar=False)
conv_colorbar(fig, pm, axs[1, 1])
draw(axs[2, 0], masked(nflh, all_ocean), f"nFLH ({L2_VERSION}, where valid)", "viridis", vmin=-0.02, vmax=0.3,
     label="nFLH [W m$^{-2}$ µm$^{-1}$ sr$^{-1}$]")
draw(axs[2, 1], masked(O["is_dark"].astype(float), all_ocean), "is_dark: black = dark open water, yellow = bright (cloud/edge)",
     ListedColormap(["#f0e442", "#222222"]), vmin=0, vmax=1, cbar=False)
fig.suptitle(f"PACE {STAMP}  O2A SVD retrievals, ALL ocean pixels (no L2-flag / dark / convergence filtering)", fontsize=11)
fig.tight_layout()
fig.savefig(OUT / "7_maps_all_pixels.png", dpi=130)
print("wrote", OUT / "5_maps_clear_sky.png", OUT / "6_maps_convergence_all.png", OUT / "7_maps_all_pixels.png")
