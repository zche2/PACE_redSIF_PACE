#!/usr/bin/env python3
"""GCHP aerosol diagnostics + 500 open-ocean columns for vSmartMOM RT.

Converted from GCHP_aerosol.ipynb.

1. Global sea-salt column maps (fine / mid / coarse / total) and a few
   open-ocean vertical profiles.
2. Rebuild `output_aerosol_profiles/gchp_ocean_columns_n500.nc` with the
   GETDP / Mie physics in `build_gchp_ocean_columns.py`.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_gchp_ocean_columns import (
    M_AIR,
    NC_PATH,
    OUT_DIR,
    OUT_NC,
    OUT_MAP,
    SPC_MW,
    main as build_ocean_columns,
)

# TOMAS SS molar mass (species_database_tomas.yml). The notebook used 31.4 g/mol,
# which is bulk GEOS-Chem SALA, not the sectional SS tracer.
MW_SS = SPC_MW["SS"]

BIN_GROUPS = {
    "fine_SS01-05": list(range(1, 6)),
    "mid_SS06-10": list(range(6, 11)),
    "coarse_SS11-15": list(range(11, 16)),
    "total_SS01-15": list(range(1, 16)),
}


def plot_ss_maps(lons, lats, col):
    def flat_lonlat(field):
        return lons.ravel(), lats.ravel(), np.asarray(field, dtype=float).ravel()

    def plot_ss_map(field, title, fname, vmax=None):
        lon, lat, val = flat_lonlat(field)
        fig = plt.figure(figsize=(11, 5))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        if vmax is None:
            vmax = np.nanpercentile(val, 99)
        sc = ax.scatter(
            lon, lat, c=val, s=12, cmap="YlGnBu",
            vmin=0.0, vmax=vmax, transform=ccrs.PlateCarree(),
            linewidths=0, rasterized=True,
        )
        ax.coastlines(linewidth=0.6)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.4)
        ax.set_global()
        ax.set_title(title)
        cb = fig.colorbar(sc, ax=ax, shrink=0.75, pad=0.02)
        cb.set_label("sea-salt column mass (kg / grid box)")
        fig.tight_layout()
        out = os.path.join(OUT_DIR, fname)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print("Wrote", out)

    plot_ss_map(col["total_SS01-15"], "GCHP sea-salt column (SS01–SS15)", "gchp_seasalt_total_map.png")
    plot_ss_map(col["fine_SS01-05"], "GCHP sea-salt fine bins (SS01–SS05)", "gchp_seasalt_fine_map.png")
    plot_ss_map(col["mid_SS06-10"], "GCHP sea-salt mid bins (SS06–SS10)", "gchp_seasalt_mid_map.png")
    plot_ss_map(col["coarse_SS11-15"], "GCHP sea-salt coarse bins (SS11–SS15)", "gchp_seasalt_coarse_map.png")
    return flat_lonlat


def select_profile_sites(lons, lats, col_total):
    lon_f, lat_f, tot_f = lons.ravel(), lats.ravel(), np.asarray(col_total, dtype=float).ravel()
    idx_flat = np.arange(tot_f.size)
    ocean = (np.abs(lat_f) < 60.0) & (tot_f > 0)
    cand = idx_flat[ocean]
    order = cand[np.argsort(tot_f[cand])[::-1]]
    picked = []
    for min_sep in (35, 25, 18, 12):
        picked = []
        for i in order:
            if len(picked) >= 6:
                break
            ok = True
            for j in picked:
                dlon = abs((lon_f[i] - lon_f[j] + 180) % 360 - 180)
                dlat = abs(lat_f[i] - lat_f[j])
                if dlon < min_sep and dlat < min_sep:
                    ok = False
                    break
            if ok:
                picked.append(int(i))
        if len(picked) >= 6:
            break

    nf, ny, nx = lons.shape
    coords = [np.unravel_index(i, (nf, ny, nx)) for i in picked]
    print("Selected ocean columns:")
    for i, (f, y, x) in enumerate(coords):
        print(
            f"  {i+1}: face={f} Y={y} X={x}  lon={lons[f,y,x]:.1f} lat={lats[f,y,x]:.1f}  "
            f"SS_total={col_total[f,y,x]:.3e} kg/box"
        )

    fig = plt.figure(figsize=(11, 5))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.scatter(
        lon_f, lat_f, c=tot_f, s=8, cmap="YlGnBu", alpha=0.35,
        transform=ccrs.PlateCarree(), linewidths=0, rasterized=True,
    )
    for k, (f, y, x) in enumerate(coords):
        ax.plot(
            lons[f, y, x], lats[f, y, x], "o", ms=8, transform=ccrs.PlateCarree(),
            label=f"{k+1}: ({lats[f,y,x]:.0f}°, {lons[f,y,x]:.0f}°)",
        )
    ax.coastlines(linewidth=0.6)
    ax.set_global()
    ax.set_title("GCHP ocean columns for SS vertical profiles")
    ax.legend(loc="lower left", fontsize=8, frameon=True)
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "gchp_seasalt_sites.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Wrote", out)
    return coords


def plot_ss_profiles(lons, lats, p_mid, ss_mass, coords):
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharey=True)
    axes = axes.ravel()
    bin_sets = [
        ("fine", list(range(1, 6))),
        ("mid", list(range(6, 11))),
        ("coarse", list(range(11, 16))),
        ("total", list(range(1, 16))),
    ]
    colors = {"fine": "C0", "mid": "C2", "coarse": "C3", "total": "k"}

    p_all = []
    for ax, (f, y, x) in zip(axes, coords):
        p = p_mid[:, f, y, x]
        p_all.append(p)
        for label, bins in bin_sets:
            prof = np.sum([ss_mass[i][:, f, y, x] for i in bins], axis=0)
            ax.plot(
                prof, p, color=colors[label],
                lw=1.6 if label != "total" else 2.0,
                ls="-" if label != "total" else "--",
                label=label,
            )
        ax.set_xlabel("SS mass (kg / layer)")
        ax.set_title(f"lon={lons[f,y,x]:.0f}° lat={lats[f,y,x]:.0f}°")
        ax.grid(True, alpha=0.3)

    pmax = max(float(np.nanmax(p)) for p in p_all)
    axes[0].set_ylim(pmax, 0.0)
    axes[0].set_ylabel("pressure (hPa)")
    axes[3].set_ylabel("pressure (hPa)")
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle("GCHP sea-salt vertical profiles (selected ocean columns)", y=1.01)
    fig.tight_layout()
    out = os.path.join(OUT_DIR, "gchp_seasalt_profiles.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Wrote", out)


def diagnose_n500():
    cols = Dataset(OUT_NC)
    print("samples", cols.dimensions["sample"].size, "lev", cols.dimensions["lev"].size)
    print("method", getattr(cols, "method", ""))
    print("tau_ref", float(cols["tau_ref"][:].min()), "..", float(cols["tau_ref"][:].max()))
    print("ss_col", float(cols["ss_col"][:].min()), "..", float(cols["ss_col"][:].max()))
    print("vars", [v for v in cols.variables if v.startswith("AOD_GCHP") or v.startswith("RadiWet")])
    cols.close()
    print("sites map:", OUT_MAP)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Reading {NC_PATH}")
    ds = Dataset(NC_PATH)
    lons = np.asarray(ds["lons"][:], dtype=float)
    lats = np.asarray(ds["lats"][:], dtype=float)
    AD = np.asarray(ds["Met_AD"][0], dtype=float)
    DELP = np.asarray(ds["Met_DELP"][0], dtype=float)
    PS = np.asarray(ds["Met_PS1WET"][0], dtype=float)

    ss_mass = {}
    for i in range(1, 16):
        vmr = np.asarray(ds[f"SpeciesConcVV_SS{i:02d}"][0], dtype=float)
        ss_mass[i] = vmr * AD * (MW_SS / M_AIR)

    def column_mass(bins):
        return np.sum([ss_mass[i].sum(axis=0) for i in bins], axis=0)

    col = {name: column_mass(bins) for name, bins in BIN_GROUPS.items()}
    for name, field in col.items():
        print(f"{name}: {np.nanmin(field):.3e} .. {np.nanmax(field):.3e} kg/box")

    nlev = DELP.shape[0]
    p_edge = np.empty((nlev + 1,) + DELP.shape[1:], dtype=float)
    p_edge[0] = PS
    np.cumsum(DELP, axis=0, out=p_edge[1:])
    p_edge[1:] = PS - p_edge[1:]
    p_edge = np.maximum(p_edge, 1e-4)
    p_edge = np.minimum.accumulate(p_edge, axis=0)
    p_mid = 0.5 * (p_edge[:-1] + p_edge[1:])

    plot_ss_maps(lons, lats, col)
    coords = select_profile_sites(lons, lats, col["total_SS01-15"])
    plot_ss_profiles(lons, lats, p_mid, ss_mass, coords)
    ds.close()

    print("\n=== rebuilding", OUT_NC, "===")
    build_ocean_columns()
    print("\n=== n500 column file ===")
    diagnose_n500()


if __name__ == "__main__":
    main()
