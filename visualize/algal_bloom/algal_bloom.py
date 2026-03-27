import earthaccess
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from datetime import datetime
import pandas as pd

# ── User-configurable settings ────────────────────────────────────────────────

# Product to visualise: "chlor_a" or "nflh"
PRODUCT = "chlor_a"

REGION_OF_INTEREST = "South Georgia"


# ── region catalogue ─────────────────────────────────────────────────────────

REGION_CONFIG = {
    "Southern Australia": {
        "bbox": (135, -36, 140, -32),
        "start_date": datetime(2025, 2, 1),
        "end_date": datetime(2025, 9, 30),
    },
    "South Georgia": {
        "bbox": (-42, -55, -35, -46),
        "start_date": datetime(2026, 1, 1),
        "end_date": datetime(2026, 1, 30),
    },
}

if REGION_OF_INTEREST not in REGION_CONFIG:
    raise ValueError(f"Unknown region '{REGION_OF_INTEREST}'. Choose from: {list(REGION_CONFIG)}")

# extract bbox, start_date, end_date
BBOX = REGION_CONFIG[REGION_OF_INTEREST]["bbox"]
start_date = REGION_CONFIG[REGION_OF_INTEREST]["start_date"]
end_date = REGION_CONFIG[REGION_OF_INTEREST]["end_date"]

# ── Product catalogue ─────────────────────────────────────────────────────────
# Each entry defines how to build the URL, which variable to read,
# and how to render the data (colormap, normalization, axis label).

PRODUCT_CONFIG = {
    "chlor_a": {
        "suite":    "CHL",
        "variable": "chlor_a",
        "cmap":     "YlGn",
        "norm":     mcolors.LogNorm(vmin=0.01, vmax=20.0),
        "cbar_label": "Chlorophyll-a [mg m⁻³]",
        "title_var":  "Chlorophyll-a",
    },
    "nflh": {
        "suite":    "FLH",
        "variable": "nflh",
        "cmap":     "RdYlGn",
        "norm":     mcolors.Normalize(vmin=-0.1, vmax=0.8),
        "cbar_label": "nFLH [W m⁻² μm⁻¹ sr⁻¹]",
        "title_var":  "nFLH",
    },
}

if PRODUCT not in PRODUCT_CONFIG:
    raise ValueError(f"Unknown product '{PRODUCT}'. Choose from: {list(PRODUCT_CONFIG)}")

cfg = PRODUCT_CONFIG[PRODUCT]

# ── Authentication & fsspec session ──────────────────────────────────────────

earthaccess.login()
fs = earthaccess.get_fsspec_https_session()

# ── Date list ─────────────────────────────────────────────────────────────────

date_list = pd.date_range(start_date, end_date, freq='D')

# ── Data access helpers ───────────────────────────────────────────────────────

def get_pace_url(date):
    """Construct PACE OCI L3m DAY URL for the selected product."""
    date_str = date.strftime('%Y%m%d')
    suite    = cfg["suite"]
    variable = cfg["variable"]
    return (
        f"https://obdaac-tea.earthdatacloud.nasa.gov/ob-cumulus-prod-public/"
        f"PACE_OCI.{date_str}.L3m.DAY.{suite}.V3_1.{variable}.0p1deg.nc"
    )


def load_data(date):
    """Load the selected product for one date, subset to BBOX."""
    try:
        url      = get_pace_url(date)
        variable = cfg["variable"]

        ds = xr.open_dataset(fs.open(url))

        if variable not in ds:
            print(f"Warning: '{variable}' not found for {date.date()}")
            ds.close()
            return None, None, None

        lon_min, lat_min, lon_max, lat_max = BBOX
        ds_sub = ds.sel(
            lon=slice(lon_min, lon_max),
            lat=slice(lat_max, lat_min),  # lat is descending in L3m files
        )

        data = ds_sub[variable].values
        lat  = ds_sub["lat"].values
        lon  = ds_sub["lon"].values

        ds.close()
        return lat, lon, data

    except Exception as e:
        print(f"Error loading data for {date.date()}: {e}")
        return None, None, None

# ── Figure setup ──────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(12, 8))
ax  = plt.axes(projection=ccrs.PlateCarree())

ax.add_feature(cfeature.LAND,      facecolor="lightgray", edgecolor="black", linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.BORDERS,   linewidth=0.3, alpha=0.5)
ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5, linestyle="--")

lon_min, lat_min, lon_max, lat_max = BBOX
ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

img        = None
title_text = ax.set_title("", fontsize=14, fontweight="bold")

sm   = plt.cm.ScalarMappable(cmap=cfg["cmap"], norm=cfg["norm"])
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation="horizontal", pad=0.08, shrink=0.8, extend="both")
cbar.set_label(cfg["cbar_label"], fontsize=11)

# ── Animation ─────────────────────────────────────────────────────────────────

def update(frame):
    global img

    date = date_list[frame]
    lat, lon, data = load_data(date)

    if lat is None:
        return

    if img is not None:
        img.remove()

    data_masked = np.ma.masked_invalid(data)

    img = ax.pcolormesh(
        lon, lat, data_masked,
        transform=ccrs.PlateCarree(),
        cmap=cfg["cmap"],
        norm=cfg["norm"],
        shading="auto",
    )

    title_text.set_text(
        f"PACE OCI: {REGION_OF_INTEREST} {cfg['title_var']}\n{date.strftime('%Y-%m-%d')}"
    )
    print(f"Frame {frame + 1}/{len(date_list)}: {date.strftime('%Y-%m-%d')}")


print("Creating animation...")
anim = FuncAnimation(fig, update, frames=len(date_list), interval=100, repeat=True)

print("Saving animation...")
out_file = (
    f"pace_{PRODUCT}_{REGION_OF_INTEREST.replace(' ', '_')}_"
    f"{start_date.strftime('%Y%m%d')}_{end_date.strftime('%Y%m%d')}.gif"
)
anim.save(out_file, writer=PillowWriter(fps=10), dpi=150)

print(f"Animation saved as '{out_file}'")
plt.close()
