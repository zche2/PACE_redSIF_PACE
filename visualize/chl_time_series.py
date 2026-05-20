import earthaccess
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from datetime import datetime
import pandas as pd

earthaccess.login()

# ── Time range ───────────────────────────────────────────────────────────────
start_date = datetime(2025, 4, 1)
end_date   = datetime(2025, 7, 1)

# ── Search for 4 km daily CHL granules ──────────────────────────────────────
results = earthaccess.search_data(
    short_name="PACE_OCI_L3M_CHL",
    temporal=(start_date.strftime("%Y-%m-%d"), end_date.strftime("%Y-%m-%d")),
)
results = [r for r in results if "4km" in r.data_links()[0]]
print(f"Found {len(results)} × 4 km CHL granules")

date_list = pd.to_datetime([
    r["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"][:10]
    for r in results
])


def load_chl_data(granule):
    """Load chlor_a from a single DataGranule; return (lat, lon, chl) or Nones."""
    try:
        fo = earthaccess.open(granules=[granule])
        with xr.open_dataset(fo[0]) as ds:
            if "chlor_a" not in ds:
                print(f"  chlor_a not found in {granule}")
                return None, None, None
            chl = ds["chlor_a"].load().values
            lat = ds["lat"].values
            lon = ds["lon"].values
        return lat, lon, chl
    except Exception as e:
        print(f"  Error: {e}")
        return None, None, None


# ── Determine global colour scale from a sample of dates ────────────────────
print("Sampling colour scale …")
sample = results[::max(1, len(results) // 10)][:10]
all_vals = []
for g in sample:
    _, _, chl = load_chl_data(g)
    if chl is not None:
        vals = chl[~np.isnan(chl)]
        if len(vals):
            all_vals.extend(np.log10(vals[vals > 0]).tolist())

if all_vals:
    vmin = np.percentile(all_vals, 2)
    vmax = np.percentile(all_vals, 98)
else:
    vmin, vmax = -2, 1          # log10 fallback: 0.01 – 10 mg m⁻³

levels = np.linspace(vmin, vmax, 20)
print(f"log₁₀(Chl-a) colour range: {vmin:.2f} → {vmax:.2f}")

# ── Set up figure ────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 8))
ax  = plt.axes(projection=ccrs.PlateCarree())
ax.add_feature(cfeature.LAND,      facecolor="lightgray", edgecolor="black", linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.BORDERS,   linewidth=0.3, alpha=0.5)
ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5, linestyle="--")
ax.set_global()

contour    = None
title_text = ax.set_title("", fontsize=14, fontweight="bold")

sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation="horizontal", pad=0.05, shrink=0.8, extend="both")
cbar.set_label("log₁₀(Chlorophyll-a)  [mg m⁻³]", fontsize=11)
cbar.set_ticks(np.linspace(vmin, vmax, 5))
cbar.set_ticklabels([f"$10^{{{v:.1f}}}$" for v in np.linspace(vmin, vmax, 5)])


# ── Animation update ─────────────────────────────────────────────────────────
def update(frame):
    global contour

    date    = date_list[frame]
    granule = results[frame]
    lat, lon, chl = load_chl_data(granule)

    if chl is None:
        return

    if contour is not None:
        try:
            contour.remove()
        except (AttributeError, ValueError):
            if hasattr(contour, "collections"):
                for c in contour.collections:
                    try:
                        c.remove()
                    except ValueError:
                        pass

    with np.errstate(divide="ignore", invalid="ignore"):
        chl_log = np.where(chl > 0, np.log10(chl), np.nan)

    chl_masked = np.ma.masked_invalid(chl_log)

    contour = ax.contourf(
        lon, lat, chl_masked,
        levels=levels,
        transform=ccrs.PlateCarree(),
        cmap="viridis",
        extend="both",
        vmin=vmin, vmax=vmax,
    )

    title_text.set_text(
        f"PACE OCI 4 km: Global Chlorophyll-a\n{date.strftime('%Y-%m-%d')}"
    )
    print(f"  Frame {frame + 1}/{len(date_list)}: {date.strftime('%Y-%m-%d')}")


# ── Build and save animation ─────────────────────────────────────────────────
print("Creating animation …")
anim = FuncAnimation(fig, update, frames=len(date_list), interval=100, repeat=True)

out = (f"pace_chl_animation"
       f"_{start_date.strftime('%Y%m%d')}"
       f"_{end_date.strftime('%Y%m%d')}.gif")

print(f"Saving → {out}")
anim.save(out, writer=PillowWriter(fps=10), dpi=150)
print("Done.")
plt.close()
