import earthaccess
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from datetime import datetime, timedelta
import pandas as pd

# Authenticate with NASA Earthdata
earthaccess.login()

# Generate date range
start_date = datetime(2025, 7, 1)
end_date = datetime(2025, 9, 1)
date_list = pd.date_range(start_date, end_date, freq='D')

# Alternatively, construct URLs manually for specific dates
def get_pace_url(date):
    """Construct PACE OCI L3b URL for a given date"""
    date_str = date.strftime('%Y%m%d')
    url = f"https://obdaac-tea.earthdatacloud.nasa.gov/ob-cumulus-prod-public/PACE_OCI.{date_str}.L3m.DAY.FLH.V3_1.nflh.0p1deg.nc"
    return url

# Function to load and process data for a single date
def load_nflh_data(date):
    """Load nFLH data for a specific date"""
    try:
        url = get_pace_url(date)
        granule = earthaccess.open(granules=[url])
        ds = xr.open_dataset(granule[0])
        
        # Extract nFLH data (adjust variable name as needed)
        if 'nflh' in ds:
            nflh = ds['nflh'].values
        elif 'chlor_a' in ds:  # Fallback to chlor_a if nflh not available
            nflh = ds['chlor_a'].values
        else:
            print(f"Warning: nflh variable not found for {date}")
            return None, None, None
        
        lat = ds['lat'].values
        lon = ds['lon'].values
        
        ds.close()
        return lat, lon, nflh
    except Exception as e:
        print(f"Error loading data for {date}: {e}")
        return None, None, None

# Create figure
fig = plt.figure(figsize=(14, 8))
ax = plt.axes(projection=ccrs.PlateCarree())

# Add map features
ax.add_feature(cfeature.LAND, facecolor='lightgray', edgecolor='black', linewidth=0.5)
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.5)
ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5, linestyle='--')
ax.set_global()

# Initialize plot elements
contour = None
title_text = ax.set_title('', fontsize=14, fontweight='bold')

# Determine global color scale by sampling a few dates
print("Determining global color scale...")
sample_dates = date_list[::30]  # Sample every 30 days
all_values = []
for date in sample_dates[:10]:  # Sample first 10
    _, _, nflh = load_nflh_data(date)
    if nflh is not None:
        all_values.extend(nflh[~np.isnan(nflh)])

vmin = np.nanpercentile(all_values, 2)
vmax = np.nanpercentile(all_values, 98)
levels = np.linspace(vmin, vmax, 20)

print(f"Color scale: {vmin:.4f} to {vmax:.4f}")

# Animation update function
def update(frame):
    global contour
    
    date = date_list[frame]
    lat, lon, nflh = load_nflh_data(date)
    
    if lat is None:
        return
    
    # Remove previous contour
    if contour is not None:
        try:
            contour.remove()
        except (AttributeError, ValueError):
            # If remove() doesn't work, try removing collections
            if hasattr(contour, 'collections'):
                for coll in contour.collections:
                    try:
                        coll.remove()
                    except ValueError:
                        pass
    
    # Mask invalid values
    nflh_masked = np.ma.masked_invalid(nflh)
    
    # Plot new data
    contour = ax.contourf(lon, lat, nflh_masked,
                          levels=levels,
                          transform=ccrs.PlateCarree(),
                          cmap='RdYlGn',
                          extend='both',
                          vmin=vmin,
                          vmax=vmax)
    
    # Update title
    title_text.set_text(f'PACE OCI: Global nFLH Distribution\n{date.strftime("%Y-%m-%d")}')
    
    print(f"Processed frame {frame+1}/{len(date_list)}: {date.strftime('%Y-%m-%d')}")

# Add colorbar (only once)
sm = plt.cm.ScalarMappable(cmap='RdYlGn', norm=plt.Normalize(vmin=vmin, vmax=vmax))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', pad=0.05, shrink=0.8, extend='both')
cbar.set_label('Normalized Fluorescence Line Height (nFLH) [W m⁻² μm⁻¹ sr⁻¹]', fontsize=11)

# Create animation
print("Creating animation...")
anim = FuncAnimation(fig, update, frames=len(date_list), interval=100, repeat=True)

# Save animation
print("Saving animation...")
writer = PillowWriter(fps=10)
anim.save(f'pace_nflh_animation_{start_date.strftime("%Y%m%d")}_{end_date.strftime("%Y%m%d")}.gif', writer=writer, dpi=150)

print(f"Animation saved as 'pace_nflh_animation_{start_date.strftime('%Y%m%d')}_{end_date.strftime('%Y%m%d')}.gif'")
plt.close()