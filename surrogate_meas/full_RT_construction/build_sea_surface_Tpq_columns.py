#!/usr/bin/env python3
"""Sample MERRA-2 T/p/q columns over open ocean for vSmartMOM RT.

Reads the same MERRA-2 reanalysis file used to build the summer transmittance
library, keeps only sea-surface columns with |lat| ≤ 60°, randomly draws
N_SAMPLES profiles, and reduces each from 72 to ~24 layers (every 3rd
half-level, matching `PROFILE_STRIDE=3` in `generate_test_spec.jl`).

Output (under this directory by default):
  merra2_sea_Tpq_columns_n1500.nc
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import interp1d

from cartopy.io import shapereader
from shapely.geometry import Point
from shapely.ops import unary_union
from shapely.prepared import prep

MERRA_PATH = os.environ.get(
    "MERRA_PATH",
    "/home/zhe2/data/MERRA2_reanalysis/MERRA2_400.inst6_3d_ana_Nv.20240705.nc4",
)
OUT_DIR = Path(__file__).resolve().parent
N_SAMPLES = int(os.environ.get("N_SAMPLES", "1500"))
SEED = int(os.environ.get("SEED", "20260919"))
LAT_MIN, LAT_MAX = -60.0, 60.0
PROFILE_STRIDE = int(os.environ.get("PROFILE_STRIDE", "3"))
OUT_NC = OUT_DIR / os.environ.get(
    "OUT_NC", f"merra2_sea_Tpq_columns_n{N_SAMPLES}.nc"
)


def downsample_half_levels(p_half: np.ndarray, stride: int = PROFILE_STRIDE) -> np.ndarray:
    """Keep every `stride`-th half-level, always including the surface."""
    stride = max(int(stride), 1)
    n = len(p_half)
    idx = list(range(0, n, stride))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    return np.asarray(p_half, dtype=np.float64)[idx]


def reduce_profile(
    T72: np.ndarray,
    q72: np.ndarray,
    p_half72: np.ndarray,
    stride: int = PROFILE_STRIDE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Map 72-layer MERRA column onto a ~24-layer half-level grid.

    Half-levels are subsampled (stride), then T/q are log-p interpolated onto
    the new layer mid-pressures (same idea as generate_test_spec.jl).
    """
    p_half72 = np.asarray(p_half72, dtype=np.float64)
    if not (p_half72[0] < p_half72[-1]):
        raise ValueError("Expected TOA→BOA half-levels (increasing p)")

    p_half = downsample_half_levels(p_half72, stride=stride)
    p_full_src = 0.5 * (p_half72[:-1] + p_half72[1:])
    p_full = 0.5 * (p_half[:-1] + p_half[1:])

    log_src = np.log(p_full_src)
    log_dst = np.log(p_full)
    T = interp1d(log_src, np.asarray(T72, dtype=np.float64), kind="linear",
                 bounds_error=False, fill_value=(float(T72[0]), float(T72[-1])))(log_dst)
    q = interp1d(log_src, np.asarray(q72, dtype=np.float64), kind="linear",
                 bounds_error=False, fill_value=(float(q72[0]), float(q72[-1])))(log_dst)
    q = np.maximum(q, 0.0)
    return p_half, p_full, T.astype(np.float64), q.astype(np.float64)


def build_ocean_mask(lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
    """Boolean (nlat, nlon) True over open ocean within LAT_MIN..LAT_MAX."""
    print("Building Natural Earth ocean mask…")
    land_shp = shapereader.natural_earth(
        resolution="110m", category="physical", name="land"
    )
    land = prep(unary_union(list(shapereader.Reader(land_shp).geometries())))

    lon2d, lat2d = np.meshgrid(lons, lats)
    is_ocean = np.zeros(lon2d.shape, dtype=bool)
    for j in range(lats.size):
        la = float(lats[j])
        if not (LAT_MIN <= la <= LAT_MAX):
            continue
        for i in range(lons.size):
            lo = float(lons[i])
            if not land.contains(Point(lo, la)):
                is_ocean[j, i] = True
    print(f"Ocean cells in [{LAT_MIN},{LAT_MAX}]: {is_ocean.sum()} / {is_ocean.size}")
    return is_ocean


def main() -> None:
    print(f"Reading {MERRA_PATH}")
    ds = Dataset(MERRA_PATH)
    lons = np.asarray(ds["lon"][:], dtype=np.float64)
    lats = np.asarray(ds["lat"][:], dtype=np.float64)
    ak = np.asarray(ds.getncattr("ak"), dtype=np.float64)  # Pa
    bk = np.asarray(ds.getncattr("bk"), dtype=np.float64)
    n_time = ds.dimensions["time"].size
    n_lev = ds.dimensions["lev"].size
    assert n_lev == 72, f"expected 72 levels, got {n_lev}"
    assert len(ak) == n_lev + 1

    is_ocean = build_ocean_mask(lons, lats)
    ocean_ji = np.argwhere(is_ocean)  # (j_lat, i_lon)
    if ocean_ji.shape[0] == 0:
        raise RuntimeError("No ocean cells found in latitude band")

    # Candidate pool: all times × ocean lat/lon
    n_ocean = ocean_ji.shape[0]
    n_cand = n_time * n_ocean
    print(f"Candidate columns: {n_cand} (= {n_time} times × {n_ocean} ocean cells)")
    if n_cand < N_SAMPLES:
        raise RuntimeError(f"Need {N_SAMPLES} samples but only {n_cand} ocean candidates")

    rng = np.random.default_rng(SEED)
    picks = rng.choice(n_cand, size=N_SAMPLES, replace=False)

    # Probe one column for reduced layer count
    j0, i0 = map(int, ocean_ji[0])
    ps0 = float(ds["PS"][0, j0, i0])  # Pa
    p_half0 = (ak + bk * ps0) / 100.0
    T0 = np.asarray(ds["T"][0, :, j0, i0], dtype=np.float64)
    q0 = np.asarray(ds["QV"][0, :, j0, i0], dtype=np.float64)
    ph_r, pf_r, _, _ = reduce_profile(T0, q0, p_half0)
    n_half = len(ph_r)
    n_layer = len(pf_r)
    print(f"Layer reduction: {n_lev} → {n_layer} (half-levels {n_lev + 1} → {n_half}, "
          f"stride={PROFILE_STRIDE})")

    T_out = np.zeros((n_layer, N_SAMPLES), dtype=np.float32)
    q_out = np.zeros((n_layer, N_SAMPLES), dtype=np.float32)
    p_mid_out = np.zeros((n_layer, N_SAMPLES), dtype=np.float32)
    p_half_out = np.zeros((n_half, N_SAMPLES), dtype=np.float32)
    ps_out = np.zeros(N_SAMPLES, dtype=np.float32)
    lat_out = np.zeros(N_SAMPLES, dtype=np.float32)
    lon_out = np.zeros(N_SAMPLES, dtype=np.float32)
    time_out = np.zeros(N_SAMPLES, dtype=np.int32)
    j_out = np.zeros(N_SAMPLES, dtype=np.int32)
    i_out = np.zeros(N_SAMPLES, dtype=np.int32)

    print(f"Sampling {N_SAMPLES} columns…")
    for s, flat in enumerate(picks):
        it = int(flat // n_ocean)
        k = int(flat % n_ocean)
        j, i = map(int, ocean_ji[k])
        ps_pa = float(ds["PS"][it, j, i])
        T72 = np.asarray(ds["T"][it, :, j, i], dtype=np.float64)
        q72 = np.asarray(ds["QV"][it, :, j, i], dtype=np.float64)
        p_half72 = (ak + bk * ps_pa) / 100.0  # hPa, TOA→BOA
        p_half, p_full, T, q = reduce_profile(T72, q72, p_half72)

        T_out[:, s] = T
        q_out[:, s] = q
        p_mid_out[:, s] = p_full
        p_half_out[:, s] = p_half
        ps_out[s] = ps_pa / 100.0
        lat_out[s] = lats[j]
        lon_out[s] = lons[i]
        time_out[s] = it
        j_out[s] = j
        i_out[s] = i
        if (s + 1) % 250 == 0 or s == 0:
            print(f"  {s + 1}/{N_SAMPLES}  lat={lat_out[s]:.1f} lon={lon_out[s]:.1f} "
                  f"ps={ps_out[s]:.1f} hPa")

    ds.close()

    OUT_NC.parent.mkdir(parents=True, exist_ok=True)
    if OUT_NC.exists():
        OUT_NC.unlink()
    print(f"Writing {OUT_NC}")
    with Dataset(OUT_NC, "w") as out:
        out.createDimension("sample", N_SAMPLES)
        out.createDimension("layer", n_layer)
        out.createDimension("half", n_half)
        out.createDimension("half_full", len(ak))  # native MERRA-2 (73)

        def put(name, data, dims, **attrs):
            var = out.createVariable(name, data.dtype, dims, zlib=True, complevel=4)
            var[:] = data
            for k, v in attrs.items():
                setattr(var, k, v)
            return var

        put("lat", lat_out, ("sample",), units="degrees_north",
            long_name="latitude")
        put("lon", lon_out, ("sample",), units="degrees_east",
            long_name="longitude")
        put("ps", ps_out, ("sample",), units="hPa",
            long_name="surface pressure")
        put("time_index", time_out, ("sample",),
            long_name="0-based MERRA-2 time index in source file")
        put("lat_index", j_out, ("sample",),
            long_name="0-based MERRA-2 lat index")
        put("lon_index", i_out, ("sample",),
            long_name="0-based MERRA-2 lon index")
        put("T", T_out, ("layer", "sample"), units="K",
            long_name="temperature on reduced mid-levels (TOA→BOA)")
        put("q", q_out, ("layer", "sample"), units="kg/kg",
            long_name="specific humidity on reduced mid-levels (TOA→BOA)")
        put("p_mid", p_mid_out, ("layer", "sample"), units="hPa",
            long_name="layer mid pressure (TOA→BOA)")
        put("p_half", p_half_out, ("half", "sample"), units="hPa",
            long_name="half-level pressure (TOA→BOA, increasing p)")

        put("ak", ak.astype(np.float64), ("half_full",), units="Pa",
            long_name="MERRA-2 hybrid A coefficient (native 73 half-levels)")
        put("bk", bk.astype(np.float64), ("half_full",), units="1",
            long_name="MERRA-2 hybrid B coefficient (native 73 half-levels)")

        out.title = "MERRA-2 open-ocean T/p/q columns (24-layer reduced)"
        out.source_file = MERRA_PATH
        out.seed = np.int32(SEED)
        out.profile_stride = np.int32(PROFILE_STRIDE)
        out.lat_min = np.float32(LAT_MIN)
        out.lat_max = np.float32(LAT_MAX)
        out.n_samples = np.int32(N_SAMPLES)
        out.vertical_order = "TOA→BOA (p_half increasing)"
        out.ocean_mask = "Natural Earth 110m land; cells not containing land"
        out.history = (
            f"Built by build_sea_surface_Tpq_columns.py from {Path(MERRA_PATH).name}; "
            f"sea-only, lat∈[{LAT_MIN},{LAT_MAX}], N={N_SAMPLES}, "
            f"stride={PROFILE_STRIDE} → {n_layer} layers"
        )

    print(f"Done: {OUT_NC}  ({N_SAMPLES} samples × {n_layer} layers)")


if __name__ == "__main__":
    main()
