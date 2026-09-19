#!/usr/bin/env python3
"""Add Met_T / Met_SPHU / p_half to gchp_ocean_columns_n500.nc.

Writes a new NetCDF (avoids HDF r+ while notebooks hold the file open), then
replaces the original. Levels remain BOA-first (same as p_mid).
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from build_gchp_ocean_columns import NC_PATH, OUT_NC

OUT = Path(OUT_NC)
TMP = OUT.with_suffix(".nc.tmp")


def column_p_half(ps, delp):
    """BOA-first half-levels [hPa] matching build_gchp_ocean_columns."""
    nlev = delp.shape[0]
    p_edge = np.empty(nlev + 1, dtype=np.float64)
    p_edge[0] = ps
    np.cumsum(delp, out=p_edge[1:])
    p_edge[1:] = ps - p_edge[1:]
    p_edge = np.maximum(p_edge, 1e-4)
    return np.minimum.accumulate(p_edge)


def main():
    print(f"Reading columns from {OUT}")
    din = Dataset(OUT, "r")
    nlev = len(din.dimensions["lev"])
    ns = len(din.dimensions["sample"])
    face = np.asarray(din["i_face"][:], dtype=np.int32)
    iy = np.asarray(din["i_y"][:], dtype=np.int32)
    ix = np.asarray(din["i_x"][:], dtype=np.int32)
    p_mid = np.asarray(din["p_mid"][:])

    print(f"Reading met from {NC_PATH}")
    ds = Dataset(NC_PATH)
    T_all = np.asarray(ds["Met_T"][0], dtype=np.float64)
    SPHU = np.asarray(ds["Met_SPHU"][0], dtype=np.float64)
    DELP = np.asarray(ds["Met_DELP"][0], dtype=np.float64)
    PS = np.asarray(ds["Met_PS1WET"][0], dtype=np.float64)
    ds.close()

    T_s = np.zeros((nlev, ns), dtype=np.float32)
    q_s = np.zeros((nlev, ns), dtype=np.float32)
    p_half_s = np.zeros((nlev + 1, ns), dtype=np.float32)
    for s in range(ns):
        f, y, x = int(face[s]), int(iy[s]), int(ix[s])
        T_s[:, s] = T_all[:, f, y, x]
        q_s[:, s] = np.maximum(SPHU[:, f, y, x], 0.0) * 1e-3
        p_edge = column_p_half(PS[f, y, x], DELP[:, f, y, x])
        p_half_s[:, s] = p_edge
        p_mid_1 = 0.5 * (p_edge[:-1] + p_edge[1:])
        if not np.allclose(p_mid[:, s], p_mid_1, rtol=1e-3, atol=0.05):
            print(f"  warn sample {s}: p_mid drift max={np.max(np.abs(p_mid[:, s] - p_mid_1)):.3f} hPa")

    if TMP.exists():
        TMP.unlink()
    print(f"Writing {TMP}")
    shutil.copy2(OUT, TMP)
    dout = Dataset(TMP, "r+")
    if "half" not in dout.dimensions:
        dout.createDimension("half", nlev + 1)

    def upsert(name, data, dims, **attrs):
        if name in dout.variables:
            dout[name][:] = np.asarray(data)
        else:
            v = dout.createVariable(name, np.float32, dims)
            v[:] = np.asarray(data)
            for k, val in attrs.items():
                v.setncattr(k, val)
        for k, val in attrs.items():
            dout[name].setncattr(k, val)

    upsert("T", T_s, ("lev", "sample"), units="K",
           long_name="GCHP Met_T (BOA-first, same as p_mid)")
    upsert("q", q_s, ("lev", "sample"), units="kg kg-1",
           long_name="GCHP specific humidity from Met_SPHU (BOA-first)")
    upsert("p_half", p_half_s, ("half", "sample"), units="hPa",
           long_name="GCHP pressure half-levels (BOA=index0 → TOA)")
    dout.met_source = NC_PATH
    # Copy global attrs already present; close handles
    dout.close()
    din.close()

    bak = OUT.with_suffix(".nc.bak")
    if bak.exists():
        bak.unlink()
    OUT.rename(bak)
    TMP.rename(OUT)
    print(f"Replaced {OUT} (backup {bak})")
    print(f"  T: {T_s.min():.1f} .. {T_s.max():.1f} K")
    print(f"  q: {q_s.min():.2e} .. {q_s.max():.2e} kg/kg")
    print(f"  p_half surface: {p_half_s[0].min():.1f} .. {p_half_s[0].max():.1f} hPa")
    print(f"  p_half top:     {p_half_s[-1].min():.3f} .. {p_half_s[-1].max():.3f} hPa")


if __name__ == "__main__":
    main()
