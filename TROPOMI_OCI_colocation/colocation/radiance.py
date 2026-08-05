"""Units and radiance helpers."""

from __future__ import annotations

import numpy as np

Na = 6.022e23
h = 6.626e-34
c_light = 2.998e8


def tropomi_mol_to_mW(radiance_mol, wavelength_nm):
    """mol s-1 m-2 nm-1 sr-1 → mW m-2 nm-1 sr-1."""
    return radiance_mol * (h * c_light * Na / wavelength_nm) * 1e9 * 1e3


def fill_nan_1d(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float).copy()
    ok = np.isfinite(y)
    if not ok.any():
        return y
    idx = np.where(ok, np.arange(y.size), 0)
    np.maximum.accumulate(idx, out=idx)
    y[~ok] = y[idx[~ok]]
    first = int(np.argmax(ok))
    y[:first] = y[first]
    return y
