"""Convert TROPOMI sif_red_pc1/pc2 coefficients to SIF at a reporting wavelength."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.linalg import svd


@dataclass(frozen=True)
class SIFBasisWeights:
    """Per detector-row weights for sif_red_pc1 and sif_red_pc2 at one wavelength."""

    ground_pixels: np.ndarray  # 0-based detector rows, shape (n_rows,)
    row_lookup: np.ndarray  # index by 0-based gp → basis row or -1
    weights: np.ndarray  # (n_rows, 2)


def _gaussian_sif_basis(
    wavelength: np.ndarray,
    center: float = 683.0,
    fwhm: float = 25.0,
) -> np.ndarray:
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    values = np.exp(-0.5 * ((wavelength - center) / sigma) ** 2)
    peak = np.max(values)
    return values / peak if peak > 0 else values


def shifted_gaussian_sif_subspace(
    wavelength: np.ndarray,
    *,
    center: float = 683.0,
    fwhm: float = 25.0,
    maximum_shift: float = 2.0,
    shift_step: float = 0.1,
    components: int = 2,
) -> np.ndarray:
    """Port of Retrieval.ShiftedGaussianSIFSubspace / sif_basis (Julia)."""
    shifts = np.arange(-maximum_shift, maximum_shift + shift_step * 0.5, shift_step)
    profiles = np.column_stack(
        [_gaussian_sif_basis(wavelength, center + shift, fwhm) for shift in shifts]
    )
    _, _, vt = svd(profiles.T, full_matrices=False)
    basis = vt[:components].T.copy()
    for j in range(basis.shape[1]):
        k = int(np.argmax(np.abs(basis[:, j])))
        basis[:, j] *= np.sign(basis[k, j])
        m = float(np.max(np.abs(basis[:, j])))
        if m > 0:
            basis[:, j] /= m
    return basis


def _interpolation_weights(wavelength: np.ndarray, target: float) -> tuple[int, int, float]:
    wl = np.asarray(wavelength, dtype=np.float64)
    if target < wl[0] or target > wl[-1]:
        raise ValueError(
            f"reporting wavelength {target} nm outside fit window [{wl[0]:.3g}, {wl[-1]:.3g}]"
        )
    upper = int(np.searchsorted(wl, target))
    if upper == 0:
        return 0, 0, 0.0
    lower = upper - 1
    if upper == lower:
        return lower, lower, 0.0
    fraction = (target - wl[lower]) / (wl[upper] - wl[lower])
    return lower, upper, float(fraction)


def load_sif_basis_weights(basis_path: str | Path, target_nm: float = 678.0) -> SIFBasisWeights:
    """Load row-dependent weights for sif_red_pc1/pc2 → SIF at target_nm."""
    basis_path = Path(basis_path)
    with h5py.File(basis_path, "r") as file:
        # Julia stores (n_basis_rows, n_channels); h5py reads (n_channels, n_basis_rows).
        wavelength = np.asarray(file["wavelength_nm"][:], dtype=np.float64)
        ground_pixels = np.asarray(file["ground_pixels"][:], dtype=np.int32)

    if wavelength.shape[0] == len(ground_pixels):
        # (n_rows, n_channels)
        row_wavelengths = [wavelength[row, :] for row in range(len(ground_pixels))]
    else:
        # (n_channels, n_rows)
        row_wavelengths = [wavelength[:, row] for row in range(len(ground_pixels))]

    weights = np.zeros((len(ground_pixels), 2), dtype=np.float64)
    for row, row_wl in enumerate(row_wavelengths):
        row_wl = np.sort(row_wl)
        fluorescence = shifted_gaussian_sif_subspace(row_wl)
        lower, upper, fraction = _interpolation_weights(row_wl, target_nm)
        weights[row, :] = (1.0 - fraction) * fluorescence[lower, :] + fraction * fluorescence[
            upper, :
        ]

    row_lookup = np.full(int(ground_pixels.max()) + 1, -1, dtype=np.int32)
    for basis_row, gp in enumerate(ground_pixels):
        row_lookup[int(gp)] = basis_row

    return SIFBasisWeights(
        ground_pixels=ground_pixels,
        row_lookup=row_lookup,
        weights=weights,
    )


def sif_from_coefficients(
    detector_row: np.ndarray,
    pc1: np.ndarray,
    pc2: np.ndarray,
    basis: SIFBasisWeights,
) -> np.ndarray:
    """Combine sif_red_pc1/pc2 with row-dependent basis weights → SIF [SI]."""
    n = pc1.size
    out = np.full(n, np.nan, dtype=np.float64)
    dr = np.asarray(detector_row, dtype=np.int32)
    for i in range(n):
        if not (np.isfinite(pc1[i]) and np.isfinite(pc2[i])):
            continue
        gp = int(dr[i])
        if gp < 0 or gp >= basis.row_lookup.size:
            continue
        row = int(basis.row_lookup[gp])
        if row < 0:
            continue
        w = basis.weights[row]
        out[i] = w[0] * pc1[i] + w[1] * pc2[i]
    return out
