#!/usr/bin/env python3
"""Sample 500 open-ocean GCHP/TOMAS columns for vSmartMOM RT.

Produces layer AOD + wet radii aggregated to the species list used by
`ocean_coxmunk_0912.yaml` / `ocean_column_aerosols.jl`:

  so4  ← SF01–15
  sala ← SS01–13   (dry r ≲ 0.5 μm, GEOS-Chem SALA cut)
  salc ← SS14–15   (dry r ≳ 0.5 μm; TOMAS15's two ×32 coarse bins)
  ocpi ← OCIL + OCOB  (hydrophilic + hydrophobic OC; YAML has one OC slot)
  bcpi ← ECIL + ECOB  (hydrophilic + hydrophobic EC; YAML has one BC slot)
  dust1..dust7 ← DUST bins grouped 01-02, 03-04, …, 13-15

Wet radius follows GEOS-Chem TOMAS GETDP (tomas_mod.F90):
  M_wet = (M_dry + 0.1875·M_SF + M_AW) / N     # kg / particle; 0.1875 = 18/96
  r = 0.5 · (6·M_wet / (π·ρ_wet))^(1/3)       # ρ_wet volume-additive

Layer AOD at λ_ref = 0.55 μm uses Mie Q_ext of the internally mixed
wet particle (volume-average refractive index), then splits τ among YAML
species by dry-volume fraction:

  τ = N · π · r² · Q_ext(x, m) · Δz
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

# Optional cartopy for ocean mask / site map
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.io import shapereader
from shapely.geometry import Point
from shapely.ops import unary_union
from shapely.prepared import prep

NC_PATH = "/kiwi-data/Data/model/GeosChem/GEOSChem.Custom.20190702_0000z.nc4"
OUT_DIR = Path(__file__).resolve().parent / "output_aerosol_profiles"
OUT_NC = OUT_DIR / "gchp_ocean_columns_n500.nc"
OUT_MAP = OUT_DIR / "gchp_ocean_columns_sites.png"

N_SAMPLES = int(os.environ.get("N_SAMPLES", "500"))
SEED = 20260916
LAT_MIN, LAT_MAX = -60.0, 60.0

# SpeciesConcVV is mol/mol dry; Met_AD is kg dry air.
# MWs must match species_database_tomas.yml (the values History used to
# write VV). TOMAS's unused MOLWT table lists SF as 98; History uses 96.
M_AIR = 28.9644e-3  # kg/mol dry air
R_DRY = 287.05      # J/(kg K)
G0 = 9.80665
LAMBDA_REF_UM = 0.550  # YAML scattering.λ_ref / aod_reference_wavelength

# kg/mol — species_database_tomas.yml
SPC_MW = {
    "SF": 96.0e-3,
    "SS": 58.5e-3,
    "ECOB": 12.01e-3,
    "ECIL": 12.01e-3,
    "OCOB": 12.01e-3,
    "OCIL": 12.01e-3,
    "DUST": 100.0e-3,
    "AW": 18.0e-3,
}
# kg/m³ — GETDP/AERODENS hard-coded densities (tomas_mod.F90)
SPC_RHO = {
    "SF": 1500.0,
    "SS": 1500.0,
    "ECOB": 2200.0,
    "ECIL": 2200.0,
    "OCOB": 1400.0,
    "OCIL": 1400.0,
    "DUST": 2650.0,
    "AW": 1000.0,
    "NH4": 1500.0,
}
# Volume-mixing RI at 550 nm (OPAC / typical GEOS-Chem optics).
# Imaginary part ≥ 0 (Bohren-Huffman n + iκ convention).
SPC_N = {
    "SF": 1.53 + 0.0j,
    "SS": 1.50 + 0.0j,
    "ECOB": 1.75 + 0.44j,
    "ECIL": 1.75 + 0.44j,
    "OCOB": 1.53 + 0.006j,
    "OCIL": 1.53 + 0.006j,
    "DUST": 1.53 + 0.008j,
    "AW": 1.333 + 0.0j,
    "NH4": 1.53 + 0.0j,
}

DRY_PREFS = ["SF", "SS", "ECOB", "ECIL", "OCOB", "OCIL", "DUST"]

# YAML species → list of (TOMAS prefix, 1-based bins)
SPECIES_GROUPS = {
    "so4": [("SF", list(range(1, 16)))],
    "sala": [("SS", list(range(1, 14)))],   # dry r ≲ 0.5 μm
    "salc": [("SS", list(range(14, 16)))],  # dry r ≳ 0.5 μm
    "ocpi": [("OCIL", list(range(1, 16))), ("OCOB", list(range(1, 16)))],
    "bcpi": [("ECIL", list(range(1, 16))), ("ECOB", list(range(1, 16)))],
    "dust1": [("DUST", list(range(1, 3)))],
    "dust2": [("DUST", list(range(3, 5)))],
    "dust3": [("DUST", list(range(5, 7)))],
    "dust4": [("DUST", list(range(7, 9)))],
    "dust5": [("DUST", list(range(9, 11)))],
    "dust6": [("DUST", list(range(11, 13)))],
    "dust7": [("DUST", list(range(13, 16)))],
}


def _read_vv(ds, prefix: str, ibin: int) -> np.ndarray:
    name = f"SpeciesConcVV_{prefix}{ibin:02d}"
    return np.asarray(ds[name][0], dtype=np.float64)  # (lev, nf, Y, X)


def number_cm3(nk_vv, air_mol, airvol_m3):
    """TOMAS NK SpeciesConcVV → #/cm³.

    NK is stored as a pseudo-VMR with MW_g = 1, so moles_NK = VV · n_air
    equals 1000 × particle number (1 particle ≡ 1 'fake kg').
    """
    vol_cm3 = airvol_m3 * 1e6
    return (nk_vv / 1000.0) * air_mol / np.maximum(vol_cm3, 1e-30)


def mass_ug_m3(vv, air_mol, airvol_m3, mw_kg_mol):
    """SpeciesConcVV (mol/mol dry) → μg/m³ using the History molar mass."""
    mol_spc = vv * air_mol
    return (mol_spc * (mw_kg_mol * 1e9)) / np.maximum(airvol_m3, 1e-30)


def _mie_qext_scalar(x: float, m: complex) -> float:
    """Bohren-Huffman Mie extinction efficiency for a sphere.

    ``x = 2πr/λ``, ``m = n + iκ`` with κ ≥ 0. Rayleigh for x < 0.05.
    """
    x = float(x)
    m = complex(m)
    if x <= 0.0 or not np.isfinite(x):
        return 0.0
    if x < 0.05:
        chi = (m * m - 1.0) / (m * m + 2.0)
        qabs = 4.0 * x * chi.imag
        qsca = (8.0 / 3.0) * x**4 * abs(chi) ** 2
        return float(max(qabs + qsca, 0.0))
    nstop = int(x + 4.0 * x ** (1.0 / 3.0) + 2.0)
    z = m * x
    nmx = int(max(nstop, abs(z)) + 15)
    d = np.zeros(nmx + 1, dtype=np.complex128)
    for n in range(nmx, 0, -1):
        d[n - 1] = n / z - 1.0 / (d[n] + n / z)
    psi0, psi1 = np.cos(x), np.sin(x)
    chi0, chi1 = np.sin(x), -np.cos(x)
    xi1 = psi1 + 1j * chi1
    qext = 0.0
    for n in range(1, nstop + 1):
        psi = (2 * n - 1) * psi1 / x - psi0
        chi = (2 * n - 1) * chi1 / x - chi0
        xi = psi + 1j * chi
        an = ((d[n] / m + n / x) * psi - psi1) / ((d[n] / m + n / x) * xi - xi1)
        bn = ((m * d[n] + n / x) * psi - psi1) / ((m * d[n] + n / x) * xi - xi1)
        qext += (2 * n + 1) * (an.real + bn.real)
        psi0, psi1 = psi1, psi
        chi0, chi1 = chi1, chi
        xi1 = psi1 + 1j * chi1
    return float(max((2.0 / x**2) * qext, 0.0))


def mie_qext(x, m):
    """Vectorized wrapper around `_mie_qext_scalar`."""
    x_arr = np.asarray(x, dtype=np.float64)
    m_arr = np.asarray(m, dtype=np.complex128)
    if x_arr.shape == ():
        return _mie_qext_scalar(float(x_arr), complex(np.broadcast_to(m_arr, (1,))[0]))
    m_arr = np.broadcast_to(m_arr, x_arr.shape)
    out = np.empty(x_arr.shape, dtype=np.float64)
    for idx in np.ndindex(x_arr.shape):
        out[idx] = _mie_qext_scalar(x_arr[idx], m_arr[idx])
    return out


def wet_particle(N_cm3, masses_ug_m3):
    """GETDP monodisperse wet radius [μm] and volume-mix n_eff for one bin."""
    N_m3 = np.maximum(N_cm3, 0.0) * 1e6
    m = {k: np.maximum(v, 0.0) * 1e-9 for k, v in masses_ug_m3.items()}  # kg/m³
    M_AW = m.get("AW", 0.0)
    M_SF = m.get("SF", 0.0)
    M_NH4 = 0.1875 * M_SF
    dry_keys = [k for k in m if k != "AW"]
    M_dry = sum(m[k] for k in dry_keys) if dry_keys else 0.0
    V = 0.0
    nV = 0.0 + 0.0j
    M_tot = M_dry + M_NH4 + M_AW
    for k in dry_keys:
        v_k = m[k] / SPC_RHO[k]
        V = V + v_k
        nV = nV + v_k * SPC_N[k]
    v_nh4 = M_NH4 / SPC_RHO["NH4"]
    v_aw = M_AW / SPC_RHO["AW"]
    V = V + v_nh4 + v_aw
    nV = nV + v_nh4 * SPC_N["NH4"] + v_aw * SPC_N["AW"]
    rho_wet = np.where(V > 0, M_tot / np.maximum(V, 1e-40), 1500.0)
    n_eff = np.asarray(np.where(V > 0, nV / np.maximum(V, 1e-40), 1.5 + 0.0j), dtype=np.complex128)
    n_eff = n_eff.real + 1j * np.maximum(n_eff.imag, 0.0)
    ok = (N_m3 > 1e-6) & (M_dry > 1e-20) & (M_tot > 0)
    M_part = np.where(ok, M_tot / N_m3, 0.0)
    d_m = np.where(ok, (6.0 * M_part / (np.pi * np.maximum(rho_wet, 1.0))) ** (1.0 / 3.0), 0.0)
    r_um = 0.5 * d_m * 1e6
    r_um = np.where(r_um > 100.0, np.nan, r_um)
    r_um = np.where(ok, r_um, np.nan)
    n_eff = np.where(ok, n_eff, 1.5 + 0.0j)
    return r_um, n_eff


def dry_volume_conc(masses_ug_m3):
    """Dry-species volume concentration (m³/m³) per prefix, shape matches mass."""
    return {k: np.maximum(masses_ug_m3[k], 0.0) * 1e-9 / SPC_RHO[k] for k in DRY_PREFS}


def layer_tau(N_cm3, r_um, qext, dz_m):
    r_m = np.nan_to_num(r_um, nan=0.0) * 1e-6
    N_m3 = np.maximum(N_cm3, 0.0) * 1e6
    q = np.nan_to_num(np.asarray(qext, dtype=np.float64), nan=0.0)
    return N_m3 * (np.pi * r_m**2) * np.maximum(q, 0.0) * np.maximum(dz_m, 0.0)


def hydrostatic_edges_and_dz(PS, DELP, T, SPHU):
    """Wet-air hydrostatic p_edge, p_mid [hPa] and Δz [m] (BOA-first).

    Uses Tv = T (1 + 0.61 q) and Δz = (Rd Tv / g) ln(p_bot / p_top),
    matching GCHPIO `_layer_thickness_m`. ``SPHU`` is g kg⁻¹.
    """
    nlev = DELP.shape[0]
    p_edge = np.empty((nlev + 1,) + DELP.shape[1:], dtype=np.float64)
    p_edge[0] = PS
    np.cumsum(DELP, axis=0, out=p_edge[1:])
    p_edge[1:] = PS - p_edge[1:]
    # DELP can slightly overshoot PS at model top → clip and keep BOA→TOA decreasing.
    p_edge = np.maximum(p_edge, 1e-4)
    p_edge = np.minimum.accumulate(p_edge, axis=0)
    p_mid = 0.5 * (p_edge[:-1] + p_edge[1:])
    q = np.maximum(SPHU, 0.0) * 1e-3  # kg/kg
    Tv = T * (1.0 + 0.61 * q)
    p_bot = p_edge[:-1]
    p_top = p_edge[1:]
    with np.errstate(divide="ignore", invalid="ignore"):
        dz = np.where(
            p_bot > p_top,
            (R_DRY * Tv / G0) * np.log(p_bot / p_top),
            0.0,
        )
    dz = np.maximum(dz, 0.0)
    return p_mid, dz


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Reading {NC_PATH}")
    ds = Dataset(NC_PATH)
    lons = np.asarray(ds["lons"][:], dtype=np.float64)
    lats = np.asarray(ds["lats"][:], dtype=np.float64)
    AD = np.asarray(ds["Met_AD"][0], dtype=np.float64)
    AIRVOL = np.asarray(ds["Met_AIRVOL"][0], dtype=np.float64)
    DELP = np.asarray(ds["Met_DELP"][0], dtype=np.float64)
    T = np.asarray(ds["Met_T"][0], dtype=np.float64)
    PS = np.asarray(ds["Met_PS1WET"][0], dtype=np.float64)
    SPHU = np.asarray(ds["Met_SPHU"][0], dtype=np.float64)
    nlev, nf, ny, nx = AD.shape
    air_mol = AD / M_AIR
    p_mid_grid, dz = hydrostatic_edges_and_dz(PS, DELP, T, SPHU)

    # Reference p_mid profile (grid-mean) for NC lev axis
    p_mid_ref = np.mean(p_mid_grid.reshape(nlev, -1), axis=1)

    print("Building ocean mask…")
    land_shp = shapereader.natural_earth(resolution="110m", category="physical", name="land")
    land = prep(unary_union(list(shapereader.Reader(land_shp).geometries())))
    flat_lon = lons.ravel()
    flat_lat = lats.ravel()
    is_ocean = np.zeros(flat_lon.size, dtype=bool)
    for i, (lo, la) in enumerate(zip(flat_lon, flat_lat)):
        if LAT_MIN <= la <= LAT_MAX and not land.contains(Point(float(lo), float(la))):
            is_ocean[i] = True
    ocean_flat = np.flatnonzero(is_ocean)
    print(f"Ocean cells in [{LAT_MIN},{LAT_MAX}]: {ocean_flat.size} / {flat_lon.size}")

    # Prefer columns with some sea-salt mass
    print("Sea-salt column mass for sampling…")
    ss_col = np.zeros((nf, ny, nx), dtype=np.float64)
    for ib in range(1, 16):
        vv = _read_vv(ds, "SS", ib)
        ss_col += mass_ug_m3(vv, air_mol, AIRVOL, SPC_MW["SS"]).sum(axis=0)
    ss_flat = ss_col.ravel()
    cand = ocean_flat[ss_flat[ocean_flat] > 1e-3]
    if cand.size < N_SAMPLES:
        cand = ocean_flat
    rng = np.random.default_rng(SEED)
    picks_flat = rng.choice(cand, size=N_SAMPLES, replace=False)
    picks = [np.unravel_index(int(i), (nf, ny, nx)) for i in picks_flat]
    print(f"Sampled {N_SAMPLES} columns")

    # Load all needed tracers once
    print("Loading TOMAS tracers…")
    prefixes = ["NK", "SF", "SS", "ECOB", "ECIL", "OCOB", "OCIL", "DUST", "AW"]
    vv = {p: [_read_vv(ds, p, ib) for ib in range(1, 16)] for p in prefixes}
    ds.close()

    # Output arrays (lev, sample)
    aod_out = {name: np.zeros((nlev, N_SAMPLES), dtype=np.float32) for name in SPECIES_GROUPS}
    rad_out = {name: np.zeros((nlev, N_SAMPLES), dtype=np.float32) for name in SPECIES_GROUPS}
    p_mid_s = np.zeros((nlev, N_SAMPLES), dtype=np.float32)
    lat_s = np.zeros(N_SAMPLES, dtype=np.float32)
    lon_s = np.zeros(N_SAMPLES, dtype=np.float32)
    i_face = np.zeros(N_SAMPLES, dtype=np.int32)
    i_y = np.zeros(N_SAMPLES, dtype=np.int32)
    i_x = np.zeros(N_SAMPLES, dtype=np.int32)

    def cell_fields(f, y, x):
        """Return N_cm3[bin], mass_ug[bin,spc], r_um[bin], tau_bin[bin] at one column."""
        am = air_mol[:, f, y, x]
        av = AIRVOL[:, f, y, x]
        dzz = dz[:, f, y, x]
        N = np.stack([number_cm3(vv["NK"][ib][:, f, y, x], am, av) for ib in range(15)], axis=0)
        mass = {}
        for pref in DRY_PREFS + ["AW"]:
            mass[pref] = np.stack(
                [mass_ug_m3(vv[pref][ib][:, f, y, x], am, av, SPC_MW[pref]) for ib in range(15)],
                axis=0,
            )
        r = np.full((15, nlev), np.nan)
        tau = np.zeros((15, nlev))
        for ib in range(15):
            masses = {k: mass[k][ib] for k in mass}
            r[ib], n_eff = wet_particle(N[ib], masses)
            x_mie = np.where(np.isfinite(r[ib]), 2.0 * np.pi * r[ib] / LAMBDA_REF_UM, 0.0)
            qext = mie_qext(x_mie, n_eff)
            tau[ib] = layer_tau(N[ib], r[ib], qext, dzz)
        return N, mass, r, tau

    print("Computing AOD + wet radius for sampled columns…")
    for s, (f, y, x) in enumerate(picks):
        i_face[s], i_y[s], i_x[s] = f, y, x
        lat_s[s] = lats[f, y, x]
        lon_s[s] = lons[f, y, x]
        p_mid_s[:, s] = p_mid_grid[:, f, y, x]
        _, mass, r, tau = cell_fields(f, y, x)
        # Dry-volume concentration per bin for apportioning mixed-particle τ.
        V_pref = dry_volume_conc(mass)  # prefix → (15, lev)
        V_dry = sum(V_pref[p] for p in DRY_PREFS)
        for name, parts in SPECIES_GROUPS.items():
            tau_g = np.zeros(nlev)
            w_r = np.zeros((0, nlev))
            r_stack = np.zeros((0, nlev))
            for pref, bins in parts:
                ibs = [b - 1 for b in bins]
                frac = np.divide(
                    V_pref[pref][ibs],
                    np.maximum(V_dry[ibs], 1e-40),
                    out=np.zeros_like(V_pref[pref][ibs]),
                    where=V_dry[ibs] > 0,
                )
                tau_g = tau_g + (tau[ibs] * frac).sum(axis=0)
                w_r = np.vstack([w_r, np.maximum(mass[pref][ibs], 0.0)])
                r_stack = np.vstack([r_stack, r[ibs]])
            wsum = w_r.sum(axis=0)
            r_mean = np.full(nlev, np.nan)
            for k in range(nlev):
                if wsum[k] > 0:
                    ww = w_r[:, k]
                    rr = r_stack[:, k]
                    ok = np.isfinite(rr) & (ww > 0)
                    if np.any(ok):
                        r_mean[k] = np.average(rr[ok], weights=ww[ok])
            aod_out[name][:, s] = tau_g.astype(np.float32)
            rad_out[name][:, s] = np.nan_to_num(r_mean, nan=0.0).astype(np.float32)
        if (s + 1) % 50 == 0:
            print(f"  {s+1}/{N_SAMPLES}")

    # Column diagnostics
    tau_ref = np.zeros(N_SAMPLES, dtype=np.float32)
    ss_col = np.zeros(N_SAMPLES, dtype=np.float32)
    for name, arr in aod_out.items():
        tau_ref += arr.sum(axis=0)
        if name in ("sala", "salc"):
            ss_col += arr.sum(axis=0)
    ss_frac = np.divide(ss_col, tau_ref, out=np.full(N_SAMPLES, np.nan, dtype=np.float32), where=tau_ref > 0)

    if OUT_NC.exists():
        OUT_NC.unlink()
    print(f"Writing {OUT_NC}")
    dout = Dataset(OUT_NC, "w")
    dout.createDimension("sample", N_SAMPLES)
    dout.createDimension("lev", nlev)

    def put(name, data, dims, **attrs):
        arr = np.asarray(data)
        dtype = np.int32 if np.issubdtype(arr.dtype, np.integer) else np.float32
        v = dout.createVariable(name, dtype, dims)
        v[:] = arr
        for k, val in attrs.items():
            v.setncattr(k, val)

    put("lat", lat_s, ("sample",), units="degrees_north")
    put("lon", lon_s, ("sample",), units="degrees_east")
    put("i_face", i_face, ("sample",), long_name="GCHP cubed-sphere face index 0-based")
    put("i_y", i_y, ("sample",), long_name="Ydim index 0-based")
    put("i_x", i_x, ("sample",), long_name="Xdim index 0-based")
    put("lev", np.arange(1, nlev + 1, dtype=np.float32), ("lev",), units="1",
        long_name="GCHP level (1=BOA, 72=TOA)")
    put("p_mid", p_mid_s, ("lev", "sample"), units="hPa",
        long_name="layer mid pressure (surface-first)")
    put("p_mid_ref", p_mid_ref.astype(np.float32), ("lev",), units="hPa",
        long_name="grid-mean mid pressure (for plotting)")

    put("tau_ref", tau_ref, ("sample",), units="1",
        long_name="column total AOD at 550 nm (Mie Qext, all YAML species)")
    put("ss_col", ss_col, ("sample",), units="1",
        long_name="column sea-salt AOD at 550 nm (sala+salc)")
    put("ss_fraction", ss_frac, ("sample",), units="1")

    # Names expected by updated ocean_coxmunk_0912.yaml
    aod_name = {
        "so4": "AOD_GCHP_SO4",
        "sala": "AOD_GCHP_SALA",
        "salc": "AOD_GCHP_SALC",
        "ocpi": "AOD_GCHP_OCPI",
        "bcpi": "AOD_GCHP_BCPI",
        "dust1": "AOD_GCHP_DUST1",
        "dust2": "AOD_GCHP_DUST2",
        "dust3": "AOD_GCHP_DUST3",
        "dust4": "AOD_GCHP_DUST4",
        "dust5": "AOD_GCHP_DUST5",
        "dust6": "AOD_GCHP_DUST6",
        "dust7": "AOD_GCHP_DUST7",
    }
    rad_name = {k: v.replace("AOD_GCHP_", "RadiWet_GCHP_") for k, v in aod_name.items()}

    for key, an in aod_name.items():
        put(an, aod_out[key], ("lev", "sample"), units="1",
            long_name=f"layer AOD at 550 nm for {key} from GCHP/TOMAS")
        put(rad_name[key], rad_out[key], ("lev", "sample"), units="um",
            long_name=f"mass-weighted wet radius for {key}")

    put("AOD_seasalt", aod_out["sala"] + aod_out["salc"], ("lev", "sample"), units="1")
    put("AOD_total_layer", sum(aod_out.values()), ("lev", "sample"), units="1")

    dout.title = "GCHP/TOMAS open-ocean aerosol columns for vSmartMOM RT"
    dout.source = NC_PATH
    dout.method = (
        "GETDP wet radius (volume-additive ρ); Mie Qext at 550 nm; "
        "dry-volume AOD split; History MW (SF=96)"
    )
    dout.n_samples = np.int32(N_SAMPLES)
    dout.seed = np.int32(SEED)
    dout.lat_min = np.float32(LAT_MIN)
    dout.lat_max = np.float32(LAT_MAX)
    dout.lambda_ref_um = np.float32(LAMBDA_REF_UM)
    group_str = ",".join(
        f"{k}:" + "+".join(f"{p}{list(b)}" for p, b in parts)
        for k, parts in SPECIES_GROUPS.items()
    )
    dout.species_groups = group_str
    dout.close()
    print(f"Wrote {OUT_NC}")
    print(f"  tau_ref: {tau_ref.min():.4f} .. {tau_ref.max():.4f}")
    print(f"  ss_col:  {ss_col.min():.4f} .. {ss_col.max():.4f}")

    fig = plt.figure(figsize=(11, 5.2))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.coastlines(linewidth=0.6)
    vmax = float(np.percentile(ss_col, 99)) if ss_col.size else 1.0
    sc = ax.scatter(
        lon_s, lat_s, c=ss_col, s=12, cmap="Blues",
        transform=ccrs.PlateCarree(), vmin=0, vmax=max(vmax, 1e-6), zorder=3,
    )
    ax.set_global()
    ax.set_title(f"{N_SAMPLES} GCHP open-ocean columns (colored by sea-salt AOD at 550 nm)")
    fig.colorbar(sc, ax=ax, pad=0.02, shrink=0.85, label="Sea-salt AOD (550 nm)")
    fig.tight_layout()
    fig.savefig(OUT_MAP, dpi=150, bbox_inches="tight")
    print("Wrote", OUT_MAP)
    plt.close(fig)


if __name__ == "__main__":
    main()
