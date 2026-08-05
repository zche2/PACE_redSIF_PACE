"""Write match-up NetCDF product."""

from __future__ import annotations

import json
from pathlib import Path

import netCDF4 as nc
import numpy as np

from .config import Config


def write_matchup_nc(
    path: Path,
    cfg: Config,
    matches: dict,
    sim: dict,
    *,
    fingerprints: dict[str, str],
    svd_meta: dict | None = None,
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()

    n_match = int(matches["dist_km"].size)
    oci_bands = np.asarray(sim["oci_bands"], dtype=np.float32)
    n_band = int(oci_bands.size)

    with nc.Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("match", n_match)
        ds.createDimension("band", n_band)
        ds.createDimension("strlen", 256)

        def _f(name, dims, data, units=None, long_name=None, dtype="f4"):
            v = ds.createVariable(name, dtype, dims, zlib=True, complevel=4)
            v[:] = data
            if units:
                v.units = units
            if long_name:
                v.long_name = long_name
            return v

        band = ds.createVariable("band", "f4", ("band",))
        band[:] = oci_bands
        band.units = "nm"
        band.long_name = "OCI band center wavelength"

        match = ds.createVariable("match", "i4", ("match",))
        match[:] = np.arange(n_match, dtype=np.int32)

        _f(
            "lt_oci_sim",
            ("match", "band"),
            sim["lt_oci"],
            units="W m-2 sr-1 um-1",
            long_name="Simulated OCI radiance from TROPOMI (denoised+convolved, no noise)",
        )
        _f(
            "lt_oci_sim_noisy",
            ("match", "band"),
            sim["lt_oci_noisy"],
            units="W m-2 sr-1 um-1",
            long_name="Simulated OCI radiance with PACE SNR white noise",
        )
        _f(
            "lt_oci_meas",
            ("match", "band"),
            sim["lt_oci_meas"],
            units="W m-2 sr-1 um-1",
            long_name="Measured PACE OCI L1B radiance at match location",
        )
        _f(
            "sigma",
            ("match", "band"),
            sim["sigma"],
            units="W m-2 sr-1 um-1",
            long_name="Noise std-dev from PACE SNR LUT",
        )

        _f("lat_pace", ("match",), matches["pace_lat"], "degrees_north", "PACE latitude")
        _f("lon_pace", ("match",), matches["pace_lon"], "degrees_east", "PACE longitude")
        _f("lat_trop", ("match",), matches["trop_lat"], "degrees_north", "TROPOMI latitude")
        _f("lon_trop", ("match",), matches["trop_lon"], "degrees_east", "TROPOMI longitude")
        _f("sza_pace", ("match",), matches["pace_sza"], "degree", "PACE solar zenith")
        _f("vza_pace", ("match",), matches["pace_vza"], "degree", "PACE sensor zenith")
        _f("sza_trop", ("match",), matches["trop_sza"], "degree", "TROPOMI solar zenith")
        _f("vza_trop", ("match",), matches["trop_vza"], "degree", "TROPOMI viewing zenith")
        _f(
            "time_pace",
            ("match",),
            matches["time_pace"],
            "seconds since 1970-01-01 00:00:00 UTC",
            "PACE scan-line time",
            dtype="f8",
        )
        _f(
            "time_trop",
            ("match",),
            matches["time_trop"],
            "seconds since 1970-01-01 00:00:00 UTC",
            "TROPOMI scan-line time",
            dtype="f8",
        )
        _f("dt_min", ("match",), matches["dt_min"], "min", "TROPOMI−PACE time difference")
        _f("dist_km", ("match",), matches["dist_km"], "km", "Great-circle distance")
        _f(
            "mu_trop_over_mu_oci",
            ("match",),
            sim["mu_trop_over_mu_oci"],
            "1",
            "cos(SZA_trop)/cos(SZA_pace)",
        )

        for name, key in [
            ("pace_scan", "pace_scan"),
            ("pace_pix", "pace_pix"),
            ("trop_scan", "trop_scan"),
            ("trop_pix", "trop_pix"),
            ("swath_id", "swath_id"),
            ("trop_gp", None),
        ]:
            data = sim["trop_gp"] if key is None else matches[key]
            v = ds.createVariable(name, "i4", ("match",), zlib=True)
            v[:] = np.asarray(data, dtype=np.int32)

        def _str_var(name: str, values) -> None:
            v = ds.createVariable(name, "S1", ("match", "strlen"))
            arr = np.zeros((n_match, 256), dtype="S1")
            for i, s in enumerate(values):
                b = str(s).encode("utf-8")[:255]
                arr[i, : len(b)] = np.frombuffer(b, dtype="S1")
            v[:] = arr

        _str_var("pace_path", matches["pace_path"])
        _str_var("tropomi_path", matches["tropomi_path"])
        _str_var("pace_name", matches["pace_name"])

        ds.title = "TROPOMI BD5 → PACE OCI co-located match-ups"
        ds.Conventions = "CF-1.8"
        ds.source = "TROPOMI_OCI_colocation pipeline"
        ds.n_pc = np.int32(sim["n_pc"])
        ds.add_noise = np.int32(1 if sim["add_noise"] else 0)
        ds.noise_seed = np.int32(sim["noise_seed"])
        ds.comment = (
            "lt units: W m-2 sr-1 um-1 numerically equal to mW m-2 nm-1 sr-1. "
            "SVD trained with loose TROPOMI filters; out-of-sample denoise via Vt projection. "
            "No BRDF/VZA path correction beyond mu0 amplitude ratio."
        )
        ds.config_json = json.dumps(cfg.to_plain_dict())
        ds.fingerprints_json = json.dumps(fingerprints)
        if svd_meta:
            ds.svd_n_train = np.int32(svd_meta.get("n_train", -1))
            ds.svd_sza_max_deg = np.float32(cfg.svd.sza_max_deg)
            ds.svd_vza_max_deg = np.float32(cfg.svd.vza_max_deg)
            ds.svd_lt_max = np.float32(cfg.svd.lt_max)

    return path
