"""Compare TOA spectra across Rayleigh optical-depth scales (output of generate_rayleigh_scaling_spectra.jl).

    python plot_rayleigh_scaling_spectra.py [path/to/rt_toa_rayleigh_scaling_*.nc] [-o OUTPUT_NAME] [-z O2-A]

Panels: (a) OCI TOA radiance, (b) radiance relative to scale 1, (c) ΔSIF = TOA SIF signal − water-leaving
SIF, (d) ΔSIF / water-leaving SIF, (e) RT-grid radiance in the --zoom band (wavelength range from ZOOM_BANDS).
The figure goes next to the NetCDF unless OUTPUT_NAME contains a directory; with no name it is <nc stem>.png.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

HERE = Path(__file__).absolute().parent
DEFAULT_NC = HERE / "output_rayleigh_scaling" / "rt_toa_rayleigh_scaling_ocean_lambertian_1010_wavnum0p05_NoRayl.nc"
DEFAULT_OUTPUT_NAME = "sunlight_on.png"   # None → <nc stem>.png
UNITS = r"W m$^{-2}$ sr$^{-1}$ µm$^{-1}$"
# RT-grid zoom windows for panel (e): key → (title label, wavelength range in nm)
ZOOM_BANDS = {
    "O2-B": ("O$_2$-B band", (660.0, 710.0)),
    "O2-A": ("O$_2$-A band", (755.0, 772.0)),
    "H2O": ("H$_2$O band", (715.0, 735.0)),
    "SIF-peak": ("SIF red peak", (670.0, 690.0)),
}
DEFAULT_ZOOM = "O2-B"


def output_path(nc_path, output_name):
    if not output_name:
        return Path(nc_path).with_suffix(".png")
    out = Path(output_name)
    return out if out.parent != Path(".") else Path(nc_path).parent / out


def main(path, output_name=DEFAULT_OUTPUT_NAME, zoom=DEFAULT_ZOOM):
    zoom_label, zoom_nm = ZOOM_BANDS[zoom]
    ds = xr.open_dataset(path).transpose(..., "scale")   # Julia writes (scale, band) in C order
    scales = ds["rayleigh_scale"].values
    lam = ds["wavelength"].values
    lam_h = ds["wavelength_hres"].values
    i_ref = int(np.argmin(np.abs(scales - 1.0)))
    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.9, len(scales)))
    lab = [f"{s:g}" for s in scales]

    fig = plt.figure(figsize=(13, 12), layout="constrained")
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.9])
    ax_r, ax_rel = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ax_sif, ax_srel = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    ax_o2 = fig.add_subplot(gs[2, :])
    axes = (ax_r, ax_rel, ax_sif, ax_srel, ax_o2)
    ref = ds["radiance_clean"].values[:, i_ref]
    sif_wl = ds["sif_waterleaving"].values
    has_sif = sif_wl > 0.02 * np.nanmax(sif_wl)          # relative ΔSIF is meaningless in the SIF tails
    sif_txt = rf"water-leaving SIF$_{{678}}$ = {ds.attrs['sif_678']:g} {UNITS}"
    for k, c in enumerate(colors):
        rad = ds["radiance_clean"].values[:, k]
        dsif = ds["sif_toa"].values[:, k] - sif_wl
        ax_r.plot(lam, rad, color=c, lw=1.4, label=lab[k])
        ax_rel.plot(lam, 100 * (rad / ref - 1), color=c, lw=1.4)
        ax_sif.plot(lam, dsif, color=c, lw=1.4)
        ax_srel.plot(lam, np.where(has_sif, 100 * dsif / sif_wl, np.nan), color=c, lw=1.4)
        m = (lam_h >= zoom_nm[0]) & (lam_h <= zoom_nm[1])
        ax_o2.plot(lam_h[m], ds["radiance_hres"].values[m, k], color=c, lw=0.5)
    sig = ds["sigma_noise"].values[:, i_ref]
    ax_rel.fill_between(lam, -100 * sig / ref, 100 * sig / ref, color="0.85", lw=0,
                        label=r"OCI noise $\pm1\sigma$ (scale 1)")
    ax_r.set_title("(a) OCI TOA radiance", loc="left")
    ax_r.set_ylabel(f"Radiance ({UNITS})")
    ax_r.legend(title="Rayleigh scale", frameon=False, fontsize=8, title_fontsize=8)
    ax_rel.set_title(f"(b) Relative to Rayleigh scale {scales[i_ref]:g}", loc="left")
    ax_rel.set_ylabel("Radiance difference (%)")
    ax_rel.axhline(0, color="0.5", lw=0.7)
    ax_rel.legend(frameon=False, fontsize=8, loc="lower right")
    ax_sif.set_title(r"(c) $\Delta$SIF = TOA SIF signal − water-leaving SIF" + f"\n({sif_txt})", loc="left")
    ax_sif.set_ylabel(rf"$\Delta$SIF ({UNITS})")
    ax_srel.set_title(r"(d) Relative $\Delta$SIF = $\Delta$SIF / water-leaving SIF" + f"\n({sif_txt})", loc="left")
    ax_srel.set_ylabel(r"$\Delta$SIF / water-leaving SIF (%)")
    for ax in (ax_sif, ax_srel):
        ax.axhline(0, color="0.5", lw=0.7)
        ax.set_xlim(lam[has_sif].min() - 5, lam[has_sif].max() + 5)
    ax_o2.set_title(f"(e) {zoom_label} ({zoom_nm[0]:g}–{zoom_nm[1]:g} nm), RT grid", loc="left")
    ax_o2.set_xlim(*zoom_nm)
    ax_o2.set_ylabel(f"Radiance ({UNITS})")
    for ax in axes:
        ax.set_xlabel("Wavelength (nm)")
        ax.spines[["top", "right"]].set_visible(False)

    aod = ds["aod"].values
    aer = "no aerosol" if not np.isfinite(aod).any() else f"GCHPIO column {ds.attrs['aer_column']}"
    fig.suptitle(f"{Path(ds.attrs['ocean_yaml']).name}: SZA {ds.attrs['sza']:g}°, VZA {ds.attrs['vza']:g}°, "
                 f"VAZ {ds.attrs['vaz']:g}°, albedo {ds.attrs['surface_albedo']:g}, "
                 f"SIF$_{{678}}$ {ds.attrs['sif_678']:g}, {aer}", fontsize=11)
    out = output_path(path, output_name)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200)
    print(f"saved {out}")

    print(f"{'scale':>8} {'τ_R(678)':>9} {'R(678)':>8} {'ΔR(678)%':>9} {'ΔSIF(678)':>10}")
    i678 = int(np.argmin(np.abs(lam - 678.0)))
    tau = ds["rayleigh_od_column"].values[i678]
    for k, s in enumerate(scales):
        r = ds["radiance_clean"].values[i678, k]
        dsif = ds["sif_toa"].values[i678, k] - sif_wl[i678]
        print(f"{s:8g} {s * tau:9.4f} {r:8.3f} {100 * (r / ref[i678] - 1):9.2f} {dsif:10.4f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("nc", nargs="?", default=DEFAULT_NC, help="input NetCDF (default: DEFAULT_NC)")
    ap.add_argument("-o", "--output", default=DEFAULT_OUTPUT_NAME,
                    help="figure file name, or a path (default: DEFAULT_OUTPUT_NAME)")
    ap.add_argument("-z", "--zoom", default=DEFAULT_ZOOM, choices=list(ZOOM_BANDS),
                    help=f"RT-grid zoom band for panel (e) (default: {DEFAULT_ZOOM})")
    args = ap.parse_args()
    main(args.nc, args.output, args.zoom)
