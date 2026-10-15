"""TOA spectra by SIF intensity from the Rayleigh × aerosol × SIF × view grid (generate_rayleigh_scaling_spectra.jl).

    python plot_rayleigh_scaling_grid.py [path/to/*_grid.nc] [--vza 30] [--rayleigh 1] [--aerosol 1] [--band 678]

Figure 1 (<nc stem>_by_sif_spectra.png), one scene picked by --vza/--rayleigh/--aerosol (nearest grid value):
  (a) OCI TOA radiance, (b) TOA SIF signal vs water-leaving SIF, (c) ΔSIF = TOA − water-leaving SIF,
  (d) TOA SIF signal / OCI noise σ, all coloured by SIF intensity.
Figure 2 (<nc stem>_by_sif_<band>nm.png), every scene at --band: ΔSIF and SIF/σ against SIF intensity,
  one column per view, colour = Rayleigh scale, line style = aerosol scale.
The RT is linear in SIF, so ΔSIF ∝ SIF intensity and ΔSIF / water-leaving SIF does not depend on it.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

HERE = Path(__file__).absolute().parent
DEFAULT_NC = (HERE / "output_rayleigh_scaling" /
              "rt_toa_rayleigh_scaling_ocean_lambertian_1010_wavnum0p05_NoRayl_sunlight_on_AerColumn_100_grid.nc")
UNITS = r"W m$^{-2}$ sr$^{-1}$ µm$^{-1}$"
LINESTYLES = ["-", "--", ":", "-."]


def nearest(values, x):
    return int(np.argmin(np.abs(np.asarray(values) - x)))


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)


def plot_spectra(ds, iv, ir, ia, out):
    sel = dict(view=iv, rayleigh=ir, aerosol=ia)
    lam = ds["wavelength"].values
    sif = ds["sif_678"].values
    rad = ds["radiance_clean"].isel(sel).transpose("band", "sif").values
    rad0 = ds["radiance_nosif"].isel(sel).values
    toa = ds["sif_toa"].isel(sel).transpose("band", "sif").values
    sig = ds["sigma_noise"].isel(sel).transpose("band", "sif").values
    wl = ds["sif_waterleaving"].transpose("band", "sif").values
    has_sif = wl[:, -1] > 0.02 * np.nanmax(wl[:, -1])
    xlim = (lam[has_sif].min() - 5, lam[has_sif].max() + 5)
    colors = plt.get_cmap("plasma")(np.linspace(0.05, 0.85, len(sif)))

    fig, axs = plt.subplots(2, 2, figsize=(13, 9), layout="constrained")
    (ax_r, ax_s), (ax_d, ax_n) = axs
    ax_r.plot(lam, rad0, color="k", lw=1.0, ls="--", label="no SIF")
    for k, c in enumerate(colors):
        ax_r.plot(lam, rad[:, k], color=c, lw=1.2, label=f"{sif[k]:g}")
        ax_s.plot(lam, toa[:, k], color=c, lw=1.4)
        ax_s.plot(lam, wl[:, k], color=c, lw=1.0, ls="--")
        ax_d.plot(lam, toa[:, k] - wl[:, k], color=c, lw=1.4)
        ax_n.plot(lam, toa[:, k] / sig[:, k], color=c, lw=1.4)
    ax_r.set_title("(a) OCI TOA radiance", loc="left")
    ax_r.set_ylabel(f"Radiance ({UNITS})")
    ax_r.legend(title=r"SIF$_{678}$", frameon=False, fontsize=8, title_fontsize=8)
    ax_s.set_title("(b) TOA SIF signal (solid) and water-leaving SIF (dashed)", loc="left")
    ax_s.set_ylabel(f"SIF ({UNITS})")
    ax_d.set_title(r"(c) $\Delta$SIF = TOA SIF signal − water-leaving SIF", loc="left")
    ax_d.set_ylabel(rf"$\Delta$SIF ({UNITS})")
    ax_n.set_title(r"(d) TOA SIF signal / OCI noise $\sigma$", loc="left")
    ax_n.set_ylabel(r"SIF / $\sigma$")
    ax_n.axhline(1, color="0.5", lw=0.7, ls=":")
    ax_d.axhline(0, color="0.5", lw=0.7)
    for ax in (ax_s, ax_d, ax_n):
        ax.set_xlim(*xlim)
    for ax in axs.flat:
        ax.set_xlabel("Wavelength (nm)")
        style(ax)
    fig.suptitle(f"{scene_text(ds)}\nVZA {ds['vza'].values[iv]:g}°, VAZ {ds['vaz'].values[iv]:g}°, "
                 f"Rayleigh ×{ds['rayleigh_scale'].values[ir]:g}, aerosol ×{ds['aerosol_scale'].values[ia]:g}",
                 fontsize=11)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"saved {out}")


def plot_band(ds, band_nm, out):
    ib = nearest(ds["wavelength"].values, band_nm)
    lam_b = ds["wavelength"].values[ib]
    sif = ds["sif_678"].values
    order = ("sif", "view", "rayleigh", "aerosol")
    toa = ds["sif_toa"].isel(band=ib).transpose(*order).values
    sig = ds["sigma_noise"].isel(band=ib).transpose(*order).values
    wl = ds["sif_waterleaving"].isel(band=ib).values
    dsif = toa - wl[:, None, None, None]
    rel = 100 * dsif[-1] / wl[-1]                                   # (view, rayleigh, aerosol)
    ray, aer = ds["rayleigh_scale"].values, ds["aerosol_scale"].values
    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.85, len(ray)))
    n_v = ds.sizes["view"]

    fig, axs = plt.subplots(2, n_v, figsize=(4.6 * n_v, 8), layout="constrained", sharex=True, sharey="row")
    axs = axs.reshape(2, n_v)
    for iv in range(n_v):
        ax_d, ax_n = axs[:, iv]
        for ir, c in enumerate(colors):
            for ia, a in enumerate(aer):
                ls = LINESTYLES[ia % len(LINESTYLES)]
                ax_d.plot(sif, dsif[:, iv, ir, ia], color=c, ls=ls, marker="o", ms=3, lw=1.3,
                          label=f"R×{ray[ir]:g}, A×{a:g} ({rel[iv, ir, ia]:+.2f}%)")
                ax_n.plot(sif, toa[:, iv, ir, ia] / sig[:, iv, ir, ia], color=c, ls=ls, marker="o", ms=3, lw=1.3)
        ax_d.set_title(f"VZA {ds['vza'].values[iv]:g}°, VAZ {ds['vaz'].values[iv]:g}°", loc="left")
        ax_d.axhline(0, color="0.5", lw=0.7)
        ax_n.axhline(1, color="0.5", lw=0.7, ls=":")
        ax_d.legend(title=r"Rayleigh, aerosol (relative $\Delta$SIF)", frameon=False, fontsize=7, title_fontsize=7)
        ax_n.set_xlabel(rf"Water-leaving SIF$_{{678}}$ ({UNITS})")
        for ax in (ax_d, ax_n):
            style(ax)
    axs[0, 0].set_ylabel(rf"$\Delta$SIF at {lam_b:.1f} nm ({UNITS})")
    axs[1, 0].set_ylabel(rf"TOA SIF / $\sigma$ at {lam_b:.1f} nm")
    aod = ds["aod"].isel(band=ib).values
    aod_txt = ", ".join(f"×{a:g}: {d:.3f}" for a, d in zip(aer, aod)) if np.isfinite(aod).any() else "none"
    fig.suptitle(f"{scene_text(ds)}\nAOD at {lam_b:.1f} nm ({aod_txt}); "
                 f"column Rayleigh OD at ×1: {float(ds['rayleigh_od_column'].values[ib]):.3f}", fontsize=11)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"saved {out}")

    print(f"Relative ΔSIF at {lam_b:.1f} nm (%), same for every SIF intensity:")
    print(f"{'VZA':>5} {'Rayleigh':>9} " + " ".join(f"{'aer×' + format(a, 'g'):>9}" for a in aer))
    for iv in range(n_v):
        for ir, r in enumerate(ray):
            print(f"{ds['vza'].values[iv]:5g} {r:9g} " + " ".join(f"{x:9.2f}" for x in rel[iv, ir]))


def scene_text(ds):
    return (f"{Path(ds.attrs['ocean_yaml']).name}: SZA {ds.attrs['sza']:g}°, "
            f"albedo {ds.attrs['surface_albedo']:g}, GCHPIO column {ds.attrs['aer_column']}, "
            f"sunlight {'on' if ds.attrs['sunlight_on'] else 'off'}")


def main(path, vza, rayleigh, aerosol, band_nm):
    path = Path(path)
    ds = xr.open_dataset(path)
    iv, ir, ia = (nearest(ds["vza"].values, vza), nearest(ds["rayleigh_scale"].values, rayleigh),
                  nearest(ds["aerosol_scale"].values, aerosol))
    plot_spectra(ds, iv, ir, ia, path.with_name(f"{path.stem}_by_sif_spectra.png"))
    plot_band(ds, band_nm, path.with_name(f"{path.stem}_by_sif_{band_nm:g}nm.png"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("nc", nargs="?", default=DEFAULT_NC, help="grid NetCDF (default: DEFAULT_NC)")
    ap.add_argument("--vza", type=float, default=30.0, help="view for figure 1 (nearest VZA, default 30)")
    ap.add_argument("--rayleigh", type=float, default=1.0, help="Rayleigh scale for figure 1 (default 1)")
    ap.add_argument("--aerosol", type=float, default=1.0, help="aerosol scale for figure 1 (default 1)")
    ap.add_argument("--band", type=float, default=678.0, help="OCI band for figure 2, nm (default 678)")
    args = ap.parse_args()
    main(args.nc, args.vza, args.rayleigh, args.aerosol, args.band)
