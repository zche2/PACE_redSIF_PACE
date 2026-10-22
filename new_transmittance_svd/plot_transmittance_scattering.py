"""Convolved gas transmittance with an added Rayleigh/Mie-like scattering term.

    python plot_transmittance_scattering.py [--n N] [--n-plot 8] [--lam0 680] [--tau0 0.01 0.3] [--alpha 0 4] [--relative]

Generates n spectra (one per distinct profile, n <= number of profiles in NC; default: all) and writes them to
NetCDF, then plots n_plot of them:
(a) one-way O2 + H2O transmittance T_gas(λ)
(b) the same spectra × exp(-m τ0 (λ/λ0)^-α), with m the profile's one-way AMF, and τ0 (log-uniform)
    and α (uniform) drawn per spectrum. τ(λ) = τ0 (λ/λ0)^-α is the Ångström law (Rayleigh α ≈ 4,
    aerosol α ≈ 0–2.5). With --relative the factor is exp(-m τ0 [(λ/λ0)^-α - 1]): 1 at λ0, so only the
    spectral shape changes (it exceeds 1 beyond λ0).
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

HERE = Path(__file__).absolute().parent
NC = Path("/home/zhe2/data/MyProjects/PACE_redSIF_PACE/convolved_transmittance/"
          "transmittance_winter_FineWvResModel_FullRange_Aug01.nc")
MAX_LEGEND = 12


def scattering_factor(lam, lam0, tau0, alpha, amf, relative=False):
    """exp(-m τ0 [(λ/λ0)^-α - c]) for arrays of samples (rows) × bands (columns)."""
    shape = (lam[None, :] / lam0) ** (-alpha[:, None]) - (1.0 if relative else 0.0)
    return np.exp(-(amf * tau0)[:, None] * shape)


def generate(ds, n, lam0, tau0_range, alpha_range, rng, relative):
    lam = ds.band.values
    T = ds.transmittance.transpose("profile", "band").values
    amf = ds.AMF.values
    idx = np.sort(rng.choice(T.shape[0], n, replace=False))
    tau0 = np.exp(rng.uniform(*np.log(tau0_range), n))
    alpha = rng.uniform(*alpha_range, n)
    S = scattering_factor(lam, lam0, tau0, alpha, amf[idx], relative)
    # (band, sample) like the source file, so Julia (load_svd_basis) reads (sample, band)
    return xr.Dataset(
        {
            "transmittance_gas": (("band", "sample"), T[idx].T, {"long_name": "one-way O2 + H2O transmittance"}),
            "scattering_factor": (("band", "sample"), S.T,
                                  {"long_name": "exp(-AMF tau0 [(lambda/lambda0)^-alpha" + (" - 1" if relative else "") + "])"}),
            "transmittance": (("band", "sample"), (T[idx] * S).T,
                              {"long_name": "gas transmittance × scattering factor"}),
            "profile_index": ("sample", idx, {"long_name": "0-based profile index in the source file"}),
            "AMF": ("sample", amf[idx], {"long_name": "one-way air mass factor"}),
            "tau0": ("sample", tau0, {"long_name": "scattering optical depth at lambda0"}),
            "alpha": ("sample", alpha, {"long_name": "Angstrom exponent"}),
        },
        coords={"band": ("band", lam, {"long_name": "Wavelength", "units": "nm"})},
        attrs={"source": str(NC), "lambda0_nm": lam0, "tau0_range": list(tau0_range), "tau0_sampling": "log-uniform",
               "alpha_range": list(alpha_range), "alpha_sampling": "uniform", "relative": int(relative)},
    )


def plot(g, k_plot, lam0, tau0_range, alpha_range, relative, out):
    lam = g.band.values
    cols = plt.get_cmap("viridis")(np.linspace(0.05, 0.9, len(k_plot)))
    legend = len(k_plot) <= MAX_LEGEND
    fig, ax = plt.subplots(1, 2, figsize=(14, 5), layout="constrained", sharey=True)
    for c, k in zip(cols, k_plot):
        s = g.isel(sample=k)
        S = s.scattering_factor.values
        ax[0].plot(lam, s.transmittance_gas, color=c, lw=1.1,
                   label=f"profile {int(s.profile_index)}, AMF {float(s.AMF):.2f}")
        ax[1].plot(lam, s.transmittance, color=c, lw=1.1,
                   label=f"τ$_0$ = {float(s.tau0):.3f}, α = {float(s.alpha):.2f}  (×{S[0]:.3f} → ×{S[-1]:.3f})")
    ax[0].set_title(f"(a) Convolved gas transmittance (O$_2$ + H$_2$O, one-way), {len(k_plot)} of {g.sizes['sample']}",
                    loc="left")
    expo = r"[(\lambda/\lambda_0)^{-\alpha}-1]" if relative else r"(\lambda/\lambda_0)^{-\alpha}"
    ax[1].set_title(rf"(b) $T_{{gas}}\,\exp(-m\,\tau_0{expo})$, " f"λ$_0$ = {lam0:g} nm, m = AMF", loc="left")
    ax[0].set_ylabel("Transmittance")
    if legend:
        ax[0].legend(fontsize=7.5, frameon=False, loc="lower left")
        ax[1].legend(fontsize=7.5, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.72))
    for a in ax:
        a.set_xlabel("Wavelength (nm)")
        a.axvline(lam0, color="0.8", lw=0.8)
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"{NC.name}: τ$_0$ log-uniform in [{tau0_range[0]:g}, {tau0_range[1]:g}], "
                 f"α uniform in [{alpha_range[0]:g}, {alpha_range[1]:g}]  (Rayleigh α≈4, aerosol α≈0–2.5)",
                 fontsize=10)
    fig.savefig(out, dpi=180)
    plt.close(fig)
    print(f"saved {out}")


def main(n, n_plot, lam0, tau0_range, alpha_range, seed, out_png, out_nc, relative=False):
    ds = xr.open_dataset(NC)
    n_prof = ds.sizes["profile"]
    n = n_prof if n is None else n
    if not 1 <= n <= n_prof:
        raise SystemExit(f"--n must be in 1..{n_prof} (profiles in {NC.name}), got {n}")
    n_plot = min(n_plot, n)

    rng = np.random.default_rng(seed)
    g = generate(ds, n, lam0, tau0_range, alpha_range, rng, relative)
    g.attrs["seed"] = seed
    out_nc = Path(out_nc or HERE / f"transmittance_scattering_n{n}{'_relative' if relative else ''}.nc")
    enc = {v: {"zlib": True, "complevel": 4} for v in g.data_vars}
    g.to_netcdf(out_nc, encoding=enc)
    print(f"saved {out_nc}  ({n} spectra × {g.sizes['band']} bands)")

    k_plot = np.sort(rng.choice(n, n_plot, replace=False))
    plot(g, k_plot, lam0, tau0_range, alpha_range, relative, out_png)
    if n_plot <= MAX_LEGEND:
        for k in k_plot:
            s = g.isel(sample=k)
            S = s.scattering_factor.values
            print(f"  profile {int(s.profile_index):5d}  AMF {float(s.AMF):.2f}  τ0 {float(s.tau0):.3f}  "
                  f"α {float(s.alpha):.2f}  factor {S[0]:.3f} at {g.band.values[0]:.0f} nm, "
                  f"{S[-1]:.3f} at {g.band.values[-1]:.0f} nm")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=None, help="spectra to generate, one per profile (default: all profiles)")
    ap.add_argument("--n-plot", type=int, default=8, help="spectra to plot (legend shown up to %d)" % MAX_LEGEND)
    ap.add_argument("--lam0", type=float, default=680.0, help="reference wavelength, nm (factor = 1 there)")
    ap.add_argument("--tau0", type=float, nargs=2, default=(0.01, 0.3), help="τ0 range at λ0 (log-uniform)")
    ap.add_argument("--alpha", type=float, nargs=2, default=(0.0, 4.0), help="Ångström exponent range")
    ap.add_argument("--seed", type=int, default=20261020)
    ap.add_argument("--relative", action="store_true", help="normalize the scattering factor to 1 at λ0")
    ap.add_argument("-o", "--output", default=str(HERE / "transmittance_lambda_alpha_scaling.png"), help="figure path")
    ap.add_argument("--nc-out", default=None, help="NetCDF path (default: transmittance_scattering_n<N>.nc here)")
    args = ap.parse_args()
    main(args.n, args.n_plot, args.lam0, args.tau0, args.alpha, args.seed, args.output, args.nc_out, args.relative)
