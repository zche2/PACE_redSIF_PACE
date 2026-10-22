"""SVD of transmittance spectra (output of plot_transmittance_scattering.py), with or without log transform.

    python svd_transmittance.py [transmittance_scattering_n8756.nc] [--window 620 780] [--n-show 6] [--n-max 30] [--no-log]

Same convention as load_svd_basis (global_svd_fit_pipeline/svd_retrieval/svd_helpers.jl): bands with
λ_min < λ < λ_max, A = log(max(T, 1e-10)) (default) or A = T − 1 (--no-log), no mean removal, SVD of A'
(band × sample), PCs = left singular vectors. Done for `transmittance` (gas × scattering) and, for comparison,
`transmittance_gas`.

Outputs next to the input: <stem>_svd_<lo>_<hi>[_lin].nc (PCs, singular values) and the matching .png:
(a) singular values (% of sum, as in the pipeline), then one panel per leading PC (gas × scattering solid,
gas only dashed, sign-aligned).
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

HERE = Path(__file__).absolute().parent
DEFAULT_NC = HERE / "transmittance_scattering_n8756.nc"


def trans_svd(T, log_transform=True):
    """T: (sample, band). SVD of log(T) or T − 1. Returns U (band, k), S (k)."""
    A = np.log(np.maximum(T, 1e-10)) if log_transform else T - 1.0
    U, S, _ = np.linalg.svd(A.T, full_matrices=False)
    return U, S


def main(path, window, n_show, n_max, log_transform=True):
    qty = "log(T)" if log_transform else "T − 1"
    path = Path(path)
    ds = xr.open_dataset(path)
    lam_all = ds.band.values
    sel = (lam_all > window[0]) & (lam_all < window[1])
    lam = lam_all[sel]
    T = {"gas × scattering": ds.transmittance.transpose("sample", "band").values[:, sel],
         "gas only": ds.transmittance_gas.transpose("sample", "band").values[:, sel]}
    print(f"{T['gas only'].shape[0]} spectra, {lam.size} bands in ({window[0]:g}, {window[1]:g}) nm "
          f"[{lam.min():.1f}–{lam.max():.1f}]")

    res = {k: trans_svd(v, log_transform) for k, v in T.items()}
    U_s, S_s = res["gas × scattering"]
    U_g, S_g = res["gas only"]
    cosines = np.linalg.svd(U_s[:, :n_show].T @ U_g[:, :n_show], compute_uv=False)
    print(f"principal angles between the first {n_show} PCs (deg):",
          np.round(np.degrees(np.arccos(np.clip(cosines, -1, 1))), 2))
    for name, (U, S) in res.items():
        var = S ** 2 / np.sum(S ** 2)
        print(f"{name:17s} S% (pipeline norm) {np.round(100 * S[:8] / S.sum(), 3)}  "
              f"cum. variance of first 1/3/6/10 PCs: " + ", ".join(f"{100 * var[:k].sum():.4f}%" for k in (1, 3, 6, 10)))

    stem = f"{path.stem}_svd_{window[0]:g}_{window[1]:g}" + ("" if log_transform else "_lin")
    out = xr.Dataset(
        {"PCs": (("band", "pc"), U_s, {"long_name": f"left singular vectors of {qty}, T = gas × scattering"}),
         "S": ("pc", S_s, {"long_name": "singular values, gas × scattering"}),
         "PCs_gas": (("band", "pc"), U_g, {"long_name": f"left singular vectors of {qty}, T = gas only"}),
         "S_gas": ("pc", S_g, {"long_name": "singular values, gas only"})},
        coords={"band": ("band", lam, {"units": "nm"}), "pc": np.arange(1, U_s.shape[1] + 1)},
        attrs={"source": str(path), "window_nm": list(window), "log_transform": int(log_transform), "mean_removed": 0,
               **{f"src_{k}": v for k, v in ds.attrs.items() if k != "source"}})
    out.to_netcdf(path.parent / f"{stem}.nc")
    print(f"saved {path.parent / f'{stem}.nc'}")

    n_max = min(n_max, lam.size)
    n_col = 3
    n_row = int(np.ceil(n_show / n_col))
    fig = plt.figure(figsize=(5 * (n_col + 1), 3.6 * n_row), layout="constrained")
    gs = fig.add_gridspec(n_row, n_col + 1)
    ax_s = fig.add_subplot(gs[:, 0])
    for (name, (U, S)), ls in zip(res.items(), ("-", "--")):
        ax_s.semilogy(np.arange(1, n_max + 1), (100 * S / S.sum())[:n_max], "o" + ls, ms=3, label=name)
    ax_s.set_xlabel("PC")
    ax_s.set_ylabel("Singular value (% of sum)")
    ax_s.set_title(f"(a) Singular values of {qty}", loc="left")
    ax_s.legend(fontsize=8, frameon=False)
    ax_s.grid(alpha=0.3)
    ax_s.spines[["top", "right"]].set_visible(False)

    for i in range(n_show):
        ax = fig.add_subplot(gs[i // n_col, 1 + i % n_col])
        sgn = np.sign(U_s[:, i] @ U_g[:, i]) or 1.0
        cos = abs(U_s[:, i] @ U_g[:, i])
        ax.plot(lam, U_s[:, i], color="C0", lw=1.3, label=f"gas × scattering ({100 * S_s[i] / S_s.sum():.3g}%)")
        ax.plot(lam, sgn * U_g[:, i], color="C1", lw=1.1, ls="--", label=f"gas only ({100 * S_g[i] / S_g.sum():.3g}%)")
        ax.axhline(0, color="0.8", lw=0.6)
        ax.set_title(f"({chr(98 + i)}) PC{i + 1}   |cos| = {cos:.3f}", loc="left", fontsize=10)
        ax.legend(fontsize=7.5, frameon=False, loc="lower left")
        ax.spines[["top", "right"]].set_visible(False)
        if i // n_col == n_row - 1:
            ax.set_xlabel("Wavelength (nm)")
    fig.suptitle(f"{path.name}: SVD of {qty}, {lam.size} bands in ({window[0]:g}, {window[1]:g}) nm, "
                 f"{T['gas only'].shape[0]} spectra, no mean removal", fontsize=11)
    fig.savefig(path.parent / f"{stem}.png", dpi=180)
    plt.close(fig)
    print(f"saved {path.parent / f'{stem}.png'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("nc", nargs="?", default=str(DEFAULT_NC))
    ap.add_argument("--window", type=float, nargs=2, default=(620.0, 780.0), help="λ_min λ_max, nm (exclusive)")
    ap.add_argument("--n-show", type=int, default=6, help="PCs plotted, one panel each")
    ap.add_argument("--n-max", type=int, default=30, help="PCs in panel (a)")
    ap.add_argument("--no-log", action="store_true", help="SVD of T − 1 instead of log(T) (svd_log_transform = false)")
    args = ap.parse_args()
    main(args.nc, args.window, args.n_show, args.n_max, not args.no_log)
