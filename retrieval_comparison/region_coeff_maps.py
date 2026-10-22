"""Region maps of every fitted state element for three O2A SVD retrievals of one
granule, styled after notebooks/figures/Figure2/Figure2_turbulent_feature.ipynb
(Gulf of Maine / Georges Bank box, L1B RGB basemap, nFLH from L2 AOP).

    python retrieval_comparison/region_coeff_maps.py

Retrievals
  npoly3     new_svd_retrieval_output_O2A          Aug01 transmittance SVD, 3 Legendre
  npoly5     new_svd_retrieval_output_npoly5       Aug01 transmittance SVD, 5 Legendre
  newSVD_npoly3  new_svd_retrieval_output_O2A_newSVD         scattering-included SVD, 3 Legendre
  newSVD_npoly2  new_svd_retrieval_output_O2A_newSVD_npoly2  scattering-included SVD, 2 Legendre

npoly3 and npoly5 share one SVD basis, so their PC coefficients are directly
comparable; newSVD uses a different basis, so its PC maps get their own colour
scale and are never differenced against the others.
"""
from datetime import datetime
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.ticker as mticker
import netCDF4 as nc
import numpy as np
import xarray as xr
from cartopy.mpl.ticker import LatitudeFormatter
from matplotlib.colors import LogNorm
from scipy.spatial import cKDTree

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

GRANULE = "20250724T170203"
NAME, CENTER = "Gulf of Maine / Georges Bank", (41.0, -68.5)
LAT_HALF, LON_HALF = 2.5, 4.0
TRANSECT = ((43.2, -67.4), (38.9, -67.4))

PACE = Path("/home/zhe2/data/PACE")
KIWI = Path("/kiwi-data/Data/satellite/PACE_OCI")
RUNS_ALL = [  # (key, label, folder, SVD basis tag); order = column order
    ("npoly2", "npoly2\n(Aug01 SVD)", "new_svd_retrieval_output_npoly2", "aug01"),
    ("npoly3", "npoly3\n(Aug01 SVD)", "new_svd_retrieval_output_O2A", "aug01"),
    ("npoly5", "npoly5\n(Aug01 SVD)", "new_svd_retrieval_output_npoly5", "aug01"),
    ("newSVD_npoly2", "newSVD, npoly2\n(scattering SVD)", "new_svd_retrieval_output_O2A_newSVD_npoly2", "new"),
    ("newSVD_npoly3", "newSVD, npoly3\n(scattering SVD)", "new_svd_retrieval_output_O2A_newSVD", "new"),
]
OUT = Path(__file__).resolve().parent / "out" / "region"  # all runs that are complete
OUT.mkdir(parents=True, exist_ok=True)
UNITS = r"W m$^{-2}$ sr$^{-1}$ µm$^{-1}$"
RGB_FLOOR_PCT, RGB_RHO_MAX, RGB_GAMMA = 1, 0.35, 0.8
PROJ = ccrs.PlateCarree()
FL_CLIM, SNR_CLIM = (-0.1, 0.5), (1.0, 50.0)
lines = []


def say(s=""):
    print(s)
    lines.append(s)


# ----------------------------------------------------------------- loading
def dated(root, name):
    t = datetime.strptime(GRANULE, "%Y%m%dT%H%M%S")
    return root / f"{t:%Y}" / f"{t:%m}" / f"{t:%d}" / name


def _ready(folder):
    f = dated(PACE / folder, f"interim_{GRANULE}_svd_retrieval.nc")
    try:
        with nc.Dataset(f) as d:
            return f.stat().st_size > 3e8 and d["x_hat"].shape[0] > 0
    except Exception:
        return False


RUNS = [r for r in RUNS_ALL if _ready(r[2])]
for r in RUNS_ALL:
    if r not in RUNS:
        print(f"[skip] {r[0]}: output not complete / not readable yet")


def dark_basemap(rgb):
    valid = np.all(np.isfinite(rgb) & (rgb > 0), axis=-1)
    floor = np.array([np.nanpercentile(rgb[..., i][valid], RGB_FLOOR_PCT) for i in range(3)])
    out = np.clip((rgb - floor) / (RGB_RHO_MAX - floor.mean()), 0.0, 1.0) ** RGB_GAMMA
    out[~valid] = 0.0
    return out


aop = dated(KIWI / "L2_AOP_V3.2", f"PACE_OCI.{GRANULE}.L2.OC_AOP.V3_2.nc")
nav = xr.open_dataset(aop, group="navigation_data")
LON = nav["longitude"].values.astype(float)
LAT = nav["latitude"].values.astype(float)
NFLH = xr.open_dataset(aop, group="geophysical_data")["nflh"].values.astype(float)

lat0, lon0 = CENTER
inbox = (np.abs(LAT - lat0) < LAT_HALF + 0.3) & (np.abs(LON - lon0) < LON_HALF + 0.5)
rr, cc = np.where(inbox)
R0, R1, C0, C1 = rr.min(), rr.max() + 1, cc.min(), cc.max() + 1
SL = (slice(R0, R1), slice(C0, C1))
say(f"region block rows {R0}:{R1}  cols {C0}:{C1}  ({R1 - R0} x {C1 - C0} pixels)")

cache = OUT / f"basemap_{GRANULE}.npy"
if cache.exists():
    BASE = np.load(cache)
else:
    l1b = dated(KIWI / "L1B_V3", f"PACE_OCI.{GRANULE}.L1B.V3.nc")
    sb = xr.open_dataset(l1b, group="sensor_band_parameters")
    obs = xr.open_dataset(l1b, group="observation_data")
    near = lambda wl, t: int(np.argmin(np.abs(wl - t)))
    bw, rw = sb["blue_wavelength"].values, sb["red_wavelength"].values
    rgb = np.stack([
        obs["rhot_red"].isel(red_bands=near(rw, 645.0)).values,
        obs["rhot_blue"].isel(blue_bands=near(bw, 555.0)).values,
        obs["rhot_blue"].isel(blue_bands=near(bw, 469.0)).values,
    ], axis=-1).astype(float)
    BASE = dark_basemap(rgb)[SL]
    np.save(cache, BASE)
BASE = BASE if BASE.shape[:2] == (R1 - R0, C1 - C0) else dark_basemap(BASE)[SL]


def load_run(folder):
    ds = xr.open_dataset(dated(PACE / folder, f"interim_{GRANULE}_svd_retrieval.nc"))
    names = ds.attrs["state_names_csv"].split(",")
    x = ds["x_hat"].values.astype(float)
    var = ds["S_posterior_diag"].values.astype(float)
    sif = ds["sif_radiance_678nm"].values.astype(float)
    ev1 = x[-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = sif / ev1
    w = float(np.nanmedian(ratio[np.isfinite(ratio) & (np.abs(ev1) > 1e-6)]))
    r = dict(names=names, x=x, sif=sif, unc=np.abs(w) * np.sqrt(np.maximum(var[-1], 0.0)),
             converged=ds["converged"].values == 1, chi2=ds["reduced_chi2"].values.astype(float),
             ocean=ds["is_ocean"].values == 1, dark=ds["is_dark"].values == 1,
             status=ds["status_code"].values)
    r["alpha"] = 10.0 / (1.0 + np.exp(-x[names.index("alpha")])) + 1.0
    return r


D = {k: load_run(f) for k, _, f, _ in RUNS}
HAS = np.isfinite(NFLH)
for k in D:
    D[k]["sif_ok"] = D[k]["ocean"] & D[k]["dark"] & D[k]["converged"] & HAS & np.isfinite(D[k]["sif"])
    D[k]["snr"] = np.where(HAS, D[k]["sif"] / D[k]["unc"], np.nan)
inreg = (np.abs(LAT - lat0) < LAT_HALF) & (np.abs(LON - lon0) < LON_HALF) & HAS


# ----------------------------------------------------------------- plotting
def cmap_clear(name):
    c = plt.get_cmap(name).copy()
    c.set_bad(alpha=0.0)
    return c


def panel(ax, field, cmap, show_ylabels=False, show_xlabels=False, **kw):
    lon, lat = LON[SL], LAT[SL]
    ax.pcolormesh(lon, lat, BASE, transform=PROJ, shading="auto", rasterized=True)
    pm = ax.pcolormesh(lon, lat, np.ma.masked_invalid(field[SL]), transform=PROJ, cmap=cmap_clear(cmap),
                       shading="auto", rasterized=True, **kw)
    ax.set_extent([lon0 - LON_HALF, lon0 + LON_HALF, lat0 - LAT_HALF, lat0 + LAT_HALF], crs=PROJ)
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=0.5, edgecolor="k")
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color="0.4", alpha=0.5, linestyle=":")
    gl.top_labels = gl.right_labels = False
    gl.left_labels, gl.bottom_labels = show_ylabels, show_xlabels
    gl.xlocator, gl.ylocator = mticker.MultipleLocator(2), mticker.MultipleLocator(1)
    gl.xlabel_style = gl.ylabel_style = {"size": 8}
    (la_a, lo_a), (la_b, lo_b) = TRANSECT
    ax.plot([lo_a, lo_b], [la_a, la_b], color="w", lw=1.6, transform=PROJ, zorder=5)
    ax.plot([lo_a, lo_b], [la_a, la_b], color="k", lw=0.7, ls="--", transform=PROJ, zorder=6)
    return pm


def robust_lim(arrs, sym=True, q=(2, 98)):
    v = np.concatenate([a[SL][inreg[SL]] for a in arrs])
    v = v[np.isfinite(v)]
    lo, hi = np.percentile(v, q)
    if sym:
        m = max(abs(lo), abs(hi))
        return -m, m
    return lo, hi


# --------------------------------------------------------------- figure 1
def tag(ax, text, size=8.5):
    ax.text(0.015, 0.975, text, transform=ax.transAxes, ha="left", va="top", fontsize=size, fontweight="bold",
            zorder=8, bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.85))


def fig_overview():
    cols = [("L2 nFLH", "nflh", "PuBuGn", dict(vmin=FL_CLIM[0], vmax=FL_CLIM[1])),
            ("SIF 678 nm", "sif", "PuBuGn", dict(vmin=FL_CLIM[0], vmax=FL_CLIM[1])),
            (r"SNR (SIF / $\sigma_{SIF}$)", "snr", "magma", dict(norm=LogNorm(*SNR_CLIM))),
            ("reduced $\\chi^2$", "chi2", "cividis", dict(norm=LogNorm(0.5, 10))),
            ("converged (green) / not (red)", "conv", "RdYlGn", dict(vmin=0, vmax=1))]
    nr = len(RUNS)
    fig = plt.figure(figsize=(17, 3.1 * nr + 1.1))
    gs = fig.add_gridspec(nr + 1, 5, height_ratios=[1] * nr + [0.07], hspace=0.1, wspace=0.07,
                          left=0.05, right=0.99, top=0.94, bottom=0.07)
    pm = {}
    for i, (key, label, _, _) in enumerate(RUNS):
        d = D[key]
        for j, (title, f, cm, kw) in enumerate(cols):
            ax = fig.add_subplot(gs[i, j], projection=PROJ)
            fld = {"nflh": NFLH, "sif": np.where(HAS, d["sif"], np.nan), "snr": d["snr"],
                   "chi2": np.where(HAS, d["chi2"], np.nan),
                   "conv": np.where(HAS, d["converged"].astype(float), np.nan)}[f]
            pm[f] = panel(ax, fld, cm, show_ylabels=j == 0, show_xlabels=i == len(RUNS) - 1, **kw)
            if i == 0:
                ax.set_title(title, fontsize=11)
            if j == 0:
                tag(ax, label.replace("\n", " "))
    for cols_, f, lab in (((0, 2), "sif", f"nFLH, SIF ({UNITS})"), ((2, 3), "snr", r"SNR (SIF / $\sigma_{SIF}$)"),
                          ((3, 4), "chi2", "reduced $\\chi^2$")):
        cax = fig.add_subplot(gs[len(RUNS), cols_[0]:cols_[1]])
        cb = fig.colorbar(pm[f], cax=cax, orientation="horizontal", extend="both")
        cb.ax.tick_params(labelsize=8)
        cb.set_label(lab, fontsize=9)
    fig.suptitle(f"{NAME}  {GRANULE[:4]}-{GRANULE[4:6]}-{GRANULE[6:8]}   (RGB where no L2 nFLH; dashed line = transect)",
                 fontsize=11)
    fig.savefig(OUT / "A_region_overview.png", dpi=130)
    plt.close(fig)


# ------------------------------------------------------------- coefficient pages
def element_field(d, name):
    if name == "alpha_phys":
        return np.where(HAS, d["alpha"], np.nan)
    if name == "sif_678":
        return np.where(HAS, d["sif"], np.nan)
    if name not in d["names"]:
        return None
    return np.where(HAS, d["x"][d["names"].index(name)], np.nan)


def _fmt(cb, lim):
    if max(abs(lim[0]), abs(lim[1])) < 2e-3:
        cb.formatter = mticker.FormatStrFormatter("%.1e")
        cb.update_ticks()


def fig_page(elements, fname, title):
    """elements: list of (state name, label, cmap, sym, share_all[, fixed (vmin, vmax)])."""
    n, nm = len(elements), len(RUNS)
    aug = [j for j, r in enumerate(RUNS) if r[3] == "aug01"]
    oth = [j for j, r in enumerate(RUNS) if r[3] != "aug01"]
    share_any = any(e[4] for e in elements)
    if share_any:   # one colour bar per row, shared by all columns
        wr, mapcol = [1] * nm + [0.045], list(range(nm))
    else:           # bar 1 for the Aug01-basis runs, bar 2 for the other basis
        wr = [1] * len(aug) + [0.085] + [1] * len(oth) + [0.085]
        mapcol = [None] * nm
        for k, j in enumerate(aug):
            mapcol[j] = k
        for k, j in enumerate(oth):
            mapcol[j] = len(aug) + 1 + k
    fig = plt.figure(figsize=(3.9 * nm + 1.2, 3.05 * n + 0.8))
    gs = fig.add_gridspec(n, len(wr), width_ratios=wr, hspace=0.1, wspace=0.12, left=0.055,
                          right=0.965 if share_any else 0.94, top=1 - 0.7 / (3.05 * n + 0.8), bottom=0.045)
    for i, el in enumerate(elements):
        name, label, cm, sym, share_all = el[:5]
        fixed = el[5] if len(el) > 5 else None
        flds = [element_field(D[k], name) for k, *_ in RUNS]
        have = [f for f in flds if f is not None]
        lim = [None] * nm
        if fixed is not None:
            lim = [fixed] * nm
        elif share_all:
            lim = [robust_lim(have, sym)] * nm
        else:
            ga = robust_lim([flds[j] for j in aug if flds[j] is not None], sym)
            for j in aug:
                lim[j] = ga
            go = robust_lim([flds[j] for j in oth if flds[j] is not None], sym) if oth else None
            for j in oth:
                lim[j] = go
        pms = [None] * nm
        for j, (k, lab, _, basis) in enumerate(RUNS):
            ax = fig.add_subplot(gs[i, mapcol[j]], projection=PROJ)
            if flds[j] is None:
                ax.set_axis_off()
                ax.text(0.5, 0.5, f"not in the {k} state vector", transform=ax.transAxes,
                        ha="center", va="center", fontsize=9, color="0.4")
            else:
                pms[j] = panel(ax, flds[j], cm, show_ylabels=j == 0, show_xlabels=i == n - 1,
                               vmin=lim[j][0], vmax=lim[j][1])
            if i == 0:
                ax.set_title(lab.replace("\n", " "), fontsize=9.5)
            if j == 0:
                tag(ax, label.replace("\n", " "), 10)
        if share_any:
            cb = fig.colorbar(next(p for p in pms if p is not None), cax=fig.add_subplot(gs[i, nm]), extend="both")
            cb.ax.tick_params(labelsize=7)
            _fmt(cb, lim[aug[0]] if lim[aug[0]] else fixed)
        else:
            cb = fig.colorbar(next(pms[j] for j in aug if pms[j] is not None),
                              cax=fig.add_subplot(gs[i, len(aug)]), extend="both")
            cb.ax.tick_params(labelsize=7)
            _fmt(cb, lim[aug[0]])
            first_oth = next((j for j in oth if pms[j] is not None), None)
            if first_oth is not None:
                cb2 = fig.colorbar(pms[first_oth], cax=fig.add_subplot(gs[i, len(wr) - 1]), extend="both")
                cb2.ax.tick_params(labelsize=7)
                _fmt(cb2, lim[first_oth])
    fig.suptitle(title, fontsize=10)
    fig.savefig(OUT / fname, dpi=120)
    plt.close(fig)


def fig_pages():
    pcs = [(f"pc_{k}", f"PC {k}", "RdBu_r", True, False) for k in range(1, 16)]
    basis_note = ("PC coefficients. 1st colour bar: npoly3 & npoly5 (same Aug01 SVD basis).\n"
                  "2nd colour bar: both newSVD runs (scattering SVD, a different basis) -> own scale, values not comparable with the 1st")
    for p in range(3):
        fig_page(pcs[5 * p:5 * p + 5], f"B{p + 1}_pc_coeffs_{5 * p + 1}-{5 * p + 5}.png",
                 f"{NAME}: {basis_note}")
    rest = [("alpha_phys", r"$\alpha$ (path-length exp.)", "viridis", False, True, (7.0, 11.0)),
            ("legendre_p0", "Legendre p0", "RdBu_r", True, True),
            ("legendre_p1", "Legendre p1", "RdBu_r", True, True),
            ("legendre_p2", "Legendre p2", "RdBu_r", True, True),
            ("legendre_p3", "Legendre p3", "RdBu_r", True, True)]
    fig_page(rest, "C_alpha_legendre_p0-p3.png",
             f"{NAME}: alpha = 10/(1+e^-raw)+1 and baseline Legendre coefficients (shared colour scale)")
    last = [("legendre_p4", "Legendre p4", "RdBu_r", True, True),
            ("legendre_p5", "Legendre p5", "RdBu_r", True, True),
            ("sif_ev1", "SIF coeff\n(sif_ev1)", "RdBu_r", True, True),
            ("sif_678", "SIF 678 nm", "PuBuGn", False, True)]
    fig_page(last, "D_legendre_p4-p5_sif.png", f"{NAME}: higher Legendre terms (npoly5 only) and SIF")


# ---------------------------------------------------------------- transect
R_EARTH = 6371.0


def xyz(lon, lat):
    lo, la = np.deg2rad(lon), np.deg2rad(lat)
    return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def transect_vals(fields, step_km=1.0, max_km=1.5):
    (la_a, lo_a), (la_b, lo_b) = TRANSECT
    pa, pb = xyz(np.array([lo_a]), np.array([la_a]))[0], xyz(np.array([lo_b]), np.array([la_b]))[0]
    om = np.arccos(np.clip(pa @ pb, -1, 1))
    t = np.linspace(0, 1, int(om * R_EARTH // step_km) + 1)
    pts = (np.sin((1 - t) * om)[:, None] * pa + np.sin(t * om)[:, None] * pb) / np.sin(om)
    ok = np.isfinite(LON) & np.isfinite(LAT) & HAS
    chord, idx = cKDTree(xyz(LON[ok], LAT[ok])).query(pts, k=1)
    near = 2 * R_EARTH * np.arcsin(np.clip(chord / 2, 0, 1)) <= max_km
    lat_p = np.rad2deg(np.arcsin(pts[:, 2]))
    return lat_p, {k: np.where(near, v[ok][idx], np.nan) for k, v in fields.items()}


def fig_transect():
    C = {"npoly2": "#e7298a", "npoly3": "#0B7A75", "npoly5": "#d95f02", "newSVD_npoly2": "#7570b3", "newSVD_npoly3": "#a6761d"}
    fields = {"nflh": NFLH}
    for k in D:
        d = D[k]
        fields[f"sif_{k}"] = np.where(d["sif_ok"], d["sif"], np.nan)
        fields[f"alpha_{k}"] = np.where(d["sif_ok"], d["alpha"], np.nan)
        fields[f"p0_{k}"] = np.where(d["sif_ok"], d["x"][d["names"].index("legendre_p0")], np.nan)
    lat_p, v = transect_vals(fields)
    fig, axs = plt.subplots(3, 1, figsize=(8.5, 8), sharex=True)
    axs[0].plot(lat_p, v["nflh"], color="0.2", lw=1.0, label="L2 nFLH")
    for k, *_ in RUNS:
        axs[0].plot(lat_p, v[f"sif_{k}"], color=C[k], lw=0.9, label=f"SIF {k}")
        axs[1].plot(lat_p, v[f"alpha_{k}"], color=C[k], lw=0.9, label=k)
        axs[2].plot(lat_p, v[f"p0_{k}"], color=C[k], lw=0.9, label=k)
    axs[0].set_ylabel(f"Radiance ({UNITS})", fontsize=9)
    axs[0].legend(ncol=3, frameon=False, fontsize=8, loc="upper right")
    axs[1].set_ylabel(r"$\alpha$", fontsize=10)
    axs[2].set_ylabel("Legendre p0", fontsize=10)
    axs[2].set_xlim(lat_p[0], lat_p[-1])
    axs[2].xaxis.set_major_formatter(LatitudeFormatter())
    axs[2].set_xlabel("Latitude (A at 43.2°N → B at 38.9°N, lon −67.4°)")
    for a in axs:
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(labelsize=8)
    fig.suptitle(f"{NAME}: transect (converged, dark, ocean pixels only)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "E_transect.png", dpi=130)
    plt.close(fig)
    for k in D:
        s = v[f"sif_{k}"]
        g = np.isfinite(s) & np.isfinite(v["nflh"])
        say(f"transect {k:7s}: {g.sum():4d}/{g.size} pts   r(nFLH,SIF)={np.corrcoef(v['nflh'][g], s[g])[0, 1]:.3f}"
            f"   SIF/nFLH slope={np.polyfit(v['nflh'][g], s[g], 1)[0]:.3f}")


# ---------------------------------------------------------------- statistics
def stats():
    keys = [r[0] for r in RUNS]
    allconv = inreg.copy()
    for k in keys:
        allconv &= D[k]["sif_ok"]
    say(f"\nruns compared: {keys}")
    say(f"region pixels with nFLH: {inreg.sum()};  valid & converged & dark in ALL runs: {allconv.sum()}")
    for k in keys:
        say(f"  {k:7s} converged (region, nFLH-valid): {100 * D[k]['converged'][inreg].mean():5.1f}%   "
            f"sif_ok: {100 * D[k]['sif_ok'][inreg].mean():5.1f}%   "
            f"median chi2 {np.nanmedian(D[k]['chi2'][inreg]):.2f}   median sigma {np.nanmedian(D[k]['unc'][inreg]):.4f}")
    for k in keys:
        a = D[k]["alpha"][allconv]
        say(f"  alpha {k:7s}: median {np.median(a):.2f}   {100 * np.mean(a > 10.5):.1f}% of pixels within 0.5 of the upper bound (11)")
    gref = {"aug01": "npoly3", "new": "newSVD_npoly3"}
    grp = {r[0]: r[3] for r in RUNS}
    say(f"\nper-element mean/std over the {allconv.sum()} common pixels, and pixel-wise r.  "
        "PCs: r vs the reference run of the same SVD basis (npoly3 for Aug01, newSVD_npoly3 for scattering SVD); "
        "other elements: r vs npoly3.  '(ref)' = reference itself, 'n/a' = element not in that state vector")
    say(f"{'element':14s}" + "".join(f"{'mean/std ' + k:>24s}" for k in keys) + "".join(f"{'r ' + k:>16s}" for k in keys))
    elems = ([f"pc_{k}" for k in range(1, 16)] + ["alpha_phys"] + [f"legendre_p{k}" for k in range(6)] + ["sif_ev1"])
    for e in elems:
        f = {k: element_field(D[k], e) for k in keys}
        if all(v is None for v in f.values()):
            continue
        v = {k: (f[k][allconv] if f[k] is not None else None) for k in keys}
        ms = "".join(f"{(f'{np.mean(v[k]):10.3e}/{np.std(v[k]):8.2e}' if v[k] is not None else 'n/a'):>24s}" for k in keys)
        rs = ""
        for k in keys:
            ref = gref[grp[k]] if e.startswith("pc_") else "npoly3"
            if ref not in v or v[ref] is None or v[k] is None:
                rs += f"{'n/a':>16s}"
            elif k == ref:
                rs += f"{'(ref)':>16s}"
            else:
                rs += f"{np.corrcoef(v[ref], v[k])[0, 1]:16.3f}"
        say(f"{e:14s}{ms}{rs}")


if __name__ == "__main__":
    fig_overview()
    fig_pages()
    fig_transect()
    stats()
    (OUT / "summary.txt").write_text("\n".join(lines) + "\n")
    print("done ->", OUT)
