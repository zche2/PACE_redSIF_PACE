#!/usr/bin/env python3
"""Aggregate vicarious-gain QC over all granules for one day.

Instead of the per-granule notebook workflow, this streams every
``vcal_gain_interim_{date}T*_svd_retrieval.nc`` under a vcal_gain directory,
accumulates statistics, and writes the same style of plots:

  1. Per-pixel mean spectrum (cyclic colorbar on pixel index) + day mean
  2. Mean spectrum binned by nflh
  3. Mean spectrum binned by chlor_a
  4. Mean spectrum binned by SZA (from paired L1B)
  5. Mean spectrum binned by VZA (from paired L1B)

Example
-------
  python vcal_gain/analyze_vcal_gain_day.py --date 20250125
  python vcal_gain/analyze_vcal_gain_day.py --date 20250125 --n-bins 10 --l1b-dir /kiwi-data/.../L1B_V3
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

_GRANULE_RE = re.compile(r"vcal_gain_interim_(\d{8}T\d{6})_svd_retrieval\.nc$")


def _list_day_files(vcal_dir: Path, date: str) -> list[Path]:
    files = sorted(vcal_dir.glob(f"vcal_gain_interim_{date}T*_svd_retrieval.nc"))
    if not files:
        raise FileNotFoundError(
            f"No files matching vcal_gain_interim_{date}T*_svd_retrieval.nc in {vcal_dir}"
        )
    return files


def _resolve_l1b(
    ds: xr.Dataset,
    granule_id: str,
    l1b_dir: Path | None,
) -> Path | None:
    """Resolve L1B path: optional --l1b-dir, else NetCDF attr, else kiwi dated layout."""
    name = f"PACE_OCI.{granule_id}.L1B.V3.nc"
    y, m, d = granule_id[:4], granule_id[4:6], granule_id[6:8]

    if l1b_dir is not None:
        for cand in (
            l1b_dir / y / m / d / name,
            l1b_dir / name,
        ):
            if cand.is_file():
                return cand

    p = Path(str(ds.attrs.get("l1b_file", "")))
    if p.is_file():
        return p

    alt = Path("/kiwi-data/Data/satellite/PACE_OCI/L1B_V3") / y / m / d / name
    return alt if alt.is_file() else None


def _as_wl_scans_pixels(da: xr.DataArray) -> np.ndarray:
    """Return float32 array with dims (wavelength, scans, pixels)."""
    need = ("wavelength", "scans", "pixels")
    missing = [d for d in need if d not in da.dims]
    if missing:
        raise ValueError(f"vcal_gain missing dims {missing}; got {da.dims}")
    return np.asarray(da.transpose(*need).values, dtype=np.float32)


def _as_scans_pixels(da: xr.DataArray, *, n_pix: int | None = None) -> np.ndarray:
    """Return float32 array with dims (scans, pixels); optionally assert pixel count."""
    need = ("scans", "pixels")
    missing = [d for d in need if d not in da.dims]
    if missing:
        raise ValueError(f"{da.name} missing dims {missing}; got {da.dims}")
    out = np.asarray(da.transpose(*need).values, dtype=np.float32)
    if n_pix is not None and out.shape[1] != n_pix:
        raise ValueError(
            f"{da.name}: expected pixels={n_pix}, got shape {out.shape} (scans, pixels)"
        )
    return out


def _retrieval_window(ds: xr.Dataset, n_scan: int, n_pix: int) -> tuple[int, int, int, int]:
    """1-based inclusive L1B window aligned to retrieval / gain grid."""
    # Prefer attrs on the paired retrieval NC
    ret_path = Path(str(ds.attrs.get("retrieval_file", "")))
    src = None
    if ret_path.is_file():
        with xr.open_dataset(ret_path) as ret:
            src = ret
            if all(k in ret.attrs for k in ("pixel_start", "pixel_end", "scan_start", "scan_end")):
                ps = int(ret.attrs["pixel_start"])
                pe = int(ret.attrs["pixel_end"])
                ss = int(ret.attrs["scan_start"])
                se = int(ret.attrs["scan_end"])
                return ps, pe, ss, se
            if "source_pixel_index" in ret and "source_scan_index" in ret:
                pix = np.asarray(ret["source_pixel_index"].values)
                scn = np.asarray(ret["source_scan_index"].values)
                return int(pix.min()), int(pix.max()), int(scn.min()), int(scn.max())

    if all(k in ds.attrs for k in ("pixel_start", "pixel_end", "scan_start", "scan_end")):
        return (
            int(ds.attrs["pixel_start"]),
            int(ds.attrs["pixel_end"]),
            int(ds.attrs["scan_start"]),
            int(ds.attrs["scan_end"]),
        )

    # Full-swath default (1-based)
    _ = src
    return 1, n_pix, 1, n_scan


def _align_l1b_angles(
    angle_full: np.ndarray,
    *,
    target_shape: tuple[int, int],
    pixel_start: int,
    pixel_end: int,
    scan_start: int,
    scan_end: int,
) -> np.ndarray | None:
    """Subset L1B (scans, pixels) angle map onto the retrieval / gain grid.

    ``pixel_*`` / ``scan_*`` are 1-based inclusive indices into the full L1B swath.
    Returns None if alignment is impossible.
    """
    n_scan, n_pix = target_shape
    # Accept either (scans, pixels) or (pixels, scans) if unambiguous
    if angle_full.ndim != 2:
        return None
    if angle_full.shape == target_shape:
        return angle_full.astype(np.float32, copy=False)
    if angle_full.shape == (n_pix, n_scan):
        angle_full = angle_full.T

    if angle_full.shape[0] < scan_end or angle_full.shape[1] < pixel_end:
        # Try full-swath match only
        if angle_full.shape == target_shape:
            return angle_full.astype(np.float32, copy=False)
        return None

    # 1-based inclusive → 0-based slices
    sub = angle_full[scan_start - 1 : scan_end, pixel_start - 1 : pixel_end]
    if sub.shape != target_shape:
        return None
    return sub.astype(np.float32, copy=False)


def _load_aligned_angles(
    l1b_path: Path,
    ds: xr.Dataset,
    target_shape: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray] | None:
    n_scan, n_pix = target_shape
    ps, pe, ss, se = _retrieval_window(ds, n_scan, n_pix)
    with xr.open_dataset(l1b_path, group="geolocation_data") as geo:
        sza_raw = np.asarray(geo["solar_zenith"].values, dtype=np.float32)
        vza_raw = np.asarray(geo["sensor_zenith"].values, dtype=np.float32)
    sza = _align_l1b_angles(
        sza_raw,
        target_shape=target_shape,
        pixel_start=ps,
        pixel_end=pe,
        scan_start=ss,
        scan_end=se,
    )
    vza = _align_l1b_angles(
        vza_raw,
        target_shape=target_shape,
        pixel_start=ps,
        pixel_end=pe,
        scan_start=ss,
        scan_end=se,
    )
    if sza is None or vza is None:
        return None
    return sza, vza


def _quantile_edges(values: np.ndarray, n_bins: int) -> np.ndarray:
    edges = np.unique(np.nanquantile(values, np.linspace(0.0, 1.0, n_bins + 1)))
    if edges.size < 3:
        edges = np.linspace(np.nanmin(values), np.nanmax(values), n_bins + 1)
    return edges


def _subsample(values: np.ndarray, max_n: int, rng: np.random.Generator) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32).ravel()
    values = values[np.isfinite(values)]
    if values.size <= max_n:
        return values
    return rng.choice(values, size=max_n, replace=False)


def _bin_index(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Return bin id in [0, n_bins-1]; -1 if out of range / non-finite."""
    idx = np.digitize(values, edges[1:-1], right=False)
    bad = ~np.isfinite(values) | (values < edges[0]) | (values > edges[-1])
    idx = idx.astype(np.int32)
    idx[bad] = -1
    return idx


def _accumulate_bins(
    sum_g: np.ndarray,
    cnt: np.ndarray,
    gain_2d: np.ndarray,
    bin_ids: np.ndarray,
) -> None:
    """Add spectra into bin accumulators. gain_2d: (λ, N), bin_ids: (N,)."""
    n_bins = sum_g.shape[0]
    for b in range(n_bins):
        m = bin_ids == b
        if not np.any(m):
            continue
        g = gain_2d[:, m]
        finite = np.isfinite(g)
        g0 = np.where(finite, g, 0.0)
        sum_g[b] += g0.sum(axis=1)
        cnt[b] += finite.sum(axis=1)


def _plot_binned(
    wavelength: np.ndarray,
    edges: np.ndarray,
    sum_g: np.ndarray,
    cnt: np.ndarray,
    *,
    title: str,
    legend_title: str,
    cmap_name: str,
    out_path: Path,
    show: bool,
    fmt: str = ".3g",
    unit: str = "",
) -> None:
    n_bins = sum_g.shape[0]
    cmap = getattr(plt.cm, cmap_name)(np.linspace(0.15, 0.9, n_bins))
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for b in range(n_bins):
        ok = cnt[b] > 0
        if not np.any(ok):
            continue
        mean_b = np.full(sum_g.shape[1], np.nan)
        mean_b[ok] = sum_g[b, ok] / cnt[b, ok]
        lo, hi = edges[b], edges[b + 1]
        n_pix = int(np.round(np.nanmean(cnt[b][ok]))) if np.any(ok) else 0
        label = f"[{lo:{fmt}}, {hi:{fmt}}){unit}  n≈{n_pix:,}"
        ax.plot(wavelength, mean_b, color=cmap[b], lw=1.4, label=label)
    ax.axhline(1.0, color="k", ls=":", lw=0.8)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Mean vcal_gain")
    ax.set_title(title)
    ax.legend(title=legend_title, fontsize=8, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)
    print(f"  wrote {out_path}")


def analyze_day(
    vcal_dir: Path,
    date: str,
    out_dir: Path,
    *,
    n_bins: int = 10,
    n_bins_angle: int = 5,
    max_pixel_lines: int = 400,
    max_files: int = 0,
    edge_samples: int = 200_000,
    l1b_dir: Path | None = None,
    seed: int = 0,
    show: bool = False,
) -> None:
    files = _list_day_files(vcal_dir, date)
    if max_files > 0:
        files = files[:max_files]
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    print(f"Found {len(files)} granules for {date} in {vcal_dir}")

    wavelength: np.ndarray | None = None
    n_wl = n_pix = 0

    sum_all: np.ndarray | None = None
    cnt_all: np.ndarray | None = None
    sum_pix: np.ndarray | None = None
    cnt_pix: np.ndarray | None = None

    # Subsampled scalars for quantile edges
    nflh_samp: list[np.ndarray] = []
    chl_samp: list[np.ndarray] = []
    sza_samp: list[np.ndarray] = []
    vza_samp: list[np.ndarray] = []
    per_file_cap = max(edge_samples // max(len(files), 1), 2_000)

    # Paths that passed pass-1 (reuse in pass-2)
    pass1_ok: list[Path] = []

    n_skip_l1b = 0
    n_skip_align = 0
    for i, path in enumerate(files):
        m = _GRANULE_RE.match(path.name)
        if not m:
            continue
        gid = m.group(1)
        try:
            with xr.open_dataset(path) as ds:
                wl = np.asarray(ds["wavelength"].values, dtype=float)
                gain = _as_wl_scans_pixels(ds["vcal_gain"])
                nflh = _as_scans_pixels(ds["nflh"], n_pix=gain.shape[2])
                chlor = _as_scans_pixels(ds["chlor_a"], n_pix=gain.shape[2])
                l1b_path = _resolve_l1b(ds, gid, l1b_dir)
                ds_attrs = dict(ds.attrs)
        except Exception as exc:
            print(f"  [skip] read failed {path.name}: {exc}")
            continue

        n_wl_g, n_scan_g, n_pix_g = gain.shape
        if wavelength is None:
            wavelength = wl
            n_wl, n_pix = n_wl_g, n_pix_g
            sum_all = np.zeros(n_wl, dtype=np.float64)
            cnt_all = np.zeros(n_wl, dtype=np.float64)
            sum_pix = np.zeros((n_wl, n_pix), dtype=np.float64)
            cnt_pix = np.zeros((n_wl, n_pix), dtype=np.float64)
            print(f"  reference grid: λ={n_wl}  pixels={n_pix}  (scans vary per granule)")
        elif wl.shape != wavelength.shape or not np.allclose(wl, wavelength, equal_nan=True):
            print(f"  [skip] wavelength mismatch: {path.name}")
            continue
        elif n_pix_g != n_pix:
            print(f"  [skip] pixel count {n_pix_g} ≠ {n_pix}: {path.name}")
            continue

        assert sum_all is not None and cnt_all is not None
        assert sum_pix is not None and cnt_pix is not None

        finite = np.isfinite(gain)
        sum_all += np.where(finite, gain, 0.0).sum(axis=(1, 2))
        cnt_all += finite.sum(axis=(1, 2))
        sum_pix += np.where(finite, gain, 0.0).sum(axis=1)  # over scans
        cnt_pix += finite.sum(axis=1)

        gain_2d = gain.reshape(n_wl, -1)
        has_gain = np.isfinite(gain_2d).any(axis=0)

        nflh_flat = nflh.ravel()
        m_n = np.isfinite(nflh_flat) & has_gain
        if np.any(m_n):
            nflh_samp.append(_subsample(nflh_flat[m_n], per_file_cap, rng))

        chl_flat = chlor.ravel()
        m_c = np.isfinite(chl_flat) & (chl_flat > 0) & has_gain
        if np.any(m_c):
            chl_samp.append(_subsample(chl_flat[m_c], per_file_cap, rng))

        if l1b_path is None:
            n_skip_l1b += 1
        else:
            tmp = xr.Dataset(attrs=ds_attrs)
            aligned = _load_aligned_angles(l1b_path, tmp, (n_scan_g, n_pix))
            if aligned is None:
                n_skip_align += 1
                print(f"  [warn] cannot align L1B angles for {gid} (L1B vs {n_scan_g}×{n_pix})")
            else:
                sza, vza = aligned
                sza_flat = sza.ravel()
                vza_flat = vza.ravel()
                m_s = np.isfinite(sza_flat) & has_gain
                m_v = np.isfinite(vza_flat) & has_gain
                if np.any(m_s):
                    sza_samp.append(_subsample(sza_flat[m_s], per_file_cap, rng))
                if np.any(m_v):
                    vza_samp.append(_subsample(vza_flat[m_v], per_file_cap, rng))

        pass1_ok.append(path)
        if (i + 1) % 20 == 0 or i + 1 == len(files):
            print(
                f"  pass1 {i + 1}/{len(files)}  ok={len(pass1_ok)}  "
                f"l1b_miss={n_skip_l1b}  align_fail={n_skip_align}"
            )

    if wavelength is None or sum_all is None or not pass1_ok:
        raise RuntimeError("No usable granules")

    n_ok = len(pass1_ok)

    # ── Plot 1: per-pixel day mean + overall mean ────────────────────────────
    mean_all = np.divide(sum_all, cnt_all, out=np.full_like(sum_all, np.nan), where=cnt_all > 0)
    mean_pix = np.divide(sum_pix, cnt_pix, out=np.full_like(sum_pix, np.nan), where=cnt_pix > 0)
    pix_ok = np.flatnonzero(np.isfinite(mean_pix).any(axis=0))
    if pix_ok.size > max_pixel_lines:
        pix_ok = pix_ok[:: max(1, pix_ok.size // max_pixel_lines)]

    cmin, cmax = 1, n_pix
    fig, ax = plt.subplots(figsize=(10, 4))
    cmap = plt.cm.twilight
    norm = plt.Normalize(vmin=cmin, vmax=cmax)
    for j in pix_ok:
        ax.plot(wavelength, mean_pix[:, j], color=cmap(norm(j + 1)), lw=0.5, alpha=0.45)
    ax.plot(wavelength, mean_all, color="silver", lw=2.0, label="day mean", zorder=5)
    ax.axhline(1.0, color="k", ls=":", lw=0.8)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("pixel index")
    cbar.set_ticks([cmin, (cmin + cmax) // 2, cmax])
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Mean vcal_gain")
    ax.set_ylim(0.994, 1.01)
    n_valid = int(np.nansum(cnt_all) / max(n_wl, 1))
    ax.set_title(
        f"{date}  ·  {n_ok} granules  ·  n_valid≈{n_valid:,}  ·  {pix_ok.size}/{n_pix} pixel lines"
    )
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    p1 = out_dir / f"vcal_gain_{date}_per_pixel.png"
    fig.savefig(p1, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)
    print(f"  wrote {p1}")

    def _cat(chunks: list[np.ndarray]) -> np.ndarray:
        if not chunks:
            return np.array([], dtype=np.float32)
        return np.concatenate(chunks)

    nflh_all = _subsample(_cat(nflh_samp), edge_samples, rng)
    chl_all = _subsample(_cat(chl_samp), edge_samples, rng)
    sza_all = _subsample(_cat(sza_samp), edge_samples, rng)
    vza_all = _subsample(_cat(vza_samp), edge_samples, rng)
    print(
        f"  edge samples (≤{edge_samples:,}): nflh={nflh_all.size:,}  chlor_a={chl_all.size:,}  "
        f"sza={sza_all.size:,}  vza={vza_all.size:,}"
    )

    edges = {
        "nflh": _quantile_edges(nflh_all, n_bins) if nflh_all.size else None,
        "chlor_a": _quantile_edges(chl_all, n_bins) if chl_all.size else None,
        "sza": _quantile_edges(sza_all, n_bins_angle) if sza_all.size else None,
        "vza": _quantile_edges(vza_all, n_bins_angle) if vza_all.size else None,
    }

    # Pass 2: only pass1_ok paths; re-check λ / pixel count
    acc: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for key, ed in edges.items():
        if ed is None:
            continue
        nb = ed.size - 1
        acc[key] = (np.zeros((nb, n_wl), dtype=np.float64), np.zeros((nb, n_wl), dtype=np.float64))

    for i, path in enumerate(pass1_ok):
        m = _GRANULE_RE.match(path.name)
        assert m is not None
        gid = m.group(1)
        with xr.open_dataset(path) as ds:
            wl = np.asarray(ds["wavelength"].values, dtype=float)
            if wl.shape != wavelength.shape or not np.allclose(wl, wavelength, equal_nan=True):
                print(f"  [pass2 skip] wavelength changed: {path.name}")
                continue
            gain = _as_wl_scans_pixels(ds["vcal_gain"])
            if gain.shape[2] != n_pix:
                print(f"  [pass2 skip] pixel count mismatch: {path.name}")
                continue
            n_scan_g = gain.shape[1]
            nflh = _as_scans_pixels(ds["nflh"], n_pix=n_pix)
            chlor = _as_scans_pixels(ds["chlor_a"], n_pix=n_pix)
            l1b_path = _resolve_l1b(ds, gid, l1b_dir)
            ds_attrs = dict(ds.attrs)

        gain_2d = gain.reshape(n_wl, -1)
        has_gain = np.isfinite(gain_2d).any(axis=0)

        if "nflh" in acc and edges["nflh"] is not None:
            ids = _bin_index(nflh.ravel(), edges["nflh"])
            ids[~has_gain] = -1
            _accumulate_bins(acc["nflh"][0], acc["nflh"][1], gain_2d, ids)

        if "chlor_a" in acc and edges["chlor_a"] is not None:
            chl = chlor.ravel()
            ids = _bin_index(chl, edges["chlor_a"])
            ids[~(has_gain & np.isfinite(chl) & (chl > 0))] = -1
            _accumulate_bins(acc["chlor_a"][0], acc["chlor_a"][1], gain_2d, ids)

        if l1b_path is not None and ("sza" in acc or "vza" in acc):
            aligned = _load_aligned_angles(
                l1b_path, xr.Dataset(attrs=ds_attrs), (n_scan_g, n_pix)
            )
            if aligned is not None:
                sza, vza = aligned
                if "sza" in acc and edges["sza"] is not None:
                    ids = _bin_index(sza.ravel(), edges["sza"])
                    ids[~has_gain] = -1
                    _accumulate_bins(acc["sza"][0], acc["sza"][1], gain_2d, ids)
                if "vza" in acc and edges["vza"] is not None:
                    ids = _bin_index(vza.ravel(), edges["vza"])
                    ids[~has_gain] = -1
                    _accumulate_bins(acc["vza"][0], acc["vza"][1], gain_2d, ids)

        if (i + 1) % 20 == 0 or i + 1 == n_ok:
            print(f"  pass2 {i + 1}/{n_ok}")

    # ── Binned plots ─────────────────────────────────────────────────────────
    if "nflh" in acc:
        _plot_binned(
            wavelength,
            edges["nflh"],
            *acc["nflh"],
            title=f"{date} · binned by nflh · {n_ok} granules",
            legend_title="nflh bin",
            cmap_name="viridis",
            out_path=out_dir / f"vcal_gain_{date}_by_nflh.png",
            show=show,
        )
    if "chlor_a" in acc:
        _plot_binned(
            wavelength,
            edges["chlor_a"],
            *acc["chlor_a"],
            title=f"{date} · binned by chlor_a · {n_ok} granules",
            legend_title="chlor_a bin (mg m⁻³)",
            cmap_name="plasma",
            out_path=out_dir / f"vcal_gain_{date}_by_chlor_a.png",
            show=show,
        )
    if "sza" in acc:
        _plot_binned(
            wavelength,
            edges["sza"],
            *acc["sza"],
            title=f"{date} · binned by SZA · {n_ok} granules",
            legend_title="SZA bin",
            cmap_name="cividis",
            out_path=out_dir / f"vcal_gain_{date}_by_sza.png",
            show=show,
            fmt=".1f",
            unit="°",
        )
    if "vza" in acc:
        _plot_binned(
            wavelength,
            edges["vza"],
            *acc["vza"],
            title=f"{date} · binned by VZA · {n_ok} granules",
            legend_title="VZA bin",
            cmap_name="magma",
            out_path=out_dir / f"vcal_gain_{date}_by_vza.png",
            show=show,
            fmt=".1f",
            unit="°",
        )

    npz_path = out_dir / f"vcal_gain_{date}_summary.npz"
    np.savez_compressed(
        npz_path,
        wavelength=wavelength,
        mean_gain=mean_all,
        mean_gain_per_pixel=mean_pix,
        n_granules=np.int32(n_ok),
        n_pixels=np.int32(n_pix),
        nflh_edges=edges["nflh"] if edges["nflh"] is not None else np.array([]),
        chlor_a_edges=edges["chlor_a"] if edges["chlor_a"] is not None else np.array([]),
        sza_edges=edges["sza"] if edges["sza"] is not None else np.array([]),
        vza_edges=edges["vza"] if edges["vza"] is not None else np.array([]),
    )
    print(f"  wrote {npz_path}")
    print("Done.")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--date", default="20250125", help="YYYYMMDD (default: 20250125)")
    p.add_argument(
        "--vcal-dir",
        type=Path,
        default=Path("/home/zhe2/data/PACE/new_svd_retrieval_output/vcal_gain"),
        help="Directory of vcal_gain_interim_*.nc files",
    )
    p.add_argument(
        "--l1b-dir",
        type=Path,
        default=None,
        help="Optional L1B root (YYYY/MM/DD/ or flat). Overrides NetCDF l1b_file when found.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Figure / npz output directory (default: <vcal-dir>/qc_<date>)",
    )
    p.add_argument("--n-bins", type=int, default=10, help="Bins for nflh / chlor_a")
    p.add_argument("--n-bins-angle", type=int, default=5, help="Bins for SZA / VZA")
    p.add_argument(
        "--edge-samples",
        type=int,
        default=200_000,
        help="Max scalar samples kept for quantile edges (subsampled)",
    )
    p.add_argument("--max-pixel-lines", type=int, default=400)
    p.add_argument("--max-files", type=int, default=0, help="If >0, only process the first N files (smoke)")
    p.add_argument("--seed", type=int, default=0, help="RNG seed for edge subsampling")
    p.add_argument("--show", action="store_true", help="Also display figures interactively")
    args = p.parse_args(argv)

    out_dir = args.out_dir or (args.vcal_dir / f"qc_{args.date}")
    analyze_day(
        args.vcal_dir,
        args.date,
        out_dir,
        n_bins=args.n_bins,
        n_bins_angle=args.n_bins_angle,
        max_pixel_lines=args.max_pixel_lines,
        max_files=args.max_files,
        edge_samples=args.edge_samples,
        l1b_dir=args.l1b_dir,
        seed=args.seed,
        show=args.show,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
