#!/usr/bin/env python3
"""
Download Sentinel-5P TROPOMI SNPP/VIIRS cloud product band 6 (S5P_OFFL_L2__NP_BD6)
via NASA earthaccess / GES DISC.

Files are placed under <output_dir>/<YYYY>/ according to each granule's
sensing-start year, e.g.:
  /kiwi-data/.../L2_VIIRS/2025/S5P_OFFL_L2__NP_BD6_20250702T….nc

Usage:
  python download_tropomi_np_bd6.py <config.toml>

Requires Python 3.11+ (stdlib tomllib) and: pip install earthaccess

Config layout (TOML):
  [download.temporal]
  start = "2025-07-02"
  end   = "2025-07-08"

  [download.spatial]               # optional
  use_bounding_box = true
  bounding_box = [-180, -90, 180, 90]   # [west, south, east, north]

  [download.collections]           # optional — defaults shown
  short_name = "S5P_OFFL_L2__NP_BD6"
  versions   = ["2", "3"]          # list of versions to search; [""] to omit filter

  [download.paths]
  output_dir = "/kiwi-data/Data/groupMembers/zhe2/MyProjects/TROPOMI_SIF/ESA/L2_VIIRS"

  [download.options]
  skip_existing  = true    # skip files already on disk (checked by filename + size)
  login_strategy = "all"   # environment → netrc → interactive
"""
from __future__ import annotations

import inspect
import os
import re
import sys
import tomllib
from collections import defaultdict
from pathlib import Path


def _die(msg: str, code: int = 1) -> None:
    print(msg, file=sys.stderr)
    raise SystemExit(code)


def _get_str(d: dict, key: str, default: str = "") -> str:
    v = d.get(key, default)
    return default if v is None else str(v).strip()


def _build_search_kwargs(
    *,
    short_name: str,
    version: str,
    temporal: tuple[str, str],
    use_bbox: bool,
    bbox: tuple[float, float, float, float] | None,
) -> dict:
    kw: dict = {"short_name": short_name, "temporal": temporal}
    v = version.strip()
    if v:
        kw["version"] = v
    if use_bbox and bbox is not None:
        kw["bounding_box"] = bbox
    return kw


def _parse_versions(col: dict) -> list[str]:
    """Return a list of version strings from the config.

    Accepts:
      versions = ["2", "3"]   # preferred — list of strings
      version  = "2"          # single string (legacy)
      versions = [""]         # empty string → omit version filter
    """
    raw = col.get("versions")
    if raw is not None:
        if not isinstance(raw, list):
            _die("[download.collections] versions must be an array, e.g. versions = [\"2\", \"3\"]")
        return [str(v).strip() for v in raw]
    # fallback to singular key
    v = str(col.get("version", "2")).strip()
    return [v]


def _download_force_kw(download_fn, force: bool) -> dict:
    sig = inspect.signature(download_fn)
    if "force" not in sig.parameters:
        return {}
    return {"force": force}


# S5P filename pattern: S5P_OFFL_L2__NP_BD6_<YYYYMMDDTHHmmss>_…
_S5P_YEAR_RE = re.compile(r"S5P_\w+_(\d{4})\d{6}T")


def _granule_year(granule) -> str | None:
    """Extract the 4-digit sensing year from the granule's data link filename."""
    links = granule.data_links() if hasattr(granule, "data_links") else []
    for link in links:
        fname = Path(link).name
        m = _S5P_YEAR_RE.search(fname)
        if m:
            return m.group(1)
    # Fallback: try the granule temporal metadata
    try:
        t = granule["umm"]["TemporalExtent"]["RangeDateTime"]["BeginningDateTime"]
        return t[:4]
    except Exception:
        return None


def _group_by_year(granules: list) -> dict[str, list]:
    """Return {year: [granule, …]} mapping."""
    by_year: dict[str, list] = defaultdict(list)
    for g in granules:
        year = _granule_year(g) or "unknown"
        by_year[year].append(g)
    return dict(by_year)


def _filter_existing(granules: list, year_dir: Path) -> tuple[list, int]:
    """Drop granules whose file already exists in *year_dir* with non-zero size."""
    if not year_dir.exists():
        return granules, 0
    existing = {
        p.name for p in year_dir.iterdir() if p.is_file() and p.stat().st_size > 0
    }
    to_download, n_skipped = [], 0
    for g in granules:
        links = g.data_links() if hasattr(g, "data_links") else []
        fname = Path(links[0]).name if links else None
        if fname and fname in existing:
            n_skipped += 1
        else:
            to_download.append(g)
    return to_download, n_skipped


def main() -> None:
    if len(sys.argv) != 2:
        _die("usage: download_tropomi_np_bd6.py <config.toml>")
    cfg_path = Path(sys.argv[1]).expanduser().resolve()
    if not cfg_path.is_file():
        _die(f"config not found: {cfg_path}")

    with open(cfg_path, "rb") as f:
        cfg = tomllib.load(f)

    dl = cfg.get("download")
    if not isinstance(dl, dict):
        _die(
            "Missing [download] table. Expected keys: "
            "[download.temporal], [download.spatial], [download.collections], "
            "[download.paths], [download.options]"
        )

    # --- temporal ---
    temporal = dl.get("temporal") or {}
    t0 = _get_str(temporal, "start")
    t1 = _get_str(temporal, "end")
    if not t0 or not t1:
        _die("[download.temporal] start and end are required (YYYY-MM-DD)")

    # --- spatial ---
    spatial = dl.get("spatial") or {}
    use_bbox = bool(spatial.get("use_bounding_box", False))
    bb_raw = spatial.get("bounding_box")
    bbox: tuple[float, float, float, float] | None = None
    if use_bbox:
        if not isinstance(bb_raw, list) or len(bb_raw) != 4:
            _die(
                "[download.spatial] bounding_box must be [west, south, east, north] "
                "when use_bounding_box=true"
            )
        bbox = tuple(float(x) for x in bb_raw)

    # --- collections ---
    col = dl.get("collections") or {}
    short_name = _get_str(col, "short_name", "S5P_OFFL_L2__NP_BD6")
    versions   = _parse_versions(col)

    # --- paths ---
    paths = dl.get("paths") or {}
    out_s = _get_str(paths, "output_dir")
    if not out_s:
        _die("[download.paths] output_dir is required")
    base_dir = Path(os.path.expanduser(out_s))

    # --- options ---
    opts = dl.get("options") or {}
    skip_existing = opts.get("skip_existing", True)
    if not isinstance(skip_existing, bool):
        _die("[download.options] skip_existing must be a boolean")
    force_redownload = not skip_existing

    login_strategy = _get_str(opts, "login_strategy", "all").lower()
    allowed_login = {"all", "environment", "netrc", "interactive"}
    if login_strategy not in allowed_login:
        _die(
            f"[download.options] login_strategy must be one of: "
            f"{', '.join(sorted(allowed_login))}"
        )

    # --- earthaccess ---
    try:
        import earthaccess
    except ImportError:
        _die("earthaccess not installed (pip install earthaccess)")

    try:
        auth = earthaccess.login(strategy=login_strategy)
    except Exception as e:
        _die(
            "earthaccess.login failed.\n"
            f"  login_strategy={login_strategy!r}\n"
            "Credentials:\n"
            "  • environment: EARTHDATA_USERNAME + EARTHDATA_PASSWORD, or EARTHDATA_TOKEN\n"
            "  • netrc: ~/.netrc (or NETRC) with machine urs.earthdata.nasa.gov\n"
            "  • all: environment → netrc → interactive\n"
            f"Details: {type(e).__name__}: {e}"
        )
    if not getattr(auth, "authenticated", False):
        _die(
            "earthaccess login did not authenticate.\n"
            "Check credentials or set [download.options] login_strategy."
        )

    # Search each requested version separately and merge results.
    granules: list = []
    ver_label = ", ".join(f"v{v}" if v else "any" for v in versions)
    print(f"[NP_BD6] Searching {short_name} [{ver_label}] from {t0} to {t1} ...")
    for ver in versions:
        kw = _build_search_kwargs(
            short_name=short_name,
            version=ver,
            temporal=(t0, t1),
            use_bbox=use_bbox,
            bbox=bbox,
        )
        try:
            found = earthaccess.search_data(**kw)
        except Exception as e:
            _die(f"[NP_BD6] search_data (version={ver!r}) failed: {e}")
        n_ver = len(found) if found else 0
        label = f"v{ver}" if ver else "any version"
        print(f"[NP_BD6]   {label}: {n_ver} granule(s)")
        if found:
            granules.extend(found)

    # Deduplicate by data URL in case version ranges overlap.
    seen: set[str] = set()
    unique: list = []
    for g in granules:
        links = g.data_links() if hasattr(g, "data_links") else []
        key = links[0] if links else id(g)
        if key not in seen:
            seen.add(key)
            unique.append(g)
    granules = unique

    print(f"[NP_BD6] Total: {len(granules)} unique granule(s)  base_dir={base_dir}")
    if not granules:
        print("Nothing to download.")
        return

    # Group by year → download into <base_dir>/<YYYY>/
    by_year = _group_by_year(granules)
    dl_kw = _download_force_kw(earthaccess.download, force_redownload)

    total_downloaded = 0
    total_skipped = 0

    for year in sorted(by_year):
        year_dir = base_dir / year
        year_dir.mkdir(parents=True, exist_ok=True)
        year_granules = by_year[year]

        if skip_existing:
            year_granules, n_skipped = _filter_existing(year_granules, year_dir)
            total_skipped += n_skipped
            if n_skipped:
                print(f"[NP_BD6] {year}/  skipping {n_skipped} already-downloaded granule(s).")

        if not year_granules:
            print(f"[NP_BD6] {year}/  all present, nothing to download.")
            continue

        print(f"[NP_BD6] {year}/  downloading {len(year_granules)} granule(s) → {year_dir}")
        try:
            files = earthaccess.download(year_granules, local_path=str(year_dir), **dl_kw)
        except Exception as e:
            _die(f"[NP_BD6] {year}/ download failed: {e}")

        if files is None:
            files = []
        # earthaccess creates files with 0o600; fix to group-readable 0o644.
        for fp in files:
            try:
                Path(fp).chmod(0o644)
            except Exception:
                pass
        total_downloaded += len(files)
        print(f"[NP_BD6] {year}/  {len(files)} file(s) written (permissions set to 644).")

    print(
        f"\n[NP_BD6] Summary: {total_downloaded} downloaded, "
        f"{total_skipped} skipped (already on disk)."
    )
    print("Done.")


if __name__ == "__main__":
    main()
