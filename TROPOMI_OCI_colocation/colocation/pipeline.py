"""Orchestrate co-location stages with caching."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import cache
from .config import Config, simulate_fingerprint_payload
from .match import run_matches_stage
from .pair import run_pairs_stage
from .rsr import run_rsr_stage
from .simulate import run_simulate
from .svd_train import run_svd_stage
from .write_nc import write_matchup_nc


def run_pipeline(
    cfg: Config,
    *,
    force: set[str] | None = None,
    from_stage: str | None = None,
    log=print,
) -> Path:
    force = force or set()
    dirs = cache.ensure_output_dirs(cfg.output_dir)
    cache.write_run_meta(cfg.output_dir, cfg.to_plain_dict(), extra={"stages": list(cache.STAGES)})

    start_i = 0
    if from_stage is not None:
        start_i = cache.stage_order_index(from_stage)

    def need(stage: str) -> bool:
        return cache.stage_order_index(stage) >= start_i

    # ── pairs ────────────────────────────────────────────────────────────────
    if need("pairs"):
        swath_list, pairs_fp = run_pairs_stage(
            cfg, dirs, force=("pairs" in force), log=log
        )
    else:
        list_path = dirs["pairs"] / "swath_list.json"
        meta_path = dirs["pairs"] / "pairs_meta.json"
        if not list_path.is_file() or not meta_path.is_file():
            raise RuntimeError("--from skipped pairs but pairs cache is missing")
        swath_list = cache.read_json(list_path)
        pairs_fp = cache.read_json(meta_path)["fingerprint"]
        log(f"[pairs] loaded existing cache (from={from_stage})")

    # ── matches ──────────────────────────────────────────────────────────────
    if need("matches"):
        matches, matches_fp = run_matches_stage(
            cfg,
            dirs,
            swath_list,
            pairs_fp,
            force=("matches" in force) or ("pairs" in force),
            log=log,
        )
    else:
        npz_path = dirs["matches"] / "matches.npz"
        meta_path = dirs["matches"] / "matches_meta.json"
        paths_json = dirs["matches"] / "match_paths.json"
        if not all(p.is_file() for p in (npz_path, meta_path, paths_json)):
            raise RuntimeError("--from skipped matches but matches cache is missing")
        matches = cache.load_npz(npz_path)
        path_info = cache.read_json(paths_json)
        matches["pace_path"] = np.asarray(path_info["pace_path"])
        matches["tropomi_path"] = np.asarray(path_info["tropomi_path"])
        matches["pace_name"] = np.asarray(path_info["pace_name"])
        matches_fp = cache.read_json(meta_path)["fingerprint"]
        log(f"[matches] loaded existing cache (from={from_stage})")

    tropomi_paths = sorted(set(str(p) for p in matches["tropomi_path"]))

    # ── svd ──────────────────────────────────────────────────────────────────
    if need("svd"):
        trop_svd, svd_fp = run_svd_stage(
            cfg,
            dirs,
            tropomi_paths,
            matches_fp,
            force=("svd" in force) or ("matches" in force) or ("pairs" in force),
            log=log,
        )
    else:
        npz_path = dirs["svd"] / "trop_svd.npz"
        meta_path = dirs["svd"] / "svd_meta.json"
        if not npz_path.is_file() or not meta_path.is_file():
            raise RuntimeError("--from skipped svd but svd cache is missing")
        arrays = cache.load_npz(npz_path)
        meta = cache.read_json(meta_path)
        trop_svd = {
            "wl": arrays["wl"],
            "mean": arrays["mean"],
            "s": arrays["s"],
            "Vt": arrays["Vt"],
            "var_frac": arrays["var_frac"],
            "n_train": int(meta["n_train"]),
            "mean_center": bool(meta["mean_center"]),
            "wl_min": float(meta["wl_min"]),
            "wl_max": float(meta["wl_max"]),
            "n_per_file": meta.get("n_per_file", {}),
        }
        svd_fp = meta["fingerprint"]
        log(f"[svd] loaded existing cache (from={from_stage})")

    # ── rsr ──────────────────────────────────────────────────────────────────
    if need("rsr"):
        rsr_bundles, rsr_fp = run_rsr_stage(
            cfg,
            dirs,
            tropomi_paths,
            matches_fp,
            force=("rsr" in force) or ("matches" in force) or ("pairs" in force),
            log=log,
        )
    else:
        # reload via stage helper with force=False and injected fingerprint match
        rsr_bundles, rsr_fp = run_rsr_stage(
            cfg, dirs, tropomi_paths, matches_fp, force=False, log=log
        )

    # ── simulate + NC (always when need simulate) ────────────────────────────
    if not need("simulate"):
        out_nc = dirs["products"] / "matchup.nc"
        if not out_nc.is_file():
            raise RuntimeError("Nothing to do: --from beyond simulate and product missing")
        log(f"[simulate] skipped; existing product {out_nc}")
        return out_nc

    # Force rebuild of product when simulate forced or upstream forced
    sim = run_simulate(cfg, matches, trop_svd, rsr_bundles, log=log)
    fingerprints = {
        "pairs": pairs_fp,
        "matches": matches_fp,
        "svd": svd_fp,
        "rsr": rsr_fp,
        "simulate": cache.fingerprint(
            {
                "svd_fp": svd_fp,
                "rsr_fp": rsr_fp,
                "matches_fp": matches_fp,
                "sim": simulate_fingerprint_payload(cfg),
            }
        ),
    }
    out_nc = dirs["products"] / "matchup.nc"
    write_matchup_nc(
        out_nc,
        cfg,
        matches,
        sim,
        fingerprints=fingerprints,
        svd_meta={"n_train": trop_svd.get("n_train")},
    )
    log(f"[products] wrote {out_nc}")
    return out_nc
