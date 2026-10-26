"""Shared filter for the nflh_vs_sif plots: L2 flag cut AND pandas query on the extracted soundings.

    --l2 zero          l2_flags == 0 (no L2 flag set at all)
    --l2 any           no L2 flag cut
    --l2 NAME,NAME     these flags must be off (names as in l2_flags:flag_meanings)
    --query "..."      any pandas query on the columns; "" for none
"""
import re
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

AOP_DIR = "/kiwi-data/Data/satellite/PACE_OCI/L2_AOP_V3.2"
DATA_DIR = "/home/zhe2/data/PACE/nflh_vs_sif_soundings/new_svd_retrieval_output_O2A"
DEFAULT_L2, DEFAULT_QUERY = "zero", "converged == 1"
UNIT = "W m$^{-2}$ µm$^{-1}$ sr$^{-1}$"


def add_filter_args(p):
    p.add_argument("--dates", nargs="+", default=["20250115", "20250415"])
    p.add_argument("--data-dir", default=DATA_DIR)
    p.add_argument("--l2", default=DEFAULT_L2)
    p.add_argument("--query", default=DEFAULT_QUERY)
    p.add_argument("--tag", default=None, help="output subfolder (default: built from the filter)")


def out_dir(here, a):
    tag = a.tag or re.sub(r"[^A-Za-z0-9_.=<>-]+", "_", f"l2{a.l2}_{a.query or 'noquery'}").strip("_")
    tag = tag.replace("<", "lt").replace(">", "gt").replace("==", "eq").replace("=", "")
    d = Path(here) / "out" / tag
    d.mkdir(parents=True, exist_ok=True)
    return d


def flag_masks(granule):
    """l2_flags name -> bit mask, read from one AOP file of the date."""
    f = f"{AOP_DIR}/{granule[:4]}/{granule[4:6]}/{granule[6:8]}/PACE_OCI.{granule}.L2.OC_AOP.V3_2.nc"
    with netCDF4.Dataset(f) as d:
        v = d["geophysical_data"]["l2_flags"]
        return dict(zip(v.flag_meanings.split(), (np.array(v.flag_masks).astype(np.int64) & 0xFFFFFFFF).tolist()))


def load_filtered(a, date, columns=None):
    """Soundings of one date after the --l2 / --query filter, and the counts after each step."""
    df = pd.read_parquet(Path(a.data_dir) / f"soundings_{date}.parquet", columns=columns)
    steps = [("soundings with L2 nflh", len(df))]
    fl = df.l2_flags.to_numpy().astype(np.int64)
    if a.l2 == "zero":
        m = fl == 0
        steps.append(("l2_flags == 0", int(m.sum())))
    elif a.l2 == "any":
        m = np.ones(len(df), bool)
    else:
        masks = flag_masks(str(df.granule.iloc[0]))
        bad = 0
        for nm in a.l2.split(","):
            bad |= masks[nm.strip()]
        m = (fl & bad) == 0
        steps.append((f"flags off: {a.l2}", int(m.sum())))
    df = df[m]
    if a.query:
        df = df.query(a.query)
        steps.append((a.query, len(df)))
    df = df[np.isfinite(df.sif_678) & np.isfinite(df.nflh)]
    steps.append(("finite SIF and nFLH", len(df)))
    return df, steps


def filter_text(a):
    return f"filter: --l2 {a.l2}   --query \"{a.query}\""
