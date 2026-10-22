"""Pick pixels for averaging-kernel (AK) analysis of one retrieved granule.

    python retrieval_comparison/select_ak_pixels.py [--run npoly3] [--granule 20250724T170203]
                                                     [--n 30] [--converged-in npoly3,newSVD_npoly3]

Writes out/ak/selected_pixels.csv for ak_selected_pixels.jl.

Selection: ocean, L2 flags ok (no cloud/glint/...), chl valid, is_dark==1 (clear open
water), plus converged in every run listed in --converged-in (default: --run), so the
same pixels can be compared across retrieval configurations.
Categories: high_chl (>=1), mid_chl (0.1-0.3), low_chl (<0.05) -- converged;
            nonconverged_lowchl -- clear, dark, low chl but NOT converged in --run.
Indices are 1-based (scan, pixel) as used by the Julia pipeline.
"""
import argparse
from pathlib import Path

import netCDF4 as nc
import numpy as np

RUN_DIRS = {
    "npoly2": "new_svd_retrieval_output_npoly2",
    "npoly3": "new_svd_retrieval_output_O2A",
    "npoly5": "new_svd_retrieval_output_npoly5",
    "newSVD_npoly2": "new_svd_retrieval_output_O2A_newSVD_npoly2",
    "newSVD_npoly3": "new_svd_retrieval_output_O2A_newSVD",
}
PACE = Path("/home/zhe2/data/PACE")
KIWI = Path("/kiwi-data/Data/satellite/PACE_OCI")
BAD_FLAGS = ["ATMFAIL", "LAND", "HIGLINT", "HILT", "HISATZEN", "STRAYLIGHT", "CLDICE", "HISOLZEN", "CHLFAIL"]

ap = argparse.ArgumentParser()
ap.add_argument("--run", default="npoly3", choices=list(RUN_DIRS))
ap.add_argument("--granule", default="20250724T170203")
ap.add_argument("--n", type=int, default=30, help="pixels drawn per category (4 categories -> 120)")
ap.add_argument("--converged-in", default=None, help="comma list of runs that must all converge")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
gid = a.granule
date = f"{gid[:4]}/{gid[4:6]}/{gid[6:8]}"
need = (a.converged_in or a.run).split(",")


def retrieval(run):
    ds = nc.Dataset(PACE / RUN_DIRS[run] / date / f"interim_{gid}_svd_retrieval.nc")
    ds.set_auto_mask(False)
    out = {k: ds[k][:] for k in ["latitude", "longitude", "converged", "is_dark", "is_ocean", "sif_ev1"]}
    ds.close()
    return out


R = {r: retrieval(r) for r in dict.fromkeys(need + [a.run])}
ref = R[a.run]
v = "V3_2"
with nc.Dataset(KIWI / "L2_BGC_V3.2" / date / f"PACE_OCI.{gid}.L2.OC_BGC.{v}.nc") as ds:
    g = ds["geophysical_data"]
    g.set_auto_mask(True)
    chl = np.ma.filled(g["chlor_a"][:].astype(float), np.nan)
    fl = g["l2_flags"]
    flags = np.ma.filled(fl[:], 0).astype(np.int64)
    names, masks = fl.getncattr("flag_meanings").split(), np.atleast_1d(fl.getncattr("flag_masks"))
with nc.Dataset(KIWI / "L2_AOP_V3.2" / date / f"PACE_OCI.{gid}.L2.OC_AOP.{v}.nc") as ds:
    g = ds["geophysical_data"]
    g.set_auto_mask(True)
    nflh = np.ma.filled(g["nflh"][:].astype(float), np.nan)
bad = sum(int(masks[names.index(n)]) for n in BAD_FLAGS if n in names)
chl[chl <= 0] = np.nan

clear = (ref["is_ocean"] == 1) & (ref["is_dark"] == 1) & ((flags & bad) == 0) & np.isfinite(chl)
conv_all = np.ones_like(clear)
for r in need:
    conv_all &= R[r]["converged"] == 1
cats = {
    "high_chl": clear & conv_all & (chl >= 1.0),
    "mid_chl": clear & conv_all & (chl >= 0.1) & (chl < 0.3),
    "low_chl": clear & conv_all & (chl < 0.05),
    "nonconverged_lowchl": clear & (ref["converged"] != 1) & (chl < 0.1),
}
rng = np.random.default_rng(a.seed)
rows = []
for cat, m in cats.items():
    idx = np.argwhere(m)
    print(f"{cat:22s} candidates: {len(idx)}")
    if len(idx) == 0:
        continue
    pick = idx[rng.choice(len(idx), size=min(a.n, len(idx)), replace=False)]
    for s, p in pick:
        rows.append((cat, s + 1, p + 1, ref["latitude"][s, p], ref["longitude"][s, p], chl[s, p], nflh[s, p],
                     ref["sif_ev1"][s, p]))
out = Path(__file__).resolve().parent / "out" / "ak"
out.mkdir(parents=True, exist_ok=True)
with open(out / "selected_pixels.csv", "w") as f:
    f.write("category,scan,pixel,lat,lon,chl,nflh,sif_ev1_stored\n")
    for r in rows:
        f.write(f"{r[0]},{r[1]},{r[2]},{r[3]:.4f},{r[4]:.4f},{r[5]:.5f},{r[6]:.5f},{r[7]:.6f}\n")
print(f"wrote {out / 'selected_pixels.csv'}  ({len(rows)} pixels, converged required in {need})")
