#!/usr/bin/env bash
# Fixed-SIF EV1 ensemble: identical SIF shape (retrieval basis EV1), varying strength only.
# Retrieval stores water-leaving SIF @ 678 nm (not SIF×T₁).
set -euo pipefail
cd /home/zhe2/FraLab/PACE_redSIF_PACE
export OUT_SUFFIX=_new
export FIXED_SIF_EV1=1
export GKSwstype=100
LOGDIR=surrogate_meas/batch_ensemble
mkdir -p "$LOGDIR"

echo "======== $(date -Is) START npoly3_fixedSIF_new (generate truth) ========"
REUSE_TRUTH=0 ZERO_SIF=0 \
  SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly3.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl \
  2>&1 | tee "$LOGDIR/run_npoly3_fixedSIF_new.log"
echo "======== $(date -Is) DONE  npoly3_fixedSIF_new ========"

echo "======== $(date -Is) START npoly5_fixedSIF_new (reuse truth) ========"
REUSE_TRUTH=1 ZERO_SIF=0 \
  SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly5.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl \
  2>&1 | tee "$LOGDIR/run_npoly5_fixedSIF_new.log"
echo "======== $(date -Is) DONE  npoly5_fixedSIF_new ========"

echo "All fixed-SIF EV1 ensemble runs finished at $(date -Is)"
