#!/usr/bin/env bash
# Re-run batch_ensemble with decay ∈ [4,12] and OUT_SUFFIX=_new
set -euo pipefail
cd /home/zhe2/FraLab/PACE_redSIF_PACE
export OUT_SUFFIX=_new
export GKSwstype=100
LOGDIR=surrogate_meas/batch_ensemble
mkdir -p "$LOGDIR"

echo "======== $(date -Is) START npoly3_new (reuse truth) ========"
REUSE_TRUTH=1 ZERO_SIF=0 \
  SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly3.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl \
  2>&1 | tee "$LOGDIR/run_npoly3_new.log"
echo "======== $(date -Is) DONE  npoly3_new ========"

echo "======== $(date -Is) START npoly5_new (reuse truth) ========"
REUSE_TRUTH=1 ZERO_SIF=0 \
  SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly5.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl \
  2>&1 | tee "$LOGDIR/run_npoly5_new.log"
echo "======== $(date -Is) DONE  npoly5_new ========"

echo "======== $(date -Is) START npoly3_zeroSIF_new (reuse truth) ========"
REUSE_TRUTH=1 ZERO_SIF=1 \
  SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly3.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl \
  2>&1 | tee "$LOGDIR/run_npoly3_zeroSIF_new.log"
echo "======== $(date -Is) DONE  npoly3_zeroSIF_new ========"

echo "======== $(date -Is) START npoly5_zeroSIF_new (reuse zeroSIF truth) ========"
REUSE_TRUTH=1 ZERO_SIF=1 \
  SVD_CONFIG=surrogate_meas/configs/svd_nPC15_npoly5.toml \
  julia --project=. -t 8 surrogate_meas/run_ensemble.jl \
  2>&1 | tee "$LOGDIR/run_npoly5_zeroSIF_new.log"
echo "======== $(date -Is) DONE  npoly5_zeroSIF_new ========"

echo "All _new ensemble runs finished at $(date -Is)"
