#!/usr/bin/env bash
# Run SVD retrieval for download.run.Sundarbans.toml; logs under repo logs/
set -euo pipefail

REPO="/home/zhe2/FraLab/PACE_redSIF_PACE"
CONFIG="${REPO}/svd_configs/download.run.Sundarbans.toml"
JULIA="${JULIA:-/home/zhe2/.juliaup/bin/julia}"
LOG_DIR="${REPO}/logs"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="${LOG_DIR}/svd_retrieval_sundarbans_${STAMP}.log"

mkdir -p "${LOG_DIR}"
cd "${REPO}"

rc=0
{
  echo "=== start $(date -Is) ==="
  echo "host=$(hostname)  user=$(whoami)  pwd=$(pwd)"
  echo "julia=$("${JULIA}" --version)"
  echo "config=${CONFIG}"
  echo "threads=auto"
  echo "---"
  set +e
  "${JULIA}" --project=. -t auto \
    global_svd_fit_pipeline/svd_retrieval/run_svd_fit.jl \
    "${CONFIG}"
  rc=$?
  set -e
  echo "=== exit ${rc} $(date -Is) ==="
} >> "${LOG}" 2>&1

echo "Log: ${LOG}"
exit "${rc}"
