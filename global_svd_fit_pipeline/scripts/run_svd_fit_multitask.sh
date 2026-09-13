#!/usr/bin/env bash
# Launch many Julia SVD retrieval processes that share ONE pipeline TOML.
#
# Each worker gets the same physics / paths; only --date or --granule differs.
# Use this instead of [svd_retrieval].parallel_granules (unsafe with NetCDF/HDF5).
#
# Usage (from repo root):
#   ./global_svd_fit_pipeline/scripts/run_svd_fit_multitask.sh path/to/pipeline.toml
#
# Environment (optional):
#   MAX_JOBS=4              # concurrent Julia processes (default: 4)
#   JULIA_THREADS=8         # -t for each process (default: 8; use 2 with GPU)
#   SPLIT=date|granule      # one job per calendar day (default) or per L1B granule
#   LOG_DIR=...             # worker logs (default: <output_dir>/multitask_logs)
#   DRY_RUN=1               # print commands only
#   JULIA_BIN=julia         # julia executable
#   JULIA_PROJECT=.         # --project path
#
# Examples:
#   MAX_JOBS=6 JULIA_THREADS=8 \
#     ./global_svd_fit_pipeline/scripts/run_svd_fit_multitask.sh \
#       svd_configs/global_fit_new/new.Jan25global_fit_pipeline.npoly3nPC15.toml
#
#   SPLIT=granule MAX_JOBS=8 JULIA_THREADS=4 \
#     ./global_svd_fit_pipeline/scripts/run_svd_fit_multitask.sh my_pipeline.toml

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${PIPE_DIR}/.." && pwd)"
RUN_JL="${PIPE_DIR}/svd_retrieval/run_svd_fit.jl"

CONFIG="${1:-}"
if [[ -z "${CONFIG}" || "${CONFIG}" == "-h" || "${CONFIG}" == "--help" ]]; then
  sed -n '2,30p' "$0"
  exit 0
fi
if [[ ! -f "${CONFIG}" ]]; then
  # allow relative to repo root
  if [[ -f "${REPO_ROOT}/${CONFIG}" ]]; then
    CONFIG="${REPO_ROOT}/${CONFIG}"
  else
    echo "ERROR: config not found: ${CONFIG}" >&2
    exit 1
  fi
fi
CONFIG="$(cd "$(dirname "${CONFIG}")" && pwd)/$(basename "${CONFIG}")"

MAX_JOBS="${MAX_JOBS:-4}"
JULIA_THREADS="${JULIA_THREADS:-8}"
SPLIT="${SPLIT:-date}"          # date | granule
DRY_RUN="${DRY_RUN:-0}"
JULIA_BIN="${JULIA_BIN:-julia}"
JULIA_PROJECT="${JULIA_PROJECT:-${REPO_ROOT}}"

cd "${REPO_ROOT}"

# ── parse dates / L1B root from TOML (stdlib tomllib on 3.11+, else julia) ──
mapfile -t TASKS < <(
  CONFIG="${CONFIG}" SPLIT="${SPLIT}" python3 - <<'PY' 2>/dev/null || true
import os, re, sys
from pathlib import Path
from datetime import datetime, timedelta

cfg_path = Path(os.environ["CONFIG"])
split = os.environ.get("SPLIT", "date")
text = cfg_path.read_text()

def table(name):
    m = re.search(rf"^\[{re.escape(name)}\]\s*$", text, re.M)
    if not m:
        return ""
    start = m.end()
    m2 = re.search(r"^\[", text[start:], re.M)
    return text[start: start + m2.start()] if m2 else text[start:]

def get_str(block, key, default=""):
    m = re.search(rf'^\s*{re.escape(key)}\s*=\s*"([^"]*)"', block, re.M)
    return m.group(1) if m else default

def get_arr(block, key):
    m = re.search(rf'^\s*{re.escape(key)}\s*=\s*\[([^\]]*)\]', block, re.M | re.S)
    if not m:
        return []
    return re.findall(r'"([^"]+)"', m.group(1))

gen = table("general")
svd = table("svd_retrieval")
l1b = get_str(gen, "L1B_dir")
out_dir = get_str(svd, "output_dir", "")

dates = get_arr(svd, "dates")
if not dates:
    ds, de = get_str(svd, "date_start"), get_str(svd, "date_end")
    if not ds:
        ds = get_str(svd, "date")
        de = ds
    if ds and de:
        step = 1
        m = re.search(r'^\s*date_interval_days\s*=\s*(\d+)', svd, re.M)
        if m:
            step = max(1, int(m.group(1)))
        d0 = datetime.strptime(ds, "%Y%m%d").date()
        d1 = datetime.strptime(de, "%Y%m%d").date()
        cur = d0
        while cur <= d1:
            dates.append(cur.strftime("%Y%m%d"))
            cur += timedelta(days=step)

if not dates:
    sys.stderr.write("ERROR: could not resolve dates from [svd_retrieval] "
                     "(need dates= [...] or date_start/date_end or date)\n")
    sys.exit(2)

print(f"#META out_dir={out_dir}", file=sys.stderr)
print(f"#META l1b_dir={l1b}", file=sys.stderr)
print(f"#META n_dates={len(dates)}", file=sys.stderr)

if split == "date":
    for d in dates:
        print(d)
    sys.exit(0)

# granule split: discover L1B under YYYY/MM/DD or flat
if not l1b or not Path(l1b).is_dir():
    sys.stderr.write(f"ERROR: [general].L1B_dir missing or not a dir: {l1b!r}\n")
    sys.exit(2)
l1b = Path(l1b)
gids = []
for d in dates:
    day = l1b / d[0:4] / d[4:6] / d[6:8]
    roots = [day] if day.is_dir() else [l1b]
    for root in roots:
        for f in sorted(root.glob(f"PACE_OCI.{d}T*.L1B.V3.nc")):
            # PACE_OCI.YYYYMMDDTHHMMSS.L1B.V3.nc
            stem = f.name
            m = re.match(r"PACE_OCI\.(.+)\.L1B\.V3\.nc$", stem)
            if m:
                gids.append(m.group(1))
if not gids:
    sys.stderr.write("ERROR: no L1B granules found for resolved dates under L1B_dir\n")
    sys.exit(2)
print(f"#META n_granules={len(gids)}", file=sys.stderr)
for g in gids:
    print(g)
PY
)

# Fallback if python helper failed entirely
if [[ ${#TASKS[@]} -eq 0 ]]; then
  echo "ERROR: failed to build task list from ${CONFIG}" >&2
  echo "Ensure Python 3 is available and [svd_retrieval] has dates or date_start/date_end." >&2
  exit 2
fi

# Filter META lines if any leaked to stdout (shouldn't)
TASK_LIST=()
for t in "${TASKS[@]}"; do
  [[ "${t}" == \#* ]] && continue
  [[ -z "${t}" ]] && continue
  TASK_LIST+=("${t}")
done
TASKS=("${TASK_LIST[@]}")
N_TASKS=${#TASKS[@]}

# Log directory
if [[ -z "${LOG_DIR:-}" ]]; then
  OUT_DIR="$(grep -E '^\s*output_dir\s*=' "${CONFIG}" | head -1 | sed -E 's/.*=\s*"([^"]+)".*/\1/' || true)"
  if [[ -n "${OUT_DIR}" ]]; then
    LOG_DIR="${OUT_DIR}/multitask_logs"
  else
    LOG_DIR="${REPO_ROOT}/svd_retrieval_multitask_logs"
  fi
fi
if ! mkdir -p "${LOG_DIR}" 2>/dev/null; then
  LOG_DIR="${REPO_ROOT}/svd_retrieval_multitask_logs"
  mkdir -p "${LOG_DIR}"
  echo "WARN: falling back to LOG_DIR=${LOG_DIR}" >&2
fi
STAMP="$(date +%Y%m%dT%H%M%S)"
MASTER_LOG="${LOG_DIR}/master_${STAMP}.log"

echo "=============================================="
echo "SVD multitask launcher"
echo "  config      : ${CONFIG}"
echo "  split       : ${SPLIT}   (n_tasks=${N_TASKS})"
echo "  max_jobs    : ${MAX_JOBS}"
echo "  julia -t    : ${JULIA_THREADS}"
echo "  log_dir     : ${LOG_DIR}"
echo "  dry_run     : ${DRY_RUN}"
echo "==============================================" | tee "${MASTER_LOG}"

run_one() {
  local task="$1"
  local tag log cmd
  if [[ "${SPLIT}" == "granule" ]]; then
    tag="granule_${task}"
    cmd=( "${JULIA_BIN}" --project="${JULIA_PROJECT}" -t "${JULIA_THREADS}"
          "${RUN_JL}" "${CONFIG}" --granule "${task}" )
  else
    tag="date_${task}"
    cmd=( "${JULIA_BIN}" --project="${JULIA_PROJECT}" -t "${JULIA_THREADS}"
          "${RUN_JL}" "${CONFIG}" --date "${task}" )
  fi
  log="${LOG_DIR}/${tag}_${STAMP}.log"
  echo "[$(date +%H:%M:%S)] START ${tag}" | tee -a "${MASTER_LOG}"
  if [[ "${DRY_RUN}" == "1" ]]; then
    printf '  DRY:'; printf ' %q' "${cmd[@]}"; echo "  > ${log}"
    return 0
  fi
  if "${cmd[@]}" >"${log}" 2>&1; then
    echo "[$(date +%H:%M:%S)] DONE  ${tag}" | tee -a "${MASTER_LOG}"
    return 0
  else
    local rc=$?
    echo "[$(date +%H:%M:%S)] FAIL  ${tag} (rc=${rc})  see ${log}" | tee -a "${MASTER_LOG}"
    return "${rc}"
  fi
}

# Simple job pool
pids=()
fails=0
for task in "${TASKS[@]}"; do
  # wait if at capacity
  while (( ${#pids[@]} >= MAX_JOBS )); do
    new_pids=()
    for pid in "${pids[@]}"; do
      if kill -0 "${pid}" 2>/dev/null; then
        new_pids+=("${pid}")
      else
        if ! wait "${pid}"; then
          fails=$((fails + 1))
        fi
      fi
    done
    pids=("${new_pids[@]}")
    (( ${#pids[@]} >= MAX_JOBS )) && sleep 2
  done
  run_one "${task}" &
  pids+=($!)
done

for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    fails=$((fails + 1))
  fi
done

echo "==============================================" | tee -a "${MASTER_LOG}"
echo "Finished: ${N_TASKS} tasks, failures=${fails}" | tee -a "${MASTER_LOG}"
echo "Master log: ${MASTER_LOG}" | tee -a "${MASTER_LOG}"
exit "${fails}"
