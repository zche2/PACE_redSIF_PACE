#!/usr/bin/env bash
# Launch generate_test_spec_lambertian_gchpio.jl in the background on the first GPU with enough
# free memory. Every setting below can be edited here or overridden from the environment, e.g.
#   N_SAMPLES=200 OUTDIR=/path/out OCEAN_YAML=/path/cfg.yaml bash generate_test_spec_lambertian_gchpio.sh
#   GPU_ID=1 N_RUN=3 OUT_NC=/tmp/test.nc bash generate_test_spec_lambertian_gchpio.sh
# Same ENSEMBLE_SEED and N_SAMPLES as an existing ensemble → same profiles / aerosol columns /
# SIF / geometry / albedos (sample k matches sample k).
set -euo pipefail

DIR=/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/full_RT_construction
DATA=/home/zhe2/data/MyProjects/PACE_redSIF_PACE/Files_in_use

# ── ensemble ────────────────────────────────────────────────────────────────
: "${N_SAMPLES:=1000}"            # design size (changing it changes all draws)
: "${N_RUN:=0}"                   # compute only the next N_RUN samples this run (0 = all)
: "${ENSEMBLE_SEED:=20260913}"    # design seed (albedos use seed + 1)
: "${ENABLE_AEROSOLS:=true}"

# ── inputs ──────────────────────────────────────────────────────────────────
: "${CONFIG_DIR:=$DIR/../configs}"
: "${OCEAN_YAML:=$CONFIG_DIR/ocean_lambertian_1010_wavnum0p05.yaml}"            # spectral grid, streams, surface
: "${MERRA_SEA_TPQ_NC:=$DIR/merra2_sea_Tpq_columns_n1500.nc}"         # T/p/q profiles
: "${OCEAN_COLS_NC:=$DIR/output_aerosol_profiles/gchp_ocean_columns_n500.nc}"  # GCHP cell index
: "${GCHP_NC:=}"                  # GCHP file; empty = the `source` attribute of OCEAN_COLS_NC
: "${PACE_RSR_NC:=$DATA/PACE_OCI_RSRs.nc}"
: "${SIF_LIB:=$DATA/SIF_singular_vector.jld2}"
: "${SNR_FILE:=$DATA/PACE_OCI_L1BLUT_baseline_SNR_1.1.txt}"

# ── outputs ─────────────────────────────────────────────────────────────────
: "${OUTDIR:=$DIR/output_new_Lambertian_gchpio_finerwavnum0p05}"
: "${OUT_NC:=$OUTDIR/rt_toa_ensemble_lambertian_gchpio_n1000.nc}"          # existing file → resume
: "${OPTICS_CACHE:=${OUT_NC%.nc}_aerosol_optics.jld2}"               # per-column GCHPIO optics
: "${LOG:=$OUTDIR/generate_test_spec_lambertian_gchpio_n${N_SAMPLES}_$(date +%Y%m%d_%H%M).log}"

# ── aerosol optics / RT ─────────────────────────────────────────────────────
: "${AER_TRUNCATION:=dbge}"       # dbge (δBGE; within 0.005 of 24 streams) | yaml | none
: "${N_NODES:=3}"                 # Mie λ-nodes (GCHPIO auto-raises for wide bands: 7 for 623–850 nm)
: "${N_QUERY:=24}"                # λ points where τ is evaluated across the RT grid
: "${MIE_ARCH:=CPU}"              # CPU (threads) | GPU
: "${ARCH:=GPU}"                  # RT architecture: GPU | CPU | DEFAULT

# ── runtime ─────────────────────────────────────────────────────────────────
: "${JULIA_ENV:=$DIR/gchpio_env}"
: "${THREADS:=16}"                # julia -t (Mie threads)
: "${GPU_ID:=}"                   # fixed GPU index; empty = wait for one with ≥ NEED_MIB free
: "${NEED_MIB:=12288}"

mkdir -p "$OUTDIR" "$(dirname "$OUT_NC")" "$(dirname "$LOG")"
for f in "$OCEAN_YAML" "$MERRA_SEA_TPQ_NC" "$PACE_RSR_NC" "$SIF_LIB" "$SNR_FILE" "$JULIA_ENV/Project.toml"; do
  [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done
[[ "$ENABLE_AEROSOLS" != true || -f "$OCEAN_COLS_NC" ]] || { echo "missing: $OCEAN_COLS_NC" >&2; exit 1; }

export N_SAMPLES N_RUN ENSEMBLE_SEED ENABLE_AEROSOLS OCEAN_YAML MERRA_SEA_TPQ_NC OCEAN_COLS_NC \
       GCHP_NC PACE_RSR_NC SIF_LIB SNR_FILE OUT_NC OPTICS_CACHE AER_TRUNCATION N_NODES N_QUERY \
       MIE_ARCH ARCH JULIA_ENV THREADS GPU_ID NEED_MIB DIR LOG

{
  echo "[$(date -Iseconds)] settings:"
  for v in N_SAMPLES N_RUN ENSEMBLE_SEED ENABLE_AEROSOLS OCEAN_YAML MERRA_SEA_TPQ_NC OCEAN_COLS_NC \
           GCHP_NC OUT_NC OPTICS_CACHE AER_TRUNCATION N_NODES N_QUERY MIE_ARCH ARCH THREADS GPU_ID; do
    printf '  %-17s %s\n' "$v" "${!v}"
  done
} | tee -a "$LOG"

nohup bash -c '
set -euo pipefail
if [[ -n "$GPU_ID" ]]; then
  idx=$GPU_ID
else
  echo "[$(date -Iseconds)] wait for GPU ≥${NEED_MIB} MiB free"
  while true; do
    mapfile -t freelist < <(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits \
      | tr -d " " | sort -t, -k2 -nr)
    best=${freelist[0]}; idx=${best%%,*}; free=${best##*,}
    (( free >= NEED_MIB )) && break
    echo "[$(date -Iseconds)] GPUs busy (${freelist[*]}); sleeping 5m"
    sleep 300
  done
fi
echo "[$(date -Iseconds)] using GPU $idx"
export CUDA_VISIBLE_DEVICES=$idx
cd "$DIR"
exec julia --project="$JULIA_ENV" -t "$THREADS" generate_test_spec_lambertian_gchpio.jl
' >>"$LOG" 2>&1 &

echo "PID=$!  LOG=$LOG  OUT_NC=$OUT_NC"
