# Launch generate_test_spec_lambertian.jl (1000 samples) on the first GPU with ≥12 GiB free.
# Same seed as generate_test_spec.sh → same profiles / aerosol columns / SIF / geometry
# as output/rt_toa_ensemble.nc, with a Lambertian surface (albedo ~ U(0.005, 0.025)).
DIR=/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/full_RT_construction
OUTDIR=$DIR/output_new_Lambertian
mkdir -p "$OUTDIR"
STAMP=$(date +%Y%m%d_%H%M)
LOG=$OUTDIR/generate_test_spec_lambertian_n1000_${STAMP}.log
OUT_NC=$OUTDIR/rt_toa_ensemble_lambertian.nc
NEED_MIB=12288   # ≥12 GiB free

nohup bash -c "
set -euo pipefail
echo \"[\$(date -Iseconds)] wait for GPU ≥${NEED_MIB} MiB free\" | tee -a '$LOG'
while true; do
  mapfile -t freelist < <(nvidia-smi --query-gpu=index,memory.free --format=csv,noheader,nounits \
    | tr -d ' ' | sort -t, -k2 -nr)
  best=\${freelist[0]}
  idx=\${best%%,*}
  free=\${best##*,}
  if (( free >= $NEED_MIB )); then
    echo \"[\$(date -Iseconds)] using GPU \$idx (free MiB from: \${freelist[*]})\" | tee -a '$LOG'
    export CUDA_VISIBLE_DEVICES=\$idx
    export ARCH=GPU
    export N_SAMPLES=1000
    export OUT_NC='$OUT_NC'
    cd '$DIR'
    exec julia --project=/home/zhe2/FraLab/vSmartMOM.jl generate_test_spec_lambertian.jl
  fi
  echo \"[\$(date -Iseconds)] GPUs busy (\${freelist[*]}); sleeping 5m\" | tee -a '$LOG'
  sleep 300
done
" >>"$LOG" 2>&1 &

echo "PID=$!  LOG=$LOG  OUT_NC=$OUT_NC"
