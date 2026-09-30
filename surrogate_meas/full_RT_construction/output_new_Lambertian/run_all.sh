#!/bin/bash
# 4 runs: {Lambertian 0.02, Cox-Munk 5 m/s} × {no aerosol, GCHP column 100}; same script, same atmosphere/SIF/geometry.
cd /home/zhe2/FraLab/PACE_redSIF_PACE/.claude/worktrees/individual-granule-notebook-ebca58
S=surrogate_meas/full_RT_construction/test_run_lambertian.jl
for Y in ocean_lambertian_0929.yaml ocean_coxmunk_0912.yaml; do
  for A in false true; do
    echo "=== $Y aerosols=$A"
    OUT_DIR=/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/full_RT_construction/output_new_Lambertian OCEAN_YAML=/home/zhe2/FraLab/PACE_redSIF_PACE/surrogate_meas/configs/$Y ENABLE_AEROSOLS=$A COLUMN_INDEX=100 julia -t 8 $S
  done
done
echo ALL_DONE
