#!/bin/bash
# Run generate_rayleigh_scaling_surface.jl (TOA + near-surface Ed, Eu, Lw, L_down, Rrs) at one or more SZAs.
#
# For each SZA that differs from the base YAML, a copy of the base YAML is written next to it with only
# `geometry: sza:` changed, named <base>_sza<SZA with . -> p>.yaml (configs/ is git-ignored, so the copy is
# local). The output NetCDF name follows from the YAML name, e.g.
#   rt_toa_rayleigh_scaling_ocean_lambertian_1010_wavnum0p05_NoRayl_sza60p5_sunlight_on_AerColumn_100_grid_surface.nc
# Existing outputs are never overwritten (the Julia script refuses).
#
# SZA 60° is avoided on purpose: with an odd number of streams (11 in the base YAML) there is a Gauss node
# at mu = 0.5, and the direct beam is double-counted in hdr / bhr_dw / bhr_uw (Lw, Ed, Eu +13.6%); TOA radiance
# is unaffected. Use 60.5° (or an even number of streams).
#
# Usage (all settings can be overridden from the environment):
#   bash surrogate_meas/full_RT_construction/generate_rayleigh_scaling_surface.sh
#   SZAS="30 45 60.5" AER_COLUMN=77 bash .../generate_rayleigh_scaling_surface.sh
#   DRY_RUN=1 bash .../generate_rayleigh_scaling_surface.sh      # only write the YAMLs and print the commands
set -euo pipefail

REPO=${REPO:-/home/zhe2/FraLab/PACE_redSIF_PACE}
RT_DIR=$REPO/surrogate_meas/full_RT_construction
BASE_YAML=${BASE_YAML:-$REPO/surrogate_meas/configs/ocean_lambertian_1010_wavnum0p05_NoRayl.yaml}
SZAS=${SZAS:-"30 60.5"}
JULIA_PROJECT_DIR=${JULIA_PROJECT_DIR:-$RT_DIR/gchpio_env}
DRY_RUN=${DRY_RUN:-0}

export AER_COLUMN=${AER_COLUMN:-100}        # 0 = no aerosol
export SUNLIGHT_ON=${SUNLIGHT_ON:-true}      # false = solar beam off
export SIF_INDEX=${SIF_INDEX:-1}
export TMPDIR=${TMPDIR:-$REPO/.claude/scratch/tmp}
mkdir -p "$TMPDIR"

base_sza=$(awk '/^geometry:/{g=1; next} g && /^[^ #]/{g=0} g && $1=="sza:"{print $2; exit}' "$BASE_YAML")
[ -n "$base_sza" ] || { echo "no geometry: sza: in $BASE_YAML" >&2; exit 1; }

for sza in $SZAS; do
    if awk -v a="$sza" 'BEGIN{exit !(a+0 == 60)}'; then
        echo "WARNING: SZA $sza sits on a Gauss node with odd nstreams; surface fields will be wrong (see header)." >&2
    fi
    if awk -v a="$sza" -v b="$base_sza" 'BEGIN{exit !(a+0 == b+0)}'; then
        yaml=$BASE_YAML
    else
        yaml=${BASE_YAML%.yaml}_sza${sza/./p}.yaml
        # change only the sza line inside the geometry: block
        awk -v s="$sza" '/^geometry:/{g=1} g && /^[^ #]/ && !/^geometry:/{g=0}
                         g && $1=="sza:" {sub(/sza:[ \t]+[0-9.]+/, "sza:     " s)} {print}' "$BASE_YAML" > "$yaml"
        echo "wrote $yaml"
    fi
    cmd=(julia --project="$JULIA_PROJECT_DIR" "$RT_DIR/generate_rayleigh_scaling_surface.jl")
    echo "SZA $sza: SCAN_YAML=$yaml ${cmd[*]}"
    [ "$DRY_RUN" = 1 ] || SCAN_YAML=$yaml "${cmd[@]}"
done
