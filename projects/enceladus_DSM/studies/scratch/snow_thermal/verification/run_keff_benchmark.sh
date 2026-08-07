#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_keff_benchmark.sh — the analytic gate on the k_eff cell-problem solver.
#
# A 50/50 ice/air layered medium (-ic_type ice_slab) has closed-form effective
# conductivity, so it decides whether the cell problem is solving the right
# equations. Two sweeps, because two different things are being tested:
#
#   S  SOLVER convergence: fix eps, refine the mesh. Both tensor components must
#      converge to the EXACT value for the diffuse profile actually present
#      (computed by 1-D quadrature in check_keff_benchmark.py, no PDE needed --
#      for a layered medium k_par is the mean of k and k_perp is the harmonic
#      mean of k). This isolates the solver.
#
#   E  PHYSICS: shrink eps at adequate mesh. k_perp then approaches the
#      SHARP-interface harmonic mean from above, because the diffuse band is a
#      high-conductivity short circuit across an otherwise series resistance.
#      This is a property of the model, not an error, and it calibrates how to
#      read the safety 0.5 vs 1.0 comparison in the packing runs.
#
# Direct LU throughout: this gate exists to validate the assembly and the
# tensor, so the linear solver must not be a variable. Step 5 revalidates
# CG+GAMG against these same numbers.
#
# Usage:  ./run_keff_benchmark.sh [output_dir]
# ---------------------------------------------------------------------------
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
PROJECT_ROOT="$PWD"
EXEC="$PROJECT_ROOT/enceladus_dsm"

OUT="${1:-$PROJECT_ROOT/studies/snow_thermal/verification}"
RAW="$OUT/raw_keff"
CSV="$OUT/keff_benchmark.csv"
mkdir -p "$RAW"

: "${PETSC_DIR:=/Users/jacksonbaglino/Petsc-3.20}"
: "${PETSC_ARCH:=MAC_Si-M2_MPICH-4.3.0_PETSC-3.20_GCC-12.4.0}"
: "${PETIGA_DIR:=/Users/jacksonbaglino/PetIGA-3.20}"
export PETSC_DIR PETSC_ARCH PETIGA_DIR

[[ -x "$EXEC" ]] || { echo "ERROR: $EXEC missing; run make." >&2; exit 1; }

LX=1.0e-3

# sweep : eps : Nx
CONFIGS=(
    # --- S: fixed eps, refining mesh -> isolates the solver ---
    "S_mesh:2.0e-5:64"
    "S_mesh:2.0e-5:128"
    "S_mesh:2.0e-5:256"
    "S_mesh:2.0e-5:512"
    # --- E: shrinking eps at adequate mesh -> the sharp-interface limit ---
    "E_eps:3.125e-5:256"
    "E_eps:1.5625e-5:256"
    "E_eps:7.8125e-6:512"
)

echo "sweep,eps,Nx,phi_bar,k00,k01,k10,k11,wall_s" > "$CSV"

for cfg in "${CONFIGS[@]}"; do
    IFS=: read -r sweep eps nx <<< "$cfg"
    tag="${sweep}_eps${eps}_N${nx}"
    log="$RAW/${tag}.log"
    folder="$RAW/$tag"; mkdir -p "$folder"; export folder

    "$EXEC" \
        -options_file "$PROJECT_ROOT/inputs/solver.opts" \
        -ic_type ice_slab -dim 2 \
        -Lx "$LX" -Ly "$LX" -Nx "$nx" -Ny "$nx" \
        -p 2 -C 1 -periodic 1 \
        -eps "$eps" -eps_temp_override 1 \
        -temp -20 -humidity 1 -t_final 1.0 -delt_t 1.0e-4 \
        -pf_output 0 -pf_monitor 0 \
        -keff 1 -keff_only 1 \
        -keff_ksp_type preonly -keff_pc_type lu \
        > "$log" 2>&1
    rc=$?
    if [[ $rc -ne 0 ]]; then echo "  FAILED (rc=$rc): $tag -- see $log" >&2; continue; fi

    phi=$(awk '/\[keff\] step/ {for(i=1;i<=NF;i++) if($i=="phi_bar") {print $(i+2); exit}}' "$log")
    wall=$(awk '/\[keff\] step/ {gsub(/[(),]/," "); for(i=1;i<=NF;i++) if($i=="s)"||$i=="s") ; print $(NF-1); exit}' "$log")
    read -r k00 k01 <<< "$(awk '/k = \[/ {print $4, $5; exit}' "$log")"
    read -r k10 k11 <<< "$(awk '/^ *\[ / && !/k = / {print $2, $3; exit}' "$log")"

    echo "$sweep,$eps,$nx,$phi,$k00,$k01,$k10,$k11,$wall" >> "$CSV"
    echo "  ok  $tag   phi_bar=$phi  k00=$k00  k11=$k11"
done

echo
echo "wrote $CSV"
echo "now run:  venv_enceladus/bin/python3 studies/snow_thermal/verification/check_keff_benchmark.py"
