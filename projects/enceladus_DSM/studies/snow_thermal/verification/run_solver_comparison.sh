#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_solver_comparison.sh — CG+GAMG against direct LU on an IDENTICAL assembly.
#
# The standalone homogenization project asserted in its README that iterative
# solvers are unreliable for this problem class, while its own solver.h
# documented CG+GAMG and its code shipped PREONLY+LU. That claim is the only
# thing standing between this study and a solver that can actually run at the
# study meshes (~8M unknowns per corrector, where a direct factorization per
# sample is not viable). So it gets tested rather than inherited.
#
# Because both runs use the same binary, the same initial condition and the same
# assembly, any difference in the resulting tensor is attributable to the linear
# solver alone.
#
# Usage:  ./run_solver_comparison.sh [output_dir]
# ---------------------------------------------------------------------------
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
PROJECT_ROOT="$PWD"
EXEC="$PROJECT_ROOT/enceladus_dsm"

OUT="${1:-$PROJECT_ROOT/studies/snow_thermal/verification}"
RAW="$OUT/raw_solver"
CSV="$OUT/solver_comparison.csv"
mkdir -p "$RAW"

: "${PETSC_DIR:=/Users/jacksonbaglino/Petsc-3.20}"
: "${PETSC_ARCH:=MAC_Si-M2_MPICH-4.3.0_PETSC-3.20_GCC-12.4.0}"
: "${PETIGA_DIR:=/Users/jacksonbaglino/PetIGA-3.20}"
export PETSC_DIR PETSC_ARCH PETIGA_DIR

[[ -x "$EXEC" ]] || { echo "ERROR: $EXEC missing; run make." >&2; exit 1; }

LX=1.0e-3
EPS=2.0e-5

# Meshes are multiples of the PROJECT RULE, Nx = ceil(sqrt(2)*Lx/eps), which is
# comp_eps.py's h = eps/sqrt(2) and gives ~8.5 elements across the phi=0.05-0.95
# band (comp_eps reports this as "~7.5 elements"). At eps=2e-5 on a 1 mm cell
# that rule is Nx=71, so 71/142/284/568 are 1x/2x/4x/8x the production
# resolution. Reporting the 5-95% band count, NOT the 1-99% count -- the two
# differ by a factor 9.2/6 = 1.53 and mixing them is how a production-resolution
# mesh gets mistaken for a heavily over-resolved one.
NBASE=71

echo "solver,Nx,mult,elem_per_band_5_95,ndof,phi_bar,k00,k01,k10,k11,its,wall_s" > "$CSV"

run_one () {                       # $1 = label, $2 = Nx, $3.. = solver flags
    local label="$1" nx="$2"; shift 2
    local tag="${label}_N${nx}"
    local log="$RAW/${tag}.log"
    folder="$RAW/$tag"; mkdir -p "$folder"; export folder

    /usr/bin/time -p "$EXEC" \
        -options_file "$PROJECT_ROOT/inputs/solver.opts" \
        -ic_type ice_slab -dim 2 \
        -Lx "$LX" -Ly "$LX" -Nx "$nx" -Ny "$nx" \
        -p 2 -C 1 -periodic 1 \
        -eps "$EPS" -eps_temp_override 1 \
        -temp -20 -humidity 1 -t_final 1.0 -delt_t 1.0e-4 \
        -pf_output 0 -pf_monitor 0 \
        -keff 1 -keff_only 1 "$@" \
        > "$log" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then echo "  FAILED (rc=$rc): $tag -- see $log" >&2; return; fi

    local phi its wall k00 k01 k10 k11
    phi=$(awk '/\[keff\] step/ {for(i=1;i<=NF;i++) if($i=="phi_bar"){print $(i+2); exit}}' "$log")
    its=$(awk '/\[keff\] step/ {gsub(/\(/," "); for(i=1;i<=NF;i++) if($i=="its,"){print $(i-1); exit}}' "$log")
    wall=$(awk '/^real/ {print $2; exit}' "$log")
    read -r k00 k01 <<< "$(awk '/k = \[/ {print $4, $5; exit}' "$log")"
    read -r k10 k11 <<< "$(awk '/^ *\[ / && !/k = / {print $2, $3; exit}' "$log")"

    local mult epb
    mult=$(awk -v n="$nx" -v b="$NBASE" 'BEGIN{printf "%.0f", n/b}')
    epb=$(awk -v e="$EPS" -v l="$LX" -v n="$nx" 'BEGIN{printf "%.1f", 6*e/(l/n)}')

    echo "$label,$nx,$mult,$epb,$((nx*nx)),$phi,$k00,$k01,$k10,$k11,$its,$wall" >> "$CSV"
    printf '  %-4s N=%-5s (%2sx rule, %5s elem/band)  k11=%s  its=%-5s %ss\n' \
        "$label" "$nx" "$mult" "$epb" "$k11" "$its" "$wall"
}

# LU and GAMG head to head, from the production resolution up to where a direct
# factorization stops being comfortable.
for nx in 71 142 284 568; do
    run_one lu   "$nx" -keff_ksp_type preonly -keff_pc_type lu
    run_one gamg "$nx"
done

# Iterative only, past LU's reach -- the regime the study meshes actually live in.
for nx in 1136 2272; do
    run_one gamg "$nx"
done

echo
echo "wrote $CSV"
