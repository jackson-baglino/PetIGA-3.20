#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_ic_resolution_study.sh — how well the layered benchmark IC is resolved,
# as a function of the interface width eps and the mesh.
#
# WHY THIS EXISTS. The 50/50 ice/air slab (-ic_type ice_slab) is the analytic
# benchmark for the k_eff homogenization: a layered two-phase medium has
# closed-form effective conductivity. But the analytic values are properties of
# the CONTINUOUS profile, while the solver sees a spline fit to nodal values.
# Before using the benchmark to judge the cell-problem solver, we need to know
# how much of any discrepancy the initial condition is responsible for on its
# own -- otherwise the IC's quadrature error gets reported as solver error.
#
# Three sweeps:
#   A  mesh follows eps, using the project rule Nx = ceil(sqrt(2)*Lx/eps).
#      This is the configuration the real study runs use, at three interface
#      widths corresponding to comp_eps.py safety factors 1.0 / 0.5 / 0.25.
#   B  eps varies at a FIXED fine mesh, isolating the effect of eps alone.
#   C  mesh refines at FIXED eps, isolating the effect of the mesh alone.
#
# Every run is a single time step on a 1 mm cell and takes well under a second.
#
# Usage:  ./run_ic_resolution_study.sh [output_dir]
# ---------------------------------------------------------------------------
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."   # -> project root
PROJECT_ROOT="$PWD"
EXEC="$PROJECT_ROOT/enceladus_dsm"

OUT="${1:-$PROJECT_ROOT/studies/snow_thermal/verification}"
RAW="$OUT/raw"
CSV="$OUT/ic_resolution_study.csv"
mkdir -p "$RAW"

: "${PETSC_DIR:=/Users/jacksonbaglino/Petsc-3.20}"
: "${PETSC_ARCH:=MAC_Si-M2_MPICH-4.3.0_PETSC-3.20_GCC-12.4.0}"
: "${PETIGA_DIR:=/Users/jacksonbaglino/PetIGA-3.20}"
export PETSC_DIR PETSC_ARCH PETIGA_DIR

if [[ ! -x "$EXEC" ]]; then
    echo "ERROR: $EXEC not found. Run 'make' in $PROJECT_ROOT first." >&2
    exit 1
fi

LX=1.0e-3

# All three sweeps are anchored on the PROJECT MESH RULE, Nx = ceil(sqrt(2)*Lx/eps)
# (comp_eps.py's h = eps/sqrt(2)), which gives ~8.5 elements across the
# phi=0.05-0.95 band. Every sweep passes through that production point rather
# than sitting to one side of it, so the sweeps bracket the resolution the real
# runs use instead of only sampling over-resolved meshes.
#
#   eps      rule Nx    (= ceil(sqrt(2) * 1e-3 / eps))
#   4.0e-5   36
#   2.0e-5   71
#   1.0e-5   142
#   5.0e-6   283
#
# sweep_label : eps : Nx
CONFIGS=(
    # --- A: mesh follows eps -- every point AT the production resolution ---
    "A_follow:4.0e-5:36"
    "A_follow:2.0e-5:71"
    "A_follow:1.0e-5:142"
    "A_follow:5.0e-6:283"
    # --- B: eps varies at fixed mesh; the rule's own Nx for eps=1e-5 ---
    "B_fixed_mesh:4.0e-5:142"
    "B_fixed_mesh:2.0e-5:142"
    "B_fixed_mesh:1.0e-5:142"
    "B_fixed_mesh:5.0e-6:142"
    # --- C: mesh refines at fixed eps, spanning 0.5x .. 4x the rule ---
    "C_fixed_eps:2.0e-5:36"
    "C_fixed_eps:2.0e-5:71"
    "C_fixed_eps:2.0e-5:142"
    "C_fixed_eps:2.0e-5:284"
)

echo "sweep,eps,Nx,dx,elem_per_band_5_95,tot_ice,tot_air,interf,phi_bar_solver,phi_bar_clone,proj_worst,phi_min,phi_max" > "$CSV"

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
        -ts_max_steps 1 -pf_output 0 \
        -keff 1 -keff_only 0 -keff_debug_phibar 1 \
        > "$log" 2>&1
    rc=$?

    if [[ $rc -ne 0 ]]; then
        echo "  FAILED (rc=$rc): $tag  -- see $log" >&2
        continue
    fi

    # Step-0 row of the monitor table: STEP TIME DT TOT_ICE TOT_AIR TEMP TOT_RHOV I-A_INTERF TOTAL_MASS
    read -r tot_ice tot_air interf <<< "$(
        awk '/^ *0 \|/ {gsub(/\|/," "); print $4, $5, $8; exit}' "$log")"

    read -r pb_clone pb_solver <<< "$(
        awk '/mean ice fraction/ {print $4, $5; exit}' "$log")"

    worst=$(awk '/worst =/ {print $3; exit}' "$log")

    read -r pmin pmax <<< "$(
        awk '/BOUNDS: phi_ice/ {gsub(/[\[\],]/," "); print $3, $4; exit}' "$log")"

    dx=$(awk -v l="$LX" -v n="$nx" 'BEGIN{printf "%.6e", l/n}')
    # Elements spanning the phi = 0.05 .. 0.95 band (6*eps). THIS is the project
    # convention: comp_eps.py's mesh rule h = eps/sqrt(2) targets ~7.5-8.5 across
    # it, and that is the number to compare a mesh against. The 1-99% band
    # (9.2*eps) gives a count 1.53x larger for the same mesh; quoting that one
    # makes a production-resolution mesh look heavily over-resolved.
    epb=$(awk -v e="$eps" -v l="$LX" -v n="$nx" 'BEGIN{printf "%.3f", 6.0*e/(l/n)}')

    echo "$sweep,$eps,$nx,$dx,$epb,$tot_ice,$tot_air,$interf,$pb_solver,$pb_clone,$worst,$pmin,$pmax" >> "$CSV"
    echo "  ok  $tag   phi_bar=$pb_solver  interf=$interf  elem/band=$epb"
done

echo
echo "wrote $CSV"
