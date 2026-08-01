#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run_t0_probe.sh — isolate the eps ARTEFACT in k_eff from the eps PHYSICS.
#
# Comparing k_eff between safety 0.5 and safety 1.0 at the end of a run mixes at
# least three things:
#   1. throat closure   -- the diffuse band is wider than a pore throat, so phi
#                          never reaches 0 there and k_eff reads high;
#   2. sharp-interface-limit error in the kinetics -- which is what the K&P
#                          safety factor is actually there to bound;
#   3. a genuinely different microstructure at t = t_final, because the two eps
#                          sinter at different rates.
#
# At t = 0 only (1) can be present. The geometries use -ic_grain_union 1, whose
# whole point is that the phi = 0.5 contour sits on the SHARP union surface for
# any eps -- so both resolutions describe the identical microstructure and differ
# only in how thick the diffuse band around it is. No time integration is needed
# at all (-keff_only 1), so this costs one cell solve per configuration.
#
# Both resolutions read the SAME grains.dat with the same -Lx and the same union
# IC; only -eps and -Nx differ. The comparison is therefore paired within a
# configuration, and packing-to-packing realization scatter cancels exactly --
# each packing is only ever compared with itself.
#
# READ THE RESULT AS AN UPPER BOUND, NOT AS THE STUDY'S ERROR BAR. t = 0 is where
# this artefact is largest, by construction: the packings are jammed with a 4 um
# design contact gap, the solid barely percolates
# (solid_largest_cluster_frac ~ 0.009), so every conduction path runs through
# those 4 um gaps -- and the diffuse band is 9.2 um (safety 0.5) or 18.4 um
# (safety 1.0), wider than the gap either way. At t = 0 the band IS the
# conduction path. Once DSM builds necks the path is solid ice and the
# sensitivity should fall. The number that matters for the study is the gap at
# the END of a run; see the README for how to measure it.
#
# Usage:  ./run_t0_probe.sh [safety05|safety10|both]
# ---------------------------------------------------------------------------
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."
PROJECT_ROOT="$PWD"
EXEC="$PROJECT_ROOT/enceladus_dsm"

WHICH="${1:-both}"
OUT="$PROJECT_ROOT/studies/snow_thermal/verification"
RAW="$OUT/raw_t0"
CSV="$OUT/t0_probe.csv"
mkdir -p "$RAW"

: "${PETSC_DIR:=/Users/jacksonbaglino/Petsc-3.20}"
: "${PETSC_ARCH:=MAC_Si-M2_MPICH-4.3.0_PETSC-3.20_GCC-12.4.0}"
: "${PETIGA_DIR:=/Users/jacksonbaglino/PetIGA-3.20}"
: "${NPROCS:=1}"
export PETSC_DIR PETSC_ARCH PETIGA_DIR

[[ -x "$EXEC" ]] || { echo "ERROR: $EXEC missing; run make." >&2; exit 1; }

PACKINGS=(phi0.250_seed1 phi0.325_seed3 phi0.400_seed5)

[[ -f "$CSV" ]] || echo "packing,safety,eps,Nx,phi_bar,k00,k01,k10,k11,k_iso,its,wall_s" > "$CSV"

probe () {                                   # $1 = packing, $2 = safety label
    local pack="$1" slabel="$2" geom
    case "$slabel" in
        0.5) geom="inputs/geometry/${pack}_T-20.opts" ;;
        1.0) geom="inputs/geometry/${pack}_T-20_s10.opts" ;;
        *)   echo "bad safety label $slabel" >&2; return 1 ;;
    esac
    [[ -f "$geom" ]] || { echo "  MISSING $geom" >&2; return 1; }

    local tag="${pack}_s${slabel}"
    local log="$RAW/${tag}.log"
    folder="$RAW/$tag"; mkdir -p "$folder"; export folder

    local eps nx
    eps=$(awk '/^-eps /{print $2; exit}' "$geom")
    nx=$(awk '/^-Nx /{print $2; exit}' "$geom")
    echo "  running $tag  (eps=$eps, Nx=$nx) ..."

    local t_start t_end
    t_start=$(date +%s)
    if [[ "$NPROCS" -gt 1 ]]; then
        "$PETSC_DIR/$PETSC_ARCH/bin/mpiexec" -np "$NPROCS" "$EXEC" \
            -options_file "$PROJECT_ROOT/inputs/solver.opts" \
            -options_file "$PROJECT_ROOT/$geom" \
            -options_file "$PROJECT_ROOT/inputs/experiment/snow_T-20_1day.opts" \
            -pf_output 0 -pf_monitor 0 -keff 1 -keff_only 1 > "$log" 2>&1
    else
        "$EXEC" \
            -options_file "$PROJECT_ROOT/inputs/solver.opts" \
            -options_file "$PROJECT_ROOT/$geom" \
            -options_file "$PROJECT_ROOT/inputs/experiment/snow_T-20_1day.opts" \
            -pf_output 0 -pf_monitor 0 -keff 1 -keff_only 1 > "$log" 2>&1
    fi
    local rc=$? ; t_end=$(date +%s)
    if [[ $rc -ne 0 ]]; then echo "    FAILED (rc=$rc) -- see $log" >&2; return 1; fi

    local phi its k00 k01 k10 k11 kiso
    phi=$(awk '/\[keff\] step/ {for(i=1;i<=NF;i++) if($i=="phi_bar"){print $(i+2); exit}}' "$log")
    kiso=$(awk '/\[keff\] step/ {for(i=1;i<=NF;i++) if($i=="k_iso"){print $(i+2); exit}}' "$log")
    its=$(awk '/\[keff\] step/ {gsub(/\(/," "); for(i=1;i<=NF;i++) if($i=="its,"){print $(i-1); exit}}' "$log")
    read -r k00 k01 <<< "$(awk '/k = \[/ {print $4, $5; exit}' "$log")"
    read -r k10 k11 <<< "$(awk '/^ *\[ / && !/k = / {print $2, $3; exit}' "$log")"

    echo "$pack,$slabel,$eps,$nx,$phi,$k00,$k01,$k10,$k11,$kiso,$its,$((t_end-t_start))" >> "$CSV"
    echo "    k_iso=$kiso  phi_bar=$phi  ($its its, $((t_end-t_start)) s)"
}

for p in "${PACKINGS[@]}"; do
    [[ "$WHICH" == "safety10" || "$WHICH" == "both" ]] && probe "$p" 1.0
    [[ "$WHICH" == "safety05" || "$WHICH" == "both" ]] && probe "$p" 0.5
done

echo
echo "wrote $CSV"
