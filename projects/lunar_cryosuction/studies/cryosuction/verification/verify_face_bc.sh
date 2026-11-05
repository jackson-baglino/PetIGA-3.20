#!/usr/bin/env bash
# =============================================================================
# studies/cryosuction/verification/verify_face_bc.sh [BASE_REF]
#
# Gates for the opt-in face-list BC options -rhovfix_faces and -Tfix_faces
# (one-sided reservoirs: an open end and a dead pore). Like verify_wall_bc.sh
# these call the binary directly: they give pass/fail numbers, not results.
#
#   F1  Regression   -- a channel run with the legacy reservoir options
#                       (-rhovfix_axis 0, -flag_BC_Tfix, walls at 60 deg) is
#                       bit-for-bit the BASE_REF binary (default HEAD).
#   F2  Equivalence  -- -rhovfix_faces x0,x1 reproduces -rhovfix_axis 0, and
#                       -Tfix_faces x0,x1 with a gradient along x reproduces
#                       the legacy gradient logic, bit for bit.
#   F3  One-sided    -- -rhovfix_faces x0 -Tfix_faces x0: the printed boundary
#                       table must show x0 open with T and rho_v pinned, and x1
#                       sealed, insulating and dphi/dn = 0 (a dead pore); the
#                       run must differ from the two-sided one.
#   F4  Guards       -- misuse exits non-zero: -rhovfix_faces without
#                       -flag_BC_rhovfix, with -rhovfix_axis, and a -Tfix_faces
#                       face parallel to -grad_temp0.
#
# Writes verify_face_bc.csv and verify_face_bc.log next to this script.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXEC="$ROOT/lunar_cryosuction"
CSV="$HERE/verify_face_bc.csv"
LOG="$HERE/verify_face_bc.log"
BASE_REF="${1:-HEAD}"
WORK="${TMPDIR:-/tmp}/verify_face_bc.$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT
: > "$LOG"
echo "gate,case,result" > "$CSV"
fail=0
say() { printf '%s\n' "$*" | tee -a "$LOG"; }
rec() { echo "$1,$2,$3" >> "$CSV"; [[ $3 == FAIL ]] && fail=1; say "  $1 $2: $3"; }

[[ -x "$EXEC" ]] || { echo "missing $EXEC -- run 'make' first" >&2; exit 2; }

# Small wall-bounded channel with a bridging disc, reservoirs off saturation so
# the vapor BC matters, a gradient along x so the T BC matters.
COMMON=(-options_file "$ROOT/inputs/solver.opts" -periodic 0 -pf_monitor 1 -pf_output 0
        -dim 2 -Nx 24 -Ny 12 -Lx 6.25e-5 -Ly 3.125e-5 -eps 1.2e-6
        -ic_type multi_grains -ice_grain_cx 3.125e-5 -ice_grain_cy 1.5625e-5
        -ice_grain_R 2.3438e-5 -temp -20 -humidity 1.0
        -t_final 2.0e2 -delt_t 1.0e-2
        -wall_faces y0,y1 -contact_angle_deg 60
        -flag_BC_rhovfix 1 -rhovfix_lo 0.99 -rhovfix_hi 1.0
        -flag_BC_Tfix 1 -grad_temp0 100.0,0.0,0.0)

runcase() {  # name binary extra-args...
  local name=$1 bin=$2; shift 2
  local dir="$WORK/$name"; mkdir -p "$dir"
  folder="$dir" "$bin" "${COMMON[@]}" "$@" > "$dir/out.txt" 2>&1
  echo $? > "$dir/rc"
}
same() { cmp -s "$WORK/$1/SSA_evo.dat" "$WORK/$2/SSA_evo.dat" && [[ -s "$WORK/$1/SSA_evo.dat" ]]; }

say "BASE_REF = $BASE_REF ($(git -C "$ROOT" rev-parse --short=8 "$BASE_REF"))"

# ---------------------------------------------------------------- F1 ----
say ""; say "=== F1: legacy reservoir options, base vs new binary ==="
SNAP="$WORK/base"; mkdir -p "$SNAP"
if git -C "$ROOT/.." archive "$BASE_REF" lunar_cryosuction | tar -x -C "$SNAP" \
   && ( cd "$SNAP/lunar_cryosuction" && make >>"$LOG" 2>&1 ); then
  runcase f1_base "$SNAP/lunar_cryosuction/lunar_cryosuction" -rhovfix_axis 0
  runcase f1_new  "$EXEC" -rhovfix_axis 0
  if same f1_base f1_new; then rec F1 legacy_bitwise PASS; else rec F1 legacy_bitwise FAIL; fi
else
  rec F1 legacy_bitwise SKIP
fi

# ---------------------------------------------------------------- F2 ----
say ""; say "=== F2: face lists reproduce the legacy options ==="
runcase f2_legacy "$EXEC" -rhovfix_axis 0
runcase f2_rhov   "$EXEC" -rhovfix_faces x0,x1
runcase f2_T      "$EXEC" -rhovfix_axis 0 -Tfix_faces x0,x1
if same f2_legacy f2_rhov; then rec F2 rhovfix_faces_x0x1 PASS; else rec F2 rhovfix_faces_x0x1 FAIL; fi
if same f2_legacy f2_T;    then rec F2 Tfix_faces_x0x1    PASS; else rec F2 Tfix_faces_x0x1    FAIL; fi

# ---------------------------------------------------------------- F3 ----
say ""; say "=== F3: one-sided reservoir ==="
runcase f3_one "$EXEC" -rhovfix_faces x0 -Tfix_faces x0
grep -E "^ +[xy] = " "$WORK/f3_one/out.txt" | tee -a "$LOG" > "$WORK/f3_table.txt"
if [[ $(cat "$WORK/f3_one/rc") == 0 ]] && ! same f2_legacy f3_one; then
  rec F3 runs_and_differs PASS
else
  rec F3 runs_and_differs FAIL
fi
x0=$(grep "x = 0 " "$WORK/f3_table.txt"); x1=$(grep "x = Lx" "$WORK/f3_table.txt")
if [[ $x0 == *"open"*"Dirichlet"*"Dirichlet"* ]]; then rec F3 x0_open_and_pinned PASS; else rec F3 x0_open_and_pinned FAIL; fi
if [[ $x1 == *"dphi/dn=0"*"insulating"*"sealed"* ]]; then rec F3 x1_dead_end PASS; else rec F3 x1_dead_end FAIL; fi

# ---------------------------------------------------------------- F4 ----
say ""; say "=== F4: guards ==="
guard() {  # label extra-args...
  local label=$1; shift
  local dir="$WORK/f4_$label"; mkdir -p "$dir"
  if folder="$dir" "$EXEC" "${COMMON[@]}" -t_final 0 "$@" > "$dir/out.txt" 2>&1; then
    rec F4 "$label" FAIL
  else
    rec F4 "$label" PASS
  fi
}
guard no_flag      -flag_BC_rhovfix 0 -rhovfix_faces x0
guard with_axis    -rhovfix_axis 0 -rhovfix_faces x0
guard T_parallel   -Tfix_faces y0

say ""
if [[ $fail -eq 0 ]]; then say "ALL GATES PASSED   (see $CSV)"; else say "FAILURES -- see $CSV"; fi
exit $fail
