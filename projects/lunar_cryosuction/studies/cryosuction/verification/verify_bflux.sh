#!/usr/bin/env bash
# =============================================================================
# studies/cryosuction/verification/verify_bflux.sh [BASE_REF]
#
# Gates for bflux.dat (BoundaryFluxMonitor, output only). Like the other
# verification drivers these call the binary directly for pass/fail numbers.
#
#   H1  Regression  -- with reservoirs on (so the monitor is active), SSA_evo.dat
#                      is bit-for-bit the BASE_REF binary (default HEAD): the
#                      monitor reads the solution and must not change it.
#   H2  Analytic    -- an ice-free channel (the one IC grain sits 1 mm outside), vapor pinned at 0.99 rho_vs on x0
#                      and rho_vs on x1, isothermal. At steady state the flux is
#                      J = D_v(T0) * 0.01 rho_vs(T0) / Lx * Ly out of x0 and the
#                      same into x1. D_v and rho_vs come from preprocess/
#                      mesh_params (the Python mirror of the C laws), not from
#                      the solver. Tolerance 1e-3 relative (the mirror's
#                      constants agree with the C to about 1e-4).
#   H3  Unscaled    -- H2 again with -xi_v 1e-2: the physical flux must not
#                      depend on xi_v (same tolerance).
#   H4  Off switch  -- -bflux_output 0 writes no bflux.dat.
#
# Writes verify_bflux.csv and verify_bflux.log next to this script.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXEC="$ROOT/lunar_cryosuction"
PY="$ROOT/venv_cryosuction/bin/python"
CSV="$HERE/verify_bflux.csv"
LOG="$HERE/verify_bflux.log"
BASE_REF="${1:-HEAD}"
WORK="${TMPDIR:-/tmp}/verify_bflux.$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT
: > "$LOG"
echo "gate,case,value,expected,result" > "$CSV"
fail=0
say() { printf '%s\n' "$*" | tee -a "$LOG"; }
rec() { echo "$1,$2,$3,$4,$5" >> "$CSV"; [[ $5 == FAIL ]] && fail=1; say "  $1 $2: $3 (expected $4) $5"; }

[[ -x "$EXEC" ]] || { echo "missing $EXEC -- run 'make' first" >&2; exit 2; }
LX=6.25e-5; LY=3.125e-5; T0=-20
BASE=(-options_file "$ROOT/inputs/solver.opts" -periodic 0 -pf_monitor 1 -pf_output 0
      -dim 2 -Nx 24 -Ny 12 -Lx $LX -Ly $LY -eps 1.2e-6 -temp $T0 -humidity 1.0)

runcase() {  # name binary args...
  local name=$1 bin=$2; shift 2
  local dir="$WORK/$name"; mkdir -p "$dir"
  folder="$dir" "$bin" "$@" > "$dir/out.txt" 2>&1
}

say "BASE_REF = $BASE_REF ($(git -C "$ROOT" rev-parse --short=8 "$BASE_REF"))"

# ---------------------------------------------------------------- H1 ----
say ""; say "=== H1: monitor does not change the solution ==="
H1=("${BASE[@]}" -ic_type multi_grains -ice_grain_cx 3.125e-5 -ice_grain_cy 1.5625e-5
    -ice_grain_R 2.3438e-5 -t_final 2.0e2 -delt_t 1.0e-2
    -wall_faces y0,y1 -contact_angle_deg 60
    -flag_BC_rhovfix 1 -rhovfix_axis 0 -rhovfix_lo 0.99 -rhovfix_hi 1.0
    -flag_BC_Tfix 1 -grad_temp0 100.0,0.0,0.0)
SNAP="$WORK/base"; mkdir -p "$SNAP"
if git -C "$ROOT/.." archive "$BASE_REF" lunar_cryosuction | tar -x -C "$SNAP" \
   && ( cd "$SNAP/lunar_cryosuction" && make >>"$LOG" 2>&1 ); then
  runcase h1_base "$SNAP/lunar_cryosuction/lunar_cryosuction" "${H1[@]}"
  runcase h1_new  "$EXEC" "${H1[@]}"
  if cmp -s "$WORK/h1_base/SSA_evo.dat" "$WORK/h1_new/SSA_evo.dat" && [[ -s "$WORK/h1_new/SSA_evo.dat" ]]; then
    rec H1 ssa_bitwise identical identical PASS
  else
    rec H1 ssa_bitwise differs identical FAIL
  fi
  [[ -s "$WORK/h1_new/bflux.dat" ]] && rec H1 bflux_written yes yes PASS || rec H1 bflux_written no yes FAIL
else
  rec H1 ssa_bitwise build_failed identical SKIP
fi

# ---------------------------------------------------------- H2, H3 ----
EXPECT=$(cd "$ROOT" && "$PY" -I -c "
import sys; sys.path.insert(0, '.')
from preprocess.mesh_params import physics as p
print(p.dv_m2_per_s($T0) * 0.01 * p.rho_vs_sat_kg_per_m3($T0) / $LX * $LY)")
say ""; say "=== H2/H3: analytic steady flux, expected J = $EXPECT kg/(m s) ==="
H2=("${BASE[@]}" -ic_type multi_grains -ice_grain_cx -1.0e-3 -ice_grain_cy 1.5625e-5
    -ice_grain_R 1.0e-9 -t_final 2.0e1 -delt_t 1.0e-2
    -flag_BC_rhovfix 1 -rhovfix_faces x0,x1 -rhovfix_lo 0.99 -rhovfix_hi 1.0
    -flag_BC_Tfix 1)
runcase h2 "$EXEC" "${H2[@]}"
runcase h3 "$EXEC" "${H2[@]}" -xi_v 1e-2
for c in h2 h3; do
  G=$(echo "$c" | tr '[:lower:]' '[:upper:]')
  f="$WORK/$c/bflux.dat"
  if [[ ! -s $f ]]; then rec "$G" bflux_missing none "$EXPECT" FAIL; continue; fi
  # columns: t step J_x0 Q_x0 J_x1 Q_x1 (x faces first; y faces follow when pinned)
  read -r jx0 jx1 <<< "$(grep -v '^#' "$f" | tail -1 | awk '{print $3, $5}')"
  for pair in "x0 $jx0 $EXPECT" "x1 $jx1 -$EXPECT"; do
    set -- $pair
    ok=$(awk -v a="$2" -v b="$3" 'BEGIN{e=(a-b)/b; if(e<0)e=-e; print (e<1e-3)?"PASS":"FAIL"}')
    rec "$G" "J_$1" "$2" "$3" "$ok"
  done
done

# ---------------------------------------------------------------- H4 ----
say ""; say "=== H4: -bflux_output 0 ==="
runcase h4 "$EXEC" "${H2[@]}" -t_final 1.0 -bflux_output 0
[[ -e "$WORK/h4/bflux.dat" ]] && rec H4 no_file present absent FAIL || rec H4 no_file absent absent PASS

say ""
if [[ $fail -eq 0 ]]; then say "ALL GATES PASSED   (see $CSV)"; else say "FAILURES -- see $CSV"; fi
exit $fail
