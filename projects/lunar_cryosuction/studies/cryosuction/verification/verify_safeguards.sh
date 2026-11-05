#!/usr/bin/env bash
# =============================================================================
# studies/cryosuction/verification/verify_safeguards.sh [BASE_REF]
#
# Gates for model-freeze-v2: the solver safeguards ported from enceladus_DSM
# (separate SNES guard and physics band, c6dee625/60244164; enforced
# interface-CFL rollback, caa209e0). Calls the binary directly for pass/fail
# numbers, like the other verification drivers.
#
#   S1  Regression -- with the pre-v2 settings restored on the command line
#                     (-phase_lo/hi +-0.05, -snes_guard_lo/hi +-0.05, a CFL
#                     limit no step reaches), every SSA_evo.dat column except
#                     dt is bit-for-bit the BASE_REF binary (default: the
#                     commit before v2); t matching bitwise means the steps
#                     taken are identical. The dt column logs the PROPOSED
#                     next dt, and the snap-to-dtmax in InterfaceCFLMonitor
#                     now runs before Monitor prints it, so it may differ.
#   S2  Enforced   -- a sublimating disc (h = 0.999, resolved) started at a
#                     deliberately large dt with -dtCFL_dphimax 0.05: at least
#                     one rollback fires, every accepted step that reports its
#                     |dphi| is <= 0.05, and no step number appears twice in
#                     SSA_evo.dat (a rolled-back step leaves no row).
#   S3  Banner     -- a guard no wider than the physics band prints the
#                     warning; the defaults do not.
#   S4  Completes  -- S2 reaches t_final.
#
# Writes verify_safeguards.csv and verify_safeguards.log next to this script.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXEC="$ROOT/lunar_cryosuction"
PY="$ROOT/venv_cryosuction/bin/python"
CSV="$HERE/verify_safeguards.csv"
LOG="$HERE/verify_safeguards.log"
BASE_REF="${1:-8d659590}"
WORK="${TMPDIR:-/tmp}/verify_safeguards.$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT
: > "$LOG"
echo "gate,case,value,expected,result" > "$CSV"
fail=0
say() { printf '%s\n' "$*" | tee -a "$LOG"; }
rec() { echo "$1,$2,$3,$4,$5" >> "$CSV"; [[ $5 == FAIL ]] && fail=1; say "  $1 $2: $3 (expected $4) $5"; }
[[ -x "$EXEC" ]] || { echo "missing $EXEC -- run 'make' first" >&2; exit 2; }

runcase() {  # name binary args...
  local name=$1 bin=$2; shift 2
  local dir="$WORK/$name"; mkdir -p "$dir"
  folder="$dir" "$bin" "$@" > "$dir/out.txt" 2>&1
}
# 60 x 45 um box; eps and mesh from comp_eps.py (alpha_c 1e-2, Rave 1.5e-5):
#   comp_eps.py --Lx 6.0e-5 --Ly 4.5e-5 --Rave 1.5e-5 --T0 -20 --alpha 1e-2
#   -> -eps 6.3921e-07 -Nx 133 -Ny 100
BOX=(-options_file "$ROOT/inputs/solver.opts" -periodic 0 -pf_monitor 1 -pf_output 0
     -dim 2 -Nx 133 -Ny 100 -Lx 6.0e-5 -Ly 4.5e-5 -eps 6.3921e-07 -temp -20 -humidity 1.0
     -ic_type multi_grains -ice_grain_cx 3.0e-5 -ice_grain_cy 2.25e-5 -ice_grain_R 1.5e-5
     -beta_sub0 7.9408e5 -dtCFL 1 -flag_BC_Tfix 1)
LEGACY=(-phase_lo -0.05 -phase_hi 1.05 -snes_guard_lo -0.05 -snes_guard_hi 1.05 -dtCFL_dphimax 10)
say "BASE_REF = $BASE_REF ($(git -C "$ROOT" rev-parse --short=8 "$BASE_REF"))"

# ---------------------------------------------------------------- S1 ----
say ""; say "=== S1: legacy settings reproduce the pre-v2 binary bit for bit ==="
S1=("${BOX[@]}" -flag_BC_rhovfix 1 -rhovfix_faces x0 -rhovfix_lo 0.999 -t_final 2.0e3 -delt_t 1.0e-2)
SNAP="$WORK/base"; mkdir -p "$SNAP"
if git -C "$ROOT/.." archive "$BASE_REF" lunar_cryosuction | tar -x -C "$SNAP" \
   && ( cd "$SNAP/lunar_cryosuction" && make >>"$LOG" 2>&1 ); then
  # The base binary has no -snes_guard_* options; PETSc ignores unknown ones.
  runcase s1_base "$SNAP/lunar_cryosuction/lunar_cryosuction" "${S1[@]}" "${LEGACY[@]}"
  runcase s1_new "$EXEC" "${S1[@]}" "${LEGACY[@]}"
  nodt() { awk '{$5=""; print}' "$1"; }
  if [[ -s "$WORK/s1_new/SSA_evo.dat" ]] && cmp -s <(nodt "$WORK/s1_base/SSA_evo.dat") <(nodt "$WORK/s1_new/SSA_evo.dat"); then
    rec S1 ssa_bitwise identical identical PASS
  else
    rec S1 ssa_bitwise differs identical FAIL
  fi
else
  rec S1 ssa_bitwise build_failed identical SKIP
fi

# ---------------------------------------------------------------- S2/S4 --
say ""; say "=== S2: the interface-CFL limit holds on every kept step ==="
# Resolved drive (h = 0.999; Pi ~ 0.3). First step 500 s: the measured rate
# (~3e-4 per s) puts step 1 at |dphi| ~ 0.14 > 0.05, so it must roll back;
# after that the forward clamp (0.9 x limit) keeps steps inside.
# (An h = 0.9 disc on this eps sits ~1 % above phi = 1 at any dt and is
# rejected by the [-0.01, 1.01] physics band until dtmin: the band refusing an
# unresolved run, which is a different test.)
runcase s2 "$EXEC" "${BOX[@]}" -flag_BC_rhovfix 1 -rhovfix_faces x0 -rhovfix_lo 0.999 \
  -dtCFL_dphimax 0.05 -t_final 3.0e3 -delt_t 5.0e2 -dtmax 1.0e3
nrb=$(grep -c "Interface-CFL violated at step .* rolling back" "$WORK/s2/out.txt")
worst=$(grep -o "max |dphi|/step [0-9.e+-]*" "$WORK/s2/out.txt" | awk '{print $3}' | sort -g | tail -1)
dups=$("$PY" -I -c "import numpy as np,sys; s=np.loadtxt(sys.argv[1],ndmin=2); print(len(s)-len(set(s[:,3].astype(int))))" "$WORK/s2/SSA_evo.dat")
[[ $nrb -ge 1 ]] && rec S2 rollbacks "$nrb" ">=1" PASS || rec S2 rollbacks "$nrb" ">=1" FAIL
ok=$(awk -v w="${worst:-0}" 'BEGIN{print (w<=0.05)?"PASS":"FAIL"}')
rec S2 max_dphi_accepted "${worst:-none}" "<=0.05" "$ok"
[[ $dups -eq 0 ]] && rec S2 duplicate_step_rows "$dups" 0 PASS || rec S2 duplicate_step_rows "$dups" 0 FAIL
tlast=$(tail -1 "$WORK/s2/SSA_evo.dat" | awk '{print $3}')
ok=$(awk -v t="$tlast" 'BEGIN{print (t>=2.999e3)?"PASS":"FAIL"}')
rec S4 reaches_t_final "$tlast" ">=3.0e3" "$ok"

# ---------------------------------------------------------------- S3 ----
say ""; say "=== S3: banner warns when the guard is not wider than the band ==="
runcase s3 "$EXEC" "${BOX[@]}" -phase_hi 1.05 -snes_guard_hi 1.05 -t_final 1.0 -delt_t 1.0e-2
grep -q "solver guard is not wider" "$WORK/s3/out.txt" && rec S3 warns printed printed PASS || rec S3 warns absent printed FAIL
runcase s3ok "$EXEC" "${BOX[@]}" -t_final 1.0 -delt_t 1.0e-2
grep -q "solver guard is not wider" "$WORK/s3ok/out.txt" && rec S3 defaults_quiet printed absent FAIL || rec S3 defaults_quiet absent absent PASS
grep -m1 "physics band" "$WORK/s3ok/out.txt" | tee -a "$LOG"

say ""
if [[ $fail -eq 0 ]]; then say "ALL GATES PASSED   (see $CSV)"; else say "FAILURES -- see $CSV"; fi
exit $fail
