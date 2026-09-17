#!/usr/bin/env bash
# =============================================================================
# verify_restart.sh — does -initial_cond resume a run seamlessly?
#
# Runs a small case to t1, restarts from its last snapshot WITH NO FLAGS beyond
# -initial_cond, and checks the continuation is indistinguishable from having
# never stopped:
#
#   1. the clock resumes at the snapshot's t          (not 0)
#   2. the step size resumes at the snapshot's dt     (not -delt_t = 1e-4)
#   3. the first NEW step lands on the same (t, dt) the uninterrupted run
#      reached at the corresponding step
#
# (2) is the one worth having a test for. A restart that begins at -delt_t has
# to climb six orders of magnitude back to the working step size through the
# NRmin/NRmax growth heuristic -- ~145 full nonlinear solves before any new
# physics happens, which on a production mesh is hours. Both values are read
# from the snapshot's own SSA_evo.dat, so neither has to be typed in and
# neither can be typed in wrong.
#
# The physics here is meaningless (48x48, 50 steps). Only the plumbing is
# under test.
#
#   ./studies/restart/verification/verify_restart.sh
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

COMMON=(-options_file inputs/solver.opts -dim 2 -Nx 48 -Ny 48 -Lx 1e-3 -Ly 1e-3
        -periodic 1 -ic_type ice_slab -eps 4e-5 -eps_temp_override 1
        -temp -20 -humidity 1 -delt_t 1.0e-4 -dtmax 5.0e1 -outp 1 -keff 0)

echo "== leg 1: run to t = 3.0e2 =="
mkdir -p "$WORK/leg1" "$WORK/leg2"
folder="$WORK/leg1" mpiexec -np 2 ./enceladus_dsm "${COMMON[@]}" \
    -t_final 3.0e2 -output_path "$WORK/leg1" > "$WORK/leg1.log" 2>&1

# Restart from the SECOND-TO-LAST snapshot, so leg 1 has a row to compare the
# continuation against.
last=$(ls "$WORK"/leg1/sol_*.dat | sort | tail -2 | head -1)
step=$(basename "$last" .dat | sed 's/sol_0*//')
read -r _ _ t0 _ dt0 _ <<< "$(awk -v s="$step" '$4==s' "$WORK/leg1/SSA_evo.dat" | tail -1)"
read -r _ _ t1 _ dt1 _ <<< "$(awk -v s="$((step+1))" '$4==s' "$WORK/leg1/SSA_evo.dat" | tail -1)"
echo "   snapshot step $step:  t = $t0  dt = $dt0"
echo "   uninterrupted next :  t = $t1  dt = $dt1"

echo "== leg 2: restart, no flags but -initial_cond =="
folder="$WORK/leg2" mpiexec -np 2 ./enceladus_dsm "${COMMON[@]}" \
    -t_final 5.0e2 -initial_cond "$last" -output_path "$WORK/leg2" \
    > "$WORK/leg2.log" 2>&1

read -r _ _ r0 _ rd0 _ <<< "$(awk '$4==0' "$WORK/leg2/SSA_evo.dat" | tail -1)"
read -r _ _ r1 _ rd1 _ <<< "$(awk '$4==1' "$WORK/leg2/SSA_evo.dat" | tail -1)"

fail=0
chk() {  # name expected actual tol
    local ok
    ok=$(awk -v a="$2" -v b="$3" -v t="$4" 'BEGIN{d=(a-b);if(d<0)d=-d;
         print (a==0 ? (d<t) : (d/(a<0?-a:a) < t)) ? "PASS" : "FAIL"}')
    printf "  [%s] %-34s expected %-14s got %s\n" "$ok" "$1" "$2" "$3"
    [ "$ok" = PASS ] || fail=1
}
echo "== gates =="
chk "clock resumed at snapshot t"  "$t0" "$r0" 1e-9
chk "dt resumed at snapshot dt"    "$dt0" "$rd0" 1e-9
chk "first new step matches t"     "$t1" "$r1" 1e-9
chk "first new step matches dt"    "$dt1" "$rd1" 1e-9
# the failure this is really guarding against
bad=$(awk -v a="$rd0" 'BEGIN{print (a < 1e-3) ? 1 : 0}')
[ "$bad" = 1 ] && { echo "  [FAIL] dt fell back to -delt_t -- the climb-back bug"; fail=1; }

cp "$WORK/leg2.log" "$HERE/restart.log" 2>/dev/null || true
echo
[ "$fail" -eq 0 ] && echo "all gates passed" || echo "FAILED"
exit $fail
