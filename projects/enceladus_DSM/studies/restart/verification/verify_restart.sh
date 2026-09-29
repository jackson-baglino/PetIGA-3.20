#!/usr/bin/env bash
# =============================================================================
# verify_restart.sh — does a resume IN PLACE look like the run never stopped?
#
# Two runs of the same small case:
#
#   A  uninterrupted, 0 -> t_end
#   B  0 -> t_stop, then resumed IN ITS OWN DIRECTORY from its second-to-last
#      snapshot -- exactly what scripts/HPC/resume_batch.sh does after a
#      timeout, including setting aside what B wrote past that snapshot --
#      and carried on to t_end
#
# and B's directory must be indistinguishable from A's:
#
#   1. SSA_evo.dat: same steps, one row each, same t and dt to 1e-9 and the
#      same totals to 1e-6 (the restart reloads U but not its time
#      derivative, so the totals may drift in the last digits -- report it)
#   2. the same sol_*.dat file names -- none overwritten, none off-cadence
#   3. the resumed dt is the working dt, not -delt_t (the climb-back bug)
#
# Before 2026-09-29 a resume restarted the step count at 0 and opened
# SSA_evo.dat for WRITE, so in place it overwrote sol_00000.. and truncated
# the first leg's log. Gates 1 and 2 are what catch that.
#
# B stops at t_stop via -t_final, so its LAST step is cut short to land on
# t_stop. Resuming from the second-to-last snapshot steps over that cut step,
# which is also why resume_batch.sh does it after a real kill.
#
# The physics here is meaningless (48x48, a few dozen steps). Only the
# plumbing is under test.
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
T_STOP=3.0e2
T_END=5.0e2
run() {  # dir log [extra...]
    local dir=$1 log=$2; shift 2
    folder="$dir" mpiexec -np 2 ./enceladus_dsm "${COMMON[@]}" "$@" \
        -output_path "$dir" > "$log" 2>&1
}

echo "== A: uninterrupted, to t = $T_END =="
mkdir -p "$WORK/A" "$WORK/B"
run "$WORK/A" "$WORK/A.log" -t_final "$T_END"

echo "== B: to t = $T_STOP, then resumed in place to t = $T_END =="
run "$WORK/B" "$WORK/B1.log" -t_final "$T_STOP"
snap=$(ls "$WORK"/B/sol_*.dat | sort | tail -2 | head -1)
step=$(basename "$snap" .dat | sed 's/sol_0*//')
read -r _ _ t0 _ dt0 _ <<< "$(awk -v s="$step" '$4==s' "$WORK/B/SSA_evo.dat" | tail -1)"
echo "   resuming from step $step:  t = $t0  dt = $dt0"
# resume_batch.sh's set-aside: later snapshots and log rows are abandoned.
for f in "$WORK"/B/sol_*.dat; do
    n=$(basename "$f" .dat | sed 's/sol_0*//'); [ "${n:-0}" -gt "$step" ] && rm "$f"
done
awk -v s="$step" 'NF<8 || $4+0 <= s+0' "$WORK/B/SSA_evo.dat" > "$WORK/ssa.tmp"
mv "$WORK/ssa.tmp" "$WORK/B/SSA_evo.dat"
run "$WORK/B" "$WORK/B2.log" -t_final "$T_END" -initial_cond "$snap"

echo "== gates =="
fail=0
PY_BIN="$ROOT/venv_enceladus/bin/python"; [ -x "$PY_BIN" ] || PY_BIN=python3
"$PY_BIN" - "$WORK/A" "$WORK/B" "$step" "$dt0" <<'PY' || fail=1
import sys, os, glob
import numpy as np
A, B, step, dt0 = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
a = np.loadtxt(os.path.join(A, "SSA_evo.dat")); b = np.loadtxt(os.path.join(B, "SSA_evo.dat"))
ok = True
def gate(name, cond, detail):
    global ok
    print(f"  [{'PASS' if cond else 'FAIL'}] {name:44s} {detail}")
    ok &= bool(cond)
sa, sb = a[:, 3].astype(int), b[:, 3].astype(int)
gate("one row per step in B", len(set(sb)) == len(sb), f"{len(sb)} rows, {len(set(sb))} steps")
gate("B has the same steps as A", np.array_equal(sa, sb), f"A {sa.min()}..{sa.max()}, B {sb.min()}..{sb.max()}")
if np.array_equal(sa, sb):
    rel = lambda c: float(np.max(np.abs(a[:, c] - b[:, c]) / np.maximum(np.abs(a[:, c]), 1e-300)))
    gate("same t on every step", rel(2) < 1e-9, f"max rel diff {rel(2):.2e}")
    gate("same dt on every step", rel(4) < 1e-9, f"max rel diff {rel(4):.2e}")
    tot = max(rel(c) for c in (1, 5, 6, 7))
    gate("same totals (ice, air, rhov, mass)", tot < 1e-6, f"max rel diff {tot:.2e}")
fa = sorted(os.path.basename(f) for f in glob.glob(os.path.join(A, "sol_*.dat")))
fb = sorted(os.path.basename(f) for f in glob.glob(os.path.join(B, "sol_*.dat")))
gate("same sol_*.dat files", fa == fb, f"A {len(fa)}, B {len(fb)}")
gate("resumed at the working dt, not -delt_t", dt0 > 1e-3, f"dt = {dt0:.3e} s")
sys.exit(0 if ok else 1)
PY

cp "$WORK/B2.log" "$HERE/restart.log" 2>/dev/null || true
echo
[ "$fail" -eq 0 ] && echo "all gates passed" || echo "FAILED"
exit $fail
