#!/usr/bin/env bash
# =============================================================================
# studies/cryosuction/verification/verify_resolution.sh [BASE_REF]
#
# Gates for the resolution checker (ResolutionMonitor / ResolutionStaticCheck,
# output only). Calls the binary directly for pass/fail numbers, like the
# other verification drivers.
#
#   K1  Regression -- SSA_evo.dat bit-for-bit the BASE_REF binary (default
#                     HEAD) on a run with the checker active.
#   K2  Speed      -- a disc under pure curvature flow: the checker's v_n
#                     (the LARGEST |dphi/dt|/|grad phi| in the band) against
#                     the mean front speed from the area, v = -(dA/dt)/(2 pi R),
#                     R = sqrt(A/pi), from SSA_evo.dat over the same step. The
#                     maximum sits on the inner side of the band (r ~ R - 2 eps),
#                     where curvature flow is faster (v ~ D/r), so the checker
#                     reads high by up to ~2 eps/R (about 10 % here): it is
#                     conservative. Pass: 1.0 <= v_chk/v_area <= 1.15.
#   K3  Warning    -- a fast sublimating disc on a coarse eps (alpha_c 1e-3,
#                     h = 0.9) must print the RESOLUTION WARNING and count it.
#   K4  Static     -- a mesh 2x coarser than eps/sqrt(2) must print the
#                     start-up RESOLUTION WARNING.
#   K5  Reservoir  -- a disc cut by a vapor-reservoir face must print the
#                     RESERVOIR-FACE WARNING; the same disc fully inside must
#                     not.
#   K6  Watch      -- the warning does not depend on the phi condition: the
#                     same cut disc with -phi_open_bc 0 (90-degree face) still
#                     prints it.
#
# Writes verify_resolution.csv and verify_resolution.log next to this script.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXEC="$ROOT/lunar_cryosuction"
PY="$ROOT/venv_cryosuction/bin/python"
CSV="$HERE/verify_resolution.csv"
LOG="$HERE/verify_resolution.log"
BASE_REF="${1:-HEAD}"
WORK="${TMPDIR:-/tmp}/verify_resolution.$$"
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
# 60 x 45 um box, eps and mesh from comp_eps.py (alpha_c 1e-2, Rave 1.5e-5):
#   comp_eps.py --Lx 6.0e-5 --Ly 4.5e-5 --Rave 1.5e-5 --T0 -20 --alpha 1e-2
#   -> -eps 6.3921e-07 -Nx 133 -Ny 100
BOX=(-options_file "$ROOT/inputs/solver.opts" -periodic 0 -pf_monitor 1 -pf_output 0
     -dim 2 -Nx 133 -Ny 100 -Lx 6.0e-5 -Ly 4.5e-5 -eps 6.3921e-07 -temp -20 -humidity 1.0
     -ic_type multi_grains -ice_grain_cy 2.25e-5 -ice_grain_R 1.5e-5
     -beta_sub0 7.9408e5 -dtCFL 1 -flag_BC_Tfix 1)
say "BASE_REF = $BASE_REF ($(git -C "$ROOT" rev-parse --short=8 "$BASE_REF"))"

# ---------------------------------------------------------------- K1 ----
say ""; say "=== K1: checker does not change the solution ==="
K1=("${BOX[@]}" -ice_grain_cx 3.0e-5 -flag_BC_rhovfix 1 -rhovfix_faces x0 -rhovfix_lo 0.999
    -t_final 2.0e3 -delt_t 1.0e-2)
SNAP="$WORK/base"; mkdir -p "$SNAP"
if git -C "$ROOT/.." archive "$BASE_REF" lunar_cryosuction | tar -x -C "$SNAP" \
   && ( cd "$SNAP/lunar_cryosuction" && make >>"$LOG" 2>&1 ); then
  runcase k1_base "$SNAP/lunar_cryosuction/lunar_cryosuction" "${K1[@]}"
  runcase k1_new "$EXEC" "${K1[@]}"
  if cmp -s "$WORK/k1_base/SSA_evo.dat" "$WORK/k1_new/SSA_evo.dat" && [[ -s "$WORK/k1_new/SSA_evo.dat" ]]; then
    rec K1 ssa_bitwise identical identical PASS
  else
    rec K1 ssa_bitwise differs identical FAIL
  fi
else
  rec K1 ssa_bitwise build_failed identical SKIP
fi

# ---------------------------------------------------------------- K2 ----
say ""; say "=== K2: measured v_n against the area-derived speed ==="
runcase k2 "$EXEC" "${BOX[@]}" -ice_grain_cx 3.0e-5 -decouple_phase_change 1 -t_final 2.0e4 -delt_t 1.0e-2 -dtmax 1.0e3
res=$("$PY" -I - "$WORK/k2" <<'EOF'
import sys, math, numpy as np
d = sys.argv[1]
s = np.loadtxt(d + "/SSA_evo.dat", ndmin=2); r = np.loadtxt(d + "/resolution.dat", ndmin=2)
t, A = s[:, 2], s[:, 1]
i = len(t) - 1                       # last accepted step
R = math.sqrt(A[i] / math.pi)
v_area = -(A[i] - A[i - 1]) / (t[i] - t[i - 1]) / (2 * math.pi * R)
j = np.argmin(abs(r[:, 0] - t[i]))
print(f"{r[j, 3]:.6e} {v_area:.6e} {r[j, 3] / v_area:.4f}")
EOF
)
read -r vchk varea ratio <<< "$res"
ok=$(awk -v e="$ratio" 'BEGIN{print (e>=1.0 && e<=1.15)?"PASS":"FAIL"}')
rec K2 vn_vs_area "$vchk (ratio $ratio)" "$varea" "$ok"

# ---------------------------------------------------------------- K3 ----
say ""; say "=== K3: fast front on a coarse eps prints the warning ==="
runcase k3 "$EXEC" "${BOX[@]}" -ice_grain_cx 3.0e-5 -beta_sub0 7.9408e6 -eps 2.5e-6 -Nx 34 -Ny 26 \
  -flag_BC_rhovfix 1 -rhovfix_faces x0 -rhovfix_lo 0.9 -humidity 0.9 -t_final 5.0e3 -delt_t 1.0e-2 -dtmax 5.0e2
if grep -q "RESOLUTION WARNING step" "$WORK/k3/out.txt"; then rec K3 warning printed printed PASS; else rec K3 warning absent printed FAIL; fi
grep -m1 "RESOLUTION SUMMARY" "$WORK/k3/out.txt" | tee -a "$LOG"

# ---------------------------------------------------------------- K4 ----
say ""; say "=== K4: coarse mesh trips the start-up check ==="
runcase k4 "$EXEC" "${BOX[@]}" -ice_grain_cx 3.0e-5 -Nx 66 -Ny 50 -t_final 1.0 -delt_t 1.0e-2
if grep -q "RESOLUTION WARNING: mesh coarser" "$WORK/k4/out.txt"; then rec K4 static printed printed PASS; else rec K4 static absent printed FAIL; fi
runcase k4ok "$EXEC" "${BOX[@]}" -ice_grain_cx 3.0e-5 -t_final 1.0 -delt_t 1.0e-2
if grep -q "RESOLUTION WARNING: mesh coarser" "$WORK/k4ok/out.txt"; then rec K4 comp_eps_mesh_quiet printed absent FAIL; else rec K4 comp_eps_mesh_quiet absent absent PASS; fi

# ---------------------------------------------------------------- K5 ----
say ""; say "=== K5: ice on a vapor-reservoir face ==="
runcase k5 "$EXEC" "${BOX[@]}" -ice_grain_cx -7.5e-6 -decouple_phase_change 1 -flag_BC_rhovfix 1 -rhovfix_faces x0 \
  -t_final 1.0e3 -delt_t 1.0e-2 -dtmax 1.0e2
if grep -q "RESERVOIR-FACE WARNING" "$WORK/k5/out.txt"; then rec K5 touching printed printed PASS; else rec K5 touching absent printed FAIL; fi
runcase k5ok "$EXEC" "${BOX[@]}" -ice_grain_cx 3.0e-5 -decouple_phase_change 1 -flag_BC_rhovfix 1 -rhovfix_faces x0 \
  -t_final 1.0e3 -delt_t 1.0e-2 -dtmax 1.0e2
if grep -q "RESERVOIR-FACE WARNING" "$WORK/k5ok/out.txt"; then rec K5 clear printed absent FAIL; else rec K5 clear absent absent PASS; fi
runcase k6 "$EXEC" "${BOX[@]}" -ice_grain_cx -7.5e-6 -decouple_phase_change 1 -flag_BC_rhovfix 1 -rhovfix_faces x0 \
  -phi_open_bc 0 -t_final 1.0e3 -delt_t 1.0e-2 -dtmax 1.0e2
if grep -q "RESERVOIR-FACE WARNING" "$WORK/k6/out.txt"; then rec K6 closed_face_touching printed printed PASS; else rec K6 closed_face_touching absent printed FAIL; fi

say ""
if [[ $fail -eq 0 ]]; then say "ALL GATES PASSED   (see $CSV)"; else say "FAILURES -- see $CSV"; fi
exit $fail
