#!/usr/bin/env bash
# =============================================================================
# verify_kinetics_scaling.sh — gate for -mob_scale / -alph_scale.
#
# WHY THIS EXISTS
# ---------------
# Under -alpha_pointwise 1 (what every Molaro production run uses) mob_sub and
# alph_sub are rebuilt at every quadrature point in SubKinetics(), so the
# pre-existing -mob_sub / -alph_sub ABSOLUTE overrides in main() are silently
# ignored — exactly like -beta_sub0. -mob_scale / -alph_scale are multipliers
# applied on BOTH paths so option 3 of the Molaro three-option study can scale
# M_0 and alph_sub independently without leaving the production code path.
#
# Three things have to hold, and only a run can show them:
#   1. DEFAULT IS A NO-OP. Omitting the flags and passing 1.0 must give
#      byte-identical kinetics, or every earlier result is invalidated.
#   2. THE FACTORS BITE, on the pointwise path. mob_sub must come out 5x and
#      alph_sub 0.01x of baseline when asked -- this is the whole point, and
#      it is what -mob_sub/-alph_sub fail to do.
#   3. THE JACOBIAN STAYS EXACT. The scaling is linear, so each derivative
#      carries its own factor and -snes_test_jacobian must still agree. If it
#      does not, a derivative was missed in SubKinetics()'s tail.
#
# Runs on a deliberately tiny mesh (32 x 16, eps coarsened to match) because
# -snes_test_jacobian builds the finite-difference Jacobian one column per DoF.
# The physics is meaningless at this resolution; the derivative algebra is not,
# and that is all this gate tests. Direct ./enceladus_dsm invocation is the
# unit-gate exception in CLAUDE.md — this script IS the reproducible record.
#
# Usage:  bash studies/molaro_2019/verification/verify_kinetics_scaling.sh
# Writes: kinetics_scaling.csv and kinetics_scaling.log next to this script.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXEC="$ROOT/enceladus_dsm"
LOG="$HERE/kinetics_scaling.log"
CSV="$HERE/kinetics_scaling.csv"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

GEOM="$ROOT/inputs/geometry/molaro/molaro_2D_L387um_eps0.87um_axisym_T-20pair_tangent.opts"
EXP="$ROOT/inputs/experiment/molaro/molaro_T-20_h0.99816_2h_a1e-1_dirichlet.opts"
SOLVER="$ROOT/inputs/solver.opts"

# 32 x 16 over Lx = 3.87e-4 gives dx = 1.21e-5; eps = dx/sqrt(2) keeps the
# interface nominally resolved so no term in the residual is degenerate.
COARSE=(-Nx 32 -Ny 16 -eps 8.55e-6 -delt_t 1.0e-4 -ts_max_steps 1
        -pf_output 0 -pf_monitor 0 -t_out_log 0 -outp 0)

[[ -x "$EXEC" ]] || { echo "FAIL: $EXEC not built. Run 'make all' first."; exit 1; }

run_case () {           # run_case <name> [extra flags...]
    local name="$1"; shift
    "$EXEC" -options_file "$SOLVER" -options_file "$GEOM" -options_file "$EXP" \
            -output_path "$WORK" "${COARSE[@]}" "$@" > "$WORK/$name.out" 2>&1
    echo "$WORK/$name.out"
}

# The pointwise banner block prints these AFTER the '--- pointwise' header;
# anchor on that so the scalar echoes higher up cannot be picked up by mistake.
kin () {                # kin <file> <mob_sub|alph_sub>
    sed -n '/--- pointwise/,$p' "$1" | grep -m1 "^   $2 " | awk '{print $3}'
}

: > "$LOG"
echo "run,mob_scale,alph_scale,mob_sub,alph_sub" > "$CSV"
fail=0

echo "=== 1. baseline (flags omitted) vs explicit 1.0 ===" | tee -a "$LOG"
f_base=$(run_case baseline)
f_one=$(run_case unity -mob_scale 1.0 -alph_scale 1.0)
mb=$(kin "$f_base" mob_sub);  ab=$(kin "$f_base" alph_sub)
mo=$(kin "$f_one"  mob_sub);  ao=$(kin "$f_one"  alph_sub)
[[ -n "$mb" && -n "$ab" ]] || { echo "FAIL: could not parse the baseline banner" | tee -a "$LOG"; cat "$f_base" >> "$LOG"; exit 1; }
echo "  baseline  mob_sub=$mb  alph_sub=$ab" | tee -a "$LOG"
echo "  unity     mob_sub=$mo  alph_sub=$ao" | tee -a "$LOG"
echo "baseline,,,$mb,$ab"   >> "$CSV"
echo "unity,1.0,1.0,$mo,$ao" >> "$CSV"
if [[ "$mb" == "$mo" && "$ab" == "$ao" ]]; then
    echo "  PASS: the default is a no-op" | tee -a "$LOG"
else
    echo "  FAIL: -mob_scale/-alph_scale 1.0 changed the kinetics" | tee -a "$LOG"; fail=1
fi

echo "=== 2. -mob_scale 5 -alph_scale 0.01 on the pointwise path ===" | tee -a "$LOG"
f_sc=$(run_case scaled -mob_scale 5.0 -alph_scale 0.01)
ms=$(kin "$f_sc" mob_sub); as=$(kin "$f_sc" alph_sub)
echo "  scaled    mob_sub=$ms  alph_sub=$as" | tee -a "$LOG"
echo "scaled,5.0,0.01,$ms,$as" >> "$CSV"
# Tolerance is 1e-3, not machine epsilon, because these numbers come off the
# banner, which prints %.4e: a ratio built from two 5-significant-figure values
# is only good to ~1e-4 relative. Anything tighter fails on print precision
# alone (5.0 came out as 4.999822 the first time this ran).
python3 - "$mb" "$ab" "$ms" "$as" <<'RATIOS' | tee -a "$LOG"
import sys
mb, ab, ms, as_ = (float(x) for x in sys.argv[1:5])
ok = True
for name, got, want in (("mob_sub", ms/mb, 5.0), ("alph_sub", as_/ab, 0.01)):
    rel = abs(got - want)/want
    ok &= rel < 1e-3
    print("  %-9s ratio %.6f vs %.4g  (rel %.2e)  %s"
          % (name, got, want, rel, "PASS" if rel < 1e-3 else "FAIL"))
sys.exit(0 if ok else 1)
RATIOS
[[ ${PIPESTATUS[0]} -eq 0 ]] || fail=1

echo "=== 3. -snes_test_jacobian with the scalings on ===" | tee -a "$LOG"
f_j=$(run_case jacobian -mob_scale 5.0 -alph_scale 0.01 -snes_test_jacobian)
# PETSc writes "||J - Jfd||_F/||J||_F = <x>, ||J - Jfd||_F = <y>" once per
# Newton iteration. Gate on the WORST of them, not the first.
ratio=$(grep -o '||J - Jfd||_F/||J||_F = [0-9.eE+-]*' "$f_j" \
        | awk '{print $NF}' | sort -g | tail -1)
if [[ -z "$ratio" ]]; then
    echo "  FAIL: -snes_test_jacobian produced no comparison" | tee -a "$LOG"
    grep -i -m5 'jacobian\|error' "$f_j" >> "$LOG"; fail=1
else
    njac=$(grep -c '||J - Jfd||_F/||J||_F' "$f_j")
    echo "  worst ||J - Jfd||_F/||J||_F = $ratio  (over $njac Newton iterations)" | tee -a "$LOG"
    echo "jacobian,5.0,0.01,,$ratio" >> "$CSV"
    # PETSc's own guidance is that O(1e-8) means "probably correct"; this
    # solver lands near 1e-9, so 1e-7 is a real gate and not a rubber stamp.
    python3 -c "import sys; sys.exit(0 if $ratio < 1e-7 else 1)" \
        && echo "  PASS: analytic Jacobian still agrees with the FD one" | tee -a "$LOG" \
        || { echo "  FAIL: $ratio exceeds 1e-7 -- a derivative in SubKinetics() is unscaled" | tee -a "$LOG"; fail=1; }
fi

echo | tee -a "$LOG"
if [[ $fail -eq 0 ]]; then echo "ALL CHECKS PASSED  -> $CSV" | tee -a "$LOG"
else echo "SOME CHECKS FAILED -> $LOG" | tee -a "$LOG"; fi
exit $fail
