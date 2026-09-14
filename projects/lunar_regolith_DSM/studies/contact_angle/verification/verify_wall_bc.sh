#!/usr/bin/env bash
# =============================================================================
# studies/contact_angle/verification/verify_wall_bc.sh
#
# Unit gates for the prescribed-contact-angle wall free-energy term. These are
# the CLAUDE.md-sanctioned exception to "always go through the run script":
# they produce no simulation result worth reading, only pass/fail numbers, and
# this script is the committed reproducible record of them.
#
#   G1  Regression       -- without -wall_faces the solver is bit-for-bit what
#                           it was before the feature existed. Compares against
#                           a binary built from the merge-base with main.
#   G2  Wall Jacobian    -- the analytic wall block matches finite differences.
#                           The wall block is ISOLATED by differencing J and F
#                           between costhet on and off, because the full system
#                           is dominated by the interior form's 3M/eps terms and
#                           is insensitive to the boundary block on its own.
#   G3  Surface measure  -- on a uniform phi=1/2 field at saturation every
#                           interior contribution to R[.][0] vanishes
#                           identically, so the assembled vector is the wall
#                           term alone and must integrate to
#                           -3*M*(1/4)*cos(theta)*|Gamma_wall|. This is the one
#                           thing that cannot be checked by reading the code:
#                           on a plain Cartesian patch PetIGA takes the
#                           detS = 1.0 branch rather than computing a geometric
#                           surface Jacobian.
#
# Writes verify_wall_bc.csv and verify_wall_bc.log next to this script.
# Exits non-zero if any gate fails.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXEC="$ROOT/lunar_regolith_dsm"
CSV="$HERE/verify_wall_bc.csv"
LOG="$HERE/verify_wall_bc.log"
WORK="${TMPDIR:-/tmp}/verify_wall_bc.$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

# Tolerances. G2 is bounded by central-difference truncation at h=1e-7, not by
# the Jacobian, so 1e-6 is loose enough to be stable and tight enough to catch
# a wrong coefficient (which shows up at O(1)).
TOL_JAC=1e-6
TOL_MEASURE=1e-10

export folder="$WORK"          # the solver reads this for SSA_evo.dat etc.
: > "$LOG"
echo "gate,case,quantity,value,tolerance,result" > "$CSV"
fail=0

say() { printf '%s\n' "$*" | tee -a "$LOG"; }
run() { "$EXEC" "$@" >>"$LOG" 2>&1; }

[[ -x "$EXEC" ]] || { echo "missing $EXEC -- run 'make' first" >&2; exit 2; }

# A disc clipped by BOTH walls: R = 0.75*Ly > Ly/2. The ice-vapour interface
# must actually cross a wall face or G2 is vacuous (the wall term is ~0 where
# h'(phi) ~ 0, i.e. in either bulk phase).
GEOM=(-options_file "$ROOT/inputs/solver.opts" -periodic 0 -pf_output 0 -pf_monitor 0
      -dim 2 -Nx 24 -Ny 12 -Lx 6.25e-5 -Ly 3.125e-5 -eps 1.2e-6
      -ic_type multi_grains -ice_grain_cx 3.125e-5 -ice_grain_cy 1.5625e-5
      -ice_grain_R 2.3438e-5 -temp -20 -humidity 1.0)

# -------------------------------------------------------------------- G2 ----
say ""; say "=== G2: wall Jacobian vs finite differences ==="
say "    (wall block isolated by differencing costhet on/off)"
for th in 0 30 45 60 90 120 150 180; do
  out="$WORK/jac_$th.txt"
  run "${GEOM[@]}" -t_final 0 -test_wall_jacobian \
      -wall_faces y0,y1 -contact_angle_deg "$th" || true
  worst=$(grep -E '^   worst' "$LOG" | tail -1 | awk '{print $3}')
  scale=$(grep -E '^     4 ' "$LOG" | tail -1 | awk '{print $4}')
  if [[ $th == 90 ]]; then
    # cos(90) = 0, so the wall block must be EXACTLY zero -- a relative error is
    # 0/0 here. This is the no-op gate: enabling the boundary form on a face and
    # then asking for 90 degrees must contribute nothing at all, which is what
    # makes theta=90 reproduce the pre-feature solver (G1).
    ok=$(awk -v s="$scale" 'BEGIN{print (s!="" && s+0<1e-25)?"PASS":"FAIL"}')
    [[ $ok == FAIL ]] && fail=1
    printf '  theta=%-4s  wall block is zero (no-op): ||Jbnd*v||=%-12s %s\n' "$th" "$scale" "$ok" | tee -a "$LOG"
    echo "G2,theta=$th,norm_wall_block_must_be_zero,$scale,1e-25,$ok" >> "$CSV"
  else
    ok=$(awk -v w="$worst" -v t="$TOL_JAC" -v s="$scale" \
          'BEGIN{print (w!="" && w+0<t+0 && s+0>0) ? "PASS" : "FAIL"}')
    [[ $ok == FAIL ]] && fail=1
    printf '  theta=%-4s  rel.err=%-14s  ||Jbnd*v||=%-14s  %s\n' "$th" "$worst" "$scale" "$ok" | tee -a "$LOG"
    echo "G2,theta=$th,rel_err_wall_block,$worst,$TOL_JAC,$ok" >> "$CSV"
  fi
done

# Same thing reached through Young's equation rather than -contact_angle_deg.
run "${GEOM[@]}" -t_final 0 -test_wall_jacobian \
    -wall_faces y0,y1 -gamma_ia 0.109 -gamma_as 0.300 -gamma_is 0.24550 || true
worst=$(grep -E '^   worst' "$LOG" | tail -1 | awk '{print $3}')
ok=$(awk -v w="$worst" -v t="$TOL_JAC" 'BEGIN{print (w!="" && w+0<t+0)?"PASS":"FAIL"}')
[[ $ok == FAIL ]] && fail=1
printf '  gammas(60deg) rel.err=%-14s  %s\n' "$worst" "$ok" | tee -a "$LOG"
echo "G2,gammas_theta60,rel_err_wall_block,$worst,$TOL_JAC,$ok" >> "$CSV"

# -------------------------------------------------------------------- G3 ----
say ""; say "=== G3: boundary surface measure ==="
MEAS=(-options_file "$ROOT/inputs/solver.opts" -periodic 0 -pf_output 0 -pf_monitor 0
      -dim 2 -Lx 6.25e-5 -Ly 3.125e-5 -eps 3.0e-6 -ic_type single_ice -RCice 1.0e-5
      -temp -20 -humidity 1.0 -t_final 0 -test_wall_measure)
# Mesh-independence, face selection, and angle scaling all in one pass.
#   Lx = 6.25e-5, Ly = 3.125e-5  =>  y-pair 1.25e-4, y-single 6.25e-5, x-pair 6.25e-5
while read -r label nx ny faces th; do
  [[ -z "$label" ]] && continue
  run "${MEAS[@]}" -Nx "$nx" -Ny "$ny" -wall_faces "$faces" -contact_angle_deg "$th" || true
  err=$(grep -E 'rel\. error' "$LOG" | tail -1 | awk '{print $4}')
  area=$(grep -E '\|Gamma_wall\|' "$LOG" | tail -1 | awk '{print $3}')
  ok=$(awk -v e="$err" -v t="$TOL_MEASURE" 'BEGIN{print (e!="" && e+0<t+0)?"PASS":"FAIL"}')
  [[ $ok == FAIL ]] && fail=1
  printf '  %-16s |Gamma|=%-16s rel.err=%-12s %s\n' "$label" "$area" "$err" "$ok" | tee -a "$LOG"
  echo "G3,$label,rel_err_measure,$err,$TOL_MEASURE,$ok" >> "$CSV"
done <<'CASES'
mesh_10x6        10  6 y0,y1       0
mesh_20x12       20 12 y0,y1       0
mesh_41x23       41 23 y0,y1       0
one_wall         20 12 y0          0
x_faces          20 12 x0,x1       0
theta60          20 12 y0,y1      60
theta135         20 12 y0,y1     135
all_faces        20 12 x0,x1,y0,y1 0
CASES

# -------------------------------------------------------------------- G4 ----
# Direction. A correct magnitude with a flipped sign would pass G2 and G3 and
# then quietly drive every validation run the wrong way, so check it explicitly.
# R = N*phi_t + ... = 0, so a NEGATIVE wall contribution gives phi_t > 0: ice
# grows at the wall and the contact line ADVANCES, which is what a wetting wall
# (cos theta > 0, i.e. gamma_as > gamma_is) must do.
say ""; say "=== G4: sign -- wetting advances the contact line ==="
while read -r th want; do
  [[ -z "$th" ]] && continue
  run "${MEAS[@]}" -Nx 20 -Ny 12 -wall_faces y0,y1 -contact_angle_deg "$th" || true
  v=$(grep -E 'sum F' "$LOG" | tail -1 | awk '{print $4}')
  ok=$(awk -v v="$v" -v w="$want" 'BEGIN{
         if (v=="") {print "FAIL"; exit}
         if (w=="neg") print (v+0 < -1e-20) ? "PASS" : "FAIL";
         else if (w=="pos") print (v+0 > 1e-20) ? "PASS" : "FAIL";
         else print (v+0 < 1e-20 && v+0 > -1e-20) ? "PASS" : "FAIL"}')
  [[ $ok == FAIL ]] && fail=1
  printf '  theta=%-4s sum_wall R = %-16s expect %-4s %s\n' "$th" "$v" "$want" "$ok" | tee -a "$LOG"
  echo "G4,theta=$th,sign_of_wall_residual,$v,$want,$ok" >> "$CSV"
done <<'SIGNS'
0   neg
60  neg
90  zero
120 pos
180 pos
SIGNS

# -------------------------------------------------------------------- G1 ----
# Build the pre-feature solver from the merge-base with main and check that a
# short run without -wall_faces is byte-identical.
say ""; say "=== G1: regression, no -wall_faces reproduces pre-feature solver ==="
BASE_REF=$(git -C "$ROOT" merge-base HEAD main 2>/dev/null || echo "")
if [[ -z "$BASE_REF" ]]; then
  say "  SKIP: no merge-base with main"
  echo "G1,regression,,,,SKIP" >> "$CSV"
else
  SNAP="$WORK/base"; mkdir -p "$SNAP"
  if git -C "$ROOT/.." archive "$BASE_REF" lunar_regolith_DSM 2>/dev/null | tar -x -C "$SNAP"; then
    if ( cd "$SNAP/lunar_regolith_DSM" && make >>"$LOG" 2>&1 ); then
      RUNOPTS=(-options_file "$ROOT/inputs/solver.opts" -periodic 0 -pf_monitor 1 -pf_output 0
               -dim 2 -Nx 24 -Ny 12 -Lx 6.25e-5 -Ly 3.125e-5 -eps 1.2e-6
               -ic_type multi_grains -ice_grain_cx 3.125e-5 -ice_grain_cy 1.5625e-5
               -ice_grain_R 2.3438e-5 -temp -20 -humidity 1.0
               -t_final 2.0e2 -delt_t 1.0e-2
               # Pin the physics on BOTH sides. This gate asks one question --
               # is the contact-angle feature inert without -wall_faces? -- and
               # a bitwise comparison answers it only if everything else is held
               # equal. -thin_iface_corr's default changed OFF -> ON on
               # 2026-09-13, so without pinning it here the gate would compare
               # two different models and fail for a reason that has nothing to
               # do with the feature under test.
               -thin_iface_corr 0)
      for which in base new; do
        folder="$WORK/g1_$which"; export folder; mkdir -p "$folder"
        bin="$EXEC"; [[ $which == base ]] && bin="$SNAP/lunar_regolith_DSM/lunar_regolith_dsm"
        "$bin" "${RUNOPTS[@]}" > "$WORK/g1_$which.out" 2>&1 || true
      done
      export folder="$WORK"
      if diff -q "$WORK/g1_base/SSA_evo.dat" "$WORK/g1_new/SSA_evo.dat" >/dev/null 2>&1; then
        say "  SSA_evo.dat identical            PASS"
        echo "G1,no_wall_faces,ssa_bitwise_identical,identical,exact,PASS" >> "$CSV"
      else
        say "  SSA_evo.dat DIFFERS              FAIL"
        diff "$WORK/g1_base/SSA_evo.dat" "$WORK/g1_new/SSA_evo.dat" | head -5 | tee -a "$LOG"
        echo "G1,no_wall_faces,ssa_bitwise_identical,differs,exact,FAIL" >> "$CSV"
        fail=1
      fi
    else
      say "  SKIP: could not build $BASE_REF"; echo "G1,regression,,,,SKIP" >> "$CSV"
    fi
  else
    say "  SKIP: could not export $BASE_REF"; echo "G1,regression,,,,SKIP" >> "$CSV"
  fi
fi

say ""
if [[ $fail -eq 0 ]]; then say "ALL GATES PASSED   (see $CSV)"; else say "FAILURES -- see $CSV"; fi
exit $fail
