#!/usr/bin/env bash
# =============================================================================
# studies/heatflux_bc/verification/verify_heatflux.sh
#
# Deterministic checks of the thermal face conditions added on branch
# heatflux-bc: -qflux_faces (prescribed heat flux) and -Twall_faces
# (prescribed wall temperature T_w(s,t), weak/Nitsche). Every case is a small
# run made through scripts/Studio/run_cryosuction.sh (one at a time, < 50k
# DOF, seconds to a minute each); check_heatflux.py compares it with a
# closed-form solution and prints PASS/FAIL rows.
#
#   J   boundary Jacobian block vs central differences (-test_wall_jacobian,
#       thermal terms switched off with the wall term to isolate the block);
#       -snes_test_jacobian on one step is reported for information
#   a   steady 1D flux: T linear with slope q/k, on two meshes
#   b   transient flux into a slab: T(0,t) - T0 = (2q/k) sqrt(alpha t/pi),
#       air and ice
#   c   erfc wall temperature on both y walls: interior follows the 1D erfc,
#       two meshes, observed order of convergence
#   c2  linear wall temperature on all faces: steady state exact
#   d   energy budget: stored heat change = boundary heat input
#   W1  bitwise: with -wall_faces (contact angle), runs that do not name a
#       thermal face reproduce the merge-base with main exactly
#   G1  is a separate gate: studies/contact_angle/verification/verify_wall_bc.sh
#
# The run script builds the solver from $PETIGA_DIR/projects/lunar_cryosuction,
# so PETIGA_DIR must point at the tree that holds this script.
# Python: $PYTHON, or venv_cryosuction/bin/python of this project.
#
# Writes verify_heatflux.csv and verify_heatflux.log next to this script.
# Exits non-zero if any check fails.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
RUN="$ROOT/scripts/Studio/run_cryosuction.sh"
CHECK="$HERE/check_heatflux.py"
CSV="$HERE/verify_heatflux.csv"
LOG="$HERE/verify_heatflux.log"
PY="${PYTHON:-$ROOT/venv_cryosuction/bin/python}"

if [[ "$(cd "${PETIGA_DIR:-/nonexistent}/projects/lunar_cryosuction" 2>/dev/null && pwd)" != "$ROOT" ]]; then
  echo "PETIGA_DIR/projects/lunar_cryosuction is not $ROOT; the run script would build another tree." >&2
  exit 2
fi
[[ -x "$PY" ]] || { echo "no python at $PY (set PYTHON)" >&2; exit 2; }

: > "$LOG"
echo "test,case,quantity,measured,expected,tolerance,result" > "$CSV"
fail=0
say() { printf '%s\n' "$*" | tee -a "$LOG"; }

# run_case <geometry> <experiment> <tag> [extra solver flags...]
# Runs one simulation through the run script; sets $RUNDIR.
run_case() {
  local g="$1" e="$2" tag="$3"; shift 3
  local out
  out=$("$RUN" "$g" "$e" "$tag" -- "$@" 2>&1)
  local rc=$?
  # The full run-script output is in the run folder (outp.txt); keep the log short.
  printf '%s\n' "$out" | grep -E "Results:|rror|FAIL" >> "$LOG"
  RUNDIR=$(printf '%s\n' "$out" | awk '/^  Results: /{print $2}' | tail -1)
  if [[ $rc -ne 0 || -z "$RUNDIR" ]]; then
    say "  run $g / $e FAILED (exit $rc)"; fail=1; return 1
  fi
  say "  run $g / $e -> $RUNDIR"
}

check() {  # check <test> --tol X [--tol-flux Y] <run dirs...>
  local out rc
  out=$("$PY" "$CHECK" "$@" 2>>"$LOG"); rc=$?
  printf '%s\n' "$out" | tee -a "$LOG" "$CSV"
  [[ $rc -ne 0 ]] && fail=1
  return 0
}

# ------------------------------------------------------------------- W1 ---
# G1 (verify_wall_bc.sh) checks the new code is bit-for-bit inert WITHOUT
# -wall_faces. The thermal terms share the boundary branch with the
# contact-angle term, which was restructured, so W1 repeats the bitwise check
# WITH -wall_faces (contact-angle channel, theta 60): the merge-base with main
# is exported, run through its own run script, and SSA_evo.dat compared.
say ""; say "=== W1: bitwise inert with -wall_faces (base = merge-base with main) ==="
W1_ARGS=(-t_final 2.0e2 -delt_t 1.0e-2 -n_out 2)
BASE_REF=$(git -C "$ROOT" merge-base HEAD main 2>/dev/null || echo "")
if [[ -z "$BASE_REF" ]]; then
  say "  SKIP: no merge-base with main"; echo "W1,wall_faces,ssa_bitwise,,,,SKIP" >> "$CSV"
else
  BT="$(mktemp -d "${TMPDIR:-/tmp}/heatflux_base.XXXXXX")"
  REPO="$(git -C "$ROOT" rev-parse --show-toplevel)"
  if git -C "$REPO" archive "$BASE_REF" include lib projects/lunar_cryosuction | tar -x -C "$BT" \
     && ln -s "$PETIGA_DIR/$PETSC_ARCH" "$BT/$PETSC_ARCH"; then
    base_out=$(PETIGA_DIR="$BT" "$BT/projects/lunar_cryosuction/scripts/Studio/run_cryosuction.sh" \
               channel_2D_H31um_eps0.86um droplet_T-20_theta60 w1_base -- "${W1_ARGS[@]}" 2>&1)
    printf '%s\n' "$base_out" | grep -E "Results:|rror|FAIL" >> "$LOG"
    BASE_RUN=$(printf '%s\n' "$base_out" | awk '/^  Results: /{print $2}' | tail -1)
    run_case channel_2D_H31um_eps0.86um droplet_T-20_theta60 w1_new "${W1_ARGS[@]}"
    if [[ -n "$BASE_RUN" && -n "${RUNDIR:-}" ]] && cmp -s "$BASE_RUN/SSA_evo.dat" "$RUNDIR/SSA_evo.dat"; then
      n=$(wc -l < "$RUNDIR/SSA_evo.dat" | tr -d ' ')
      say "  SSA_evo.dat identical ($n rows, base $BASE_REF)   PASS"
      echo "W1,wall_faces_theta60_base_${BASE_REF:0:8},ssa_bitwise_identical,$n rows,identical,exact,PASS" >> "$CSV"
    else
      say "  SSA_evo.dat DIFFERS or a run failed   FAIL"; fail=1
      echo "W1,wall_faces_theta60_base_${BASE_REF:0:8},ssa_bitwise_identical,differs,identical,exact,FAIL" >> "$CSV"
    fi
  else
    say "  SKIP: could not export $BASE_REF"; echo "W1,wall_faces,ssa_bitwise,,,,SKIP" >> "$CSV"
  fi
  rm -rf "$BT"
fi

# ------------------------------------------------------------------- J ----
# Tolerance 1e-6: the isolated-block error is bounded by central-difference
# truncation at h = 1e-7 (same as verify_wall_bc.sh G2); a wrong coefficient
# shows up at O(1).
say ""; say "=== J: thermal boundary Jacobian block ==="
for variant in "wall+thermal::" "thermal_only::-contact_angle_deg 90"; do
  name="${variant%%::*}"; extra="${variant#*::}"
  # shellcheck disable=SC2086
  run_case hf_jac_box62x31um_disc hf_jacobian "jac_${name/+/_}" -test_wall_jacobian $extra || continue
  worst=$(awk '/^   worst/{print $3}' "$RUNDIR/outp.txt" | tail -1)
  scale=$(awk '/^     4 /{print $4}' "$RUNDIR/outp.txt" | tail -1)
  ok=$(awk -v w="$worst" -v s="$scale" 'BEGIN{print (w!="" && w+0<1e-6 && s+0>0)?"PASS":"FAIL"}')
  [[ $ok == FAIL ]] && fail=1
  say "  $name  rel.err(boundary block) = $worst  ||Jbnd v|| = $scale  $ok"
  echo "J,$name,rel_err_boundary_block,$worst,0,1.0e-06,$ok" >> "$CSV"
done
run_case hf_jac_box62x31um_disc hf_jacobian sntj -t_final 1.0e-3 -delt_t 1.0e-3 -ts_max_steps 1 -snes_test_jacobian
if [[ -n "${RUNDIR:-}" ]]; then
  r=$(awk -F'= ' '/\|\|J - Jfd\|\|_F\/\|\|J\|\|_F =/{split($2,a,","); if (a[1]+0>m) m=a[1]+0} END{print m}' "$RUNDIR/outp.txt")
  say "  -snes_test_jacobian worst ||J-Jfd||/||J|| over the Newton iterations = $r (information: dominated by the interior form's FD noise)"
  echo "J,snes_test_jacobian,max_rel_FD_diff_full_system,$r,0,,INFO" >> "$CSV"
fi

# ------------------------------------------------------------------- a ----
# Tolerance 1e-8 relative to the rise qL/k: the linear steady state lies in
# the spline space, so the only errors are the SNES tolerance and the decayed
# transient (exp(-pi^2 alpha t/(4L^2)) < 1e-30 at t_final).
say ""; say "=== a: steady 1D flux ==="
A=()
run_case hf_air_box200x20um_N40 hf_steady_q100 verify && A+=("$RUNDIR")
run_case hf_air_box200x20um_N80 hf_steady_q100 verify && A+=("$RUNDIR")
[[ ${#A[@]} -eq 2 ]] && check a --tol 1e-8 --tol-flux 1e-8 "${A[@]}"

# ------------------------------------------------------------------- b ----
# Tolerance 1e-3 of the surface rise, on the surface value and on the whole
# profile (far face by images). A-priori size: spatial (h/sqrt(alpha t))^2 ~
# (2 um / 60 um)^2 ~ 1e-3 at the earliest time checked, smaller later; the
# time error of the flux switching on at t = 0 decays like dt/t.
say ""; say "=== b: transient flux into a slab ==="
B=()
run_case hf_air_box600x20um_N300 hf_slab_q100_air verify && B+=("$RUNDIR")
run_case hf_ice_box600x20um_N300 hf_slab_q1000_ice verify && B+=("$RUNDIR")
[[ ${#B[@]} -gt 0 ]] && check b --tol 1e-3 "${B[@]}"

# ------------------------------------------------------------------- c ----
say ""; say "=== c: erfc wall temperature (Nitsche) ==="
C=()
run_case hf_air_box600x100um_N75x12  hf_twall_erfc_air verify && C+=("$RUNDIR")
run_case hf_air_box600x100um_N150x25 hf_twall_erfc_air verify && C+=("$RUNDIR")
run_case hf_air_box600x100um_N75x12  hf_twall_erfc_air verify_dt4 -delt_t 2.5e-8 -dtmax 5.0e-7 && C+=("$RUNDIR")
# Tolerance 1e-3 of dT at every time on every run (the worst, at t = t_final/64
# on the coarse mesh, is 4e-4). Order h -> h/2 must exceed 1.5 (p = 2: 2
# expected in this norm) at the times where the dt/4 run shows the error is
# spatial (changes < 5 %); earlier, the t = 0 jump's time error dominates.
[[ ${#C[@]} -eq 3 ]] && check c --tol 1e-3 "${C[@]}"
run_case hf_air_box200x100um_N20x10 hf_twall_linear_air verify && check c2 --tol 1e-8 "$RUNDIR"

# ------------------------------------------------------------------- d ----
# Flux runs: the input xi_T q Ly t is exact and the stored heat is integrated
# exactly from the control values. The only error is a fixed offset from the
# integrator start (the first step stores 1.2 q Ly dt0, reported as INFO);
# tolerance 1e-4 of the final input (measured 1e-5).
# Wall-temperature run (fine mesh): the input is the time integral (trapezoid
# over every step) of -k dT/dn on the wall faces from bflux.dat; tolerance
# 2e-3 (the Nitsche face also exchanges beta (T - T_w), which bflux omits, and
# the trapezoid rule on a 1/sqrt(t) flux).
[[ ${#B[@]} -gt 0 ]] && check d --tol 1e-4 "${B[@]}"
[[ ${#C[@]} -ge 2 ]] && check d --tol 2e-3 "${C[1]}"

say ""
if [[ $fail -eq 0 ]]; then say "ALL CHECKS PASSED   (see $CSV)"; else say "FAILURES -- see $CSV"; fi
exit $fail
