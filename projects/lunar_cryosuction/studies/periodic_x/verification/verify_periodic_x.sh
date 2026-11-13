#!/usr/bin/env bash
# =============================================================================
# studies/periodic_x/verification/verify_periodic_x.sh
#
# Gates for -periodic_x (solver 44f8e40f) and the periodic layouts of
# preprocess/build_ice_layout.py. Every simulation goes through
# scripts/Studio/run_cryosuction.sh (staged, reproducible run folders), one at
# a time, chained with &&: 9 runs of 13.7 k DOF, 2 simulated days each, a few
# minutes in all. The checks are in check_periodic_x.py (integrals at Gauss
# points from the sol_*.dat files; see its docstring).
#
#   a  Mesh/IC periodicity (no run): periodic knots = PetIGA's, wrapped control
#      columns = first ones + Lx, walls and slope continuous at the seam, phi
#      file = the N periodic nodes, seam layouts = mid layouts shifted N/2.
#   b  Flat periodic channel (P_mid): lens center drift ~ 0, water Q conserved
#      as well as in the sealed open channel (O_mid, reference).
#   c  Seam crossing: c1 P_seam vs P_mid, c2 CD_shift (throat on the seam, lens
#      across it) vs CD: area and center mod Lx equal at every snapshot.
#   d  Converging-diverging cell (CD): the lens center moves toward the throat;
#      Q conserved.
#   f  Built-in rectangle (-geom_file "" -Nx 114 -Ny 38) = the flat periodic
#      geom_file mesh (RECT_mid vs P_mid).
#   g  MPI: CD on 4 ranks (TARGET_DOFS_PER_CORE=3500) = CD on 1 rank.
#   h  CD and CD_shift with a tightened Newton solve: the water change and the
#      c2 differences shrink, so at default tolerances they are the solver
#      tolerance, not the seam.
#   e  G1 is studies/contact_angle/verification/verify_wall_bc.sh (run it
#      separately; it builds the merge-base solver).
#
# The run script takes its project root from $PETIGA_DIR, so PETIGA_DIR must
# be the checkout that holds this script.
#
#   bash studies/periodic_x/verification/verify_periodic_x.sh            # all
#   bash studies/periodic_x/verification/verify_periodic_x.sh --reuse F  # checks only,
#        on the runs listed in F (lines "label run_dir"; the last run's list is
#        kept as verify_periodic_x_runs.txt)
#
# Writes verify_periodic_x.csv and .log next to this script; exits non-zero on
# any FAIL.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
PY="$ROOT/venv_cryosuction/bin/python"
RUN="$ROOT/scripts/Studio/run_cryosuction.sh"
EXP=px_T-20_iso_theta60_2d_ac1e-2
CSV="$HERE/verify_periodic_x.csv"; LOG="$HERE/verify_periodic_x.log"
RUNS="$HERE/verify_periodic_x_runs.txt"
: > "$LOG"
say() { printf '%s\n' "$*" | tee -a "$LOG"; }

if [[ "$ROOT" != "$(cd "${PETIGA_DIR:-/nonexistent}/projects/lunar_cryosuction" 2>/dev/null && pwd)" ]]; then
  echo "PETIGA_DIR must point at the checkout holding this script ($ROOT)" >&2; exit 2
fi
cd "$ROOT" || exit 2

say "=== inputs (studies/periodic_x/scripts/make_px_inputs.py) ==="
"$PY" studies/periodic_x/scripts/make_px_inputs.py >> "$LOG" 2>&1 || { say "input build failed"; exit 2; }

echo "gate,case,quantity,measured,expected,tolerance,result" > "$CSV"
say ""; say "=== a: mesh and phi-file periodicity ==="
"$PY" "$HERE/check_periodic_x.py" mesh inputs/geometry/periodic_x | tee -a "$LOG" | grep -v '^#' >> "$CSV"
fa=${PIPESTATUS[0]}

if [[ "${1:-}" == "--reuse" ]]; then
  RUNS="$2"
else
  : > "$RUNS"
  # run_one LABEL GEOMETRY [extra solver flags...]; env TARGET_DOFS_PER_CORE passes through
  run_one() {
    local lab="$1" geo="$2"; shift 2
    local out="$HERE/.run_$lab.out"
    say "  run $lab: $geo $*"
    "$RUN" "$geo" "$EXP" "verify_$lab" -- "$@" > "$out" 2>&1 || { say "  $lab FAILED (see $out)"; return 1; }
    local dir; dir=$(sed -n 's/^  Results: //p' "$out" | tail -1)
    echo "$lab $dir" >> "$RUNS"; rm -f "$out"
  }
  TIGHT=(-snes_rtol 1e-10 -snes_atol 1e-30 -snes_stol 1e-14)
  say ""; say "=== runs (run_cryosuction.sh, one at a time) ==="
  run_one P_mid    px_flat120x40_lens_mid_eps1.5um &&
  run_one O_mid    px_open_flat120x40_lens_mid_eps1.5um &&
  run_one P_seam   px_flat120x40_lens_seam_eps1.5um &&
  run_one CD       px_cd120_W40_wt16_lens_eps1.5um &&
  run_one CD_shift px_cd120_W40_wt16_shift_lens_eps1.5um &&
  run_one RECT_mid px_flat120x40_lens_mid_eps1.5um -geom_file "" -Nx 114 -Ny 38 &&
  TARGET_DOFS_PER_CORE=3500 run_one CD_mpi px_cd120_W40_wt16_lens_eps1.5um &&
  run_one CD_tight       px_cd120_W40_wt16_lens_eps1.5um       "${TIGHT[@]}" &&
  run_one CD_shift_tight px_cd120_W40_wt16_shift_lens_eps1.5um "${TIGHT[@]}" ||
  { say "a run failed; checks skipped"; exit 1; }
fi

say ""; say "=== b, c, d, f, g, h: runs ==="
"$PY" "$HERE/check_periodic_x.py" runs "$RUNS" | tee -a "$LOG" | grep -v '^#' >> "$CSV"
fr=${PIPESTATUS[0]}

say ""
if [[ $fa -eq 0 && $fr -eq 0 ]]; then say "verify_periodic_x: ALL PASS (see $CSV)"; exit 0; fi
say "verify_periodic_x: FAIL (see $CSV)"; exit 1
