#!/usr/bin/env bash
# =============================================================================
# studies/contact_angle/verification/dtmax_study.sh
#
# How large can -dtmax be before the ANSWER changes?
#
# There is no reliable closed-form estimate for dtmax here, and the tau_sub/10
# rule this replaces was calibrated on three points whose failure turned out to
# have a different cause (the wall-term sign flip, since clamped). The solver is
# implicit, so there is no explicit stability limit; what actually constrains dt
# is temporal accuracy, and the honest way to find that is to vary dtmax and
# watch the measured quantity.
#
# Note the code ALREADY has the right dynamic guard: InterfaceCFLMonitor caps dt
# so no DOF changes by more than -dtCFL_dphimax (0.2) per step, measured from
# the previous accepted step. It never fired in batch 2026-09-12, meaning dt was
# nowhere near that limit. dtmax should therefore be a loose ceiling, with the
# CFL limiter disposing -- this study finds where "loose" stops being safe.
#
# Runs on the cheapest geometry (eps = 3.00 um, 142x48) so the whole sweep is
# minutes, and reports theta_inf, whether the CFL limiter fired, and the phi
# bounds for each dtmax.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
GEOM=channel_2D_H100um_eps3.00um
EXP=relax_T-20_theta60
CSV="$HERE/dtmax_study.csv"
PY="$ROOT/venv_lunar/bin/python"

echo "dtmax,dtmax_over_tau,steps,cfl_caps,phi_min,phi_max,theta_inf,theta_err,wall_s" > "$CSV"
echo
printf "%-10s %-11s %-8s %-9s %-12s %-11s %-9s\n" dtmax dtmax/tau steps cfl_caps phi_min theta_inf wall[s]

for DT in 2.0e3 4.0e3 8.0e3 1.6e4 3.2e4; do
  TAG="dtmax${DT}"
  S=$(date +%s)
  "$ROOT/scripts/Studio/run_lunar.sh" "$GEOM" "$EXP" "$TAG" -- -dtmax "$DT" \
      > "$HERE/dtmax_${DT}.log" 2>&1
  E=$(date +%s)
  D=$(ls -td "$HOME/SimulationResults/lunar_regolith_DSM/scratch/$GEOM"/*"$TAG" 2>/dev/null | head -1)
  [[ -z "$D" ]] && { echo "  $DT  FAILED (no run dir)"; continue; }

  TAU=$(grep -oE 'tau_sub  =  [0-9.e+-]+' "$D/outp.txt" | head -1 | awk '{print $3}')
  RAT=$(awk -v a="$DT" -v b="$TAU" 'BEGIN{printf "%.4f", a/b}')
  CFL=$(grep -c "Interface-CFL cap" "$D/outp.txt" 2>/dev/null | tr -d "\n" || echo 0)
  STEPS=$(tail -1 "$D/SSA_evo.dat" 2>/dev/null | awk '{print $4}')

  read -r PMIN PMAX TH <<< "$("$PY" - "$D" <<'PYEOF'
import sys,os,glob,numpy as np,csv,re
from igakit.io import PetIGA
D=sys.argv[1]; io=PetIGA(); nrb=io.read(os.path.join(D,"igasol.dat"))
lo,hi=0.0,1.0
for f in sorted(glob.glob(D+"/sol_*.dat")):
    p=io.read_vec(f,nrb)[...,0]; lo=min(lo,p.min()); hi=max(hi,p.max())
th="nan"
c=os.path.join(D,"contact_angle.csv")
if os.path.exists(c):
    m=re.search(r'theta_inf = ([-\d.]+)', open(c).read())
    if m: th=m.group(1)
print(f"{lo:.3e} {hi:.6f} {th}")
PYEOF
)"
  ERR=$(awk -v t="$TH" 'BEGIN{ if (t=="nan") print "nan"; else printf "%+.3f", t-60.0 }')
  printf "%-10s %-11s %-8s %-9s %-12s %-11s %-9s\n" "$DT" "$RAT" "$STEPS" "$CFL" "$PMIN" "$TH" "$((E-S))"
  echo "$DT,$RAT,$STEPS,$CFL,$PMIN,$PMAX,$TH,$ERR,$((E-S))" >> "$CSV"
done
echo
echo "-> $CSV"
