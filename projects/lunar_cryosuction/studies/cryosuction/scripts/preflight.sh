#!/usr/bin/env bash
# studies/cryosuction/scripts/preflight.sh <tests-file>
#
# Parse check before an HPC submission: for every geom:exp line, resolve both
# names, then start the solver with -t_final 0 on one rank (set-up, IC and the
# header only; no time step is taken) and check that it exits cleanly and that
# the printed boundary table and tau_sub are what the files claim. This is a
# set-up check like the verification drivers, not a simulation; it writes only
# to a temporary folder. It first runs check_open_bc.sh (plan tool T5), which
# fails any pair whose open vapor-reservoir face still lets ice through
# (-phi_open_bc not 0). For a geometry read from a mesh file (-geom_file) it
# also runs preprocess/check_walls.py: walls must be smooth, with no valley or
# tip of radius below 20 eps and no corner (user, 2026-11-13).
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
source "$ROOT/scripts/lib/opts.sh"
EXEC="$ROOT/lunar_cryosuction"
[[ -x $EXEC ]] || { echo "missing $EXEC -- run make" >&2; exit 2; }
WORK="${TMPDIR:-/tmp}/cryo_preflight.$$"; mkdir -p "$WORK"; trap 'rm -rf "$WORK"' EXIT
fail=0
echo "-- open-boundary check (check_open_bc.sh)"
"$ROOT/studies/cryosuction/scripts/check_open_bc.sh" "$1" || fail=1
echo "-- solver set-up check"
while IFS= read -r line || [[ -n $line ]]; do
  line="${line%%#*}"; line="$(echo "$line" | tr -d '[:space:]')"
  [[ -z $line ]] && continue
  g="${line%%:*}"; e="${line#*:}"
  gf=$(resolve_opts "$ROOT/inputs/geometry" "$g") || { echo "FAIL $line: geometry not found"; fail=1; continue; }
  ef=$(resolve_opts "$ROOT/inputs/experiment" "$e") || { echo "FAIL $line: experiment not found"; fail=1; continue; }
  d="$WORK/$g.$e"; mkdir -p "$d"
  if ( cd "$ROOT" && folder="$d" "$EXEC" -options_file inputs/solver.opts -options_file "$gf" \
         -options_file "$ef" -t_final 0 -pf_output 0 > "$d/out.txt" 2>&1 ); then
    tau=$(grep -m1 -E "tau_sub *=" "$d/out.txt" | sed 's/^ *//')
    dtr=$(grep -m1 -E "dtmax/tau_sub" "$d/out.txt" | sed 's/^ *//')
    echo "OK   $line"
    echo "     $tau   $dtr"
    grep -E "h_[xyz] = L/N|h_[xyz] = mapped|eps / smallest|only valid for front|RESOLUTION WARNING" "$d/out.txt" | sed 's/^ */     /'
    grep -q "RESOLUTION WARNING" "$d/out.txt" && fail=1
    grep -E "^ +[xy] = " "$d/out.txt" | sed 's/^/    /'
    mesh=$(awk '$1=="-geom_file"{print $2}' "$gf" | head -n1)
    if [[ -n $mesh ]]; then
      PY="$ROOT/venv_cryosuction/bin/python"
      if [[ -x $PY ]]; then
        epsv=$(awk '$1=="-eps"{print $2}' "$gf" | head -n1)
        "$PY" "$ROOT/preprocess/check_walls.py" "$ROOT/$mesh" ${epsv:+--eps $epsv} | sed 's/^/   /' || fail=1
      else
        echo "     (wall check skipped: no venv_cryosuction python)"
      fi
    fi
  else
    echo "FAIL $line"; tail -5 "$d/out.txt" | sed 's/^/     /'; fail=1
  fi
done < "$1"
exit $fail
