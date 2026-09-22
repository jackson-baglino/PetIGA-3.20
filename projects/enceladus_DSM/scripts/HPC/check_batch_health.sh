#!/usr/bin/env bash
# =============================================================================
# check_batch_health.sh — did these runs actually run cleanly?
#
# Run on the HPC, over a directory of run folders. Greps the SLURM .e/.o files
# and outp.txt for the failures that are silent in the results:
#
#   Stale file handle / make Error   two jobs compiling into the same obj/.
#                                    run_enceladus.sh:162-165 documents it. A
#                                    job can survive this and still have linked
#                                    against a half-written object.
#   Bus error / core dumped          the same race seen from the other side --
#                                    the executable relinked underneath ranks
#                                    that had it mapped.
#   DIVERGED / SNES / PETSC ERROR    solver failure.
#   PHASE GUARD TRIPPED              phi left [0,1] by more than the guard band;
#                                    numerics, not physics.
#
# A run that hit the compile race is NOT necessarily wrong, but it is not
# trustworthy either, and the only safe response is to rerun it from a build
# that was not being written at the time.
#
#   ./scripts/HPC/check_batch_health.sh /resnick/scratch/$USER/enceladus_DSM
# =============================================================================
set -uo pipefail
root="${1:?usage: $0 <dir of run folders>}"
bad=0
for d in "$root"/*/; do
    [ -d "$d" ] || continue
    name=$(basename "$d")
    hits=""
    for pat in "Stale file handle" "Bus error" "core dumped" "Error 1" \
               "DIVERGED" "PETSC ERROR" "PHASE GUARD TRIPPED"; do
        n=$(grep -rslF "$pat" "$d" --include='*.e*' --include='*.o*' \
              --include='outp.txt' 2>/dev/null | wc -l | tr -d ' ')
        [ "$n" -gt 0 ] && hits="${hits}${pat} (${n}); "
    done
    steps=$( [ -f "$d/SSA_evo.dat" ] && wc -l < "$d/SSA_evo.dat" | tr -d ' ' || echo 0 )
    keff=$( [ -f "$d/k_eff.csv" ] && echo $(( $(wc -l < "$d/k_eff.csv") - 1 )) || echo 0 )
    if [ -n "$hits" ]; then
        printf "  FLAG  %-58s steps=%-6s keff=%-4s %s\n" "${name:0:58}" "$steps" "$keff" "$hits"
        bad=$((bad+1))
    else
        printf "  ok    %-58s steps=%-6s keff=%s\n" "${name:0:58}" "$steps" "$keff"
    fi
done
echo ""
[ "$bad" -eq 0 ] && echo "  no failures found" \
                 || echo "  $bad run(s) flagged -- rerun them through submit_batch.sh"
exit 0
