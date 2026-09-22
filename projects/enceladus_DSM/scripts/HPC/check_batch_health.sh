#!/usr/bin/env bash
# =============================================================================
# check_batch_health.sh — did these runs actually run cleanly?
#
#   ./scripts/HPC/check_batch_health.sh <dir> [--t-final <s>]
#
# COMPLETION IS THE PRIMARY SIGNAL, log greps are secondary. A run either
# produced a full trajectory or it did not, and that is checkable directly
# from SSA_evo.dat and k_eff.csv. Log strings are a supporting hint, and a
# noisy one -- see below.
#
# THE FIRST VERSION OF THIS SCRIPT FLAGGED EVERY RUN, including four that had
# already been analysed in detail and were known good. Three separate bugs,
# all worth stating because they are easy to repeat:
#
#   1. It searched `--include='*.o*'` for SLURM stdout. That glob also matches
#      *.opts -- and inputs/solver.opts contains the word DIVERGED four times
#      in a comment about SNES convergence. Every run stages a copy, so every
#      run "failed". The SLURM pattern is *.o<digits>, not *.o*.
#   2. It grepped for DIVERGED at all. Rejected steps are NORMAL here: the
#      interface-CFL limiter is designed to reject and retry, so transient
#      divergence is the machinery working. Only terminal failures count.
#   3. It looked for SSA_evo.dat directly inside each subdirectory, but runs
#      live one level deeper (<geom>/<timestamp>/ on scratch), so every run
#      reported steps=0 -- including runs with 6647 steps.
#
# It also reported stale failures: a run directory accumulates the .e/.o of
# every attempt, so a Bus error from a previous failed submission is still
# there after a successful rerun. Only the NEWEST log pair is read.
#
# What is treated as a hard failure:
#   Stale file handle / make: ***   two jobs compiling into the same obj/.
#                                   run_enceladus.sh:162-165 documents it; 2+
#                                   jobs must go through submit_batch.sh.
#   Bus error / core dumped         the same race from the other side -- the
#                                   executable relinked under mapped ranks.
#   [ABORT]                         the solver's own terminal abort.
#
# BUT A LOG HIT ON A RUN THAT FINISHED IS A STALE LOG, NOT A FAILURE. A run
# directory keeps the .e/.o of every attempt, and "newest by mtime" does not
# identify the successful one: a clean .e keeps its creation time while a
# crashing job writes to its .e at crash time, so the FAILED log is often the
# newer. Taking log hits at face value flagged four L/R_ave = 64 runs that had
# each reached 100% of t_final with 52 k_eff samples. So completion wins: if
# the run finished, a log hit is reported as a stale note and nothing else.
#
# t_final IS READ FROM THE RUN'S OWN .opts, not passed in. One --t-final for a
# whole scratch directory is wrong the moment it holds more than one campaign:
# it reported the 15-hour molaro runs as "0% of t_final" against a 30-day
# target and flagged them as truncated.
#
# k_eff IS ONLY EXPECTED IF IT WAS ASKED FOR. The curvature-calibration and
# molaro runs do not pass -keff, so k_eff.csv is absent by design; requiring
# it flagged every one of them.
# =============================================================================
set -uo pipefail
root="${1:?usage: $0 <dir> [--t-final <s>]}"; shift || true
tfinal=""
while [ $# -gt 0 ]; do
    case "$1" in
        --t-final) tfinal="$2"; shift 2 ;;
        *) shift ;;
    esac
done

bad=0; n=0
# Discover by CONTENT: a directory holding SSA_evo.dat is a run, wherever it
# sits. Name-based discovery has broken on every layout change so far.
while IFS= read -r ssa; do
    d=$(dirname "$ssa"); n=$((n+1))
    name=$(realpath --relative-to="$root" "$d" 2>/dev/null || basename "$d")

    steps=$(wc -l < "$ssa" | tr -d ' ')
    keff="-"; keff_asked=0
    if [ -f "$d/k_eff.csv" ]; then
        keff_asked=1; keff=$(( $(wc -l < "$d/k_eff.csv") - 1 ))
    fi
    tend=$(awk 'END{printf "%.4g", $3}' "$ssa" 2>/dev/null)

    # t_final from THIS run's staged opts; --t-final is only a fallback
    tf=$(grep -hE '^-t_final' "$d"/*.opts "$d"/inputs/*.opts 2>/dev/null \
         | awk '{print $2}' | tail -1)
    [ -z "$tf" ] && tf="$tfinal"
    frac=""; pct=0
    if [ -n "$tf" ]; then
        pct=$(awk -v a="$tend" -v b="$tf" 'BEGIN{printf "%.0f", 100*a/b}')
        frac=$(printf " (%s%% of t_final)" "$pct")
    fi

    # newest SLURM log pair only -- a run dir keeps every attempt's .e/.o
    newest=$(ls -t "$d"/*.e[0-9]* "$d"/*.o[0-9]* 2>/dev/null | head -2)
    hits=""
    if [ -n "$newest" ]; then
        for pat in "Stale file handle" "Bus error" "core dumped" "make: ***" "[ABORT]"; do
            if grep -qlF "$pat" $newest 2>/dev/null; then hits="${hits}${pat}; "; fi
        done
    fi
    # outp.txt is the run's own record, and only terminal markers count there
    [ -f "$d/outp.txt" ] && grep -qF "[ABORT]" "$d/outp.txt" 2>/dev/null && hits="${hits}ABORT in outp; "

    # Completion is the verdict. A log hit on a finished run is a stale note.
    done_ok=0
    [ -n "$tf" ] && [ "$pct" -ge 99 ] && done_ok=1
    incomplete=0
    [ "$steps" -lt 2 ] && incomplete=1
    [ "$keff_asked" -eq 1 ] && [ "$keff" -lt 1 ] && incomplete=1
    [ -n "$tf" ] && [ "$done_ok" -eq 0 ] && incomplete=1

    if [ "$incomplete" -eq 1 ]; then
        printf "  FLAG  %-56s steps=%-6s keff=%-4s t=%s%s  %s\n" \
               "${name:0:56}" "$steps" "$keff" "$tend" "$frac" "$hits"
        bad=$((bad+1))
    elif [ -n "$hits" ]; then
        printf "  ok*   %-56s steps=%-6s keff=%-4s t=%s%s  [stale log: %s]\n" \
               "${name:0:56}" "$steps" "$keff" "$tend" "$frac" "$hits"
    else
        printf "  ok    %-56s steps=%-6s keff=%-4s t=%s%s\n" \
               "${name:0:56}" "$steps" "$keff" "$tend" "$frac"
    fi
done < <(find "$root" -name SSA_evo.dat 2>/dev/null | sort)

echo ""
[ "$n" -eq 0 ] && { echo "  no runs found under $root (looked for SSA_evo.dat)"; exit 0; }
[ "$bad" -eq 0 ] && echo "  $n run(s), all complete  (ok* = finished, but an older attempt's log is still present)" \
                 || echo "  $bad of $n run(s) INCOMPLETE -- rerun those through submit_batch.sh"
exit 0
