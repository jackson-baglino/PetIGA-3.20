#!/usr/bin/env bash
# =============================================================================
# select_snapshots.sh — mark the four snapshots worth downloading from a run.
#
#   bash scripts/lib/select_snapshots.sh <run_dir> [<run_dir> ...]
#
# Writes <run_dir>/.rsync-snapshots, an rsync filter file listing
#   the opening frame   first sol_*.dat with 1 s <= t <= 1 h (pplib's t = 0)
#   t_final / 3         the snapshot nearest one third of the run
#   2 t_final / 3       the snapshot nearest two thirds
#   the last snapshot
# plus igasol.dat. scripts/HPC/fetch_stage.sh reads it through rsync's
# per-directory merge (dir-merge), so a stage download brings exactly these
# four files per run (~0.8 GB at L/R 40) instead of all ~44 (~9 GB), and
# postprocess/render_snapshots.py turns them into PNGs locally.
#
# Pure bash + awk on SSA_evo.dat (time = column 3, step = column 4): it runs
# at the end of every job (scripts/HPC/run_enceladus.sh) in about a second,
# and on the login node for runs that finished before it existed.
# Never fails the caller: a run without snapshots or SSA_evo.dat is skipped.
# =============================================================================
for run in "$@"; do
    run="${run%/}"
    ssa="$run/SSA_evo.dat"
    [[ -f "$ssa" ]] || { echo "  select_snapshots: no SSA_evo.dat in $run" >&2; continue; }
    steps=$(ls "$run" 2>/dev/null | sed -n 's/^sol_0*\([0-9][0-9]*\)\.dat$/\1/p; s/^sol_0\{5\}\.dat$/0/p' | sort -n -u)
    [[ -n "$steps" ]] || { echo "  select_snapshots: no sol_*.dat in $run" >&2; continue; }
    pick=$(awk -v snaps="$(echo $steps)" '
        BEGIN { n = split(snaps, S, " "); for (i = 1; i <= n; i++) want[S[i]] = 1 }
        { t[$4 + 0] = $3 + 0; if ($3 + 0 > tend) tend = $3 + 0 }
        END {
            t[0] = 0
            open = ""
            for (i = 1; i <= n; i++) { s = S[i]; if ((s in t) && t[s] >= 1 && t[s] <= 3600) { open = s; break } }
            if (open == "") open = S[1]
            t0 = t[open]
            printf "%s\n", open
            for (f = 1; f <= 2; f++) {
                target = t0 + f / 3 * (tend - t0); best = ""; bd = -1
                for (i = 1; i <= n; i++) { s = S[i]; if (!(s in t)) continue
                    d = t[s] - target; if (d < 0) d = -d
                    if (bd < 0 || d < bd) { bd = d; best = s } }
                if (best != "") printf "%s\n", best
            }
            printf "%s\n", S[n]
        }' "$ssa" | sort -n -u)
    {
        echo "# written by scripts/lib/select_snapshots.sh -- the snapshots fetch_stage.sh downloads"
        echo "+ igasol.dat"
        for s in $pick; do printf '+ sol_%05d.dat\n' "$s"; done
    } > "$run/.rsync-snapshots"
    echo "  select_snapshots: $(basename "$run"): steps $(echo $pick)"
done
exit 0
