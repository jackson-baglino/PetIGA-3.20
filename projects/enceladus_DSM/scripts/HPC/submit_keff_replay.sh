#!/usr/bin/env bash
# =============================================================================
# submit_keff_replay.sh — recompute k_eff on finished runs under a different
# band-interpolation law, without repeating the sintering.
#
# WHY. -keff_interp defaulted to `arith` until 2026-09-23 and defaults to
# `tensor` now, so every run before that date carries the arithmetic law's
# O(eps) normal-flux bias. The phase field on disk is unaffected -- the law
# enters only the steady-state cell problem -- so the fix is a replay, not a
# rerun: -keff_replay reads each sol_*.dat, re-solves the correctors, and
# writes k_eff_<law>.csv alongside the run. The in-line k_eff.csv is left
# alone (src/keff.c:377-381), so the old and new numbers sit side by side.
#
# WHAT IT DISCOVERS, AND WHY IT DOES NOT TAKE PATHS. A batch run writes to
#   $SCRATCH/enceladus_DSM/batch_<ts>[_<tag>]/<geom>__<exp>/
# but a run RESUMED after a walltime kill goes through submit_enceladus.sh,
# which has no BATCH_OUT_DIR and therefore lands in the single-run tree
#   $SCRATCH/enceladus_DSM/<geom>/<ts>_<exp>_resume_job<id>/
# So one logical run's snapshots can be split across two parents, and the
# merged form that postprocess/merge_restart_legs.py builds exists only after
# download. This script therefore SCANS for directories that actually contain
# sol_*.dat and submits one job per directory -- legs included. Merge the
# resulting per-leg CSVs afterwards, the same way the legs themselves are
# merged.
#
# Geometry and experiment are read from the .opts files run_enceladus.sh
# staged into each run directory, not parsed out of the directory name: the
# batch and resume directory layouts differ, the staged files do not.
#
# USAGE
#   ./scripts/HPC/submit_keff_replay.sh --roots <dir> [<dir>...] [options]
#
#   --roots <dirs>    where to scan (required). Give both the batch parent and
#                     $SCRATCH/enceladus_DSM/<geom> if the run was resumed.
#   --laws "<list>"   band interpolations to replay (default: tensor).
#                     e.g. --laws "tensor sharp" for the cross-check.
#   --match <glob>    only run directories whose name matches (default: *)
#   --target-samples N  aim for ~N samples per logical run (default 60), and
#                     derive the stride from it. Legs are strided TOGETHER from
#                     the group's total, so both halves of a resumed
#                     trajectory get the same cadence.
#   --stride N        fixed stride, overriding --target-samples.
#   --max-cost D      refuse to submit above this estimate (default 100).
#   --tag <tag>       batch tag (default: keff_replay)
#   --dry-run         print the discovered runs, the cost estimate and the
#                     tests-file, and submit nothing.
#
# YOU ALMOST NEVER WANT arith HERE. Each run's existing in-line k_eff.csv IS
# the arithmetic result -- it is what the run wrote while it was integrating.
# Replay arith only to check the replay path itself against that file.
#
# WHY THE DEFAULT IS A SAMPLE COUNT AND NOT stride 1. The snapshot count on
# $SCRATCH is not the count in a downloaded copy. The 2026-09-16 pilot held
# ~670 snapshots per seed on the cluster and ~55 after thin_snapshots.py ran
# on the download. Sizing the replay from the thinned number and running at
# stride 1 cost 12x what it should have. 60 samples over 30 days matches the
# cadence of the in-line k_eff.csv being compared against.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
source "$PROJECT_ROOT/scripts/lib/provenance.sh"
cd "$PROJECT_ROOT"

GEOMETRY_DIR="$PROJECT_ROOT/inputs/geometry"
EXPERIMENT_DIR="$PROJECT_ROOT/inputs/experiment"

# Per-sample cost model, from the 2026-09-16 pilot's own k_eff.csv wall_s
# column (36-64 s per sample at 241 ranks under arith) scaled by the ~2.5x CG
# iteration count the tensor law shows on studies/keff_sharp_limit/disk.
# Indicative only -- it is here so nobody submits without a number.
SEC_PER_SAMPLE=100
RATE_PER_CORE_HOUR=0.012

roots=(); laws="tensor"; match='*'; stride=0; target=60; max_cost=100
tag="keff_replay"; dry=0
while [ $# -gt 0 ]; do
    case "$1" in
        --roots)   shift; while [ $# -gt 0 ] && [[ "$1" != --* ]]; do roots+=("$1"); shift; done ;;
        --laws)    laws="$2"; shift 2 ;;
        --match)   match="$2"; shift 2 ;;
        --stride)  stride="$2"; shift 2 ;;
        --target-samples) target="$2"; shift 2 ;;
        --max-cost) max_cost="$2"; shift 2 ;;
        --tag)     tag="$2"; shift 2 ;;
        --dry-run) dry=1; shift ;;
        -h|--help) sed -n '2,50p' "${BASH_SOURCE[0]}"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

[ "${#roots[@]}" -gt 0 ] || { echo "❌ --roots is required" >&2; exit 2; }

# Resolve a run directory's geometry/experiment from the .opts it staged.
# Returns "<geom_stem> <exp_stem>" or nothing.
resolve_opts() {
    local run="$1" geom="" exp="" f stem
    for f in "$run"/*.opts; do
        [ -e "$f" ] || continue
        stem="$(basename "$f" .opts)"
        [ "$stem" = "solver" ] && continue
        if   find "$GEOMETRY_DIR"   -name "${stem}.opts" -print -quit | grep -q .; then geom="$stem"
        elif find "$EXPERIMENT_DIR" -name "${stem}.opts" -print -quit | grep -q .; then exp="$stem"
        fi
    done
    [ -n "$geom" ] && [ -n "$exp" ] && echo "$geom $exp"
}

TESTS_FILE="$(mktemp -t keff_replay_tests.XXXXXX)"
: > "$TESTS_FILE"

echo "=== scanning for finished runs with snapshots ==="
RUNS=(); GEOMS=(); EXPS=(); NSNAPS=(); NSAMPS=(); STRIDES=()
for root in "${roots[@]}"; do
    [ -d "$root" ] || { echo "  ⚠ not a directory, skipping: $root"; continue; }
    # maxdepth 2: catches both <batch>/<geom>__<exp>/ and <geom>/<ts>_..._resume/.
    while IFS= read -r run; do
        case "$(basename "$run")" in $match) ;; *) continue ;; esac
        nsnap=$(find "$run" -maxdepth 1 -name 'sol_*.dat' | wc -l | tr -d ' ')
        [ "$nsnap" -gt 0 ] || continue

        geom=""; exp=""
        read -r geom exp < <(resolve_opts "$run") || true
        if [ -z "$geom" ] || [ -z "$exp" ]; then
            echo "  ⚠ SKIP $(basename "$run") — could not resolve its .opts"; continue
        fi

        RUNS+=("$run"); GEOMS+=("$geom"); EXPS+=("$exp"); NSNAPS+=("$nsnap")
        echo "  $(basename "$run")"
        echo "      $nsnap snapshots"
    done < <(find "$root" -mindepth 1 -maxdepth 2 -type d | sort)
done

n_runs=${#RUNS[@]}
[ "$n_runs" -gt 0 ] || { echo "❌ no replayable runs found under the given roots" >&2; exit 1; }

# A run's output folder is <geom>__<exp>__<label>. The label is the law alone
# unless two directories share a geometry AND an experiment -- which is exactly
# what a resumed leg does -- in which case they would collide and each needs a
# discriminator too. Count first, decorate only where it is needed.
# ---------------------------------------------------------------------------
# Stride. A run's LEGS are strided together, from the group's total snapshot
# count, so both legs land on the same cadence -- leg 1 typically holds 3x the
# snapshots of its resume, and per-directory striding would sample the second
# half of every trajectory three times as densely as the first.
#
# The default targets a sample COUNT rather than taking stride 1, because the
# snapshot count on $SCRATCH is not the count in a downloaded copy: the
# downloads here had been through thin_snapshots.py and held ~55 per seed
# where the originals held ~670. Pricing the replay off the thinned number and
# running at stride 1 was a 12x cost error.
# ---------------------------------------------------------------------------
for i in "${!RUNS[@]}"; do
    if [ "$stride" -gt 0 ]; then
        STRIDES+=("$stride")
    else
        group_snaps=0
        for j in "${!RUNS[@]}"; do
            if [ "${GEOMS[$i]}" = "${GEOMS[$j]}" ] && [ "${EXPS[$i]}" = "${EXPS[$j]}" ]; then
                group_snaps=$((group_snaps + ${NSNAPS[$j]}))
            fi
        done
        st=$(( (group_snaps + target - 1) / target ))
        [ "$st" -lt 1 ] && st=1
        STRIDES+=("$st")
    fi
    NSAMPS+=($(( (${NSNAPS[$i]} + ${STRIDES[$i]} - 1) / ${STRIDES[$i]} )))
done

n_jobs=0; total_samples=0
for i in "${!RUNS[@]}"; do
    dup=0
    for j in "${!RUNS[@]}"; do
        [ "$i" = "$j" ] && continue
        [ "${GEOMS[$i]}" = "${GEOMS[$j]}" ] && [ "${EXPS[$i]}" = "${EXPS[$j]}" ] && dup=1
    done

    runtok=""
    if [ "$dup" -eq 1 ]; then
        base="$(basename "${RUNS[$i]}")"
        # The SLURM job id is in a resume directory's name and is the shortest
        # thing that is actually unique; fall back to the position.
        runtok="$(printf '%s' "$base" | grep -oE 'job[0-9]+' | tail -1 || true)"
        [ -n "$runtok" ] || runtok="r$(printf '%02d' "$i")"
    fi

    for law in $laws; do
        [ "$law" = "${laws%% *}" ] && printf "  %-58s stride %-3s -> %s samples\n" \
            "$(basename "${RUNS[$i]}")" "${STRIDES[$i]}" "${NSAMPS[$i]}"
        if [ -f "${RUNS[$i]}/k_eff_${law}.csv" ]; then
            echo "  ⚠ ${RUNS[$i]##*/}/k_eff_${law}.csv already exists — it will be overwritten"
        fi
        printf '%s:%s:--label %s -keff 1 -keff_interp %s -keff_replay %s -keff_replay_stride %s\n' \
            "${GEOMS[$i]}" "${EXPS[$i]}" "${law}${runtok:+_$runtok}" \
            "$law" "${RUNS[$i]}" "${STRIDES[$i]}" >> "$TESTS_FILE"
        n_jobs=$((n_jobs + 1))
        total_samples=$((total_samples + ${NSAMPS[$i]}))
    done
done

[ "$n_jobs" -gt 0 ] || { echo "❌ no replayable runs found under the given roots" >&2; exit 1; }

# Cost. nprocs is per-geometry, so read it back from the first spec's geometry
# rather than assuming the pilot's 241.
nprocs=$(awk -F: 'NR==1{print $1}' "$TESTS_FILE" | while read -r g; do
    f=$(find "$GEOMETRY_DIR" -name "${g}.opts" -print -quit)
    nx=$(grep -E '^[[:space:]]*-Nx[[:space:]]' "$f" | awk '{print $2}' | head -1)
    ny=$(grep -E '^[[:space:]]*-Ny[[:space:]]' "$f" | awk '{print $2}' | head -1)
    python3 -c "import math;print(math.ceil(3*${nx:-1}*${ny:-1}/100000))"
done)
core_hours=$(python3 -c "print(f'{$total_samples*$SEC_PER_SAMPLE*${nprocs:-1}/3600:.1f}')")
cost=$(python3 -c "print(f'{$total_samples*$SEC_PER_SAMPLE*${nprocs:-1}/3600*$RATE_PER_CORE_HOUR:.2f}')")

echo ""
echo "============================================================"
print_repo_provenance "$PROJECT_ROOT"
echo "  run directories : $n_runs"
echo "  laws            : $laws"
echo "  jobs            : $n_jobs"
echo "  samples total   : $total_samples"
echo "  ranks per job   : ${nprocs:-?}"
echo "  estimated       : ~${core_hours} core-hours  ~\$${cost} at \$${RATE_PER_CORE_HOUR}/core-hour"
echo "                    (indicative: ${SEC_PER_SAMPLE} s/sample, tier-1 rate)"
echo "============================================================"
echo ""

# Refuse anything above --max-cost. The estimate is printed either way, but a
# printed number does not stop a submission and this one did not: a stride-1
# run of the un-thinned pilot was priced at $435 and went through.
over=$(python3 -c "print(1 if $cost > $max_cost else 0)")
if [ "$over" -eq 1 ] && [ "$dry" -eq 0 ]; then
    rm -f "$TESTS_FILE"
    cat >&2 <<MSG
❌ estimated \$${cost} exceeds --max-cost \$${max_cost}; nothing submitted.

   Lower the sample count with --target-samples (default $target per run,
   legs included) or raise the ceiling with --max-cost if this is intended.
   The existing in-line k_eff.csv of each run is the ARITH result already --
   there is no need to replay arith to obtain it.
MSG
    exit 3
fi

if [ "$dry" -eq 1 ]; then
    echo "=== tests file (dry run — nothing submitted) ==="
    cat "$TESTS_FILE"
    rm -f "$TESTS_FILE"
    exit 0
fi

# One submit_batch.sh call: it builds once on the login node and exports
# SKIP_COMPILE=1, which is what keeps the jobs out of each other's obj/.
exec ./scripts/HPC/submit_batch.sh --tag "$tag" --tests-file "$TESTS_FILE"
