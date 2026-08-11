#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# scripts/lib/opts.sh — single source of truth for resolving a bare .opts name
# to a path.
#
# Geometry and experiment .opts live in per-family / per-campaign
# subdirectories (inputs/geometry/sinter/, inputs/experiment/molaro/, ...), but
# every script's CLI takes a BARE NAME. Resolving one to the other is therefore
# something six scripts need, and until 2026-08-10 only two of them had it:
# submit_enceladus.sh and Studio/run_enceladus.sh carried private copies of
# resolve_opts, while submit_batch.sh, HPC/run_enceladus.sh and
# run_batch_tests.sh still did a flat "$dir/$name.opts".
#
# The failure mode was quiet and expensive. submit_batch.sh reported
# "geometry file not found" for every job in a batch and exited 0 with
# "Submitted: 0" -- and because submit_enceladus.sh resolved the path only for
# its OWN allocation sizing and then handed the bare name to
# HPC/run_enceladus.sh, single submits of a subdirectory geometry would have
# died inside the SLURM job instead, after queueing.
#
# Same lesson as alloc.sh: put it HERE and nowhere else.
# ---------------------------------------------------------------------------

# resolve_opts <dir> <bare-name>
#   Prints the full path to <dir>/**/<bare-name>.opts on stdout, returns 0.
#   Returns 1 if there is no match. Returns 1 and prints the candidates to
#   stderr if the name is ambiguous -- guessing between two files that differ
#   only by campaign is worse than stopping.
#   A top-level hit always wins, so an explicit "subdir/name" also works.
resolve_opts() {
    local dir="$1" name="$2" hit
    if [ -f "$dir/${name}.opts" ]; then printf '%s' "$dir/${name}.opts"; return 0; fi
    hit=$(find "$dir" -type f -name "${name}.opts" -print 2>/dev/null)
    if [ "$(printf '%s' "$hit" | grep -c .)" -gt 1 ]; then
        echo "❌ Ambiguous name '${name}':" >&2
        printf '%s\n' "$hit" | sed 's|^|     |' >&2
        return 1
    fi
    [ -n "$hit" ] && printf '%s' "$hit" && return 0
    return 1
}

# require_opts <dir> <bare-name> <label>
#   resolve_opts, but on failure prints a diagnostic naming what was searched
#   and what exists, then returns 1. Callers that must abort should do so
#   themselves so batch drivers can skip one spec and continue.
require_opts() {
    local dir="$1" name="$2" label="$3" path
    if path="$(resolve_opts "$dir" "$name")"; then printf '%s' "$path"; return 0; fi
    {
        echo "❌ ${label} '${name}' not found under ${dir}/"
        echo "   Searched recursively for ${name}.opts. Available ${label} names:"
        find "$dir" -type f -name '*.opts' -exec basename {} .opts \; 2>/dev/null \
            | sort | sed 's|^|     |' | head -40
    } >&2
    return 1
}
