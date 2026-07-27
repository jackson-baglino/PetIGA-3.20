#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# check_ic_types.sh — guard against dead -ic_type / mesh-sizing regressions.
#
# The C solver (src/permafrost2.c) dispatches a fixed set of -ic_type values
# and SETERRQs on anything else in 2D/3D. Config files drift out of sync with
# that list every time the model changes (the 2026-06-13 two-phase fork left
# 21 .opts files that aborted at startup and silently broke the regression
# suite). This script catches that class of failure before a run is launched.
#
# Checks, per inputs/**/*.opts:
#   1. Every 2D/3D file's -ic_type is in the live dispatch list.
#      (1D falls through to a default, so any value is accepted there.)
#   2. Any file setting -geom_file also carries a "# DOF_GRID: nx ny [nz]"
#      comment — run_permafrost.sh parses it for rank sizing; without it the
#      run silently drops to 1 rank.
#
# Exit 0 if clean, 1 if any problem is found. Zero external dependencies.
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUTS="$ROOT/inputs"
SRC="$ROOT/src/permafrost2.c"

# Derive the valid 2D/3D ic_type list directly from the source dispatch, so
# this script never goes stale relative to the code it guards.
VALID=$(grep -oE 'strcmp\(ic_type, "[a-z_0-9]+"\)' "$SRC" \
        | sed -E 's/.*"([a-z_0-9]+)".*/\1/' | sort -u)
if [[ -z "$VALID" ]]; then
    echo "check_ic_types: could not extract valid ic_type list from $SRC" >&2
    exit 1
fi

fail=0
while IFS= read -r f; do
    dim=$(grep -oE '^-dim[[:space:]]+[0-9]' "$f" | grep -oE '[0-9]' || true)
    ic=$(grep -oE '^-ic_type[[:space:]]+[A-Za-z_0-9]+' "$f" | awk '{print $2}' || true)

    # ic_type validity (2D/3D only; 1D has a fall-through default)
    if [[ -n "$ic" && "$dim" != "1" ]]; then
        if ! grep -qx "$ic" <<< "$VALID"; then
            echo "DEAD ic_type '$ic' (dim=${dim:-?}): ${f#"$ROOT"/}" >&2
            fail=1
        fi
    fi

    # DOF_GRID comment required whenever -geom_file is set
    if grep -qE '^-geom_file[[:space:]]' "$f"; then
        if ! grep -qE '^#[[:space:]]*DOF_GRID:[[:space:]]*[0-9]' "$f"; then
            echo "MISSING '# DOF_GRID: nx ny [nz]' (uses -geom_file): ${f#"$ROOT"/}" >&2
            fail=1
        fi
    fi
done < <(find "$INPUTS" -name '*.opts' | sort)

if [[ $fail -eq 0 ]]; then
    echo "check_ic_types: OK — all .opts pass (valid ic_types: $(echo $VALID | tr '\n' ' '))"
fi
exit $fail
