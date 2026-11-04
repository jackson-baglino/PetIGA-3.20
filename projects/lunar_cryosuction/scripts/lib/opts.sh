#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# scripts/lib/opts.sh — single source of truth for .opts name resolution.
#
# inputs/geometry/ and inputs/experiment/ are organised into family
# sub-directories (geometry/wedge/, experiment/tgrad/, ...), so the bare name a
# user types is NOT a path. Every script that turns a name into a file needs
# the same search, and the copies drifted: on 2026-08-25 the two RUN scripts
# had it while both BATCH drivers and the HPC runner did not, so a whole
# four-arm batch was skipped with "geometry file not found" pointing at the
# flat path -- even though the runner it would have called resolves that name
# fine. Same failure mode alloc.sh was created to end. Change it HERE only.
# ---------------------------------------------------------------------------

# resolve_opts <dir> <name>  ->  echoes the path, or returns 1
#
# Exact match at the top level wins, so a name can always be pinned by placing
# the file directly in <dir>. Otherwise search the tree. More than one hit is a
# hard error rather than a silent pick -- two families holding the same name
# means the name is not a unique handle and the caller must disambiguate.
resolve_opts() {
    local dir="$1" name="$2" hit
    if [[ -f "$dir/${name}.opts" ]]; then printf '%s' "$dir/${name}.opts"; return 0; fi
    hit=$(find "$dir" -type f -name "${name}.opts" -print 2>/dev/null)
    if [[ "$(printf '%s' "$hit" | grep -c .)" -gt 1 ]]; then
        echo "❌ Ambiguous name '${name}':" >&2
        printf '%s\n' "$hit" | sed 's|^|     |' >&2
        return 1
    fi
    [[ -n "$hit" ]] && printf '%s' "$hit" && return 0
    return 1
}

# list_opts <dir>  ->  every available name, relative and without .opts
# Use in "not found" messages so the user sees what IS available.
list_opts() {
    local dir="$1"
    [[ -d "$dir" ]] || return 0
    (cd "$dir" && find . -type f -name '*.opts' | sed 's|^\./||; s/\.opts$//' | sort)
}
