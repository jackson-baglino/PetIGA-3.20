#!/usr/bin/env bash
# =============================================================================
# provenance.sh — print which commit of this repo a submit script is running.
#
# WHY. A cost estimate printed on the cluster is only as current as the
# checkout that printed it. On 2026-09-23 a replay was priced at $19 locally
# and quoted $217 on the cluster; both numbers were arithmetically right and
# the difference was entirely that the cluster had not pulled the commit which
# changed the default stride. Nothing in the output said so.
#
# Printing the commit, its date, and whether the tree is dirty or behind its
# upstream makes that visible in the same banner as the estimate.
# =============================================================================

print_repo_provenance() {
    local root="${1:-.}"
    command -v git >/dev/null 2>&1 || return 0
    git -C "$root" rev-parse --git-dir >/dev/null 2>&1 || return 0

    local sha date dirty behind upstream
    sha="$(git -C "$root" rev-parse --short HEAD 2>/dev/null)"
    date="$(git -C "$root" log -1 --format=%cd --date=format:'%Y-%m-%d %H:%M' 2>/dev/null)"
    dirty=""
    [ -n "$(git -C "$root" status --porcelain 2>/dev/null)" ] && dirty="  [uncommitted changes]"

    behind=""
    upstream="$(git -C "$root" rev-parse --abbrev-ref --symbolic-full-name '@{u}' 2>/dev/null)"
    if [ -n "$upstream" ]; then
        # Count without fetching: this reports what the last fetch knew, which
        # is the honest thing to say on a login node that may have no network
        # credentials. A non-zero count is a hard signal; zero is not a
        # guarantee the remote has not moved since.
        local n
        n="$(git -C "$root" rev-list --count "HEAD..$upstream" 2>/dev/null)"
        if [ -n "$n" ] && [ "$n" -gt 0 ]; then
            behind="  ⚠ $n commit(s) behind $upstream -- git pull before trusting this"
        fi
    fi

    echo "  repo            : $sha  ($date)$dirty$behind"
}
