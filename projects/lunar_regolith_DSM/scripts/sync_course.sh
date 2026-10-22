#!/bin/bash
# Mirror projects/lunar_regolith_DSM to the course repository (remote "course",
# citgitfr/BaglinoJackson_IceSedimentInteractions): main and every feat/* branch
# that contains this project, as subtree-split history. Pushes are fast-forward
# only; archive/* tags and branches are never pushed.
#   Usage (from anywhere in the repo):  scripts/sync_course.sh [--dry-run]
set -euo pipefail
PREFIX=projects/lunar_regolith_DSM
DRY=${1:-}
cd "$(git rev-parse --show-toplevel)"
git remote get-url course >/dev/null 2>&1 || {
    echo "no 'course' remote; add it with:"
    echo "  git remote add course https://github.com/citgitfr/BaglinoJackson_IceSedimentInteractions.git"
    exit 1; }

push() {  # <local branch> <course branch>
    local sha
    sha=$(git subtree split --prefix="$PREFIX" "$1" 2>/dev/null | tail -1)
    [ -n "$sha" ] || { echo "skip $1: no history under $PREFIX"; return; }
    echo "$1 -> course/$2  ($(git log -1 --format='%h %s' "$sha"))"
    [ "$DRY" = "--dry-run" ] || git push course "${sha}:refs/heads/$2"
}

push main main
for b in $(git for-each-ref --format='%(refname:short)' 'refs/heads/feat/*'); do
    if git cat-file -e "$b:$PREFIX" 2>/dev/null; then push "$b" "$b"
    else echo "skip $b: does not contain $PREFIX (push it by hand if wanted)"; fi
done
