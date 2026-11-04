#!/usr/bin/env bash
# check_freeze.sh — is the paper-2 model still the frozen model?
#
# The solver here is frozen at git tag $TAG (see FREEZE.md). This script
# compares the guarded files -- src/, include/ and inputs/solver.opts -- with
# the tag and checks that every commit since the split that touched them has a
# row in the FREEZE.md ledger. It reads git only; it builds and runs nothing.
#
#   Usage (from anywhere in the repo):  scripts/check_freeze.sh
#
# Fails (exit 1) when:
#   - inputs/solver.opts differs from the tag (solver defaults are frozen);
#   - a commit after $BASELINE touches a guarded path and its short hash
#     (first 8 characters) is not in the FREEZE.md ledger;
#   - a guarded path has uncommitted changes (commit them, then log them).
# Prints, but does not fail on, the per-file differences of src/ and include/
# from the tag: allowed setup changes (new ICs, opt-in options, output, bug
# fixes) show up there once they are logged.
set -uo pipefail

TAG=lunar_cryosuction/model-freeze-v1
OLD=projects/lunar_regolith_DSM        # where the frozen tree lives in the tag
BASELINE=52e7ba16                      # last split commit (rename and rewiring)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$(git -C "$HERE" rev-parse --show-toplevel)"
NEW="${HERE#"$REPO"/}"                 # projects/lunar_cryosuction
LEDGER="$HERE/FREEZE.md"
GUARDED=("$NEW/src" "$NEW/include" "$NEW/inputs/solver.opts")
fail=0

echo "freeze tag : $TAG ($(git -C "$REPO" rev-parse --short "$TAG^{commit}"))"
echo "baseline   : $BASELINE   project: $NEW"

# --- 1. per-file comparison with the tag ------------------------------------
# The main file was renamed lunar_main.c -> cryosuction_main.c at the split.
map_new() { local f=$1; [[ $f == src/lunar_main.c ]] && f=src/cryosuction_main.c; echo "$f"; }
echo; echo "== guarded files vs the tag (working tree)"
changed=0
while IFS= read -r path; do
    rel=${path#"$OLD"/}; new=$(map_new "$rel")
    if [[ ! -f "$HERE/$new" ]]; then echo "  removed   $new"; changed=1
    elif ! git -C "$REPO" show "$TAG:$path" | cmp -s - "$HERE/$new"; then
        echo "  modified  $new"; changed=1
        if [[ $new == inputs/solver.opts ]]; then
            echo "  FAIL: inputs/solver.opts holds the frozen solver defaults; it must match the tag"
            fail=1
        fi
    fi
done < <(git -C "$REPO" ls-tree -r --name-only "$TAG" -- "$OLD/src" "$OLD/include" "$OLD/inputs/solver.opts")
while IFS= read -r f; do   # files added since the tag
    rel=${f#"$NEW"/}; [[ $rel == src/cryosuction_main.c ]] && continue
    git -C "$REPO" cat-file -e "$TAG:$OLD/$rel" 2>/dev/null || { echo "  added     $rel"; changed=1; }
done < <(git -C "$REPO" ls-files -- "$NEW/src" "$NEW/include")
[[ $changed == 0 ]] && echo "  none: src/, include/ and inputs/solver.opts match the tag"

# --- 2. uncommitted changes in guarded paths ---------------------------------
dirty=$(git -C "$REPO" status --porcelain -- "${GUARDED[@]}")
if [[ -n $dirty ]]; then
    echo; echo "== FAIL: uncommitted changes in guarded paths (commit, then add a ledger row)"
    echo "$dirty" | sed 's/^/  /'; fail=1
fi

# --- 3. every guarded commit since the baseline is in the ledger -------------
echo; echo "== commits since $BASELINE touching guarded paths"
n=0
while IFS= read -r line; do
    [[ -z $line ]] && continue
    n=$((n + 1)); h=${line%% *}
    if grep -q "$h" "$LEDGER"; then echo "  logged    $line"
    else echo "  MISSING   $line"; fail=1; fi
done < <(git -C "$REPO" log --format='%h %s' --abbrev=8 "$BASELINE"..HEAD -- "${GUARDED[@]}")
[[ $n == 0 ]] && echo "  none"

echo
if [[ $fail == 0 ]]; then echo "check_freeze: OK"; else
    echo "check_freeze: FAIL (see FREEZE.md for what is allowed and how to log it)"; fi
exit $fail
