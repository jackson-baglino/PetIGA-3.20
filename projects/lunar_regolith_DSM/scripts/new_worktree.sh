#!/bin/bash
# Create a git worktree for parallel work (for example a feat/* physics branch)
# that is ready to build, run and lint, then print how to start Claude in it.
#
#   Usage (from anywhere in the repo):
#     scripts/new_worktree.sh <name> [branch] [--no-make]
#
#   <name>    folder name: the worktree goes to <repo>/.claude/worktrees/<name>/
#   [branch]  an existing branch (local or on origin) is checked out as is; a new
#             name creates that branch from the current local HEAD. Default:
#             new branch worktree-<name> from the current local HEAD.
#   --no-make skip building the solver in the worktree.
#
# Why not plain `claude --worktree <name>`: by default it branches from
# origin/main, so local commits that are not pushed yet are missing, and a fresh
# checkout has none of the git-ignored files the project needs. This script
# starts from local HEAD (or the branch you name) and links, not copies, the
# shared git-ignored files from the main checkout:
#   venv_lunar/, venv_workflows/  one shared environment (pins in requirements*.txt;
#                                 give a worktree its own venv only if its
#                                 branch changes the requirements)
#   CLAUDE.local.md               personal rules
# The solver build (obj/, lunar_regolith_dsm) is per worktree on purpose, so a
# physics branch never overwrites main's binary.
#
# Remove a worktree when its branch is merged or archived (CLAUDE.md):
#   git worktree remove .claude/worktrees/<name>      (from the repo root)
# This deletes the links, not the shared venvs they point to.
set -euo pipefail

usage() { sed -n '5,12p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }

name="" branch="" make_it=1
for a in "$@"; do
    case "$a" in
        --no-make) make_it=0 ;;
        -h|--help) usage ;;
        -*) echo "unknown option: $a" >&2; usage ;;
        *) if [ -z "$name" ]; then name="$a"
           elif [ -z "$branch" ]; then branch="$a"
           else echo "too many arguments" >&2; usage; fi ;;
    esac
done
[ -n "$name" ] || usage
case "$name" in */*|.*|*" "*) echo "name must be a plain folder name: $name" >&2; exit 2 ;; esac

ROOT="$(git rev-parse --show-toplevel)"
PROJECT_REL=projects/lunar_regolith_DSM
MAIN_PROJECT="$ROOT/$PROJECT_REL"
WT="$ROOT/.claude/worktrees/$name"
WT_PROJECT="$WT/$PROJECT_REL"

[ -d "$MAIN_PROJECT" ] || { echo "not found: $MAIN_PROJECT" >&2; exit 1; }
[ ! -e "$WT" ] || { echo "already exists: $WT" >&2; exit 1; }
[ ! -L "$ROOT/.claude" ] && [ ! -L "$ROOT/.claude/worktrees" ] ||
    { echo ".claude or .claude/worktrees is a symlink; refusing" >&2; exit 1; }

branch="${branch:-worktree-$name}"
git check-ref-format --branch "$branch" >/dev/null ||
    { echo "not a valid branch name: $branch" >&2; exit 2; }

mkdir -p "$ROOT/.claude/worktrees"
if git -C "$ROOT" show-ref --verify --quiet "refs/heads/$branch"; then
    echo "checking out existing branch $branch"
    git -C "$ROOT" worktree add "$WT" "$branch"
elif git -C "$ROOT" show-ref --verify --quiet "refs/remotes/origin/$branch"; then
    echo "checking out origin/$branch as local branch $branch"
    git -C "$ROOT" worktree add --track -b "$branch" "$WT" "origin/$branch"
else
    echo "new branch $branch from local HEAD ($(git -C "$ROOT" log -1 --format='%h %s'))"
    git -C "$ROOT" worktree add -b "$branch" "$WT" HEAD
fi

for f in venv_lunar venv_workflows CLAUDE.local.md; do
    if [ -e "$MAIN_PROJECT/$f" ]; then
        ln -s "$MAIN_PROJECT/$f" "$WT_PROJECT/$f"
        echo "linked $f -> main checkout"
    else
        echo "skipped $f (not in the main checkout)"
    fi
done

# The links must be ignored by git, or they could be committed by accident.
untracked=$(git -C "$WT_PROJECT" status --porcelain --untracked-files=all -- venv_lunar venv_workflows CLAUDE.local.md)
if [ -n "$untracked" ]; then
    echo "warning: git does not ignore these links in the worktree:" >&2
    echo "$untracked" >&2
fi

if [ "$make_it" -eq 1 ]; then
    echo "building the solver in the worktree (make)..."
    make -C "$WT_PROJECT" >/dev/null
    echo "built $WT_PROJECT/lunar_regolith_dsm"
fi

cat <<EOF

Worktree ready: $WT_PROJECT  (branch $branch)
Start Claude in it:
  cd "$WT_PROJECT" && claude
Remove it when the branch is merged or archived:
  git -C "$ROOT" worktree remove "$WT"
EOF
