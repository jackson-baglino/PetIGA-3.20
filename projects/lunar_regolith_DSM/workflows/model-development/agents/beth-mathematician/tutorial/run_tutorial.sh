#!/usr/bin/env bash
# The beth-mathematician tutorial, end to end: the agent rederives the current phase-field model blind,
# one equation at a time (T04 heat, T05 vapor, T06 free energy, T07 phase field), each answer is checked by a
# grader that needs no answer key, and the same four work orders are run as single prompts with no tools
# (the baseline). Walk-through: README.md in this folder.
#
#   bash workflows/model-development/agents/beth-mathematician/tutorial/run_tutorial.sh            # live
#   bash workflows/model-development/agents/beth-mathematician/tutorial/run_tutorial.sh --replay   # no model
#
# live     (needs Claude Code, logged in, and macOS for the sandbox) runs, in order:
#            1. the harness and grader tests (no model; each grader proves itself on planted errors)
#            2. the red-team escape probes, which the launcher requires before any real call
#            3. T04 to T07 with the agent           (about 70 s each)
#            4. T04 to T07 as single-prompt baselines (--no-tools)
#            5. a summary table of all eight runs
# --replay (no model, no login; only venv_workflows) runs step 1, then regrades the recorded answers of the
#          runs listed in tutorial/runs.txt and re-renders their derivation.md into a temporary folder.
# The solver comparison (the derived equations written as C code and run against the current model) needs
# PETSc and PetIGA; it is not run here. Its commands and results: studies/derived_model/README.md on branch
# feat/derived-residual.
#
# Every run folder is under runs/model-development/beth-mathematician/<timestamp>_<task>[_baseline]/; read its
# PROVENANCE.md first, then workspace/derivation.md, then answer.json and grade.json.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
cd "$ROOT"
PY=venv_workflows/bin/python
AGENT=beth-mathematician
HERE=workflows/model-development/agents/$AGENT/tutorial
TASKS=(T04_heat_latent T05_vapor_balance T06_free_energy T07_phase_field)
MODE=${1:-live}

step() { printf '\n==== %s\n' "$*"; }

step "1. Harness and grader tests (no model is called)"
$PY -m unittest discover -s scripts/agentkit/tests 2>&1 | tail -3
$PY -m unittest discover -s workflows/model-development/agents/$AGENT/tests 2>&1 | tail -3
echo "(One known failure while the user decides it: test_exchange_counted_twice, the T05 grader's blind spot;"
echo " see docs/planning/30_rederivation_review.md section 1.)"

if [ "$MODE" = "--replay" ]; then
    step "2. Regrade the recorded answers (no model is called)"
    RUNS=()
    while IFS= read -r r; do RUNS+=("$r"); done < <(grep -v '^#' "$HERE/runs.txt" | sed '/^$/d')
    $PY "$HERE/summarize.py" --regrade "${RUNS[@]}"
    OUT=$(mktemp -d)
    for r in "${RUNS[@]}"; do
        $PY scripts/agentkit/render_derivation.py "$r" --stdout > "$OUT/$(basename "$r").md" 2>/dev/null
    done
    echo "derivation.md re-rendered from each answer.json into $OUT"
    exit 0
fi

command -v claude >/dev/null || { echo "Claude Code (claude) is not on PATH; use --replay"; exit 2; }
LOG="runs/model-development/$AGENT/tutorial_$(date +%Y%m%d-%H%M%S)"
mkdir -p "$LOG"
RUNS=()
launch() {   # launch.py prints "<VERDICT>: ... run folder <path>"
    local line
    line=$($PY scripts/agentkit/launch.py "$AGENT" "$@" 2>&1 | tail -1)
    echo "$line" | tee -a "$LOG/launch.log"
    RUNS+=("$(echo "$line" | sed -n 's/.*run folder //p')")
}

step "2. Red-team escape probes (required before any real call on the current harness)"
$PY scripts/agentkit/launch.py "$AGENT" --task T01_heat_weak_form --redteam 2>&1 | tail -1 | tee -a "$LOG/launch.log"

step "3. The agent rederives the model, one equation per call"
for t in "${TASKS[@]}"; do launch --task "$t"; done

step "4. Baseline: the same work orders as single prompts, no tools"
for t in "${TASKS[@]}"; do launch --task "$t" --no-tools; done

step "5. Summary"
printf '%s\n' "${RUNS[@]}" > "$LOG/runs.txt"
$PY "$HERE/summarize.py" --json "$LOG/summary.json" "${RUNS[@]}" | tee "$LOG/summary.txt"
echo
echo "Run folders: $LOG/runs.txt. In each, read PROVENANCE.md, then workspace/derivation.md, then grade.json."
