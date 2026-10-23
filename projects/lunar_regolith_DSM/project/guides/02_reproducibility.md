# Reproducibility

The course grades reproducibility at 20 %: another student must be able to
install and run the system. For an LLM workflow, reproducible means a rerun
reaches the same verified conclusions and passes the same checks. It does not
mean identical wording; models are not deterministic, and current Claude
models do not accept a temperature setting.

## What must be pinned

| Item | Where |
|---|---|
| Model per session | `.claude/settings.json` → `"model"` |
| Model per agent | full ID in the agent's frontmatter (`model: claude-sonnet-5-5`), never an alias |
| Claude Code version | `claude --version`, written to each run manifest |
| Python packages | `requirements.txt` (exact versions) |
| Solver build | git commit; PETSc/PetIGA versions in `README.md` |
| Agents, skills, prompts, workflows, settings | git; a run uses only committed versions |

## What each run records

`runs/<workflow>/<run_id>/`:

- `manifest.json`: run ID, date, git commit, Claude Code version, model ID per
  agent, sha256 of every prompt file and input file, check verdicts.
  Committed.
- `transcript.jsonl`: `claude -p "..." --output-format stream-json`. Kept
  locally; selected traces are copied into reports.
- `runs.db`: one row per agent call (see guide 01). Kept locally.

One machine-written manifest replaces hand-written provenance notes. Prose
goes in the report, once.

## Habits that make runs repeatable

- Numbers come from scripts with explicit inputs, never from model memory.
- Every check that can be a script is a script. Its verdict is identical on
  every run; an agent's is not.
- Run workflow stages headless (`claude -p ... --permission-mode dontAsk`),
  so the permission set is the same every run and nothing depends on what was
  clicked.
- Commit before a run; the manifest's commit must contain everything used.
- Test stability: run a stage several times on the same input and count how
  often the verified conclusion changes. Twenty clean runs bound the failure
  rate at about 15 % (rule of three: 3/n).

## What the final repository needs (Project Guideline)

- Source, README with run instructions, `requirements.txt`, pinned model ID.
- Representative configuration and the tests or evaluation scripts.
- Prompts and system instructions (the agent, skill and workflow files),
  selected traces, outputs, plots and data provenance.
- No secrets, restricted data or paper full texts.
