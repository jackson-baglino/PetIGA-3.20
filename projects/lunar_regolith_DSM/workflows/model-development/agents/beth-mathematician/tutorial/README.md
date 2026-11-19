# Tutorial: the agent rederives our phase-field model, blind

`beth-mathematician` (Claude Opus 5.5, run headless) derives the current two-phase model of ice, water vapor and heat one
equation per call. It is never shown the answer: each work order states the physics and the symbols, and a grader it
cannot see checks the result with physical tests that need no answer key. The same work orders, run as single prompts
with no tools, are the baseline.

| Task | Derives | Grader checks (no answer key) |
|---|---|---|
| T04 | heat equation with latent heat, for temperature $T$ | units, equilibrium, latent-heat sign, conduction, conservation |
| T05 | water-vapor balance, for vapor density $\rho_v$ | units, equilibrium, mass exchange at one ice density, diffusion, conservation |
| T06 | interface free energy $F[\phi]$ | equal wells, isotropic gradient energy, interface energy $=\gamma$ |
| T07 | phase-field equation, for the ice fraction $\phi$ | units, pure phases stay put, growth sign, relaxation lowers $F$ |

## Run it

```bash
bash workflows/model-development/agents/beth-mathematician/tutorial/run_tutorial.sh            # live, about 12 min
bash workflows/model-development/agents/beth-mathematician/tutorial/run_tutorial.sh --replay   # no model, no login
```

Live mode needs Claude Code, logged in, on macOS (the sandbox is macOS Seatbelt). It runs:

1. the harness and grader tests, with no model;
2. the red-team escape probes, which the launcher requires before any real call;
3. T04–T07 with the agent;
4. the same four as baselines (`--no-tools`);
5. a summary table.

The replay mode needs only `venv_workflows`. It regrades the recorded answers of the runs in `runs.txt` and re-renders
their `derivation.md`.

## What happens in one call

1. **Work order.** `tasks/<task>/task.json` holds the statement, the symbols with units and LaTeX, the given equations
   and the names the answer may use. There is no answer key in it.
2. **Launcher** (`scripts/agentkit/launch.py`). It first refuses the call unless a red-team run passed on the current
   harness (the drift gate). It then copies the work order into a fresh `workspace/` and starts `claude -p` in the
   sandbox, which can reach the workspace only and has no network.
3. **Agent loop.** The agent reads `task.json`, writes its checks to `scratch/check.py` (SymPy for algebra and
   calculus, pint for units), runs them with `python scratch/check.py` (the only command it may run; a guard hook
   blocks everything else), reads what they print, and revises. A Stop hook refuses to let it finish before
   `scratch/check.py` exists.
4. **Answer.** A JSON object checked against `schemas/result.schema.json`. It holds:
   - the result as a SymPy expression;
   - each step's words, LaTeX equation and check;
   - the symbols, with units and typical values;
   - the neglected effects, with reasons.
5. **Rendering.** `scripts/agentkit/render_derivation.py` writes `workspace/derivation.md` in one fixed layout and
   lints it.
6. **Grading.** `tasks/<task>/grade.py`, unseen by the agent, runs the physical checks on the agent's own symbols at
   physical sizes. It then lines the result up, term by term, with the current model's code (`src/assembly.c`). Each
   grader first proves itself on planted errors (`key_check`, `tests/test_grade_*.py`).
7. **Record.** `PROVENANCE.md` (hash-chained) and `manifest.json`.

## Read a run

Each run is in `runs/model-development/beth-mathematician/<timestamp>_<task>[_baseline]/`. Read, in order:

1. `PROVENANCE.md`: what ran, the barrier events, the answer and the verdict;
2. `workspace/derivation.md`: the derivation, step by step;
3. `grade.json`: each check and the comparison with the code;
4. `transcript.jsonl`: every tool call and its output.

## Not run here

The derived equations, written as a C residual and Jacobian and run against the current model, need PETSc and PetIGA.
The code, commands and results are in `studies/derived_model/README.md` on branch `feat/derived-residual`.
