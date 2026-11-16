# WORKFLOW: model-development (the model line)

Status: **design approved by the user 2026-11-16, with one change** (the Fable
check is not part of this line). S3:Build is next. Nothing below
is built except `library_search.py` and `add_paper.py` (S1:Library). The
flowchart is `flowchart.html` in this folder.

Background, findings and the answers to the user's questions:
`docs/planning/25_model_line_plan.md`.

## 1. Purpose

**Question.** What model of ice, water vapor and heat in a regolith pore space
should this project solve? How exactly is it derived, discretized and
implemented, and does the implementation do what the math says?

**A finished output** is four things, tied together by one run folder:
1. a derivation with every step sourced or derived, and every claim checked by
   a script someone other than the author wrote;
2. a frozen specification of the weak form, residual, Jacobian, auxiliary
   functions and time-step rule;
3. C code behind a runtime option, default off, on a `feat/<model>` branch,
   passing a locked test list whose tolerances the user set;
4. a linear teaching document that explains the model, the derivation and the
   code to the user.

**Acceptance test.** The output counts as correct when all of these hold:
- every verifier check passes;
- the locked tests pass, including:
  - the regression to today's model;
  - the neck-width comparison of two sintering grains with no sediment;
  - manufactured solutions at the expected order;
  - the full Jacobian test;
- the mutation check catches the planted bugs;
- the user has read the teaching document and accepted it at G:Write-up.

**Baselines (non-agent), used for comparison:**
- **Derivation:** doc 18, the F2 pore-space derivation of 2026-10-18. It is
  hidden from the agents, and a script diffs the agent's strong form against
  it term by term.
- **Tests:** a fixed rule table ("Jacobian test fails → code") is the baseline
  for the foreman's triage.
- **Solver settings:** the current `inputs/solver.opts` is the baseline for
  M8:Performance.

## 2. Stages

Setup steps S1–S3 are done once. (The daily Fable access check is a general
session hook, `~/.claude/hooks/fable_check.sh`, not part of this line; user, 2026-11-16.) Line stages M1–M9 run through the line script
`scripts/model_development/run_md.py`. Each stage ends in a file that something
other than its producer checks.

| # | Stage | Agent or script | Input | Output | Check | Stop / escalate |
|---|---|---|---|---|---|---|
| M1 | Scope | beth-modeler | work order; `src/`, `include/`; library | `derivation/spec_problem.md` | verifier: `checks/m1_*.py` (Pint units, scaling numbers recomputed); `library_search.py quote` for each citation | gate G:Scope |
| M2 | Selection | beth-modeler | `spec_problem.md` | `derivation/selection.md`, `derivation/citations.json` | quote check on every citation | gate G:Selection |
| M3 | Derivation | beth-modeler | M1 and M2 outputs | `derivation/ledger.json`, `derivation/derivation.md` | verifier: `checks/m3_*.py` (reduction to `Residual_A1`, conservation, energy, units, each named operation, unlisted assumptions); `compare_doc18.py` (evaluation only) | gate G:Derivation |
| M4 | Discretization | beth-modeler | ledger | `derivation/spec.md`, `derivation/code_map.json`, `derivation/mms/` | verifier: `checks/m4_*.py` (weak ↔ strong, symbolic Jacobian vs differentiated residual, dimensionless groups) | gate G:Discretization (spec frozen by hash) |
| M5 | Test lock | beth-verifier | `spec.md`, `code_map.json`, `main` | `studies/<model>/tests/`, `tests.json` | every test runs; regressions pass on `main`; each test inside the run caps | gate G:Test lock (tests frozen by hash) |
| M6 | Implementation | beth-modeler ⇄ tests; beth-model-foreman on failure | spec, test list | C code on `feat/<model>` | build; full Jacobian test; locked tests | loop limits (§3.3); dispute → user |
| M7 | Mutation | script `mutate.py` | branch, test list | `mutation.json` | each planted bug caught by at least one test | missed bug → verifier adds a test → user re-locks |
| M8 | Performance | beth-modeler, measured by `perf.py` | branch | `perf.csv`, settings | all tests still pass after each change | review gate |
| M9 | Write-up | beth-writer (+ student review) | ledger, spec, check and test results, `perf.csv` | `docs/model/<model>.md` | `check_doc.py`: every symbol defined before first use, every number traced to a file; citation check | gate G:Write-up |

Each stage is described in full (what, why, how, done when, then) in the
approved plan; the essentials are below.

### M1: Scope (spec and scaling)

- **What:** `spec_problem.md`, the contract every later stage reads. It holds:
  - the physical problem and the questions the model must answer;
  - the mechanisms in and out, each with a scaling estimate:
    - vapor diffusion against Knudsen flow in vacuum;
    - surface diffusion;
    - radiation in pores;
    - grain conduction;
  - the constitutive laws and their validity ranges for lunar temperatures:
    ρ_vs, D, P_sat, k_ice, c;
  - mixed conductivity, by the Nicoli et al. (2011) arithmetic–harmonic blend;
  - regolith k and c_p, from the literature for the lunar composition or from
    a proxy material;
  - the notation lock: ψ, φ, φ_i = ψφ, φ_a = ψ(1 − φ), φ_s = 1 − ψ,
    ψ̃ = ψ + υ, eps, eps_s, plus the sign of the normal and the units;
  - what is out of scope.
- **Why:** today's constitutive laws are Earth-air laws, and a derivation
  without a locked scope drifts.
- **Done when:**
  - the verifier's checks pass;
  - every citation passes the quote check and is listed for the user;
  - **G:Scope:** the user chooses:
    - the scope;
    - conducting or insulating grains;
    - the source of the regolith properties;
    - air or lunar laws for the first version.
- **Then:** M2.

### M2: Selection (model type)

- **What:** at least four candidates, scored against the M1 requirements:
  - body-fitted two-phase (today's model);
  - diffuse domain or smoothed boundary;
  - three-phase phase field;
  - sharp interface or level set.

  The criteria:
  - fidelity: contact line, no spurious film, conservation;
  - geometry flexibility;
  - cost in DOFs;
  - fit to PetIGA against another code;
  - whether it can be built and verified by December.
- **Why:** the method is not fixed in advance. The three-phase model's
  spurious air film is a known failure, and any candidate must avoid it.
- **Done when:** every citation passes the quote check, and **G:Selection:**
  the user picks the type and names the model (the `<model>` in
  `feat/<model>`).
- **Then:** M3.

### M3: Derivation

- **What:** `ledger.json` and `derivation.md`. They cover:
  - the strong form;
  - the scales and dimensionless groups, with a validity criterion for ξ;
  - the floors (smooth, on products);
  - the wall or contact condition, starting from today's validated boundary
    term;
  - whether the thin-interface correction applies;
  - the sharp-interface limit to first order.

  Every ledger step is one of:
  - **cited:** paper, equation, page and quote;
  - **derived:** from steps i and j by a named operation;
  - **assumption** A#;
  - **choice** C#.
- **Why:** everything after rests on these equations, and the sharp-interface
  limit is the evidence that the model means what we say.
- **How:**
  - the modeler derives;
  - the verifier writes the SymPy checks from the ledger's *claims*, not from
    its reasoning;
  - `compare_doc18.py` produces the term-by-term diff against doc 18 for the
    user.
- **Done when:**
  - all checks pass;
  - no step lacks a source;
  - **G:Derivation:** the user answers each C# and confirms the cited claims.
- **Then:** M4.

### M4: Discretization (spec freeze)

- **What:** `spec.md` holds:
  - the weak form;
  - the IGA spatial discretization and generalized-alpha in time;
  - the residual and Jacobian entries;
  - the auxiliary functions: the indicator and its gradient, the initial
    conditions and their extension into grains, the post-processing masks;
  - the a-priori time scales: interface relaxation, tau_kin, and vapor and
    thermal diffusion at the mesh and domain sizes;
  - the dt rule: an interface-displacement limit (≤ c·eps per step) combined
    with an error-based controller;
  - the manufactured-solution sources.

  `code_map.json` maps each term to its `assembly.c` location and its runtime
  option.
- **Why:** the coder and the tester must work from one fixed document.
  Today's dt rule uses no scale a priori.
- **How:**
  - the modeler first checks in the PETSc 3.20 source whether TSAdapt's error
    estimate works with TSALPHA, or whether the `//AM` patch must stay;
  - the verifier checks:
    - weak ↔ strong;
    - the symbolic Jacobian against differentiation of the symbolic residual;
    - that the groups are consistent.
- **Done when:** the checks pass and **G:Discretization:** the user freezes the
  spec; its hash is recorded.
- **Then:** M5.

### M5: Test lock

- **What:** the test list (`tests.json`) and its scripts:

  | Test | What it checks |
  |---|---|
  | Regression | Option off is identical to `main`; ψ ≡ 1 agrees with today's model. |
  | Full Jacobian | The analytic Jacobian, including the shift terms. |
  | Manufactured solutions | Convergence order in h and in dt. |
  | Closed-box water | Total water is conserved. |
  | Energy | The free energy decays. |
  | 1D front speed | Speed of a front in a tapered pore. |
  | 90° contact | The contact angle at a wall. |
  | Separate refinement | eps and the grain width, refined one at a time. |
  | ξ sensitivity | The effect of the time-scaling factor. |
  | Neck width | Two sintering ice grains with ψ ≡ 1, new model against today's model (same mesh, eps, parameters), neck width against time. The neck measure itself is tested on a known shape first. |
  | K&P lamellae and bubble | Regressions, if those cases are confirmed in the paper. |
  | Cost budget | A cost limit per test. |
- **Why:**
  - the tests are written from the spec, before the code, by an agent that
    never sees the code;
  - the tolerances are set by the user, so the coder cannot loosen them.
- **Done when:**
  - every test runs;
  - the regressions pass on `main`;
  - **G:Test lock:** the user sets the tolerances, and the tests are hashed and
    locked by a hook.
- **Then:** M6.

### M6: Implementation (in layers)

- **What:** the C code, layer by layer: first "infra" (a manufactured-solution
  forcing option and a full Jacobian test that includes the shift terms; the
  existing `-test_wall_jacobian` checks dR/dU only), then heat, vapor, phase
  field, coupling, wall term. Each layer must pass before the next.
- **The loop** (§3.3):
  1. **Each layer:** build, run the full Jacobian test, run the layer's tests,
     then run all earlier tests.
  2. **A failure** goes back to the modeler, with at most 3 fix rounds.
  3. **The foreman then triages:**
     - a code bug goes to the modeler;
     - a spec error goes to M3 or M4, with a gate;
     - a test bug goes to the verifier;
     - a tolerance or environment problem goes to the user.
  4. **Every passing run logs** wall time, Newton and Krylov iterations and dt.
- **Done when:** the locked list is green and the last layer is committed.
- **Then:** M7.

### M7: Mutation

- **What:** about 10 planted bugs, each applied in a scratch worktree. For
  example:
  - the sign of the wall term;
  - a missing weight;
  - a dropped Jacobian term;
  - a wrong floor;
  - eps_s below the mesh resolution.

  The record says which test caught each bug.
- **Why:** this measures the verifier. A green suite that catches nothing
  proves nothing.
- **Done when:** every mutation is caught, or each miss has a new test and the
  user has re-locked the list (G:Test lock reopened; blocking).
- **Then:** M8.

### M8: Performance

- **What:** settings chosen by measurement, recorded in `perf.csv` and committed
  as the model's defaults on the branch.
- **How:**
  - check that the nondimensional solve makes the tolerances act per field;
  - on two fixed cases, compare:
    - BiCGStab or GMRES;
    - ASM with ILU(k), block Jacobi, field split (ice | T | ρ_v);
    - MUMPS for 2D, if the PETSc build has it;
  - check the dt rule against a fixed small-dt reference;
  - run strong scaling at 20k, 50k and 100k DOFs per rank, which sets this
    project's local target in `alloc.sh` behind a flag.

  All tests stay green after every change.
- **Then:** M9.

### M9: Write-up

- **What:** `docs/model/<model>.md`, in this order. Each part uses only what
  came before.
  1. the problem;
  2. unknowns and assumptions;
  3. the strong form, term by term, with meaning and source;
  4. the nondimensional form;
  5. the sharp-interface limit and why it justifies the model;
  6. the weak form;
  7. discretization;
  8. the residual and Jacobian;
  9. the dt rule;
  10. the code map;
  11. each test and what it proves;
  12. performance choices;
  13. limitations.
- **How:**
  - beth-writer writes it from the files, not from memory;
  - `check_doc.py` checks that each symbol is defined before its first use and
    traces every number to a file;
  - a fresh-context "student" review reads only the document. It must redo the
    reduction at ψ ≡ 1 and one Jacobian entry, and it lists every step it
    could not follow.
- **Done when:** the scripts pass, the student's list is resolved, and
  **G:Write-up:** the user has read it.

## 3. Agents

All three new agents are defined in `workflows/model-development/agents/`.
`beth-writer` is reused from `.claude/agents/`.
- **Model:** all four run on `claude-opus-5-5`; Fable is not available
  (doc 25 §1.3).
- **How they run:** headless (`claude -p`) through `scripts/agentlib.py`, with
  fresh context per call.
- **Permissions:** `permissionMode: dontAsk`, plus a per-agent settings file in
  `.claude/permissions/model-development/`.
- **Output:** each call returns one JSON object in its schema
  (`schemas/result.schema.json`): `status` (done | blocked | dispute), the
  output paths, the evidence, and the open questions for the user.

```
name:         beth-modeler
model:        claude-opus-5-5
tools:        Read, Grep, Glob, Write, Edit, Bash
skills:       library-search
MCP servers:  none
reads:        src/, include/, inputs/, memory/semantic/, docs/interface_width_conventions.md,
              library/papers/** (through library_search.py and the page files it names),
              its run folder; in M6 also the feature worktree and test results
writes:       <run>/derivation/**; in M6 and M8 src/ and include/ in the feature worktree only
bash:         venv_workflows/bin/python <run>/derivation/scratch/*.py (its own SymPy work);
              library_search.py search|quote; make (worktree); run_test.sh <test id>;
              git status|diff (worktree)
may not:      read docs/planning/** (independence from doc 18); edit tests, checks, tests.json
              or tolerances; touch main; use the web; run anything but run_test.sh
context:      fresh per task; state passes through the run folder
```

```
name:         beth-verifier
model:        claude-opus-5-5
tools:        Read, Grep, Glob, Write, Edit, Bash
skills:       none
reads:        <run>/derivation/ outputs (spec_problem, ledger, derivation.md, spec, code_map),
              src/ on main (for regressions), memory/semantic/units.md, its own checks
writes:       <run>/checks/**; studies/<model>/tests/** and tests.json (M5 only, and when the
              foreman routes a test bug to it)
bash:         venv_workflows/bin/python <run>/checks/*.py; run_test.sh <test id> (M5: on main)
may not:      read the feature branch's src/ or include/ or its diff; read docs/planning/**;
              edit src/ or include/; use the web
context:      fresh per task
```

```
name:         beth-model-foreman
model:        claude-opus-5-5
tools:        Read, Grep, Glob, Write
reads:        the run folder: test logs, the failing test, the spec section, the diff the
              script saved, earlier triage records
writes:       <run>/triage/<n>.json, or <run>/gates/<gate>.md (the gate summary)
may not:      run anything; edit code, tests, specs or definitions
context:      fresh per call; one decision per call
```

`beth-writer` (existing general helper, `.claude/agents/beth-writer.md`, tools Read,
Grep, Glob, Write) is called twice in M9, each time fresh: once to write the
document, and once in review mode as the "student", given only the finished
document.

The library is a search script plus a skill (`library-search`), not an agent.
The modeler opens the page itself before relying on a passage.

### 3.1 The line

```
user (owner)            approves the design, every gate, the test tolerances, the final document
  run_md.py             walks the stage graph in flowchart.json; refuses a stage whose inputs are
                        not done; keeps runs.db; stops at gates and on a usage-limit error
    station scripts     one per stage: write the work order from a fixed template, call the
                        agent through agentlib.py, audit the transcript, run the checks,
                        retry once on a malformed result, write report.json
      agents            beth-modeler, beth-verifier, beth-model-foreman, beth-writer
```

**Gates** (all blocking): G:Scope, G:Selection, G:Derivation, G:Discretization,
G:Test lock, G:Write-up. At each gate the station calls `beth-model-foreman`
once, with fresh context, to write `gates/<gate>.md`: what was produced, what the
checks found, and the decisions the user must make, each as a direct question.
A test added after an M7 mutation miss reopens G:Test lock, which blocks M8
until the user re-locks the list.

**Review items** (queued; the line goes on):
- the doc 18 comparison;
- the M8 settings.

### 3.2 Keeping agents on task

| Mechanism | What it does |
|---|---|
| The script owns the loop | No agent waits, polls or remembers. State lives in `runs.db` and files. |
| One work order per call | One deliverable, a written done criterion, `maxTurns`, and a timeout (agentlib's default 1200 s; 3600 s for M3, M4 and M6). |
| Stop hook in each agent's frontmatter | Runs the stage's structural pre-check, for example ledger schema valid and no unsourced step. It blocks the end of the turn while the check fails, unless the agent returns `blocked` with a reason. Claude Code's 8-block cap and `stop_hook_active` stop a loop. |
| Schema output | `--json-schema`: a result without evidence or paths is rejected and retried once. |
| Transcript audit (`agentlib.py`) | Checks the pinned model, the allowed tools only, reads inside bounds, and no permission denials. |
| Path guards (PreToolUse hooks) | The modeler cannot write tests or checks; the verifier cannot write `src/`; no agent writes to `main` or outside its paths. |
| Usage-limit handling | A limit error is recorded, the script stops with the resume command, and it never retries in a loop. |
| Status | `run_md.py status` prints one line; `notify.sh` alerts at a gate. |

### 3.3 The iterate loop (M6)

1. The modeler gets one failing test, its log and the spec section, and has at
   most 3 fix rounds for that failure.
2. If the test still fails, the foreman gets the evidence. It returns one route,
   with its reason:
   - `code`: back to the modeler;
   - `spec`: back to M3 or M4, which reopens that gate;
   - `test`: to the verifier;
   - `tolerance` or `environment`: to the user.
3. The modeler may dispute a test, with evidence only. The verifier answers in
   one call, and the foreman rules. An unresolved dispute goes to the user.
4. The foreman's routes are compared, at the end, with the fixed rule table
   (§1 baseline).

## 4. Checks

| Output | Check | Kind | Pass criterion |
|---|---|---|---|
| any citation | `library_search.py quote` | deterministic | quote on the cited page (whitespace and case normalized only) |
| scaling numbers (M1) | `checks/m1_*.py` with Pint | deterministic, independent author | recomputed value within 1 % and the same units |
| ledger (M3) | schema check plus `checks/m3_*.py` (SymPy) | deterministic, independent author | every check prints PASS; no step without a source |
| strong form vs doc 18 | `compare_doc18.py` | deterministic diff | review item for the user (evaluation, not a gate) |
| spec (M4) | `checks/m4_*.py` | deterministic, independent author | weak ↔ strong identity; Jacobian entries equal the derivatives |
| code | build; full Jacobian test; locked tests | deterministic | the tolerances the user set at G:Test lock |
| tests | `mutate.py` | deterministic | every planted bug caught |
| triage | rule table at the end | deterministic baseline | agreement reported, not gated |
| document | `check_doc.py`; citation check; student review | deterministic plus independent agent | no undefined symbol, no untraced number; the student's list resolved |

**Measuring the checker** (TEMPLATE §4): M7 measures the verifier's tests with
planted bugs. The verifier's SymPy checks are measured the same way: before
M3's checks count, three planted errors in a copy of the ledger (a sign, a
missing factor, a wrong exponent) must each make a check fail.

## 5. Human decision points

| Gate | The user is shown | The user decides | What breaks if skipped |
|---|---|---|---|
| G:Scope | `spec_problem.md`, the verifier's results, the citations | scope; grain heat; regolith property source; air or lunar laws | the derivation drifts or uses laws that do not hold in vacuum |
| G:Selection | `selection.md` with the sourced scores | model type and name | the line derives a model the user did not choose |
| G:Derivation | `derivation.md`, the check results, the C# list, the doc 18 diff | each choice; the cited claims | an unchecked equation or claim becomes the spec |
| G:Discretization | `spec.md`, check results | freeze | code and tests could target different equations |
| G:Test lock | the test list with the proposed tolerances and reasons | the tolerances | the coder could pass tests the user would not accept |
| G:Write-up | the document and the student's list | accept or send back | the user does not learn the model |

**User actions outside the gates:**
- deciding the three B1 claims the model rests on (KP c001, c012; MF c011);
- running `add_paper.py` and the index build, since the library is closed to
  normal sessions.

## 6. Reproducibility

Each run writes `runs/model-development/<run_id>/`:
- `manifest.json` (written by `run_md.py init`):
  - the git commit and dirty flag;
  - the Claude Code version;
  - the model ID per agent;
  - the hashes of agent files, prompt templates, schemas and check scripts;
  - the library index hash (`library_fts.sqlite`, `files_sha256`);
  - the frozen spec and test hashes, once locked;
  - the check verdicts.
- `runs.db` (SQLite): one row per agent call and per check (stage, agent,
  model, input hash, output path, verdict, timestamp), so a run can be resumed
  and audited.
- `transcripts/<stage>/<label>_try<k>.jsonl` from `--output-format
  stream-json`.
- `derivation/`, `checks/`, `triage/`, `gates/`, `perf.csv`, `mutation.json`.

Code and tests are on `feat/<model>`. Their commits are recorded in the
manifest at each gate. Solver runs go through `run_test.sh` →
`scripts/Studio/run_lunar.sh`, which stages a reproducible copy of every run.

Rerunning a stage with the same inputs must reach the same verified
conclusions, though not identical wording.

## 7. Phases and freeze

- **Develop:** S3:Build and the toy dry run. Runs are tagged `dev` and are not
  results.
- **Freeze:** before the first real run, the user approves
  `freeze/md-v1.json` (the hash of every definition) and a git tag. The phase
  becomes execute.
- **Execute:** `run_md.py` refuses to start if a frozen hash differs, and a hook
  asks before any edit to a frozen file. Ideas go to `BACKLOG.md`.
- **The one way back:** a real error found by a check stops the run. The phase
  returns to develop, the fix is recorded, and `md-v2` is approved before the
  run resumes.

The spec and the tests are products of a run, not definitions. They are
frozen by their own gates (G:Discretization, G:Test lock), with the hash in the
manifest.

## 8. Risks and failure modes

| Risk | How it would show | What catches it |
|---|---|---|
| The modeler's derivation copies a known answer instead of deriving | matches doc 18 with no derivation steps | ledger check: every step derived or cited; doc 18 is hidden by permission |
| The verifier's checks are too weak | everything passes | planted ledger errors (§4); M7 mutation |
| A cited passage is misread | wrong equation used | quote check, plus the user confirms cited claims at the gates |
| Paper text carries instructions | an agent acts on text from a paper | paper text is data; agents have no web or push; the audit flags out-of-bounds actions |
| Tests written to the code, not the spec | tests pass on a buggy branch | the verifier cannot read the branch's `src/`; M7 |
| Runaway runs | CPU saturated; one run blocks others | `run_test.sh`: lock file, ≤ 8 ranks, 20k–100k DOFs per rank, 30 min timeout, locked list only |
| Floor slows Newton or breaks the per-field test | many Newton iterations; dt collapses | logged per run in M6; M8 solver study |
| eps_s under-resolved | wall error grows | separate refinement test in M5 |
| Usage limit mid-stage | partial outputs | the stage is atomic: outputs are kept only when the call returns `done` and the checks pass; resume restarts the stage |
| The agent picks a model the user would not | — | G:Selection |
