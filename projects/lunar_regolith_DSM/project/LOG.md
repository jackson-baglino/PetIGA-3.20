# Running log

Newest first. One entry per session: done, remaining, bottlenecks.

## 2026-11-04 (logged 2026-11-13): Execute run 20261105-pilot2 started; C0 passed

- Session: f83f89f4-9803-416c-ac96-bc7d58c2305c (a background job from the foreman session; it ran as plain
  Claude Code, not as beth-foreman)
- Done:
  - Initialized execute run `runs/literature-review/b1/20261105-pilot2` under b1-v2 (3296e375). It has nine
    papers, commit e0b3c296 and a clean tree.
  - The user approved U0 and ran the line. C0 (held-out calibration, 2353.7 s) passed on its first attempt:
    rule kinds 14/14, sign or exponent 6/6, E4 plants 30/31 (lower bound 85.6 %), controls flagged 3/40. The
    one miss was a dropped condition.
  - Then A2 ran for three papers and P1 for one.
- Remaining (decide in the beth-foreman session, not a background job):
  - The C0 gate, "accept the frozen calibration result".
  - The A2 gate for li_2009: the copy kind is "unknown" and no DOI was found on pages 1–2.
  - The A2 gate for lervag_lowengrub_2015: the copy is an arXiv version and no DOI was found on pages 1–2.
    Check both copy kinds against the contract's Copy kinds paragraph.
  - P1 yu_chen_thornton_2012 FAILED (exit 1, 0.1 s). The log is
    `untrusted/logs/P1_yu_chen_thornton_2012_20261104-195827.log`; read it in review mode. Its A2 was clean:
    publisher copy, DOI found.
  - Ask the foreman why A2 and P1 ran while the C0 gate was still waiting.
  - Then rerun the line: `run_b1.py run --run runs/literature-review/b1/20261105-pilot2`.
- Bottlenecks:
  - `run_b1.py --run` takes a path. A bare name made a stray `./20261105-pilot2/`, which was deleted before
    any stage ran.
  - Background jobs started from the foreman session do not inherit `--agent beth-foreman`, so this
    session's setup steps (init, deleting the stray folder, a commit) ran without the foreman's limits.
    Candidate backlog item: block line commands outside beth-foreman.

## 2026-11-02/03: Pilot line closed out; literature review v2 designed; isolation fix

- Session: bd0b31b7-3240-4887-9a00-1fe467db9be6 (continued after compaction; still open)
- Done, 2026-11-02 (B1 line, develop phase):
  - Build tests passed:
    - T-trifecta (9dcf731d, 89c9aa89);
    - T-applicability;
    - T-restore (491 archived, 208 restored, 0 mismatches).
  - Line fixes:
    - a waiting gate holds only its item (2513649b);
    - X1 calls go to the ledger (518f4297);
    - status shown in dependency order (e9d4f737);
    - C0 with fewer than 19 plants is reported, not failed (54c3e90f).
  - Tuning diagnostic `calib-dev6`: 32/32 plants caught, **9/40 controls flagged**, against the rule of at most 4
    (d2296ac2).
  - The user's U3 decisions now reach the library:
    - new stage L1b (142657d2);
    - `u3_review.py --redo` (35e56cfd).
  - Dev run `dev-20261101-kp` complete end to end, U4 approved (c73554b1).
- Done, 2026-11-03:
  - **v2 design agreed** with the user over four rounds.
    - Planning notes `docs/planning/19`–`22` are uncommitted on purpose, because they discuss course examples
      (Lecture 10, Scenario C).
    - Result: `workflows/literature-review/STORY.md` (the user's version) and `SPEC.md` (the agent spec; replaces
      WORKFLOW.md at Build 10) (7fccc572, d8ce3290).
    - `TODO.md`: six phases from today to a working library, with "⟂ depends" branches (27367234).
  - Key decisions:
    - Claude Code on the user's login only (no API key, no Agent SDK);
    - judgments as single schema-checked calls;
    - only the scout loops, run through agentlib;
    - usage cap at 90 % of the five-hour window (no weekly limit);
    - two library layers;
    - paragraph IDs;
    - condition labels: necessary, scope, unnecessary;
    - stage codes F/G/R/L/M/K/C, written "code: name";
    - the foreman stays as the run interface, with a planned guard hook (H8).
  - **Marker removed**; `venv_workflows` rebuilt (90 → 22 packages, `anthropic` 2.4.0); weights cache deleted
    (0eff0689).
  - **Worker isolation canary test:**
    - our flags kept CLAUDE.md out but loaded the repo-root auto memory index;
    - fixed with `CLAUDE_CODE_DISABLE_AUTO_MEMORY=1` in agentlib and tres (61446b41);
    - `--safe-mode` tested as an equivalent; `--bare` needs an API key.
  - **c006:**
    - our premise was wrong: the checker had called the statement supported;
    - `u3_review` now shows why a claim is disputed (3892163f);
    - the user kept "right" (the quote supports only the stated half; the other half is an inference);
    - L1b republished (34437223);
    - backlog item: a form for inferred facts (6b13e453).
  - Memories added (auto memory):
    - long answers go to markdown;
    - account limits;
    - v2 design;
    - descriptive stage names;
    - say why for user tasks.
- Remaining:
  1. **User: TODO 1.3**, the nine control rulings (G001 G004 G007 G014 G018 G031 G032 G036 G039) in
     `runs/literature-review/calib-dev6/untrusted/control_review_tuning.md`. Judge whether the paper establishes
     the statement (the c006 lesson), not whether it is quoted word for word. Then TODO 1.4–1.11 (gold fixes or
     prompt change, tuning rerun, backlog decision, close dev run `dev-20261101-a2`, pin Claude Code 2.1.317, the
     freeze), then Phase 2, the nine-paper pilot run.
  2. The v2 build (TODO Phase 3) after the pilot. The `.claude/` hook (Build 4) and three new agent definitions need
     the user's approval.
  3. `project/PILLARS.md` Pillar 2/3 framing is still old (from 2026-11-01).
- **For a slides session (class share, 2026-11-04).** Sources, verified this session:
  - `STORY.md` (the whole design in plain words) and `SPEC.md` §3 (roles table) and §4 (stage table);
  - `flowchart.html` (the v1 line);
  - the Lab Report 2 figures in the course folder.

  Numbers you can cite, each with its source:
  - one-paper trial: about 12 min machine time per paper (`runs/literature-review/b1/dev-20261101-kp/runs.db`);
  - calibration tuning: 32/32 plants, 9/40 controls (d2296ac2);
  - T-restore: 0 mismatches;
  - canary test: control leaked everything; our flags leaked auto memory; fixed (WORKFLOW §12);
  - the D/R audit and the four reasons for a separate agent (STORY §2).

  Course material and slides stay in the course folder, never committed.
- Bottlenecks:
  - Untrusted run files can't be read in a normal session (by design), so c006's failed test had to be found from
    code. The review tool's display was the real gap.
  - The auto-mode classifier denies reading `untrusted/` even for label-only scripts.

## 2026-11-01 — Lab Report 2; beth-writer; README; pushed

- Session: bd0b31b7-3240-4887-9a00-1fe467db9be6
- Done: Lab Report 2 written and built in the course folder
  (`ME295/LabReport2/`, not committed): report.md, facts.py (every number
  with its source; refuses hand-typed numbers), build_report.py (HTML ->
  docx/pdf via soffice), four simplified Graphviz flowcharts, motivation
  from Johnson et al. 2025 (Zotero 7VPUC8G4), generic agent names, Pillar 2
  = derivation, Pillar 3 = implementation and testing. New agent
  `.claude/agents/beth-writer.md` (Opus; Read/Grep/Glob/Write; approved,
  d9dddd20) reviewed the report twice (34 then 17 findings; all applied
  except those left to the user); reviews in `LabReport2/review/`. README
  updated with the workflows, layout and tests (d6f6f8f7). Pushed origin and
  synced the course repository (course main 482a1760); the 35 unpushed
  commits are now on both.
- Remaining:
  1. User: fill the "[AUTHOR: ...]" line in the report's AI Use
     Declaration (or ask to have the proposed text inserted), read the PDF,
     submit; then update the deadline line in CLAUDE.local.md.
  2. User decisions left open by the report: an abstention rule and
     temperature note in WORKFLOW.md §5; adding the script-only checker
     baseline to §5 before the freeze; related-work claims from read papers.
  3. project/PILLARS.md still has the old Pillar 2/3 framing; update it to
     Pillar 2 = model derivation, Pillar 3 = implementation and testing
     (needs the user's OK).
  4. Freeze hashes omit the parser version and venv_docparse has no lock
     file (found by the review).
  5. Resume the B1 build from the 2026-11-01 line session
     (47a4c410-b4fb-46dc-9c25-6c04ed2b31e1): its Remaining list above.
- Bottlenecks: headless agents' file writes need an `Edit(//abs/path/**)`
  allow rule; `Write(...)` rules were denied twice (findings recovered from
  the denied calls in run JSON).

## 2026-11-01 — Literature review as a line: foreman, stations, hooks, sandbox

- Session: 47a4c410-b4fb-46dc-9c25-6c04ed2b31e1
- Done: plan to keep agents on task (plan file in ~/.claude/plans/), then
  built it (8d4ae5d5..HEAD, ~27 commits, 112 tests). Develop/freeze/execute
  phases and "the line" in workflows/TEMPLATE.md and WORKFLOW.md §0;
  .claude/rules/ (workflows, solver-c, inputs), CLAUDE.md trimmed;
  scripts/agentlib.py shared runner (audit, schema, retry, runs.db); all
  literature-review agents in workflows/literature-review/agents/ (6 workers
  + beth-foreman, started by scripts/literature_review/foreman.py);
  run_b1.py (stage graph = run blocks in flowchart.json; gates; freeze/thaw;
  --dev --papers --skip), lr_status.py, lr_backlog.py, BACKLOG.md; hooks
  H1-H7 and skills lr-status, lr-freeze, lr-backlog, new-agent; stations
  check_copy (A2), P1 link in read_check, revise (E5), curator (L1),
  library-based Y1/Y2/Y3 with Q2/Q3/Leads and the one Y3 revision,
  blind_read + compare_reads (X1), u3_review, backup + restore_test;
  wf_launch (sandbox per station; agent stations need a deny-list profile,
  a write allow-list broke Claude Code's login and raised a keychain
  dialog) and the review session (Read/Grep/Glob only, tested live).
  A2 found 5 arXiv copies; the user put published copies in Zotero and
  ran zotero_pull. Dev run runs/literature-review/b1/dev-20261101-kp
  (Kaempfer & Plapp, C0 skipped): 23 claims, E3 18/23, E5 10 confirmed,
  12 revised, 1 disputed (c006), 22 published; memo v1 fails Y2 (9 uncited
  sentences), as Y2's revision step expects.
- Remaining:
  1. Answer: Moure & Fu, published copy or arXiv v1? Then update the
     B1_CONTRACT.md Papers table to the new Zotero copies (paths and hashes
     in Literature/incoming/b1/zotero_pull.json; Lervag arXiv only; Li and
     Folch & Plapp published) and rerun A2.
  2. Dev run: the U3a gate opened before X1 ran (backlog item: make U3a
     require X1). Fix that, then resume: run_b1.py run --run
     runs/literature-review/b1/dev-20261101-kp (or the foreman); U3 with
     u3_review.py in a terminal of the user's own; then Y2, Y3, U4, BK.
  3. Open question: should U3 "missed statements" become claims through a
     fresh reader, E3 and E4? LaTeX-source comparison dropped (user).
  4. Step 6 build tests (T-det, T-units, T-applicability, T-perm,
     T-sandbox, T-trifecta); step 7 calibration diagnostic rerun (prompts
     changed form) and /lr-freeze; step 8 the nine-paper run with T-repeat.
  5. Leftover test folder runs/literature-review/b1/_review_test/ (two
     dummy files) for the user to delete.
- Bottlenecks: token budget before the 2026-11-01 deadlines (Lab Report 2,
  weekly homework); agent stations share the session's usage limit.

## 2026-11-01 — Prompt audit of instruction files, hooks, memory and worktrees

- Session: 8272b16e-7e7f-4c96-bb27-dd1f38bb621c
- Done: /claude-api prompt-audit of every instruction file that loads into
  sessions here, then of the hooks, the auto-memory index and the
  literature-review prompts. Fixes (11 commits, a7af16ed..0260e544):
  -t_out (no such option) -> -t_interv in memory/MEMORY.md and the
  monitoring.c comment; Proposals/ read only on request (MEMORY.md old line
  marked superseded; CLAUDE.local.md matches); morty-unit-reader lost the
  "zoom by reading again" line, units calibration C0_nozoom on run
  20261022-141826 passes as before (S3 8/8, S4 8/8, S4m 2/2, 0/4 control
  false alarms); hooks fixed for worktrees (ruff hook says "skipped" when
  venv_lunar is missing; session log always goes to the main checkout's
  sessions.tsv, both tested in a real worktree); flag_classifier.md removed
  (unused); skeptic_omissions.md takes {n_papers} from the claims (75 tests
  pass); permafrost references replaced where the code moved
  (src/permafrost2.c -> src/lunar_main.c, CLAUDE.md, ruff.toml); planning
  docs 16 and 18 copied back into docs/planning/ with a provenance note;
  scripts/new_worktree.sh (worktree from local HEAD or a named branch,
  linked venvs and CLAUDE.local.md, per-worktree make; tested create, lint,
  venv, remove) and .gitignore venv_*/ -> venv_*; effort decision added to
  memory/MEMORY.md. Not in git: CLAUDE.local.md (account entry removed,
  artifact spelling), auto-memory index cut from 78 to 41 lines (enceladus
  carry-overs and practice-round entries dropped from the index, files kept;
  push rule now "ask first"; seven files reworded off permafrost; backup of
  the old folder in the session scratchpad).
- Remaining: (1) CLAUDE.md "Sessions and housekeeping" still says parallel
  work uses `claude -w <name>`; propose pointing it at
  scripts/new_worktree.sh (plain `claude -w` branches from origin/main and
  misses unpushed commits). Needs user approval (CLAUDE.md edit).
  (2) Push: main is 13 commits ahead of origin; ask the user, then
  scripts/sync_course.sh. (3) -thin_iface_corr help text in
  src/lunar_main.c says "Default 0" but the code default is on; offer to fix
  the string. (4) docs/model_description.md and docs/material_parameters.md
  describe the retired three-phase model (gamma_iv/im/mv not in src/); mark
  or archive them. (5) ~/.claude/skills/.trash kept for now (user's call).
  (6) Lab Report 2 is written tonight; update the deadline line in
  CLAUDE.local.md once submitted.
- Bottlenecks: cp and rm are aliased with -i, so they stall in tool calls;
  use `command cp -f` / `command rm`. Reading runs/literature-review/*/untrusted
  is denied outside review mode, as designed.

## 2026-11-01 — CLAUDE.md review, linting, plugins, plot_fields fix

- Session: 8a61e2cf-a3f4-46e5-9347-3c087f8a51c9
- Done: CLAUDE.md gained commands, an architecture section and 8 fixes
  from an /init review (old test/ suite, HPC --out-root, literature-review
  agents and review mode, edits that need approval, memory format, tracked
  settings.local.json). CLAUDE_CODE_NEW_INIT=1 in user settings (mirrored
  in the committed copy). Linters: ruff 0.17.3 (ruff.toml, check only, now
  zero findings; legacy igakit folder excluded) with a report-only
  PostToolUse hook (.claude/hooks/ruff_check.sh); clang-tidy 23.1.4 from
  the pip wheel via scripts/lint_c.sh (zero findings); --compile-flags
  writes the git-ignored compile_flags.txt for clangd. Plugins installed at
  user scope: skill-creator, pyright-lsp (pyright 1.1.414 as a uv tool),
  math-proof; math-olympiad skipped (competition-tuned, auto-triggers on
  "verify this proof"). Solver lint: dead -outp test for SSA_evo.dat
  removed (every-step rows kept, as the user wants), strict strtod parsing
  of grain files, NOLINT on restart-log sscanf, seed cast. Before/after
  runs (channel_2D_H31um_eps0.86um, dtlad20_..._dt10tau, tags
  lintC_before/lintC_after): all outputs byte-identical. plot_fields.py,
  plot_fields_highres.py and scripts/plot_fields.py matched solver_evo.dat
  as a snapshot since 2026-10-01 (no pf.pvd, post-processing exit 1);
  pattern now sol_[0-9]*.dat. Pushed to origin and the course repo.
- Remaining: triage pyright's 145 errors / 118 warnings in active code
  (pyright from the project root; much is numpy/matplotlib typing noise).
  Run folders made before today keep the buggy staged plot_fields.py; to
  get pf.pvd, replace the folder's postprocess/ and rerun
  run_postprocess.sh. math-proof needs its env block in a dedicated
  per-problem folder (not global: it turns background tasks off) once
  model development starts. This log was last updated 2026-10-22; decide
  whether the B1 contract's Status section replaces it.
- Bottlenecks: none. Lessons: rm is aliased to rm -i and ls to eza (use
  `command rm`, `command ls`); run_postprocess.sh always uses the run
  folder's staged scripts.

## 2026-10-22 (5) — Units-verification workflow built and run; units.md verified

- Session: cd3a2de5-fd76-4f3c-8357-c1b9d1164930
- Done: built the workflow in three reviewed stages
  (scripts/units_verification/, 48 tests): S0 locate, S1 render, S1b Marker,
  S3 compare (units as dimensions; printed symbols), agents
  morty-unit-reader and beth-unit-checker (approved, merged, pushed), S4
  adjudication, S4m meaning review with a second pass on definition pages,
  C0 calibration (passed twice: 8/8, 8/8, 2/2; 0/4 false alarms on the
  final design), reader stability, U1 inspection set, S5 apply + recheck.
  Run 20261022-155335: user decided 22 rows; units.md now 54 verified
  literature rows, value column, UNRELIABLE header removed. Dropped j_T,
  j_v, t_n, T_cold, T_hot, beta_sol. Page images and Marker text stay
  local; DOIs and hashes make the run reproducible.
- Remaining: push the commits since the merge (asked); optional second
  pass for the solver's own parameters (M0, alpha) with a code-reading
  checker.
- Bottlenecks: none. Lesson: cp is aliased to cp -i; use `command cp -f`.

## 2026-10-22 (4) — Units workflow designed; flowchart agent; Marker; wrap-up

- Session: db191e62-cff4-48e1-a6d2-099d62fee556 (started before the
  SessionStart hook existed; ID from the transcript file name)
- Done: remaining hard-coded d0 values replaced by the shared function.
  Units-verification workflow designed (workflows/units-verification/:
  WORKFLOW.md, flowchart.json/.html). Flowchart renderer with overlap lint
  (scripts/flowchart/) and the morty-flowchart agent (approved). Sources
  for the units workflow copied from Zotero into Literature/; the two
  Moure & Fu SI copies shown to be the same document (one resaved by
  Preview). Marker 1.10.2 added (stage S1b) in venv_workflows;
  requirements-workflows.txt/.lock. Old remote branch deleted; wheeler_2010
  DOI confirmed by the user.
- Remaining:
  - Build the units-verification workflow in stages (S0/S1/S1b/S3 scripts
    with tests → agents + C0 calibration → real run + U1 inspection). User
    answers WORKFLOW.md §8 in the next prompt.
  - User decisions still open: delete the duplicate SI
    (Literature/cg4c00539_si_001.pdf)?; prune the repo-root
    PetIGA-3.20/.claude/settings.local.json as proposed?; regenerate the
    wedge/two-throat geometries (their EPS = 8.5840e-7 used the old d0,
    0.14 % large)?
  - Discuss the contact-angle free-flow condition (user, later).
  - Docs synthesis: deferred by the user (not a project priority).
- Bottlenecks: none.

## 2026-10-22 (3) — comp_eps d0 aligned; course repository mirrored

- Done: comp_eps.py computes d0 exactly as the solver (−0.14 %; golden
  values regenerated). Course repository now mirrors this project:
  old export history kept as tag archive/2026-10-18-export, main replaced
  by the subtree-split history (one approved force-push),
  feat/diag-air-film pushed (split from projects/sublimation_pf).
  scripts/sync_course.sh mirrors main and feat/* from now on.
- Remaining: user checks wheeler_2010 in library/references.bib; decide on
  the visibility of jackson-baglino/PetIGA-3.20 (personal repos cannot be
  "internal"); prune the repo-root settings.local.json; units verification
  workflow; docs-synthesis pilot design.
- Bottlenecks: none.

## 2026-10-22 (2) — d0 computed once; settings, library DOIs, units, memory split

- Done: d0 computed once at start-up (K&P Eq. 13) and used everywhere;
  -d0_sub0 removed from the solver and 106 opts files (b927b982). Simulation
  rule loosened to at most one run without instructions. settings.local.json
  pruned 197 → 21 rules, global settings copied in, committed. Removed a
  malformed rule from the repo-root settings.local.json (startup warning).
  library/references.bib: 85 DOIs, all resolve via doi.org; one title
  mismatch (wheeler_2010). Units table cleaned and flagged unreliable.
  Lunar auto-memory split into its own folder with a sectioned index.
  Contact-angle status corrected (static angle only).
- Remaining: decide how the course repo is synced (separate history);
  prune the repo-root settings.local.json; align comp_eps d0 with the
  solver (-0.14 %)?; units verification workflow; docs-synthesis pilot design.
- Bottlenecks: the course repo's history is unrelated to this repo's.

## 2026-10-22 — Project memory, permissions, skills, guides; branches cleaned

- Done: merged feat/d0-kaempfer-plapp (d0 from K&P Eq. 13), verified it
  compiles. Deleted the 7 practice-round/merged branches and all 7 sibling
  worktrees; feat/diag-air-film kept as a branch. Restored the comp_eps engine
  (broken by the restart) as preprocess/mesh_params with its 10 tests.
  Added memory/ (README decision table, MEMORY.md, USER.md,
  semantic/units.md from the archived units registry), PERMISSIONS.md,
  .claude/settings.json (pinned model, deny rules, session-ID hook), skills
  comp-eps and wrap-up, workflows/TEMPLATE.md, project/guides/01–05.
  Moved Library into library/ (ignored); archived the iCloud Planning, Review
  and Conventions folders.
- Remaining: user review of the auto-memory pruning list and the
  settings.local.json pruning; delete the old remote branch; design the docs
  synthesis pilot workflow.
- Bottlenecks: none.

## 2026-10-21 — Restart: clean slate for the agentic workflow

- Done: tagged `main` and all 8 branches under `archive/2026-10-21/`. Moved the
  practice-round material (agentic/, .agent_state*, .claude agents, settings
  and activity log, outputs/, the course/research symlinks, stale
  NOTES.md/TESTS.md and the three sediment-era docs) to
  `scratch/2026-10-21_practice_round/`. Rewrote CLAUDE.md and
  CLAUDE.local.md. Wrote `docs/contact_angle/README.md` (synthesis) and
  `project/PILLARS.md`. Stopped tracking LaTeX aux files.
- Remaining: user review of the new CLAUDE files; decide whether `project/` is
  committed; delete old branches and worktrees after confirmation; design the
  Pillar 1 workflow.
- Bottlenecks: none yet.
