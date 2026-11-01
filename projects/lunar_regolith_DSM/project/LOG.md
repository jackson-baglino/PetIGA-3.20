# Running log

Newest first. One entry per session: done, remaining, bottlenecks.

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
