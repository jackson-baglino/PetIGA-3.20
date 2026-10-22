# Running log

Newest first. One entry per session: done, remaining, bottlenecks.

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
