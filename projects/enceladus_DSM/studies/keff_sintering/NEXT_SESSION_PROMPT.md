# Prompt for the next session

Copy everything below the line into a fresh Claude Code session started in
`/Users/jacksonbaglino/PetIGA-3.20/projects/enceladus_DSM`.

---

I'm starting the production campaign for the sintering → effective thermal
conductivity study. Read `studies/keff_sintering/audit_unresolved/README.md`
first — it is the audit that killed the previous attempt and it defines what
this campaign has to do differently.

## Where things stand

**The old results are dead.** Nine 28-day runs in
`~/SimulationResults/HPC_results/enceladus_DSM/unresolved_results/` are not
usable. `beta_sub0` was hardcoded at `1.4e5` for every temperature, which makes
the implied `alpha_c` *rise* 8.4× as `rho_vs` falls 8.0× — their product is flat
to 4.2%, so the "k_eff evolution is temperature-independent" result was forced by
construction. Don't cite them. The porosity sweep survives only as a qualitative
"the trend is monotone and steep."

**Decisions already made — do not relitigate these:**

- `alpha_c = 1.0e-3`, **constant at every temperature.** The anchored Arrhenius
  `alpha_c(T)` is abandoned: it is too restrictive to implement and we showed it
  doesn't do what was intended. If a constant is an approximation, it is at
  least the *same* approximation at every point of the sweep, which is exactly
  what the old runs got wrong. Pass `--alpha-c 1.0e-3` to
  `generate_study_opts.py`; do not pass `--alpha_arrhenius` to `comp_eps.py`.
- **`alpha_c` does not affect the mesh.** With `--vn_feature`, `comp_eps.py`
  derives `v_n` from `R_feat`, which makes the K&P Eq.(45) kinetic bound exactly
  equal `R_feat`. It is binding at every `alpha_c` from 1e-5 to 3e-2, and
  `eps = safety·R_feat` regardless. Mesh cost is set by `R_feat` alone.
- **`dtmax` is derived per temperature** as `1.09·tau_sub`, never hand-set and
  never tied to the output cadence. At `eps = 4.5e-7`, `alpha_c = 1e-3`:
  1753 s (−20 °C) → 14692 s (−40 °C). `generate_study_opts.py` does this already.
- **`-eps_valid_temp`** must be set on every run. One ε per temperature. This is
  the guard that would have caught the old bug.
- Geometry: `phi = 0.325`, `R_ave = 50 µm`, `L/R_ave = 40`, fully periodic (`xy`).
  Packings already exist: `inputs/packings/pilot_LR40/` (4 seeds) and
  `inputs/packings/rev_LR64/` (4 seeds).
- REV: `L/R_ave ≥ 40`, and it must hold at **`t_final`**, not `t = 0` —
  coarsening makes it fall during the run.

## Two things are unresolved and one of them gates the campaign

1. **GATING — the `k_eff` interpolation coefficient.** `KeffPointCond` in
   `src/keff_cell.c:33` uses an arithmetic mix, which gives a +17.8% rise where
   the sharp-threshold and tensorial forms both give +57% / +54.6% (they agree
   with each other to 0.7%, and with the old standalone code's ~+50%). That is a
   **3× correction on the headline quantity.** Fix this and re-verify with
   `-keff_replay` on existing snapshots *before* any production run is
   submitted. Do not spend money producing numbers in the wrong coefficient.
2. **Not gating — `R_feat/R_ave` = 1/25 vs 1/50.** 1/25 is the coarsest
   defensible mesh (`Nx = 2829`, ~24M DOF, ~301 cores). The de-risking plan was
   to run the first condition at both, 3 seeds each, ~$100–150. Decide whether
   that is worth it given (1) has to happen first anyway.

## What I want from you, in this order

**Step 1 — fix the gating item.** Change the `k_eff` coefficient to
sharp/tensorial, verify via `-keff_replay` against snapshots we already have,
and show me the before/after on a run whose old number I can compare to.

**Step 2 — generate the first batch and hand me the command.** A *small* first
batch: I want 3–4 runs, not the campaign. My instinct is the temperature axis at
fixed packing (one seed), since that is the axis the old work destroyed — but
argue for something else if there's a better first cut. Build it through
`preprocess/generate_study_opts.py --alpha-c 1.0e-3` and hand me a
`scripts/HPC/submit_batch.sh` invocation. Tell me the expected cost and
wall-clock before I submit.

**Step 3 — write the campaign plan** to `studies/keff_sintering/CAMPAIGN.md`,
with explicit **checkpoints**: after each batch, what number do we read, and what
does it have to show for the next batch to be worth submitting? These runs are
expensive; I want to stop early if the signal isn't there. Include a cost
estimate per stage and a running total. Also flag which stages produce a figure
my advisor can look at.

## How I work

- **Never launch a simulation yourself.** Give me the command; I run it.
- **Two or more HPC jobs go through `scripts/HPC/submit_batch.sh`** — chained
  single submits race in `obj/` and die with "Stale file handle".
- **Push to GitHub after committing**, before I submit anything, so I can pull on
  the cluster.
- Always go through `scripts/Studio/run_enceladus.sh` or
  `scripts/HPC/submit_enceladus.sh`, never bare `./enceladus_dsm` or `mpiexec`.
- Results must land in `$RESULTS_BASE` or in the repo under
  `studies/<study>/`, and must come with plots. Say where they are.
- Commit autonomously as you go; never `--no-verify`, never force-push `main`.
- Before writing any new input file, `find inputs/scratch -name '*<something>*'`
  and `git mv` it back out if it already exists.

Start with Step 1.
