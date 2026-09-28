# k_eff sintering campaign — to-do list

The working list for the manuscript runs. **Read this first in a new session.**
Tick items as they close, add the date and where the result lives, and commit.
Background and reasoning: `CAMPAIGN.md` (plan and results by stage),
`studies/rve_anisotropy/README.md` (domain size and packing bias),
`.claude/ACTIVITY_LOG.md` (what happened, session by session).

**Standing rules**
- Every manuscript run goes through `scripts/HPC/submit_keff_production.sh`.
  That script owns every option; stage files list runs only. Never hand-type
  `--extra-opts` for these runs.
- Submit in stages. A stage's check must pass before the next is submitted.
- Jackson launches every run; the assistant only prepares commands.
- Push before submitting. The script refuses to run from an uncommitted or
  unpushed checkout.

---

## Now: batch 3, the porosity × temperature matrix (125 runs)

- **Packings:** `inputs/packings/keff_LR40/`, φ {0.275, 0.325, 0.375, 0.425,
  0.475} × 5.
- **Temperatures:** T {−5, −10, −20, −30, −40} °C.
- **Index:** `batch3_phi_T.txt` (never submitted whole).

- [ ] **3a shakedown**: `batch3a_shakedown.txt`, 7 runs, ~$60–100.
  - Submit: `./scripts/HPC/submit_keff_production.sh studies/keff_sintering/batch3a_shakedown.txt`
  - Checks, all must pass:
    - [ ] all 7 reach 30 d; `health_check.py` clean (mass, φ bounds, k_eff KSP)
    - [ ] `k_eff.csv` is written (not `k_eff_tensor.csv`): one row per 5 steps,
      and every step at −40 °C
    - [ ] a `sol_*` at t ≥ 1 s exists (the movie opening frame), 50 log
      snapshots, ~9 GB/run
    - [ ] record s/step at 401 ranks against batch 2's ~20 s at 241, then
      decide: keep 60k DoF/core or revert to 100k (`scripts/lib/alloc.sh`)
    - [ ] −5 °C wall time comfortably inside 24 h
    - [ ] φ 0.475 runs stably; k_eff is monotone in φ
    - [ ] −5 / −20 / −40 on the 0.325 packing collapse on k(SSA) and on t/τ_sub
- [ ] **3b**: rest of −20 °C, `batch3b_T-20.txt`, 20 runs.
  - Check: seed scatter of the k_iso rise per φ (~3–4% expected at 0.325);
    the φ trend is larger than the scatter.
  - Check: k_xx/k_yy per φ follows the contact fabric (1.03 at 0.275 →
    ~1.1 at 0.425).
- [ ] **3c**: cold end, `batch3c_cold.txt`, 49 runs (−30, −40).
  - Check: every φ collapses onto its −20 °C curve. If one doesn't, the
    time-rescaling result depends on φ; rethink before buying the warm end.
- [ ] **3d**: −10 °C, `batch3d_T-10.txt`, 25 runs.
- [ ] **3e**: −5 °C, `batch3e_T-5.txt`, 24 runs. The most expensive stage;
  last.
- [ ] **After each stage:**
  - [ ] download the tables (see the `rsync` recipe in `ACTIVITY_LOG.md`, 2026-09-26)
  - [ ] run `postprocess/plot_keff.py` per run
  - [ ] run `postprocess/compare_keff.py` per batch
  - [ ] run `health_check.py`
  - [ ] record the measured cost here

## Open questions and follow-ups

- [ ] **k_eff solve cost.** 50–84 s per sample in batch 2, against 4 s on
  rev64 for the same solve at similar iteration counts, with 10× swings
  within a run. Looks environmental; unexplained.
  - Proposed: a small replay benchmark (tolerance 1e-7, GAMG threshold,
    hypre, pc_freeze) with `-log_view`. A few dollars, no new simulation.
  - Awaiting a go-ahead.
- [ ] **rev64 seeds 2–4** (2026-09-21, jobs 3264515–3264517). Are they on
  scratch? If so, fetch `k_eff.csv` + `SSA_evo.dat`: they test the
  sintered-state RVE directly.
- [ ] **Seam relaxation reach.** It moves 86 grains by more than 0.1 R (up
  to 1.2 R), reaching ~30% of L from the seam; seen in the deposition movie.
  - The generator's comment says only grains near the seam move. Correct the
    comment; decide whether it matters.
- [ ] **φ 0.475 connectivity at the band.** The generator's percolation test
  uses the sharp geometry, but the solver joins gaps smaller than 9.2·eps.
  Measure connectivity at the band before writing about 0.475.
- [ ] **Kink near day 7.5** in seed 1, −20 °C (batch 2): the time step is
  flat there, so it's likely a topology event. Look at the snapshots.
- [ ] **Anisotropy direction vs real snow.** Our packings give k_xx > k_yy.
  - Check against measured snow anisotropy before claiming anything.
  - Cite only closely matching work, with DOIs verified via Crossref.
- [ ] Correct `studies/packing_design/README.md`: "vertical load paths"
  predicted k_yy > k_xx; the data say the opposite.
- [ ] `effective_thermal_cond/docs/tensor_conductivity_law.tex`: cite Nicoli,
  Plapp & Henry 2011 (prior art for the tensor law).
- [ ] Optional: an SSA-triggered k_eff cadence (sample when SSA drops by δ),
  giving evenly spaced points on the k–SSA curve. About 30 lines in
  `KeffDue` (`src/keff_sample.c`). Not needed for batch 3.

## Later (CAMPAIGN.md stages 6–7)

- [ ] Sensitivity arms from the stage-5 centre point, 3 seeds each:
  - α_c {1e-4, 1e-3, 1e-2}: may change the kinetic regime, and with it the
    T collapse;
  - R_ave ×0.5 / ×2 at fixed L/R;
  - σ_ln 0.2;
  - eps ×2.
- [ ] Analysis:
  - partial correlation of k_eff with SSA at fixed porosity;
  - anisotropy from seed-mean k_xx/k_yy (never per seed);
  - per-snapshot metrics (Euler characteristic, chord lengths, D_eff);
  - lever ranking with seed error bars.

## Manuscript obligations (things the text must say)

- [ ] **Necks are under-resolved**, by design: floor r/R = √(12·eps/R) ≈ 0.49.
  No neck-growth exponent is claimed.
- [ ] **In this model, temperature is a pure time rescaling by τ_sub**
  (attachment-limited, α_c = 1e-3). State it as a model property, not a
  result about snow.
- [ ] **Anisotropy comes from the drop-and-roll deposition rule.** Report it
  as a seed mean; k_xy is a zero-mean fluctuation that shrinks as 1/L.
- [ ] **The y-seam gate** (≥ 0.76 of interior contact density) and its
  residual bound: k_yy at most ~1–1.5% low at L/R 40.
- [ ] **φ 0.475 is past 2D solid percolation**: it marks the limit of a 2D
  snow analogue.
- [ ] **The baseline is t = 1 d**, not t = 0 (IC relaxation). Compare
  temperatures at matched t/τ_sub or matched SSA.
- [ ] **Pore connectivity:** fragmented at t = 0, D_eff ≤ 0.0014 D_v.
  Long-range ripening is suppressed in 2D.

## Done (most recent first)

- 2026-09-28 — Production submit script `scripts/HPC/submit_keff_production.sh`;
  this list.
- 2026-09-28 — Batch 3 staged into 3a–3e; 125 opts generated.
- 2026-09-28 — Deposition movie + storyboard
  (`postprocess/make_deposition_movie.py`); outputs in
  `~/SimulationResults/HPC_results/enceladus_DSM/GrainPackingSintering/presentation/`.
- 2026-09-28 — Packings: φ 0.275–0.475 by 0.05 × 5, unique seeds, y-seam
  gate, porosity-salted RNG.
- 2026-09-28 — RVE and anisotropy check from existing data
  (`studies/rve_anisotropy/`): L/R 40 is adequate; the seam is the only
  generator artifact.
- 2026-09-27 — Batch 2 (warm end, 9 runs) analysed: the dt check passes, T is
  a time rescaling, it cost $342; k_eff every 5 steps is enough.
- 2026-09-26 — k_eff plots (absolute, normalized, off-diagonal), compare
  plots; `k_eff.csv` naming fix; t ≥ 1 s opening frame.
