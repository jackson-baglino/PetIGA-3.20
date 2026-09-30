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

- [x] **3a shakedown** — ran 2026-09-28 (commit f4179a9, 401 ranks), analysed
  2026-09-29. Results: `~/SimulationResults/HPC_results/enceladus_DSM/GrainPackingSintering/keff_b3a_shakedown/`
  (`compare/` holds the overlays). Cost **$195** (estimated $60–100).
  - [x] all 7 reach 30 d; health clean (mass drift ≤ 5e-7, 0 KSP failures,
    symmetry ≤ 1e-11, 0 bound trips)
  - [x] `k_eff.csv` every 5 steps; every step at −40 °C (108 samples)
  - [x] opening frame at t = 1.26 s; 45 snapshots, 8.1 GB/run
  - [x] **s/step at 401 ranks is no faster** (15–28 s vs ~20 at 241), and the
    k_eff solve is SLOWER: 0.85 s/iteration vs 0.49 at 241 ranks. So 60k
    DoF/core costs 1.66× for nothing.
    → **DECISION PENDING: revert `scripts/lib/alloc.sh` to 100k**
  - [x] −5 °C used 16.6 h of 24 (≈11–13 h expected at 241 ranks)
  - [x] φ 0.475 stable. k_iso falls monotonically with φ (0.86 → 0.31 W/m/K at
    30 d). Rise 1→30 d: +25.9 / +25.0 / +24.7 / +22.6 / +19.7 %.
  - [x] T collapse on the 0.325 packing: k_iso at matched SSA agrees to 0.06%;
    the speed-up equals the τ_sub ratio to ±1.5% (−5: 3.83×; −40: 0.128×)
  - Finding: at φ 0.475 (ice percolates in y only) k_xx rises +6% against
    +26% for k_yy; k_xx/k_yy goes 0.50 → 0.37.
- [x] **k_eff cadence → SSA trigger** (2026-09-29).
  - The 3a curves had visible corners at every 5 steps (0.49 d = 5.4 τ_sub
    apart). The new setting samples every 0.1% drop in SSA after 11 τ_sub,
    every 5 steps before that, with a 20 τ_sub backstop.
  - The same option at every temperature; the −40 rule is gone.
  - Max error 0.02% of the plotted range, 60× below the smallest real kink
    (the SSA ≈ 15,300 event on seed 301, ≥ 1.3%).
  - Samples per run ~308 / 179 / 46 at −5 / −20 / −40.
  - Code: `src/keff_sample.c` (KeffDue). `studies/keff_sintering/predict_cadence.py`
    mirrors it and reproduced all 7 3a schedules on the old path; run
    `--check` on the rerun.
- [ ] **Rerun 3a** with the final options once the two tests below are done.
  The first 3a is kept as the shakedown record, not pooled.
- [ ] **Before 3b** (cost, not correctness), in this order:
  - [x] **DoF/core scaling test**, ran 2026-09-29 (`scaling/README.md`).
    - The fewest ranks is cheapest: 198k DoF/core (121 ranks) is 3–5×
      cheaper than 60k.
    - k_eff is 15–30× slower at ≥ 8 nodes, from communication; the suspect
      is the mixed node-type constraint.
    - [ ] follow-up (~$5): bracketed `--constraint` at 241 ranks; 300k/400k
      targets.
    - [ ] set `scripts/lib/alloc.sh` (Jackson's call; 200k recommended now)
  - [ ] k_eff solver benchmark: `solver_benchmark/` (~$10, replay only).
    Now LOW value: at ≤ 6 nodes a sample is ~5 s, ~15 min per −20 °C run.
    Optional. Adopt a faster setting only if max |Δk|/k_iso < 1e-5, and record
    it in the production script before 3b.
  - Projected cost of all 125 runs at 200k DoF/core with the SSA cadence:
    roughly $0.9k. That's per-temperature $/run from `scaling/README.md`
    (−30/−40 °C are cheaper still) × 25 runs each.
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
- [ ] Figures normalized by k_eff(t = 0), e.g. the AGU snapshot set, are a
  post-processing choice (user, 2026-09-29), not a problem.
  - k(0) is the unrelaxed analytic IC: the first ~3.7 τ_sub are relaxation,
    worth ~+35–50% in k.
  - Say in the caption that the normalization includes that relaxation.
  - Before comparing porosities that way, check whether k(11 τ)/k(0) differs
    with φ.
- [ ] HPC support: job 3606286 (3a, −5 °C) stuck in COMPLETING on hpc-34-37
  since ~06:30 2026-09-29; billing should have stopped at job end (check `sacct`).
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

- 2026-09-30 — One normalization for every k_eff figure: `plot_keff.py`,
  `compare_keff.py` and `plot_keff_snapshots.py` open on the first k_eff
  sample with t >= 1 s, call it t = 0, and divide by its measured values
  (k_eff,0, SSA_0). The 11 τ_sub interpolated baseline ("b") and the greyed
  relaxation are gone; seed means cover only the span every seed does.
  Opening at 1.7 s vs the snapshots' 126 s differs by 0.12 % in k_iso.
  Colours: single runs cmocean `deep` (k_xx dashed, k_yy dotted, k_iso
  solid), temperature `thermal`, porosity amp-to-black.
  `coefficient_fix/compare_laws.py` keeps its own `--baseline-days` (a study
  record, not a figure script).

- 2026-09-29 — 3a analysed. Fixed on the way:
  - `health_check.py` keyed runs by seed, so it skipped the same packing at
    other temperatures;
  - `compare_keff.py` and `plot_keff.py` baselines are now 11 τ_sub, and
    interpolated. With −40 °C present, "1 d" was 1.4 τ_sub, inside the
    relaxation;
  - `run_batch_postprocess.sh` takes POSTPROCESS_DIR/PYTHON overrides.

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
