# k_eff sintering campaign — record of work

Everything done and still planned, in order. Each finished item gives where its
data lives and one sentence on what it was and what came out of it. Written for
the manuscript: the methods and supplement should be reconstructible from this
page alone.

**Keep it current.** Update this file whenever a task finishes or a plan
changes, in the same commit. `TODO.md` stays the working checklist (checks,
thresholds, open questions); this page is the history.

**Path shorthand**

| | path |
|---|---|
| `GRP` | `/resnick/groups/rubyfu/jbaglino/simulation_outputs/enceladus_DSM` (permanent group storage) |
| `CAMP` | `GRP/keff_sintering_campaign` (every manuscript run since 2026-10-01; one folder per run `<geom>__<exp>/`, submission records in `stages/`) |
| `SCR` | `/resnick/scratch/jbaglino/enceladus_DSM` (**purged** — copy anything still needed to `GRP`) |
| `LOC` | `~/SimulationResults/HPC_results/enceladus_DSM` (local copies; `LOC/GPS` = `LOC/GrainPackingSintering`) |
| `REPO` | `projects/enceladus_DSM` |

---

## Done

### Model and method

- [x] **2026-09-22 — Audit of the 2025-09 temperature sweep.**
  Data: `LOC/unresolved_results/` (local only). Audit: `REPO/studies/keff_sintering/audit_unresolved/`.
  Nine old runs re-examined from their staged source. **Finding:** a hard-coded `beta_sub0` made the implied α_c rise 8.4× as ρ_vs fell 8×, so their "k_eff is temperature-independent" result was forced by construction. The runs were discarded and α_c is held constant at 1e-3 from here on.

- [x] **2026-09-22/23 — Tensor conductivity law (`-keff_interp tensor`).**
  Verification: `REPO/studies/keff_sharp_limit/verification/` (planar slab) and `REPO/studies/keff_sharp_limit/disk/` (cylinder array). Write-up: `effective_thermal_cond/docs/tensor_conductivity_law.tex`.
  The law is arithmetic along the interface and harmonic across it. **Finding:** it removes the diffuse band's O(eps) bias (slab error 59% → 4e-7; cylinder +1.2% → +0.006%). It is the default from 2026-09-23. Prior art: Nicoli, Plapp & Henry 2011, still to be cited in the .tex.

- [x] **2026-09-23 — Periodic IC stretch fixed.** Code only (`FillIC1D/2D`).
  The benchmark ICs were stretched by (N+1)/N on periodic axes. The packings were never affected.

- [x] **2026-09-24 — Disk f-sweep (36 runs).**
  HPC: `SCR/batch_2026-09-24__14.23.15_keff_disk_fsweep/` (location not verified). Local: `LOC/GPS/batch_2026-09-24__14.23.15_keff_disk_fsweep/`. Analysis: `REPO/studies/keff_sharp_limit/disk/analyze_fsweep.py`, `fsweep.png`.
  Single ice cylinders, 6 radii × 3 eps/R × 2 laws. **Finding:** the arithmetic-law error grows with solid fraction (+1.9% at f = 0.05 → +38.6% at f = 0.50), while the tensor law stays under 1%.

- [x] **2026-10-01 — Grain-size scaling (no runs).** `REPO/studies/grain_size_scaling/` (script, CSV, figure, README).
  An R_ave sweep from 0.05 to 500 µm through the campaign's sizing code. **Finding:** the mesh is scale-free (Nx = 2829 at every size, since eps = R/50 and L = 40R), and so is the cost (~$3–4 per run) while the run stays attachment-limited (R ≪ L* = 139 µm). Simulated time scales as R² (30 d → 7.6 s at 0.1 µm). Below ~1 µm the model's physics fails (Knudsen vapour transport and gas conduction, dominant surface diffusion), not its cost.

### Pilot, replay and REV checks (φ 0.325, −20 °C)

- [x] **2026-09-16 — Pilot: 4 seeds, L/R 40, 30 d, arithmetic law.**
  HPC: `SCR/batch_2026-09-16__13.16.11_pilot_keff/`, plus the resume legs in `SCR/packing_2D_pilot_phi0.325_…_seed<N>_…/`. Local: `LOC/GPS/batch_2026-09-16__13.16.11_pilot_keff/`.
  First full sintering → k_eff trajectories. **Findings:**
  - Leg 1 ran at the old flat `-dtmax 200`, which inflated the step count about 18×. This led to the derived `dtmax = 1.09·τ_sub` and to the corrected cost model.
  - Seed 2 lost ranks to a node fault (`Slurmd could not connect IO`) and has no 1-day baseline.
  - The data exist only on scratch.

- [x] **2026-09-18/21 — rev64: L/R 64 REV check, 4 seeds.**
  HPC: `SCR/packing_2D_rev_phi0.325_Rave50um_LR64_seed<N>_L3.2mm_…/` (jobs 3015847, 3264514–3264517). Local: `LOC/GPS/rev64/` (all 4 seeds, with `k_eff.csv`) and `LOC/rev64/`.
  These were the first runs on the derived dtmax throughout. **Findings:**
  - 30 d took 371 steps, 2 h 19 min at 615 ranks, ≈ $17: the basis for the corrected cost model.
  - A k_eff sample took 4 s.
  - Seeds 2–4 have not been analysed yet (open item).

- [x] **2026-09-24 — Stage 1: pilot replay under the tensor law (8 jobs, $43).**
  HPC: `SCR/batch_2026-09-24__12.18.18_keff_replay/` (location not verified); the tensor CSVs were merged into the pilot run folders. Local: `LOC/GPS/batch_2026-09-24__12.18.18_keff_replay/`. Analysis: `REPO/studies/keff_sintering/coefficient_fix/`.
  The stored pilot snapshots were re-solved for k_eff; no new simulation. **Finding:** the sintering-driven k_eff rise from day 1 goes from +15.8% (arithmetic) to +28.8% (tensor), ×1.83. The arithmetic law read k_eff(0) 71% high on tangent grains.

- [x] **2026-09-28 — RVE and anisotropy check from existing data** (no new runs).
  Results: `REPO/studies/rve_anisotropy/` (README, CSVs, figures).
  Contact fabric and y-seam analysis of all 42 packings on disk. **Findings:**
  - L/R 40 is adequate for k_iso and SSA.
  - k_xx > k_yy (≈ 1.09) comes from the drop-and-roll rule, not the domain size. The `packing_design` README has it backwards.
  - k_xy is a zero-mean finite-size fluctuation.
  - The y-seam is a contact-poor generator artifact. This led to the seam gate, and the 3×3-superdomain idea was rejected.

### Warm-end temperature batch

- [x] **2026-09-25/27 — Batch 2: φ 0.325, seeds 1/3/4 × −20/−15/−10 °C, 9 runs, k_eff every step, $342.**
  HPC: `GRP/batch_2026-09-25__18.06.50_keff_T_warm/`. Local: `LOC/GPS/keff_T_warm_phi0.325/` (overlays in `compare/`).
  **Findings:**
  - The time-step check passes: the dtmax-200 pilot and the dtmax-8526 rerun agree to 0.01%.
  - **Temperature is a pure time rescaling** in this attachment-limited model. k(SSA) is the same at every T to < 0.5%, and the rate equals the τ_sub ratio.
  - Sampling k_eff at every step was the cost driver (50–84 s per sample).

### Packings

- [x] **2026-09-28 — Seam gate and porosity-salted seeds in the generator.** Code: `REPO/preprocess/generate_packing_gravity.py`.
  `--min-seam-contact 0.76` rejects packings with a contact-poor y-seam. Seeds are salted by porosity so streams are independent.

- [x] **2026-09-28 — First production packings (gated), φ 0.275–0.475 × 5.**
  `REPO/inputs/packings/keff_LR40_gated/` (seeds 301–905).
  Built with every homogeneity gate on. **Superseded** on 2026-10-01 (option B, below). Kept because the first 3a ran on them and batch_rve measures the gate effect with them.

- [x] **2026-09-28 — Deposition movie and storyboard.**
  Local: `LOC/GPS/presentation/`. Script: `REPO/postprocess/make_deposition_movie.py`.
  A replay of one packing's build, verified against `grains.dat` to 5e-16 m. **Finding:** the seam relaxation moves grains up to 1.2 R, reaching about 30% of L (open item: correct the generator's comment).

- [x] **2026-09-30 — Convergence-study packings, L/R 20–80, ungated.**
  `REPO/inputs/packings/rve_phi0.325/` (26 packings). Script: `REPO/studies/keff_sintering/make_rve_packings.sh`.
  **Finding:** the first build switched void gates with size and confounded the two (z_band jumped at the switch). The ungated rebuild shows z_band converged from L/R 30 (3.47–3.53), and that the gated production set sat about 5% low (3.32 vs 3.48).

- [x] **2026-10-01 — Production packings rebuilt ungated (option B), φ 0.275–0.475 × 5.**
  `REPO/inputs/packings/keff_LR40/` (seeds 1601–2005; `packings_summary.csv`). Script: `REPO/studies/keff_sintering/make_packings.sh`.
  Only the seam and percolation gates are on (percolation is off at 0.475). **Findings:**
  - z_band at φ 0.325 is 3.46, matching the converged ungated value.
  - Seed-mean z_band rose 0.14–0.20 at φ ≤ 0.375 and was unchanged at 0.425/0.475.
  - φ 0.475 now percolates in both axes in 4 of 5 packings (1 of 5 gated): the old gates also favoured disconnected packings there.

### Shakedown, cadence and allocation

- [x] **2026-09-28/29 — Batch 3a shakedown (first): 7 runs on the gated packings, $195.**
  HPC: `GRP/batch_2026-09-28__13.50.25_keff_3a_shakedown/`. Local: `LOC/GPS/keff_b3a_shakedown/` (overlays in `compare/`).
  One packing per φ at −20 °C, plus φ 0.325 at −5 and −40 °C. **Findings:**
  - Every correctness check passed: health, opening frame, mass.
  - The T collapse holds to 0.06%.
  - k_iso falls monotonically with φ (0.86 → 0.31 W/m/K at 30 d).
  - At φ 0.475 (ice percolating in y only) k_xx rose +6% against +26% for k_yy.
  - Every-5-step k_eff sampling gave visible corners, which led to the SSA trigger.
  - 401 ranks was no faster, which led to the scaling test.
  - Kept as the shakedown record; not pooled.

- [x] **2026-09-29 — SSA-triggered k_eff sampling.** Code: `REPO/src/keff_sample.c`. Mirror and checker: `REPO/studies/keff_sintering/predict_cadence.py`.
  k_eff is sampled at every 0.1% drop in SSA after 11 τ_sub, every 5 steps before that, and at least every 20 τ_sub. **Finding:** the maximum error is 0.02% of the plotted range, 60× below the smallest real kink, with about 308/179/46 samples per run at −5/−20/−40 °C.

- [x] **2026-09-29/10-01 — DoF/core scaling tests.**
  HPC: test 1 `GRP/batch_2026-09-29__*_scaling_*/`; test 2 `GRP/batch_2026-09-30__*_scaling_*/`. Local: `LOC/batch_2026-09-29__scaling/`. Results: `REPO/studies/keff_sintering/scaling/README.md`, `scaling.png`.
  Restarts of a 3a run from day 10 at 20k–400k DoF/core. **Findings:**
  - The phase-field step does not strong-scale, so fewer ranks are cheaper.
  - The k_eff solve is 15–30× slower at ≥ 8 nodes (communication).
  - Cost is flat (±10%) from 61 to 161 ranks.
  - **Chosen: 200k DoF/core = 121 ranks = 4 nodes** (≈ $3.85 per −20 °C run, ≈ $12 per −5 °C run).
  - Jobs are capped at 6 nodes (`MAX_NODES_PER_JOB`).

- [x] **2026-09-30 — Memory sizing.** Code: `REPO/scripts/lib/alloc.sh` `mem_per_cpu`.
  The 61-rank jobs were OOM-killed at 1 GB per core. The fitted model is peak = 0.20 GB + 2.35 GB/MDoF per rank + 8 B × total DoF, requested at 1.5× that. That gives 2G for production runs and 4G for L/R 80. Memory is not billed.

### Infrastructure

- [x] **2026-09-28 — Production submit script and stage files.** `REPO/scripts/HPC/submit_keff_production.sh`, `REPO/studies/keff_sintering/batch3*.txt`, `batch_rve.txt`.
  One script owns every run option; stage files list runs only. It refuses an uncommitted or unpushed tree.
- [x] **2026-10-01 — Per-run time limits.** `submit_keff_production.sh` (`time_limit_for`), `submit_batch.sh` (per-job `--time`).
  The jobs had sat on Priority for over 24 h. **Finding:** `sshare` showed `rubyfu` at 1.5× its fair share (`LevelFS` 0.66), 98% of it this campaign's usage. Jobs now ask about 2× their predicted wall time instead of 24 h, so the scheduler can backfill them: −5 18 h, −10 12 h, −20 6 h, −30/−40 4 h, scaled up for the L/R 56/80 meshes (8 h / 16 h).
- [x] **2026-10-01 — One campaign folder and per-stage download.** `CAMP/`; `REPO/scripts/HPC/fetch_stage.sh <stage file>` downloads a stage in one rsync.
- [x] **2026-09-26/30 — Figure pipeline.**
  - `plot_keff.py` and `compare_keff.py` make the per-run and per-batch k_eff plots.
  - `plot_keff_snapshots.py` makes the manuscript snapshot figures.
  - Movies open on the t ≥ 1 s frame.
  - One normalization everywhere: k_eff,0 and SSA_0 at the first sample with t ≥ 1 s.

- [x] **2026-10-02 — Timestep ceiling raised from 1.09·τ_sub to 2·τ_sub.**
  Data: `~/SimulationResults/lunar_regolith_DSM/scratch/batch_2026-10-01__{14.25.24_dtmax_small,14.52.38_dtmax_sinter,15.49.53_dtmax_ripen}` (local only). Tables: `projects/lunar_regolith_DSM/studies/contact_angle/velocity_plan_2026-10-01/README.md`.
  A six-rung ladder (0.8–40·τ_sub) in the lunar solver (same `src/`) on channel growth, a sintering pair and a ripening pair, all at −20 °C, α_c = 1e-3. **Finding:** 2·τ_sub agrees with 0.8·τ_sub to ≤ 0.1 % in every measure; 5 is within 1 %; 10 is 2–4 % slow; 20–40 are 10–30 % slow, always stable. `snow_T*_30d.opts` and `generate_study_opts.py` now use 2·τ_sub. **Caveat:** runs finished before this date used 1.09·τ_sub; the ladder says the two are indistinguishable but it was not run on a packing. The early k_eff cadence (`-keff_freq 5`, every 5 steps) is step-based, so it samples about half as often in time before 11·τ_sub.

- [x] **2026-10-04 — dtmax check on a packing: PASS.** `REPO/studies/keff_sintering/dtmax_check/`. The gated seed301 packing at −20 °C, run at 1.09·τ_sub (first shakedown) and at 2·τ_sub (3a rerun, job 3774294). **Finding:** SSA(t) agrees to within 0.21%, k_iso(t) to 0.13% and k_iso(SSA) to 0.08%, in 240 steps instead of 367. All 2·τ_sub runs stand.

- [x] **2026-10-04 — dtmax capped at 0.2 d for −30/−40 °C.** `REPO/preprocess/generate_study_opts.py` (`DTMAX_CAP_S`), `snow_T-30/-40_h1.00_30d.opts`.
  The 3a −40 °C k_eff(t) curve had visible corners. Two causes: k_eff every 5 steps before 11 τ_sub (= 7.7 d at −40 °C), and steps of about 1.1 d (2·τ_sub = 1.4 d). At −20 °C steps are 0.18 d. **Change:** dtmax = min(2·τ_sub, 0.2 d), which binds only at −30/−40 °C and gives them the −20 °C time resolution, for about +$1 per run. Seed1701 at −40 °C is rerun in 3c, and its 3a copy is moved to `superseded_*`.

- [x] **2026-10-04 — Porosity-series figures.** `REPO/studies/keff_sintering/phi_summary.py` → `LOC/keff_sintering_campaign/compare/phi_summary/T-20/` (`phi_trends.png`, `seeds_by_phi.png`, `anisotropy_time.png`, `phi_trends.csv`).
  **Finding:** k_iso scales as SSA^−0.80 to SSA^−0.83 for φ ≤ 0.325 after 11 τ_sub; the exponent weakens to −0.67 at φ 0.475. The anisotropy is set almost entirely by the initial packing: it shifts once during the first day, then drifts only slowly. Seed1601 (φ 0.275) is a low outlier in absolute k. Also added `rve_convergence/kxy_vs_L.py`: k_xy is a zero-mean fluctuation whose RMS falls from 7.9% at L/R 20 to 3.1% at L/R 80.

- [x] **2026-10-04 — Four snapshot PNGs for every run.** `REPO/scripts/lib/select_snapshots.sh` (each job marks t = 0, t_final/3, 2 t_final/3 and the last snapshot in `.rsync-snapshots`), `fetch_stage.sh` (downloads those four `sol_*.dat` by default; `--no-snapshots`), `REPO/postprocess/render_snapshots.py` (renders locally after the fetch into `<run>/plots/snapshots/`; `--no-render`).
  The PNGs are about 3080 px at one pixel per element, 2–3 MB each (L/R 80: ~11 MB). The download is about 0.8 GB per L/R 40 run. Rendering on the HPC was considered: it would mean a 1-core follow-on job, since SLURM cannot release a running job's cores, or holding 121 cores for a minute in-job (3–5% of a run). Local rendering was chosen because it needs no HPC Python.

### Molaro 2019 grain-pair validation (manuscript Fig. 2; separate study)

- [x] −20 °C round 2: local `LOC/GrainPairSintering/batch_2026-09-08__17.20.46_molaro_T-20_round2/`. −5 °C: `LOC/GrainPairSintering/batch_2026-09-29__10.01.06_molaro_T-5_round2/` and `LOC/GrainPairSintering/molaro_2D_…_T-5pair_…_h0.99674_2h_…/`. HPC: scratch, exact path not recorded. Study: `REPO/studies/molaro_2019/`, figures in `studies/molaro_2019/manuscript/`.
  Two grains sintering, compared with the Molaro cryostage data. **Finding:** the model reproduces the neck-growth exponent and about 50% of the rate (the vapour share); the remainder is surface diffusion, which the model doesn't include.
  - A −5 °C run timed out and exposed a destructive restart, which is now fixed (resume in place).

---

## In progress

- [ ] **Batch 3a rerun** — first submitted 2026-10-01 at commit d7452d5 and cancelled while pending, because dtmax changed to 2·τ_sub on 2026-10-02. Resubmitted at the 2·τ_sub commit.
  - **What:** 8 runs on the option-B packings: one per φ at −20 °C (seeds 1601/1701/1801/1901/2001), plus 1701 at −5 and −40 °C, plus gated seed301 at −20 °C as the paired dtmax check against the first shakedown (1.09·τ_sub, same packing). Final options, 121 ranks, per-T time limits.
  - **Data:** HPC `CAMP/<geom>__<exp>/`, record in `CAMP/stages/batch3a_shakedown__<ts>/`. Local `LOC/keff_sintering_campaign/` (`fetch_stage.sh studies/keff_sintering/batch3a_shakedown.txt`).
  - **Check:** the 3a checks, plus `predict_cadence.py --check` on every run.
  - *Result (2026-10-04, all 8):* all reached 30 d (to within the final-step quirk below); health clean; all k_eff schedules match. dtmax pair PASS (see the dtmax check entry). φ 0.475 is stable. The seed1701 −5 °C run fits the T collapse: k at matched SSA agrees to 0.04%, with a time ratio of 3.85 against the τ_sub ratio of 3.79; 4 h 43 min, $6.85. The 3a + 3b + convergence round cost $164 in total.
  - *Partial result (2026-10-03, 3 of 8: φ 0.275 and 0.325 at −20 °C, 0.325 at −40 °C; jobs 3774288, 3774289, 3774293).* Local: `LOC/keff_sintering_campaign/` (per-run `plots/keff/`, overlays in `compare/batch3a/`, snapshot figure in the 1701 −20 run's `plots/keff/snapshots/`).
    - All reached 30 d. Health is clean, and the k_eff schedule matches the solver exactly on every run.
    - **2·τ_sub delivered:** 232 steps at −20 °C (368 at 1.09) in 1 h 42 min for $2.48; −40 °C took 52 min for $1.27.
    - **T collapse holds at 2·τ_sub:** k at matched SSA, −20 vs −40, agrees to 0.06%. The time ratio is 7.86 against the τ_sub ratio of 7.72, the same 1.8% offset the old 1.09 runs show (7.84). So the dtmax change moved the kinetics by < 0.3%. The direct same-packing dtmax check (gated 301) is still to come.
    - Rise from 11 τ_sub to 30 d: +24.6% (0.275) and +25.1% (0.325), against +25.9% / +25.0% on the gated packings.
    - Absolute k_iso at 30 d for 1701 is 5% below gated 301. 1701 has the lowest contact count of the new 0.325 set (z_band 3.22 vs 3.37), so this is one seed, not the recipe.
    - Anisotropy varies a lot by seed: k_xx/k_yy is 1.20 on 1601 against 0.86 on the gated 601, and 1601 has k_xy ≈ 12% of k_xx. Report anisotropy only as a seed mean.
    - **Found:** every-5-step sampling before 11 τ_sub left 3 samples across the first-day rise (corners in the k–time plots), so `-keff_freq 1` from 3b on. Also, the production manifest's run list was never written (a `grep` under `pipefail`); fixed.
    - At −40 °C only 19 steps fall after 11 τ_sub, all sampled. That's enough, because its k(SSA) path is the −20 °C path.

## Planned

- [x] **Convergence study** `batch_rve.txt` — done 2026-10-04. HPC `CAMP/` (L/R 20–80 runs plus gated 302–305). Local `LOC/keff_sintering_campaign/`. Results: `REPO/studies/keff_sintering/rve_convergence/README.md`.
  - **Finding:** the sintering rise is size-independent: 27.2 ± 1.0% at L/R 40 against 27.9 ± 1.3% at L/R 80. Absolute k_iso at L/R 40 is within the standard error of L/R 80 (+0.6 to +1.3% ± 4–5%). L/R 30 sits 13% high in absolute k, read as realization statistics; its rise is normal. The gated vs ungated packings differ by +1.6 ± 1.4 points in rise, which is not significant.
- [ ] (original plan) **Convergence study** `batch_rve.txt`: 26 runs, about $155.
  - **What:** φ 0.325, −20 °C, ungated L/R 20/30/56/80, plus the 5 gated L/R 40 packings. The L/R 40 point is production seeds 1701–1705, from 3a and 3b.
  - **Data:** HPC `CAMP/`. Analysis: `rve_convergence/analyze_rve.py`.
  - **Question:** is L/R 40 within the standard error of L/R 80? And how much did the old gates shift k_eff?
  - *Result:* —
- [x] **3b**, the rest of −20 °C — done 2026-10-04. HPC `CAMP/`, local `LOC/keff_sintering_campaign/` (per-run `plots/keff/`; overlays `compare/batch3ab/`; table `compare/summary/stage_summary.csv`; health `compare/health_check.txt`).
  - **Findings (−20 °C, 5 seeds per φ):**

    | φ | rise from day 1 to 30 d | k_xx/k_yy |
    |---|---|---|
    | 0.275 | 27.4 ± 2.3% | 1.02 |
    | 0.325 | 27.2 ± 2.3% | 1.02 |
    | 0.375 | 24.5 ± 1.6% | 0.91 |
    | 0.425 | 21.2 ± 2.5% | 0.79 |
    | 0.475 | 19.0 ± 3.2% | 0.64 |

  - The rise falls with porosity above 0.325, by more than the seed scatter.
  - **The anisotropy reverses** (k_yy > k_xx) above φ ≈ 0.35, against the contact fabric, and its scatter grows toward percolation. It is not explained by contact orientation or by ice chord length (`REPO/studies/rve_anisotropy/chord_anisotropy.py`, README). The amplification-near-percolation idea is open.
  - Health: 1–2 phase-bound rollbacks (retried at half dt) in 9 runs. One run had 5.3% CFL rollbacks.
  - **Quirk:** when the final step is rejected, the run ends up to 0.25 d early and the 30-day snapshot is lost. Data are unaffected.
  - Full snapshots checked for φ 0.475 seed2001, φ 0.425 seed1902 and L/R 80 seed1401: snapshot figures look right; vertical ice columns at φ 0.475.
- [ ] (original plan) **3b**, the rest of −20 °C (`batch3b_T-20.txt`): 20 runs, about $80.
  - **Question:** is the seed scatter smaller than the φ trend, and does k_xx/k_yy follow the contact fabric?
  - *Result:* —
- [ ] **Matrix reordered (2026-10-03):** 3c–3e run seeds 1–3 of each φ; seeds 4–5 at the non-(−20 °C) temperatures moved to an optional final stage 3f. Temperature is a per-packing time rescaling, so the −20 °C column (5 seeds) carries the seed means.
- [ ] **3c**, −30 and −40 °C (`batch3c_cold.txt`): 29 runs (seeds 1–3).
  - **Question:** does every φ collapse onto its −20 °C curve?
  - *Result:* —
- [ ] **3d**, −10 °C (`batch3d_T-10.txt`): 15 runs (seeds 1–3).
  - *Result:* —
- [ ] **3e**, −5 °C (`batch3e_T-5.txt`): 14 runs (seeds 1–3), the most expensive temperature.
  - *Result:* —
- [ ] **3f (optional)**, seeds 4–5 at −5/−10/−30/−40 °C (`batch3f_seeds45.txt`): 40 runs. Run only if a porosity fails the collapse test or a reviewer asks for five packings in every cell.
  - *Result:* —
- [ ] **Pre-manuscript checks** (no new simulations):
  - rev64 seeds 2–4 analysis (already local);
  - φ 0.475 connectivity at the diffuse band;
  - the day-7.5 kink (batch 2, seed 1);
  - our anisotropy direction against measured snow;
  - the seam-relaxation comment;
  - the `packing_design` README correction;
  - the Nicoli 2011 citation;
  - the t = 0 normalization caption note.
- [ ] **Sensitivity: eps ×2 only** (`batch_eps2.txt`): φ 0.325, seeds 1701–1703, −20 °C, eps = 2 µm with the production experiment (dtmax unchanged in seconds). Asks whether the k_eff trajectory depends on the interface width. Decided 2026-10-04: α_c, R_ave and σ_ln arms dropped as giving the manuscript nothing. R_ave is pure time rescaling (grain-size study), and the α_c regime limit is stated, not run.
  - *Result:* —
- [ ] **Analysis** (stage 7):
  - k_eff vs SSA at fixed φ;
  - anisotropy from seed means;
  - per-snapshot metrics (Euler characteristic, chord lengths, D_eff);
  - ranking which factors matter most.
- [ ] **Optional:**
  - the k_eff solver benchmark (`solver_benchmark/`, replay only, low value now);
  - a portable `-march` build if a job hangs again.

## Not run / withdrawn

- **Batch 1** (`batch1_T.txt`, cold end on the pilot packings): no results recorded anywhere, so treated as never run. It was superseded by the batch 3 matrix on the gated and then option-B packings.
- **R_feat/R_ave = 1/50 and the eps-correction runs** (PLAN/ROADMAP): withdrawn, because they were sized against the arithmetic-law bias.
