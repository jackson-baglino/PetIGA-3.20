# Campaign: how sintering changes k_eff

> **Working to-do list: [`TODO.md`](TODO.md).** Submit every manuscript run
> with `scripts/HPC/submit_keff_production.sh <stage file>` -- it owns all
> options.

The executable plan, with a **stop rule at every stage**. `PLAN.md` holds the
study design and `ROADMAP.md` the earlier staging; this file supersedes both
wherever the tensor conductivity law changed the answer, and says so each time.

**The deliverable is comparative** — how much sintering changes `k_eff`, and
which levers matter most — not absolute `k_eff` in W/m/K.

**Scope is steady-state.** The cell problem, the surface excesses and the
tensor law concern the corrector problem on a frozen `φ`. Nothing here may be
used to justify an interpolation change in the transient metamorphism model
(`ρc`, latent-heat placement). If something there looks analogous, it is a
separate question.

## How to read a checkpoint

Each stage names **one number**, the threshold it must clear, and what happens
if it does not. Costs are at the tier-1 Resnick rate, $0.012/core-hour
(`scripts/HPC/hpc_cost.sh`); cost ≈ ranks × wall-hours × rate.

> **Cost model corrected 2026-09-25 — earlier estimates were ~15–20× high.**
> They took the pilot's ~6 650 steps at −20 °C as the step count and scaled it
> by ρ_vs. That count is an artefact: pilot leg 1 ran at the old flat
> `-dtmax 200` (**97% of its steps at that cap**, not "2.4%"), and after the
> derived `dtmax = 1.09·τ_sub` (a7dda40) arrived for leg 2, days 15–30 took
> **165 steps**. The rev64 run, on the derived `dtmax` throughout, took
> **371 steps in 2 h 19 min at 615 ranks, ≈ $17** for 30 days (79% of steps at
> the cap). At 241 ranks: ~22.5 s/step, ~12 s per tensor `k_eff` sample, and
> steps ≈ `t_final/dtmax` + ~70 early CFL-limited steps. Every 30-day run
> fits one 24 h job — no restart legs at any temperature in −40…−10 °C.

| stage | what | cost | running | figure |
|---|---|---|---|---|
| 0 | tensor law + both ladders | done | — | ✔ committed |
| 1 | pilot replay under tensor (no new simulation) | **$43 actual** | $43 | ✔ **done** |
| 2 | batch 1 — cold end of the T axis | < ~$30 | ~$73 | ✔ |
| 3 | batch 2 — warm end, −20/−15/−10 × seeds 1,3,4 | **$342 actual** | ~$415 | ✔ **done** |
| 4 | production packings | local | ~$223 | — |
| 5 | main matrix | measured at stage 2 | — | ✔ **SLIDES** |
| 6 | sensitivity arms | ~12 runs | — | — |
| 7 | analysis | — | — | ✔ **SLIDES** |

---

## Stage 0 — the conductivity law · DONE

`-keff_interp` defaults to `tensor` since 2026-09-23: arithmetic along the
interface, harmonic across, `n = ∇φ/|∇φ|`. Both branches satisfy `K(0)=K_air`
and `K(1)=K_ice`, so the law changes nothing outside the band and cancels the
band's first-order normal-flux bias.

| ladder | arith | tensor |
|---|---|---|
| planar slab, `studies/keff_sharp_limit/verification/` | `k_⊥` +59.4% → +3.8% | flat at the sharp value to **4e−7** |
| single ice cylinder (Rayleigh), `studies/keff_sharp_limit/disk/` | +18.19% → +1.93%, matching the no-free-constant dipole prediction to −1.9% | +1.212% → **+0.006%**, order **≈ 2.5** |

Writeup: `effective_thermal_cond/docs/tensor_conductivity_law.tex`.

**Every run before 2026-09-23 used `arith`**, the 2026-09-16 pilot included.

---

## Stage 1 — replay the pilot · GATES EVERYTHING

**Nothing is simulated here.** The four pilot trajectories finished in
September. A replay re-opens their stored snapshots, reads `φ`, and re-solves
the steady-state corrector problem under the tensor law. Full detail, including
why four runs become eight jobs, is in
[`coefficient_fix/README.md`](coefficient_fix/README.md).

Neither verification ladder has a **contact**, and the contact is the point:
the packing's grains are seated in exact tangency, so the path between two
grains starts as a single point and the arithmetic band bridges it with
conducting material. The single-cylinder benchmark's ~+8% at the pilot's
`eps/R_ave = 0.02` is the isolated-interface part of that error with the
contact effect excluded by construction — a lower bound, not an estimate.

```bash
./scripts/HPC/submit_keff_replay.sh --dry-run \
    --roots $SCRATCH/enceladus_DSM/batch_2026-09-16__13.16.11_pilot_keff \
            $SCRATCH/enceladus_DSM/packing_2D_pilot_phi0.325_Rave50um_LR40_seed*_L2mm_eps1000nm_perxy_T-20
venv_enceladus/bin/python studies/keff_sintering/coefficient_fix/compare_laws.py <batch>
```

Give both roots — all four seeds were resumed, and a resume leg lands in the
single-run tree, not the batch parent. 8 jobs, ~$19, tensor only: each run's
existing in-line `k_eff.csv` **is** the arith result, so there is nothing to
buy on that side.

### RESULT (2026-09-24) — the gate passes, continue to stage 2

8 jobs, stride 12, ~$5.4 each. Measured, from the 1-day baseline:

| seed | arith | tensor |
|---|---|---|
| 1 | +13.9% | **+27.4%** |
| 3 | +18.3% | **+31.9%** |
| 4 | +15.0% | **+27.2%** |
| ensemble | **+15.8%** (sd 2.3) | **+28.8%** (sd 2.6) |

**The headline is ×1.83, from +15.8% to +28.8%.** Seed 2 is excluded: its
leg-1 job lost 30 MPI ranks to `Slurmd could not connect IO` on `hpc-19-16`,
so only its second half exists and it has no 1-day baseline.

The mechanism is confirmed directly. On seed 1 `k_eff(0)` falls from **0.6312
(arith) to 0.3695 (tensor)** — the arithmetic law was reading **+71% high** on
tangent grains, exactly the pre-welded-contact effect. By day ~10 the two
absolute curves cross: once necks are real, the band correction matters less.

Independently corroborated by `studies/keff_sharp_limit/disk/` f-sweep: the
arith error grows with solid fraction (+1.9% at `f`=0.05 to **+38.6%** at
`f`=0.50, at the campaign's `eps/R`=0.02) while the tensor error stays **under
1%** across the whole range. Two unrelated geometries, one mechanism.

Costs, against what was predicted: **~100 s/sample was right** (measured
83–110 s). `ksp_its` went 106–130 (arith) → 136–180 (tensor), **~1.4×**, not
the ~2.5× predicted from the cylinder benchmark — so the tensor operator is
cheaper on the packing than feared.

**Consequence for stage 2:** the cold-end threshold is now against the tensor
seed scatter, **sd 2.6 points**, not the arith 2.1.

### The baseline is t = 1 day

The initial condition is an analytic clamped sum of `tanh` profiles, not an
equilibrated phase field, so the first hours of every run are the field
relaxing rather than sintering — and that is where the error is largest.
Measured on pilot seed 1, `d ln SSA / d ln t` runs 0.000 → −0.075 over the
first ~8 h and is flat at ~−0.077 after. **Nothing is read at `t = 0`.**

It matters: the pilot ensemble rise is +19.6% from `t = 0` and **+15.9% from
1 day**. Every earlier `+19.6%` in this campaign came from `t = 0` and is
superseded. `ROADMAP.md:260`'s ban on `k_eff(t=0)` therefore stands — for this
reason, which is about the initial condition and applies under any
conductivity law, not only the one it was written for.

> **Withdrawn here:** `PLAN.md:71-79`'s `R_feat/R_ave = 1/50` recommendation.
> Its case was a ~33%-vs-~15% `k_eff` eps bias, which was an *arithmetic-law*
> property and is gone. See "Necks are not resolved" below.
>
> **Also withdrawn:** `ROADMAP.md` Phase 2 / `PLAN.md` Stage 2's eps-correction
> runs were sized to measure the `O(eps)` `k_eff` bias the tensor law nulls.
> What survives is whether the *evolution* differs at smaller eps — which a
> replay cannot answer and which needs one real run, not a ladder. Deferred to
> stage 6.

---

## Necks are not resolved, by design

**The only length-scale ratio this campaign imposes is the RVE one,
`L/R_ave = 40`**, which is measured: seed scatter in `k_eff` falls from 12–13%
at `L/R_ave` 10–20 to 4.3% at 40. It must hold at **`t_final`**, not `t = 0` —
coarsening lowers it during the run.

`R_feat/R_ave` is **not** a resolution requirement. It is the cost knob: it
sets `eps = safety·R_feat`, hence `Nx`, hence the bill, and nothing else
depends on it. The campaign uses **1/25** (`eps` = 1 µm, `Nx` = 2829, 24 M DOF,
241 ranks) because that is what is affordable across ~80 runs.

At that mesh a neck is only represented above `r/R = √(12·eps/R) = 0.49`, so
**most of every trajectory's neck is below the floor.** That is accepted, not
fixed: resolving necks at this domain size is not something these simulations
can be asked for, and buying it would cost 4× the campaign for a quantity the
manuscript does not report. `k_eff`, SSA and porosity are the reported
quantities and none of them needs a resolved neck.

**Manuscript obligation.** State the floor as a measured number and its
consequence — necks are under-resolved, so no neck-growth exponent is claimed
from these runs, and the sintering *rate* inherits an error of unquantified
sign at early times. Do not present it as a caveat on `k_eff` itself; the
conductivity is a property of the field as it stands, whatever produced it.

Do not reopen this. If a later question genuinely needs neck radius, it needs
a different, smaller-domain experiment — not a finer mesh on an 80-run matrix.

---

## Stage 2 — batch 1, the cold end · < ~$30 (re-priced 2026-09-25)

The axis the 2025-09 sweep destroyed, rebuilt with `alpha_c = 1e-3` constant
and `beta_sub0` therefore spanning 8.4× across the axis, one `eps` and one
`-eps_valid_temp` per temperature.

```bash
./scripts/HPC/submit_batch.sh --tag keff_T_batch1 \
    --tests-file studies/keff_sintering/batch1_T.txt \
    --extra-opts "-keff 1 -keff_step0 1 -keff_t_interv 5.0e4 \
                  -keff_ksp_type cg -keff_pc_type gamg"
```

−40 °C and −30 °C, seeds 1 and 2 — **two temperatures × two seeds, not three
× one**. The pilot already fixes the scatter at −20 °C; what is missing is
scatter at the cold end, which is what makes a temperature difference a result
rather than a realization.

~~The table that stood here (~825 / ~2 450 steps, ~$16 / ~$47 each) scaled the
pilot's dtmax-200 step count and is withdrawn — see "Cost model corrected"
above.~~ Colder means a larger derived `dtmax`, hence fewer steps than the
−20 °C run's ~370: each cold-end run costs less than a −20 run, ≲ $7.
**Record the measured s/step here.**

**Read:** the `k_eff` rise at −40 and −30 against the −20 ensemble, all on the
tensor scale.

- **cold-end rise separates from −20 by more than the seed scatter** →
  temperature is a real lever. Continue to stage 3.
- **it does not** → 30 days is too short at −40 °C for the effect to develop
  (ρ_vs is 8× lower, so the same wall-clock buys 8× less sintering). Do **not**
  buy the matrix. Either extend `t_final` at the cold end, or reformulate the
  axis at matched sintering state, and re-read this checkpoint.

---

## Stage 3 — batch 2, the warm end · **$342 actual** · **SLIDES 1**

φ = 0.325, T ∈ {−20, −15, −10} °C, seeds 1, 3, 4, 30 d. Batch file and the
full rationale: [`batch2_T_warm.txt`](batch2_T_warm.txt). ~~"−10 °C is ~16 700
steps, ~111 h and five 24 h restart legs"~~ — withdrawn, same error as stage 2:
at `dtmax` 3449 s it is ~850 steps, ~8 h, one job.

```bash
./scripts/HPC/submit_batch.sh --tag keff_T_warm \
    --tests-file studies/keff_sintering/batch2_T_warm.txt \
    --out-root /resnick/groups/rubyfu/jbaglino/simulation_outputs \
    --extra-opts "-keff 1 -keff_step0 1 -keff_freq 1 -t_out_log 150 -t_out_log_t0 60 -keff_ksp_type cg -keff_pc_type gamg"
```

| T | dtmax | est. steps | est. wall @241 (k_eff every step) | est. cost |
|---|---|---|---|---|
| −20 | 8526 s | ~370 | ~3.5 h | ~$10 |
| −15 | 5376 s | ~560 | ~5.5 h | ~$16 |
| −10 | 3449 s | ~850 | ~8 h | ~$24 |

- **Output.** `k_eff` at every accepted step (`-keff_freq 1`), because the
  deliverable plot is `k_eff` vs SSA and `SSA_evo.dat` is per step too: ~300
  samples in the 1–30 d window at −20 °C. The old `-keff_t_interv 5e4` put only
  ~6 samples in 1–3 d, where SSA moves most. Snapshots: `-t_out_log 150` from
  60 s, plus step 1 (~29 GB/run). Without it `-outp 1` writes **every** step,
  since `t_final/dtmax` < 1000 never trips the 1000-snapshot cap.
- **Written to the group directory** (`--out-root`), not `$SCRATCH`: paper runs.
- **−20 °C is re-run, not reused.** The pilot mixes dtmax 200 (0–15 d) and 8526
  (15–30 d), and d ln SSA/d ln t jumps −0.087 → −0.096 exactly at the switch.
  Seeds 1/3/4 match the pilot's valid seeds, so 0–15 d is a free dt check.

**Read, first:** new −20 vs pilot, same seed, relative `k_eff` rise at 15 d.
If they differ by more than the seed sd (2.6 points), `dtmax = 1.09·τ_sub` is
too loose — lower `--dtmax-over-tau` and re-run before reading −15/−10.
**Then:** `k_eff` rise vs T with seed error bars, and `k_eff` vs SSA per T.

### RESULT (2026-09-27) — temperature is a pure time rescaling

All 9 runs reached 30 d. Local copy:
`~/SimulationResults/HPC_results/enceladus_DSM/GrainPackingSintering/keff_T_warm_phi0.325/`
(cluster: `.../simulation_outputs/enceladus_DSM/batch_2026-09-25__18.06.50_keff_T_warm`).
Overlays in its `compare/`, made by `postprocess/compare_keff.py`.

- **dt check passes.** New −20 (dtmax 8526 throughout) vs pilot leg 1 (dtmax
  200), same packings: k_iso agrees to 0.01% from day 5; rise 1→15 d differs
  by 0.13–0.38 points against the 2.6-point stop threshold.
- **At matched SSA, k_iso is the same at every T** — within 0.5% (mean
  0.01–0.10%) for each packing, over SSA 13 200–17 900 m⁻¹. Temperature moves a
  packing along one k(SSA) path; it does not change the path.
- **The speed along it is exactly τ_sub(−20)/τ_sub(T)**: the time to reach a
  given SSA is 1.59× shorter at −15 and 2.47–2.50× shorter at −10, every seed,
  against τ_sub ratios 1.59 and 2.47. With every run normalized at the same
  t/τ_sub (11), the k/k_b curves of all three temperatures coincide.
- **So "k_eff rise from 1 d" is flat in T** (+30.0 / +29.7 / +30.5%, sd 2.3–3.3,
  paired differences within ±2.1 points) **only because a fixed-time baseline
  is a different sintering state at each T.** Absolute k_iso at 30 d is 8–10%
  higher at −10 than at −20 for the same packing, because it got further.
- **This is a property of the model regime, not a finding about snow.** At
  α_c = 1e-3 the problem is attachment-limited (L* ≫ R, see
  studies/packing_design §4), and temperature then enters only through the
  kinetic prefactor. The collapse would break where the regime changes:
  vapour-diffusion-limited kinetics (larger α_c or grains), or an imposed
  temperature gradient.
- **Anisotropy** k_xx/k_yy at 30 d: 1.092 / 1.089 / 1.079 (sd 0.08–0.11) —
  set by the packing, not by T. k_xy ranges from −7% to +10% of k_xx by seed,
  with random sign.

**Consequence for stage 5.** Along T, the matrix buys nothing but τ_sub(T):
any temperature's trajectory is the −20 trajectory with time rescaled. Unless
a regime change is added on purpose, one T per porosity is enough, and the T
axis is reported as a rate law.

**Cost, measured: $342, not ~$150.** Step counts and s/step were as predicted
(368/545/814 steps, ~20 s/step at 241 ranks). The overrun is **k_eff at every
step**: 50–84 s per sample (median), 67–78% of wall time. The rev64 run solved
the same corrector in 4 s at 615 ranks with the same iteration counts; the time
per CG iteration swings 10× (0.06–0.55 s) within single runs here, so it looks
environmental (node contention), not algorithmic — unexplained. Seed 1 at −10
used 23 h 56 m of its 24 h limit. **For future runs: sample k_eff every ~5–10
steps, or by SSA change, not every step.**

> **SLIDES 1 — "the temperature lever is real and we can resolve it"**
> `k_eff` rise vs T with seed error bars; the audit table showing why the old
> sweep could not have seen it; the conductivity-law correction from stage 1.
> The ask: sign-off on the condition matrix before stage 5.

---

## Stage 4+5 (revised 2026-09-28) — packings built; batch 3 is the matrix

Superseded plan below kept for the record. What was actually done:

- **Packings.** `inputs/packings/keff_LR40/`: φ {0.275, 0.325, 0.375, 0.425,
  0.475} × 5, one unique seed each, y-seam gate on, and
  `make_packings.sh` as the record.
  - φ = 0.475 has the percolation gate off: the ice spans both axes in 1 of 5.
    It is kept to show where 2D stops working, not as a snow condition.
  - The RVE and anisotropy basis is in `studies/rve_anisotropy/README.md`.
- **Batch 3.** [`batch3_phi_T.txt`](batch3_phi_T.txt): the 25 packings ×
  T {−5, −10, −20, −30, −40} °C = 125 runs.
  - k_eff every 5 steps (every step at −40).
  - 50 log-spaced snapshots per run.
  - **Submitted in five checked stages**, never whole: 3a shakedown (7 runs),
    3b −20 °C rest (20), 3c −30/−40 (49), 3d −10 (25), 3e −5 (24). Each
    `batch3?_*.txt` lists the check that must pass before the next.
  - ~$1.7k on the batch-2 cost basis; up to ~$2.8k at 60k DoF/core if steps
    don't speed up.
- **The weekend runs (batch 2) are not rerun.**
  - Their conclusions (the dt check, T as a time rescaling) are paired
    comparisons on the same packings, and hold.
  - They are not pooled with batch 3: those packings predate the seam gate.

## Stage 4 — production packings · local, minutes

φ ∈ {0.25, 0.30, 0.35, 0.40} × 5 seeds, `L/R_ave` ≥ 40 at **`t_final`**, not
`t = 0` — coarsening lowers it during the run. Use `--max-void-per-L`, never
the mean-radii gate, when domain size varies. Porosity ceiling is **0.40**: the
solid stops percolating between 0.40 and 0.45, above which `k_eff` is not the
conductivity of a connected medium. Gate: solid percolates both axes on every
packing.

---

## Stage 5 — the main matrix · **SLIDES 2**

φ {0.25, 0.30, 0.35, 0.40} × T {−40, −30, −20, −10} × 5 seeds, minus what
stages 2–3 already cover. **Price it from stage 2's measured s/step**, not from
this file's estimates.

**Read:** between-condition difference against within-condition SEM. If the
matrix is unaffordable, cut temperatures or porosities — **never seeds**. The
seed count is what makes a trend a result rather than an anecdote, and the SEM
is already the limit on what can be claimed.

---

## Stage 6 — sensitivity arms · ~12 runs

One at a time from the stage 5 centre point, 3 seeds each:

- **`alpha_c`** {1e-4, 1e-3, 1e-2} — the one genuine free parameter, and it
  moves `L*` across the grain scale, so it may change the *mechanism*.
- **`R_ave`** {×0.5, ×1, ×2} **at fixed `L/R_ave`** — the packings are then
  geometrically identical and only the physics changes.
- **`σ_ln`** {0.2, 0.5} — how much the suppressed-ripening 2D artifact matters.
- **`eps` ×2 on one condition** — what survives of the withdrawn stage-1
  eps work: whether the *evolution*, not `k_eff`, differs at smaller eps.

---

## Stage 7 — analysis · **SLIDES 3**

1. **Make the SSA claim falsifiable.** Regress `k_eff` on porosity alone, then
   porosity + SSA, and report the **partial** correlation — pooled *and within
   each condition*. A relation that appears only when porosity varies is weaker
   evidence than one that holds at fixed porosity.
2. **Anisotropy from the tensor eigenvalues**, not `k_yy/k_xx`: `|k01|/k00`
   reaches 0.081, so the principal axes are not the grid axes. It is flat in
   time and seed-dependent — inherited from the packing, not made by sintering.
3. **Per-snapshot metrics**: Euler characteristic, chord-length distributions
   (what Calonne/Torquato models are built on, so the sharpest test of 1),
   `D_eff`. **Not neck radius** — it is below its resolution floor for most of
   every trajectory at this mesh; see "Necks are not resolved".
4. **Lever ranking** with seed-scatter error bars — the actual deliverable.

The strongest result already in hand: `k_eff` rises while the ice fraction does
not move at all (`phi_bar` constant to seven figures in every pilot run). Mass
is conserved in a closed box, so density is constant *by construction* and a
density-only parameterisation predicts a flat line. It is not that density
fails to explain the rise — density **cannot** be the explanation.

---

## Standing rules

- Runs are launched by Jackson, never by the assistant.
- 2+ HPC jobs go through `submit_batch.sh`; chained single submits race in
  `obj/` and give "Stale file handle".
- Push before submitting, so the cluster checkout matches.
- `merge` → `thin` → `health_check` on every batch, in that order.
- Production runs are 2D. The tensor law is verified in 2D; 3D is untested.
