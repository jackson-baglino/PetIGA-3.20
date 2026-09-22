# Campaign roadmap: sintering → k_eff

From the jobs running now to the finalised simulation set. Ordered by
dependency, not by preference: each gate's outcome changes what the next step
should be, and the three marked **SLIDES** are the points where that choice is
worth putting in front of an advisor before spending the compute.

**The deliverable is relative/comparative** — how much sintering changes k_eff,
and which levers matter most — not absolute k_eff in W/m/K. That choice is
already supported: the L/R=64 run reproduced the relative response to −0.3 sd
while sitting +5% in absolute.

## Triage (2026-09-22)

The checking was hitting diminishing returns, so the remaining items were
sorted by what skipping each actually costs. Two are not worth running.

| item | verdict | why |
|---|---|---|
| L/R=64 seed-1 rebuild | **skip** | The decision is already made — relative rise agrees across sizes (−0.7 sd), and we claim ratios. Costs one caveat sentence on a slide. |
| curvature at ε/4 | **run in parallel** | 2 tiny runs (Nx 2524). Gates nothing — a uniform 2D→3D rate offset cancels in every cross-condition comparison. Only needed for a timing claim and the discussion. |
| sharp/tensorial k_eff coefficient | **do first** | Not a test, a fix: the headline reads +17.8% arithmetic vs +57.0% sharp. A few lines in `keff_cell.c` plus `-keff_replay`, no new runs. Skipping means 80 runs measure something 3× too small. |
| resolution `R_feat/R_ave` | **decide explicitly** | The one genuinely open question, and the expensive one. See below. |

### The resolution decision is the real one

At the pilot's `R_feat/R_ave = 1/25`, a neck is only represented above
`r/R = √(12ε/R) = 0.49` — half the grain radius. Much of every trajectory is
below that floor, and the campaign's headline is the *rise from t = 0*.

| `R_feat/R_ave` | ε | r/R floor | Nx | DOF | cores |
|---|---|---|---|---|---|
| 1/25 (pilot) | 0.96 µm | **0.49** | 2829 | 24M | 301 |
| 1/50 | 0.48 µm | 0.35 | 5657 | 96M | 1201 |
| 1/100 | 0.24 µm | 0.25 | 11314 | 384M | 4801 |

Going to 1/50 is **4× the whole campaign**. Guessing costs either 4× compute or
80 redone runs.

**De-risk it inside the campaign instead of ahead of it.** Run the first
condition at both 1/25 and 1/50, 3 seeds each (~4.5k core-allocations, ~19% of
an 80-run campaign at 1/25), and compare the *relative* rise:

- agrees → run everything at 1/25, having saved 4× on the other 79 runs;
- differs → found on 6 runs instead of 80.

| | what | cost | gate |
|---|---|---|---|
| **0** | REV ensemble + curvature calibration | done | — |
| **1** | sharp/tensorial coefficient + replay | code change, no runs | **required** |
| **2** | first condition at 1/25 **and** 1/50, 3 seeds each | ~4.5k core-alloc | **decides resolution** |
| **3** | production packings | minutes, local | — |
| **4** | main campaign | see Phase 5 | **SLIDES 2** |
| **5** | sensitivity arms | ~12–16 runs | — |
| **6** | analysis and writeup | — | **SLIDES 3** |
| *par* | curvature at ε/4 | 2 tiny runs | alongside, gates nothing |

---

## Phase 0 — in flight

```
./scripts/HPC/submit_rev_check.sh         # 4 seeds, L/R_ave = 64, ~769 cores each
./scripts/HPC/submit_curvature_calib.sh   # axisym vs planar two-grain pair
```

Download tables only (`k_eff.csv`, `SSA_evo.dat`, `outp.txt`, `*.opts`) — the
REV question needs no snapshots, and at L/R=64 each is ~2.6× larger.

---

## Phase 1 — GATE: domain size, and what 2D can claim

### 1a. REV

```bash
venv_enceladus/bin/python studies/keff_sintering/compare_rev64.py <pilot> <rev64>
```

Compare **ensemble to ensemble**, not to the single earlier seed.

| outcome | meaning | consequence |
|---|---|---|
| L/R=64 mean within ~2 SEM of L/R=40 | 40 is converged | **use L/R=40** — the campaign costs 2.56× less |
| outside, relative response still agrees | absolute biased, ratios fine | **use L/R=40** and report only ratios; state the offset |
| outside, relative response differs too | 40 is not usable | **use L/R=64** and cut the matrix to fit |

The middle row is the likely one and the cheapest. It is also the row that
makes the relative-only framing load-bearing rather than a convenience.

### 1b. Curvature

```bash
venv_enceladus/bin/python postprocess/neck_width.py <curvcal_axisym_run>
venv_enceladus/bin/python postprocess/neck_width.py <curvcal_planar_run>
```

Ratio of dr/dt is the 2D→3D rate correction. Planar should be **faster** —
a 3D neck is a saddle, `H ≈ 1/r − 1/ρ`, and the `+1/r` term that partially
cancels is absent in 2D.

> ### SLIDES 1 — "can this question be answered in 2D?"
> The deck that buys in (or does not buy in) to the whole approach. Have:
> - The three pore diagnostics (`pore_scale.png`): L/ξ_pore 59 → 52.5, largest
>   pore 3.3% → 4.0% of L, **pore never percolates** — homogenization holds,
>   tested three ways, because each is blind to a different failure.
> - The 2D limits, stated as measured numbers rather than caveats: pore
>   fragmented into ~230 pieces from t=0 (and **not worsening** — 243 → 224,
>   usable pore 78.8% → 85.0%); long-range vapour transport suppressed ≥100×;
>   curvature correction from 1b.
> - `rev64_compare.png`: absolute offset vs identical relative response — the
>   evidence for the relative-only framing.
> - The ask: agreement that the claim is comparative, and the domain size.

---

## Phase 2 — GATE: does the headline number need an ε correction?

**Do this before the campaign.** It is not polish; it decides whether the
headline relative number is right.

The ε bias on k_eff is **not constant in time**. It scales with `eps × SSA`,
and SSA falls 37% over a run, so the bias shrinks — and therefore does *not*
cancel in the ratio `k(t)/k(0)`. Propagating the measured ladder
(`studies/packing_design/bias.csv`) through the pilot:

| | measured | implied true | bias |
|---|---|---|---|
| t = 0 | 0.6331 | 0.5136 | +23.3% |
| t_end | 0.7570 | 0.6443 | +17.5% |
| **relative rise** | **+19.6%** | **+25.4%** | — |

So the headline may understate by ~30% *of itself*. The prefactor was fitted on
a different packing, so treat the size as indicative — but the **direction is
structural**: SSA falls, the bias shrinks, the measured rise under-states.

**Run:** one condition, two seeds, at `eps` and `eps/2` (`R_feat/R_ave` 1/25 →
1/50, Nx 4526 → 9052 if at L/R=64, or 2829 → 5657 at 40). Compare *relative*
rises, not levels.

- Rises agree within seed scatter → quote the measured number, no correction.
- Rises differ → extrapolate in ε and quote the corrected number with the
  ladder as evidence.

**This also settles `R_feat`**, which nothing else has. The neck-resolution
floor is `r/R = √(12·eps/R)`: 0.49 at 1/25, 0.35 at 1/50. If neck radius is a
*reported* metric, 1/25 is not usable — most of the trajectory is below its
floor.

---

## Phase 3 — production packings

Once domain size and `eps` are fixed:

```bash
for phi in 0.25 0.30 0.35 0.40; do for s in 1 2 3 4 5; do
  venv_enceladus/bin/python preprocess/generate_packing_gravity.py \
    --Lx <L> --porosity $phi --mean-r 50e-6 --sigma-ln 0.5 --periodic xy \
    --seed $s --max-void-per-L 0.0335 --band-per-mean-r <9.2*eps/R_ave> \
    --out inputs/packings/prod/phi${phi}_seed${s}
done; done
venv_enceladus/bin/python preprocess/generate_study_opts.py \
  --packings-dir inputs/packings/prod --out-dir inputs \
  --R-feat <R_ave/N> --alpha-c 1e-3 --Rave 50e-6 --temps -40 -30 -20 -10
```

Use `--max-void-per-L`, never the mean-radii gate, whenever domain size varies.
Porosity ceiling is **0.40**: the solid stops percolating between 0.40 and 0.45,
and above that k_eff is not the conductivity of a connected medium at all.

---

## Phase 4 — main campaign, first condition

Run **one** condition (φ = 0.325, T = −20 °C) × 5 seeds to completion before
launching the rest. It confirms cost, `t_final`, and that the seed count clears
the between-condition differences you need to resolve.

> ### SLIDES 2 — "the signal is real and we can resolve it"
> - `pilot_keff.png`: k_eff rises ~19% while **ice fraction does not move at
>   all** (0.675798 at every sample, drift exactly 0.0). Mass is conserved in a
>   closed box, so density is constant *by construction* and a density-only
>   parameterisation predicts a flat line. This is the strongest result you
>   have: it is not that density fails to explain the rise, it is that density
>   **cannot** be the explanation.
> - Regression: ice fraction alone R² = 0.13; + SSA → 0.85.
> - Seed scatter and the SEM, so the resolvable difference is explicit.
> - The ask: sign-off on the condition matrix before spending Phase 5.

---

## Phase 5 — main campaign, remainder

**porosity {0.25, 0.30, 0.35, 0.40} × temperature {−40, −30, −20, −10 °C} ×
5 seeds = 80 runs**, minus the condition already done.

If that is unaffordable at the Phase 1 domain size, cut **temperatures or
porosities, never seeds** — the seed count is what makes a trend a result
rather than an anecdote, and the SEM is already the limit on what can be
claimed.

Submit via `scripts/HPC/submit_batch.sh` (2+ jobs go through the batch
submitter). After download: `merge_restart_legs.py` → `thin_snapshots.py` →
`health_check.py` on every batch, in that order.

---

## Phase 6 — sensitivity arms

One at a time from the Phase 4 centre point, 3 seeds each:

- **α_c** {1e-4, 1e-3, 1e-2} — the one genuine free parameter, and it moves
  `L*` across the grain scale, so it may change the *mechanism* rather than
  only the rate.
- **Grain size R_ave** {×0.5, ×1, ×2} **at fixed L/R_ave** — the packings are
  then geometrically identical and only the physics changes, which isolates
  the grain-size effect from a domain-size effect.
- **σ_ln** {0.2, 0.5} — a direct test of how much the suppressed-ripening
  artifact matters. Near-monodisperse barely ripens, so if results are
  insensitive to σ_ln the 2D artifact is not biting.

---

## Phase 7 — analysis

1. **Make the SSA claim falsifiable.** Regress k_eff on porosity alone, then
   porosity + SSA, and report the *partial* correlation — pooled **and within
   each condition**. A relation that only appears when porosity varies is
   weaker evidence than one that holds at fixed porosity.
2. **Anisotropy from the tensor eigenvalues**, not `k_yy/k_xx`: off-diagonal
   `|k01|/k00` reaches 0.081, so the principal axes are not the grid axes.
   It is flat in time and seed-dependent (0.86–1.00) — inherited from the
   packing, not made by sintering.
3. **Per-snapshot metrics** worth adding: Euler characteristic (connectivity
   once grains merge and coordination stops being definable), chord-length
   distributions (what Calonne/Torquato models are built on, so the sharpest
   test of whether SSA is cause or proxy), `D_eff` (the 2D transport limit as
   a function of time).
4. **Lever ranking** — the actual deliverable. Effect size per lever with
   seed-scatter error bars.

> ### SLIDES 3 — results
> - The lever ranking, with error bars, as the headline.
> - One worked trajectory for mechanism **plus** the ensemble for the claim.
> - Limitations as measured numbers, not hedges: the ε bias and its correction
>   (Phase 2), the curvature factor (Phase 1b), no long-range ripening, both
>   phases cannot percolate in 2D, anisotropy is deposition-induced.

---

## Standing rules

- Runs are launched by Jackson, never by the assistant.
- 2+ HPC jobs go through `submit_batch.sh`; chained single submits race in
  `obj/` and give "Stale file handle".
- Push before submitting, so the cluster checkout matches.
- `merge` → `thin` → `health_check` on every batch, in that order. Thinning
  before merging frees nothing, because the snapshots are hardlinks.
- Never use `k_eff(t=0)` as a baseline: at t=0 the diffuse band is standing in
  for necks that have not grown, so it is an artifact by construction.
