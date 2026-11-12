# Run progress: what paper 2 / AGU needs, and what has run

One page per glance: every simulation the plan calls for, by type, with its
state. Detail for each run (folder, model version, full outcome) is in
`RUNS.md`; the day-to-day project state is in `STATUS.md`; the plan itself is
`SimulationPlan_v2.docx` (iCloud, `.../SimulationProposalsSetup/`).

**Last updated:** 2026-11-11

**Full run plan (setup, schematics, every matrix):** `RunPlan_2026-11-11.docx`
in the iCloud plan folder, built by `docs/run_plan/build_runplan.py`. It
proposes production gradients G = 1, 3, 10 K/m (30 as one above-threshold
point) instead of 100/300, pending the user and the S7 face-stability scan;
types 1, 2 (gradient runs), 6 and 7C/7D change accordingly. Type 5 is
unaffected.

States: **done** (ran and analyzed) · **ran** (ran, analysis open) ·
**running/queued** (on the HPC) · **ready** (inputs built and preflight OK, not
submitted) · **to build** (inputs not written yet).

## Summary

| Type | What it asks | Runs | State | AGU? (D4) |
|---|---|---|---|---|
| 0 | Shakedown: does the solver survive harvesting runs? | stage 0, 0b, 0c, v2check | **done** 2026-11-05 | background |
| I1 | Is the 0.86 um mesh accurate enough? (issue 1) | 1 fine-mesh run | **done** 2026-11-07: within 3 % | methods |
| 1 | One wall droplet (cap): drift and growth vs gradient, vapor path, size | 19 | **ready** (`type1_tests.txt`) | certain |
| 2 | Several droplets: ripening vs gradient, shielding vs spacing | 14 | **ready** (`type2_tests.txt`) | certain |
| 3 | Lens in a flat pore: does its speed follow 2 d0 cos(theta)/H in series with diffusion? | 10 | **done** 2026-11-07: law holds to 1 % | trend |
| 4 | Wedge: does a taper move ice to the narrow end? | 7 | **done** 2026-11-07: yes, every wedge | trend |
| 5 | Throat: where does a lens come to rest vs gradient? | 12 | **ready** (`type5_tests.txt`) | trend |
| 6 | Constriction (dog bone): how much does a neck slow a lens? | 7 | **ran** 2026-11-09/11: 5 finished, 2 hit the 48 h limit (below) | certain |
| 7 | Realistic rough pore, mixed ice, mining conditions | 7A + 7B, 7C-1/2, 7D-1/2 (6) | 7A **ready**; 7B-7D **to build** | certain |
| 8 | Lunar physics | - | paper 3, not this model (D6) | no |
| S2 | Straight pore vs the same pore with one throat: throat as an equivalent length | 2 | **to build** | if time |
| S5 | Two-throat pore, warm top, G 300 and 1000 K/m, h_open 0.99 and 0.995: recapture and harvest cutoff | 4 | **to build** (1000 K/m needs a finer eps) | if time |
| later | "Too long" pore where ice gets stuck (about 150 cores) | 1 | after the other runs are final (user, 2026-11-09) | no |

S1 (eps pairs and a time-step halving per new geometry), S3, S4 and S6 are in
the plan (section 9) but not scheduled.

## Per type

### Done

- **Type 3** (10 runs, `batch_2026-11-07__08.00.43_type3`, report
  `REPORT_type3.md/.pdf` there): the lens speed follows the series law within
  1 % with nothing fitted; the closed control (3C-1) does not move.
- **Type 4** (7 runs, `batch_2026-11-07__08.00.47_type4`, `REPORT_type4`):
  ice moves to the narrow end in every wedge; series law with theta -/+ omega
  to 3-5 %. Every theta-60 lens reached a vapor face, so later geometries got
  longer margins.
- **I1** (`batch_2026-11-07__08.00.52_issue1_eps_pair`, `REPORT_issue1`): fine
  mesh within 3 % (area 0.31 %, vanishing time 3.0 %); production stays at
  eps 0.8584 um.

### Type 6: ran 2026-11-09 to 11-11 (`batch_2026-11-09__15.09.21_type6`)

Geometry: 1100 x 150 um pore, neck 300 um long (x = 400-700 um) of width wc;
the lens starts at x = 855-1045 um (cold end); warm end at x = 0; G 100 or
300 K/m; local-saturation vapor at both ends. wc 150 um is the straight
control (no neck).

| Run | Finished? | Lead edge enters neck / leaves neck / all ice past neck / reaches x = 0 (days) | Peak ice | Resolution (Pi > 1, share of steps) |
|---|---|---|---|---|
| 6A wc 150, G 100 | yes, 60 d, 29 h | 6 / 14 / 24 / 26 | +12.5 % | 0.5 % |
| 6A wc 80, G 100 | yes, 60 d, 27 h | 8 / 14 / 24 / 28 | +11.2 % | 1.2 % |
| 6A wc 40, G 100 | yes, 60 d, 35 h | 13 / 18 / 28 / 40 | +10.7 % | 4.1 % |
| 6A wc 150, G 300 | yes, 30 d, 30 h | 2 / 5 / 9 / 9.5 | +15.3 % | 1.6 % |
| 6A wc 80, G 300 | yes, 30 d, 35 h | 2.5 / 5 / 9 / 10.5 | +13.7 % | 17 % |
| 6A wc 40, G 300 | **no**: day 9.6 of 30 at 46 h | 4 / 6.5 / - / - | +12.3 % | **65 %** |
| 6B wc 40, theta 90, G 300 | **no**: day 11.7 of 30 at 46 h | 5 / 7.5 / 11 / - | +13.7 % | **57 %** |

Times are from snapshots about 0.5-1 d apart (ice_bodies.py), so +/- 1 d.

### Ready to submit (pushed 2026-11-11)

Types 1, 2, 5 (45 runs) and 7A. From the project root on the HPC after
`git pull`:

```
TARGET_DOFS_PER_CORE=50000 ./scripts/HPC/submit_batch.sh --tag type1 --tests-file studies/cryosuction/type1_tests.txt --out-root /resnick/groups/rubyfu/jbaglino/simulation_outputs -- --time=2-00:00:00
```

and the same with `type2` and `type5`. Type 7A (2.41 M DOF, about 48 cores
at 50k/core) gets a pre-check first (rule: jobs above 30-40 cores).

### To build

- 7B (insulated dead end), 7C-1/2 (warm top, G 100 and 300), 7D-1/2 (warm dead
  end, G 100 and 300): experiment files on layout A.
- S2, S5: geometries and experiments.

## Decisions waiting on the user

- D1 theta_ref = 60 deg; D4 which types go to AGU (STATUS.md).
- Type 7 ice fraction: 24.6 % now, or 30 %.
- Droplet angle: 40 deg is a placeholder (user, 2026-11-09).
- Type 6: what to do with the two runs that did not finish (see the reply of
  2026-11-11 and RUNS.md).
