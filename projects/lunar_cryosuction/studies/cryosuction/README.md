# Cryosuction study (paper 2)

How ice in a regolith fracture migrates when one end is open to a dry
reservoir and the other is a dead pore. Knobs: geometry (pore throats),
open-end humidity, thermal BCs, and contact angle. Everything runs at −20 °C
as a lab analog; lunar temperatures enter later by scaling (decision of
2026-11-04).

## Setup

- **Open end, x = 0.** T and ρ_v are pinned there with `-Tfix_faces x0` and
  `-rhovfix_faces x0` (opt-in, FREEZE ledger 60bd33a5).
  - φ is open on that face (`-phi_open_bc`, default on): nothing is imposed
    on it, so the ice passes through.
- **Dead pore, x = Lx.** Zero vapor flux. It is regolith (`-wall_faces …,x1`)
  with the same contact angle as the walls.
  - T is either fixed (`-Tfix_faces x0,x1`) or insulated (`-Tfix_faces x0`).
- **Outflow.** `bflux.dat` (ledger 770a9bb5) is the physical vapor and heat
  flux through each pinned face. It uses the unscaled D_v and k, so it does
  not depend on ξ_v.

## Stage 0 (shakedown, `stage0_tests.txt`, 12 runs)

**First submission (2026-11-04) was cancelled.**
- It reused eps 0.86 µm with α_c = 1e-3.
- At its front speeds, `comp_eps.py` called for eps 17–170× smaller (K&P
  Eq. 45 kinetic bound).
- The local A runs from that submission found the open-face film (see
  `figures/stage0_A_area_law.png`).

**Requeue settings:**
- α_c = 1e-2, so `-beta_sub0` is 7.9408e5 s/m.
- eps and mesh come from `comp_eps.py`, evaluated at each run's expected front
  speed. The speed is the series-law maximum, `scripts/front_speed.py`.
  Every geometry file records the commands used.
- dtmax = 2τ_sub sits in each geometry file. τ_sub is the value the solver
  prints.
- Every run has the resolution checker (FREEZE ledger 8c96ab09). It prints a
  boxed warning if the measured front speed breaks the bound or ice touches an
  open face, and it writes `resolution.dat`.

| Test | Geometry (eps) | Experiment | Checks |
|---|---|---|---|
| A1 | `openface_2D_box60um_disc_offcenter` (0.64 µm) | `openface_T-20_curvflow_22h_ac1e-2` | open face: the cut disc against the full-disc law |
| A2 | same | `…_closed` | 90° control |
| A3 | `openface_2D_box60um_disc_interior` | `openface_T-20_curvflow_22h_ac1e-2` | D_eff |
| A4 | `openface_2D_box60x30um_cap` | same | wall/open-face corner |
| B1–B3 | `channel_2D_H125um_eps0.39um_deadend` (0.39 µm) | `cryo_T-20_h0.990_iso_3d_ac1e-2_xiv{1e-3,1e-4,1e-2}` | ξ_v invariance; IC transient |
| B4 | same | `…_xiv1e-3_xiT1e-3` | ξ_T = ξ_v |
| C1 | `throat_2D_L850um_eps0.86um_narrowL_deadend` (0.86 µm) | `cryo_T-20_h0.990_iso_5d_ac1e-2` | two plugged throats, isothermal |
| C2 | `regolithpore_2D_L814um_eps0.86um_throat_deadend` (0.86 µm) | `cryo_T-20_h0.998_iso_5d_ac1e-2` | regolith pore. h = 0.998 because h = 0.99 needs eps 0.20 µm (34 M DOF). |
| C3 | two-throat | `cryo_T-20_h0.990_G+100_warmdead_5d_ac1e-2` | warm dead end, +100 K/m |
| C4 | two-throat | `cryo_T-20_h0.990_insdead_5d_ac1e-2` | insulated dead end |

C3 uses +100 K/m, not +1000 K/m. At h = 0.99, +1000 K/m would make the cold
open end supersaturated (local h = 1.031), so the vapor would flow the wrong
way.

**Before submitting:** run `scripts/preflight.sh stage0_tests.txt`. It
prints τ_sub, dtmax/τ_sub and the start-up resolution check for each pair,
and fails on a resolution warning.

**Submit (HPC):**

```bash
scripts/HPC/submit_batch.sh --tag cryo_stage0b \
    --tests-file studies/cryosuction/stage0_tests.txt \
    --out-root /resnick/groups/rubyfu/jbaglino/simulation_outputs
```

**Analyze** (after download):

```bash
venv_cryosuction/bin/python studies/cryosuction/scripts/analyze_stage0.py --batch <batch folder> --fig-a <png>
```

## Run ledger

Every run and batch, queued or finished, with its question, what to look
for, and the outcome: `RUNS.md`.

## Scripts

- `scripts/supersat_compare.py`: compares the imposed undersaturation h − 1
  with the Gibbs–Thomson terms of a grain and of a contact-angle meniscus.
  It also runs the K&P Eq. 48 and quasi-steady checks.
- `scripts/front_speed.py`: expected front speed (series law) for `comp_eps.py --vn`.
- `scripts/preflight.sh`: the set-up check described above.
- `scripts/analyze_stage0.py`: the pass/fail numbers for A and B, and the
  outflow and mass balance for C.
- `verification/verify_face_bc.sh`, `verify_bflux.sh`, `verify_resolution.sh`:
  the gates for the three solver additions.
