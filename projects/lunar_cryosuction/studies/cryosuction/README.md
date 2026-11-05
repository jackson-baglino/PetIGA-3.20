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

| Test | Geometry | Experiment | Checks |
|---|---|---|---|
| A1 | `openface_2D_box60um_disc_offcenter` | `openface_T-20_curvflow_7d` | open face: a cut disc shrinks as the full disc would |
| A2 | same | `…_closed` (`-phi_open_bc 0`) | the 90° face misses that law, so A1 is a real test |
| A3 | `openface_2D_box60um_disc_interior` | `openface_T-20_curvflow_7d` | D_eff from dA/dt = −2πD |
| A4 | `openface_2D_box60x30um_cap` | `openface_T-20_curvflow_7d` | contact line at the wall/open-face corner |
| B1–B3 | `channel_2D_H125um_eps0.86um_deadend` | `cryo_T-20_h0.990_iso_3d_xiv{1e-3,1e-4,1e-2}` | ξ_v invariance (K&P 2009, Sec. II D); IC transient |
| B4 | same | `…_xiv1e-3_xiT1e-3` | K&P choice ξ_T = ξ_v |
| C1 | `throat_2D_L850um_eps0.86um_narrowL_deadend` | `cryo_T-20_h0.900_iso_5d` | pilot: two plugged throats, isothermal |
| C2 | `regolithpore_2D_L814um_eps0.86um_throat_deadend` | `cryo_T-20_h0.900_iso_5d` | pilot: regolith pore, central bridge plus wall caps |
| C3 | two-throat | `cryo_T-20_h0.900_G+1000_warmdead_5d` | warm dead end, +1000 K/m |
| C4 | two-throat | `cryo_T-20_h0.900_insdead_5d` | insulated dead end (latent-heat effect against C1) |

**Notes on the inputs:**
- All A–C runs use α_c = 1e-3: β_sub0 = 7.9408e6 s/m, τ_sub = 5.8604e3 s, and
  dtmax = 2τ_sub = 1.17e4 s, set in the experiment file.
- The two pilot geometries come from the deterministic builders:
  - `build_geometry_two_throat.py --orientation narrow-left`
  - `build_geometry_regolith_pore.py --ice-placement throat_bridge --tag _throat`
  - Their parameters match the `inputs/scratch/` originals line for line.
  - The mesh files are in `inputs/geometry/meshes/`.

**Before submitting:** `scripts/preflight.sh stage0_tests.txt` resolves every
pair. It then starts the solver with `-t_final 0` and prints τ_sub and the
boundary table. No time step is taken.

**Submit (HPC):**

```bash
scripts/HPC/submit_batch.sh --tag cryo_stage0 \
    --tests-file studies/cryosuction/stage0_tests.txt \
    --out-root /resnick/groups/rubyfu/jbaglino/simulation_outputs
```

**Analyze** (after download):

```bash
venv_cryosuction/bin/python studies/cryosuction/scripts/analyze_stage0.py --batch <batch folder>
```

## Scripts

- `scripts/supersat_compare.py`: compares the imposed undersaturation h − 1
  with the Gibbs–Thomson terms of a grain and of a contact-angle meniscus.
  It also runs the K&P Eq. 48 and quasi-steady checks.
- `scripts/preflight.sh`: the set-up check described above.
- `scripts/analyze_stage0.py`: the pass/fail numbers for A and B, and the
  outflow and mass balance for C.
- `verification/verify_face_bc.sh` and `verification/verify_bflux.sh`: the
  gates for the two solver additions.
