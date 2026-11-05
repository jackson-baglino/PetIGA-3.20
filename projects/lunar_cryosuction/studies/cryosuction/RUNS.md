# Run ledger: cryosuction study

One row per run or batch: where it ran, which model version, what question it
answers, what to look for in the results, its status, and the outcome. Add the
row when the run is queued, and update Status and Outcome when it ends.

- **Model** is the freeze tag the run's binary belongs to (`FREEZE.md`,
  Versions table).
- **HPC batches** live under `/resnick/groups/rubyfu/jbaglino/simulation_outputs/`
  and, once downloaded, under `~/SimulationResults/HPC_results/lunar_cryosuction/`.
- **Local runs** live under `~/SimulationResults/lunar_cryosuction/scratch/<geometry>/`.
- **Paper?** says whether the result can go into paper 2 as it stands.

## Active and queued

| ID | Queued | Where / folder | Model | Geometry : experiment | Question | Look for | Status |
|---|---|---|---|---|---|---|---|
| 0c-B1 | 2026-11-05 | HPC, `--tag cryo_stage0c` (`stage0c_tests.txt`) | v2 | channel_2D_H125um_eps0.39um_deadend : cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-3 | Baseline harvesting run: an ice bridge in a rectangular channel draining to a reservoir at h = 0.99 over 3 days. | Completes; bflux outflow / ice lost ≈ 1; outflow J_x0(t) and front position against the series law v = abs(sigma)/(beta0 + K l); Pi ≤ 0.6 after start-up; no RESERVOIR-FACE warning. | submitted (binary at ae3ddbd2: no start-up window in its checker; analyze_stage0.py applies it) |
| 0c-C3 | 2026-11-05 | same batch | v2 | throat… : cryo_T-20_h0.990_G+100_warmdead_5d_ac1e-2 | Same diagnostic, with a warm dead end (v1 aborted at 4.3e4 s). | As 0c-C1. | submitted (binary at ae3ddbd2: no start-up window in its checker; analyze_stage0.py applies it) |
| 0c-C4 | 2026-11-05 | same batch | v2 | throat… : cryo_T-20_h0.990_insdead_5d_ac1e-2 | Same diagnostic, insulated dead end (v1 completed but overshot to phi 1.11). | As 0c-C1; does max phi_ice now stay ≤ 1.01? | submitted (binary at ae3ddbd2: no start-up window in its checker; analyze_stage0.py applies it) |
| 0c-C2 | 2026-11-05 | same batch | v2 | regolithpore_2D_L814um_eps0.86um_throat_deadend : cryo_T-20_h0.998_iso_5d_ac1e-2 | Regolith pore with wall caps and a central grain at h = 0.998: which ice goes first? Rerun of a clean v1 run so the result is v2. | Order in which the caps and the grain disappear; outflow / ice lost ≈ 1; agreement with the v1 run (0b-C2) where v1 had no rollbacks. | submitted (binary at ae3ddbd2: no start-up window in its checker; analyze_stage0.py applies it) |

## Finished

| ID | Ran | Folder | Model | Geometry : experiment | Question | Outcome | Paper? |
|---|---|---|---|---|---|---|---|
| v2check-HPC | 2026-11-05 | `~/SimulationResults/HPC_results/lunar_cryosuction/batch_2026-11-05__07.55.32_cryo_v2check/` | v2 | channel_2D_H125um_eps0.39um_deadend : cryo_T-20_h0.990_iso_20ks_ac1e-2 | Same as v2check-local, on the HPC | **Matches v2check-local**: 530 rows, same steps, ice within 5e-7 relative (round-off from a different rank layout); 0 rejections, 0 rollbacks. | no (2e4 s check) |
| 0c-C1 | 2026-11-05 | `~/SimulationResults/HPC_results/lunar_cryosuction/batch_2026-11-05__08.46.44_cryo_stage0c/throat_2D_L850um_eps0.86um_narrowL_deadend__cryo_T-20_h0.990_iso_5d_ac1e-2/` (job 6262619) | v2 (binary ae3ddbd2, no start-up window) | throat… : cryo_T-20_h0.990_iso_5d_ac1e-2 | Does v2 stop the distorted-mesh crash (v1 aborted at 3.0e4 s)? | **Yes.** Completed 5 days: 1283 steps, 1.00 Newton/step, 0 rejections, 0 CFL or band rollbacks, phi_ice in [-4.8e-19, 1.0000], h in [0.990, 1.000] throughout, mass balance 1.0000. "CHECK THIS RUN" comes from Pi alone: (1) t < 1e3 s start-up, Pi 25–30 (the binary predates the window); (2) 1e3–2e4 s, Pi median 1.03 (max 1.31): the narrow-plug front sits at the K&P Eq. 45 bound on this mesh; (3) 5.3–6.0e4 s: the narrow plug pinches off at its walls and the two leftover wall caps vanish (Pi up to 1e3, v_n 1.5e-6 m/s at the final collapse; a real curvature-driven event, no rollbacks); (4) after 1e5 s, Pi 0.17. Science: the wide plug's open-side face stays at x ≈ 510 µm until the narrow plug is gone, then retreats to ≈ 595 µm by day 5; its dead-end face stays at ≈ 695 µm; 48 % of the ice left at day 5. Snapshots: `~/SimulationResults/review_scratch/2026-11-05_cryo_stage0b/0c-C1_snapshots.png`. | no: narrow-plug phase is at the Eq. 45 bound (needs the resized mesh) |
| v2check-local | 2026-11-05 07:57 | local, 10 ranks, `~/SimulationResults/lunar_cryosuction/scratch/channel_2D_H125um_eps0.39um_deadend/2026-11-05__07.57.41_cryo_T-20_h0.990_iso_20ks_ac1e-2_v2check` | v2 | channel_2D_H125um_eps0.39um_deadend : cryo_T-20_h0.990_iso_20ks_ac1e-2 | Do the v2 safeguards alone carry stage-0b B1 through 4.5–8.9e3 s? | **Yes.** Reached 2.0e4 s in 529 steps, 1.00 Newton/step, 0 rejections, 0 CFL rollbacks, 0 band rollbacks; phi_ice stayed in [-4.8e-19, 1.0000] (v1 overshot from 4.7e3 s). Pi steady 0.56 after start-up (44 warnings, all before t = 95 s). Ice matches v1 to 0.01–0.08 % of the change, so v2 only removed the defect. Mass balance 0.998. Dead-end ice never reached the reservoir. | no (2e4 s check) |
| 0-* | 2026-11-04 | HPC `cryo_stage0` (cancelled) and local A runs, `~/SimulationResults/review_scratch/2026-11-04_cryo_stage0_openface` | v1 | 12 stage-0 runs at eps 0.86 µm, alpha_c 1e-3 | Shakedown | Cancelled: eps 17–170× too coarse for the drive (comp_eps). The local A runs first showed the open-face film. | no |
| 0b-A1 | 2026-11-04 | `batch_2026-11-04__19.34.44_cryo_stage0b/openface_2D_box60um_disc_offcenter__openface_T-20_curvflow_22h_ac1e-2` | v1 | disc cut by an open face (curvature flow) | Does ice pass through the do-nothing open face? | No. The film along the whole face is 1.66 eps thick at both eps 0.86 and 0.64 µm; ice is not removed. The condition is ill-posed (BC review 2026-11-05). Led to `-phi_open_bc 0`. | no |
| 0b-A2 | 2026-11-04 | `…/openface_2D_box60um_disc_offcenter__openface_T-20_curvflow_22h_ac1e-2_closed` | v1 | same, 90° face | Control | Matches the half-disc law: empties at 0.45 days against 0.447 predicted. | no |
| 0b-A3 | 2026-11-04 | `…/openface_2D_box60um_disc_interior__openface_T-20_curvflow_22h_ac1e-2` | v1 | interior disc | D_eff of curvature flow | D_eff 1.1401e-15 against eps²/tau_sub 1.1365e-15 m²/s (0.3 %). | no |
| 0b-A4 | 2026-11-04 | `…/openface_2D_box60x30um_cap__openface_T-20_curvflow_22h_ac1e-2` | v1 | wall cap at the wall/open-face corner | Contact line at the corner | Film climbs the open face; ice area grows 75 → 93 µm². | no |
| 0b-B1 | 2026-11-04 | `…/channel_2D_H125um_eps0.39um_deadend__cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-3` | v1 | rectangular channel, xi_v 1e-3 | Baseline; IC transient | Aborted at 1.37e5 s (min dt; stuck at the shared ±0.05 guard edge). IC transient 0.14 % of the ice, so no pre-solve is needed. | no |
| 0b-B2 | 2026-11-04 | `…_xiv1e-4` | v1 | same, xi_v 1e-4 | xi_v invariance | Ice change within 0.6 % of B1 over 8.8e4 s; aborted at 8.8e4 s (phi bounds). | no |
| 0b-B3 | 2026-11-04 | `…_xiv1e-2` | v1 | same, xi_v 1e-2 | xi_v invariance | Within 0.04 % of B1 over 1.37e5 s; completed (ice fully gone at about 2.5e5 s) but overshot phi to 1.18. | no |
| 0b-B4 | 2026-11-04 | `…_xiv1e-3_xiT1e-3` | v1 | same, xi_T = xi_v = 1e-3 | xi_T | Aborted at 5.1e3 s. xi settled at Moure & Fu values (user, 2026-11-05). | no |
| 0b-C1 | 2026-11-04 | `…/throat_2D_L850um_eps0.86um_narrowL_deadend__cryo_T-20_h0.990_iso_5d_ac1e-2` | v1 | two-throat, isothermal | Plug order | Aborted at 3.0e4 s (phi bounds). Mass balance 0.9986. | no |
| 0b-C2 | 2026-11-04 | `…/regolithpore_2D_L814um_eps0.86um_throat_deadend__cryo_T-20_h0.998_iso_5d_ac1e-2` | v1 | regolith pore, h 0.998 | Which ice goes first | Completed with no rollbacks; 85 % of the ice gone in 5 days; mass balance 0.9991. | no (v1) |
| 0b-C3 | 2026-11-04 | `…__cryo_T-20_h0.990_G+100_warmdead_5d_ac1e-2` | v1 | two-throat, warm dead end | Thermal BC | Aborted at 4.3e4 s (min dt). Mass balance 0.9969. | no |
| 0b-C4 | 2026-11-04 | `…__cryo_T-20_h0.990_insdead_5d_ac1e-2` | v1 | two-throat, insulated dead end | Thermal BC | Completed; 44 % of the ice gone; the narrow plug went first, then the wide plug retreated from its open side. Overshot to phi 1.11. T moved ≤ 5e-5 K, so sublimation cooling is negligible at −20 °C. | no |
