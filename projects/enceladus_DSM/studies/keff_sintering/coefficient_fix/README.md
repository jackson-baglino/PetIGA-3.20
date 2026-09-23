# The band-interpolation law on a real packing

`-keff_interp` selects how the conductivity is interpolated across the diffuse
band in the steady-state cell problem. It defaulted to `arith` until
2026-09-23 and defaults to `tensor` now. **Every run before that date, the
2026-09-16 pilot included, carries the arithmetic law.**

This directory measures what the change is worth on the packing the campaign
actually runs.

## Why the existing ladders do not already answer this

`studies/keff_sharp_limit/` verifies the tensor law where the answer is known
in closed form:

| ladder | geometry | arith | tensor |
|---|---|---|---|
| `verification/` | planar slab | `k_⊥` +59.4% → +3.8% down the `eps` ladder | flat at the sharp value to **4e−7** |
| `disk/` | isolated disk array, Rayleigh | +18.19% → +1.93% | +1.212% → **+0.006%**, order ≈ 2.5 |

**Neither has a contact.** The packing's grains are seated in exact tangency,
so at `t = 0` the conduction path between two grains is a single point. That
is where the arithmetic law does its worst: the diffuse band bridges the
contact with partially-conducting material, inflating `k_eff(0)` for a
structure that ought to barely conduct — and thereby suppressing the
*relative* rise, which is the campaign's headline.

At the pilot's `eps/R_ave = 0.02` the disk ladder puts the arith bias at about
**+8%**. That is a **lower bound** on the packing, not an estimate of it: it is
the isolated-interface part of the error with the contact effect excluded by
construction.

## The "before"

Measured from the four pilot runs' own in-line `k_eff.csv`
(`~/SimulationResults/HPC_results/enceladus_DSM/GrainPackingSintering/batch_2026-09-16__13.16.11_pilot_keff/`;
φ=0.325, `L/R_ave`=40, `R_ave`=50 µm, `eps`=1 µm, `Nx=Ny=2829`, 241 ranks,
`alpha_c`=1e-3, T=−20 °C, 30 d):

| seed | `k_iso(0)` | `k_iso(end)` | rise |
|---|---|---|---|
| 1 | 0.6312 | 0.7408 | +17.4% |
| 2 | 0.6289 | 0.7533 | +19.8% |
| 3 | 0.6462 | 0.7905 | +22.3% |
| 4 | 0.6262 | 0.7438 | +18.8% |

Ensemble **+19.6%, sd 2.1 points, SEM 1.1**.

A numpy finite-volume twin (`studies/packing_design/cell_solve.py`, commit
`0a7dac1`) put the sharp-law rise at +57.0% on seed 1, but it ran on the 945²
`.vts` output grid rather than the 2829² solve grid and left no driver. Treat
it as indicative of direction and scale only — the replay below supersedes it.

## Producing the "after"

The phase field on disk is unaffected by the law, so this is a replay, not a
rerun:

```bash
# 1. on the cluster — always dry-run first, the snapshot count sets the bill
./scripts/HPC/submit_keff_replay.sh --dry-run --laws "tensor sharp" \
    --roots $SCRATCH/enceladus_DSM/batch_2026-09-16__13.16.11_pilot_keff \
            $SCRATCH/enceladus_DSM/packing_2D_pilot_phi0.325_Rave50um_LR40_seed*_L2mm_eps1000nm_perxy_T-20
# 2. drop --dry-run to submit
# 3. download, then:
venv_enceladus/bin/python studies/keff_sintering/coefficient_fix/compare_laws.py <batch_dir>
```

Give **both** roots, and note the `seed*` glob: a run resumed after a walltime
kill leaves its later snapshots in the single-run tree `<geom>/<ts>_..._resume_job<id>/`,
not in the batch parent, and **all four** pilot seeds were resumed. Naming one
seed there would silently replay half of the other three. Read the `--dry-run`
listing and check it shows eight directories, not four. The script scans for `sol_*.dat` and submits one job per
directory it finds; `compare_laws.py` concatenates the per-leg CSVs by time.

`k_eff_tensor.csv` and `k_eff_sharp.csv` are written *alongside* the original
`k_eff.csv`, which is left untouched — that is what makes the before/after a
comparison rather than a replacement.

## What to read from it

1. **The rise**, `arith → tensor`. This is the correction to the campaign's
   headline and it re-scales every threshold downstream.
2. **`k_eff(0)` specifically.** `ROADMAP.md:260` forbids using `k_eff(t=0)` as
   a baseline *because* arith pre-welded tangent grains. If the tensor `k(0)`
   drops substantially, that rule is obsolete and `CAMPAIGN.md` should say so.
3. **`tensor` vs `sharp`.** Two completely different routes to removing the
   same bias — one analytic, one by thresholding the geometry. Agreement is
   the evidence that either is right. A gap above ~5% is a result in its own
   right and everything downstream waits on it. The driver prints this check.
4. **`ksp_its`.** The tensor operator is locally anisotropic inside the band
   (`K_∥/K_⊥ ≈ 29` at `φ = 0.5`) and costs ~2.5× the CG iterations on the disk
   ladder. A larger rise on the packing is a real campaign cost signal.

## Files

| | |
|---|---|
| `compare_laws.py` | discovery, leg merge, table, cross-check, figure |
| `compare_laws.csv` | one row per (seed, law) |
| `compare_laws.png` | (a) absolute `k_eff(t)`, (b) normalised to each law's own `t = 0` |

Until the replay lands, the CSV and figure hold the arith baseline only, and
the driver says so on stdout.
