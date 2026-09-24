# k_eff disk array: ice-fraction sweep

The disk ladder (`../disk/`) compares the arithmetic and tensor conductivity laws
against Rayleigh's sharp solution at a single ice fraction, `f = 0.196`. This
sweep repeats that comparison at six disk radii, so each interface width
becomes a curve `k_eff(f)`.

| | |
|---|---|
| radii | `R = 125, 175, 250, 300, 350, 400 µm` in a 1 mm periodic cell, giving `f = 0.049 … 0.503` |
| widths | `eps/R = 0.04, 0.02, 0.01`. **0.02 is production**: `eps = 1 µm` on `R_ave = 50 µm`. |
| laws | arithmetic (`-keff_interp arith`) and tensor |
| mesh | `eps/h = 5.12`, as in both ladders, giving `N = 320 … 4096` |
| jobs | 36. Each is `-keff_only` with no time integration. |

## Why sweep f

A single point can show that the bias exists. It cannot show that the bias
**depends on the microstructure**. The first-order theory
(`../disk/keff_disk_analytic.py`) predicts that at `eps/R = 0.02` the
arithmetic law reads +1.8% high at `f = 0.05` and +24% high at `f = 0.50`. If
the bias were a fixed fraction, it would cancel in `k(t)/k(0)`. Because it
varies with the geometry, it does not cancel.

The `R = 250 µm` rows repeat the ladder rungs L/100, L/200 and L/400, so they
double as a consistency check against `../disk/keff_disk.csv`.

## Running it (HPC)

Before submitting, push from the Mac and pull on the cluster.

```bash
cd $HOME/PetIGA-3.20/projects/enceladus_DSM      # wherever the repo lives on the cluster
venv_enceladus/bin/python studies/keff_sharp_limit/disk_sweep/make_sweep_spec.py
./scripts/HPC/submit_batch.sh --tag keff_disk_fsweep \
    --tests-file studies/keff_sharp_limit/disk_sweep/sweep_tests.txt -- --time=0-01:00:00
```

The run folder name still says `R250um`, because every job uses that geometry
file with `-RCice` overridden per job. The `--label`
(`arithmetic_R125_e0.02`, …) is what identifies each run.

`submit_batch.sh` now sizes each job from its own `-Nx/-Ny`. Before this
change it used the geometry file's 256², which would have given every job 2
ranks. A `-keff_only` job is sized on 1 DoF per node, because the cost is the
scalar corrector solve. The largest job (4096²) gets 168 ranks.

## After it finishes

Download the batch folder. Only the small CSVs are needed. Then run:

```bash
venv_enceladus/bin/python studies/keff_sharp_limit/disk_sweep/collect_sweep.py <batch_dir>
venv_enceladus/bin/python docs/keff_tensor_slides/make_figures.py
```

`collect_sweep.py` writes `keff_disk_sweep.csv` and exits non-zero if any run
is missing, fails square symmetry, or built a different disk than the one
requested (`phi_bar` against `f(1 + π²ε²/3R²)`).
