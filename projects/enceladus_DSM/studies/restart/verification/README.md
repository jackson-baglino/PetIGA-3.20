# Restart verification

Does `-initial_cond` resume a run seamlessly?

```bash
./studies/restart/verification/verify_restart.sh
```

## What restart does

```bash
./enceladus_dsm <same opts as leg 1> -initial_cond <rundir>/sol_NNNNN.dat
```

**That is the whole invocation.** The clock and the time step are both read
from the snapshot's own `SSA_evo.dat`, which carries one row per step:

```
ssa/eps   tot_ice   t   step   dt   tot_air   tot_rhov   tot_mass
```

so `NNNNN` from the filename indexes straight into it. `-t_start` and
`-dt_start` override, but should not normally be needed — and a value that is
never typed is a value that cannot be typed wrong.

## Why dt is adopted, not reset

A restart that begins at `-delt_t` (1e-4 s by default) has to climb six orders
of magnitude back to the working step size through the NRmin/NRmax growth
heuristic — roughly 145 full nonlinear solves before any new physics happens.
On a production mesh that is hours of compute to get back to where the run
already was. Adopting the snapshot's dt makes a continuation cost what
continuing would have cost.

## Step count, output schedule and the log (2026-09-29)

A resume carries on **counting steps from the snapshot's own number**
(`TSSetStepNumber`), resumes the **output schedule** at the next time the
uninterrupted run would have written, and **appends** to `SSA_evo.dat`
(skipping the resumed step's row, which is already there). Before this, a
resume restarted at step 0 and opened `SSA_evo.dat` for write -- so resuming
in place (`RESUME_INTO`, the default in `resume_batch.sh`) overwrote the
first leg's `sol_00000.dat`, `sol_00001.dat`, ... and truncated its log.

`resume_batch.sh` moves whatever the first leg wrote past the resume point
into `<run>/abandoned_after_step<N>_<timestamp>/` first: the continuation
recomputes that stretch, and keeping both would give two states per step.

## Gates

`verify_restart.sh` runs the case uninterrupted (A) and stopped-then-resumed
in place (B), and requires B's directory to match A's:

| gate | checks |
|---|---|
| one row per step, same steps as A | no step-0 restart, no duplicated resume row |
| same t and dt on every step (1e-9) | the continuation follows the same trajectory |
| same totals (1e-6) | ice, air, vapour, mass -- U is reloaded, its time derivative is not |
| same `sol_*.dat` names | nothing overwritten, nothing written off-cadence |
| resumed dt is the working dt | not `-delt_t` -- the climb-back bug below |

## Restarting from a killed job

Use the **second-to-last** `sol_*.dat`, not the last. `SSA_evo.dat` is flushed
every step so it is always intact, but a snapshot the scheduler interrupted
mid-write can be truncated, and `IGAReadVec` would either error or — worse —
load a short vector. One extra step is cheap insurance.

```bash
ls <rundir>/sol_*.dat | sort | tail -2 | head -1
```

The restart leg must use the **same geometry** (mesh, `-p`, `-C`, `-dof`,
domain) as the leg that wrote the snapshot. `IGAReadVec` checks the vector
length, which catches a wrong geometry file. It does **not** check `-eps`, the
kinetics, or `-dtmax`, all of which may legitimately change on restart.
