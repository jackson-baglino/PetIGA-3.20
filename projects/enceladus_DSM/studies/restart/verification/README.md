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

## Gates

| gate | checks |
|---|---|
| clock resumed at snapshot `t` | not 0 — otherwise the two legs overlap and anything measured across the join is wrong |
| dt resumed at snapshot `dt` | not `-delt_t` — the climb-back bug above |
| first new step matches `t` | the continuation lands where the uninterrupted run did |
| first new step matches `dt` | likewise |

The comparison is against the run's *own* next step: leg 1 is run past the
restart point, the snapshot is taken from the second-to-last output, and leg 2
must reproduce leg 1's following row. Measured, it does so to all printed
digits.

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
