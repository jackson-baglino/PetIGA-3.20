# MEMORY.md — curated durable facts

## Project

- Paper 2 of the thesis (paper 1 is `enceladus_DSM`, paper 3 is `lunar_regolith_DSM`). Copied from `lunar_regolith_DSM` on 2026-11-04 with the model frozen at tag `lunar_cryosuction/model-freeze-v1` (commit 147a49f0); physics frozen, setup open, see `FREEZE.md`.  (2026-11-04 · user)

## Solver (copied from lunar_regolith_DSM 2026-11-04)

- d0 = γ·a³/(k_B·T) (Kaempfer and Plapp 2009, Eq. 13) is computed once at start-up from temp0 and used for the whole run, including monitoring.c; it is not an input (`-d0_sub0` removed).  (2026-10-22 · user; b927b982)
- comp_eps.py computes d0 with the solver's formula and constants (1.0152e-9 m at −20 °C). Before 2026-10-22 it used γ·V_m/(R·T) = 1.0166e-9 m, so geometry files made earlier carry an eps about 0.14 % larger than comp_eps now gives (kinetic-bound cases).  (2026-10-22 · user)
- The wall free-energy term recovers the static (Young) contact angle. It is not fully validated: the implementation does not apply something like a free-flow boundary condition at the contact line; to be discussed. Status in `docs/contact_angle/README.md`.  (2026-10-22 · user)
- `SSA_evo.dat` gets a row every time step by design; `-outp` / `-t_interv` control only the `sol_*.dat` snapshots.  (2026-11-01 · user)
