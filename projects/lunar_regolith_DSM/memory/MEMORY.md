# MEMORY.md — curated durable facts

<!-- Short, durable, sourced. Not a log (that is project/LOG.md).
     Format: - fact.  (date · source) -->

## Project

- The agentic workflow was restarted from a clean slate; the first attempt is archived under `scratch/2026-10-21_practice_round/` and the `archive/2026-10-21/*` tags.  (2026-10-21 · user)
- Research-folder use is limited to `Literature/` and `Proposals/` in iCloud; everything else is saved in this project or the course folder.  (2026-10-22 · user)
- `memory/semantic/units.md` verified 2026-10-22 by the units-verification workflow (run 20261022-155335): 54 literature rows checked against page images; the value column is the printed value.  (2026-10-22 · workflow; user-approved)

## Solver

- d0 = γ·a³/(k_B·T) (Kaempfer and Plapp 2009, Eq. 13) is computed once at start-up from temp0 and used for the whole run, including monitoring.c; it is not an input (`-d0_sub0` removed).  (2026-10-22 · user; b927b982)
- comp_eps.py computes d0 with the solver's formula and constants (1.0152e-9 m at −20 °C). Before 2026-10-22 it used γ·V_m/(R·T) = 1.0166e-9 m, so geometry files made earlier carry an eps about 0.14 % larger than comp_eps now gives (kinetic-bound cases).  (2026-10-22 · user)
- The wall free-energy term recovers the static (Young) contact angle. It is not fully validated: the implementation does not apply something like a free-flow boundary condition at the contact line; to be discussed. Status in `docs/contact_angle/README.md`.  (2026-10-22 · user)
- `src/assembly.c` and `src/material_properties.c` differ from `enceladus_DSM`; only `src/snes_convergence.c` is still identical.  (2026-10-22 · cmp)
