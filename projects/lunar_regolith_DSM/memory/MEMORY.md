# MEMORY.md — curated durable facts

<!-- Short, durable, sourced. Not a log (that is project/LOG.md).
     Format: - fact.  (date · source) -->

## Project

- The agentic workflow was restarted from a clean slate; the first attempt is archived under `scratch/2026-10-21_practice_round/` and the `archive/2026-10-21/*` tags.  (2026-10-21 · user)
- Research-folder use is limited to `Literature/` and `Proposals/` in iCloud; everything else is saved in this project or the course folder.  (2026-10-22 · user)

## Solver

- d0 = γ·a³/(k_B·T) (Kaempfer and Plapp 2009, Eq. 13) is computed once at start-up from temp0 and used for the whole run, including monitoring.c; it is not an input (`-d0_sub0` removed).  (2026-10-22 · user; b927b982)
- comp_eps.py still sizes eps with d0 = γ·V_m/(R·T), 0.14 % above the solver's value; whether to align it is open.  (2026-10-22 · comp_eps.py check `d0_matches_solver`)
- The wall free-energy term recovers the static (Young) contact angle. It is not fully validated: the implementation does not apply something like a free-flow boundary condition at the contact line; to be discussed. Status in `docs/contact_angle/README.md`.  (2026-10-22 · user)
- `src/assembly.c` and `src/material_properties.c` differ from `enceladus_DSM`; only `src/snes_convergence.c` is still identical.  (2026-10-22 · cmp)
