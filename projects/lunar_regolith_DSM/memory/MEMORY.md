# MEMORY.md — curated durable facts

<!-- Short, durable, sourced. Not a log (that is project/LOG.md).
     Format: - fact.  (date · source) -->

## Project

- The agentic workflow was restarted from a clean slate; the first attempt is archived under `scratch/2026-10-21_practice_round/` and the `archive/2026-10-21/*` tags.  (2026-10-21 · user)
- Research-folder use is limited to `Literature/` and `Proposals/` in iCloud; everything else is saved in this project or the course folder.  (2026-10-22 · user)

## Solver

- d0 is computed from Kaempfer and Plapp (2009) Eq. 13, d0 = γ·a³/(k_B·T) at temp0; `-d0_sub0` overrides and is announced. Merged into main in a982048b.  (2026-10-22 · user; src/lunar_main.c)
- Open: `monitoring.c` hard-codes d0 = 2.548e-7/T_K, 1.0 % off the Eq. 13 value; comp_eps.py flags this.  (2026-10-22 · comp_eps.py check `d0_matches_solver`)
- The prescribed contact angle (wall free energy) is validated and carried forward; status and open items in `docs/contact_angle/README.md`.  (2026-10-21 · docs/contact_angle)
- `src/assembly.c` and `src/material_properties.c` differ from `enceladus_DSM`; only `src/snes_convergence.c` is still identical.  (2026-10-22 · cmp)
