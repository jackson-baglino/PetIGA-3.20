# _trash — staged for deletion, pending review

Files here have been `git mv`'d out of the tree so the move is reviewable in a
single commit and fully recoverable from git history. Nothing here is deleted
yet. Once reviewed, a follow-up commit runs `git rm -r _trash/`.

To rescue a file: `git mv _trash/<path> <original path>`.

---

## Batch 1 — dead `-ic_type` configs (2026-07-21)

Dropped in the 2026-06-13 two-phase fork, which removed the sediment degree of
freedom. The two-phase solver (`src/permafrost2.c`) dispatches only
`two_ice_grains_boundary`, `ice_slab`, `single_ice`, `multi_grains` in 2D/3D
and `SETERRQ`s on anything else. These files broke the regression suite: every
`DEFAULT_TESTS` entry that named one aborted at IC setup.

### `inputs/geometry/` — 2D, abort at startup (invalid ic_type)
| file | ic_type |
|---|---|
| `2D_ice_sed_pair.opts` / `_hires` | `ice_sed_pair` |
| `2D_separated_grains.opts` / `_hires` | `enclosed` |
| `2D_touching_grains.opts` / `_hires` | `enclosed` |
| `2D_single_sed.opts` | `single_sed` |

### `inputs/geometry/` — 1D, fall through but semantically obsolete
`1D_ice_sed_pair.opts` / `_hires` (`ice_sed_pair`), `1D_separated_grains.opts`
/ `_hires` (`enclosed`), `1D_touching_grains.opts` / `_hires` (`enclosed`).
These do not abort (1D has a default fall-through) but describe sediment-era
grain configurations that no longer have meaning in the two-phase model.

### `inputs/tests/` — entire directory (8 files, all abort)
`test1_IceCap` (`ice_cap`), `test2_CapillaryBridge` (`capillary`),
`test4_ContactSedPair` (`contact_sed`), `test5_SlabAndGrains` (`slab_and_grains`),
each with a `_hires` variant. All request removed ic_types.

### `tmp/advisor_slides.tests`
2026-05-25 `<geom>:<exp>` list for advisor-slide batches; 5 of 6 entries name
dead geometries and no script consumes the file.

**Companion edits (not in _trash — applied in place):**
- `scripts/Studio/run_batch_tests.sh` and `scripts/HPC/submit_full_suite.sh`:
  `DEFAULT_TESTS` and `ADVISOR_SLIDES_GEOMS` repaired to reference only live
  geometries.
- `inputs/README.md`: rewritten for the two-phase model.
- `scripts/check_ic_types.sh`: new guard that validates every `.opts` ic_type
  against the source dispatch list (and checks the `# DOF_GRID:` comment).

---

## Batch 2 — Gibbs-Thomson / d0_GT removal (2026-07-21)

GT curvature is no longer pursued. The code was removed in a preceding commit
(flat-interface `rho_vs` everywhere; `Curvature()` and the phi Hessian gone).
`-d0_GT` is now an unregistered option that PETSc silently ignores.

### `preprocess/` — GT-only tooling
`plot_alpha_gt_sweep.py`, `plot_gt_temperature_dependence.py`.

### `inputs/experiment/` — 21 GT-specific experiments
Everything matching `*_GT`, `*_GTphys`, `*_GTamp`, `*_nopen_GT*`: the GT
parameter-sweep diagnostics (`1day_T-20_nopen_GT*`) and the long-duration
GT-ripening runs (`{2wk,6mo,90day,365day,3yr,5yr,10yr,100yr}…_GTphys`,
`30day_T-5_h1.00_GTamp`, `30day_T-20_h1.00_nopen_GT*`). Without a curvature
term these have no Ostwald-ripening driver, so they are moot as written.

**Companion edits (in place):** the inert `-d0_GT` option line was stripped
from 30 substantive experiment files (Molaro, Tgrad, ripening, collapse,
sinter, dt studies) that set `-d0_GT 0` but are not GT-specific. Their physics
comments referencing capillarity are left as documentation.
