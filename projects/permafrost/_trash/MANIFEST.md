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
