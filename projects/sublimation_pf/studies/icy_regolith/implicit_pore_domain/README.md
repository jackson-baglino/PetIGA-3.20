# Implicit pore-domain approach (Effort 1)

**Question.** How does ice adhere to and metamorphose within lunar-regolith pore
space, as a function of pore geometry? Here the regolith is modelled
*implicitly* — it is the deformed top/bottom boundary of the domain, not a
simulated field. Ice / vapour / temperature evolve in the pore space between the
regolith walls (3 DOF, no solver changes). Interior regolith grains are not
representable in a single tensor-product patch; that is the explicit approach
(`../explicit_sediment_phase/`, Effort 2).

## Geometry

A two-sided pore channel whose walls are lined with regolith grains drawn from a
grain-size distribution (median radius ~50 µm, `docs/material_parameters.md`
§2.3), with a central **throat pinch** (the defining pore constriction). Built by:

    preprocess/build_geometry_regolith_pore.py

It writes one shared mesh `inputs/geometry/regolith_pore.dat` (the walls) and a
per-strategy `.opts` (the ice). Deterministic; asserts every geometric
constraint (throat gap, no buried/overlapping ice, `eps/R` resolution) before
writing. See the header of that script for all knobs (`R_MED`, `N_BUMPS`,
throat multipliers, wall seeds).

## Ice placement is the study variable

How ice sits in the pore is exactly what we're probing, so the initial ice
configuration is a swappable `--ice-placement` strategy, each testing a
different adhesion hypothesis:

| Strategy (`.opts`) | grains | Hypothesis it probes |
|---|---|---|
| `flank_caps` → `2D_regolith_pore.opts` | 6 | Ice caps on grain tops; throat bare. Baseline. |
| `throat_bridge` → `2D_regolith_pore_throat.opts` | 7 | Ice suspended in the throat — does it persist/adhere at the constriction? |
| `pore_lining` → `2D_regolith_pore_lining.opts` | 15 | Small ice grains seated in wall troughs (reentrant pore corners) — distributed adhesion. |

Preview each config (walls + ice) at `preprocess/regolith_pore{,_throat,_lining}.png`.
Add a strategy by registering a function in `ICE_STRATEGIES`.

## Running

    # from the project root
    ./scripts/Studio/run_permafrost.sh 2D_regolith_pore        30day_T-20_h1.00_arrh
    ./scripts/Studio/run_permafrost.sh 2D_regolith_pore_throat 30day_T-20_h1.00_arrh
    ./scripts/Studio/run_permafrost.sh 2D_regolith_pore_lining 30day_T-20_h1.00_arrh

(3 DOF; `run_permafrost.sh` sizes the rank count from the `# DOF_GRID` comment.)

## Regenerating / varying the geometry

    python3 preprocess/build_geometry_regolith_pore.py --ice-placement flank_caps
    python3 preprocess/build_geometry_regolith_pore.py --ice-placement throat_bridge --tag throat
    python3 preprocess/build_geometry_regolith_pore.py --ice-placement pore_lining  --tag lining

To vary the *pore* geometry (wall roughness, throat tightness, GSD), edit the
module constants and re-run; the `.dat` and all `.opts` regenerate.

## Status & caveats

- **`eps` is loose** (reused ~8.58e-7 at this scale). **Recompute with
  `preprocess/comp_eps.py`** for the actual run temperature before any
  production run — the builder prints this reminder.
- First gate: a short run that passes mass conservation (`postprocess/plot_mass.py`).
- Then: sweep ice placement × throat tightness, and quantify the resulting
  adhesion morphology (which configs coarsen toward the walls vs. the throat,
  SSA evolution, ice–wall contact length).
