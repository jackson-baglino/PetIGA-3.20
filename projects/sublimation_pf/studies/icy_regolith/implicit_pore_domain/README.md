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

## Parameters — eps is proper, but temperature-locked

`eps = 8.5840e-7` is **not** an arbitrary/loose value: it is exactly what
`preprocess/comp_eps.py` recommends for **T = −20 °C, α_c = 1.341e-2**, and it
reproduces the validated `2D_ripening_two_sided` reference run bit-for-bit
(same `beta_sub0=5.9216e5`, `d0_sub0=1.0166e-9`). The binding K&P constraint is
the temperature-dependent **kinetic bound**, so eps is grain-size independent
here — all three strategies pass `eps/R_ave < 5%` (3.0%, 3.0%, 4.8%).

**Consequence:** these geometries are valid **only at −20 °C**, paired with a
−20 °C experiment (`30day_T-20_h1.00_arrh`, which sets the matching kinetics).
For any other temperature — notably lunar PSR temperatures (40–120 K,
`docs/material_parameters.md`) — **recompute eps** and regenerate:

    python3 preprocess/comp_eps.py --Lx 8.14e-4 --Ly 2.6933e-4 \
        --Rave <smallest_ice_R> --T0 <degC> --alpha <alpha_c>
    # then set T0_C / ALPHA_C / EPS in build_geometry_regolith_pore.py and re-run

## Status & first gate

- First gate: a short −20 °C run that passes mass conservation
  (`postprocess/plot_mass.py`), confirming the geometry + IC behave.
- Then: sweep ice placement × throat tightness, and quantify the resulting
  adhesion morphology (which configs coarsen toward the walls vs. the throat,
  SSA evolution, ice–wall contact length).
