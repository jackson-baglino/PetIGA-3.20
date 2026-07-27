# Study: Icy lunar regolith (Paper 1)

**Scientific question.** How does dry-snow-style metamorphism (sublimation /
vapor redistribution of ice) proceed when the ice occupies the pore space of an
inert granular regolith? What controls ice redistribution, neck formation, and
loss rates in an icy lunar-regolith analog?

**Approach: implicit pore domain** (`implicit_pore_domain/`). The regolith is
represented by the **domain geometry** — baked into the mesh boundary — so the
solver stays at 3 DOF (ice, T, vapor) and needs no model changes.
Branch: `exp/regolith-implicit-pore-domain`.

The competing **explicit sediment phase** approach (a 4th phase-field DOF φ_s
with ∂φ_s/∂t = 0 and a triple-well potential) was set aside on 2026-07-24 and
retired on 2026-07-27. It failed because the variational-inequality bound
constrains the ice DOF, but air is a *derived* field (1 − φ_i − φ_s), so it
overshoots at the triple junction — `DIVERGED_FUNCTION_DOMAIN` at step 26.
Making it work needs a constrained-sum / obstacle multiphase solver. The
three-phase formulation is now gone from both sibling projects; its last state
is at tag `archive/dry_snow_metamorphism-legacy` and in this repo's history.

**Shared solver.** Consumes the model in `../../src` + `../../preprocess` +
`../../postprocess`, which is byte-identical to the sibling `enceladus_DSM`.

**Parameter regime.** Sub-mm to few-mm domains; lunar surface/near-surface
temperatures. ε and mesh sized via `../../preprocess/comp_eps.py`.
