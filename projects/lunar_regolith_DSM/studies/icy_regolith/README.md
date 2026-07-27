# Study: Icy lunar regolith (Paper 1)

**Scientific question.** How does dry-snow-style metamorphism (sublimation /
vapor redistribution of ice) proceed when the ice occupies the pore space of an
inert granular regolith? What controls ice redistribution, neck formation, and
loss rates in an icy lunar-regolith analog?

**Two competing approaches** to representing the regolith. Both target the same
result; one will be chosen for the paper.

| Approach | Directory | Regolith representation | DOF |
|---|---|---|---|
| Implicit | `implicit_pore_domain/` | Pore-space **domain geometry** (regolith baked into the mesh boundary) | 3 (ice, T, vapor) |
| Explicit | `explicit_sediment_phase/` | 4th **phase-field DOF** φ_s with ∂φ_s/∂t = 0, triple-well potential | 4 |

**Shared solver.** Both consume the master model in `../../src` +
`../../preprocess` + `../../postprocess`. The implicit approach needs no solver
changes (branch `exp/regolith-implicit-pore-domain`); the explicit approach adds
a 4th field and a triple well (branch `exp/regolith-explicit-sediment-phase`).

**Parameter regime.** Sub-mm to few-mm domains; lunar surface/near-surface
temperatures. ε and mesh sized via `../../preprocess/comp_eps.py`.

**Status (2026-07-21).** Not started. Effort 2's surface energies are the known
prior blocker — parameterize by contact angle θ (only γ_is > (γ_iv/2)(1−cos θ)
is a real constraint); the beta-eliminated equations are still to be supplied.
Start explicit from the simplest geometry (one slab, one grain) and gate on
3-phase mass conservation before adding complexity.
