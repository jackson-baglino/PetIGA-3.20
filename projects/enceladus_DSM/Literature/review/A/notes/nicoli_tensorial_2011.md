# nicoli_tensorial_2011 — Nicoli, Plapp & Henry 2011, Phys. Rev. E 84:046707 (doi 10.1103/PhysRevE.84.046707)

Whole paper read (pp. 046707-1 to -6).

**Summary.** Steady two-phase transport j = −M(φ)∇V, ∇·j = 0 through a fixed diffuse-interface structure. Shows two thin-interface effects of order W: a surface excess current (removed by direct/arithmetic interpolation, Eq. 9) and an interface (Kapitza) resistance (removed by inverse/harmonic interpolation, Eq. 16). Both are removed only by a tensor mobility M = M⊥ n⊗n + M∥(1 − n⊗n) (Eq. 17) with M∥ arithmetic and M⊥ harmonic. Validated in 2D: disk inclusion (M1 = 1, M2 = 1/2) and two half-disks (M2/M1 = 0.1), finite differences on MAC grid.

**Findings relevant to the manuscript**
- p.1 abstract: "both effects can be eliminated simultaneously only if the coefficient of proportionality between flux and potential gradient (mobility) is allowed to become a tensor in the interfaces."
- p.2: assumes "the profile of φ is odd with respect to the point φ = 1/2" (φ(x) = 1 − φ(−x)), "true in all standard phase-field models".
- p.3: valid for curved interfaces "with a local radius of curvature R satisfying R ≫ W ... at least up to second order in ε = W/R".
- p.4 Fig. 2: "the direct and inverse interpolations converge only linearly with ε to this limiting value, whereas the tensorial interpolation suppress linear thin-interface effects, which leads to a convergence that is almost quadratic in ε."
- Table I: first-order coefficients c1 = 0.0739 (direct), −0.152 (inverse), −0.00135 (tensorial).
- p.1: lists "conduction of electrical current or heat through a two-phase material with different conductivities" and tomographic structures as targets.

**Claims it supports / contradicts / qualifies**
- Draft Eqs. (16)–(17) (L404–409), tensor law with K∥ arithmetic and 1/K⊥ harmonic: **supported** (Nicoli Eqs. 9, 16, 17).
- Draft L418 "excess resistance normal to it and excess conductance tangent to it, both of order ε for a scalar conductivity law ... The tensorial law removes both at the first order": **supported** (p.2–3; Fig. 2).
- Outline L783–786 "define the law and show the cancellation for a planar interface at contrasts of 2 and 10": **supported** (analysis planar p.2–3; numerical tests M2/M1 = 1/2 and 0.1, p.4).
- Caveat to add (Limitations or Sec 2.3): Nicoli's result holds for R ≫ W (p.3); in the sintered packings the throats/necks are not much wider than the interface (cf. project memory: band 9.2–18.4 μm vs throats 4–12.6 μm), so first-order cancellation is not guaranteed there. One clause.
- Nicoli's analysis is for a fixed (stationary) φ — consistent with outline L790 "cell problem is solved on frozen snapshots of φ".

**Candidate citations / claims**
- Sec 2.3: "the cancellation holds where the interface is flat on the scale of its width (Nicoli et al., 2011); at the narrowest throats it is approximate."
