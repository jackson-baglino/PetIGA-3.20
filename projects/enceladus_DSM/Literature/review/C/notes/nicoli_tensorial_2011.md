# nicoli_tensorial_2011 — Nicoli, Plapp & Henry, Phys. Rev. E 84, 046707 (2011)
DOI (PDF p.1): 10.1103/PhysRevE.84.046707. Read in full (6 pp).

## Summary
Analysis of steady two-phase transport j = −M(φ)∇V, ∇·j = 0 through a fixed diffuse-interface structure. For a planar interface they show two O(W) thin-interface effects: an excess surface current (eliminated by direct/arithmetic interpolation, Eq. 9) and an interface (Kapitza) resistance (eliminated by inverse/harmonic interpolation, Eq. 16). A tensor mobility M = M⊥ n⊗n + M∥(1 − n⊗n) with M⊥ harmonic and M∥ arithmetic (Eq. 17) removes both. Numerical tests in 2D (finite differences, MAC grid): a disk inclusion (M1 = 1, M2 = 1/2, R = 0.25) and two half-disks (M2/M1 = 0.1).

## Findings relevant to the manuscript
- "both effects can be eliminated simultaneously only if the coefficient of proportionality between flux and potential gradient (mobility) is allowed to become a tensor in the interfaces" (abstract p.1/046707-1).
- Each scalar interpolation removes one effect: direct (Eq. 9) cancels the surface conductivity; inverse (Eq. 16) cancels the interface resistance; "Since these interpolations are mutually exclusive, it seems as if necessarily one of the two effects must remain nonzero" (p.3/046707-3).
- Convergence: "the direct and inverse interpolations converge only linearly with ε to this limiting value, whereas the tensorial interpolation suppress linear thin-interface effects, which leads to a convergence that is almost quadratic in ε" (p.4/046707-4); Table I: c1 = 0.0739 (direct), −0.152 (inverse), −0.00135 (tensorial).
- Validity: planar result extends to curved interfaces "with a local radius of curvature R satisfying R ≫ W ... at least up to second order in ε = W/R" (p.3/046707-3).
- Motivation explicitly includes heat conduction through two-phase materials and tomographic voxel data smoothed by a phase field (p.1/046707-1).
- Not covered: moving interfaces ("still not capable of removing all thin-interface effects from phase-field models with a nonstationary φ", p.5/046707-5) — consistent with the paper solving the cell problem on frozen φ snapshots.

## Claims supported / contradicted / qualified
- Draft L400–418, Eqs. (16)–(17) (K∥ arithmetic, K⊥ harmonic, tensor): supported (Eqs. 9, 16, 17).
- Draft L418 "excess resistance normal to it and excess conductance tangent to it, both of order ε for a scalar conductivity law": partially. Each effect is O(W) (p.2–3), but a scalar law can cancel one of them (arithmetic cancels the tangential excess; harmonic cancels the normal resistance); a scalar law cannot cancel both. Suggested: "a scalar law can remove one of the two but not both".
- Draft L418 "The tensorial law removes both at the first order": supported (p.3, Fig. 2, Table I).
- Outline L783–786 "show the cancellation for a planar interface at contrasts of 2 and 10": partially. The cancellation is derived analytically for a planar interface (any contrast); the numerical checks are a 2D disk inclusion (contrast 2) and two half-disks (contrast 10). Rephrase: "derive the law for a planar interface and test it on a disk inclusion and two half-disks (contrasts 2 and 10, in 2D)".
- Outline L1715–1718: bib key still old style: confirmed key `nicoli_tensorial_2011`.

## Candidate additions
- Section 2 (tensor law): the law holds to O(ε) for R ≫ W (p.3); at sintering necks the curvature radius can approach the interface width, so the first-order cancellation is only guaranteed where R ≫ W. One clause in Limitations.
