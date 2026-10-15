# nicoli_tensorial_2011 — Nicoli, Plapp & Henry (2011), Phys. Rev. E 84, 046707

**Summary.** Analysis of steady two-phase transport (flux = −M(φ)∇V, ∇·j = 0) through a fixed diffuse-interface structure. For a planar interface with an odd-symmetric profile, a scalar mobility produces a surface excess current (zero only for arithmetic interpolation) and a potential jump/Kapitza resistance (zero only for harmonic interpolation). A tensorial mobility M = M⊥ n⊗n + M∥(1−n⊗n) with harmonic M⊥ and arithmetic M∥ removes both. Validated numerically (finite differences, MAC grid, SOR) on a disk inclusion and on two half-disks.

**Findings relevant to the manuscript**
- p.1 abstract: "the finite thickness of the interfaces induces two effects ... a surface excess current and a potential jump at the interfaces ... both effects can be eliminated simultaneously only if the ... mobility is allowed to become a tensor in the interfaces."
- Surface current ∝ W, zero for direct interpolation M = M1φ + M2(1−φ) (Eqs. 7–9, p.2); interface resistance ∝ W, zero for inverse interpolation 1/M = φ/M1 + (1−φ)/M2 (Eqs. 15–16, p.3); tensor Eq. (17), p.3.
- Validity, p.3: "for a sufficiently smooth interface (... R ≫ W), a local curvilinear coordinate system can be defined in which the above relations remain valid at least up to second order in ε = W/R."
- Numerical tests, p.4: disk "M1 = 1, M2 = 1/2" (contrast 2); half-disks "M2/M1 = 0.1" (contrast 10). Fig. 2: "the direct and inverse interpolations converge only linearly with ε ..., whereas the tensorial interpolation suppress linear thin-interface effects, which leads to a convergence that is almost quadratic in ε." Table I: c1 = 7.4×10⁻², −1.5×10⁻¹, −1.35×10⁻³ (direct, inverse, tensorial).
- Applicability to tomographic data, p.1: smoothed phase-field representations of µCT data could be used directly for transport.
- Moving interfaces not solved: p.5 "this methodology is still not capable of removing all thin-interface effects from phase-field models with a nonstationary φ".

**Claims supported / qualified**
- Draft L400–410 (Following Nicoli et al., tensor K with K∥ arithmetic, 1/K⊥ harmonic): **supported** — identical construction (Eq. 17, p.3).
- Draft L418 "excess resistance normal to it and excess conductance tangent to it, both of order ε for a scalar conductivity law ... The tensorial law removes both at the first order": **partially** — each is O(W), but a scalar law can remove one of them (arithmetic removes the tangential excess, harmonic the normal resistance); the point is that no scalar law removes both. Tensorial removal at first order: supported (Fig. 2, Table I).
- Draft L418 "The diffuse interface formulation tends to the sharp interface model as the interface width ε → 0": supported (Fig. 1, all interpolations converge to J0). Note Nicoli's W is the interface width; the manuscript's ε is the interface parameter — fine, but don't call ε the "width".
- Outline L783–786 "Nicoli, Plapp & Henry (2011) define the law and show the cancellation for a planar interface at contrasts of 2 and 10": **partially** — the cancellation is derived analytically for a planar interface (any contrast); the numerical demonstrations at contrasts 2 and 10 are curved geometries (a disk; two half-disks), with validity requiring R ≫ W.
- Outline §2.3 "a frozen microstructure" (cell problem steady): matches Nicoli's scope (fixed, immobile structure, p.2: "Since we consider a fixed and immobile two-phase structure, the phase field is independent of time").

**Candidate additions.** In §2.3 state the condition R ≫ ε·(band factor) for the tensor law's accuracy (Nicoli p.3), and connect it to the paper's own limitation that pore throats narrower than the 9.2 µm band are unresolved (outline L1368–1369) — this is exactly where Nicoli's curvature condition fails.
