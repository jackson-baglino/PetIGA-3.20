# fourteauImpactWaterVapor2021 — Fourteau, Domine & Hagenmuller (2021b), The Cryosphere 15, 2739–2755

**Summary.** Theory of coupled heat conduction and latent-heat transport by vapor in snow under a temperature gradient, in the slow- and fast-kinetics limits; steady-state FE simulations (Elmer) on 34 µCT snow microstructures; effective conductivity K_eff = K_cond + latent contribution; temperature dependence.

**Findings relevant to the manuscript**
- Abstract: "Numerical simulations on measured three-dimensional snow microstructures reveal that the effective thermal conductivity of snow can be significantly larger, by up to about 50 % for low-density snow, than if water vapor transport is neglected" (fast kinetics); "the fast kinetics hypothesis could be a reasonable assumption for modeling heat and mass transport in snow".
- Slow kinetics: "Kslow = 0" latent contribution (Sect. 2) — the conductive-only result is the slow-kinetics limit.
- Intro: "Calonne et al. (2011) 'recommend purely conductive effects ... to be considered separately from non-conductive processes'"; conductivity "plays an integral role in the energy budgets of the ground ..., ice caps and glaciers ... and sea ice"; "variations of thermal conductivity between snow layers impact the temperature gradients at the layer scale and thus in part govern snow metamorphism (Vionnet et al., 2012)".
- Temperature: with slow kinetics K_eff "only barely decreases over the range of temperature studied, consistent with ... Calonne et al. (2011)"; with fast kinetics low-density snow K_eff increases strongly with T (Fig. 5).
- Microscale treated as steady state (Sect. 2): "the timescale governing the microscopic scale is much shorter than the macroscopic timescale".

**Claims supported / qualified**
- Outline L565–568 ("vapor transport under a gradient adds to the apparent conductivity; the conductivity here is conductive only, which is exact for isothermal sintering"): **partially**. Supported: latent heat adds up to ~50 % (low density, fast kinetics) and is zero in the slow-kinetics limit. Qualification: "exact for isothermal sintering" conflates the state of the sample with the property: any heat flux through the deposit requires a gradient, and then the latent contribution appears whenever kinetics are fast. Better: "the conductivity computed here is the conductive part (K_cond); with α_c = 10⁻³ the latent contribution is expected to be small (slow-kinetics limit), but with fast kinetics it can add up to ~50 % in low-density snow (Fourteau et al. 2021b)". Also note K_eff here is 2-D.
- Outline §2.3 "the cell problem is a steady problem on a frozen microstructure": consistent with Fourteau's steady microscale treatment.
- Outline Intro ¶2 add (conductivity sets heat flow and temperature profile; input to snowpack models): Fourteau intro is a closely matching citation (energy budgets; layer-scale gradients govern metamorphism, citing Vionnet 2012).
- Outline §5.6 Limitations (Calonne 2019 on temperature-dependent conductivities): Fourteau Fig. 3/5 shows T-dependence of ki, ka and of K_eff — supports listing "constant k_ice, k_air" as a limitation; in the slow-kinetics case the effect is small.

**Candidate additions.** Rephrase L565–568 as above; cite Fourteau (2021b) in §5.6 for the T-dependence of K_eff through k_ice(T) being small in the slow-kinetics case.
