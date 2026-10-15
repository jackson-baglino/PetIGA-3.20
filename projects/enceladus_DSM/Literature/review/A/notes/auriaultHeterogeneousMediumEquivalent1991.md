# auriaultHeterogeneousMediumEquivalent1991 — Auriault 1991, Int. J. Engng Sci. 29(7):785–795

**Summary.** Methodological paper: when does a finely heterogeneous medium under an excitation admit an intrinsic macroscopic description? Answer: when scales separate (l/L = ε ≪ 1) for both geometry and excitation; double-scale asymptotic expansion with periodic (or stationary) cell problems is the self-consistent test. Two worked examples only: Stokes flow → Darcy law (Part 3) and elastodynamics of composites (Part 4). **No heat-conduction / thermal-conductivity example.**

**Findings relevant to the manuscript**
- Scale separation: p.2 (786) "the macroscopic domain must contain a (very) large number of heterogeneities ... there exists an elementary representative volume (ERV)"; "l/L = ε ≪ 1". Uses "elementary representative volume (ERV)" — draft writes "representative element volume (REV)".
- Periodicity as convenience: p.4 (788) "If we limit ourselves to the discovery of the structure of the macroscopic description, periodic and random media are equivalent, if homogenization is possible. It is therefore convenient to assume that the medium is periodic." Also p.5, p.7: local stationariness of a random medium "would lead to identical results".
- Cell problem / compatibility condition: p.4 (788) eqs. of form div_y Φ^(i) + div_x Φ^(i-1) = 0; zero-average source condition gives the macroscopic law.
- Effective coefficient as volume average of a cell solution: Darcy K = ⟨k⟩ (eq. 17, p.6) and c = ⟨a e_y(ξ + I)⟩ (p.10) — structurally the same as K_eff = ⟨K(∇t + I)⟩, but the paper never writes the conduction cell problem.
- No vector named **t**; elastic corrector is "ξ" (third-order tensor).

**Claims it supports / contradicts / qualifies**
- Draft L282 "homogenization technique initially outlined by Auriault (1991)": partially — the asymptotic two-scale method is outlined here (but originates in Bensoussan et al. 1978 and Sanchez-Palencia 1980, refs [1,2], p.1/p.11); the conduction problem itself is not in this paper.
- Draft L353 "We further assume, without less of generality, that the medium is periodic ... as described in Calonne 2011 and Auriault 1991": Auriault actually supports that periodicity is a convenience for obtaining the *structure* of the macroscopic law ("periodic and random media are equivalent, if homogenization is possible", p.4). This qualifies the outline's \tC that "periodicity is an assumption, not a convenience": for the structure of the law it is a convenience; for the *value* of K_eff of a finite random sample it is an approximation (cell must be an REV). Suggest wording: "We assume the medium is periodic, which leaves the form of the macroscopic law unchanged (Auriault, 1991) provided the cell is a representative volume."
- Draft L362 "As is shown in Auriault (1991) and Calonne (2011), the corresponding macroscopic description is ∇·(K_eff∇T)=0" and the K_eff formula: partially — Auriault 1991 shows the method on Darcy/elastic examples only; heat conduction is shown in Calonne 2011 (and the 2009 book). Better cite Auriault, Boutin & Geindreau 2009 (book) for the conduction cell problem.

**Candidate citations / claims**
- Bensoussan, Lions & Papanicolaou 1978 / Sanchez-Palencia 1980 are the origins of periodic homogenization (refs [1,2]); not in library, only mention if wanted.
- Prefer the 2009 book (Auriault, Boutin & Geindreau) for the heat-conduction cell problem and the name of the vector t (check in book notes).
