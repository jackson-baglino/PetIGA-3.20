# Gate G:Selection (M2): toy 1D heat conduction (dry run only)

## 1. What the stage produced
- `derivation/selection.md`: three candidates and a recommendation.
  - **FD-BE:** finite differences with backward Euler.
  - **FEM-Galerkin:** P1 elements or quadratic B-splines, with backward Euler.
  - **Series:** the closed-form eigenfunction series.
  - The file also has the modeler's derivations D1–D9, a comparison table, the proposed verification tests V1–V4 and a proposed reference discretization.
- `derivation/citations.json`: an empty list. No papers are cited.
- Check results: `checks/results_M2.json`.

## 2. What the checks found
- `all_pass: true`. The only check recorded is `"citations": "OK citations"`. **No failures.**
- **Coverage:** no verifier check of the derivations (D1–D9) or the numbers is recorded for M2. The modeler also says it could not run scratch Python (a guard blocked it), so every number is a hand computation.
- **What I checked by hand.** These agree with the file:
  - D5: (π/64)²/12 = 2.008·10⁻⁴.
  - D6: about 2.34·10⁻⁴ peak time error in mode 1 for Δτ = 10⁻³.
  - D9: b_k = 4/(kπ)³ for odd k, and the sum at ξ = 1/2 equals 1/8.
  - The end-flux error is h/2, and the half-cell fix gives 1/2.
  - The truncation tail for K = 11 is about 2.7·10⁻⁴.
  - D8: the off-diagonal sign condition is Δτ ≥ h²/6.
  - Every spec line cited in the file (L3–5, 42, 46–47, 50–52, 98–100, 108–111, 129–130) says what the file claims.
- **A point that passed but affects V1's tolerance.** V1 calls for agreement "to roundoff level," but at τ_end = 1 a residual transient remains.
  - The modeler estimates it with the exact decay: 0.129·e^(−λ₁) ≈ 6.7·10⁻⁶.
  - With the backward-Euler factor (1 + λ₁Δτ)^(−1000), my estimate is about 7.0·10⁻⁶.
  - Either way, V1 at τ_end = 1 cannot pass at roundoff unless a direct steady solve is added.

## 3. Claims from papers to confirm
None. `citations.json` is empty, and every claim comes from the spec or is derived in the file.

## 4. Decisions for you
- **Which discretization should the toy use?**
  - Options: (a) FD-BE, with Series as the reference solution in the transient test; (b) FEM-Galerkin, which reopens the brief's choice of FD with backward Euler.
  - **Modeler's recommendation: (a).**
  - Reasons:
    - The brief already fixes FD with backward Euler in NumPy (spec L129–130).
    - FD-BE is stable for every Δτ and keeps 0 ≤ θ ≤ 1/8 (D3, D4).
    - Lumped P1 FEM gives the same scheme in 1D (D8b).
    - Comparing with the series (V2) covers the spatial order, which the steady test cannot catch.
  - Proposed model name: `heat1d-fd-be`, on the branch `feat/heat1d-fd-be`.
- **Do you approve the reference discretization N = 64, Δτ = 10⁻³, τ_end = 1?**
  - Options: approve, or give other values.
  - **Modeler's recommendation: approve.** Reason: τ_end = 1 is well above τ₉₉ ≈ 0.47, which is the margin G:Scope asked for. The cost is negligible: 63 unknowns over 1000 steps.
- **How should V1 (the steady-state test) handle the residual transient?**
  - Options: (a) compare at τ_end = 1 with a tolerance above about 7·10⁻⁶; (b) add a direct steady solve and compare at roundoff level; (c) both.
  - **No recommendation from the modeler.** The file raises the option of a direct solve but does not choose.
- **Do you approve the G:Scope decisions C1 and C2?** In `decisions.json`, the main session approved them "for plumbing only (not a user decision)":
  - C1 = (a), the unit case;
  - C2 = any real q in the code, tested with q > 0.
  - This stage relies on both.

## 5. What I could not verify
- **Brief line citations:** the brief `workflows/model-development/problems/toy_heat1d.md` is outside my read roots, so I could not check the citation "brief lines 14–16."
- **Unchecked derivation claims:** I did not recheck these by hand:
  - D2 (eigenvalues of D_h);
  - D4 (M-matrix bounds);
  - D7 and D7b (FEM nodal exactness);
  - the Crank–Nicolson sign-flip claim.
  - No verifier check covers them either.
- **Unshown claim (D8):** D8 shows only that a sufficient condition is lost for Δτ < h²/6. Whether negative values actually occur with a consistent mass matrix is not shown. This matters only if you choose FEM-Galerkin.
