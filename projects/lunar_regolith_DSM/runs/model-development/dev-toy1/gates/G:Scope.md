# Gate G:Scope (M1): toy 1D heat conduction (dry run only)

## 1. What the stage produced
- `derivation/spec_problem.md`: the problem, mechanisms (3 in, 3 out), scales, constitutive laws, notation, out-of-scope list, and two choices (C1, C2).
- `derivation/citations.json`: an empty list. No papers are cited.
- `derivation/scratch/m1_scales.py`: the modeler's scale computations.
- Verifier checks: `checks/m1_scales.py`, `m1_nondim.py`, `m1_units.py`, `m1_validity.py`, `m1_inputs.json`; results in `checks/results_M1.json`.

## 2. What the checks found
- All 4 scripts exited 0 with empty stderr. `all_pass: true`. **No failures.**
- Count: I count 58 claim entries with `pass: true` (nondim 15, scales 16, units 20, validity 7). The verifier's summary says 56. The difference does not change the outcome, but the two counts do not match.
- Points to note (these passed, but you should know about them):
  - **Steady-state time:** at τ = t_99 ≈ 0.4666, the whole-solution error is 1.03 % of the peak (0.00129 against 0.125), not 1 % as stated. The spec rounds τ_steady to 0.47; the computed value is 0.4666 (difference 0.73 %, within the 1 % check tolerance). When you set t_end for the tests, allow some margin above 0.47.
  - **Constant k:** β ΔT_max = 0.125·m for a power law k ∝ T^m. That is below 1 only for |m| < 8. No β is sourced; the brief accepts the assumption.
  - **T_0 = 1 K in the unit case:** the verifier marks continuum validity at 1 K as "not assessable: no material" (informational only). This does not affect the dimensionless problem.
  - **Swapped parameters:** every reference parameter is 1, so the numeric checks cannot catch swapped parameters. The symbolic and unit checks cover those formulas.

## 3. Claims from papers to confirm
None. `citations.json` is empty. Fourier's law and the caloric law are attributed to the brief (line 8) as standard laws.

## 4. Decisions for you
`decisions.json` is empty. The spec raises these two choices:
- **C1: Which reference parameter values should the line use?** Options: (a) the unit case L = k = ρ = c_p = q = T_0 = 1 (SI units), which gives Q = 1; (b) a physical material, such as a metal rod, whose values would need a source; (c) values you specify. **Modeler's recommendation: (a).** Reason: the equation in θ has no parameters, so the unit case tests the plumbing fully, gives a clean analytic peak of 0.125 K, and needs no external source. All the numbers in the spec assume (a).
- **C2: Should the code allow negative q (a heat sink)?** Options: allow any real q, or restrict to q ≥ 0. **Modeler's recommendation:** allow any real q in the code and test only with q > 0. Reason: the parabola test holds for any sign of q; only the stated range of T assumes q ≥ 0.
- The spec also asks you to approve the locked notation table ("pending user approval").

## 5. What I could not verify
- **Brief line citations:** the brief `workflows/model-development/problems/toy_heat1d.md` is outside my read roots. I could not check any of the line citations to the brief (lines 3–4, 8–9, 13–17, 21–23, 26).
- **Scratch script output:** I did not check that the output of `derivation/scratch/m1_scales.py` matches the spec. The spec says it does. The verifier's own `checks/m1_scales.py` reproduces the values.
