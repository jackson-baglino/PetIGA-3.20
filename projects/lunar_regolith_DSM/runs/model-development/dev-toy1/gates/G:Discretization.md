# Gate G:Discretization (M4), model `heat1d-fd-be`: DEV TOY DRY RUN

This is a plumbing test only. It is not a result and must not be cited (spec.md, line 3).

## 1. What M4 produced
- `derivation/spec.md`: weak form (M4-1, M4-2), FD rows (M4-3), backward Euler (M4-4), direct steady rows (M4-5), residual and term table T1–T5 (M4-6, M4-6b), Jacobian (M4-7), dt rule (M4-8), manufactured solution (M4-9), the layers and proposed tests (infra-1..4, heat-1..8), and a read-only report on the anchor solver's TSALPHA/TSAdapt.
- `derivation/code_map.json`: 20 entries. They map T1–T5, their Jacobian rows, the Jacobian check, MMS, the scales, the BE step, the Thomas solve, the steady solve, the end fluxes and the dt helper to functions in `toy/heat1d.py`.
- `derivation/mms/mms_sources.py`: a SymPy script that prints f, checks the data and corners, and prints samples.

## 2. What the checks found (`checks/results_M4.json`)
- 4 check scripts plus 1 shared parser. All exited 0, with empty stderr. **32 of 32 claims pass, and there are 0 failures.**
  - `m4_jacobian`: 6 of 6 pass. The full-J relative error is ≤ 1.1e-13, and ∂F/∂U̇ alone is ≤ 3.5e-11, against a tolerance of 1e-8. When the check was run on a deliberately broken Jacobian (a planted error), it gave errors of 1.0e-3 with the shift term dropped and 0.27 with an off-diagonal sign flipped, so it catches both errors.
  - `m4_mms`: 8 of 8 pass. The observed orders are 2.008 in space and 1.007 in time. The reference error is 4.3056e-5.
  - `m4_scales`: 7 of 7 pass. At τ = 1 the discrete floor is 7.0172e-6 in θ, which leaves a margin of only 1.43× below 1e-5. The half-cell end flux converges at order 2.0.
  - `m4_weak_strong`: 11 of 11 pass.
- Two sets of numbers differ, but neither is a failure. The modeler's own scratch Jacobian check reports larger errors than the verifier's: 7e-11 to 1e-10 overall, ∂F/∂U̇ up to 6.3e-8, and a drop-shift error of 9.6e-2. The two checks are separate scripts at different states. Both are far inside their tolerance, and both catch a broken Jacobian.
- The verifier says it stopped an earlier, slow version of `m4_weak_strong.py` and did not restart it. The results above come from the replacement version.

## 3. Paper claims to confirm
None. `derivation/citations.json` is empty, and the spec cites no paper.

## 4. Decisions for you
1. **Where should the new discretization symbols (N, h, Δτ, U, shift, J, f, …) be recorded?**
   - (a) Keep them only in the spec.md symbol table.
   - (b) Add them to ledger.json when the ledger is next opened.

   The modeler gave no recommendation. M4 is not allowed to write to the ledger.
2. **What should `solve` do when t_end/dt is not an integer to within 1e-10 relative?**
   - (a) Raise `ValueError`.
   - (b) Shrink dt silently so the steps fit.

   The modeler recommends (a), because the M5 time-order test needs the step to be exactly the dt the caller passed.
3. **Do you accept the proposed solver settings?** These are: fixed-step BE; one Newton step per time step, followed by a residual guard; the Thomas solve; float64; and the same solve with shift = 0 for the steady case.
   - (a) Accept them.
   - (b) Change them.

   The modeler recommends (a). R is affine in U, so one Newton step is exact up to roundoff, and the code then calls the same `residual` and `jacobian` that the Jacobian test checks. The constant c in the guard is left for M5 or the verifier to set; the modeler would start at 10³.
4. **Reminder, already decided (D-d):** the heat-4 tolerance on τ_end against S44 must be looser than 1e-5 in θ. M5 sets the value. The floor is 7.0e-6.

Earlier gates were approved by the main session as a dry run, not by you (decisions.json).

## 5. What I could not verify from the files
- `toy/heat1d.py` does not exist yet, so none of the code_map locations can be checked.
- `mms/mms_sources.py` was never run as itself. The modeler ran a copy, and the verifier's MM6 compared its expressions with the spec.
- The dt rule M4-8 is a first-order estimate for mode 1, checked numerically to within 3 %. It is not proved to be a bound.
- No check confirms the file and line citations in the TSALPHA section (alpha1.c, lunar_main.c, monitoring.c, tsadapt.c). J5 checks only the algebra. The patched header `tsimpl.h` was not read, so the declarations and defaults of `adap`, `NRmin`, `NRmax`, `factor`, `dtmin`, `dtmax` and `maxdt_dyn` are unknown.
