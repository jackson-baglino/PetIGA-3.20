# Gate G:Test lock (M5), model `heat1d-fd-be`: DEV TOY DRY RUN

This is a plumbing test only. It is not a result and must not be cited. Earlier gates were approved by the main session, not by you (decisions.json).

## 1. What M5 produced
- `studies/heat1d-fd-be/tests/tests.json` (in the worktree): 18 tests. Tests infra-1..4 and heat-1..8 come from the spec's Layers. Tests heat-9..14 are the verifier's additions: transient scaling P vs C1, the dt rule M4-8, the ValueError rule, a dense M4-4 oracle, energy decay, and the discrete heat balance. All use `"baseline": "none"`, runner python, 1 rank, and 2 to 10 min wall time.
- Scripts `test_infra_1..4.py` and `test_heat_1..14.py`, plus `INTERFACE.md` and `reference/README.md` (listed in `gates/M5_verifier.result.json`).
- `dev-toy1/tests/_verifier_ref/`: the verifier's own reference module `ok/` and six planted-error copies (mut_ftime, mut_kappa, mut_gR, mut_shift, mut_flux, mut_mmsQ). None of them has been run.
- `tests.json` also has a "not_applicable" list: the regression run on main, psi = 1 against today's model, neck width, and the closed-box water test.

## 2. What the checks found
- **`checks/results_M5.json` does not exist.** The only M5 evidence is the verifier's report and the draft run logs.
- Draft run on the worktree: 18 of 18 tests produced a result JSON and a log. **18 of 18 fail**, each with `ModuleNotFoundError: No module named 'toy'` (log line 2 in every log). This failure is expected before M6. It shows that each script compiles, imports from `--project`, and prints a verdict. It does not show anything more.
- **The test bodies have never run.** The guard blocked running the tests against `_verifier_ref/ok` and the mutants. So nobody has shown that any test passes on a correct implementation or fails on a planted error.
- Spec numbers the tests depend on were not reproduced independently: e_16 = 5.2340e-4 and e_32 = 1.3016e-4 (heat-1), e(1e-2) = 7.5869e-4 and e(5e-3) = 3.7746e-4 (heat-2), the floor 7.0172e-6 (heat-4), and the ratio of about 0.9993 at K = 1000 (heat-13).

## 3. Paper claims to confirm
None. The tests cite only spec items (M4-x, S-x, D-x).

## 4. Decisions for you
1. **Do you approve the call signatures in `INTERFACE.md`?** The spec names these functions but does not give their arguments:
   - `residual(tau, U, Udot, n, params)`
   - `jacobian_tridiag(shift, n, params) -> (sub, diag, sup)`
   - `solve_steady(n, params) -> (x, T)`
   - `end_fluxes(..., f_end=(f0, f1)) -> (F'_0, F'_1)`
   - `check_jacobian` returns a dict with `'pass'` and looks up `jacobian_tridiag` by name at call time (infra-4 depends on this).

   Options: (a) approve; (b) change them. M6 needs an answer before it starts. The verifier asks for approval.
2. **Should the line allow a calibration run before the tests are locked?** The run would be `run_test.sh ... --project <run>/tests/_verifier_ref/<name> --out <run>/tests/_verifier_ref/results --draft`. Options: (a) allow it, and hold the lock until each test passes on `ok` and fails on its mutant; (b) lock now without calibration. The verifier asks for (a). Its reason: the checking code has never executed, so a test that always fails or always passes would not show up until M6, where it would look like a code fault.
3. **Do you accept the tolerances?**
   - heat-1/2: |p − 2| ≤ 0.1 and |p − 1| ≤ 0.1, with magnitudes within 2 %.
   - heat-4: d ≤ 2e-5 (this meets D-d), plus the floor 7.0172e-6 within 1 %.
   - heat-3/7: 1e-11 in θ′.
   - heat-12: 1e-12.
   - heat-13: ratio ≥ 0.99.
   - heat-14: 1e-9·max(1, |Q|).
   - infra-3: 1e-8 relative, with an FD step of 1e-3.

   Options: (a) accept them all; (b) keep only the order and bound claims (drop the 2 % and 1 % magnitude bands) until calibration reproduces the spec numbers. The verifier gives no recommendation. The magnitude bands are the only checks that catch forcing evaluated at τⁿ and a wrong forcing coefficient (heat-12 also catches the first of these).
4. **Is it enough for heat-5 to check D4 and monotonicity at sampled step counts?** Options: (a) accept sampling; (b) add `step_be(U, tau_new, dtau, n, params)` to the interface so heat-5 can check every step. The verifier gives no recommendation.

## 5. What I could not verify from the files
- I did not read the test scripts, `INTERFACE.md`, or `reference/README.md`. The guard blocks reads of the worktree for this stage, apart from `tests.json`. The signatures above are as the verifier reports them.
- There is no M5 modeler result file, only `M5_verifier.result.json`.
- The infra-1 claim in `tests.json` says T3 is "−g/h²" in rows 1 and N−1. The verifier's report says "−1/h²". These agree only if g = 1.
- `state.json` still lists M6, M7 and M9 as "pending after reopened M3", even though M3 was approved at 18:46. The note looks stale.
- The verifier reports that it found no text in the files it read that looked like instructions. I found none in the files I read.
