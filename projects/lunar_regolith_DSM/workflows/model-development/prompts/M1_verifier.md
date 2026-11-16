# Work order: M1:Scope checks (beth-verifier)

Run folder: {run}
Read `{run}/derivation/spec_problem.md` (the claims to test). Do not read anything under docs/planning.
Write check scripts `{run}/checks/m1_*.py`, one per group of claims, as your instructions describe.
Test: every scaling estimate and dimensionless group recomputed from the stated inputs (Pint for units,
within 1 %); every unit in the notation table consistent with the formulas that use it; every constitutive
law's stated validity range consistent with the problem's temperature range as stated in the file.
Run them with `{py} {run}/checks/<file>.py`.

Done criterion: every numeric claim in `## Mechanisms` and `## Scales` has a check that ran; the scripts
print one JSON line per claim and exit 0 only when all pass. Report failures as findings in `evidence`.
