# Work order: M4:Discretization checks (beth-verifier)

Run folder: {run}
Read `{run}/derivation/spec.md`, `code_map.json`, `mms/` and the ledger. Do not read docs/planning.
Write `{run}/checks/m4_*.py` (with `--spec`/`--ledger` arguments as in M3). Test: integrating the weak form
by parts gives back the strong form of the ledger and its boundary conditions; each Jacobian entry equals
the derivative of its residual entry, including the time-derivative (shift) terms; the manufactured-solution
sources satisfy the strong form; the dimensionless groups are consistent with M1.
Run them with `{py} {run}/checks/<file>.py`.

Done criterion: every residual and Jacobian entry is covered by a check that ran.
