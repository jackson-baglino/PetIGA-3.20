# Work order: M3:Derivation checks (beth-verifier)

Run folder: {run}
Read `{run}/derivation/ledger.json` (the claims). You may read `{run}/derivation/derivation.md` for the
statements but test each claim independently. Read `src/assembly.c` on main ({root}/src) for the reduction
check. Do not read docs/planning.
Write `{run}/checks/m3_*.py`. Every script takes `--ledger <path>` (default `{run}/derivation/ledger.json`)
so the line can run it against a ledger with planted errors; it must read the statements it tests from the
ledger, not copy them into the script.
Test every derived step (the named operation really gives it), the reduction to today's equations,
conservation of total water in a closed box, free-energy decay without the source, units of every term, and
flag any step that needs an assumption not in the ledger.
Run them with `{py} {run}/checks/<file>.py`.

Done criterion: every derived step is covered by a check that ran. Report failures as findings.
