# Work order: M4:Discretization (beth-modeler)

Run folder: {run}
Read the ledger and `derivation.md` in `{run}/derivation/`, `{run}/decisions.json`, and in {root}:
`src/assembly.c`, `src/lunar_main.c` (time stepping), `src/monitoring.c` (time-step control),
`inputs/solver.opts`. Check in the PETSc source ({petsc_dir}/src/ts/impls/implicit/alpha/) how TSALPHA
interacts with TSAdapt, and say what you found with file and line.
{feedback}
Write `{run}/derivation/spec.md`, `{run}/derivation/code_map.json` and `{run}/derivation/mms/` (manufactured
solutions: a script that prints the source terms).

`spec.md` sections, in order: `## Weak form`, `## Discretization`, `## Residual`, `## Jacobian`,
`## Auxiliary functions` (indicator and gradient, initial conditions, extension into grains, output masks),
`## Time scales and time step` (a priori scales with values; the dt rule), `## Solver settings` (proposed,
with reasons), `## Manufactured solutions`, `## Layers` (the implementation order, each layer with the tests
that must pass before the next).
`code_map.json`: a list of {{"term", "row", "where" (file and function), "option", "layer"}}.

Done criterion: all sections exist; every residual term has its Jacobian entries; every term is in
code_map.json; the dt rule names its inputs and formula.
