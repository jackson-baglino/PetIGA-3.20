# Work order: M8:Performance (beth-modeler)

Run folder: {run}
Feature worktree: {worktree}
{feedback}
Measure, one run at a time through run_test.sh, the solver settings, the dt rule and the parallel sizing
listed in `## Solver settings` and `## Time scales and time step` of the spec, on the two cases the work
order names: {perf_cases}. Every change is behind an option; every locked test must still pass after the
change you keep. Write `{run}/perf.csv` (case, setting, ranks, dofs_per_rank, wall_s, newton_its, krylov_its,
steps, result) from the run logs only, and `{run}/perf.md` with the choices and their measured reasons.
Done criterion: each choice has a measured reason; the locked tests pass with the chosen settings.
