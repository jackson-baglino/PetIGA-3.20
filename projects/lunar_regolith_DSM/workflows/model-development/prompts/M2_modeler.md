# Work order: M2:Selection (beth-modeler)

Run folder: {run}
Read `{run}/derivation/spec_problem.md` (approved by the user at G:Scope) and `{run}/decisions.json`.
{feedback}
Write `{run}/derivation/selection.md` and add any new citations to `{run}/derivation/citations.json`.

`selection.md`: one `## Candidate: <name>` section for each of at least {min_candidates} model types{candidates_hint}.
In each: the equations in outline, how the grains enter, how the contact line and the no-flux conditions
are imposed, what can go wrong (with sources), the cost in DOFs for a reference case, how it fits this
solver, and how it would be verified. Then `## Comparison` (a table of the criteria: fidelity, geometry
flexibility, cost, fit to the solver, buildable and verifiable by December), then `## Recommendation` with
your choice, the reason, and a proposed short model name (letters, digits and hyphens) for the branch
`feat/<name>`.

Done criterion: at least {min_candidates} candidate sections, the comparison table and the recommendation
exist; every quote passes the quote check. The user makes the choice at G:Selection; you recommend.
