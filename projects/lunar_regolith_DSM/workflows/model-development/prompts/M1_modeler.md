# Work order: M1:Scope (beth-modeler)

Run folder: {run}
Problem brief: {problem}
User decisions so far: {run}/decisions.json (may be empty)
{feedback}
Write `{run}/derivation/spec_problem.md` and `{run}/derivation/citations.json`.

`spec_problem.md` must have these sections, as level-2 headings, in this order:
1. `## Problem`: the physical system and the questions the model must answer.
2. `## Mechanisms`: each mechanism considered, in or out, with a scaling estimate (an order of magnitude with
   the numbers and their sources) that justifies the decision.
3. `## Scales`: the characteristic scales (length, time, temperature, vapor density) and every dimensionless
   group with its value for a reference case; say which terms dominate.
4. `## Constitutive laws`: each law, its formula, its source, its validity range, and whether it holds at
   the temperatures and pressures of the problem.
5. `## Notation`: every symbol, its meaning and units; the sign of the normal; the range of each field.
6. `## Out of scope`: what is left out and why.
7. `## Choices for the user`: each decision that is the user's (C1, C2, ...), with the options, your
   recommendation and the reason.

`citations.json` is a list of objects {{"id", "paper_id", "page_pdf", "quote", "used_for"}}, one per cited
passage; each quote verbatim and checked with
`{py} {root}/scripts/model_development/library_search.py quote "<quote>" --paper <id> --page <n>`
(search with `... library_search.py search "<terms>"`). An empty list is allowed when nothing is cited.

Done criterion: both files exist; the seven sections exist in order; every quote passes the quote check;
every number in `## Scales` and `## Mechanisms` names its source (a citation id, a file and line, or a
computation shown in the text). Put your own scratch scripts in `{run}/derivation/scratch/`.
