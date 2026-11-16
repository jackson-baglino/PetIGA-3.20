# Work order: M9:Write-up (beth-writer)

Write `{doc}`: a teaching document that explains the model, its derivation and its implementation to the
user, a PhD student who will build on it. Sources (read them; do not write from memory):
`{run}/derivation/` (spec_problem.md, selection.md, ledger.json, derivation.md, spec.md, code_map.json),
`{run}/decisions.json`, `{run}/checks/` results, `{run}/tests/` results, `{run}/mutation.json`,
`{run}/perf.md`.
{feedback}
Sections, in this order, each using only what came before: Problem; Unknowns and assumptions; Strong form
(term by term, meaning and source); Nondimensional form; Sharp-interface limit; Weak form; Discretization;
Residual and Jacobian; Time step; Code map; Tests and what each proves; Performance choices; Limitations.
Start with a `## Symbols` table (symbol, meaning, units, where defined). Every equation carries the ledger
step id; every number names its file. Plain American English, concise and linear.
Done criterion: the file exists with every section; every symbol in the table is defined before its first
use in the text.
