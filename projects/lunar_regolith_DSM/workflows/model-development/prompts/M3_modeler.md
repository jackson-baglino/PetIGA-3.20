# Work order: M3:Derivation (beth-modeler)

Run folder: {run}
Read `{run}/derivation/spec_problem.md`, `{run}/derivation/selection.md` and `{run}/decisions.json` (the
user's choices so far, including the model type chosen at G:Selection). Read `src/assembly.c`
(`Residual_A1`) and `src/material_properties.c` in {root}.
{feedback}
Write `{run}/derivation/ledger.json` and `{run}/derivation/derivation.md`.

`ledger.json`: {{"symbols": [{{"symbol", "meaning", "units"}}], "steps": [ ... ]}}. Each step:
{{"id": "S1", "type": "cited"|"derived"|"assumption"|"choice", "statement": "<LaTeX>", ...}} with
- cited: "paper_id", "page_pdf", "quote", "equation" (printed number or "");
- derived: "from": [earlier step ids], "operation": "<named operation>";
- assumption: "label": "A1", "if_wrong": "<what changes>";
- choice: "label": "C1", "options": [...], "recommendation": "...".
Step ids are unique; a derived step cites only earlier steps.

The derivation covers, in order: the strong form of every equation; the reduction to today's equations
where the grains are absent, term by term; the nondimensional form with each group; the boundary and wall
conditions (start from today's validated wall term); whether the thin-interface correction applies; the
sharp-interface limit to first order in the interface widths. `derivation.md` presents the ledger in order,
each displayed equation tagged with its step id.

Done criterion: every step is typed and sourced; every cited quote passes the quote check; every symbol is
in "symbols" before its first use; `derivation.md` mentions every step id. Stop with `blocked` at any choice
that is the user's.
