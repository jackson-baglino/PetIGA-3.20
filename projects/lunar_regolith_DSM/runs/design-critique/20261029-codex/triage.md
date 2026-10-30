# Triage of the Codex review (2026-10-29)

Codex (gpt-6.1-sol, clean CODEX_HOME, read-only, no web searches in its log)
reviewed revision 9 plus STATUS.md. Findings: `findings.md`. Recommendations
below are the main session's; the user decides.

| # | Finding | Checked | Recommendation |
|---|---|---|---|
| 1 | Build a small B1 path with a stopping point; defer LaTeX/JATS, whole-paper tables and figures, embeddings, desk, discovery | fits the model-first rule | **Accept.** Becomes revision 10. Keep P1 as built. |
| 2 | Blind transcription cannot check prose claims (conditions, attribution) | yes, X1 outputs are LaTeX or value+unit only | **Accept.** X1 for equations and values; prose claims get a semantic check and user approval when they support a model option. |
| 3 | Checker only sees condition pages the reader chose | yes (E4 inputs) | **Accept.** Checker writes an applicability checklist from the methods section before seeing claims; test it with the user's dropped-condition plants. |
| 4 | Two table extractions can agree on a wrong cell | design point | **Accept.** Defer automated tables; check memo-used cells against the tile image. |
| 5 | LaTeX canonicalizer may erase bold/vector distinctions | design point, not built | **Accept.** Keep font commands; paired must-match / must-differ tests when built. |
| 6 | Dimensional checks cannot catch signs, coefficients, scaling | correct | **Accept.** Defer the general operator layer; a user-checked symbol and scaling table for the equations that enter the model. |
| 7 | Calibration mix tests transcription more than judgment | partly | **Accept.** Use the existing 40 statements and 68 plants as a diagnostic; write no more plants for B1. |
| 8 | Per-kind gate impossible at 8–18 items | **verified**: n of n gives 0.05^(1/n); 8/8 = 68.8 %, 18/18 = 84.7 %, 19/19 = 85.4 % | **Accept.** Remove the 8-item rule; report per kind with intervals. |
| 9 | Fresh held-out sets after failures will run out; plants from one statement are correlated | correct | **Accept.** One frozen final evaluation for B1; report by statement and paper too. |
| 10 | Too many roles (E2, X1, two memo checkers, revisions) | design point | **Accept.** One reader, fresh checker, no routine E2, one skeptical memo checker. |
| 11 | Counts-only main session makes evidence review awkward | real trade-off | **User decides.** Option: a dedicated review session that may read selected evidence but has no network, writes or broad file access; defer the desk. |
| 12 | Security tests check exfiltration, not a steered false claim | correct | **Accept.** One injected-instruction fixture whose scientific claim must stay right. |
| 13 | No single B1 roster or run contract | partly: STATUS.md said "9 seed papers"; B1 has 7 seeds (§9 item 38) plus Kaempfer & Plapp and Moure & Fu = 9 papers | **Accept.** Write the B1 run contract (papers, hashes, questions, outputs, stage exceptions). |
| 14 | Immutable P1 content vs later symbols/units; memo verification order | correct | **Accept.** Analysis records stored separately; claims → draft memo → verify cited claims → approved memo. |
| 15 | No restore test; no public fixture | correct | **Accept in part.** One restore test; defer the content-addressed archive. Public fixture depends on the course requirement (user). |

Questions for the user (from Codex, still open):
- Which one or two modeling decisions must the B1 memo support?
- How many hours can go to pipeline building and evidence review before the conference?
- Does any B1 paper need a scan, supplement or correction?
- What reproducibility evidence does the course need?
