# BACKLOG: literature-review

Ideas for improving the workflow that came up while doing planned work
(`.claude/rules/workflows.md`). They are read and decided at the start of the
next develop phase, not acted on when written. One entry each: date, source
(session, stage or run), the idea, and why.

## Open
- 2026-11-01 · dev-20261101-kp, foreman report · run_b1.py status lists stages in flowchart node order (A2 after U4), which reads as if A2 ran late; list them in dependency order instead. Why: The foreman took the listing for a stage-order problem; the line itself follows the requires links and ran A2 before P1.
- 2026-11-01 · dev-20261101-kp · prepare.py refuses to write content that already exists, so a P1 interrupted mid-paper cannot resume; write to a staging folder and rename when complete. Why: Stopping the line during P1 would otherwise need a manual cleanup before it can run again.
- 2026-11-01 · beth-foreman run dev-20261101-kp 2026-11-01 · run_b1.py status should name the stage and paper behind each queued review item, and print the what text of a waiting gate Why: The foreman must tell the user what to decide, but status shows only counts, so the foreman cannot say which check raised the item without opening run files it may not read
- 2026-11-01 · dev-20261101-kp · Make U3a require X1 in flowchart.json, so X1's disagreements are in the user's U3 pass; in node order the U3a gate opened before X1 ran. Why: u3_review.py reviews X1 disagreements together with the model-option claims; with U3a first, the user would review twice.

## Decided
