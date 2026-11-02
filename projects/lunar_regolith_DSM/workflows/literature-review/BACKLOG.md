# BACKLOG: literature-review

Ideas for improving the workflow that came up while doing planned work
(`.claude/rules/workflows.md`). They are read and decided at the start of the
next develop phase, not acted on when written. One entry each: date, source
(session, stage or run), the idea, and why.

## Open
- 2026-11-01 · dev-20261101-kp, foreman report · run_b1.py status lists stages in flowchart node order (A2 after U4), which reads as if A2 ran late; list them in dependency order instead. Why: The foreman took the listing for a stage-order problem; the line itself follows the requires links and ran A2 before P1.
- 2026-11-01 · dev-20261101-kp · prepare.py refuses to write content that already exists, so a P1 interrupted mid-paper cannot resume; write to a staging folder and rename when complete. Why: Stopping the line during P1 would otherwise need a manual cleanup before it can run again.
- 2026-11-01 · beth-foreman run dev-20261101-kp 2026-11-01 · run_b1.py status should name the stage and paper behind each queued review item, and print the what text of a waiting gate Why: The foreman must tell the user what to decide, but status shows only counts, so the foreman cannot say which check raised the item without opening run files it may not read

## Decided
- 2026-11-01 · DONE 2026-11-01 (U3a requires X1; edge X1→U3a; test in test_run_b1.py) · dev-20261101-kp · Make U3a require X1 in flowchart.json, so X1's disagreements are in the user's U3 pass; in node order the U3a gate opened before X1 ran. Why: u3_review.py reviews X1 disagreements together with the model-option claims; with U3a first, the user would review twice.
- 2026-11-01 · DONE 2026-11-02 (run_b1.py holds the gated item, continues with the rest; test) · dev-20261101-a2, A2 via run_b1.py · A waiting per-paper gate (A2 li_2009) stops the whole line, so the other papers' A2 checks never run; let a per-paper blocking gate hold only that paper's downstream stages while the line continues with the other papers Why: In the nine-paper run two A2 gates (Li, Lervag) would make the user approve them one at a time before the other seven copies are even checked
- 2026-11-01 · DONE 2026-11-02 (ledger is the run's calls.jsonl; summary_x1_reads counts calls and reused reads; a read of an older memo version is redone; test) · dev-20261101-kp, X1 via run_b1.py (8.5 s, no call in the run ledger) · blind_read.py: log calls to the run's calls.jsonl (not untrusted/x1/calls.jsonl), and print reused vs new reads in summary_x1 (it skips a paper whose blind read already exists, so a hand test can stand in for the line's read) Why: The session cannot tell from counts whether X1's 6/6 agreement came from a read made by the line or from an earlier hand test; every check in a frozen run must trace to a recorded call
