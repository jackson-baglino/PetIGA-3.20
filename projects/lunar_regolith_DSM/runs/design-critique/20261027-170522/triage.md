# Triage: design critique 20261027-170522 (round 3)

Design: `workflows/literature-review/WORKFLOW.md`, revision 6 (commit
cd64e493). Critic: `beth-critic` (claude-opus-5-5), fresh, with the extra
questions in `runs/design-critique/focus/rev6_scripts_permissions.md`:
deterministic replacements, agent count, permissions. Result: 25 findings
(2 critical, 12 major, 11 minor) and 9 unclear points
(`critiques/beth-critic.json`). The transcript audit passed on the first
attempt; 3 findings cite a web source.

Checks made by the main session, 2026-10-27:
- **F5**, against the Codex docs (learn.chatgpt.com/docs/config-file/config-basic,
  redirected from developers.openai.com): "Codex enables web search by default
  for local chats and serves results from a web search cache"; it is turned
  off with `web_search = "disabled"`. **Confirmed.** The docs page does not
  say what the read-only sandbox can read; that is to be tested in T-codex.
- **F6**, against code.claude.com/docs/en/sandboxing: "By default, if the
  sandbox can't start … Claude Code runs commands without sandboxing", and
  `sandbox.failIfUnavailable: true` makes it exit instead. Setting
  `allowUnsandboxedCommands: false` makes Claude Code ignore
  `dangerouslyDisableSandbox`. Write paths are set per session in settings,
  not per command. **Confirmed.**
- **F13** (SymPy's LaTeX parser and ∇): not checked. It is plausible, and
  T-units will measure the coverage on the B1 seeds.

Recommendation key: **accept**; **modify** (accept with the change stated);
**user** (needs a decision only the user can make).

| # | Sev. | Finding (short) | Rec. | Change |
|---|---|---|---|---|
| F1 | crit | A checker's own correction is published unchecked | accept | an E4 correction is a new claim version and goes through E3, E3b, X1 (if important) and a fresh E4 |
| F2 | crit | Recorded facts and their inputs exist only on the local disk | user | a private, dated archive (tarball) of the P1 outputs, library and transcripts, with hashes in the manifest and checked by L2; never in git. Where it goes is the user's choice |
| F3 | major | Claims have no structured value or unit | accept | `parameter: {quantity, value, unit_printed, uncertainty?}` on claims and verdicts; X1, E3b and Y2 compare fields |
| F4 | major | LaTeX and XML text has no pages | accept | page membership is checked on the text of the imaged PDF (Marker or the text layer, used only to locate the quote); wording comes from the source; alignment failure is its own flag |
| F5 | major | Codex can read the whole disk and search the web | user (accept) | `web_search="disabled"` by `-c`; Codex runs isolated, either as a separate macOS user or in a container that mounts only the call folder (user picks); Codex's exec log is audited; T-codex plants both violations |
| F6 | major | Sandbox fails open; no per-script limits | accept | `failIfUnavailable: true`, `allowUnsandboxedCommands: false`; offline scripts run through a wrapper with a per-command profile (no network, writes to the run folder); every run records that the sandbox was active |
| F7 | major | The main session holds the lethal trifecta | modify | drivers give the main session counts, IDs and verdicts, not quotes, flags, scout reasons or passages; those go to the user in the reading view; flags that look like instructions go straight to the user; no WebFetch for the main session while a wave runs |
| F8 | major | Memos are checked by the same kind of model, and numbers without units | accept | Y2 compares value and unit through the F3 fields; for model-critical questions Y3 also runs on Codex, and disagreements go to U4 |
| F9 | major | Repeated U3b samples inflate the bound | accept | population = the wave's passed claims; after a failure the fresh 60 comes only from claims rechecked under the new prompt; the combined bound over all samples is reported |
| F10 | major | Prompts tuned and scored on the same plants; n≈4 per kind | user (modify) | split plants into tuning and held-out sets and score only the held-out set; pass on the overall lower bound and on the mechanical kinds; other kinds reported with their intervals, with a target of 8 per kind over time from the growing set. Costs more plants (mostly script and Codex) |
| F11 | major | Importance depends on the reader's chosen kind | accept | a script also marks as important any claim with a number and unit, a table cell or LaTeX; importance is recomputed after E4 |
| F12 | major | Title matching merges Part I and Part II | accept | years within ±1 and no series token differing; title-only merges go to U1 |
| F13 | major | SymPy cannot parse most phase-field equations | accept | an operator layer maps ∇, ∂ and · to their dimensions before SymPy; *not parsed* sends the equation to X1 and E4 and is not a units pass; coverage measured in T-units |
| F14 | major | Challenge sends unpublished ideas to web search | accept | a script turns the statement into keyword queries, the user approves them, and the scout sees only the terms and results |
| F15 | minor | Hashing a rebuilt SQLite file is fragile | accept | hash a canonical dump (sorted rows, embeddings rounded); one thread; a version change is a planned rebuild |
| F16 | minor | The reader gets every tile up front | accept | E1 gets the text, page images and disagreements; it reads tiles on demand; the audit logs which |
| F17 | minor | E4 resends the key pages for every page | accept | E4 runs per group of 3 to 5 pages, key pages loaded once, verdicts per claim |
| F18 | minor | R1 sends every claim in the library | accept | a script prefilters by tags, aliases and search score; the reader gets the top N plus claims on the pages found |
| F19 | minor | SC3 adjudication can be a rule | modify | rule: any include → include; both exclude on the same criterion → exclude; SC3 removed. SC0 bands skipping the screeners only after its agreement is measured over two waves |
| F20 | minor | Tagger and librarian do not justify separate agents | modify | T2 becomes rules in T1 (merge proposals above a similarity threshold), decided at U4; the reference desk becomes a skill the main session uses (`library_query.py` plus `check_citations.py`; Lab 5's skill), returning hits verbatim with an optional summary labeled as one; query logging recorded as an allowed write |
| F21 | minor | The 30-claim recheck measures repeatability, not correctness | accept | relabeled as repeatability |
| F22 | minor | A/B merge can join two different claims | accept | merge also requires similar claim text or matching parameter fields; variants are never shown to Y1 |
| F23 | minor | Marker and Surya fetch weights themselves | accept | all weights in a pinned local cache, offline variables set, hashes computed from loaded files and compared with a lock file |
| F24 | minor | An empty web search is weak evidence of a gap | accept | gaps also run D1's API queries and a P2 search; memos word a gap as "not found by queries …" with the queries listed |
| F25 | minor | The practice-library baseline measures disagreement | user (accept) | the user rules on a random 30 of the rechecked practice claims; error rates of both libraries reported against the user's verdicts |

Unclear points, answered in revision 7:
1. **What E4 sees.** E4 sees the claim's fields (it must check them), never
   the reader's reasoning or Marker text. Anchoring on the reader's values is
   why X1 is blind.
2. **U3b population.** The wave's passed claims (F9).
3. **SC0.** It uses the index's embedding model on a separate store of
   candidate abstracts.
4. ***Not parsed*.** Not a units pass (F13).
5. **X1 tiles.** A script picks the tile from the claim's page and position;
   if the quantity is near a tile edge, both tiles are given; *not found*
   goes to E4, not counted as a disagreement.
6. **Tier S.** U1 is still shown (short), A2 runs, E2 does not; the seed rule
   does not apply.
7. **Nothing goes to OpenAI until decision 21 is resolved.** C0 core starts
   with the user's and script plants; the Codex parts wait.
8. **Main session.** Its script runs are sandboxed (F6); `acquire.py` is an
   excluded command and asks each time it runs.
9. **Sections for long papers.** P1 splits by Marker's headings; claims
   spanning sections are merged by E3.

## Decisions (2026-10-27)

- **Decided by the user** (all four as recommended):
  - F2: the archive goes to the iCloud course folder, `ME295/library_archive/`.
  - F5: Codex runs as a separate macOS user.
  - F10: plants are split into tuning and held-out sets.
  - F25: the user rules on 30 practice-library claims.
- **The other 21 findings and the 9 unclear points** were applied as
  recommended in revision 7, for the user's review.
