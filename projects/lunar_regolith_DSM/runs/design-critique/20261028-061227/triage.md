# Triage: design critique 20261028-061227 (round 4)

Design: `workflows/literature-review/WORKFLOW.md`, revision 7 with the edit of
2026-10-28 (whole paper in each library record), commit eb710362. Critic:
`beth-critic` (claude-opus-5-5), fresh, with the questions in
`runs/design-critique/focus/rev7_whole_plan.md`: weakest points, review mode,
the two-layer library, calibration statistics, simplicity. Result: 25
findings (0 critical, 19 major, 6 minor) and 13 unclear points
(`critiques/beth-critic.json`).

Checks made by the main session, 2026-10-28:
- **F20**, against platform.claude.com/docs/en/build-with-claude/vision: the
  2,576 px tier applies to "Claude 4.7 and later models", which includes the
  pinned Sonnet 5.5 reader and Opus 5.5 checker. **Refuted.** The same page
  gives a limit the design misses: a request with more than 20 images caps
  each image at 2,000 px ("resize each image so that neither dimension exceeds
  2000 px, or keep the request to 20 or fewer image and document blocks").
  Tiles (1,650 px) pass; full page images at 2,576 px would be rejected once
  an agent has read more than 20 images.
- **F1** (a sandboxed process cannot apply another macOS sandbox profile): a
  nested `sandbox-exec` ran in this session, but this session has no sandbox
  enabled, so the test says nothing. **Not confirmed; plausible.** T-sandbox
  tests it at build.
- **F16** (Codex sees images only through `-i`): consistent with the
  `codex exec --help` read on 2026-10-27. Accepted without a further check.

Recommendation key: **accept**; **modify** (accept with the change stated);
**user** (needs a decision only the user can make); **reject** (with reason).

| # | Sev. | Finding (short) | Rec. | Change |
|---|---|---|---|---|
| F1 | major | `run_offline.sh` cannot nest a sandbox; excluding it lets arbitrary commands escape | accept | a fixed launcher that takes only a script name from an allow-list (no arbitrary argv), runs outside the Claude sandbox with its own profile, writes stdout and stderr to the run folder, and returns only the schema-checked summary; T-trifecta adds a planted traceback and a `cat` through the launcher |
| F2 | major | Review mode's read-deny list misses many stores of untrusted text | accept | invert it: in review mode, Read and Bash reads are allowed only on a short list (summaries, scripts, settings, briefs, WORKFLOW.md); everything else is denied. T-trifecta plants an instruction in each untrusted store |
| F3 | major | Outside waves the main session has all three legs (desk, U4, triage) | user | see Decisions below |
| F4 | major | The checker is anchored on the reader's (often Marker's) transcription | accept | for equation and parameter claims, E4 first transcribes blind (given only the quantity and its position), then sees the claim and compares; the blind transcription is the checker's field in X1 |
| F5 | major | Branch calibration with ~10 plants can never meet the 85 % bound | modify | the branch's light calibration is reported (all caught or not, with its interval) and is not called a pass; the gate is the core calibration. A branch becomes a gate only once its set holds ≥ 19 held-out plants |
| F6 | major | Pooling script-caught plants with judgment plants measures E3, not E4 | accept | kinds E3 catches by rule must all be caught; the bound is computed only on plants that reach E4. Held-out size fixed in advance (≥ 19 for an 85 % bound with no misses) |
| F7 | major | The held-out set is reused after each prompt revision | accept | a held-out set is scored once per prompt version, and its misses are never shown to whoever revises the prompts; after a failure a fresh held-out set is drawn from the growing pool, and the number of attempts is reported with the bound |
| F8 | major | Adding E4-found errors to the set inflates E4's catch rate | accept | only errors found independently of E4 (X1 resolved by the user, U3a, U3b, Y3) enter the held-out pool; E4-found errors go to the tuning set |
| F9 | major | Core calibration has no scanned or older papers | modify | plants record their copy kind; the core set includes at least one scanned paper (script and Codex plants, user-approved); a branch whose papers are mostly scans adds plants of that kind |
| F10 | major | U3b's uniform sample barely touches the claims memos rely on | accept | U3b gets two strata: every claim cited by a model-critical memo statement (a census, or a separate 60 if there are more), and the uniform 60 for the rest; a bound is reported per stratum |
| F11 | major | Y2 matches a memo number anywhere in the cited claim | accept | each memo number carries its claim ID and quantity in a structured inline tag; Y2 matches value, unit and quantity against that claim's parameter field |
| F12 | major | The memo checker can judge support from unchecked passages | accept | Y3 is two calls: support is judged only against the cited claims; passages are used only to list omissions and contradictions. The schema has no field for passage-based support |
| F13 | major | Nearly every B1 and B5 claim is "important", so X1 and disputes grow with papers | user | see Decisions below |
| F14 | major | "Normalized LaTeX" equality is undefined; formatting differences flood U3 | accept | a token-level canonicalizer with a tested list of equivalences, with its false-dispute rate measured in T-res; differences only in formatting go to E4, not to the user |
| F15 | major | The Codex binary and login sit in the main user's home folder | accept | a pinned Codex binary in a root-owned shared path, one login as the separate user, and a root-owned launcher allowed in sudoers for that one command; T-codex tests it |
| F16 | major | Codex sees images only when attached with `-i` | accept | tiles attached with `codex exec -i`; the audit allows zero shell commands; T-codex gets an X1 case with a known answer |
| F17 | major | The symbol list is in the checked layer, but E4 checks only symbols claims use | accept | symbols move to the content layer (unchecked) unless a claim uses them; E3b takes reference units only from the verified units table, never from inferred dimensions |
| F18 | major | Supplementary material and corrigenda are not stored | accept | A1 also fetches linked supplements and corrigenda (Crossref and JATS relations) into `source/`; P1 runs on them under the same paper ID; the card lists any that could not be acquired |
| F19 | major | Claims on a newly flagged retraction stay usable until the user rules | accept | Y1 excludes, and Y2 fails, any claim on a flagged paper until the user rules |
| F20 | minor | Readers may be capped at 1,568 px | reject | refuted by the docs (above). Added instead: every image an agent reads is at most 2,000 px on its long edge, because agents read more than 20 images per call; page images are downscaled for agents, and tiles carry the detail |
| F21 | minor | Full-library archives and backups grow quadratically | accept | one content-addressed store: each paper copy archived once by its hash, plus a per-wave manifest of hashes; backup and archive merge into it. Also answers unclear point 12 |
| F22 | minor | Claude Code is not version-locked | accept | Claude Code pinned with auto-update off; the version recorded per call in `runs.db`; a mismatch stops the wave |
| F23 | minor | Pinned models will be retired | accept | the design states that agent verdicts are reproducible through their archived transcripts; dated model IDs where offered; each model's deprecation date recorded |
| F24 | minor | The second table extractor adds nothing | reject | with the whole paper stored, tables not cited by any claim are content that is never checked by E4; the two-extractor comparison is the only check on them, and it marks the cells that disagree |
| F25 | minor | Several diagnostics cannot block publication | modify | the 30-claim repeatability sample runs on the first two waves and after any prompt or model change, not every wave; D2c is already on demand and Z1 deferred. Both screeners stay until SC0's agreement is measured (as decided in round 3) |

Unclear points, answered in revision 8:
1. **E4 and Marker text.** The disagreement list tells E4 where to look (page,
   position, kind), not what Marker read.
2. **Quotes from LaTeX or XML.** `pages/p001.md` is the source text split at
   the page boundaries found by aligning it with the PDF's text (F4 of round
   3); a quote is matched on that page's text, and an alignment failure is
   its own flag.
3. **Memo sentences without claims.** Model options, *library silent* and gap
   statements are marked as such and exempt from Y2's citation rule; gaps
   list their queries (D2g).
4. **Mechanical plants before any wave.** They are made from the user's gold
   statements, which are verified by the user.
5. **Kinds with 8 or more plants.** The kind then must meet the same bound as
   the overall set.
6. **acquire.py.** It is an excluded command and asks each wave; its URLs
   come only from OpenAlex and Crossref records, never from the main session.
7. **Drivers in review mode.** Drivers run through the launcher (F1) outside
   the session's sandbox, with their own profile allowing only the Anthropic
   and OpenAI hosts and their run folder.
8. **U1's sample.** It includes both rule exclusions and screener
   exclusions, sampled in proportion.
9. **Tier S.** Branch C0, X1, U3a and U3b still apply.
10. **Several copies.** Every copy is indexed; search returns the copy the
    card marks as current, and lists the others.
11. **Repeatability failure.** E4 reruns on the wave's passed claims under
    the revised prompt, after C0 passes it.
12. **Backups without images.** Answered by F21: the content-addressed store
    holds every copy's images once.
13. **Units table.** `memory/semantic/units.md` came from the
    units-verification workflow (54 verified rows); its hash goes in the
    manifest and E3b refuses to run on a mismatch.

## Decisions for the user

- **F3: the main session outside waves.** Between waves, the reference desk
  returns verbatim passages into the main session, which has web access. That
  is all three legs. Options:
  - (a) *recommended:* review mode whenever a session reads library or run
    text, and the desk runs as a headless call with no network, writing its
    hits to a file the user opens; the main session gets only IDs and pages.
    The main session can then not discuss passage text with you.
  - (b) review mode for every session in this project; web access only in a
    separate session that never opens library files.
  - (c) accept the risk outside waves, and keep it in the residual list.
- **F13: which claims get the Codex read (X1).** Options:
  - (a) *recommended:* X1 on claims a memo cites and claims on model-critical
    questions, plus a random 20 % of the other important claims to measure
    what is skipped; decide after the shakedown wave's disagreement rate.
  - (b) keep X1 on every important claim.

## Decisions (2026-10-28)

- **Decided by the user:** F3 option (a) and F13 option (a), as recommended.
- **The other findings and the 13 unclear points** were approved and applied
  as recommended in revision 8 (commit 557c851b).
