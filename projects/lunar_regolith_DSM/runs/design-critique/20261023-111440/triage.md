# Triage: design critique 20261023-111440 (round 2)

Design: `workflows/literature-review/WORKFLOW.md` revision 3 (commit 70e61f0e;
sha256 in `manifest.json`). Critic: `beth-critic` (claude-opus-5-5), fresh,
with no knowledge of round 1: 22 findings (2 critical, 11 major, 9 minor) and
9 unclear points (`critiques/beth-critic.json`). Transcript audit passed: one
Read of WORKFLOW.md, seven web searches.

Checks made by the main session, 2026-10-23:
- **F1**, on the Crossref REST API, with the retracted DOI 10.1021/am300292v:
  the original's record has no `update-to`; its `updated-by` lists the
  retraction (source `retraction-watch`), and `filter=updates:<DOI>` finds
  the notice 10.1021/acsami.9b11759, whose `update-to` points back.
  **Confirmed**: the design names the wrong field.
- **F10**, by test: a headless Haiku session in a clean folder (`--tools Read
  --permission-mode dontAsk --setting-sources ""`) asked to Read
  WORKFLOW.md by absolute path was refused ("Permission to use Read has been
  denied"). Also, `scripts/units_verification/agents.py` `_audit` already
  fails any Read of a path not handed to the agent. **Refuted**; the design
  text did not say this, so it is clarified.
- **F18**, on help.openalex.org/api-reference/authentication (updated
  2026-10-15): "basic queries with no key at all"; a free key gives "$1 a
  day"; headers X-RateLimit-Remaining, -Limit, -Credits-Used. A keyed request
  returned `x-ratelimit-limit: 10000`, `x-ratelimit-credits-used: 1` for a DOI
  filter. The design's statement is **correct**; the critic's snippets were
  outdated. Its monitoring advice is kept.
- Other sources were search results only, not opened.

Recommendation: A = accept, M = accept with a change, R = reject, D = the
user decides. The user's decision is the last column.

| ID | Sev | Finding (short) | Rec | Proposed change | User |
|---|---|---|---|---|---|
| F1 | crit | Retraction status read from the wrong Crossref field | A | Read `updated-by` on each DOI and query `filter=updates:<DOI>`; test on 10.1021/am300292v; plant a retracted paper in C0 || accepted (2026-10-23) |
| F2 | crit | Claims and memos exist only on one local disk; no backup | D | An archived snapshot of the library after every publish, hashed in the event log, copied to a private off-machine location the user chooses; L2 checks the latest snapshot's hash. (Committing claims and memos was declined on 2026-10-23.) || backup to the ME295 course folder in iCloud (2026-10-23) |
| F3 | maj | `questions[]` and tags decide what Y1/Y3 see, and nobody checks them | A | E4 confirms `questions[]` against the brief; Y3 also gets a second retrieval by full-text and tag search and reports claims missing from the question-ID set || accepted (2026-10-23) |
| F4 | maj | Papers from earlier waves are never re-tagged for new or changed questions | A | Versioned question IDs (`Q2@v1`); a re-tag pass over existing claims whenever a question is added or changed (agent proposes, E4 confirms), recorded in the manifest || accepted (2026-10-23) |
| F5 | maj | Memos keep citing claims later corrected or retracted | A | Citations carry the claim version (`[id@v1]`); L2 re-runs `check_citations.py` on every memo and returns to Y1 any memo citing a superseded or retracted version || accepted (2026-10-23) |
| F6 | maj | C0's pooled 90 % can hide a kind that is mostly missed | M | Per kind at least 3 of 4 caught, and all fabricated quotes, wrong pages and sign/exponent errors; overall ≥ 90 %; false alarms on controls ≤ 10 %; lower bound reported || accepted (2026-10-23) |
| F7 | maj | U1 and U3b report bounds but gate nothing; the remedy re-uses the checker that failed | M | U3b passes at k = 0 of 60 (bound 4.9 %). k ≥ 1: revise E4's prompt with the errors found, re-run C0 including them, re-check that kind for the wave, then draw a fresh 60. U1: any relevant paper in the exclusion sample → revise screening and re-screen | accepted as recommended (2026-10-23) |
| F8 | maj | E3's ground truth is OCR; a correct quote from the image can fail it | A | A quote corrected by E4 passes E3 by fuzzy match to the Marker span (edit distance limit), both strings recorded; §1 states the script matches OCR text || accepted (2026-10-23) |
| F9 | maj | Conditions on method or results pages are invisible to E4 | A | Reader records the page of each condition; E4 also gets the method pages and every condition page; C0 plants dropped conditions outside the abstract, model and conclusions pages || accepted (2026-10-23) |
| F10 | maj | Read not confined to the agent's folder | R | Tested: refused in dontAsk; audit checks paths. Design text now states both || rejected (test) (2026-10-23) |
| F11 | maj | Y1/Y3 get every claim in one context, which overflows for large reviews | A | Token limit per call; above it, split by sub-question or tag facet; a script checks every claim was presented to some chunk || accepted (2026-10-23) |
| F12 | maj | No failure or resume rules for most agent stages | A | Every driver checks each input item has exactly one output (partial = failure); one retry, then the user; a call counts as done only after the check; L1 refusal counts shown at U4 || accepted (2026-10-23) |
| F13 | maj | Marker/Surya output may differ across machines | A | Test determinism at build (same machine twice; second machine if available); record the backend and inference settings; the stored Marker pages and their hashes are the reference for E3 || accepted (2026-10-23) |
| F14 | min | Number matching in Y2 unspecified | A | Define the extractor, normalisation (1e-3 = 10^{-3}, %), exemptions (IDs, years, equation numbers), and test cases || accepted (2026-10-23) |
| F15 | min | Same statement with different quotes becomes two claims | A | Second matching pass by claim-text similarity on the same page; likely duplicates go to E4 to merge; omission rates before and after || accepted (2026-10-23) |
| F16 | min | Screening calibration has no borderline cases | A | 10–20 borderline relevant abstracts chosen by the user; recall threshold on them | accepted; revisit if too costly (2026-10-23) |
| F17 | min | Recalibration misses Marker/Surya/mutool and Claude Code minor updates | A | Add P1 tool or weight changes; on any Claude Code version change run a 10-item smoke subset of C0, full C0 if it fails || accepted (2026-10-23) |
| F18 | min | OpenAlex limits "disputed" | M | Design's facts verified; add: read X-RateLimit-Remaining, estimate a wave's cost before D1, prefer filter and ID lookups, fail fast if the Keychain read fails || accepted (2026-10-23) |
| F19 | min | `--setting-sources ""` drops project deny and sandbox rules for headless agents | A | Pass deny rules (and sandbox config if adopted) with `--settings <file>` on every agent call; test that the librarian cannot write to `library/` || accepted (2026-10-23) |
| F20 | min | Re-check agreement has no threshold | A | Pass ≥ 90 % agreement on the 30 re-checked claims; below → re-run E4 on the wave and re-calibrate || accepted (2026-10-23) |
| F21 | min | Preprint and journal version can both be read in one wave | A | D3 merges by OpenAlex work ID, Crossref `relation` (has-preprint / is-preprint-of), and fuzzy title plus first author || accepted (2026-10-23) |
| F22 | min | Gold recall set is one paper, 10–15 items | A | 2–3 papers of different types (article, long report, thesis), about 30 statements; bound reported | accepted; revisit if too costly (2026-10-23) |

Unclear points, proposed resolutions:

1. Screener framing B may exclude only when the candidate bears on no
   question *and* an exclusion criterion applies; otherwise context.
2. User-authored claims (U3a "missed") go through E3 and E4 like any other,
   so the acceptance test holds.
3. Disputed claims are shown at U3b in addition to the 60; "U3b is empty"
   applies to its random part only.
4. E4 cannot search for the right page; its verdict is "not on the cited
   page", and E5 (which has the whole paper) finds the page, then E3/E4.
5. Long papers read per section: a script merges the section maps; E4's
   context pages come from the merged map; a claim spanning sections cites
   its main page and lists the other pages under conditions.
6. Title similarity: `difflib.SequenceMatcher` ratio on titles normalised
   (NFKC, lower case, punctuation removed, whitespace collapsed) ≥ 0.9.
7. Memos are built when the user asks at U4 or when a review closes, from
   claims published by L1 of the latest wave.
8. E5 gets the claim, the checker's reasons, the cited page, the condition
   pages and their Marker text; the re-check is a fresh E4 that does not see
   the revision history.
9. `search/` snapshots are always committed, gzip-compressed (bibliographic
   metadata; OpenAlex data is CC0). Abstracts from Crossref may carry
   publisher copyright, so abstracts are stored separately and not committed;
   their hashes are.
