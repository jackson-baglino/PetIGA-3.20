# Triage: design critique 20261023-103742

Design: `workflows/literature-review/WORKFLOW.md` (sha256 in `manifest.json`, commit 7686e258).
Critics: `beth-critic` (claude-opus-5-5) returned 23 findings and 8 unclear points
(`critiques/beth-critic.json`). `rick-critic` (claude-fable-5-1) failed on both
attempts: "Usage credits are required for this model." (the account has no
access to Fable).

Checks made by the main session, 2026-10-23:
- Statistics (F10, F11), with scipy: the Clopper-Pearson two-sided 95 % upper
  bound for 1 error in 60 is 0.0894 (one-sided 0.0766); the 95 % CI for 5/6 is
  0.359–0.996. The critic's figures are correct.
- F15, on help.openalex.org/how-to-use-the-api/rate-limits-and-authentication
  (updated 2026-10-15): basic queries work without a key; a free key gives "10×
  the keyless budget" ($1 a day); an exhausted budget returns 429. Confirmed.
- F4, on crossref.org/documentation/retrieve-metadata/retraction-watch/:
  retractions are in the REST API `update-to` field, with source `publisher` or
  `retraction-watch`; corrections are included "but these are not as
  comprehensive as retractions". Confirmed for retractions, partly for
  corrections.
- The critic's other sources are search results only (it had WebSearch, not
  WebFetch); they were not opened.
- Transcript audit passed: pinned model, one Read of WORKFLOW.md, six web
  searches.

Recommendation: A = accept, M = accept with a change, D = the user decides.
The user's decision is the last column.

| ID | Sev | Finding (short) | Rec | Proposed change | User |
|---|---|---|---|---|---|
| F1 | crit | Reader copies Marker's misread equation; E3 agrees; image never consulted | A | E4 checks every equation, number and quote of every claim against the image; the checker gets the image only, never Marker text | |
| F3 | crit | basis, kind, tags, questions never checked | A | Add basis and kind to the E4 verdict; plant wrong-basis errors in C0; show them at U3 | |
| F4 | crit | Retractions and errata checked only for non-publisher copies | A | Query Crossref `update-to` for every DOI at D3; record on the card; re-query at L2 | |
| F2 | maj | Sonnet, Opus and Fable are all Claude: not "different families" | M | Reword to "different Claude models"; state that the U3 sample is the only fully independent check; another vendor's checker is a separate decision | D |
| F5 | maj | Librarian agent chooses the synthesizer's claims, unchecked, not reproducible | A | Y1 and Y3 get claims from `library_query.py` with fixed filters; query and returned IDs in the manifest; librarian agent only answers ad hoc questions | |
| F6 | maj | Y3 sees only cited claims, so cannot find omitted contradictions | A | Y3 gets the full claim set for the question | |
| F7 | maj | "Same claims on re-run" is not measurable | A | Reproducibility = identical hashes for deterministic stages + re-checking recorded claims gives the same verdicts at a stated agreement rate | |
| F8 | maj | No rule for matching reader A's and B's claims | A | Match by page and overlapping quote span; unmatched claims go to E4 as candidates; merged set defined before L1 | |
| F9 | maj | Checker sees only the cited page; conditions stated elsewhere are lost | A | E4 also gets the abstract, assumptions and conclusions pages; C0 plants "dropped condition" errors | |
| F10 | maj | Stratified U3 sample breaks the 3/n bound; no bound for k > 0 | A | Simple random sample of 60; Clopper-Pearson bound for any k; inspect all when a wave has ≤ 60 passed claims | |
| F11 | maj | C0 too small to measure 90 % | A | 30–50 planted errors over ≥ 3 papers, some written by the user; report the rate with its CI | D (user's time) |
| F12 | maj | Short or generic quotes pass E3 | A | Quote ≥ 8 words or a full sentence; flag quotes that occur more than once in the paper | |
| F13 | maj | Many closed-access papers have no abstract in OpenAlex | A | Record abstract presence; never exclude a candidate without one; fetch from Semantic Scholar or arXiv, else show at U1 | |
| F15 | maj | OpenAlex keyless budget is small; 429 truncates D1 | A | Free OpenAlex key stored outside the repository (needs an account); 429 → pause and resume | D (account) |
| F16 | maj | Failed API → wave continues silently | A | Mark the wave "discovery incomplete", show at U1, re-run failed queries | |
| F17 | maj | Curator crash leaves half a publish | A | Stage in a temporary folder; event appended last as the commit marker; `recover` mode | |
| F14 | min | Two Haiku screeners share errors | M | Keep two, with different framings; exclude only on an explicit exclusion criterion | |
| F18 | min | Marker output depends on version and weights | A | Pin Marker, Surya and weight hashes and flags in the manifest (as units-verification) | |
| F19 | min | Context-paper cards unchecked | A | Cards hold bibliographic fields only; any statement becomes a claim and is checked | |
| F20 | min | Y2 number matching is string-based | A | Number must occur in a cited claim's quote; Y3 checks it is the same quantity | |
| F21 | min | IDs tied to DOIs; preprint → journal not linked | A | Stable internal ID with DOI aliases; A1 links a new version to its preprint and triggers the re-check | |
| F22 | min | Reports and theses without Crossref records dropped as invented | A | Also resolve via arXiv, NTRS, OSTI, ADS; unresolved items shown at U1, not dropped | |
| F23 | min | U1 10 % of exclusions grows with candidates | A | Fixed random sample of 60 exclusions | |

Unclear points, proposed resolutions:

1. Memo closing with unsupported statements: allowed only if the user accepts
   each, and the memo marks them.
2. INFERENCE sentences cite the claims they are inferred from; Y3 judges
   whether the inference follows.
3. §5 "below": refers to the C0 paragraph; reword.
4. Checker and Marker text: the checker gets images only (see F1).
5. `page`: record both the printed page number and the PDF page index.
6. E4 batching: per page, all claims on that page; drop "fresh per item" and
   "other claims' verdicts" wording.
7. Fable fallback: Y3 on Opus is a weaker check; say so in §5 (see F2).
8. Recalibration triggers: any change to an agent file, prompt, schema, model
   ID, or Claude Code major version.
