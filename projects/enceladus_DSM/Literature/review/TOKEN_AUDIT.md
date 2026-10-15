# Token audit: literature fact-check, first iteration (2026-10-14/15)

The figures come from the subagent transcripts, which were parsed by `audit.py` and `coverage.py` in the session scratchpad.
Reading was stopped on the user's request before the cross-check and the outline update.

## Where the tokens went

| Agent | API calls | Cache read | Cache write | Output | Mean context per call | Max context |
|---|---|---|---|---|---|---|
| Collector | 33 | 2.3 M | 0.11 M | 0.8 k | 73 k | 106 k |
| Checker | 64 | 7.9 M | 0.18 M | 5 k | 126 k | 194 k |
| Reviewer A | 188 | 72.5 M | 2.7 M | 61 k | 400 k | 967 k |
| Reviewer B | 250 | 98.3 M | 2.5 M | 20 k | 403 k | 957 k |
| Reviewer C | 208 | 82.6 M | 2.8 M | 53 k | 411 k | 963 k |
| **Total** | 743 | **263.7 M** | 8.2 M | 140 k | | |

The three reviewers account for **97 %** of the tokens. Collecting and checking the library cost about 4 %.

## Bottlenecks, ranked

1. **Context build-up in long-lived reader agents (the main cost).** Every API call re-reads the whole
   conversation so far. Once a paper has been read, its text stays in context and is billed again on
   every later call. As a result, each reviewer saw only about 2–2.7 M characters of paper text
   (about 0.6–0.8 M tokens), yet consumed 75–100 M tokens. That is an amplification of **about 100–130×**.
   Mean context ran at about 400 k tokens per call and reached about 960 k before each compaction.
2. **Reading the library three times.** Independence matters for the *verdicts*. It does not matter
   for extracting what each paper says. Three full reads tripled the dominant cost.
3. **"Read every paper in full" was not achieved anyway.**

   | Reviewer | Papers with ≥ 80 % of text seen | Papers with < 20 % seen | Stub notes (< 700 B) |
   |---|---|---|---|
   | A | 36 | 70 | 53 |
   | B | 25 | 71 | 1 (stopped at 87 of 111 papers) |
   | C | 38 | 63 | 43 |

   Many stub notes say "scanned" (for example `A/notes/hemingwayInteriorEnceladus2018.md` and
   `chungTimeIntegrationAlgorithm1993.md`), but the transcripts show little or none of those papers'
   text was read. The brief's "read in full" instruction lost out to the context pressure.
   **Treat the notes for uncited papers as title- or abstract-level only.** The notes for the cited
   papers, which came first in every reading order, are much better grounded.
4. **Usage-limit stops (three times).** Three Opus agents running in parallel at about 400 k context each hit the
   session limit within 1–2 hours each time. Every resume rebuilt the context.
5. **Re-reading shared inputs.** The brief, outline, draft, MANIFEST and CHECK_COLLECTION were re-read
   after each compaction or resume: 0.2–0.6 M characters per reviewer. Bash outputs were often 20–30 k characters
   (whole papers via `cat`, or concatenated notes).
6. **Unreliable text layer (minor).** Dropped minus signs, μ and × forced page renders, but only about
   17–21 images were rendered per reviewer, so this was cheap.

## Recommendations for the next iteration

- **Use one short-lived extractor per paper (or per batch of about 3–5 papers) instead of long-lived readers.**
  Give each extractor the claims list (compact, about 10 k tokens) plus one paper's text, and have it write a
  structured `notes/<key>.md` (claims supported or contradicted, numbers with page locators). Context stays at
  about 30–80 k, so the cost scales with the library size rather than its square. The expected saving is **about 10–30×**.
- **Extract once, judge three times.** Run extraction a single time. Run the three independent reviewers only on
  the verdict step, using the claims table, the notes and targeted page checks, so their contexts stay small.
- **Make the work claims-first.** One agent builds the claim and citation list from the outline and draft
  before any reading starts. Cited papers get a full read. Uncited papers get a screen of abstract,
  conclusions and figures for "suggested additions", plus a full read only when the screen flags them.
- **Choose the model by step.** Use a cheaper model (Sonnet or Haiku) for the per-paper extraction. Keep Opus for the
  verdicts, the cross-check and the outline edits.
- **Check reading coverage automatically.** Run `coverage.py` against the transcripts before trusting any notes,
  and reject stub notes.
- **Limit concurrency and output size.** Run at most 2 large agents at once. Cap Bash output, for example
  with `head -c 15000`, and split papers by page range rather than using `cat` on whole files.
- **Clean the text once.** Pre-render the table and equation pages of the papers flagged in CHECK_COLLECTION, and
  store the corrected values in one shared file so no agent has to re-derive them.

## State at stop

- Notes: A 111 (53 stubs), B 87 (1 stub), C 111 (43 stubs).
- Draft claims tables: `A/claims_table.md` and `C/claims_table.md`. B has none.
- **Not done:** the `manuscript_outline_new_{A,B,C}.tex` files, the cross-check and the consolidated outline.
- **Before any outline update:** ask the user for the updated manuscript, as the user requested on 2026-10-15.
