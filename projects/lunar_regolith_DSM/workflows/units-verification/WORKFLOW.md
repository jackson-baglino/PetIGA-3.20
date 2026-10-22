# WORKFLOW: units-verification

Status: **design, for review** (2026-10-22). Nothing is built yet.
Flowchart: `flowchart.html` in this folder, rendered from `flowchart.json` by `scripts/flowchart/render_flowchart.py`.

## 1. Purpose

- **Question.** Is each row of `memory/semantic/units.md` (62 symbols from
  Kaempfer & Plapp 2009, Moure & Fu 2024 and its SI) correct: the unit, the
  meaning, the page, and the printed value where there is one?
- **Finished output.** `units.md` with every row either `verified` (page
  image, run ID, date) or `disputed` with the reason, and the UNRELIABLE
  header removed only when no row is unresolved.
- **Acceptance test.**
  - Every row has a final status.
  - The planted-error calibration (§4) passes before the real run.
  - The user has inspected every disputed row and a random sample of
    agreed rows.
- **Baseline.** The practice-round registry: one extractor plus one image
  check. Its rows still contain merged meanings and mangled quotes, so the
  number of rows this run corrects measures what the design adds.

## 2. Stages

| # | Stage | Agent or script | Input | Output | Check | Stop / escalate |
|---|---|---|---|---|---|---|
| U0 | Supply sources | user | — | PDFs of the three sources in `Literature/` | S0 finds each by DOI | — |
| S0 | Locate sources | `locate_sources.py` | `units.md` source ids, `Literature/` | `manifest.json`: path, sha256, DOI per source | first-page text names the paper; files whose rendered pages are identical count as one source | missing or ambiguous source → stop, ask user |
| S1 | Render pages | `render_pages.py` (`mutool draw -r 200`) | manifest; cited page per row, plus the page before and after | `pages/<source>_p<N>.png`, hashes | file exists, hash recorded | render failure → stop |
| C0 | Calibrate (once, before S2 on real rows) | S2 → S3 → S4 on a copy of 8 rows with planted errors | 8 altered rows | calibration report | S3 flags 8/8; S4 rejects ≥ 7/8 planted units | below threshold → stop, revise prompts |
| S2 | Blind read | `morty-unit-reader` | one page image + the symbol names expected on it (no units, no values) | `reads/<page>.json` per symbol: printed form, meaning, value as printed, unit as printed, where on the page, or not found | JSON schema; every requested symbol answered | schema failure → re-ask once, then log |
| S3 | Compare | `compare_units.py` | table row, blind read | AGREE / DISAGREE / NOT_FOUND per row | units parsed to dimension vector and scale (J m⁻³ K⁻¹ ≡ J/(m³ K); °C difference ≡ K); NOT_FOUND re-read on neighbour pages | — |
| S4 | Adjudicate | `beth-unit-checker` | for each DISAGREE/NOT_FOUND row: page image, table row, blind read | verdict: table right / read right / both wrong (with correction) / cannot tell, plus location on the page | verdict schema; a correction must re-pass S3 | "cannot tell" → user |
| U1 | Inspect | user | every row from S4 + a random 10 % (6 rows) of AGREE rows, each with its page crop | approve or correct | — | user decides |
| S5 | Apply | `apply_verdicts.py` | verdicts, user decisions | updated `units.md`, run manifest | every row has a status; header removed only if none unresolved | — |

## 3. Agents

```
name:         morty-unit-reader
model:        claude-sonnet-5-5
tools:        Read (the one page image it is given)
skills:       none
MCP servers:  none
reads:        one page PNG; list of symbol names
writes:       JSON array, one object per requested symbol (schema in §6)
may not:      see the table's units or values; use the PDF text layer; guess a
              symbol that is not on the page (return found=false)
context:      fresh per page
```

```
name:         beth-unit-checker
model:        claude-opus-5-5
tools:        Read (page image)
skills:       none
MCP servers:  none
reads:        page PNG, the table row, the blind read
writes:       {symbol, verdict, unit, value, meaning, location, reason}
may not:      fix rows it was not given; answer from memory of the paper
context:      fresh per disputed row
```

Each agent call runs headless (`claude -p`, `--permission-mode dontAsk`,
Read only) from a driver script; the exact invocation is fixed at build time.

## 4. Checks

| Output | Check | Kind | Pass criterion |
|---|---|---|---|
| source files | DOI on first page, sha256 | deterministic | one file per source |
| blind reads | schema; S3 unit comparison | deterministic | parses; dimensions and scale equal |
| disputed rows | `beth-unit-checker` | independent agent (different model, sees the image) | verdict with a location |
| the whole pipeline | planted errors (C0): 8 rows with a flipped exponent sign, a wrong prefix, a swapped unit, a wrong page | measured | S3 8/8; S4 ≥ 7/8 |
| reader stability | read 5 pages twice | measured | report the disagreement count |
| final table | user inspection of all disputed + 10 % of agreed rows | human | user approves |

The reader is blind to the claimed answer, so agreement with the table is
evidence rather than echo. Comparison is a script, so the same reads always
give the same verdicts.

## 5. Human decision points

- **U0:** confirm which file is each source when S0 cannot decide (see §8).
- **U1:** each disputed row and the agreed sample, shown with the page crop,
  the table row and the read side by side.

## 6. Reproducibility

`runs/units-verification/<run_id>/`:

- `manifest.json`: git commit, Claude Code version, model IDs, sha256 of the
  agent files, every PDF and every page PNG, check verdicts.
- `reads/*.json`, `verdicts.json`, `calibration.json`: committed.
- `runs.db`: one row per agent call (stage, agent, model, input hash, output
  path, verdict, time); local.
- `transcript.jsonl`: local.

Read schema (one object per symbol):
`{symbol, found, printed_as, meaning, value_as_printed, unit_as_printed, location}`.

## 7. Risks and failure modes

| Risk | Caught by |
|---|---|
| Text layer drops exponents and minus signs | agents get images only |
| Reader anchors on the claimed unit | blind read (no units or values given) |
| Wrong page in the table | neighbour pages rendered; NOT_FOUND re-read |
| Subscripts and Greek letters misread | S3 compares dimensions, not strings; S4 sees the image |
| "Verified" rows that were never re-checked | every row passes through S2–S3 again |
| Two agents share a blind spot | S4 uses a different model; user inspects a random sample |

## 8. Open decisions for the user

1. **Sources:** resolved 2026-10-22. Kaempfer & Plapp (2009,
   10.1103/PhysRevE.79.031502) and Moure & Fu (2024, 10.1021/acs.cgd.4c00539)
   were copied from the user's Zotero library into `Literature/` as
   `kaempfer_plapp_2009.pdf` and `moure_fu_2024.pdf`; the DOI on each first
   page matches.
   (Resolved: the SI's two copies render pixel-identically on all 12 pages;
   `cg4c00539_si_001.pdf` was re-saved by macOS Preview on 2026-10-18, which
   changed only its file metadata. The version of record is the untouched
   `MoureFu2024_SI/moure_fu_2024_SI.pdf`.)
2. **Value column:** also record the printed value (read from the image)? Recommended.
3. **Scope:** literature symbols only, or also the solver's own parameters
   (units from `lunar_main.c` option help text)? Recommended: literature
   first, solver parameters as a second pass with a code-reading checker.
4. **Models:** Sonnet reader and Opus adjudicator, as above?
5. **Sample size** for U1: 10 % of agreed rows?
