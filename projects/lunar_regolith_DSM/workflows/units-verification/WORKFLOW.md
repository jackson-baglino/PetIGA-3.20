# WORKFLOW: units-verification

Status: **being built in stages** (approved 2026-10-22). Stage 1 (S0, S1, S1b, S3 with tests) built;
stage 2 (agents, S2, S4, S4m, C0) built on `feat/units-agents`, awaiting the user's approval of the agents;
stage 3 (real run, U1, S5) not yet.
Flowchart: `flowchart.html` in this folder, rendered from `flowchart.json` by `scripts/flowchart/render_flowchart.py`.

## 1. Purpose

- **Question.** Is each row of `memory/semantic/units.md` (58 literature symbols
  from Kaempfer & Plapp 2009, Moure & Fu 2024 and its SI) correct: the unit,
  the printed symbol, the meaning, the page, and the printed value where there
  is one? (`j_T` and `j_v` were removed on 2026-10-22: post-processing fluxes,
  not model parameters.)
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
| S1b | Marker reading | `marker_extract.py` (Marker 1.10.2, forced OCR, no LLM) | every page of each source | `marker/<source>_p<N>.md`: page as Markdown with LaTeX math; a symbol → page index | Markdown written for every page | Marker failure on a page → that page has no Marker reading; logged |
| C0 | Calibrate (once, before S2 on real rows) | S2 → S3 → S4 → S4m on a copy of the table with 10 planted errors and 4 untouched control rows | altered rows | `calibration.json` | S3 flags 8/8 unit and page errors; S4 rejects ≥ 7/8; S4m catches 2/2 meaning errors; false alarms on the controls reported | below threshold → stop, revise prompts |
| S2 | Blind read | `morty-unit-reader` | one page image + the symbol names expected on it (no units, no values) | `reads/<page>.json` per symbol: printed form, meaning, value as printed, unit as printed, where on the page, or not found | JSON schema; every requested symbol answered | schema failure → re-ask once, then log |
| S3 | Compare | `compare_units.py` | table row, blind read, Marker reading | AGREE / DISAGREE / NOT_FOUND per row | units parsed to dimension vector and scale (J m⁻³ K⁻¹ ≡ J/(m³ K); °C difference ≡ K); printed symbols compared after dropping sub/superscript marks (ρi ≡ ρ_i; ρv ≠ ρ_vs); AGREE only if the table, the blind read and Marker all agree; NOT_FOUND re-read on the page Marker's index points to | — |
| S4 | Adjudicate | `beth-unit-checker` | for each DISAGREE/NOT_FOUND row: page image, table row, blind read, Marker reading | verdict: table right / read right / both wrong (with correction) / cannot tell, plus location on the page | verdict schema; a correction must re-pass S3 | "cannot tell" → user |
| S4m | Meaning review (after S3, on every row, independent of S4) | `beth-unit-checker` | pass 1, per page: the page image and, for every row cited on it, the table's printed symbol and meaning and the blind read's meaning; pass 2, per row that pass 1 could not settle ("cannot tell"): the cited page plus up to 3 pages where Marker finds the symbol, earliest first, including the main paper for SI rows (where the SI's symbols are defined) | per row: meaning matches / table wrong (with the page's wording and the page it is defined on) / cannot tell, plus location | verdict schema | a pass-1 "cannot tell" goes to pass 2 (or straight to U1 if Marker finds the symbol on no other page); anything but "matches" after that → the row is disputed and goes to U1 |
| U1 | Inspect | user | every row from S4 or S4m + a random 10 % of AGREE rows (rounded up), each with its page crop | approve or correct | — | user decides |
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
reads:        S4: page PNG, the table row, the blind read, the Marker reading
              S4m: page PNG, each row's printed symbol and meaning (table and read)
                   (pass 2: also the pages where the symbol is defined)
writes:       S4: {symbol, verdict, unit, value, meaning, location, reason}
              S4m: per row {symbol, meaning_verdict, page_wording, location, reason}
may not:      fix rows it was not given; answer from memory of the paper
context:      fresh per disputed row (S4); fresh per page (S4m pass 1) or per row (pass 2)
```

Each agent call runs headless from `scripts/units_verification/agents.py`,
in a fresh empty folder outside the project that holds only the page images
(so no CLAUDE.md, memory or units table is in reach):
`claude -p <task> --model <pinned> --tools Read --permission-mode dontAsk
--setting-sources "" --strict-mcp-config --no-session-persistence
--disable-slash-commands --system-prompt <agent body> --output-format
stream-json --verbose --json-schema <schema>`. In dontAsk mode the harness
refuses a Read outside that folder (tested 2026-10-22). Every transcript is
audited by script: the session model must be the pinned one and every tool
call a Read of a given image; the output is validated against the schema;
one retry, then the failure is logged.

## 4. Checks

| Output | Check | Kind | Pass criterion |
|---|---|---|---|
| source files | DOI on first page, sha256 | deterministic | one file per source |
| blind reads | schema; S3 unit comparison | deterministic | parses; dimensions and scale equal |
| Marker reading | a second, independent instrument (OCR model, not Claude) | deterministic given its pinned version | agrees with the blind read and the table |
| disputed rows | `beth-unit-checker` | independent agent (different model, sees the image) | verdict with a location |
| meanings, every row | `beth-unit-checker` (S4m) | independent agent | verdict with a location |
| the whole pipeline | planted errors (C0): 8 rows with a flipped exponent sign, a wrong prefix, a swapped unit, a wrong page; 2 rows with a swapped meaning; 4 untouched controls | measured | S3 8/8; S4 ≥ 7/8; S4m 2/2; control false alarms reported |
| reader stability | read 5 pages twice | measured | report the disagreement count |
| final table | user inspection of all disputed + 10 % of agreed rows | human | user approves |

The reader is blind to the claimed answer, so agreement with the table is
evidence rather than echo. Marker is a different kind of instrument (a
local OCR model), so the blind read and Marker do not share a model's
blind spots. Tested on Kaempfer & Plapp p. 4: it recovered every value the
text layer mangles (1.8 × 10⁶ J m⁻³ K⁻¹, 2.178 × 10⁻⁵ m² s⁻¹, …). Comparison is a script, so the same reads always
give the same verdicts.

## 5. Human decision points

- **U0:** confirm which file is each source when S0 cannot decide (see §8).
- **U1:** each disputed row and the agreed sample, shown with the page crop,
  the table row and the read side by side.

## 6. Reproducibility

`runs/units-verification/<run_id>/`:

- `manifest.json`: one entry per stage, each with its git commit, tool
  versions (Claude Code, mutool, Marker, surya, torch, model-cache folders),
  input hashes and verdict; sha256 of the agent files, every PDF, every page
  PNG and every Marker page. Committed.
- `marker_index.json`, `compare*.json`, `reads/*.json`, `verdicts.json`,
  `calibration.json`: committed (short extracted facts, not page text).
- `pages/*.png`, `marker/*.md`, `crops/`: **not committed** (user decision,
  2026-10-22). They are page images and page text of the papers. The run is
  reproducible from the DOIs and PDF hashes in the manifest: anyone with the
  papers re-renders the pages and re-runs Marker and checks the hashes.
- Environment: `requirements-workflows.txt` / `.lock` (`venv_workflows`).
- `runs.db`: one row per agent call (stage, agent, model, input hash, output
  path, verdict, time); local.
- `transcript.jsonl`: local.

Scripts live in `scripts/units_verification/`; each refuses to run while the
workflow's definitions (this folder, the scripts, `.claude/agents/`) have
uncommitted changes, and prints the inputs it used. Run from the project root:

```bash
P=venv_workflows/bin/python; S=scripts/units_verification; R=<run_id>
$P $S/locate_sources.py --run-id $R        # S0
$P $S/render_pages.py   --run-id $R        # S1   (--add SOURCE:PAGE for re-reads)
$P $S/marker_extract.py --run-id $R        # S1b  (about 6 min per source on CPU)
$P $S/calibrate.py      --run-id $R        # C0   (needs S0, S1, S1b)
$P $S/read_pages.py     --run-id $R        # S2   (--rereads after S3; --repeat 5 --out-dir reads_repeat)
$P $S/compare_units.py  --run-id $R        # S3   (--marker-only: Marker side alone)
$P $S/check_rows.py     --run-id $R --task S4    # S4, then --task S4m
$P $S/read_pages.py     --run-id $R --repeat 5 --out-dir reads_repeat --stage S2-repeat
$P $S/reader_stability.py --run-id $R      # reader stability: 5 pages read twice
$P $S/inspect_set.py    --run-id $R        # U1: inspection.json + local inspection/index.html
$P -m unittest discover -s $S/tests        # tests
```

Read schema (one object per symbol):
`{symbol, found, printed_as, meaning, value_as_printed, unit_as_printed, location}`.

## 7. Risks and failure modes

| Risk | Caught by |
|---|---|
| Text layer drops exponents and minus signs | agents get images only |
| Reader anchors on the claimed unit | blind read (no units or values given) |
| Wrong page in the table | Marker's symbol → page index; NOT_FOUND re-read on that page |
| Marker or the reader misreads Greek letters or subscripts | units are compared as dimensions, not strings; printed symbols after normalisation; S4 and S4m see the image |
| "Verified" rows that were never re-checked | every row passes through S2–S3 again |
| Two agents share a blind spot | S4 uses a different model; user inspects a random sample |
| Wrong meaning on a row whose unit agrees | S4m reviews every row's meaning against the image; S3 compares printed symbols |

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
2. **Value column:** resolved 2026-10-22: record the printed value too.
   S3 checks the read value against Marker (the table has no value column yet).
3. **Scope:** resolved 2026-10-22: literature symbols only (58 rows after removing j_T and j_v); the
   solver's parameters (`M0`, `alpha`) are out of this pass.
4. **Models:** resolved 2026-10-22: Sonnet reader, Opus adjudicator.
5. **Sample size** for U1: resolved 2026-10-22: 10 % of agreed rows.
