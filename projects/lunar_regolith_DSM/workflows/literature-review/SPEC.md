# SPEC: literature-review, version 2

The design for agents and scripts, agreed with the user on 2026-11-03 (`docs/planning/19`–`22`). The plain-language
version is `STORY.md`.

**Status: designed, not built.** `WORKFLOW.md` (v1) describes the pilot line as it is built and frozen; it governs
the nine-paper pilot run. This file replaces `WORKFLOW.md` when build step B10 (§13) lands, and `WORKFLOW.md` is then
archived. Until then, this file is the target and must not be cited as describing current behavior.

Contents:
- §1 purpose
- §2 rules
- §3 roles
- §4 stages
- §5 schemas
- §6 checks and calibration
- §7 stop and pause rules
- §8 library and index
- §9 runner, gates and ledgers
- §10 isolation and security
- §11 reproducibility and freeze
- §12 tests
- §13 build order
- §14 open and deferred
- Appendix: v1 to v2 stage codes

---

## 1. Purpose

- **Goal.** A library of literature facts the project can trust, and a checked memo per discipline. The library
  serves:
  - model choice (Pillar 2);
  - implementation checks (Pillar 3);
  - the manuscript.
- **Unit of work: a review.** One discipline, with its own brief. The order of disciplines:
  1. numerical methods;
  2. planetary observations and setting;
  3. experiments (sublimation; ice–sediment; ice sintering and metamorphism);
  4. sublimation physics;
  5. material properties (on demand: a property is looked up when the model needs it, not swept).
- **Products.**
  1. **Layer 1, the full index:** every paragraph, equation, table and caption of every collected paper, searchable
     by keyword and by meaning. No model is involved.
  2. **Layer 2, verified facts:** every fact of the core papers, and of any paragraph a later question needs, with
     conditions, quote, paragraph ID and independent checks. Only layer 2 facts may be cited.
  3. **A memo per discipline** answering the brief's questions; every sentence cites fact IDs.
- **Scope.** The review is one part of the project. It must never block model work: it runs in the background,
  with stops only at the user gates of §4.

## 2. Rules

1. **Decide or run.** Every task is D (needs judgment) or R (a script can do it). An R task is never given to a
   model.
2. **Justify every model step.** Each D step names its reason for a separate context. The reasons are:
   - **information:** it must see, or must not see, specific inputs;
   - **permissions:** it needs different rights;
   - **judgment:** it needs a distinct kind of judgment;
   - **breadth:** the work splits into independent parts.

   No reason means it merges into another call or becomes a script.
3. **No loop without need.** A D step is a **single call** (one request, script-chosen inputs, schema-checked
   answer, no tools) unless its next input depends on its own earlier findings. Only the scout (F4) qualifies.
4. **Code owns control flow.** The stage runner decides what runs next, assigns every item, pauses and stops. No
   agent assigns work or decides order.
5. **One writer per artifact.** Each output file has exactly one stage that writes it. Library records are
   immutable; a correction is a new version.
6. **Checks live in the harness.** A rule matters only if a script or hook enforces it. The prompt states rules; the
   harness checks them.
7. **The producer is not the only check.** Every fact is checked by a script and by a context that did not produce
   it.
8. **Untrusted text is data.** Paper and web text never acts as instructions (§10).
9. **Your seat only.** All model calls run through Claude Code on the user's enterprise login, which has a
   five-hour window and no weekly limit. Nothing is billed elsewhere, and no API key is used (§7.1).

## 3. Roles

| Role | Kind | Model | Reason (rule 2) | Loop | Rights |
|---|---|---|---|---|---|
| user | human | — | — | — | approves F1, F6, L2, M4, C1, the freeze; settles G2 mismatches |
| main session | interactive Claude Code | session default | designs and builds | yes | broad; **cannot start a non-dev run** (guard, §9) |
| `beth-foreman` | interactive Claude Code agent | `claude-opus-5-5` | permissions (no edits in execute); information (clean context) | yes | runs `run_*.py` (init, status, run, verify), `lr_status.py`, `lr_backlog.py`, `git status/log`; no edits, no approvals, no untrusted text |
| stage runner | script | — | — | — | starts stations through `wf_launch.py` sandbox profiles |
| `beth-searcher` | single call | `claude-opus-5-5` | judgment (search strategy) | no | none (no tools) |
| `morty-screener` | single call | `claude-sonnet-5-5` | judgment (relevance); breadth (per discipline) | no | none |
| `beth-scout` | **agent** | `claude-opus-5-5` | permissions (web); breadth | **yes** | WebSearch, WebFetch only; no Read/Write/Bash |
| `morty-reader` | single call | `claude-sonnet-5-5` | judgment (extraction) | no | none |
| `beth-checker` | single call | `claude-opus-5-5` | **information** (never sees the reader's reasoning); judgment | no | none |
| `morty-blind-reader` | single call | `claude-sonnet-5-5` | **information** (never sees the fact) | no | none |
| `beth-synthesizer` | single call | `claude-opus-5-5` | judgment (writing) | no | none |
| `beth-skeptic` | single call | `claude-opus-5-5` | information (sees memo and facts, not the drafting); judgment | no | none |
| `beth-recall-judge` | single call | `claude-opus-5-5` | judgment (does a fact match a gold statement) | no | none |

Names follow the house rule (`<character>-<purpose>`, prefix set by the model family). New definitions
(`beth-searcher`, `morty-screener`, `beth-scout`) go on `main` only after the user approves them. Models are pinned
in the model map (`MODELS`, `agent_call.py`), never as aliases. Revisions (R4) use `morty-reader` with the revise
task.

**Single-call mechanism (back end B).** `agentlib.call_agent()` runs `claude -p` with:
- `--tools ""` (only the StructuredOutput tool remains) and `--json-schema <schema>`;
- the isolation flags of §10;
- the inputs given **in the call**, not read by the model.

Inputs are given in the call through `--input-format stream-json` with image content blocks, *if* build test T-input
(§12) passes. Otherwise the model gets Read, limited to the call's clean folder, and its transcript is audited for
Reads outside that folder. The ledger records which mechanism ran.

**Back end A** (the Messages API with a Console key) is documented and not built: the user has no key. `anthropic`
2.4.0 is installed in `venv_workflows` for when that changes.

## 4. Stages

Per-unit:
- **once:** per review;
- **cand:** per candidate paper;
- **paper:** per included paper;
- **group:** per page group of at most 4 cited pages.

Gate:
- **B:** blocking; only that item and what depends on it wait;
- **R:** review; queued, and the line goes on;
- **—:** none.

| Code: name | D/R | Done by | Per | Input | Output | Checked by | Gate |
|---|---|---|---|---|---|---|---|
| **C1: calibrate checker** | D+R | `check_plants.py` + `beth-checker`; scored by `calibrate.py` | once per freeze | held-out plants and controls (§6) | score report | pass rules (§6) | B (user accepts) |
| **F1: approve discipline brief** | user | user | once | questions, seeds, scope, inclusion rules, recall list (about 20 must-find papers) | `brief.json` v*n* | schema | B |
| **F2: write queries** | D | `beth-searcher` | once | brief | `queries.json` | schema; queries are logged and shown at F6 | — |
| **F3: search and chain** | R | `search.py` | once, hop by hop | queries, seeds | `candidates` rows (source, query or hop, IDs) | dedupe by DOI, OpenAlex ID, normalized title; stop rule (§7.3) | — |
| **F3b: choose citations to follow** | D | `beth-searcher` | paper (per hop) | one included paper's reference list | `follow[]` with reasons | stop rule | — |
| **F4: scout the web** | D (loop) | `beth-scout` via `agentlib` | once | brief, current candidate list (titles, DOIs) | `proposals[]` (title, DOI or URL, reason) | each proposal goes through F3's dedupe and G2's identity check | — |
| **F5: screen relevance** | D | `morty-screener` | cand | title, abstract, brief | `{decision: include/exclude/other_discipline, discipline, reason, core_rank}` | coverage ledger: every candidate has exactly one decision | — |
| **F6: approve shortlist** | user | user | once | shortlist; a random sample of 20 rejects; queries and their hits; recall test result | approved list; the core set (seeds + top-ranked) | recall ≥ the brief's target, or the user accepts | B |
| **G1: fetch** | R | `fetch.py` | paper | DOI | PDF (+ LaTeX or XML source where open) | order: Zotero (`zotero_pull.py`), Unpaywall/OpenAlex open access, arXiv only if nowhere else; else the user's fetch list | R (list) |
| **G2: check copy** | R | `check_copy.py` | paper | PDF, DOI | identity record | Crossref title and authors on page 1; retraction status; hash | B on mismatch |
| **G3: prepare paper** | R | `prepare.py` (MinerU, `venv_docparse`, offline sandbox) | paper | PDF | page images (300 dpi), text, LaTeX equations with printed numbers, `layout.json`, **paragraph IDs** | coverage of ink by layout blocks; equation-number pass | — |
| **G4: index paragraphs** | R | `index.py` | paper | G3 text with paragraph IDs | layer 1 index (keyword: SQLite FTS5; meaning: local embedding model, chosen by T-embed) | rebuild gives identical hashes | — |
| **R1: extract facts** | D | `morty-reader` | group (core papers) | page images + numbered paragraphs of the group; brief (for tags only) | `facts[]` (§5) | schema | — |
| **R2: quote check** | R | `quote_check.py` | paper | facts, G3 text | per-fact pass or fail | quote ≥ 8 words, found **inside its paragraph ID**, unique; the same for each condition's quote | R (failures) |
| **R3a: assumption list** | D | `beth-checker` | paper | the whole paper (images + paragraphs), no facts | `assumptions[]` with paragraph IDs | quotes checked as in R2 | — |
| **R3: check facts** | D | `beth-checker` | group | the group's facts (statements, quotes, conditions; no reader reasoning), the group's pages, the assumption list | per fact: verdict; per condition: label (§6.2); missing conditions with quote and paragraph ID | missing-condition quotes checked as in R2 | — |
| **R4: revise facts** | D | `morty-reader` (revise task) | fact | the failed fact, the checker's report, cited and condition pages | revised fact or withdrawal | back through R2 and R3; at most 2 rounds (§7.2) | R (after 2 rounds) |
| **R5: blind reread** | D+R | `morty-blind-reader` + `compare_reads.py` | paper | page images + item pointers only ("equation 13, page 4"; "the value of α_c, page 6") | transcriptions | compared by script with the fact (`compare_reads.py`: canonical LaTeX; values normalized) | R (disagreements) |
| **L1: publish** | R | `curator.py` | paper | facts passing R2, R3 and R5 | layer 2 records, v1 | schema; one writer | — |
| **L2: your fact review** | user | `review.py` (renamed `u3_review.py`) in the user's terminal | once | disputed facts, R5 disagreements, R4 failures after 2 rounds; each with paper, page, **paragraph**, quote, conditions, checker report | decisions right/wrong/unsure (with `--redo`) | — | B |
| **L3: publish your decisions** | R | `curator.py --after-review` | paper | L2 decisions | new versions; rejected facts marked `rejected_by_user`, never deleted | — | — |
| **M1: draft memo** | D | `beth-synthesizer` | once | layer 2 facts for the brief's questions (JSON) | memo with fact-ID citations | — | — |
| **M2: citation check** | R | `check_citations.py` | once | memo, facts | pass or fail per sentence | every sentence cites (except marked options, gaps, *library silent*); numbers and units match the cited fact | — |
| **M3: skeptic review** | D | `beth-skeptic` | once | memo, cited facts | objections[] | revision then M2 again; leftovers go to the user | R |
| **M4: approve memo** | user | user | once | memo, M2 and M3 reports | approval | — | B |
| **K1: backup** | R | `backup.py` + `restore_test.py` | once | library, memo | dated archive + hash manifest | restore into an empty folder: 0 mismatches; M2 on the restored memo passes | — |

Order: C1 → F1 → F2 → F3 (with F3b per hop) ∥ F4 → F5 → F6 → per paper G1 → G2 → G3 → G4 → (core papers) R1 → R2 →
R3a → R3 → R4 → R5 → L1 → M1 → L2 → L3 → M2 → M3 → M4 → K1. M1 runs before L2 so that the user reviews the facts the
memo actually uses. An L2 decision that adds or removes a cited fact sends the memo back to M1 before M2 (closes the
v1 backlog item of 2026-11-02).

**A new question later:** search layer 1, then run R1–L3 on the paragraphs it finds that layer 2 doesn't cover,
then M1–M4 for that question.

## 5. Schemas

Every D step's answer is checked against a JSON Schema (`schemas/<name>.json`, `additionalProperties: false`, all
fields required). The schema guarantees shape, not truth; truth is §6. Main schemas:

- **brief:**
  - `review`, `discipline`, `questions[]` (id, text, model_critical)
  - `seeds[]` (DOI), `scope`, `include_rules[]`, `exclude_rules[]`
  - `recall_list[]` (DOI), `recall_target`
  - `version`
- **candidate** (a ledger row, not a model answer):
  - `doi`, `openalex_id`, `title_norm`
  - `found_by[]` (query or hop or scout)
  - `decision`, `decided_by`, `reason`
- **screen:**
  - `decision` (include | exclude | other_discipline), `discipline`, `reason`
  - `core_rank` (1–5)
- **fact:**
  - `id`, `paper`, `statement`
  - `kind` (equation | value | relation | method | observation)
  - `basis` (this_paper | cited_work, with the cited reference)
  - `quote`, `paragraph` (`p<page>.<n>`)
  - `equation` (`latex`, `printed_number`) or null
  - `value` (`quantity`, `value`, `unit`) or null
  - `conditions[]` (`text`, `quote`, `paragraph`)
  - `tags` (`discipline`, `questions[]`)
- **check** (R3), per fact:
  - `verdict` (supported | not_supported | partly)
  - `conditions[]` (`label`: necessary | scope | unnecessary, `reason`, `quote`, `paragraph` for an unnecessary
    label)
  - `missing_conditions[]` (`text`, `quote`, `paragraph`)
  - `meaning` (`established`, `by_whom`, `assumptions`, `evidence`), for relation and method facts
- **blind reading:** `items[]` (`pointer`, `latex` or `value` + `unit`, `not_found`)
- **memo:** Markdown with inline tags `[fact_id]` and, for numbers, `[fact_id quantity]`

Paragraph IDs: `p7.3` means page 7 (PDF page), text block 3 in reading order from `layout.json`. They are displayed
as "page 7, paragraph 3".

## 6. Checks and calibration

### 6.1 Checks per output

| Output | Check | Kind |
|---|---|---|
| copy | identity, retraction, hash (G2) | script |
| quote | inside its paragraph, ≥ 8 words, unique (R2) | script |
| fact | R3 verdict from a context that never saw the reader's reasoning | model, different input |
| conditions | R3 labels against the paper-wide assumption list; missing ones quote-checked | model + script |
| equation, value | R5 blind transcription compared by script | model + script |
| memo | M2 citations, numbers, units | script |
| memo | M3 skeptic | model |
| records | rehash; restore test (K1) | script |
| search | recall test against the brief's list; coverage ledger complete | script |

### 6.2 The condition rule

A **condition** is a statement that restricts when the paper says the fact holds: geometry, regime, parameter range,
material, model assumption, approximation. The checker labels each one:
- **necessary:** the paper states or shows the fact fails, or changes, without it;
- **scope:** the paper asserts the fact only under it and is silent beyond it;
- **unnecessary:** the paper itself states or shows the fact holds without it. It must quote the passage.

Necessary and scope conditions are kept. Unnecessary ones are flagged, not deleted. A fact is **complete** when
every necessary or scope condition the paper attaches to it is listed, wherever the paper states it: the same
paragraph, an earlier section's assumptions, a figure caption, a table header. The assumption list (R3a) is the
reference for completeness. A missing condition is reported with a quote and paragraph ID and confirmed by script.

### 6.3 Calibration (C1)

Carried over from v1 (`WORKFLOW.md` §5), with these changes:
- the condition labels are scored;
- the material is relabeled after the user rules on the paused control flags.

Material:
- the user's 40 gold statements (controls);
- the user's 11 planted errors;
- the approved script plants and the six dropped-condition plants;
- a tuning split and a held-out split (`calibrate.py split`).

Pass rules, on the one frozen held-out evaluation:
- every E3-type plant (fabricated quote, wrong page or paragraph) is caught;
- every sign or exponent plant is caught;
- **at most 10 % of controls flagged**;
- the catch rate on plants reaching R3 is reported with its one-sided 95 % Clopper–Pearson lower bound. It is a gate
  (≥ 85 %) only with 19 or more such plants.

An abstention counts against the checker. User plants and script plants are reported separately. The script-only
baseline (`baseline.py`) is reported alongside.

C1 is rerun after any change to an agent file, prompt, schema, model ID, parser version, **engine (back end or
single-call mechanism)** or Claude Code version. The effort level is Claude Code's default, recorded in the
manifest; it changes only as a diagnostic (MEMORY.md, 2026-11-01).

### 6.4 Measurements (reported, not gates)

- **Recall:** gold statements found by a fresh R1 run.
- **Repeatability:** R1–R3 run twice on one paper, compared fact by fact.
- **Hand check:** the user checks 30 random layer 2 facts against the PDFs after each discipline's 40-paper pilot.
- **Baselines:**
  - the v1 pilot (same nine papers, same calibration sets);
  - a same-effort single-agent run on the first discipline, against the team's result.

## 7. Stop and pause rules

### 7.1 Usage cap

Every `claude -p` call made with `--output-format stream-json` emits `rate_limit_event` messages. Fields:
`rateLimitType` (`five_hour`), `utilization` (0–1), `resetsAt`, `status`. This was seen in the v1 transcripts on
2026-11-03; it is typed in the Agent SDK and undocumented for `-p`. `agentlib` records each event in the ledger.

- **At utilization ≥ 0.90** (a setting, recorded in the manifest), the runner starts no new item and lets running
  items finish. It writes `pause_<time>.json` and prints the pause report:
  - what ran since the last pause;
  - every item's state;
  - what remains;
  - `resetsAt`.
- **Default mode:** stop, with exit code PAUSED.
- **`--continuous` mode:** sleep until `resetsAt` + 5 min, then resume. Each pause is still reported.
- **If no utilization is reported:** a per-run call cap (manifest setting) applies.
- **If a call fails with a limit message** ("You've hit your session limit"): the item is recorded as *not done*
  (not failed), and the runner pauses.
- There is no weekly limit on this account, so none is monitored.

### 7.2 Repair rounds

A failed fact gets at most 2 revise-and-recheck rounds (R4), then goes to L2. A failed call is retried once, then
the item fails and the runner stops with the log path.

### 7.3 Search stops

- **Citation chaining stops** when a hop adds < 5 % new included papers (setting), or after 3 hops.
- **The search passes** when it finds the recall list at or above the brief's target. Otherwise the user decides at
  F6: add queries and rerun F3, or accept.

### 7.4 Scale stops

Each discipline starts with a pilot of about 40 papers and the hand check (§6.4) before it grows.

## 8. Library and index

- **Truth is immutable files**, one per object: paper content (G3 output, written once), facts (versioned), user
  decisions, memos. Catalog and index are built from them and can be rebuilt.
- **Layer 1:** `index/` holds SQLite FTS5 over paragraphs, plus embeddings from a local model, offline and pinned.
  Every entry carries paper, paragraph ID and copy hash.
- **Layer 2:** `library/facts/` holds one file per fact version, with `events.jsonl` as the append-only log.
- **Rejected facts** stay, with status `rejected_by_user`.
- **A normal session never reads `library/`, `index/` or `untrusted/`.** Scripts print counts and IDs only. Evidence
  is read in review mode.

## 9. Runner, gates and ledgers

- **The runner** (`run_b1.py`, to be renamed `run_review.py` with `--review <name>`) walks the stage graph in
  `flowchart.json` `run` blocks:
  - `per`, `requires`, `station`, `gate`, `steps`, `done`;
  - it runs each step through `wf_launch.py` (sandbox profile, allow-listed script);
  - it checks the `done` files, never a script's own report;
  - it records every item in `runs.db`.

  Commands: `init`, `status`, `run [--parallel N] [--continuous] [--until X] [--max N]`, `approve`, `freeze`,
  `thaw`, `verify`.
- **Parallelism:** per-paper and per-group items run up to N at a time (default 4). G3 may run ahead of R stages.
- **Gates:** a blocking gate holds its item and dependents only. Approvals are recorded only by `approve`, typed by
  the user.
- **Ledgers:**
  - `runs.db`: line items, calls and usage events;
  - `calls.jsonl`: every model call with its input hash, model, mechanism, usage and stop reason;
  - the coverage ledger (`candidates` table);
  - the review queue.
- **Foreman guard (hook H8, to build; needs user approval to add under `.claude/`):**
  - `foreman.py` sets `LR_ROLE=foreman` in the foreman session's environment;
  - a PreToolUse Bash hook blocks `run … run|init|approve` on a non-dev run when `LR_ROLE` is not `foreman`, and
    names the command to start the foreman;
  - dev runs stay allowed in the main session during develop;
  - the session-start message names the session's role.
- **Status** comes from `lr_status.py`, generated from `runs.db` and the output files, never written by hand.

## 10. Isolation and security

- **Worker isolation:**
  - flags: `--setting-sources "" --strict-mcp-config --no-session-persistence --disable-slash-commands
    --permission-mode dontAsk --system-prompt <body>`;
  - environment: `CLAUDE_CODE_DISABLE_AUTO_MEMORY=1` (since 2026-11-03);
  - a clean temporary working folder holding only the call's inputs.

  The canary test of 2026-11-03 (`WORKFLOW.md` §12) showed these keep CLAUDE.md, CLAUDE.local.md and auto memory
  out. `--safe-mode` is an equivalent, tested alternative that keeps the login and permission files. `--bare` is
  not usable, because it needs an API key.
- **Single calls have no tools**, so injected text can at most corrupt an answer, which R2, R3 and R5 catch.
- **The scout** has WebSearch and WebFetch only: no Read, Write or Bash. It runs in the `anthropic` sandbox profile,
  and its proposals are data checked by F3/G2. Web text is untrusted.
- **Stations** run under `wf_launch.py` sandbox profiles (offline for G3, network only where needed). The allow-list
  is edited only by the user.
- **Hooks H1–H7 carry over** (phase guard, line guard, foreman Bash guard, stale-status stop). H8 is the foreman
  guard.

## 11. Reproducibility and freeze

- **The manifest per run** records:
  - git commit, freeze version;
  - Claude Code version, model IDs, effort;
  - back end and single-call mechanism;
  - prompt, schema and input hashes;
  - parser lock and model snapshot;
  - embedding model;
  - the usage-cap setting;
  - check verdicts.
- **Freeze** (`/lr-freeze`, user approval) hashes:
  - agent files, prompts, schemas;
  - station scripts, the model map;
  - parser and index environments;
  - `requirements-workflows.lock`.

  In execute, the runner refuses changed frozen files, and a hook asks before an edit. Ideas go to `BACKLOG.md`.
- **The only way back:** a real error found by L2, R5 or M3 stops the run. The phase returns to develop, the fix is
  made, C1 is rerun, and the user approves a new freeze.

## 12. Tests at build time

Each test prints PASS or FAIL from files and transcripts, never from an agent's report. Results go into §12 of the
active design file with dates.

| Test | What | Pass |
|---|---|---|
| T-isolation | canary words in the CLAUDE.md files and every auto memory folder; a worker reports what it sees; control without isolation | the control sees them; the worker sees none (passed 2026-11-03 for v1 flags + env) |
| T-input | single call with page images given in the call (`--input-format stream-json`), no Read tool | the model transcribes a planted value from an image correctly; no tool call in the transcript |
| T-cap | runner against a fake `claude` reporting utilization 0.5, 0.91 | starts no new item after 0.91; running items finish; pause report written; `--continuous` resumes after `resetsAt` |
| T-guard | main session tries `run` on a non-dev run; the foreman tries the same | blocked; allowed |
| T-paragraph | paragraph IDs on the three calibration papers; quotes moved to a neighboring paragraph | IDs stable across reruns; moved quotes fail R2 |
| T-recall | F3 on numerical methods with the brief's recall list | target met, or the user accepts at F6 |
| T-embed | candidate local embedding models on held-out query → paragraph pairs written by the user | the chosen model's recall at 10 reported; offline; pinned |
| T-coverage | inject a candidate with no decision | coverage script reports it |
| T-pilot | v2 on the nine pilot papers and the pilot's calibration sets | C1 passes; layer 2 facts cover the pilot's published claims or each gap is explained; per-paper time and usage reported |
| T-single | same-effort single agent on the first discipline's pilot | reported next to the team's result |
| T-repeat, T-restore, T-applicability, T-trifecta | as in v1 (`WORKFLOW.md` §12) | as in v1 |

## 13. Build order

Before B1–B10, the **pilot** (v1 line, `WORKFLOW.md`) runs to completion:
1. the c006 redo (user);
2. the 9 control rulings (user);
3. the gold fixes or prompt change, then the tuning rerun;
4. the freeze (user);
5. the live hook checks;
6. the nine-paper run;
7. the methods memo.

The target is about 8 November. Model development (Pillar 2) starts from that memo.

The v2 build runs alongside model work:
- **B1** call module for single calls in `agentlib` (mechanism from T-input); usage events in the ledger; T-input.
- **B2** usage cap and pause report in the runner, `--parallel`, `--continuous`; T-cap.
- **B3** paragraph IDs in `prepare.py`; R2 per paragraph; T-paragraph.
- **B4** foreman guard H8 and the session role message (user approves the `.claude/` change); T-guard.
- **B5** stage codes and names in `flowchart.json`, the runner and the status display (v2 codes; Appendix).
- **B6** R stages as single calls: reader, checker (R3a and R3 with condition labels), blind reader, reviser; the
  schemas of §5; C1 rerun.
- **B7** Get: `fetch.py` (G1), G2 extended to retraction status; `index.py` with FTS5 + embedding (G4); T-embed.
- **B8** Find:
  - `search.py` (F3, OpenAlex and Semantic Scholar, chaining and stop rule);
  - `beth-searcher`, `morty-screener`, `beth-scout` definitions (user approves);
  - the coverage ledger;
  - the F6 report;
  - T-recall, T-coverage.
- **B9** T-pilot: v2 on the nine pilot papers, compared with the pilot.
- **B10** this file replaces `WORKFLOW.md` (archived); `TEMPLATE.md`, the contract format and the foreman's
  definition are updated; freeze v2. Then numerical methods: a 40-paper pilot, the hand check, T-single, scale-up.

## 14. Open and deferred

- **Back end A** (Messages API) is unbuilt and needs a key the user doesn't have.
- **Agent SDK:** not used, by user decision on 2026-11-03. `agentlib` serves every worker, the scout included.
- **Lower effort for screening** is not adopted. The default-effort rule stands.
- **LangGraph:** not for this pipeline; possibly a side project.
- **Units workflow:** Marker was removed on 2026-11-03, so that workflow needs a parser before it is used again.
- **Material properties** are on demand; their brief format is to be decided when first needed.

## Appendix: v1 to v2 stage codes

| v1 | v2 | | v1 | v2 |
|---|---|---|---|---|
| U0 | F1: approve discipline brief | | E5 | R4: revise facts |
| — | F2–F6 (new) | | X1 | R5: blind reread |
| A2 | G2: check copy | | L1 | L1: publish |
| P1 | G3: prepare paper | | U3a | L2: your fact review |
| — | G1, G4 (new) | | L1b | L3: publish your decisions |
| E1 | R1: extract facts | | Y1, Y2, Y3, U4 | M1–M4 |
| E3 | R2: quote check | | BK | K1: backup |
| E4 | R3a + R3: check facts | | C0 | C1: calibrate checker |

The reviews are named by discipline: "methods review" (was B1), "planetary review", and so on.
