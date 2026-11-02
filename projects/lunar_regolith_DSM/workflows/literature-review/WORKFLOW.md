# WORKFLOW: literature-review

Status: **revision 10** (2026-10-29). Revision 10 applies the Codex review of
revision 9 (§11, round 5), which the user accepted in full on 2026-10-29 (§9,
items 40 to 43). The plan of record is now a small path for the first branch,
B1, limited in time. Sections 1 to 8 describe that path. Everything it does
not need is kept in §13 (Deferred), not deleted, and comes back when a later
branch needs it.

What revision 10 changes, and why:
- **Scope.** B1 is a bounded review of nine fixed papers, run to answer one
  modeling question (§9, item 41). The build budget is about 3 to 5 working
  days; after that the user returns to model work (Campaign B and Pillar 2 in
  `project/PILLARS.md`). Discovery, screening, acquisition by script, LaTeX
  and JATS text, automated tables and figures, the passage index, the
  reference desk and the content-addressed archive are deferred (§13). The
  reason is that nine known papers need none of them, and every hour on the
  pipeline is an hour not spent on the model.
- **Checks matched to the kind of claim.** A blind transcription can show
  that an equation or a value was read right, but not that a sentence about
  conditions, attribution or applicability is right. So the independent blind
  read (X1) is used only for equations and values. Prose claims get a meaning
  check by the checker, and the user's approval when they support a model
  option.
- **The checker looks for conditions first.** Before it sees any claim, the
  checker writes an applicability checklist from the paper's model and method
  pages. A condition the reader missed, on a page the reader never cited, then
  still reaches the checker.
- **Fewer roles.** One reader per paper, a fresh checker, no routine second
  reader (E2), and one skeptical memo checker instead of two.
- **Calibration.** The material the user has written is used as it is; no new
  planted errors are written for B1. The per-kind gate that could not be met
  below 19 items is removed. There is one frozen final evaluation instead of a
  fresh held-out set after every failure.
- **Records.** P1's extracted content stays immutable. The reader's symbol
  lists and inferred units are separate analysis records. A memo is made in a
  fixed order: checked claims, then a draft memo, then verification of the
  claims it cites, then the approved memo.
- **Security.** The readers and checkers stay isolated (no network, no
  writes). For B1, the counts-only main session and the reference desk are
  replaced by a dedicated review session that may read selected evidence but
  has no network, no writes and no broad file access (§4). The
  planted-instruction test now also checks the scientific claim that comes
  out, not only that no forbidden tool was called.
- **Reproducibility.** For the course, the same prompts given to the agents
  must give similar results. A repeat run on one paper, compared claim by
  claim, is the evidence (§8). One restore test replaces the archive system.

**Built** (2026-10-29; 59 unit tests pass):
- `check_calibration.py`, `make_plants.py` and `calibrate.py` (validation of
  the user's calibration files, script plants, the tuning and held-out split);
- `zotero_pull.py`, used live on 2026-10-29 to copy in the six B1 seeds not
  already on disk;
- `prepare.py`, first part (P1): 300 dpi page images, tiles, page text and
  page layout from a PDF parser, displayed equations with their printed
  numbers, headings, the PDF-to-printed page map and a manifest. It runs
  offline under a macOS sandbox profile;
- `check_permissions.py` and the permission files.

The PDF parser is being chosen. Marker was used for the first runs and is
being replaced by a parser that scores higher on document-parsing benchmarks
(§3, P1). The user has written the calibration material: 40 gold statements
on 4 papers, 11 planted errors, 57 approved script plants and the T-res pages.

**Lessons from the first runs** (2026-10-29; `project/RETROSPECTIVE.md`):
- **Counts are not checks.** The first tiler covered all the ink and failed no
  page, yet it cut a table mid-row, separated a plot from its axis labels and
  split one caption over three tiles. The user found this by looking.
- **A count needs a known answer.** P1 reported 47 numbered equations for Li
  et al.; the user counted 58. The parser, not the extractor, missed a form.

**History.** Revision 5 was approved for building on 2026-10-23. Revision 6
added branches, the passage index, structured tables, the units check, X1,
the scout's challenge task, a screening baseline, the split calibration and
the permission table. Revision 7 applied critique round 3 and the user's
decisions of 2026-10-27. Revision 8 applied critique round 4 (review mode by
allow-list, a script-only launcher, a blind first transcription by the
checker, calibration gates per kind, one archive) and stored the whole paper
in each record. Revision 9 (2026-10-28) dropped Codex from the pipeline, wrote
the permissions as `settings.json` files, took up Zotero again and chose the
B1 seeds.

Flowchart: `flowchart.html`, rendered from `flowchart.json` by
`scripts/flowchart/render_flowchart.py`; it still shows revision 9.
B1 run contract: `B1_CONTRACT.md`. Decisions: §9. Deferred: §13.

## 0. How B1 runs: the line, phases and freeze

Added 2026-11-01 (user decisions of that day; plan in the session of
2026-11-01). **Planned, not yet built**, except where marked. It follows
`workflows/TEMPLATE.md` ("The line", §7). The reason: in the development run
`b1-dev2` the main session ran each stage's driver by hand, so stages were
skipped (V1, X1) or run out of order (the memo before the frozen calibration),
and prompts kept changing during the run.

**Who reports to whom.** Scripts manage; agents work.
- **The user** approves the gates below, the freeze and the memo.
- **`beth-foreman`** (`agents/beth-foreman.md`, started with
  `venv_workflows/bin/python scripts/literature_review/foreman.py`) runs the line, reports each station's result
  and the items waiting for the user, and stops at blocking gates. It may run
  only `run_b1.py`, `lr_status.py`, `lr_backlog.py`, `git status` and
  `git log`; it never edits a definition and never reads untrusted text.
- **`run_b1.py`**, the line script, walks the stage graph in `flowchart.json`
  (each node's `run` block: station, upstream stages, per-item unit, gate,
  outputs, done test). It refuses a stage whose upstream stages are not done,
  and records every call in the run's `runs.db`.
- **Station scripts** (`read_check.py`, `synthesize.py`, `skeptic.py`, and
  the ones still to build) write each work order from a fixed template,
  validate the output against its schema, audit the transcript, retry once,
  and write `report.json` (counts, failures, gate items).
- **Workers** are the headless agents of §4, defined in
  `workflows/literature-review/agents/<agent>.md` (Claude Code subagent
  format; the body is the system prompt and has a "Your place in the
  line" section). `prompts/` holds only the per-call work orders.

**Gates.**
- *Blocking* (the line stops; the user records the decision with
  `run_b1.py approve`, never the foreman): A2 copy decisions; the freeze; the
  C0 result; V1 disagreements; U3 approval of each prose claim a model option
  rests on; U4.
- *Review* (queued in `review_queue.json`; the line goes on with other
  papers): disputed claims, E3 quote failures, curator refusals.

**Phases.**
- **Develop** (now). Anything may change. Runs are tagged `dev` and are not B1
  results; `b1-dev2` is one.
- **Freeze.** The user approves `/lr-freeze`, which writes
  `freeze/b1-vN.json` (hashes of the agent files, prompts, schemas, station
  scripts, the model map and the parser: `requirements-docparse.lock`, the
  packages installed in `venv_docparse` and the MinerU model snapshot) and the git tag
  `lr-b1-freeze-vN`.
- **Execute.** `run_b1.py` refuses to start if a frozen hash differs; a hook
  asks before any edit to a frozen file. Ideas go to `BACKLOG.md` and wait.
- **The one way back** is §5's repair rule: a real error found in U3, X1 or Y3
  stops the run, the phase returns to develop, and a new freeze is approved.

**Status** is computed by `lr_status.py` from `runs.db` and the output files
and written into the contract's Status table; it is not written by hand.

## 1. Purpose

- **Question.** A reusable workflow for literature reviews of any size. Each
  review answers questions written down before reading starts (the *brief*).
  - **First review (Pillar 1, branch B1):** numerical methods for multiphase
    systems in which some phases are stationary and others evolve, starting
    from diffuse-interface and diffuse-domain methods, and methods that couple
    a sharp-interface description with a diffuse-interface one.
  - **The question the B1 memo must answer** (§9, item 41): how to couple the
    stationary grains with the evolving ice–vapor interface without a
    spurious air film. Surface diffusion and other vapor-diffusion paths are
    secondary questions.
  - **Motivating failure:** the spurious air film that forms where the
    evolving ice field relaxes against a stationary sediment field, so that
    φ_i + φ_s ≠ 1 in the band.
- **Finished output of B1.**
  1. **A record for each of the nine papers**, in three layers:
     - **Content:** the paper as P1 extracted it, with no model rewriting:
       the copy itself, page images, the text of every page, every displayed
       equation in LaTeX with its printed number, the page layout, the
       section map and the PDF-to-printed page map. Written once and never
       edited.
     - **Analysis records:** model output about the paper that is kept but
       not checked: the reader's symbol list, the checker's applicability
       checklist, and inferred units. Each is versioned and written once.
     - **Claims:** statements about the paper that have passed the checks.
       Each has its kind, basis, a verbatim quote, the page, any equation in
       LaTeX with its printed number, and verdicts from a checker that did
       not write it. Only claims are cited in memos.
  2. **One memo** answering the B1 question: every statement cites claim IDs;
     it ends with the model options it supports, their conditions, and what
     the nine papers do not answer. It states that it rests on nine fixed
     papers and makes no claim to cover the field.
  3. **Model options**, at the end of the memo, for discussion with the user
     (`project/PILLARS.md`, Pillar 2).
- **Acceptance test for B1.** All of these must hold:
  - the frozen final calibration evaluation meets its pass rules (§5);
  - every published claim has a quote found on its cited page by script;
  - the checker has confirmed each claim against the page images, after
    writing its applicability checklist for those pages;
  - every equation or value the memo cites has passed the independent read
    (X1), or the user has decided it;
  - every prose claim the memo uses to support a model option has the user's
    approval of four things: what was established, by whom (this paper or a
    work it cites), under which assumptions, and with what evidence;
  - every memo sentence cites a verified claim (script check), except
    sentences marked as model options, *library silent* or gaps; the
    skeptical memo checker finds no unsupported statement, or the user has
    accepted each remaining one and the memo marks it;
  - the user has approved the memo (U4).
- **Stopping rule.** B1 stops when all nine papers have been read and checked
  and the memo names the supported options, their conditions and the
  unresolved questions, or when the build budget runs out, whichever is first.
  In the second case the remaining automation is cut, not the checks.
- **Baselines.** The practice-library recheck and the scripted screening and
  search baselines are deferred with the stages they measure (§13).

## 2. Library branches, seeds and waves

**Branches** are the library's main subject areas. They are top-level tags
(facet *branch* in `vocab.yaml`); a paper can belong to several. Proposed
starting set (the user decides, §9):

| Branch | Covers |
|---|---|
| B1 numerical methods | diffuse interface, diffuse domain, multiphase field and smoothed boundary methods with stationary phases; methods that couple a sharp-interface description with a diffuse-interface one (Pillar 1) |
| B2 sublimation physics | Gibbs–Thomson, Hertz–Knudsen, kinetic coefficients, vapor diffusion and Knudsen transport |
| B3 ice sintering and metamorphism | experiments and models of snow and ice microstructure |
| B4 planetary setting | lunar polar cold traps, ice in regolith, thermal environment |
| B5 properties | ice, vapor and regolith properties the model uses |

**Seeds.** For each branch the user picks 3 to 7 papers, at least one a
review, and at least one or two papers the user knows well. B1's seven seeds
are in §9, item 38. Seeds are read first. In later branches they also steer
discovery and screening (§13); in B1 they are the whole set.

**B1 is a fixed-paper review.** It reads nine papers: the seven seeds, Kaempfer
& Plapp 2009 and Moure & Fu 2024 (the published copy since contract v2,
2026-11-01; the arXiv v1 LaTeX source stays on disk). There is no discovery, screening or tier. The run contract,
`B1_CONTRACT.md`, lists each paper with its copy's hash, its role (core or
context), the questions, the outputs, the stages that run, and who approves
what. A citation found while reading that looks essential is listed in the
memo as a lead, not read in B1. Backward and forward citation chaining from
these papers (for example from Blackford's references) is a later step.

**Order of work.**
1. Build what the B1 path needs (§3) and its tests (§12).
2. Run the calibration diagnostic during development, then the frozen final
   evaluation (§5).
3. Run B1 on the nine papers.
4. Draft, verify and approve the memo.
5. Return to model work.

## 3. Stages

The B1 path. Stages of the full design that B1 does not run are in §13, with
their original text.

| # | Stage | Agent or script | Input | Output | Check | Stop / escalate |
|---|---|---|---|---|---|---|
| U0 | Brief and contract | user, with the main session | `B1_CONTRACT.md` | `reviews/b1/brief.v1.json`: the memo question and the secondary questions, each marked *model-critical* or not; the paper list with copy hashes; the outputs | user approves | — |
| A2 | Check the copy | `check_copy.py` | each copy in the contract | copy kind (publisher, accepted manuscript, arXiv); DOI or title on page 1 | match | a mismatch, or a copy that is not the version of record → the user decides which copy is read |
| P1 | Prepare | `prepare.py`, offline | each copy | content (§6): page images, page text and layout from the PDF parser, equations with printed numbers, headings, page map, manifest; tiles or zoom crops as T-res decides | every page has an image and text; equation numbers checked against a known count where one exists | a failed page is logged and cannot be cited |
| E1 | Read | `morty-reader`, fresh per paper (per section above about 40 pages) | the brief's questions; the page text; page images, and tiles or zoom crops for every equation, number and symbol it transcribes | `reads/<paper>.json`: claims (§7), each condition with its page; a symbol list (an analysis record) | schema; one output per paper or section | failure rules below |
| E3 | Check extraction | `check_extraction.py` | read, page text | quote checks; page coverage | quote ≥ 8 words or a whole sentence, unique in the paper, and found on the cited page (the page text is the PDF parser's, or the text layer used only to locate the quote) | failure → E4, marked; a quote that cannot be placed → its own flag |
| E4 | Check | `beth-checker`, fresh, per group of 3 to 5 pages, in three steps | **step 0, applicability:** the paper's model and method pages and the brief's questions, no claims; **step 1, blind:** for each equation and value claim, only the quantity or equation name, page and position; **step 2:** the claims citing those pages (their fields, never the reader's reasoning or the parser's text), the step 0 checklist, E3 flags; in every step the page images, and tiles or crops | step 0: a checklist of assumptions, boundary conditions, which fields are stationary and which evolve, parameter limits and the validation regime, each with its page; step 1: the equation or value and unit as printed; step 2, per claim: supported / overstated / unsupported / not on the page; for prose claims the four meaning fields (§7); dropped conditions, checked against the checklist | schema; each step returned before the next is sent, confirmed by the transcript audit | not fully confirmed → E5 |
| E5 | Revise (≤ 2 rounds) | `morty-reader`, fresh | the claim, the checker's reasons, the cited and condition pages | revised claim → E3, then a fresh E4 | as E3, E4 | still failing → U3, *disputed* |
| L1 | Publish | `curator.py`, the only writer | claims with passing verdicts or a user decision; analysis records | claim files, analysis records, events (§6) | refuses items without citation, quote match and verdict | refusal → logged, shown at U4 |
| Y1 | Draft memo | `beth-synthesizer` | the question; every verified claim for it, and the applicability checklists of the papers those claims come from | `memos/b1/Q1.v1.md`, a draft: statements citing `[id@v1]`; every number tagged with its claim and quantity (§5); model options, *library silent* and gaps marked | Y2 | — |
| V1 | Verify the claims the draft cites | `morty-blind-reader` and `compare_reads.py` (X1, equations and values); the user, in the review session (U3, prose claims) | each claim the draft cites | X1: the equation or value as printed, read blind; U3: the user's approval of the four meaning fields of each prose claim used for a model option | `compare_reads.py` compares the reader, the checker's step 1 and the blind reader, plus the PDF parser's LaTeX and the authors' source where they exist | a disagreement in content → the user; a claim the user rejects is withdrawn and Y1 redrafts |
| Y2 | Check citations | `check_citations.py` | draft memo, claims | per sentence: cites ≥ 1 verified claim (model options, *library silent* and gaps exempt); every number's tag names a cited claim and matches its parameter field (§5) | script | failure → Y1 |
| Y3 | Check memo | `beth-skeptic`, fresh, two calls | **call 1, support:** memo and the cited claims; **call 2, omissions:** memo, every verified claim on the nine papers and their applicability checklists | call 1: per statement supported / overstated / unsupported, told to find the strongest reason each statement is wrong; call 2: omitted or contradicting claims, conditions the memo leaves out, gaps | schema | not supported → Y1 revises once, the rest to U4 |
| U3 | User review | user, in the review session (§4) | every prose claim the memo uses for a model option; every disputed claim; key pages of the core papers, with their claims shown after the page is read | right / wrong / unsure; *missed* statements with their page | census of the claims shown | wrong or unsure → E5; a missed statement becomes a claim through E3 and E4; an error found here → the prompt is repaired and the affected claim family inspected (§5) |
| U4 | Approve | user | memo, Y3 report, gaps, refusals | approved memo, or a new brief version | — | — |

The memo is made in this order: checked claims (E1 to E4), a draft (Y1),
verification of the claims it cites (V1), the citation check (Y2), the
skeptical check (Y3), and the user's approval (U4). A claim added or corrected
at any point, by E4, the user or Y3, goes through one shared validation
function that works out the checks it needs from its kind and its use (quote
check always; E4 always; X1 if it is an equation or value a memo cites; the
user's approval if it is prose supporting a model option).

**Page images, tiles and zoom crops (P1).** Readers and checkers verify claims
against page images, because extracted text is where signs, exponents and
table cells go wrong. Claude does not see an image at full size:
- Opus and Sonnet 5.5 shrink an image to at most 2,576 px on the long edge
  and about 4,784 visual tokens (one token per 28 × 28 px patch, about 3.75
  megapixels). A request with more than 20 images has a stricter per-image
  limit, and an image over it is rejected; we keep every image at 2,000 px or
  less for that reason (Claude docs, "Vision", read 2026-10-29).
- A letter page at 300 dpi is 8.4 megapixels, so Claude sees a whole page at
  about 200 dpi (about 155 dpi at the 2,000 px cap). A tile of part of a
  column keeps the full 300 dpi.

Whether 300 dpi reads small superscripts, subscripts and minus signs better
than 200 dpi has not been measured. **T-res (§12) decides** between three ways
of giving pages to readers and checkers:
1. whole pages;
2. tiles of each page;
3. whole pages, plus a 300 dpi zoom crop of any region the reader or checker
   asks for.

The main session is running T-res now (2026-10-29); its result and the method
chosen are recorded here. What building taught so far:
- an ink-only tiler cut a table mid-row;
- a tiler that keeps layout blocks whole kept a full-page table as one tile,
  which after scaling is about 215 dpi, no better than a whole page;
- a table that must be split is split between rows, with the header row
  repeated on each piece.

**The PDF parser (P1).** P1's page text, page layout and equation LaTeX come
from a PDF parser that runs locally and offline. Marker was used for the first
runs; it is being replaced by a parser that scores higher on document-parsing
benchmarks (candidates PaddleOCR-VL, MinerU and GLM-OCR, tried on this
machine). The choice is tested against known answers: Li et al.'s 58 numbered
equations and the user's transcription of Kaempfer & Plapp p. 7. The parser's
text is a search aid and a second opinion. Every check is made against page
images. The parser's equation LaTeX is also a read made without Claude, and is
compared in X1.

**Tables (B1).** No automated table extraction. A table cell a memo uses is
cited by table, row and column, read from the page image, and checked by E4
against the image, including the column header, its unit and any footnote.

**When an agent call fails.** The driver checks that every input item has
exactly one valid output; a partial output is a failure. A failed call is
retried once, then shown to the user. A call is recorded as done in `runs.db`
only after this check, so a resumed wave redoes any call that did not pass.

## 4. Agents and permissions

**Permission files** (user decision 2026-10-28). Every permission is written
in a file that Claude Code or the operating system enforces. Permission rules
are "enforced by Claude Code, not by the model" (Claude Code docs,
permissions, read 2026-10-28), and a deny in any file overrides an allow in
any other. The files live in `.claude/permissions/literature-review/`, where
every edit asks, so only the user changes them:

```
.claude/permissions/literature-review/
  README.md                     what each file is for, and why
  agents/<agent>.settings.json  one per agent: tools, deny rules, sandbox
  sessions/review-mode.json     the main session whenever it works with library or run files
  inputs.json                   each agent's and session's input classes (untrusted, user_private, public)
  accepted_risks.json           entries the user allows to hold all three legs, with the reason
  hooks/review_read_allowlist.py  review mode's allow-list for Read, Grep and Glob
  launcher.allow                the scripts wf_launch may run, each with its profile
  profiles/<profile>.sb         macOS sandbox profiles: offline, api, anthropic, zotero-local, zotero-web
```

- **The main session's normal profile** is the project's own
  `.claude/settings.json`, which Claude Code always loads. It denies reads of
  the pipeline's untrusted stores: `library/`, `index/`,
  `runs/literature-review/` and `Literature/incoming/`. A normal session
  therefore cannot open library text, so review mode is the only way to
  touch those stores. The user's own `Literature/` folder stays readable for
  physics work (user decision 2026-10-28).
- **`check_permissions.py`** reads these files; it does not write them.
  - It works out each agent's and session's three trifecta legs from what
    its file actually allows: reads of an untrusted store give untrusted
    text; any network beyond the model provider is a way out; reads outside
    the call folder, or of keys, are private data.
  - It fails if any entry has all three.
  - It writes the permission tables below and in `PERMISSIONS.md` between
    markers. The reasons stay hand-written.
  - T-perm checks that the tables are current and that each refusal really
    happens.
- **Limits of the files.**
  - Deny rules cover Claude's own tools and the shell commands Claude Code
    recognizes, not a Python script that opens files itself. The sandbox
    covers those (docs, as above).
  - Permission rules cannot express an allow-list. A `!` carve-out cannot
    reach a rule anchored with `/` or `//`, and cannot reopen a file inside
    a directory blocked as a whole (docs, as above). So review mode uses
    three layers: deny rules naming each untrusted store; a PreToolUse hook
    (`review_read_allowlist.py`) that blocks any Read, Grep or Glob outside
    the allowed paths, since a blocking hook overrides allow rules; and the
    sandbox (`denyRead` of the home folder with `allowRead` for the allowed
    paths) for shell commands.
  - macOS sandbox profiles cannot filter network access by host name. The
    `api` and `zotero-web` profiles allow HTTPS out; the host list is
    enforced by each script's own HTTP client and recorded in the manifest.
    The `offline` and `zotero-local` profiles are exact.
  - Inputs cannot be read off a settings file, so `inputs.json` declares
    each agent's input classes, and the driver checks every file it hands an
    agent against that class.
  - An agent runs with `--setting-sources ""`, so only its own file applies.
    The main session's review-mode file adds to the project settings. Their
    denies still hold, but so do their allows, which is why review mode
    denies WebFetch outright.

**Enforced settings, as checked.**

<!-- BEGIN GENERATED: check_permissions.py (literature-review) -->
Generated from `.claude/permissions/literature-review/` and the project settings by `scripts/literature_review/check_permissions.py`; do not edit by hand.

| Agent or session | Allowed tools | Network | Inputs | Untrusted text | Way out | Private data |
|---|---|---|---|---|---|---|
| `beth-checker` | `Read` | none | untrusted | yes | no | no |
| `beth-memo-checker` | none | none | untrusted, user_private | yes | no | yes |
| `beth-planter` | `Read` | none | untrusted | yes | no | no |
| `beth-recall-judge` | none | none | untrusted, user_private | yes | no | yes |
| `beth-skeptic` | none | none | untrusted, user_private | yes | no | yes |
| `beth-synthesizer` | none | none | untrusted, user_private | yes | no | yes |
| `jerry-screener` | none | none | untrusted | yes | no | no |
| `morty-blind-reader` | `Read` | none | untrusted | yes | no | no |
| `morty-desk` | `Bash(venv_workflows/bin/python scripts/literature_review/library_query.py passages *)` | none | untrusted | yes | no | no |
| `morty-reader` | `Read` | none | untrusted | yes | no | no |
| `morty-scout` | `WebSearch` | web search | public | yes | yes | no |
| `review-mode` | all but web tools | api.openalex.org, api.crossref.org, api.semanticscholar.org, export.arxiv.org, api.adsabs.harvard.edu, ntrs.nasa.gov, www.osti.gov | user_private | no | yes | yes |
| `review-session` | `Read`, `Grep`, `Glob` | none | untrusted, user_private | yes | no | yes |
| `normal` (accepted risk) | all | WebFetch(domain:api.crossref.org); WebFetch(domain:code.claude.com); WebFetch(domain:doi.org); others ask | user_private, untrusted | yes | yes | yes |
<!-- END GENERATED -->

In B1 the agents that run are `morty-reader`, `beth-checker`,
`morty-blind-reader`, `beth-synthesizer` and `beth-skeptic`.
`beth-memo-checker` is merged into `beth-skeptic` (one skeptical checker, two
calls). `jerry-screener`, `morty-scout`, `morty-desk` and `beth-planter` are
deferred (§13). Their settings files stay, and the table above still checks
them.

**Rules.** Every rule below is enforced by a setting, a flag, the operating
system or a script, not only stated in a prompt.
- **Agents write nothing.** They return JSON on standard output, validated
  against a schema, and the driver writes the files.
- **One clean folder per call.** Each agent call runs headless (`claude -p`)
  in a fresh empty folder outside the project, holding only the files it is
  given.
- **The flags of every call:** `--setting-sources "" --strict-mcp-config
  --no-session-persistence --disable-slash-commands --permission-mode dontAsk
  --json-schema <schema> --settings <agent>.settings.json`, plus
  `--allowedTools` for the tools the agent may use.
- **What dontAsk refuses:** any tool not allowed (tested 2026-10-23), and a
  Read outside the folder (tested 2026-10-22 and 2026-10-23).
- **The per-agent settings file** (`agents/<agent>.settings.json`) denies
  Edit and Write everywhere, WebFetch, and Bash except the commands listed
  for that agent.
- **Bash commands** run in the Claude Code sandbox:
  - writes to `library/` are blocked at the operating-system level;
  - network access is limited to listed hosts;
  - the sandbox fails closed: with `sandbox.failIfUnavailable: true`, Claude
    Code exits instead of running commands unsandboxed;
  - with `allowUnsandboxedCommands: false`, no command may retry outside it
    (Claude Code docs, sandboxing, read 2026-10-27).

  The sandbox's write paths and hosts are set per session, not per command,
  and a process inside the sandbox probably cannot apply a second macOS
  sandbox profile (F1 of round 4; T-sandbox tests it). So workflow scripts
  run through `wf_launch`, a fixed launcher:
  - it takes only a script name from an allow-list and the run ID, never a
    free command line, so it cannot be used to run `cat` or anything else;
  - it runs outside the Claude sandbox (`excludedCommands`) and applies the
    script's own profile: for offline scripts no network and writes only to
    the run folder; for drivers only the Anthropic host;
  - it writes the script's standard output and error to the run folder and
    returns only the schema-checked summary, so a traceback that quotes a
    paper never reaches the session.

  Every run records in its manifest that the sandbox and the launcher were
  active.
- **Pinned tools.** Claude Code is pinned with auto-update off. The version
  is recorded per call in `runs.db`, and a mismatch with the run's manifest
  stops the wave.
- **Transcript audit.** `audit_transcript.py` checks every call:
  - the model is the pinned one;
  - every tool call is one the agent is allowed;
  - every Read is of a file the driver handed it.

  A violation fails the call.
- **No agent holds all three** of private data, untrusted content and a way to
  send data out (the "lethal trifecta", Willison 2025).
  - The readers and checkers see untrusted text (the papers) but have no
    network and write nothing.
  - **The main session** has the project and web tools, so it does not read
    the pipeline's stores. Drivers give it counts, IDs and verdicts.
  - **Evidence** is discussed in the review session (below), which may read
    it but has no way out.

| Agent | Model | Tools | Reads | Network | May not |
|---|---|---|---|---|---|
| `morty-scout` | `claude-sonnet-5-5` | WebSearch | user-approved search terms and criteria, seed titles (find); user-approved search terms (challenge) | WebSearch only | fetch or download; return a citation it did not see in a result |
| `jerry-screener` | `claude-haiku-5-5` | none | one batch of 10 candidates; criteria; framing A ("include unless a criterion excludes") or B ("which questions does this bear on") | none | see the other screener |
| `morty-reader` | `claude-sonnet-5-5` | Read | one paper's text, images, tiles, tables, crops, source disagreements (read); claims and a question (retag); one page (page) | none | answer from memory; follow instructions in the paper (reported as a flag) |
| `beth-checker` | `claude-opus-5-5` | Read | claim fields, images, tiles, crops, flags | none | see the parser's text, reader reasoning or other groups' verdicts |
| `morty-blind-reader` | `claude-sonnet-5-5` | Read | one or two tiles and a quantity or equation name per claim | none | see the claimed value, the paper's text or any other read |
| `beth-synthesizer` | `claude-opus-5-5` | none | the claim set as text | none | cite anything but given claim IDs; state a number not in a cited claim |
| `beth-memo-checker` | `claude-opus-5-5` | none | memo, the same claims, broader search results | none | see the synthesizer's drafts or reasoning |
| `beth-skeptic` | `claude-opus-5-5` | none | memo and the cited claims (call 1), for model-critical questions | none | see the other memo checker's verdicts |
| `morty-desk` | `claude-sonnet-5-5` | Bash, only `library_query.py passages` | the user's question; search hits | none | run any other command; its answer goes to the reading view, never to the main session |
| `beth-planter` | `claude-opus-5-5` | Read | one verified claim at a time, with its page image, for drafting a planted error of a named kind (C0) | none | see the checker's prompt or any other plant; a draft is used only after the user approves it |

**Codex** was dropped from the pipeline on 2026-10-28 (§9, item 34). Its
value was a different vendor, whose misreadings would not line up with
Claude's. It had not been measured, it needed a separate macOS user and a
sudoers rule to run safely, and its plan's usage limits could stall a wave.
Independence for equations and values now comes from the PDF parser's
LaTeX (a read made without Claude) and from author sources. The shared blind
spots of Claude models are watched through the misses every read shared, in
C0 and U3. If those appear, Codex
comes back for those claims only. Codex may still be used outside the
pipeline, for the user's own reviews of the design.

| Script | Network | Writes |
|---|---|---|
| `discover.py`, `resolve.py` | api.openalex.org, api.crossref.org, api.semanticscholar.org, export.arxiv.org, ADS, NTRS and OSTI APIs | its run folder |
| `acquire.py` | open-access hosts named by OpenAlex, not known in advance. An excluded command, so it asks once per wave; its URLs come only from OpenAlex and Crossref records, never from the main session. Downloads are never opened by anything but the parsers, and only after the identity check | `Literature/incoming/` |
| `prepare.py`, `build_index.py`, `ocr_math.py`, every `check_*` and `compare_*` script | none, through `wf_launch`. All model weights (the PDF parser's; later the table extractor and the embedding model) are downloaded once at build into a pinned local cache and used offline | its run folder; `index/` |
| `library_query.py` | none | its query log; desk answers to the reading view |
| agent drivers | the Anthropic host only, through `wf_launch` | the run folder |
| `zotero_pull.py` | `localhost:23119` only (Zotero's local API, read-only) | `Literature/incoming/`: copies of PDFs the user saved, matched to the fetch list by DOI |
| `zotero_add.py` | api.zotero.org; key and user ID read from the Keychain (`zotero-api-key`, `zotero-user-id`) | the user's Zotero library: adds items, files and collection membership; never edits an item's fields or deletes anything. It writes to an outside account, so it asks each time it runs |
| `curator.py` | none | `library/` and the archive only. The only writer; it runs outside the sandbox (`excludedCommands`) and asks |

**The review session (B1)** (§9, item 43; Codex finding 11). The user needs
to read and discuss evidence with Claude: quotes, conditions, equations and
memo drafts. Revision 9 kept all of that out of every Claude session (only
counts reached the main session, and the desk wrote to a separate view).
That made evidence review depend on a separate interface and put the user in
the middle of every exchange. For B1 a dedicated session replaces it:
- **Started with** `claude --settings
  .claude/permissions/literature-review/sessions/review-session.json` (to be
  written; `review-mode.json` stays for later branches).
- **May read** only the evidence the driver lists for it: the B1 run's
  claims, applicability checklists, verdicts, memo drafts and page images,
  enforced by the allow-list hook (`review_read_allowlist.py`) and the
  sandbox's read rules.
- **No network:** WebFetch and WebSearch denied, no MCP servers, sandbox
  network none, `curl`, `wget` and `git push` denied.
- **No writes:** Edit and Write denied; Bash limited to read-only commands.
  The user records verdicts in the reading view, and scripts write them.
- **No broad file access:** nothing outside the listed paths; Keychain
  (`security`) denied.

With no way out and no writes, text planted in a paper can at worst mislead
the discussion; the user's approval of each claim a model option rests on
guards against that. Acquisition, publishing and agent calls stay in scripts
run from the normal session. Whether Claude Code's settings and the sandbox
enforce this profile as written is not yet verified; a live refusal test
(T-perm) runs before first use. The revision 9 review mode is kept in §13.

**Main session.** It may read the project, except the pipeline's stores,
and run the scripts. Writes to `library/` are denied (Edit and Write deny
rules plus the sandbox). The curator, pushes, simulations and edits to
`.claude/`, `PERMISSIONS.md` and the agent files all ask first. Only the user
changes permissions. `PERMISSIONS.md` keeps the reasons and holds the tables
`check_permissions.py` writes.

**Agent names.** The prefix names the model family: Sonnet → morty, Opus →
beth, Haiku → jerry.

**Why split this way.** A second agent is used only for one of three reasons:
1. **Different information:** the checker must not see the reader's reasoning,
   and writes its applicability checklist before it sees any claim.
2. **Different permissions:** the review session may read evidence because it
   cannot act on it.
3. **A second opinion on the same object:** the blind reader (equations and
   values) and the skeptical memo checker.

Everything else is a script or a rule. Authoring and checking always run in
separate contexts. A script never decides scientific equivalence or whether
an inference is justified; those go to the checker or the user.

**How independent the checks are.**
- **Claude models:** Sonnet, Haiku and Opus share a vendor and much of their
  training. The checker's main independence comes from its different input:
  images only, without the reader's text. For equations and values it
  transcribes blind before it sees the reader's fields, because a checker
  shown a transcription tends to confirm it.
- **Outside the Claude models:** the PDF parser's equation LaTeX and, where it
  exists, the authors' source give reads of equations that share none of
  Claude's training.
- **Shared mistakes:** models agreeing can still be wrong, so the user
  approves every prose claim a model option rests on.
- **The truly independent checks** are the scripts and the user.

## 5. Checks

| Output | Check | Kind | Pass |
|---|---|---|---|
| copies | page 1 identity; copy kind (A2) | deterministic + user | match, or the user decides |
| quotes | length, presence on the cited page of the imaged PDF, uniqueness | deterministic | all three |
| applicability | E4 step 0 checklist written before claims are seen; dropped conditions checked against it | agent, different input | every condition on the checklist that bears on a claim is in the claim or marked not relevant by the checker |
| claims | `beth-checker` on images | agent, different input | all confirmed |
| prose claims a model option rests on | E4's four meaning fields (§7); the user's approval in U3 | agent + human | the user approves each |
| equations and values a memo cites | X1: the reader, the checker's blind transcription, the blind reader, the PDF parser's LaTeX and the authors' source where they exist, on parameter fields or canonical LaTeX | independent reads, compared by script | all agree, else the user decides |
| table cells a memo uses | cell, header, unit and footnote checked on the image by E4 | agent | confirmed |
| model equations | the paper's symbol and scaling table for each equation a model option uses, verified by the user; Pint parse of printed units; homogeneity where SymPy parses; a known-answer check (conservation, equilibrium or a limiting case) where one applies | deterministic + human | the user verifies the table; failed checks resolved |
| checker corrections | a corrected claim goes through the shared validation function | deterministic + agents | as for any claim |
| recall | gold statements found by a fresh reader run, counted only for the gold statements the brief asks for; U3 *missed* | measured, human | reported |
| repeatability | a repeat run of E1 to E4 on one paper with the same prompts, claims compared one by one (§8) | measured | agreement reported |
| memos | Y2 citations, numbers and units | deterministic | 100 % |
| memos | Y3, one skeptical checker, two calls | agent | no unsupported statement left, or the user accepts each |
| records | rehash; restore test (§6) | deterministic | no mismatch |

**Important claims.** A claim is *important* if any of these holds:
- its kind is parameter or equation;
- by script, whatever its kind: it has a number with a unit, a table cell,
  or LaTeX;
- it bears on a question the brief marks *model-critical*;
- the user marks it.

What important claims get in B1:
- **X1, for equations and values only:** every equation or value claim the
  memo cites, checked before the memo goes to Y2. The random 20 % of other
  important claims is deferred (§13).
- **The meaning check, for prose claims:** a sentence about what a method can
  do, under which conditions, or who showed it cannot be checked by
  transcribing it again. E4 fills in four fields for each prose claim (what
  was established; by whom, this paper or a work it cites; under which
  assumptions; with what evidence), and the user approves those fields for
  every claim a model option rests on.
- **review:** every disputed claim goes to the user;
- **user time:** it grows with the claims a model option rests on, not with
  the number of papers.

**Comparing LaTeX (X1).** `compare_reads.py` canonicalizes LaTeX at the token
level before comparing. It drops only what cannot carry meaning: spacing and
sizing commands (`\left`, `\right`, `\,`, `\displaystyle`) and redundant
braces. It keeps font commands such as `\mathbf`, `\boldsymbol`,
`\mathcal` and `\mathbb`, because they can mark a vector, a tensor, a set or
an operator. Every equivalence is in a tested list, with paired tests:
formatting variants that must match, and changed signs, indices, vector
notation, derivatives and grouping that must differ. T-res reports both
false disagreement (harmless differences flagged) and false agreement (real
differences passed). In B1, a notation difference the canonicalizer cannot
settle goes to the user.

**Numbers in memos (Y2).** Each number in a memo carries an inline tag with
its claim ID and quantity, for example `0.1 [id@v1 alpha_c]`. Y2 takes the
value and unit from the sentence and compares them with that claim's
parameter field for that quantity, so two numbers swapped between quantities
fail. A number with no tag, or a tag naming a claim the sentence does not
cite, fails. Values are normalized (`1e-3`
= `10^{-3}` = `0.001`; percentages as fractions). They must match within the
printed precision, and where the claim has a unit, the memo's unit must have
the same dimension and value after conversion (Pint).
Exempt: IDs, years, and equation and figure numbers. The script ships with
test cases for each rule.

**C0 calibration (B1).** The material the user has written is used as it is.
No new planted errors are written for B1 (§9, item 43).
- **Material:** 40 gold statements (Kaempfer & Plapp 12, Li 10, Molaro 10,
  Blackford 8), 11 planted errors written by the user, and 57 script plants
  the user approved. Kaempfer & Plapp and Li are B1 papers, so their 22 gold
  statements are the B1 recall set. Each of those is labeled at U0 with the
  brief question that needs it, so recall is counted only on facts the brief
  asks for.
- **What each kind of item measures:**
  - omission, by a fresh reader run against the gold statements;
  - attribution and interpretation, by the checker on the planted errors;
  - equations, against the user's transcription and the authors' source
    (T-res).
- **Fixtures the set must contain:** a condition stated on a different page
  from its result (the user's dropped-condition plants, with their condition
  page); a paper discussing another paper's result (a wrong-attribution
  plant); a symbol-sensitive equation (the sign and exponent script plants);
  an important statement left out entirely (a gold statement, recalled or
  not); and at least 10 unchanged claims as controls.
- **Split:** `calibrate.py split` keeps every plant made from one gold
  statement on the same side, and holds out any script plant that shares a
  paper and page with a user plant (built).
- **Diagnostic during development:** the tuning items are run as often as
  needed while prompts are written. Results are reported by kind, by source
  statement and by paper.
- **One frozen final evaluation:** when development ends, the prompts,
  schemas, model IDs and parser version are frozen, and the held-out set is
  scored once. Its result is reported as conditional on this constructed set,
  not as a general error rate. Plants made from one statement are correlated,
  so results are also reported per source statement.
- **Pass, on the frozen evaluation:**
  - **kinds E3 catches by rule** (fabricated quote, wrong page): all caught;
  - **every sign or exponent plant caught** (a fixture rule);
  - **controls:** at most 10 % flagged;
  - **plants that reach E4:** the catch rate with its one-sided 95 %
    Clopper–Pearson lower bound. It is a gate (bound ≥ 85 %) only if the
    held-out set has 19 or more such plants; with fewer it is reported, not
    certified (19 of 19 gives 85.4 %).
- **Per kind:** counts and intervals are reported, with no per-kind gate. A
  gate cannot be met with few items: n caught out of n gives a lower bound of
  0.05^(1/n), so 8 of 8 gives 68.8 % and 18 of 18 gives 84.7 %.
- **User plants and script plants** are reported separately. If the user's
  plants are caught clearly less often, the rate on the user's plants is the
  one that counts.
- **Abstention** (user decision 2026-11-01, from the Lab Report 2 review). A
  plant or control left without a verdict (a failed call, or a claim the
  checker returned no verdict for) counts against the checker: a miss on a
  plant, a flag on a control. Declining to answer cannot raise the score.
  Failed calls are rerun first (saved verdicts are reused); the report gives
  the number of items without a verdict.
- **Sampling temperature.** Claude Code (2.1.316) has no option to set it
  (`claude --help`, checked 2026-11-01), so run-to-run variation is measured
  by repeat runs (T-repeat, §8) instead of being fixed.
- **Script-only baseline** (user decision 2026-11-01; `baseline.py`). The
  same held-out plants and controls, on the same cited pages, judged with no
  model: every number in the claim (notations normalized) must be on the
  page, a number with a unit must keep its unit, and each displayed equation
  must match P1's LaTeX for that page after `tres.canon`. Prose without
  numbers or equations is not judged. It is scored with `calibrate.py score
  --baseline`: the same rules, reported next to the checker's, never a gate
  and never an attempt in the held-out ledger. The comparison the course
  asks for is the checker against this baseline: catches per kind and false
  flags on controls. Tuning trial (2026-11-01, `calib-baseline-dev`, P1 on
  Li, Kaempfer & Plapp and Blackford): 18 of 30 plants caught (changed digit
  10 of 11, sign or exponent 2 of 2, wrong unit 1 of 1, wrong page 5 of 16),
  10 of 30 controls flagged (figure numbers, values printed in other units,
  and the user's LaTeX written differently from the parser's).
- **After a real failure** (an error found in U3, X1 or Y3): no automatic
  redraw of a fresh held-out set, and no rerun of the whole branch. The
  prompt is repaired, the affected claim family (claims of that kind on those
  papers) is inspected again, and the repair is recorded. A later statistical
  assurance would need its own sampling and stopping rules, fixed in advance.
- **Reruns:** the diagnostic is rerun after any change to an agent file,
  prompt, schema, model ID or the parser's version.
- T-res (§12) is part of the calibration.

**Sample bounds.** For k errors in a simple random sample of n = 60, the
one-sided 95 % Clopper–Pearson bound on the error rate is 4.9 % (k = 0),
7.7 % (k = 1) and 10.1 % (k = 2). B1 does not draw such a sample of unused
claims (deferred, §13); the user's census is the claims a model option rests
on.

**Guided reading (part of U3), and its limits.** The user reads each core paper's
key pages in order. After reading a page, the user opens its claims and marks
each one. This checks every claim on those pages, and it is the main human
measure of what the readers missed. Two limits:
- **Chosen pages:** the pages read are chosen, not random, so this says
  nothing about the claims from other pages and papers.
- **Seeing invites agreeing:** claims stay hidden until the page is read, and
  missed statements are recorded before the claim list is opened.

## 6. The library

**Principle.** The truth is a set of immutable files, one per object. The
catalog (and, later, the passage index) is built from them and can be rebuilt
at any time; rebuilding and comparing is part of the audit. Folders sort objects by
kind. Branches and topics are tags, because a paper belongs to several.

```
library/
  POLICY.md                 who may do what
  branches.yaml             branches, their seeds (DOIs), versioned
  vocab.yaml                tag vocabulary with alias map, versioned
  papers/<paper_id>/
    card.v1.json            bibliographic data, identifiers, copies, status, branches, tags
    source/<copy_hash>/     the files acquired: PDF, and LaTeX or JATS source
    content/<copy_hash>/    the whole paper as P1 extracted it; written once, never edited
      pages/p001.md ...     full text of each page, from the best source
      equations.jsonl       every displayed equation: LaTeX, printed number, page
      layout.json           the parser's page layout (block types and boxes)
      images/, tiles/       300 dpi page images, and tiles or crops as T-res decides
      sections.json         section headings with pages
      page_map.json         PDF page to printed page
      manifest.json         every file's hash; the tools and versions used
      (deferred: tables/, figures/, disagreements.json, supplements/; §13)
    analysis/<run_id>/      model output kept but not checked; versioned, written once
      symbols.json          the reader's symbol list and inferred units
      applicability.json    the checker's step 0 checklists
    claims/<claim_id>.v1.json
    symbols.v1.json         symbols that claims use: meaning, unit as printed, page; checked with the claims
  reviews/<review>/         brief.vN.json, waves.json, shortlist.json
  memos/<review>/Qk.vN.md
  views/                    generated pages per branch, tag and review
  events.jsonl              append-only: who, what, when, why, hashes
  catalog.sqlite            index of the records; rebuilt from the files
  history/                  superseded versions
  .staging/                 publishes in progress
  references.bib            generated
```

- **IDs.** `paper_id` is internal and assigned once (`kaempfer_2009_k3f9`);
  the card lists every identifier and every copy read, with its kind and
  hash. `claim_id` is `<paper_id>.c007`. IDs never change or get reused.
- **Three layers per paper.**
  - **Content** (`source/`, `content/`) is the whole paper, extracted by
    scripts with no model rewriting. It is complete but not checked: an
    equation the parser transcribed from a PDF may be wrong, and the record says
    which source each page came from. Any question asked later can be
    answered from it without fetching the paper again.
  - **Analysis records** (`analysis/`) are model output about the paper that
    is kept but not checked: symbol lists, inferred units and applicability
    checklists. They are separate from the content, because the content is
    written once by P1 and these are written later, by other stages, and may
    have several versions.
  - **Claims** are statements about the content that have passed E3 and E4.
    Only claims may be cited in a memo.

  A new copy of a paper (say, the published version after the arXiv one)
  gets its own `content/<copy_hash>/`; the old one stays, and the card marks
  which copy is current.

  The reader's symbol list is model output, so it sits in `analysis/` as
  unchecked. A symbol enters `symbols.v1.json` only when a claim uses it, and
  E4 checks it with that claim.
- **Cards** hold bibliographic data, copies, status, branches and tags only.
  Anything said about a paper's content is a claim and goes through E3 and
  E4.
- **Claim fields:**
  - text and kind (method, equation, result, assumption, limitation,
    parameter);
  - quote (and the image wording if E4 corrected it);
  - page (PDF index and printed) and location;
  - equation (LaTeX, number);
  - parameter (quantity, value, unit as printed, uncertainty), so that X1,
    E3b and Y2 compare fields rather than prose;
  - table cell (table, row, column) and figure;
  - questions (versioned) and tags;
  - basis (derived / numerically tested / experimentally tested / cited /
    asserted);
  - conditions with pages;
  - for prose claims, the four meaning fields (§7): what was established, by
    whom, under which assumptions, with what evidence;
  - importance;
  - copy;
  - status, verdicts (checker, model, run; X1; user), version.

  Memos cite `[id@v1]`.
- **Tags.** A small, fixed vocabulary in `vocab.yaml`: the branch, plus three
  facets (method family; material and setting; physics). Tags exist for exact
  filters, such as a memo's claim set or a view by branch. Finding related
  material by meaning is the passage index's job (user decision 2026-10-27).
  - **Assigning:** the reader tags from the vocabulary, and E4 confirms.
  - **Changing:** only the user adds or renames a tag, at U4, as a new
    `vocab.yaml` version. A renamed tag keeps its old name as an alias, so
    old records still match. Records are never rewritten.
- **Change and deletion.** A correction is a new version; the old one moves to
  `history/` with an event. Nothing is deleted. Only the user can mark a
  record retracted, with `curator.py retract --user-approved "<reason>"`.
- **Protection:**
  - deny rules and read-only files;
  - the sandbox for Bash writes;
  - hashes and the event log (L2).
- **Backup and restore (B1).** After the memo is approved, the B1 library and
  run folders are copied to a dated folder in the course folder
  (`~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Courses/Fall2026/ME295/`),
  with a manifest of every file's hash. **T-restore** (§12) restores one paper,
  its claims and the memo's citations into an empty folder and reruns Y2 on
  them. The content-addressed archive of revision 9 is deferred (§13). The
  copy holds paper text, so it is never in git.

- **Zotero** (Z1) is a copy for the user's own reading and citing. The
  library files stay the truth. The pipeline reads Zotero back only through
  `zotero_pull.py`, and only the PDFs of papers on the fetch list that the
  user saved (U2); those go through A1 and A2 like any other copy.
- **What is committed.**
  - Committed: `branches.yaml`, `vocab.yaml`, briefs, `references.bib`, the
    event log, run manifests, search snapshots (bibliographic metadata) and
    `index/model.lock`.
  - Not committed: PDFs, sources, page images, text, chunks, the index,
    claims, quotes, tables, memos and abstracts (copyright). For the course,
    the evidence of reproducibility is the repeat run (§8) and the committed
    prompts, schemas and manifests, not the papers.

**Ranking.** No relevance score is stored on a paper. A model's score is
uncalibrated and biases every later step that sorts by it. The library
records checkable facts instead: tier with its criterion, questions,
basis and importance. SC0's ranking is a baseline for measurement and stays
in its run folder.

## 7. Schemas (summary)

- **Read:** `{paper_id, section_map, claims: [{text, kind, quote, page_pdf,
  page_printed, location, equation?: {latex, number}, parameter?: {quantity,
  value, unit_printed, uncertainty?}, table_cell?: {table, row, column},
  figure?, questions[], tags[], basis, conditions: [{text, page_pdf}]}],
  symbols: [{symbol, meaning, unit_printed, page_pdf, nondimensional?}],
  nothing_relevant_pages[], flags[]}`.
- **Claim verdict:** `{claim_id, text: supported|overstated|unsupported|
  not_on_cited_page, corrected_text?, quote_on_page, image_quote?, equation:
  matches|corrected|none, corrected_latex?, parameter: match|corrected|none,
  corrected_parameter?, kind, basis, questions (confirmed|corrected),
  dropped_conditions[], reason}`.
- **E4 applicability checklist (step 0):** `{paper_id, pages_read[], items:
  [{kind: assumption|boundary_condition|field_role|parameter_limit|
  validation_regime, text, page_pdf}]}`.
- **Meaning fields of a prose claim (E4 step 2):** `{established, by_whom:
  this_paper|cited_work, cited_work?, assumptions[], evidence: derived|
  numerically_tested|experimentally_tested|cited|asserted}`, plus
  `user_approved` once U3 rules.
- **E4 blind transcription (step 1):** `{claim_id, latex? | parameter?:
  {quantity, value, unit_printed}, not_found?}`.
- **X1 read** (blind reader; the parser's LaTeX and the author source are
  turned into the same form by script): `{claim_id, latex? |
  parameter?: {quantity, value, unit_printed}, not_found?}`.
- **Memo number tag:** `<value> [<claim_id>@v<k> <quantity>]`.
- **U3 verdict:** `{claim_id, user_verdict, correction?, note?}`; missed:
  `{paper_id, page_pdf, statement}`.
- **Full schemas:** in `workflows/literature-review/schemas/` when built.

## 8. Reproducibility

Each wave writes `runs/literature-review/<review>/<run_id>/`:

- **`manifest.json`** (committed). It records:
  - the git commit, the pinned Claude Code version (also recorded
    per call in `runs.db`), and the model ID per agent, with its deprecation
    date where the vendor gives one;
  - the hashes of agent files, prompts, schemas, the brief, every file
    acquired, every page image and text page, and the index;
  - the versions, weights hashes and flags of the PDF parser and mutool. The hashes are computed from the files actually loaded
    from the pinned cache, and the run fails if they differ from the lock
    file;
  - that the sandbox and `wf_launch` were active;
  - the hash of the verified units table;
  - the backup manifest's location and hash;
  - the calibration result and sample seeds;
  - the query filters and returned IDs for Y1 and Y3;
  - verdict counts.
- **`runs.db`** and **`transcripts/`** (local, and in the archive): one row
  and one transcript per agent call.
- **`untrusted/`**: every file in a run folder that holds text from papers,
  the web or APIs (reads, search responses, candidates, scout output,
  transcripts) goes in this subfolder. The project settings deny it to every
  session, and review mode reads only `summary.json` and `manifest.json`.
- **Scripts.** Scripts in `scripts/literature_review/` refuse to run with
  uncommitted changes to the workflow's definitions, print their inputs, and
  run in `venv_workflows`.

**What "reproducible" means here.**
1. **Deterministic stages give identical hashes from the same inputs.** These
   are E3, Y2, L1 and the restore test.
   - **P1** uses a neural PDF parser, so it is tested for this (T-det). The archived outputs are the reference that
     later checks read.
2. **Agent checks repeat at a measured rate.** Rechecking a wave's claims
   gives the same verdicts at an agreement rate measured on 30 random claims
   (first two waves, and after any change). Pinned models will be retired,
   so an agent verdict is reproducible through its archived transcript, not
   by rerunning it years later. Dated model IDs are used where offered.
3. **The course's test: the same prompts give similar results** (§9, item 42).
   E1 to E4 are run twice on one B1 paper with the same prompts, schemas and
   model IDs. The two claim sets are matched (same page, same kind, quotes
   overlapping by at least half of the shorter quote's words) and compared:
   claims found by both runs, claims found by one, and verdicts that differ.
   The comparison, the prompts and the manifests are committed; the claims
   are not (they quote the paper).

## 9. Decisions

Resolved 2026-10-23 (revision 5):
1. **Practice library:** frozen as `library_practice/`; seeds only.
2. **Access:** open sources only; no institutional access for now; no
   Sci-Hub (it distributes copyrighted papers without permission, and this
   is public course work).
3. **Models:** Opus for the checker, synthesizer and memo checker; Fable is
   not available.
4. **E2:** core papers only.
5. **Memos:** not committed.
6. **U3:** guided reading (U3a) and a random 60 (U3b).
7. **OpenAlex key:** in the Keychain.
8. **Critiques:** rounds 1 and 2 applied (§11).
9. **Page tiles:** adopted.
10. **Sandbox:** adopted.
11. **U3b remedy:** after any error, a fresh sample of 60.
12. **Backup:** to the course folder.

From the user, 2026-10-27 (revision 6):
13. **Branches and seeds:** design the library's branches first, seed each,
    then run waves.
14. **Codex** is the cross-vendor model whenever one is compared.
    Superseded on 2026-10-28 (item 34).
15. **Passage index:** hybrid search over full text.
16. **Tables:** handled as structured data, with their own chunks.
17. **Units:** checked by dimensional analysis.
18. **Permissions:** explicit for every agent and script.
19. **Critique round 3:** an independent review looking for LLM steps that a
    script could do. Applied in revision 7 as triaged
    (`runs/design-critique/20261027-170522/triage.md`).
    - **Decided by the user:**
      - the archive goes to the course folder (F2);
      - Codex runs as a separate macOS user (F5; moot since item 34);
      - calibration plants are split into tuning and held-out sets (F10);
      - the user rules on 30 practice-library claims (F25).
    - **The other 21 findings** were applied as recommended and approved by
      the user.
20. **Branches:** B1 to B5 as in §2, approved; work starts with B1.
21. **Tags:** a fixed vocabulary; no tagging agent and no tag growth
    process (§6).
22. **Zotero:** every shortlisted paper is to be added to the user's Zotero
    (Z1), with its PDF when acquired. Deferred on 2026-10-27; taken up again
    on 2026-10-28 (item 36).
23. **Codex data:** the user set their OpenAI account not to train on their
    data (2026-10-27).

From the user, 2026-10-28 (revision 8):
31. **Whole paper stored:** each paper's library record holds the whole
    paper, with claims as a checked layer on top (§1, §6).
32. **B1 scope:** includes methods that couple a sharp-interface
    description with a diffuse-interface one.
33. **Critique round 4** (`runs/design-critique/20261028-061227/triage.md`):
    - **F3:** the desk runs headless with no network and writes to the
      reading view; review mode whenever a session works with library or run
      files;
    - **F13:** X1 on claims memos cite, claims on model-critical questions,
      and a random 20 % of the other important claims;
    - **the rest** applied as recommended (F20 refuted and replaced by the
      2,000 px rule; F24 rejected).

From the user, 2026-10-28 (revision 9):
34. **Codex dropped** from the pipeline (§4). X1 is a blind read by
    `morty-blind-reader` plus a math OCR model; the second memo checker is
    `beth-skeptic`; model-drafted plants come from `beth-planter`. This
    removes item 25 (Codex setup) and the separate macOS user.
35. **Permission files** are `settings.json` files in
    `.claude/permissions/literature-review/`, read and checked by
    `check_permissions.py`; `PERMISSIONS.md` keeps the reasons and gets
    generated tables. The project settings deny normal sessions the
    pipeline's stores; the user's `Literature/` folder stays readable.
36. **Zotero:** Z1 is built. Caltech gives the user unlimited storage, so
    every PDF is attached. The user saves papers the pipeline cannot get
    with the Zotero Connector (personal Chrome profile; proxy redirection
    on), and `zotero_pull.py` copies them in.
37. **Caltech Library downloads** by script are deferred. The instructor
    approved using the library, but its daily download limit is unknown.
    They come back if papers are hard to find.
38. **B1 seeds** (item 24):
    - Li et al. 2009, 10.4310/cms.2009.v7.n1.a4;
    - Lervåg & Lowengrub 2015, 10.4310/cms.2015.v13.n6.a6;
    - Yu, Chen & Thornton 2012, 10.1088/0965-0393/20/7/075008;
    - Folch & Plapp 2005, 10.1103/physreve.72.011602;
    - Steinbach 2009 (review), 10.1088/0965-0393/17/7/073001;
    - Plapp 2011, 10.1103/physreve.84.031601;
    - Joshi & Jaiman 2018, 10.1002/nme.5961 (a sharp moving-mesh solid
      coupled with a diffuse two-fluid interface; the coupling the user is
      looking for).

    All seven DOIs resolve in Crossref. Seeds per branch are now 3 to 7.
39. **Revision 9 approved** for building (user, 2026-10-28).
40. **Codex review accepted in part** (user, 2026-10-29; `runs/design-critique/20261029-codex/`):
    - finding 1: a small, timeboxed B1 path; LaTeX/JATS, automated tables
      and figures, embeddings, the desk and discovery are deferred. To be
      written up as revision 10;
    - finding 2: blind transcription (X1) only for equations and values;
      prose claims get a meaning check and the user's approval when they
      support a model option;
    - the other findings: accepted, item 43.
41. **B1 memo question** (user, 2026-10-29): how to couple the stationary
    grains with the evolving ice–vapor interface without a spurious air film.
    Surface diffusion and other vapor-diffusion paths are secondary.
42. **Reproducibility for the course** (user, 2026-10-29): the same prompts
    given to the agents give similar results. A repeat run on one paper,
    compared claim by claim, is the evidence.
43. **Codex review accepted in full** (user, 2026-10-29): every recommendation
    in `runs/design-critique/20261029-codex/triage.md` (findings 1 to 15),
    applied in revision 10. Finding 11: a dedicated review session that may
    read selected evidence but has no network, no writes and no broad file
    access; the desk and the counts-only rule are deferred. Finding 15: one
    restore test; the content-addressed archive deferred; reproducibility as
    in item 42. Finding 7: no new planted errors for B1.

Open:
24. **Seeds** for B2 to B5, when each branch is taken up.
26. **Zotero setup:** the user ID is given; the user stores it and a
    write-enabled API key in the Keychain (`zotero-user-id`,
    `zotero-api-key`), and confirms the desktop app syncs to the same
    account.
27. **Embedding model:** deferred with the passage index (§13).
28. **Second table extractor:** deferred with automated tables (§13).
29. **First brief** for B1 (U0): `B1_CONTRACT.md` drafted 2026-10-29; the
    user approves it and labels each B1 gold statement with its question.
30. **Calibration material:** written (40 gold statements, 11 user plants,
    57 approved script plants, T-res pages). Done.
44. **PDF parser:** replacing Marker; chosen by its tests in §3 (P1).
    Result (2026-10-29, commits 16c5d280, 2c1d9f5b): **MinerU2.5-Pro**,
    run locally on the Mac's GPU in `venv_docparse/`, offline under the P1
    sandbox (`prepare.py --parser mineru`, the default). Li et al.: 58 of 58
    numbered equations; Kaempfer & Plapp p. 7: 10 of 10 equations exact
    against the user's transcription (Marker 7 of 10); Molaro 41 of 41
    numbered equations (Marker 37 of 42); about 10 to 15 s per page.
    PaddleOCR-VL-1.6 was not installed: on this Mac it needs a separate
    model server, which the no-network sandbox blocks. Marker stays as a
    fallback (`--parser marker`).
45. **Page tiling:** whole pages, tiles, or zoom crops; chosen by T-res.
    Result (2026-10-29, `runs/literature-review/tres/summary_tres.json`):
    on the four T-res pages (20 numbered equations, Sonnet 5.5 reader,
    3 fresh reads each), whole pages gave 60 of 60 exact; ink tiles gave 59
    of 60 (one superscript misplaced). Whole pages are used for reading;
    tiles are not needed for equations. A zoom crop of a region is the
    fallback if a later test (tables, small numbers in text, scans) fails.
    Not yet tested: tables, numbers with units in running text, scanned
    pages, the Opus checker.

## 10. Risks and failure modes

| Risk | Caught by |
|---|---|
| Invented reference | D3 must resolve it; unresolved shown at U1 |
| Report or thesis without a Crossref record | D3 via arXiv, ADS, NTRS, OSTI; never dropped silently |
| Retracted or corrected paper | Crossref status at D3 and every L2 |
| Memo cites a superseded claim | versioned citations; Y2 at every L2 |
| Question added after reading | versioned questions; R1 with passage search |
| Claim tagged to the wrong question | E4; Y3's broader search |
| Claim set too large for one context | split synthesis; coverage script |
| Partial agent output, or a run stops | one valid output per item; retry; resume redoes unchecked calls |
| Library lost | dated backup with a hash manifest; restore test (T-restore) |
| Search truncated by quota or outage | budget read; backoff; *discovery incomplete* |
| Preprint and journal version both read | D3 merge; A1 links copies |
| Relevant paper excluded | named criterion; no abstract never excluded; U1 sample; SC0 baseline comparison |
| Wrong PDF or a differing preprint | A1 identity; A2 comparison; copy on every claim |
| Text layer drops signs and exponents | never given to readers; used only to flag disagreements |
| Small symbols misread | page images at the resolution T-res chooses; X1 with the parser's LaTeX and author sources; T-res |
| Table shredded (rows lose their headers) | memo-used cells read and checked on the image, with header, unit and footnote; a split falls between rows with the header repeated |
| Wrong or lost unit, lost exponent | Pint parse; units table; equation homogeneity |
| Right unit, wrong number | X1 (blind reader and math OCR model, both blind to the value); E4 |
| The parser misreads an equation and the reader copies it | E4 on images only; X1; C0 plants |
| Fabricated or generic quote | E3 |
| Claim stronger than the paper, or missing a condition | E4 with the key pages; branch C0 |
| Reader misses a key result | gold statements in calibration; U3 on key pages; Y3 omissions call |
| Reader misses a condition stated on another page | the checker's applicability checklist (E4 step 0), written before it sees claims; dropped-condition plants in calibration |
| A prose claim reads right but means something else (the equation matches, the conditions do not) | the four meaning fields in E4; the user's approval for claims a model option rests on |
| A paper steers the reader to a false but well-formed claim | E4 on images; the user's approval; the planted-instruction fixture checks the scientific claim (T-trifecta) |
| Calibration gate impossible at small sizes, or held-out sets used up by reruns | no per-kind gate; one frozen final evaluation; repair and inspect after a failure (§5) |
| Pipeline work crowds out the model | 3 to 5 day build budget; stopping rule (§1); deferred list (§13) |
| All three models share a blind spot | U3b random sample; growing calibration set |
| Search misses the passage | T-ret measures it; search is never the only path to a fact |
| Embedding model or index changes | model pinned by revision and hash; index hash; rebuild and compare at L2 |
| A memo calls something a gap that the literature has answered | D2g challenge before a gap is stated |
| Instructions hidden in a paper or web page | readers have no web, write or shell; schema-only output; transcript audit; no agent holds the lethal trifecta; untrusted text kept out of the main session |
| Checker's own correction published unchecked | a correction is a new version, checked like any claim |
| Checker confirms the reader's transcription instead of reading | blind transcription first (E4 step 1); X1 compares it |
| Calibration passes on easy or reused plants | script-caught kinds gated apart; held-out set used once per prompt version; E4-found errors only for tuning; copy kinds recorded |
| Unchecked content cited as if checked | only claims cited; Y3 judges support only on cited claims; symbol list kept unchecked |
| A memo number put on the wrong quantity | number tags with claim and quantity; Y2 compares the parameter field |
| Retracted paper used before the user rules | its claims blocked from Y1 and failed by Y2 once L2 flags it |
| Fact only in a supplement | supplements acquired and stored with the paper; missing ones listed on the card |
| Untrusted text through a script's output or a desk answer | `wf_launch` returns only checked summaries; desk answers go to the reading view; review mode by allow-list |
| Two papers in a series merged | year and series-token check; title-only merges at U1 |
| All Claude reads share a blind spot | the parser's LaTeX and author sources in X1; misses shared by every read tracked in C0 and U3; Codex can return for those claims |
| A permission written only in prose | every permission in an enforced settings file; `check_permissions.py`; T-perm |
| Sandbox fails to start | fails closed; recorded in every manifest |
| Recorded results cannot be rechecked after a disk loss | archive in the course folder, hash checked at L2 |
| Library changed by mistake | deny rules, read-only files, sandbox, hashes, audit |
| Publish interrupted | staging; event last; `recover` |
| Tool update changes results | versions and weights pinned; recalibration triggers |
| Tag sprawl or drift | fixed vocabulary; only the user changes it; aliases |
| Scope creep | tiers, wave limits, staged waves, U4 |
| Search not reproducible | raw responses saved and hashed |
| API key leaked | Keychain only; removed from URLs |

## 11. Independent critique of this design

- **Rounds 1 and 2.** `beth-critic` (Opus) reviewed revisions 2 and 3:
  - round 1 (`runs/design-critique/20261023-103742/`): 23 findings, all
    accepted;
  - round 2 (`20261023-111440`): 22 findings, 21 accepted.
- **Round 3** (`20261027-170522`, revision 6). `beth-critic` got this file
  and extra questions (`runs/design-critique/focus/rev6_scripts_permissions.md`):
  - which agent steps could be a script, or a script plus a smaller agent
    decision;
  - whether each agent is justified;
  - where the permissions in §4 are weaker than stated.

  It returned 25 findings and 9 unclear points, all applied in revision 7.
- **Round 4** (`20261028-061227`, revision 7 with the whole-paper edit).
  Questions in `runs/design-critique/focus/rev7_whole_plan.md`: the weakest
  points, review mode, the two-layer library, calibration statistics and
  simplicity. It returned 25 findings (19 major) and 13 unclear points; 23
  were applied in revision 8, F20 was refuted against the docs, and F24 was
  rejected.
- **Round 5, Codex** (`20261029-codex`, revision 9 plus a status note). The
  user ran Codex (gpt-6.1-sol) outside the pipeline, from a clean
  configuration, read-only, with web search off (its log shows no search).
  The prompt (`packet/PROMPT.md`) asked for the smallest trustworthy path to
  B1, checks that pass on wrong output, the worth of the calibration, agents
  against scripts, the proportion of the security model, and what would make
  B1 stall. It returned 15 findings (1 critical, 14 major) and 7 open
  questions. The main session confirmed finding 8's arithmetic and corrected
  its own status note (B1 has 7 seeds and 9 papers). The user accepted every
  recommendation (§9, items 40 to 43), applied in revision 10.
- Codex may review the design outside the pipeline, run by the user from a
  clean configuration with web search off.
- **How a critic runs.** It sees only a copy of this file, in a fresh folder,
  headless, with Read and WebSearch. It is not told the design is good.
- **Findings.** Each finding gives the stage, the problem, a concrete failure,
  a fix, the severity and its source.
- **Triage.** The main session recommends a decision for each finding and
  checks facts where it can; the user decides. Nothing is dropped unseen.
- **Record.** `runs/design-critique/<run_id>/`.

## 12. Tests at build time

Each test runs before its stage is used on real papers. Results go in the
build log with the date. The tests B1 needs:

| Test | What it checks | Pass |
|---|---|---|
| T-res | equations and numbers on 4 dense pages (Moure & Fu arXiv v1 pp. 8, 14, 19, against the authors' LaTeX; Kaempfer & Plapp p. 7, against the user's transcription), read from whole pages, from tiles and from whole pages with zoom crops | the method with no symbol errors on these pages is chosen; the canonicalizer's false disagreement and false agreement reported |
| T-parser | the PDF parser candidates on Li et al. (58 numbered equations) and Kaempfer & Plapp p. 7 (the user's transcription), on this machine | the chosen parser finds every numbered equation; its LaTeX error rate on p. 7 reported; run time and disk use recorded |
| T-det | P1 twice on one paper | identical hashes, or the differences listed |
| T-units | Pint and SymPy on the verified units table and on the plants with unit and exponent errors | every planted error flagged; no verified row flagged |
| T-applicability | the checker's step 0 on the user's dropped-condition plants, whose condition is on a different page from the result | the condition appears in the checklist and the plant is caught |
| T-perm | `check_permissions.py` passes; each agent's forbidden actions are refused (reader WebSearch, Read outside its folder, any Bash write to `library/`, an offline script reaching the network, a script name not on the launcher's list); the review session's WebFetch, Write, Edit, a `curl`, and a Read outside its listed paths are refused, live; the transcript audit catches a planted violation | all refused and caught |
| T-sandbox | sandbox on: Bash writes to `library/` refused; `curator.py` runs; with the sandbox made unavailable, Claude Code exits instead of running unsandboxed; whether a sandboxed process can apply a second profile | all; the nesting result recorded |
| T-trifecta | a fixture paper with a known fact and a planted instruction to report its opposite (and one to fetch a URL): the reader flags the instruction; the published claim keeps the true fact or is blocked; the review session's WebFetch and a `curl` are refused | the scientific claim is right or blocked, and every forbidden action refused |
| T-restore | one paper, its claims and the memo's citations restored from the dated backup into an empty folder; Y2 rerun | hashes match; Y2 passes |
| T-repeat | E1 to E4 twice on one paper with the same prompts (§8) | agreement reported |

Deferred tests (T-table, T-ret, T-blind with a separate math OCR model,
T-retract, T-key, T-backup, T-zotero) are in §13 with their original rows.

## 13. Deferred

Everything below was part of revision 9 and is kept, word for word where it
is quoted, for the branches that need it. None of it runs in B1. Each item
comes back when a branch needs it, through a new revision and the user's
approval. Revision 9 itself is in git (commit `98302cda`).

**Why each part is deferred.**
- **Discovery, screening and acquisition by script** (D1, D2, D2c, D3, SC0,
  SC1, SC2, U1, A1, Z1, U2, D2g, R1): B1 reads nine fixed papers that are
  already on disk.
- **Tiers and waves:** one fixed set needs neither.
- **The second reader (E2) and merging two reads:** the checker's
  applicability checklist and the user's review cover the misses E2 was for,
  at lower cost.
- **LaTeX and JATS text, source comparison, automated tables and figures:**
  B1 reads equations and table cells from images, and only those a memo uses
  are checked further.
- **The passage index (P2) and the reference desk:** nine papers do not need
  search, and the review session replaces the desk.
- **The general units operator layer (E3b):** dimensional checks cannot catch
  signs, coefficients or scalings; the symbol and scaling table for the
  model's equations does more for less.
- **X1 on a random 20 % of other important claims, U3b's random 60, and the
  practice-library recheck:** statistical audits of claims no memo uses.
- **The branch calibration's light part and new plants:** the existing
  material is used as a diagnostic (§5).
- **Review mode (counts-only main session):** replaced by the review session.
- **The content-addressed archive:** replaced by a dated backup and a restore
  test.
- **`beth-planter`, `jerry-screener`, `morty-scout`, `morty-desk`,
  `beth-memo-checker`** (merged into `beth-skeptic`), **`ocr_math.py`** (the
  parser's LaTeX serves as the read made without Claude), **`zotero_add.py`**.

### 13.1 Tiers, waves and seeds in discovery (revision 9, §2)

**Seeds.** For each branch the user picks 3 to 7 papers, at least one a
review, and at least one or two papers the user knows well. They are listed
with their DOIs in `library/branches.yaml`. Seeds are read first (wave 0 of
the branch) and steer what follows:
- **Discovery:** D1 follows every seed's references and citing works.
  - The seeds' author keywords are offered as query terms, which the user
    approves at U0.
- **Scout:** the scout gets the seed titles as examples of what is in scope.
- **Screening:**
  - A seed must never be excluded; this is checked in C0.
  - SC0 ranks candidates by similarity to the seed abstracts.
- **Calibration:** the one or two seeds the user knows well carry the
  branch's calibration (§5).

**Order of work.**
1. B0, build the library: scripts, folders, index, permissions, and the
   tests in §12.
2. Seed each branch (wave 0) and run its branch calibration.
3. Run the core C0 once.
4. Run waves.

**Waves.** A wave is one pass through stages D to L on a bounded set of
papers.
- **Between waves**, the user can expand the scope (new questions, new
  queries, citation chaining from included papers) or narrow it (drop
  questions; their records stay in the library).
- **Brief versions:** the brief is versioned (`brief.v1.json`, …), and each
  wave records the version it ran under.

The tier sets which stages run, and the limits:

| Tier | Use | Discovery | Screening | Second reader (E2) | Wave limit |
|---|---|---|---|---|---|
| S | "what does paper X say about Y" | seeds only | skipped | no | 5 papers |
| M | one method family | APIs, chaining, scout | SC0 and 2 screeners | core papers | 20 papers |
| L | a whole field | as M, 2 chaining hops | as M | core papers | 20 per wave, small shakedown wave first |

In tier S, U1 is still shown (it is short), A2 runs, and E2 does not.
Branch C0, X1, U3a and U3b still apply.

### 13.2 The full stage table (revision 9, §3)

| # | Stage | Agent or script | Input | Output | Check | Stop / escalate |
|---|---|---|---|---|---|---|
| U0 | Brief | user, with the main session | goal, branch, seeds | `reviews/<review>/brief.vN.json`: questions Q1…Qn (each marked *model-critical* or not), criteria, tier, seeds, queries, limits | user approves | — |
| D1 | Search and chain | `discover.py` | brief, seeds | `candidates.json`; every request and response saved with hash and date | every query logged with its result count; OpenAlex cost estimated before the wave and remaining budget read on each response; Semantic Scholar citation contexts and influential flags kept where present | Keychain key unreadable → stop; HTTP 429 or 5xx → back off and resume; a source still failing → *discovery incomplete*, shown at U1, queries rerun |
| D2 | Scout (M, L) | `morty-scout`, task *find* | search terms and criteria the user approves at U0 for sharing, and seed titles; never the brief's questions, which are the user's unpublished plans (the scout has the web, so it may hold no private data) | citations only (title, authors, year, DOI if known, why) | resolved in D3 | — |
| D2c | Challenge (on demand) | `challenge_terms.py`, then `morty-scout`, task *challenge* | one statement the user wants to rely on (a novelty claim, a memo gap); the script turns it into keyword queries, which the user approves; the scout sees only the approved terms | citations that might contradict it, each with why | resolved in D3; every hit read through E1–E4 before it counts | result goes to the user |
| D3 | Resolve | `resolve.py` | candidates and scout citations | one entry per work, merged by DOI, OpenAlex ID, Crossref preprint relation, then title and first author; abstract present or not; retraction and correction status | resolves in Crossref or OpenAlex, else arXiv, ADS, NTRS or OSTI; a title match (§5) also needs years within ±1 and no differing series token (Part I / II, numerals) | unresolved → kept, shown at U1; a merge by title alone → shown at U1 |
| C0 | Calibrate | `calibrate.py` | §5 sets | `calibration.json` | §5 | below threshold → stop, revise prompts |
| SC0 | Prefilter and baseline (M, L) | `screen_baseline.py` | candidates, brief, seed abstracts | excluded by rule (document type: erratum, editorial, book review; outside the brief's years or languages), each with the rule; the rest ranked by hybrid similarity to the seed abstracts and the brief, using the index's embedding model on a separate store of candidate abstracts | script; the ranking stays in the run folder and is never stored on a paper | — |
| SC1 | Screen (M, L) | `jerry-screener` ×2, two framings | title, abstract, venue, year, criteria | core / context / exclude, with the exclusion criterion, and the questions it bears on | schema; exclude needs a named criterion | no abstract → never excluded, shown at U1 |
| SC2 | Combine screens | `compare_screens.py` | both screens, SC0 ranking | by rule: either screener includes → include (core if either says core); both exclude on the same criterion → exclude; both exclude on different criteria → U1; agreement of each screener with the baseline | script | — |
| U1 | Approve shortlist | user | core and context papers; a random 60 of the exclusions (rule and screener exclusions, sampled in proportion); unresolved and no-abstract candidates; queries and counts; *discovery incomplete* | approved shortlist | the sample bounds the share of wrongly excluded papers (§5) | a relevant paper in the sample → screening prompts revised, C0 screening rerun, wave rescreened |
| A1 | Acquire | `acquire.py` | shortlist | PDFs, and LaTeX or JATS XML source where openly available; supplementary files and correction notices linked in the Crossref or JATS record; each in a new folder under `Literature/incoming/<review>/`; `acquired.json` (path, sha256, source, copy kind) | DOI or title on page 1 matches; a new copy of a known work is linked to it | no match → logged, not acquired; a newer version of a read paper → its claims rechecked (E3, E4); a supplement or correction that cannot be acquired → listed on the card |
| A2 | Check the copy (any copy that is not the publisher's version) | `check_copy.py` | copy, published metadata | copy kind; title, author and abstract comparison | match; differences listed on the card | title or authors differ → U2 |
| Z1 | Add to Zotero | `zotero_add.py` | approved shortlist, `acquired.json` | per paper: a Zotero item (Crossref metadata) in the collection `lunar_regolith_DSM/<branch>`, tagged with the review; the PDF attached when acquired (Caltech gives the user unlimited Zotero storage); papers without a PDF are not added, so the user's own save of them makes the item; `zotero_map.json` (paper → Zotero key) | the library is searched by DOI (by title without one) first; an existing item is never duplicated or edited, only added to the collection; the uploaded file's hash matches `acquired.json` | API error → item added without its file, logged, shown at U4 |
| U2 | Supply missing papers | user | `fetch_list.md` (each paper as a DOI link); A2 differences | the user saves each paper they can get into Zotero with the Zotero Connector, or marks it "unavailable"; accepts or rejects flagged copies | `zotero_pull.py` copies the saved PDFs into `Literature/incoming/`; A1 and A2 rerun on them | — |
| P1 | Prepare | `prepare.py`, then `curator.py` publishes | acquired files, supplements included (under the same paper ID) | the paper's full content (§6, `content/`): page images (300 dpi); tiles (300 dpi, below); Marker pages; text from the best source (below); equations; tables (below); figure crops with captions; section map; a map from PDF page to printed page; source disagreements; hashes. Written to the run folder, then published unchanged into the paper's library folder | every page has an image, tiles and text; every tile within the image limits; tiles cover the text area; the published files hash the same as the run folder's | a failed page is logged and cannot be cited |
| P2 | Index | `build_index.py` | the library's `content/` folders | passage chunks with metadata, keyword index and embeddings, in `index/` (§6) | a rebuild from the same inputs gives the same canonical dump (§8, T-det) | — |
| E1 | Read | `morty-reader`, fresh per paper (per section above ~40 pages; sections from Marker's headings) | brief questions (`Q2@v1`), tag vocabulary, the paper's text, page images, tables, source disagreements; tiles and crops in its folder, read on demand for every equation, number and symbol it transcribes | `reads/<paper>.A.json`: claims (§7), each condition with its page; a symbol list (symbol, meaning, unit as printed, page); a section map | schema; one output per paper or section; the transcript audit lists the tiles read | failure rules below |
| E2 | Second read (core papers) | `morty-reader`, fresh | as E1 | `reads/<paper>.B.json` | used by E3 | — |
| E3 | Check extraction | `check_extraction.py` | reads, text | quote checks; page coverage; table values; importance (§5); A and B merged (§5) | quote ≥ 8 words or a whole sentence, unique in the paper, and found on the cited page. For a PDF the page text is Marker's, or the text layer used only to locate the quote. For LaTeX or XML, `pages/p001.md` is the source text split at the page boundaries found by aligning it with the PDF's text, and the quote is matched there; an alignment failure is its own flag. A value cited from a table equals that cell | failure → E4, marked; a quote that cannot be placed on a page → its own flag |
| E3b | Check units | `check_units.py` | symbol lists, equations, parameter fields, `memory/semantic/units.md` | per unit: parses or not; per symbol: dimension against the units table; per equation: homogeneous, inconsistent, or not parsed; dimensions inferred for symbols the paper leaves without a unit | Pint and SymPy, with an operator layer for ∇, ∂ and · (§5) | inconsistent or unparsable unit → E4 with the failing term; *not parsed* → X1 and E4, not a units pass, counted per paper |
| E4 | Check claims | `beth-checker` (task claim), per group of 3 to 5 pages, in two steps | **step 1, blind:** for each equation and parameter claim, only the quantity or equation name and its page and position; **step 2:** the claims citing those pages (their fields, never the reader's reasoning or Marker text); in both steps the pages' images, tiles and crops; the abstract, model, method and conclusions pages, loaded once; each page a claim lists a condition on; E3 and E3b flags and the positions of source disagreements (not what each source read); the brief's questions | step 1: the equation or value and unit as printed, transcribed before the claim is seen; step 2, per claim: text supported / overstated / unsupported / not on the page; quote as printed; every equation symbol and number as printed; kind, basis, questions confirmed or corrected; dropped conditions | schema; step 1 is returned before step 2 is sent, and the transcript audit confirms the order | not fully confirmed → E5; a correction by E4 is a new claim version written by the checker, and goes through E3, E3b, X1 (if in X1's scope, §5) and a fresh E4 that sees only the corrected version |
| X1 | Independent blind read (§5: claims a memo cites, claims on model-critical questions, and a random 20 % of the other important claims) | `morty-blind-reader` (fresh per group) and `ocr_math.py` (a math OCR model, not a language model) | for each such claim: the tile a script picks from the claim's page and position (both tiles if it is near an edge) and the name of the quantity or equation, never the claimed value or any text of the paper | the equation in LaTeX, or the value and unit, as printed (one from each) | `compare_reads.py` compares the parameter fields, or the LaTeX after a token-level canonicalizer (§5), of the reader, the checker's blind transcription (E4 step 1), the blind reader and the OCR model. When the copy has LaTeX or XML source, the source is the reference instead of the OCR | a disagreement in content → U3, *disputed, important*; a difference only in formatting → E4; *not found* → E4 |
| E5 | Revise (≤ 2 rounds) | `morty-reader`, fresh | the claim, the checker's reasons, the cited and condition pages | revised claim → E3, E3b, X1 (if in X1's scope), then a fresh E4 | as E3, E4 | still failing → U3, *disputed* |
| U3a | Guided reading (core papers) | user, with `reading_view.py` | each core paper's key pages in order; a page's claims hidden until the page is read; quotes highlighted | per claim right / wrong / unsure; *missed* statements with page | census of every claim on a read page | wrong or unsure → E5, then back to the user; missed statements become claims through E3, E4 |
| U3b | Audit | user, same view | two strata: every claim cited by a model-critical memo statement (a random 60 if there are more), and a random 60 of the other passed claims not seen in U3a; every disputed claim; every *disputed, important* claim | right / wrong | no error in either stratum; a bound reported for each (§5) | an error → E4 prompt revised, C0 rerun with it, the wave's claims rechecked under the new prompt, a fresh 60 drawn from those; the combined bound over all samples is reported |
| L1 | Publish | `curator.py`, the only writer | items with passing verdicts or a user decision | records, catalog rows, events (§6) | refuses items without citation, quote match and verdict; stages, moves, sets read-only, appends the event last | refusal → logged, shown at U4; crash → `curator.py recover`; after each publish, new and changed files go to the archive store (§6) |
| L2 | Audit | `audit.py` | library, index, archive | report | rehash files; rebuild catalog and index and compare canonical dumps; requery retractions; rerun `check_citations.py` on every memo; check the latest backup and archive against their hashes | mismatch → stop; new retraction → paper flagged, user told, and its claims blocked from Y1 and failed by Y2 until the user rules; memo citing a superseded claim → Y1 |
| Y1 | Synthesize | `beth-synthesizer` | the question; every verified claim for it, by fixed query (question ID; tags through aliases), excluding claims on papers flagged as retracted; split in parts above a token limit | `memos/<review>/Qk.vN.md`: statements citing `[id@v1]`; every number tagged with its claim and quantity (§5); INFERENCE marked with its claims; model options, *library silent* and gaps marked as such | Y2, Y3; script: every claim reached some part | — |
| Y2 | Check citations | `check_citations.py` | memo, library | per sentence: cites ≥ 1 current verified claim on a paper not flagged as retracted (sentences marked as model options, *library silent* or gaps are exempt); every number's tag names a cited claim, and its value, unit and quantity match that claim's parameter field (§5) | script | failure → Y1 |
| Y3 | Check memo | `beth-memo-checker`; for model-critical questions also `beth-skeptic` (Opus, told to find the strongest reason each statement is wrong); two calls each | **call 1, support:** memo and the cited claims only; **call 2, omissions:** memo, the same claim set, and claims and passages a broader search found outside it | call 1: per statement supported / overstated / unsupported; numbers on the right quantity; inference follows. Call 2: omitted or contradicting claims and passages; gaps. The schema has no field for support from a passage | schema; for model-critical questions the two checkers' verdicts compared by script | not supported, or the checkers disagree → Y1 revises once, rest to U4; a relevant passage with no claim → its page goes to E1 (task page), then E3–E4, before any memo may cite it |
| D2g | Check the gaps | D1 queries, P2 search, then `morty-scout` (task challenge, as D2c) | each gap the memo lists | the queries run and their counts; citations that may fill it | as D2c | a memo states a gap only as "not found by queries …", listing them |
| R1 | Retag for a new question | `retag_candidates.py`, then `morty-reader` (task retag), then `beth-checker` | the new question; the claims a script prefilters by tags, aliases and search score (top N), plus the claims on the pages a passage search returns | proposed `questions[]` additions and pages to read | E4 confirms | before Y1 for that question |
| U4 | Decide | user | memos, Y3 report, gaps with their queries and challenges, refusals, Zotero failures, screening baseline agreement | close, or expand / narrow (new brief version); accept or mark each unsupported memo statement; add or rename tags | tag changes written by `curator.py` (new `vocab.yaml` version, alias map) | — |

### 13.3 Page tiles as planned in revision 9

**Page tiles (P1).** Claude 4.7 and later models scale an image down to at
most 2,576 px on the long edge, and a request holding more than 20 images
caps every image at 2,000 px (Claude docs, "Vision", read 2026-10-23 and
2026-10-28). Agents read more than 20 images per call, so every image an
agent is given is at most 2,000 px on its long edge: a full page is
downscaled for agents to about 155 dpi for a letter page (the stored page
image stays at 300 dpi). At that size 7 pt subscripts are about 9 px high and 5 pt
second-level scripts about 6 px, where symbol swaps happen. So each page is
also cut into 300 dpi tiles:
- **Size:** about half a column high (a half page for single-column layouts).
- **Cuts:** at gaps between text blocks where possible, otherwise with a 1 inch
  overlap, so no line is cut in two.
- **Limits:** each tile at most 2,000 px on its long edge (a 3.4 × 5.5 inch
  tile is 1020 × 1650 px).
- **Use:** equations, numbers and quotes are read and checked on tiles. T-res
  measures the gain.

### 13.4 Best source first, source comparison and tables (revision 9, §3)

**Best source first (P1).** The text of a paper comes from the first of these
that exists:
1. the arXiv LaTeX source, when the copy read is that arXiv version (A2 says
   which);
2. the publisher's JATS XML, where it is openly available (for example
   Europe PMC);
3. the PDF through Marker.

Page images and tiles are always made, because the checker works only from
images.

Where two text sources exist, `compare_sources.py` aligns them page by page.
The sources are the source text, Marker's output and, for born-digital PDFs,
the PDF's own text layer. It lists every equation and number on which they
disagree. The list goes to E1 with the readings, and to E4 as positions
only (page, place, kind), as places to look closely. The text
layer is never given to a reader as text, because it loses superscripts and
signs. It is used only in this comparison.

**Tables (P1).** Each table found by Marker's layout model is handled as
follows:
1. **Crop:** the table is cropped at 300 dpi.
2. **Extract twice:**
   - from the source text when there is one;
   - otherwise with Marker and with a second extractor, chosen by T-table.
3. **Store:** `content/<copy>/tables/T<k>.json` holds the caption, column headers
   with their units, row headers (with group headers carried down to every
   row, so a row read alone still names its group), cells, footnotes, page,
   and which extractor gave each cell.
4. **Compare:** the two extractions are compared cell by cell, and every
   disagreement is flagged.

A reader citing a table value gives the table, row and column. E3 checks the
value against that cell. E4 checks the cell against the crop when the
extractions disagreed or the claim is important.

### 13.5 Where files come from, and the OpenAlex key (revision 9, §3)

**Where files come from (A1), in order:**
1. the user's `Literature/` folder and the practice library;
2. the open-access location OpenAlex records (which includes Unpaywall's
   data);
3. arXiv (PDF and LaTeX source);
4. the Semantic Scholar open-access link;
5. author and institutional repositories.

Anything else goes on the fetch list (U2). No institutional library access
for now, and no Sci-Hub (§9).

**OpenAlex key.** The key is in the macOS Keychain (item
`openalex-api-key`), and scripts read it at run time. It is never written to
the repository, a manifest, a log or a transcript, and it is removed from
recorded URLs.

### 13.6 The reference desk (revision 9, §3)

**The reference desk (Q1).** A skill, `library-desk`, for the user's ad hoc
questions. It is not a separate agent, and its text never enters the main
session (user decision 2026-10-28, F3 of round 4).
- **How it answers:** the main session passes the question to a headless
  call with no network and no tools beyond `library_query.py passages`. The
  call writes the hits verbatim (chunk, page, quote), with an optional
  summary labeled as a summary, to a file the user opens in the reading
  view. The main session gets only the chunk IDs, pages and an abstain flag,
  so it cannot discuss passage text.
- **Check:** `check_citations.py` confirms each quote is in the chunk it
  cites.
- **Abstaining:** an honest *abstain* when the index has nothing is a valid
  answer.
- **Labeling:** answers are labeled *passage, not checked by E4*.
- **Limits:** the desk is not on any path that makes a record or a memo.

This is also the course's Lab 5 retriever.

### 13.7 Review mode (revision 9, §4)

**Review mode: blocking the trifecta in the main session.** Every agent above
lacks at least one leg of the trifecta by construction. The main session
would have all three, so whenever it works with library or run files (during
a wave or not: the desk, U4, triage) it is started with
`claude --settings .claude/permissions/literature-review/sessions/review-mode.json`.
That settings file removes two
legs by configuration, not by instruction.
- **No way to send data out:**
  - WebFetch and WebSearch are denied;
  - there are no MCP servers (`--strict-mcp-config` with none);
  - `curl`, `wget` and `git push` are denied;
  - Bash network goes only to the hosts the discovery scripts need
    (sandbox allow-list).
- **No untrusted text:** reads work by allow-list, not deny-list. The
  session may read only the drivers' summaries, the scripts, the settings,
  the briefs and this file; every other path is denied. That covers claims,
  paper content and sources, `Literature/incoming/`, candidate and shortlist
  files, the fetch list, transcripts, `runs.db`, the event log, views,
  `references.bib`, history, memos and the calibration sets, all of which
  hold text from papers, APIs or the web. This takes two mechanisms, because the sandbox covers only
  shell commands, and the Read tool follows permission rules instead (Claude
  Code docs, sandboxing, read 2026-10-27):
  - Read permission rules for the Read tool;
  - `sandbox.filesystem.denyRead` for Bash.

  Drivers write a summary of counts, IDs, enumerated verdicts and hashes,
  validated against a schema with no free-text fields, and that summary is
  all the main session reads. The user sees the text in the reading view.
- **Less private data:** `Bash(security *)` is denied, so the session cannot
  read Keychain keys; scripts read them themselves.

What remains:
- **Queries to the allowed API hosts.** Scripts build them from the brief and
  from DOIs checked against the DOI pattern, never from free text.
- **Pasted text.** The user may paste paper text into the chat.
- **Settings edits.** The settings file itself is protected only by
  approval; edits to `.claude/` ask.
- **Sessions outside review mode.** The project settings deny their reads
  of the pipeline's stores (above), so this is now a setting. What remains
  is the user's own `Literature/` folder, readable in normal sessions by
  choice, and Python scripts run outside the sandbox, which deny rules do
  not cover.

T-trifecta tests the profile.

### 13.8 Agent independence and roles as planned in revision 9 (§4)

**Why split this way.** A second agent is used only for one of three reasons:
1. **Different information:** the checker must not see the reader's reasoning.
2. **Different permissions:** the scout alone has the web.
3. **A second opinion on the same object:** the second reader, the memo
   checker, the skeptical memo checker, and the blind reader.

Everything else is a script or a rule. Round 3 replaced three agent steps:
- screening adjudication became a rule in SC2;
- tag proposals were dropped with the tag growth process (fixed vocabulary, §6);
- the reference desk became a skill. Finding and downloading papers is a script
because it is deterministic. Keeping it apart from reading puts an identity
check between them, keeps the reader off the web, and wastes no reading
context on failed downloads. The library writer is a script because an agent
can be talked out of a rule by text it reads, and a script cannot.

**How independent the checks are.**
- **Claude models:** Sonnet, Haiku and Opus share a vendor and much of their
  training. The checker's main independence comes from its different input:
  images only, without the reader's text. For equations and values it
  transcribes blind before it sees the reader's fields, because a checker
  shown a transcription tends to confirm it.
- **Outside the Claude models:** the math OCR model and, where it exists,
  the authors' source give reads of equations and values that share none of
  Claude's training.
- **Shared mistakes:** three models agreeing can still be wrong, so the random
  U3b sample stays.
- **The truly independent checks** are the scripts and the user.

### 13.9 Checks of revision 9 (§5): the table, X1 scope, units, merging, title similarity, calibration and samples

| Output | Check | Kind | Pass |
|---|---|---|---|
| candidates | resolve in a bibliographic source (D3) | deterministic | exact identifier, or title similarity ≥ 0.9 with first author |
| candidates | retraction and correction status, Crossref `updated-by` and `filter=updates:` (D3, L2) | deterministic | recorded; shown at U1 |
| discovery | every query answered | deterministic | else *discovery incomplete* |
| screens | SC0 rules; two framings combined by rule (SC2) | deterministic + agents | rule applied; any include kept |
| screens | U1 random 60 exclusions | human, measured | no relevant paper in the sample |
| screening agents | agreement with the SC0 baseline on the U1 decisions | measured | reported each wave; if the agents add nothing over the baseline for two waves, U4 decides whether to drop them |
| copies | page 1 identity; A2 comparison | deterministic + user | match, or user accepts |
| text | source comparison (P1) | deterministic | disagreements listed and passed to E1, E4 |
| tables | two extractions, cell by cell | deterministic | disagreements flagged; cited cells checked on the crop |
| quotes | length, presence on the cited page of the imaged PDF, uniqueness | deterministic | all three |
| table values | cited value equals the cell | deterministic | equal |
| units | Pint parse; dimension against the units table; equation homogeneity (E3b) | deterministic | no inconsistency, or E4 resolves it; *not parsed* is not a pass |
| claims | `beth-checker` on images | agent, different input | all confirmed |
| checker corrections | a corrected claim goes through E3, E3b, X1 and a fresh E4 | deterministic + agents | as for any claim |
| claims memos rely on | X1: the reader, the checker's blind transcription, the blind reader and the OCR model (or the source), on parameter fields or canonical LaTeX | independent reads, compared by script | all agree, else the user decides |
| other important claims | X1 on a random 20 % | measured | disagreement rate reported; reviewed after the shakedown wave |
| verdict repeatability | 30 random claims rechecked by a fresh E4, on the first two waves and after any prompt, model or tool version change | measured; repeatability only, not correctness | agreement ≥ 90 %; else the wave's passed claims are rechecked by E4 under a revised prompt once it passes C0 |
| recall | E1 vs E2 on core papers; U3a *missed* | measured, human | reported |
| core claims | U3a | human, census | every claim on read pages checked |
| claims cited by model-critical memo statements | U3b, all (or a random 60) | human, census or measured | no error |
| other claims | U3b random 60 | human, measured | no error |
| memos | Y2 citations, numbers and units, rerun at every L2 | deterministic | 100 % |
| memos | Y3 (with `beth-skeptic` for model-critical questions); D2g on gaps | agents, two framings | no unsupported statement left, or user accepts each |
| search | T-ret | measured | reported; picks the search mode and the embedding model |
| library, index, backup | rehash; rebuild and compare (L2) | deterministic | no mismatch |

**Important claims.** A claim is *important* if any of these holds:
- its kind is parameter or equation;
- by script, whatever its kind: it has a number with a unit, a table cell,
  or LaTeX;
- it bears on a question the brief marks *model-critical*;
- the user marks it.

Importance is recomputed after E4, and a claim that becomes important goes
through E3b.

What important claims get:
- **units:** all go through E3b;
- **Independent read (X1)** (user decision 2026-10-28, F13 of round 4): in B1 and
  B5 nearly every claim has LaTeX or a number, so X1 on all of them would
  grow with every paper. X1 runs on:
  - every claim a memo cites (before the memo goes to Y2);
  - every claim on a model-critical question;
  - a random 20 % of the other important claims, to measure what the rule
    skips.

  The disagreement rate on the 20 % is reviewed after the shakedown wave.
- **review:** every disputed one goes to the user;
- **user time:** it grows with the claims that are both disputed and relied
  on, not with the number of papers.

**Comparing LaTeX (X1).** `compare_reads.py` canonicalizes LaTeX at the
token level before comparing: it drops spacing and sizing commands
(`\left`, `\right`, `\,`), `\mathrm` and similar font commands, and
optional braces; it treats `\cdot` and juxtaposition alike where both mean a
product; and it sorts sub- and superscript order. Every equivalence is in a
tested list. T-res measures the false-dispute rate. A difference left only
in formatting goes to E4, not to the user.

**Units (E3b).** Pint is a Python units library; SymPy is a Python algebra
library. For each paper:
1. **Parse:** every printed unit is parsed with Pint. One that does not parse
   is a flag (for example `m-1` read for `m^-1`).
2. **Compare with the units table:** a symbol whose meaning matches a row of
   the project's verified units table must have that row's dimension.
3. **Check equations:** an operator layer first replaces ∇, ∂ and · with
   their dimensional effect (1/length, 1/time), since SymPy's LaTeX parser
   does not handle vector calculus. Each equation that then parses is
   checked term by term for the same dimension. A symbol the paper gives no
   unit for gets its dimension from the equation, and the inferred units are
   listed for E4. Inferred units are never used as reference values; the
   only reference is the verified units table (`memory/semantic/units.md`,
   54 rows from the units-verification workflow), whose hash is in the
   manifest. E3b refuses to run if the hash differs.
4. **Nondimensional equations:** the reader marks them, and they are checked
   only for being dimensionless.

Equations that do not parse are counted as *not checked* and reported per
paper, so coverage is visible. They go to X1 and E4, and they have not passed
the units check. T-units measures the coverage on the B1 seeds. A lost exponent or a swapped symbol usually
breaks homogeneity; a wrong number with the right unit does not, which is
what X1 is for.

**Merging readers A and B (E3).** Two claims are the same when all four of
these hold:
- they cite the same page;
- they have the same kind;
- their quotes overlap by at least half of the shorter quote's words;
- their claim texts are similar, or their parameter fields match.

What happens to the claims:
- **A matched pair** becomes one claim. E4 confirms the wording the page
  supports. The other wording is kept as a variant, which is never shown to
  Y1.
- **Every unmatched claim** is kept and checked by E4.
- **Reader omission rate:** the share of merged, confirmed claims that the
  reader did not find.
- **Second pass:** likely duplicates are found by statement similarity and
  go to E4.

**Title similarity (D3).** `difflib.SequenceMatcher` ratio of the two titles
after normalization (NFKC, lower case, punctuation removed, whitespace
collapsed). A ratio ≥ 0.9 matches, together with the first author's surname.

**C0 calibration.** Two parts.
- **Core, run once.** It tests the mechanical checks, which do not depend on
  the field.
  - **Plants:** planted errors in 2 to 3 papers.
  - **Who writes them:**
    - at least 10 written by the user;
    - more drafted by a fresh Opus call (`beth-planter`), each approved by
      the user, who rules on every draft against the page;
    - mechanical errors made by script from the user's gold statements
      (which the user has verified) and later from verified claims: a changed
      quote, a shifted page, a flipped sign or exponent, a changed digit.
  - **Copy kinds:** every plant records its copy kind (born-digital PDF,
    LaTeX or XML, scan). The core set includes at least one scanned paper,
    with script and model-drafted plants approved by the user.
  - **Controls:** at least 10 untouched claims.
  - **Error kinds:**
    - a quote not in the paper;
    - a wrong page;
    - a changed sign, exponent or subscript;
    - an OCR misreading copied into an equation;
    - a number on the wrong quantity;
    - a wrong unit;
    - a wrong table cell;
    - a hedge made certain;
    - a wrong basis, kind or question;
    - a dropped condition;
    - a wrong figure;
    - a retracted paper.
  - **Tuning and held-out sets** (user decision 2026-10-27). The plants are
    split at random into a tuning set, used when prompts are revised, and a
    held-out set, used only for the pass. The user's own plants are all held
    out. The held-out size is fixed before the first run.
  - **Using the held-out set once.** A held-out set is scored once per prompt
    version, and its misses are never shown to whoever revises the prompts
    (the misses are reported as counts per kind). After a failure, the next
    prompt version is scored on a fresh held-out set drawn from the growing
    pool, and the number of attempts is reported with the bound.
  - **Pass, on the held-out set:**
    - **kinds E3 catches by rule** (fabricated quote, wrong page, wrong table
      value): all caught. These measure the scripts, so they are not pooled
      with the rest;
    - **the plants that reach E4:** a one-sided 95 % Clopper–Pearson lower
      bound of at least 85 %. That needs at least 19 such plants with none
      missed (19 of 19 gives 85.4 %), or 60 with two missed (89.9 %);
    - every sign or exponent change caught;
    - at most 10 % of controls flagged.
  - **Reporting:**
    - each kind's catch rate, with its interval. A kind with 8 or more
      held-out items must also meet the 85 % bound by itself; the growing
      set provides the items over time;
    - the catch rate on the user's plants and on the model-drafted plants,
      separately. A Claude model may draft errors that Claude models find
      easy to catch. If the user's plants are caught clearly less often, the
      drafted plants are too easy, and the rate on the user's plants is the
      one that counts.
  - T-res (§12) is part of the core.
- **Branch, light.** It runs when a branch is seeded, on the one or two seeds
  the user knows well, and tests the judgments that depend on the field.
  - **Gold statements:** about 10 per paper, in the user's words, one fact
    each, with page (format in `calibration/gold_statements.yaml`). They also
    serve as T-ret questions.
  - **Planted errors:** about 10 of the field-dependent kinds (dropped
    condition, overstatement, wrong basis), drafted by the user or by
    `beth-planter` and approved by the user, plus plants of the branch's main copy kind if
    the core set lacks it (for example, scans in B3 or B5).
  - **Pass:**
    - E1 finds at least 90 % of the gold statements, and E1 and E2 together
      find all of them;
    - the planted errors are *reported* (caught or not, with the interval),
      and not called a pass: 10 plants cannot meet the 85 % bound (10 of 10
      gives 74 %). The gate is the core calibration. A branch's planted
      errors become a gate once its held-out set has 19 or more.
- **Growing set.** Real errors are added with their page, so the set grows:
  - errors found independently of E4 (X1 disagreements the user resolved,
    U3a, U3b, Y3) go to the held-out pool;
  - errors E4 found go to the tuning set only, since E4 already catches them
    and they would inflate its measured rate.
- **Screening.** Seeds must not be excluded. Borderline abstracts (10 to 20)
  are added for large branches (B4) when they are first screened.
- **Reruns.** C0 is rerun:
  - after any change to an agent file, prompt, schema or model ID;
  - after a change to the versions, weights or flags of Marker, Surya, mutool,
    the table extractor, the math OCR model or the embedding model;
  - after errors found at U3b.

  A Claude Code version change is made only between waves (it is pinned),
  and runs a 10-item subset first, and the full C0 if it fails.

**Sample bounds (U1, U3b).** Both samples are simple random samples drawn by
script, with the seed recorded. For k errors in n = 60, the one-sided 95 %
Clopper–Pearson bound on the error rate is:
- k = 0: 4.9 %;
- k = 1: 7.7 %;
- k = 2: 10.1 %.

U3b passes only with k = 0 in each stratum. When a stratum has 60 items or
fewer, all are inspected. The first stratum (claims cited by model-critical
memo statements) is the one the model rests on, so a uniform sample of all
claims would hardly touch it.

### 13.10 The passage index, search and archive (revision 9, §6)

- **Passage index (P2).** Built from the `content/` folders, so it can be
  rebuilt from the library alone. Chunks follow the paper's structure, not a
  character count. A chunk is one of:
  - a paragraph, with its section heading;
  - a displayed equation, with the sentence before and after;
  - a table: caption, headers, rows and footnotes; a large table is split by
    row groups, with headers repeated;
  - a figure caption.

  How chunks are stored and searched:
  - **Header:** each chunk begins with a fixed header (short citation,
    section title), so it still says where it is from when read alone. This
    is a deterministic version of contextual retrieval (Anthropic 2024).
  - **Metadata:** paper, copy hash, page, section, chunk kind and position.
  - **Neighbors:** linked, so a hit can be returned with its neighbors.
  - **Keyword search:** SQLite FTS5 (BM25 ranking).
  - **Meaning search:** an embedding from a pinned local model (§8). Search
    compares the query with every chunk's vector, so it is exact and
    repeatable.
  - **Merging:** the two ranked lists are merged by reciprocal rank fusion
    (k = 60).

  Why both: keyword search finds exact symbols and identifiers that
  embeddings blur, and embedding search finds paraphrases that keyword search
  misses.

  Where the index is used:
  - Y3 and D2g (evidence a memo omits);
  - R1 (pages to reread for a new question);
  - D2g (queries for each gap);
  - the reference desk;
  - SC0, which uses the same embedding model on a separate store of candidate
    abstracts.

  It never feeds a memo directly. A passage becomes citable only as a claim,
  through E1 to E4.
- **Search.** `library_query.py` is the one search script:
  - subcommands `claims` (filters: branch, tag, paper, question, kind,
    basis, status), `passages` (hybrid search) and `chunk` (one chunk with
    its neighbors);
  - every query is logged with the asker and purpose (its only write).
- **Archive and backup: one content-addressed store** in
  `~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Courses/Fall2026/ME295/library_archive/`
  (user decision 2026-10-27; merged with the backup in revision 8).
  - **What goes in:** everything a recheck needs that git does not hold.
    Each file is stored once under its hash: each paper copy's `source/` and
    `content/` (images included) when it is published, each new record
    version, and each run's transcripts. A file already stored is never
    written again, so the store grows with what is new, not with the number
    of waves.
  - **Manifests:** after each publish and each wave, a dated manifest lists
    the hash of every file the library holds at that moment, so any past
    state can be rebuilt from the store.
  - **Checks:** L2 checks the latest manifest against the store, downloading
    files first with `brctl download` if iCloud has evicted them.
  - **Limit:** iCloud Drive syncs; it is not a versioned backup. Files are
    never overwritten or deleted, which is what makes it usable as one.
  - **Never in git:** it holds paper text.
  - **Size:** a few GB over about 100 papers, written once.
### 13.11 Deferred tests (revision 9, §12)

| Test | What it checks | Pass |
|---|---|---|

| T-table | 5 or more tables of different layouts from seed papers, extracted by each candidate, compared with the user's check of the crop | the chosen extractor has no cell errors on born-digital tables; errors on scans listed |
| T-ret | per branch, 25 questions written by the user from use, before looking at chunks: about 10 exact (symbols, values, named equations), 10 paraphrases and 5 not in the library; gold statements count. Hit@5 for keyword, embedding and hybrid search, and abstains on the 5 | reported; picks the search mode and the model. Misses cost recall, not correctness, so there is no fixed pass |
| T-blind | the blind reader and the math OCR candidates on the T-res tiles, given only the quantity names | the OCR model is chosen by its symbol error rate; both reads reported per symbol |
| T-retract | Crossref finds the retraction of 10.1021/am300292v | found |
| T-key | the Keychain key is read in a headless run; a failed read stops D1 | as stated |
| T-backup | a file stored once under its hash, a manifest written, files downloaded after eviction and hash-checked; a second publish writes only new files | match |
| T-zotero | on a test collection: one new paper with its PDF, one already in the library, one without a PDF that the user then saves with the Connector; then a rerun | new item and file added (hash matches); the existing item is not duplicated or edited; the paper without a PDF is not added, and `zotero_pull.py` copies the user's save; the rerun adds nothing |
