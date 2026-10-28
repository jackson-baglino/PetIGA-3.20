# Retrospective: struggles, pivots and redos

For the final presentation's "what went well, what didn't". Newest first. One
entry per pivot or redo:
- **What happened:** what we tried and why it failed or changed.
- **Pivot:** what we did instead.
- **Lesson:** the general rule we took from it.
- **Evidence:** a commit, run or file.

Added to whenever a plan is abandoned, a result is redone, or a check catches
our own mistake.

## Pivots and redos

### 2026-10-28 — Two permission mechanisms in the design did not work as written
- **What happened:** building the permission files showed:
  - **Carve-outs:** revision 9's review mode was "one deny list with `!`
    carve-outs". The docs say a carve-out cannot reach an anchored rule, and
    cannot reopen a file inside a blocked directory, so that allow-list
    would not have worked.
  - **`Write(...)` rules:** the first live test printed warnings that
    path-scoped `Write(...)` deny rules are ignored. Only `Edit(...)` covers
    file edits.
  - **The scout:** the trifecta check found that the scout's inputs (the
    brief's questions) plus its web search made all three legs.
- **Pivot:**
  - review mode uses deny rules, a PreToolUse hook as the allow-list, and
    the sandbox;
  - the `Write(...)` rules were removed;
  - the scout gets only search terms the user approves for sharing.
- **Lesson:**
  - Run a permission file against the real tool before trusting it.
  - Read the tool's warnings.
  - A check that computes the legs finds what prose review missed.
- **Evidence:** `.claude/permissions/literature-review/README.md`;
  `scripts/literature_review/check_permissions.py`.

### 2026-10-28 — Codex dropped; permissions moved from prose to settings files
- **What happened:**
  - **Codex:** it was the cross-vendor reader, so that its misreadings would
    not line up with Claude's. Making it safe needed a separate macOS user,
    a shared binary and a sudoers rule. A configuration check also found
    that the user's own interactive Codex held all three trifecta legs:
    browser and computer-use plugins, web search on, and home-folder
    access. Its benefit had never been measured, and its usage limits could
    stall a wave.
  - **Permissions:** they were described in tables and prose in
    WORKFLOW.md, enforced only where a matching setting happened to exist.
    A proposed YAML file would have enforced nothing by itself.
  - **Caltech Library downloads:** the daily limit is unknown, so scripted
    downloads were deferred.
- **Pivot:**
  - **Independence without Codex:** X1 uses a fresh Claude blind read plus a
    math OCR model, with author sources as the reference.
  - **Permissions:** written as the `settings.json` files Claude Code
    enforces, with a script that checks them for the trifecta and writes
    the tables.
  - **Zotero:** taken up again, with the user's own saves pulled back in.
- **Lesson:**
  - Weigh an unmeasured safeguard against its setup risk.
  - A permission counts only if a file the tool enforces says so.
  - Check the configuration of every tool, not only its documentation.
- **Evidence:** WORKFLOW.md revision 9, §4 and §9 items 34 to 39.

### 2026-10-28 — Fourth critique: our statistics and our security had holes
- **What happened:** round 4 found 19 major findings, among them:
  - **Calibration:** the branch calibration could never pass (10 of 10
    caught gives a 74 % lower bound, below the 85 % rule). The pooled score
    mostly measured the scripts, not the checker. The held-out set was
    reused after each prompt fix, and errors the checker found were fed back
    to inflate its own score.
  - **Review mode:** it blocked reads by deny-list and missed many folders
    holding paper or API text. It also applied only during waves, so the
    reference desk handed passages to a session with web access.
  - **Independence:** the checker was shown the reader's transcription, so
    the reader, the checker and Codex agreeing could still rest on one
    Marker misread.
- **Pivot:** revision 8. Gates are split by what they measure, and each
  held-out set is used once. Review mode works by allow-list whenever a
  session touches library files, and the desk writes to the reading view.
  The checker transcribes blind first. One finding (image limits) was
  refuted against the docs, and checking it turned up a limit we had missed
  (more than 20 images caps each at 2,000 px).
- **Lesson:**
  - Compute what a pass rule can actually reach before adopting it.
  - Allow-lists, not deny-lists, for what a privileged session may read.
  - A checker shown the answer tends to confirm it.
- **Evidence:** `runs/design-critique/20261028-061227/triage.md`; WORKFLOW.md
  revision 8.

### 2026-10-28 — The library still stored only a description of each paper
- **What happened:** revision 6 added a full-text index, but the text itself
  stayed in P1's run folder and the archive. The library record of a paper
  was still its card, claims, symbols and tables, so the design summary read
  as a set of key facts, not the paper. The user caught it.
- **Pivot:** each paper's library folder holds the whole paper (source files,
  full text per page, equations, tables, figures, page images), written once
  and never edited. Claims are a checked layer on top, and the index is
  rebuilt from the library.
- **Lesson:** adding a capability is not the same as changing the record.
  Check that the summary a user reads says what the design stores.
- **Evidence:** WORKFLOW.md §1 and §6, edit of 2026-10-28.

### 2026-10-27 — Codex is an auto-updating alpha
- **What happened:** the Codex installed for cross-vendor checks is the
  command-line tool bundled inside the VS Code extension. It is an alpha
  (0.164.0-alpha.11.1) that changes whenever the extension updates, and its
  event log does not name the model it used.
- **Pivot:** a pinned standalone install, and the model named explicitly with
  `-m` (open).
- **Lesson:** "installed" is not "pinned". Check every tool for
  reproducibility before it enters a workflow.
- **Evidence:** WORKFLOW.md §9 item 25.

### 2026-10-27 — A user-facing file format broke on a colon
- **What happened:** a calibration statement containing ": " made
  `gold_statements.yaml` invalid YAML; the file would have failed when loaded.
- **Pivot:** the statement was quoted, and a quoting rule was added to the
  file's instructions.
- **Lesson:** validate files people fill in by hand as soon as they are
  written, not when the pipeline first reads them.
- **Evidence:** commit 4c9d7b61.

### 2026-10-27 — Tagging process dropped
- **What happened:** the design had a growing tag vocabulary: readers proposed
  tags, a tagging agent proposed merges, and an alias map handled synonyms.
  Once embedding search arrived, finding related material by meaning no
  longer needed tags.
- **Pivot:** a small fixed vocabulary, kept only for exact filters (branch,
  method family, material, physics). No tagging agent, and no growth
  process.
- **Lesson:** a new capability can make an older mechanism redundant.
  Re-ask what each part is for.
- **Evidence:** commit 8709b11d.

### 2026-10-27 — Third independent critique: 2 critical, 12 major findings
- **What happened:**
  - **Critical:** a checker's own corrections were published with no second
    check, so the checker became the author.
  - **Critical:** recorded results could not be rechecked after a disk loss,
    because page images, text and claims were not archived anywhere.
  - **Main-session exposure:** the main session held the "lethal trifecta".
  - **Sandbox:** the Claude Code sandbox fails open by default.
  - **Codex:** it searches the web by default and can read the whole disk.
- **Pivot:**
  - all 25 findings were applied (revision 7);
  - three agent steps became a rule, a rule and a skill: screening
    adjudication, tag proposals and the reference desk.
- **Lesson:** an independent critic asked "what could be a script?" finds
  agents that are pipeline stages in disguise (the lecture's "cover-up
  test"). Facts about tools must be checked in their docs; the critic's two
  web-sourced claims were confirmed.
- **Evidence:**
  - `runs/design-critique/20261027-170522/triage.md`
  - commits 6564797a and 384c57c5.

### 2026-10-26/27 — The library indexed only extracted claims
- **What happened:** the design kept only the claims readers extracted. A
  question nobody asked at reading time could not be answered later without
  rereading the paper. Lecture 8 (retrieval) and the user's notes made the
  gap plain.
- **Pivot:** a full-text passage index with keyword plus embedding search,
  structured tables, and best-source-first ingestion (revision 6).
  Passages never feed a memo directly.
- **Lesson:** decide at question time what matters, not at reading time,
  but keep the verified-claim layer between the text and anything written.
- **Evidence:** commit cd64e493.

### 2026-10-23 to 27 — Calibration size
- **What happened:** the first calibration plan was large: 40 or more
  planted errors, 30 gold statements, 10 to 20 borderline abstracts, all
  written by the user. That raised whether it must be redone for every
  subtopic, and whether another model could write it.
- **Pivot:**
  - **Core plus branch:** a core calibration (mechanical checks, run once)
    and a light one per branch (judgment that depends on the field).
  - **Codex drafts:** Codex drafts plants, and the user approves each one.
  - **Control:** the user's own plants are held out as a control on whether
    Codex's plants are too easy.
  - **Growing set:** real errors are added to the set as they are found.
- **Lesson:**
  - Separate what a test measures (pipeline mechanics, or field judgment)
    before deciding how often to repeat it.
  - Never let model-written answer keys go unreviewed.
- **Evidence:** WORKFLOW.md §5 (revisions 6 and 7).

### 2026-10-23 — Models downscale page images
- **What happened:** a full page sent to Claude is scaled to about 199 dpi,
  however finely it was rendered. At that size small subscripts are 6 to
  9 px, where symbol swaps happen.
- **Pivot:** 300 dpi tiles for every equation, number and quote, with a
  resolution test (T-res).
- **Lesson:** read the model's input limits; the rendering setting is not what
  the model sees.
- **Evidence:** WORKFLOW.md §3 "Page tiles".

### 2026-10-23 — Two design critiques before building
- **What happened:**
  - Round 1: 23 findings, all accepted.
  - Round 2: 22 findings, 21 accepted, and one refuted by a test.
  - Among them, the design named the wrong Crossref field for retractions
    (`update-to` instead of `updated-by`). Testing on a known retracted DOI
    confirmed the critic was right.
- **Pivot:** revisions 3 and 4.
- **Lesson:**
  - Critique the design before writing code.
  - Test a factual finding instead of arguing with it.
- **Evidence:** `runs/design-critique/20261023-103742/`, `20261023-111440/`.

### 2026-10-23 — Planned second critic unavailable
- **What happened:** the design used a Fable critic (`rick-critic`) for
  independence from the Opus models, but the account has no Fable access.
- **Pivot:** Opus critic only, with the weaker independence stated in the
  design. Later, Codex was added as a cross-vendor model.
- **Lesson:** check model access before designing around a model.
- **Evidence:** WORKFLOW.md §4, §11.

### 2026-10-23 — Headless agents silently refused web search
- **What happened:** in `--permission-mode dontAsk`, a tool that needs
  permission (WebSearch) is refused unless it is preapproved, so the critic
  ran without web access.
- **Pivot:** `call_agent(..., allow=...)` passes `--allowedTools`.
- **Lesson:** test what an unattended agent can actually do, not what its
  definition lists.
- **Evidence:** `scripts/units_verification/agents.py`.

### 2026-10-23 — Flowcharts unreadable
- **What happened:** low contrast, and edge labels drawn over nodes and other
  edges.
- **Pivot:** a theme file, plus a renderer that places labels itself and
  lints every overlap. A chart that fails the lint is not accepted.
- **Lesson:** make visual quality checkable by script.
- **Evidence:** `scripts/flowchart/`, `test_render_flowchart.py`.

### 2026-10-22 — The units table was unreliable
- **What happened:** the units registry carried over from the practice round
  had wrong or unsupported rows.
- **Pivot:** the units-verification workflow:
  - a blind read of each page image;
  - an OCR comparison by script;
  - an Opus adjudicator;
  - the user's inspection.

  54 rows were verified, and 6 symbols were dropped.
- **Lesson:** re-verify inherited "facts"; a table is only as good as its
  provenance.
- **Evidence:** run `20261022-155335`; commits cd59dcec and 2b9ed052.

### 2026-10-22 — Two formulas for one constant
- **What happened:** `comp_eps.py` computed d0 as γV_m/(RT), while the solver
  used γa³/(k_BT) (Kaempfer & Plapp Eq. 13). They differ by 0.14 %.
  Geometry files made earlier carry the old value.
- **Pivot:** one shared formula, computed once at start-up; the input
  option was removed.
- **Lesson:** one source for each number.
- **Evidence:** commit b927b982; `memory/MEMORY.md`.

### 2026-10-22 — Course repository history
- **What happened:** the course repository's history was unrelated to this
  project's, so it could not simply be mirrored.
- **Pivot:** one approved force-push of a subtree-split history. The old
  history is kept as a tag, and `scripts/sync_course.sh` mirrors from then on.
- **Evidence:** `project/LOG.md`, 2026-10-22 (3).

### 2026-10-21 — Practice round abandoned; restart from a clean slate
- **What happened:** the first attempt (2026-10-12 to 21) built agents
  without a written design.
  - **The literature fact-check:** it used three long-lived Opus reviewers
    that each read about 110 papers. They used about 272 million tokens, hit
    the usage limit three times, and left 60 to 70 papers per reviewer
    unread (audit in enceladus_DSM, 2026-10-15).
  - **The resulting library:** 85 papers and 2,532 claims, 66 papers still
    flagged, with one reader and one checker per paper and no calibration.
- **Pivot:**
  - everything was archived (`archive/2026-10-21/*` tags,
    `scratch/2026-10-21_practice_round/`);
  - every workflow is now designed with the user, with a flowchart, before
    it is built;
  - short-lived per-paper agents replace long-lived readers.
- **Lesson:**
  - **Design first.**
  - **Agent context is the cost driver.** Fresh, small contexts per item
    beat long-lived agents.
- **Evidence:** commit 6ae0b12d; `project/PILLARS.md`.

## What went well

- **Units verification.** The first workflow built under the new rules ran
  end to end. It produced 54 verified rows, and its calibration passed twice
  (`project/LOG.md`, 2026-10-22 (5)).
- **Design first, critique before code.** Four critique rounds found 95
  findings before any literature-review code existed, so none of them cost a
  rebuild.
- **Deterministic checks catch real errors.** Examples are the Crossref
  retraction field, the d0 mismatch, the YAML colon, and the sandbox and
  Codex defaults.
- **Retrieval and the course line up.** The library's passage index is also
  the course's Lab 5.
