# Agentic workflows: what each agent does, how it is specialised, and how the pieces fit

Updated 2026-10-19. Agent instruction files are in `.claude/agents/`; tools and pipelines are listed in `agentic/tools/REGISTRY.md`; the library rules are in `agentic/library/POLICY.md`. Older design notes (rationale, course links) are in the research folder's `Planning/` (02, 07, 09, 10, 11, 15, 20, 21).

## 1. How the agents are built and specialised
- **One narrow job per agent**, written as one instruction file (name, model tier, description, rules). The instruction states exactly what the agent must return.
- **Fixed output shape.** Most agents return JSON checked against a schema, with a short reasoning field first (at most about 80 words). A wrong shape is rejected and retried once.
- **Fresh context, minimal input.** A worker is started per task with only the text it needs; it has no memory of earlier work. Headless workers run with no tools (`claude -p --tools ""`), so they cannot read or change files; the image checker alone may read one page image.
- **Model tier by role.** Haiku scores short items (triage). Sonnet does the bulk work (reading papers, drafting derivations and reports). Opus judges other agents' work (claim checks, image checks, adversary, red teams, critics, auditing) because a different, stronger tier should check the author. All are Claude models, so errors can correlate; the programs below catch what models could share.
- **Deterministic checks the agents cannot override**: quote matching against the page, a symbolic checker, numerical solutions, a fact ledger that fills every number in a report, the library's hash check.
- **Human gates**: you approve the reading list, read flagged claims and memos, and read each derivation before it is relied on. No agent submits a job, merges code or edits the library.

## 2. Roster (18 agents)
| Agent | Tier | Workflow | What it does | Specialised by |
|---|---|---|---|---|
| query-planner | Sonnet | literature | turns review questions into search queries | outputs queries only |
| triage-scorer | Haiku | literature | scores one candidate from title and abstract, 0 to 3 per topic R1 to R11 | sees no full text |
| paper-extractor | Sonnet | literature | reads one whole paper, returns about 30 claims each with page and exact quote, answers topic questions, lists equations | one paper per worker, quote and page required |
| claim-checker | Opus | literature | judges whether one claim is supported by its quote and the cited page | sees only claim, quote, page |
| image-checker | Opus | literature | reads a page image to confirm or correct equations and numbers the text layer garbled | may read only that page image |
| topic-synthesizer | Opus | literature | writes a gate memo for one topic from checked records only | records via the desk, never the PDFs |
| librarian-curator | Sonnet | library | reviews metadata and extraction before publication; proposes tags, duplicates, hold | proposes only; a program writes |
| librarian-reference | Sonnet | library | answers requests from the passages the desk retrieved, with citations; says "not found" | cannot search or write |
| librarian-auditor | Opus | library | explains the integrity and access report in plain language | cannot change anything |
| math-spec-writer | Sonnet | derivation | states the problem: equations, assumptions, expansion order, what to match | does not derive |
| math-deriver | Sonnet | derivation | derives step by step, each step carrying a machine-checkable claim | must reference trusted definitions (@name) |
| math-adversary | Opus | derivation | sees only spec and derivation; looks for missing terms, hidden assumptions, wrong limits | does not see the checker's report |
| math-scribe | Sonnet | derivation | typesets an approved derivation, adds nothing | runs only after the gate |
| formulation-reviewer | Opus | review | fresh red team of a derivation, plan and memo (physics, numerics or logic role) | no stake in the work |
| report-outliner | Sonnet | writing | plans a report from rubric, ledger and brief | never writes prose |
| report-drafter | Sonnet | writing | writes the report in the author's voice; numbers only as ledger placeholders | style profile and ledger rules |
| report-critic | Opus | writing | flags unclear, jargon-heavy or hollow passages and AI-sounding prose | sees only the rendered draft |
| report-rubric-auditor | Opus | writing | checks the draft against the course guideline and lists what only the author can confirm | guideline text |
(The four math agents live on branch `agent/math-lane` until it is merged; `formulation-reviewer` is on both.)

## 3. The workflows
**Literature.** query-planner and discovery (OpenAlex) -> triage-scorer -> **you approve the reading list** -> PDFs fetched (or you download paywalled ones) -> text and page images prepared -> one paper-extractor per paper -> program checks every quote on its page -> claim-checker on all number/equation claims plus a quarter of the rest -> image-checker on flagged pages -> a second independent reading of about a fifth of the papers -> records published to the library -> topic-synthesizer writes a gate memo from library records -> citation audit of the memo (claim-checker lines) -> **you read flagged items**. Planned: a panel of several different checkers with an adjudicator.
**Derivation.** math-spec-writer -> math-deriver -> symbolic checker recomputes every step and requires the final result to appear in a passed check -> math-adversary -> optional blind second derivation by a different tier, compared by the symbolic checker -> numerical audit against an independent solution -> **you read the summary** -> math-scribe. Known-answer tasks and planted errors calibrate the workflow before it is used on new physics.
**Review.** formulation-reviewer (physics, numerics, logic roles) reads the package; every objection carries a test; we run what we can and record the outcome.
**Writing.** evidence ledger (facts with sources) -> report-outliner -> report-drafter -> deterministic checks (length, long sentences, raw numbers, unexplained jargon, voice, sections) -> report-critic and report-rubric-auditor -> revise, keeping the best version -> renderer (Word, PDF) -> a "before you submit" checklist for you.
**Library.** The curator program publishes records, texts, memos, search sets and the candidate table (hashed, versioned, read-only) after the intake checks; agents read only through the reference desk, which needs a stated purpose, applies a who-may-ask-for-what table, verifies hashes and logs every request; the auditor verifies everything and explains problems to you.

## 4. Where things live
- Agent definitions `.claude/agents/`; code `agentic/lit`, `agentic/math`, `agentic/writing`, `agentic/library`, `agentic/tools`; tests `agentic/tests`.
- Workspace (rebuildable): `.agent_state/` (registry, extracted text, token ledger); working outputs in `research/Review/`, `research/Literature/`.
- **Library (the record of what we know):** `research/Library/` (iCloud, not in git).
- Run-by-run costs and calls: `.agent_state/ledger.jsonl` (tier, tokens, time per call).

## 5. Honest limits
All models are one family; one checker per claim so far (the panel is the planned fix); model versions are not pinned in the ledger; library protection is strong against mistakes and against agents following instructions, not against a deliberate shell command by the same user (see POLICY.md); headless workers cannot write, but a driver bug could.
