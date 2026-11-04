# TODO: from today to a working library

Written 2026-11-03, from `STORY.md`, `SPEC.md` and the line's status. Tick items off here as they are done, with the
date and commit. Ideas that come up along the way go to `BACKLOG.md`, not here.

**How to read it**
- **Who:** **you** (only you can do it or decide it), **me** (a Claude session builds or runs it), **foreman** (the
  foreman session runs the line while you watch).
- **⟂ depends:** what to do next depends on how an earlier item turns out. The branch is spelled out.
- Each item says **why** it is there.
- Stage codes: v1 codes (U3, E4, …) until Phase 3 renames them; v2 codes (F1, R3, …) after.
- Dates are proposals. The one fixed point is that the methods memo should exist before model work for Lab Report
  3 (15 Nov) needs it.

**What "a fully working library" means here** (the finish line, Phase 6):
- **Layer 1:** every paragraph of every included paper of disciplines 1–4 is indexed and searchable.
- **Layer 2:** every fact of each discipline's core papers is verified and published.
- **Memos:** an approved memo for each discipline.
- **New questions:** a later question can be answered through the on-demand path (search layer 1, extract, check,
  publish).
- **Backup:** the restore test passes.

---

## Phase 0: already done (for orientation)

- [x] The v1 line is built: 17 stages, the stage runner, the foreman, hooks H1–H7, the dev run `dev-20261101-kp` end
      to end.
- [x] Build tests passed: T-det, T-sandbox, T-perm, T-trifecta, T-restore, T-applicability, worker isolation (auto
      memory fix, 61446b41).
- [x] Tuning calibration diagnostic (`calib-dev6`): 32/32 plants caught, **9/40 controls flagged** (rule: at most
      10 %, which is 4 of 40).
- [x] Marker removed; `venv_workflows` rebuilt with `anthropic` 2.4.0 (0eff0689).
- [x] v2 design agreed and written (`STORY.md`, `SPEC.md`, 7fccc572).

---

## Phase 1: finish development of the pilot line (v1)

Goal: a line that passes its own rules, frozen, so that the nine-paper result belongs to one known version.

- [ ] **1.1 You: look at claim c006 again.** Why: on 2026-11-02 you marked it "right". We then assumed the checker had
      shown it wrong, but its last verdict calls the statement *supported*. The claim was disputed for another reason,
      which the review tool didn't show until 2026-11-03 (it now lists what failed). In your own terminal:
      `venv_workflows/bin/python scripts/literature_review/u3_review.py --run runs/literature-review/b1/dev-20261101-kp --redo kaempfer_plapp_2009-c006`
      Read "Why it is disputed", then answer for the claim as a whole. **⟂ depends:** if it is still right, answer
      `r` (nothing changes in the library); if a failed point makes it wrong, answer `w`.
- [ ] **1.2 Me: run the line on the dev run** so the stage that publishes your decisions applies your answer
      (marks c006 `rejected_by_user` if you answered `w`), and confirm it with the curator's counts.
- [ ] **1.3 You: rule on the nine flagged control statements** (G001 G004 G007 G014 G018 G031 G032 G036 G039) in
      `runs/literature-review/calib-dev6/untrusted/control_review_tuning.md`. For each, write "gold wrong" or "checker
      wrong" on its `Ruling:` line. Why: the checker must flag at most 4 of 40 correct statements. Until each flag is
      settled, we don't know whether the 9 flags are checker errors or errors in your gold statements. Use the new
      condition rule: a condition is kept if it is **necessary** or **scope** (`SPEC.md` §6.2). Note that G018 is a
      statement you ruled "checker wrong" on 30 October.
- [ ] **1.4 Me: read the rulings** (`control_rulings.py`, counts and IDs only) and act on them.
      **⟂ depends on 1.3:**
      - all or most "gold wrong": you fix those gold statements (I prepare the edits; you approve them), then 1.5;
      - several "checker wrong": I change the checker's prompt, aiming at over-flagging conditions, then 1.5;
      - a mix: both.
- [ ] **1.5 Me: rerun the tuning diagnostic** (`check_plants.py --set tuning`). **⟂ depends:** if controls flagged
      ≤ 4/40 and the plant rules hold, go on; if not, back to 1.4 (a prompt change), at most two more rounds, then we
      decide together whether the rule or the checker gives way.
- [ ] **1.6 You: decide the open backlog item** (a claim review decision that changes what the memo cites should
      send the memo back to be redrafted). Why: today your "right" on a disputed claim reaches the library but not
      the memo. v2 adopts the redraft (`SPEC.md` §4). **⟂ depends on your choice:** build it into v1 now (I add a
      runner rule, small), or accept the gap for the pilot and check the memo by hand at memo approval.
- [ ] **1.7 Me: close the second dev run** `dev-20261101-a2` (paused at the copy check for Li 2009). Why: it was a
      copy-check test; the real run repeats that gate. I mark it abandoned in the status, without deleting it.
- [ ] **1.8 Me: pin the Claude Code version.** The lock says 2.1.308; this machine runs 2.1.317. Why: the freeze
      records the version, and a version change means rerunning the calibration. I record 2.1.317 if the tuning rerun
      (1.5) ran on it. **⟂ depends:** if Claude Code updates again before the freeze, rerun 1.5.
- [ ] **1.9 Me: final development checks.** All tests pass, `check_permissions.py`, `lr_status.py --write`, commit.
- [ ] **1.10 You: approve the freeze** (`/lr-freeze`, v1). Why: from here on nothing that affects results can change
      during the run.
- [ ] **1.11 Me: live hook checks** right after the freeze. These can only be tested in the execute phase: an edit
      to a frozen file asks first (H3); a station run by hand is blocked (H4); a stop with a stale status is blocked
      once (H5); the foreman can't `cat` a file (H6).

## Phase 2: the pilot run (nine papers, v1). Target: about 8 November

Goal: the methods memo that model development needs, and the baseline v2 must beat.

- [ ] **2.1 Foreman: start the execute run** (`run_b1.py init`), then the frozen calibration (C0) on the held-out
      set, scored once.
- [ ] **2.2 You: accept or reject the calibration result.** **⟂ depends:**
      - passes: go on;
      - fails: the run stops, the line goes back to development (thaw), we fix, recalibrate on the tuning set, and
        you approve a new freeze. The held-out set is not redrawn.
- [ ] **2.3 You: copy decisions** for Li 2009 (no DOI printed) and Lervåg (arXiv only): read this copy, or supply
      the version of record.
- [ ] **2.4 Foreman: run the line** (prepare, extract, quote check, check, revise, publish, draft memo, blind
      reread). It pauses at your gates; I report usage.
- [ ] **2.5 You: claim review** (`u3_review.py` in your terminal). Why: these are the claims where independent checks
      disagree, plus every prose claim a model option rests on.
- [ ] **2.6 Foreman: citation check and skeptic review.** Then **you: approve the memo.** Then the backup.
- [ ] **2.7 Me: repeat test** on one paper (extract and check twice, compared claim by claim). Why: it measures how
      much answers vary between runs.
- [ ] **2.8 ⟂ depends on what 2.4–2.6 find:** a real error found by your review, the blind reread or the skeptic
      stops the run (thaw, fix, recalibrate, new freeze, rerun the affected claim family). Otherwise the pilot is
      done.
- [ ] **2.9 You and me: take the memo's model options into Pillar 2** (model development, Lab Report 3).

## Phase 3: build v2 (alongside model work). Proposed: 9 to 22 November

Goal: the new engine, its safety features and the Find and Get steps. Details per step: `SPEC.md` §13. Each step
ends with its test passing and a commit.

- [ ] **3.1 Build 1: single-call module** in `agentlib` (no tools, schema-checked answer, usage events in the
      ledger), with test T-input: can page images be given inside the call? **⟂ depends on T-input:**
      - passes: inputs go in the call;
      - fails: the model gets Read, limited to the call's folder, and the transcript is audited.

      Either way the stations above it don't change.
- [ ] **3.2 Build 2: usage cap** (pause at 90 % of the five-hour window, pause report, `--continuous`,
      `--parallel N`), with T-cap on a fake `claude`.
- [ ] **3.3 Build 3: paragraph numbers** in the paper-preparation step; the quote check works per paragraph;
      T-paragraph.
- [ ] **3.4 You: approve the foreman guard** (a new hook under `.claude/`). **Me: Build 4**, the guard and the
      session-role message, with T-guard.
- [ ] **3.5 Build 5: the new stage codes and names** in the stage graph, the runner and the status display.
- [ ] **3.6 Build 6: reading and checking as single calls**: extract facts (R1), assumption list (R3a), check facts
      with condition labels (R3), revise (R4), blind reread (R5); the new schemas (`SPEC.md` §5).
- [ ] **3.7 Calibrate the new checker** (C1, tuning set). **⟂ depends:**
      - meets the rules: go on;
      - doesn't: prompt work as in 1.4–1.5. A different engine is a different checker, so this can take several
        rounds.
- [ ] **3.8 Build 7: Get**: fetch (G1), copy check with retraction status (G2), the paragraph index (G4).
      **⟂ depends on T-embed:** you write about 20 query-to-paragraph pairs; the local embedding model with the best
      recall is chosen and pinned.
- [ ] **3.9 You: approve three new agent definitions**: `beth-searcher` (queries, citations to follow),
      `morty-screener` (relevance) and `beth-scout` (web scout). **Me: Build 8**, Find: search and citation chaining
      (F3), screening (F5), scout (F4, on `agentlib`), the coverage ledger, the shortlist report (F6); T-coverage.
- [ ] **3.10 Build 9: v2 on the nine pilot papers** (T-pilot), compared with the pilot. **⟂ depends:**
      - v2 matches or beats the pilot (calibration passes, the pilot's published claims covered or each gap
        explained): go on;
      - worse on something: fix that stage and rerun T-pilot.
- [ ] **3.11 Build 10:** `SPEC.md` replaces `WORKFLOW.md` (archived); update `TEMPLATE.md`, the contract format and
      the foreman's definition; tests, status.
- [ ] **3.12 You: approve the freeze** (v2). Then the live hook checks again, plus the foreman guard live.

## Phase 4: first discipline, numerical methods. Proposed: from about 23 November

Goal: the first discipline-scale review, tested small before it grows.

- [ ] **4.1 You: write the discipline brief (F1).** It holds the questions, 3–7 seeds (at least one review, at least
      one you know well), scope, inclusion rules, and **about 20 must-find papers** for the recall test. I draft the
      format and you fill it in. Why: everything after this depends on the questions and scope.
- [ ] **4.2 Foreman: queries, search, chaining, scout, screening** (F2–F5).
- [ ] **4.3 You: approve the shortlist** (F6). You see the shortlist, a random sample of 20 rejects, the queries and
      their hits, and the recall result. **⟂ depends on recall:**
      - target met: go on;
      - missed: add queries and rerun the search, or accept with a note.
- [ ] **4.4 Foreman: a 40-paper pilot**: fetch, copy check, prepare, index, then read and check the core papers among
      them.
- [ ] **4.5 You: fetch the paywalled papers** the script couldn't get (Zotero Connector), from its list. **⟂
      depends:** on how many there are.
- [ ] **4.6 You: hand check** 30 random verified facts against the PDFs. Why: the only check of the whole chain by a
      person, on a random sample. **⟂ depends:** errors found go to the repair rule (thaw, fix, recalibrate,
      refreeze) before scaling.
- [ ] **4.7 Me: same-effort single-agent comparison** (T-single) on the same 40 papers. Why: the course's standing
      rule; the multi-step design has to earn its extra parts. **⟂ depends:** if one agent does as well, we simplify.
- [ ] **4.8 Foreman: scale up** to the full shortlist; usage pauses as needed (`--continuous` if you choose).
- [ ] **4.9 You: fact review** (L2) of the disputed facts.
- [ ] **4.10 Memo** (M1–M3), then **you: approve it** (M4). Then the backup and restore test (K1).

## Phase 5: the other disciplines

Each repeats Phase 4 (brief, search, shortlist, 40-paper pilot, hand check, scale-up, review, memo, backup), in this
order:

- [ ] **5.1 Planetary observations and setting.**
- [ ] **5.2 Experiments**: sublimation; ice–sediment; ice sintering and metamorphism. **⟂ depends:** one review or
      three, decided at its brief.
- [ ] **5.3 Sublimation physics.**
- [ ] **5.4 Material properties: on demand.** No sweep. When the model needs a property, it goes through the
      new-question path (Phase 6). **⟂ depends:** the brief format is decided the first time it is needed.

Calibration (C1) is rerun before each discipline's frozen run if anything changed: prompt, model, engine or Claude
Code version.

## Phase 6: proving the library works (the finish line)

- [ ] **6.1 The new-question path, end to end.** You ask a question that no memo answers. The system searches layer
      1, extracts and checks the paragraphs found, publishes them, and writes a short memo. Why: this is the "add
      questions later" promise; it must work, not just exist.
- [ ] **6.2 A full restore** of the library from backup into an empty folder: 0 hash mismatches, and the memos'
      citation checks pass on the restored copy.
- [ ] **6.3 Coverage report:** every candidate in every discipline has exactly one recorded decision; every included
      paper is in layer 1; every core paper has its facts in layer 2.
- [ ] **6.4 You: sign off** that the library is fit for model work and the manuscript.

---

## Running notes

- **Usage:** all model work uses your five-hour window; the runner pauses at about 90 %. Phases 2 and 4–5 are the
  heavy ones.
- **Your time:** the main asks are 1.3, 2.2–2.6, 3.4, 3.9, 4.1, 4.3, 4.5–4.6, 4.9–4.10, and each discipline's
  repeat. Everything else runs without you.
- **Rules that stay in force:** no edits to frozen files during a run; agent definitions and `.claude/` changes need
  your approval; nothing from the course or research folders is committed.
