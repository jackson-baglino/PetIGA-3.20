# Pillar 1 redesign: round 4 (Marker, sign-in, isolation test, foreman guard)

2026-11-03. This answers your responses to `21_pillar1_redesign_followup.md`. Two things are already done and
committed (sections 1 and 3). Section 6 has the last two questions before the plan documents.

---

## 1. Marker: you were right, and it is gone

The paper parser has been **MinerU2.5-Pro** since 29 October. It runs in its own environment, `venv_docparse`, and
it read 10 of 10 equations on Kaempfer and Plapp p. 7, against Marker's 7. Marker stayed in `venv_workflows` for
two reasons only: as a fallback path in the paper-preparation step (P1), and for the practice units workflow.

Done (commit `0eff0689`):
- **Marker removed:**
  - the fallback path in `prepare.py`;
  - `marker_run.py`;
  - `marker-pdf`, `surya-ocr` and `torch` from the environment.
- **`venv_workflows` rebuilt** from what the scripts actually import:
  - `numpy`, `scipy`, `pillow` and `PyYAML`, at the same versions as before (they used to arrive as Marker's
    dependencies);
  - `markdown2` and `ruff`;
  - **`anthropic` 2.4.0**, the current release.

  That is 22 packages instead of 90, and 243 MB instead of 1.3 GB. The requirements file and lock are updated.
- **Tests:** all workflow tests pass in the new environment (literature review 134, units 48, flowchart 6).
- **Still to do:** the next real stage run is the first one in the new environment. The units workflow still calls
  Marker and will need it reinstalled or replaced when it is used again, as you said.

One leftover: Marker's model weights are still on disk, 3.2 GB in `~/Library/Caches/datalab`. Nothing uses them
now. Shall I delete them (question 2 in section 6)?

## 2. No Console key: the plan uses your login only

Without a Console key there is one back end, **B**: `claude -p` on your login, drawing on the five-hour window.
- The call module is built for B only.
- The A path stays a documented option that we don't build.
- `anthropic` 2.4.0 is installed but unused for now, which costs nothing.

**The Agent SDK without a key.** Your classmates skipping the API key step fits how the SDK works. It runs the same
`claude` program, and when no key is set, that program uses the login it already has. The docs don't promise this,
but it is consistent with everything we read.

**Your condition for the SDK:** if it automatically read CLAUDE.md and MEMORY.md, we'd drop it. The isolation test
(section 3) answers this for the `claude` program the SDK runs. Without isolation it reads both. With our flags plus
the new setting, or with `--safe-mode`, it reads neither. So the SDK *can* be isolated. Your condition doesn't rule
it out.

Your notes point both ways, so here are the two options. Question 1 in section 6 asks you to choose.
- **(a) Claude Code agents in `.md` files, all run through `agentlib.py`, the scout included.** Recommended.
  - It is already built, tested and audited: model pinned, tool calls checked against the allowed list, every
    attempt in the ledger.
  - The scout is the only tool-loop worker, so the SDK would add a second way of running agents for one agent.
- **(b) The SDK for the scout.** You would learn the current standard tool. It needs one sign-in test and its own
  checks, which repeat what `agentlib.py` already does.

**Google's ADK** (the repository your TAs shared) is Google's agent kit. It is built around Gemini, and other
models can be plugged in through adapters. It does the same job as the Agent SDK for a different vendor, so there is
nothing in it for this project.

## 3. The `--safe-mode` test, and what it found

You asked to test whether a worker can see memory. I ran a **canary test**:
- I planted unique marker words in all three places Claude Code keeps auto memory for this project:
  - the folder named in the project settings;
  - the folder named after the project path;
  - the folder named after the repository root (`PetIGA-3.20`).
- I asked headless workers, with no tools, which markers and instructions they could see.
- A worker with **no isolation** was the control, to prove the test can detect a leak.

| Worker | CLAUDE.md seen? | Memory marker seen? |
|---|---|---|
| no isolation (control) | **yes**: it named the Fable character "rick" and the Lab Report 3 date, and listed CLAUDE.md, CLAUDE.local.md, USER.md and both MEMORY.md files | **yes** |
| our flags (`--setting-sources ""` and the rest) | no | **yes: the repository-root memory index** (14 kB) |
| our flags + `CLAUDE_CODE_DISABLE_AUTO_MEMORY=1` | no | no |
| `--safe-mode` | no | no |
| `--safe-mode` + a `--settings` file denying `Read(./CLAUDE.md)` | the Read was **denied**; the allowed README read worked | — |

What it means:
- **Our flags kept CLAUDE.md out but let auto memory in**, as the docs said.
- **Both fixes work.**
- **`--safe-mode` keeps your login and still enforces permission files.** So it is a valid replacement for
  `--bare`, which can't use your login.

**Fixed** (commit `61446b41`): `agentlib.py` and `tres.py` now start every worker with auto memory off, and a test
checks it. The backlog item is closed, and the result is recorded in WORKFLOW §12.

**Earlier runs were probably not affected.** Workers start in a temporary folder outside the repository, and the
repository-root memory loads only for a worker started inside it, as in this test. That is inferred: I did not
retest the old runs. The calibration reruns anyway after the redesign.

**Your memory question.** We don't need to restructure memory. The fix belongs to the worker: auto memory is turned
off for workers. The main session keeps its memory as it is, and so does the foreman, which is a normal Claude Code
session.

A side note: the repository-root memory folder is an older memory set from sessions started at `PetIGA-3.20`
itself. It is not this project's memory and was left untouched.

## 4. A guard so the line can't be run from the main session

Agreed. The design, for the plan:
1. `foreman.py`, which starts the foreman session, sets a marker in that session's environment (`LR_ROLE=foreman`).
   Hooks see the session's environment.
2. **A new hook** (PreToolUse on Bash, in the project settings). In the **execute phase** it blocks
   `run_b1.py run`, `init` and `approve` when that marker is missing. The message names the command to start the
   foreman. Status and verify stay allowed everywhere, because they only read.
3. **In the develop phase,** the main session may still run *dev* runs (manifest `"dev": true`), as now. Those are
   build tests, not the real review. The guard blocks a non-dev run.
4. **The session-start message** names which session you are in ("main session: design and build" or "foreman:
   runs"), so the difference is visible from the first line.

The hook lives in `.claude/`, so I'll ask before adding it. It goes into the plan as a build step with its own test:
the main session tries `run`, gets blocked, and the foreman's attempt goes through.

## 5. Where the design now stands

| Topic | Decision |
|---|---|
| Parser | MinerU only; Marker removed |
| Engine for decision steps | single calls through `claude -p` on your login (back end B): no tools, one turn, schema-checked answer |
| Tool-loop agents | the scout only, on `agentlib.py` (or the SDK, question 1) |
| Worker isolation | our flags + auto memory off (done); `--safe-mode` tested as the alternative |
| Interactive agents | `.md` agents: the foreman for runs (with the guard), helpers such as the writer and the critic |
| Control | the stage runner; usage cap at about 90 % of the five-hour window, then pause or wait for the reset |
| Agent justification | each model step names its reason (information, permissions, judgment, breadth); no loop unless the next input depends on its findings |
| Stage codes, disciplines, condition rule, memo, two-layer extraction | as agreed in notes 20 and 21 |

## 6. Last questions before the plan documents

1. **Scout:** (a) `agentlib.py` (recommended) or (b) the Agent SDK?
2. **Marker's weights** (3.2 GB, `~/Library/Caches/datalab`): delete them?

When you answer, I'll write the two documents: the story version for you, and the agent specification that
replaces WORKFLOW.md. The paused items come first in the plan's order: the misjudged claim c006, the nine control
rulings, then the freeze.
