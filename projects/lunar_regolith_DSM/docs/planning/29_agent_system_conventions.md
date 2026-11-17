# 29. One way to build, start, contain and test every agent

Date: 2026-11-17. Status: the user agreed the direction. Sections 2–5 are proposals; nothing in them is built
yet. Today's build is only the mathematician (doc 28). The move of the older agents and the Python index are for
2026-11-18.

## 1. Words to use

- **Claude Code is the agent harness.** It runs the model loop and applies the tools, permissions, hooks and
  sandbox.
- **Our launcher** is a Python script (`scripts/agentkit/`) that starts Claude Code **headless** (`claude -p`,
  print mode). It passes the agent's definition, settings and answer schema, then checks the result.
- **The harness, in the course's sense,** is the launcher plus the hooks, schemas and graders: everything outside
  the model that decides what counts.

In one sentence: "Each agent runs headless through one launcher, which applies its definition, permissions and
sandbox, and grades its answer."

## 2. Sandboxes at four levels

A sandbox is set in a **settings file or an OS profile**, never in CLAUDE.md, because CLAUDE.md is prompt. One
fact shapes the whole design: our workers start with `--setting-sources ""`, so **the user and project settings
never reach a worker.** Only what the launcher passes reaches it. Each level therefore has its own carrier:

| Level | Protects | Carrier | Contents (proposal) |
|---|---|---|---|
| 0. You, every project | your interactive sessions | `~/.claude/settings.json` | read denies for secrets (`~/.ssh`, `~/.aws`, keychains). A full sandbox only after testing that git, ssh to the HPC and the run scripts still work. |
| 1. This project | my interactive sessions here | `.claude/settings.json` | today's denies (library, untrusted stores), plus the Claude Code sandbox for Bash once tested |
| 2. A line run | the whole launcher process: the claude process, its hooks, the Read and Write tools, every child | OS profile, `sandbox-exec -f <profile>.sb` (the literature-review design, `wf_launch.py`) | writes only to the run folder and temp; reads of secrets denied; network only HTTPS for the model call |
| 3. One agent | the agent's Bash commands and every process they start | the agent's `settings.json`, passed with `--settings` | the Claude Code sandbox (`failIfUnavailable: true`, `allowUnsandboxedCommands: false`, an empty network allowlist), permission rules, guard and Stop hooks |
| 4. The solver | an agent-built binary | `run_test.sh` wraps it in the offline profile | no network; writes only to the test's run folder |

**Rule: inner levels may only narrow.** The launcher builds each agent's settings from a project base plus the
agent's own file. A test fails if an agent's file widens the base (an extra write path, a network domain, a
removed deny).

**Should the literature-review sandbox become the project sandbox?** Its profiles should become the shared
level-2 profiles, moved to `scripts/agentkit/profiles/`:
- `offline`: no network; writes to the run folder and temp;
- `offline-plus`: the same, plus one named place;
- `agent`: today's `anthropic.sb`, for the process that calls the model.

The Zotero profiles stay with their scripts. They need two fixes first [checked, 2026-11-17]:
- **Reads are open.** Every profile starts with `(allow default)`, so `~/.ssh` and `docs/planning` are readable.
- **Network is too broad.** `anthropic.sb` allows HTTPS to any host. Level 3 closes this for commands the agent
  runs (no network at all), but the profile should also deny reads of secrets.

## 3. Asked every time an agent is built

The `new-agent` skill becomes a checklist that asks for each item. A scaffold script writes the folder with the
**strictest defaults**, so any widening is a decision with a written reason:

1. purpose, and the known-answer task it must pass first;
2. model;
3. tools: the minimum, default none;
4. skills: default none;
5. reads and writes: default the run folder only;
6. sandbox level and profile: default `offline`, no network;
7. environment variables: default PATH, HOME, TMPDIR;
8. hooks: guard and Stop, always;
9. work-order and answer schemas;
10. grader;
11. an enforcer for every "Never";
12. tests.

The reasons for any widening go in `tests/enforcers.json`, which the tests read.

## 4. One way to start any agent

- `scripts/agentkit/launch.py` provides one function: `run_agent(name, task, run_dir)`. In order, it:
  1. finds the agent by name (it fails on a duplicate);
  2. validates the work order;
  3. composes the settings (base plus agent) and builds the environment from `env.allow`;
  4. wraps the call in the level-2 profile and calls `claude -p`;
  5. audits the transcript;
  6. validates the answer;
  7. runs the grader;
  8. writes the run record.
- **A test fails if any other file calls `claude -p`.**
- **Today, eight scripts call it:** `agentlib.py`, `foreman.py`, `review_session.py`, `run_b1.py`,
  `t_security.py`, `tres.py`, `run_md.py` and `uvlib.py`.
- **Order of migration:**
  1. the model line, which moves to agentkit as we build it;
  2. design critique, which is small;
  3. units verification, which is finished, so archive it rather than migrate it;
  4. the literature review, which stays frozen until we come back to it.

  The two interactive launchers (`foreman.py`, `review_session.py`) start your own sessions, not workers. They
  stay apart from agentkit, under the same tests.

## 5. Testing the convention

Three kinds of test, each with a clear job:

1. **Structure** (fast, runs before every real call):
   - every agent folder has its files;
   - the name prefix matches the model, and the model is a full ID;
   - the frontmatter tools match the settings;
   - no name appears twice;
   - no agent widens the base;
   - every "Never" has an `enforcers.json` entry;
   - only agentkit calls `claude -p`.
2. **Barrier probes.** No model is needed. For each enforcer, a probe tries what it forbids:
   - write outside the run folder;
   - read `docs/planning` or `~/.ssh`;
   - open a socket;
   - `git push --dry-run`;
   - a guard given broken input.

   Each probe must be blocked, logged and stop the call.
3. **Tests of the tests,** which show the tests cover the full breadth:
   - every rule has an ID and at least one test;
   - every test has a **planted violation**: a copy of an agent folder with one thing broken, such as a removed
     deny, an extra write path or an alias model;
   - the right test must fail on it;
   - a meta-test fails if any rule lacks a test or a planted violation.

   An independent reviewer (beth-critic, in a clean folder) may then flag rules the list misses. It may flag;
   it never passes anything.

## 6. The existing agents

- **beth-writer:** writes and reviews documents, checking every claim, number and citation against its source.
  In review mode it reports findings and never edits. I start it for reports, and the model line's write-up
  stage (M9) uses it. **General: keep it in `.claude/agents/`.**
- **beth-critic:** reads one workflow design and reports problems with concrete failure scenarios, checking
  factual assumptions on the web. `scripts/design_critique/run_critique.py` runs it headless in an empty folder.
  **General: keep it.**
- **beth-unit-checker and morty-unit-reader** belong to the units-verification workflow of 2026-10-22, which
  checked the 54 rows of `memory/semantic/units.md` against page images (run 20261022-155335).
  - They are in `.claude/agents/` because that workflow came before the `workflows/<name>/agents/` convention
    (`scripts/agentlib.py:5`).
  - The two-agent design was for independence: the reader saw only the page, never the claimed value, and the
    checker settled disagreements.
  - The job is done, so neither is a general agent. **Proposal:** move both to
    `workflows/units-verification/agents/` and mark that workflow archived.
  - A future "check this number against the page" job belongs to the paper reader we will build.
- **Five permission files with no agent:** `beth-memo-checker`, `beth-planter`, `jerry-screener`, `morty-desk`
  and `morty-scout` (in `.claude/permissions/literature-review/agents/`). They were designed in the
  literature-review `WORKFLOW.md` and never built. **Proposal:** delete them in the move.

## 7. The Python index (2026-11-18)

There are 196 tracked Python files [checked]:

| Folder | Files |
|---|---|
| `scripts/literature_review` | 61 |
| `postprocess` | 33 |
| `preprocess` | 30 |
| `scripts/units_verification` | 17 |
| `scripts/model_development` | 15 |
| `.claude/hooks` | 10 |
| other | 30 |

190 of them open with a docstring.

**Proposal:**
- A script writes `docs/python_index.md` from the docstrings, grouped by workflow, then by role (launcher, station,
  hook, check, test), then alphabetically. Because it is generated, it cannot drift from the code, and a test
  fails on a file with no docstring.
- I then add a keep / archive / delete column, each with evidence: what calls the file, its last commit, and
  whether a live line uses it.
- You decide. Files to archive go under an `archive/*` tag, not straight to deletion.

## 8. CLAUDE.md

- **Worth doing, and done in a small way.** Project CLAUDE.md now says that building an agent goes through the
  `new-agent` skill, that every rule needs an enforcer outside the prompt, and that one agent must pass a
  known-answer task before the next is added.
- **It is still a prompt.** It reminds me; the scaffold and the tests are what make the rule hold.
- **Global `~/.claude/CLAUDE.md`:** two general lines would carry the same habits to the enceladus and cryosuction
  projects: "every agent rule needs an enforcer outside the prompt" and "one agent on a known answer first".
  Needs your OK.

## 9. Order

- **Today:** mathematician steps 0–2 (doc 28), then Lab Report 3 at about 6 pm.
- **2026-11-18:**
  - the move of the older agents and the five orphan files;
  - agentkit as the one launcher for the model line;
  - the level-2 profiles with read denies;
  - the convention tests;
  - the Python index.
- **After that:** mathematician steps 3–5, then the next agent.
