# Pillar 1 redesign: answers to your responses on note 19

2026-11-03. This answers your point-by-point responses to `19_pillar1_redesign_comparison.md`. The sections follow
yours. Section 5 lists what is agreed and section 6 the few questions left before I write the two plan documents.

Where a fact was checked, its source is given: a file and line, a run folder, or a documentation page. Anything I
inferred is labeled as inferred.

---

## 1. Phase 1 and Phase 2

### 1.1 Is fetching a script? What needs a decision?

Yes. Fetching is deterministic, and in our design it is a script. Given a DOI, the script tries these sources in
order:
1. your Zotero library (`zotero_pull.py`, built);
2. the open-access link from Unpaywall or OpenAlex;
3. arXiv, only when the paper exists nowhere else.

If none of these works, the script adds the paper to a list for you to fetch with the Zotero Connector. A script
then checks the copy it got: the DOI's Crossref metadata against the title and authors on page 1, the file hash, and
retraction status. This check is built today as "check copy" (A2).

Searching is also a script. It queries OpenAlex and Semantic Scholar and follows references and citing papers.

Three things in Phase 1 need judgment. None of them needs an agent with tools:

| Decision | Why a script can't make it | What makes it |
|---|---|---|
| **What to search for** (queries per discipline) | choosing words that find the right literature is judgment | one model call per discipline: brief in, list of queries out |
| **Is this paper relevant?** (title and abstract) | relevance to our questions is judgment | one model call per paper: abstract and discipline brief in, `{decision, reason, discipline}` out |
| **Which citations deserve another hop** | a reference list mixes core papers with background | one model call per paper: its reference list in, the references worth following out |

Two kinds of discovery fall outside the databases: NASA technical reports, theses, and conference abstracts. Only
that open-ended search, the "scout," is a real agent: it decides what to look at next based on what it found. Even
then, a script checks and fetches everything the scout proposes.

### 1.2 One team per discipline

Agreed. In Lecture 10's terms this is a split for **breadth**, and it is legitimate. I suggest building it like this:
- **A discipline brief** for each discipline: a JSON file with its questions, seed papers, scope boundaries and
  inclusion rules. You approve each one.
- The script runs the three decision calls above for every discipline **in parallel**, each call using its own
  discipline brief.
- **Duplicates are handled by script.** One paper found by two disciplines is read once and tagged with both.
  Workers never choose papers, so two of them cannot end up reviewing the same one (see 3.5).

The parallelism comes from the script running many calls at once. Having one brief per discipline is what makes
each call judge with that discipline's criteria.

### 1.3 Extraction: more information, not less

Your concern was blind spots. The two-layer design answers it like this:
- **Layer 1, the full index.** Every paragraph, equation, table and figure caption of every paper we collect is
  stored and searchable by keyword and by meaning. No model is involved. Nothing is thrown away, so a question we
  think of later can still find text we never extracted.
- **Layer 2, verified facts.** For core papers, the reader extracts **every** fact in the paper, not only the facts
  that answer our questions. Every fact gets its conditions, its quote and its location, and every fact goes
  through the checker.

A new question then needs two steps:
1. search Layer 1;
2. extract and check the paragraphs it finds that Layer 2 doesn't already cover.

This is the "add questions later" path you described.

Which papers count as "core" is a decision for you: the discipline's seed papers and the papers its relevance
screen ranks highest. Every other paper is in Layer 1 only until a question needs it.

### 1.4 Paragraph numbers

Paragraph numbers will be added. The parser already writes `layout.json`, with the blocks on each page and their
type and position (`scripts/literature_review/prepare.py:8`). A script will number the text blocks in reading
order and give each one an ID such as `p7.¶3` (page 7, paragraph 3).

Each extracted fact must then cite a paragraph ID. The quote check (E3) also gets stricter: the quote must be found
**inside that paragraph**, not just somewhere on the page. In your claim review you will see the paper, the page,
the paragraph number and the quote.

### 1.5 Copying what Scenario C does

All five, as you said:
1. a coverage ledger;
2. a pilot before scaling;
3. caps set in the run command;
4. a same-budget single-agent baseline;
5. typed JSON briefs.

Section 1.6 changes the third one, as you asked.

### 1.6 A cap based on the usage limit, and a "keep running" mode

**This is possible with what we already run, and I checked it on your account.** Every `claude -p` call made with
`--output-format stream-json`, which is how our workers run, reports a `rate_limit_event`. It gives:
- `rateLimitType` (`five_hour`, `seven_day`);
- `utilization`, a fraction from 0 to 1;
- `resetsAt`, the reset time;
- `status` (`allowed`, `allowed_warning`, `rejected`).

The trial run's transcripts already contain four of these events, all for the five-hour window, with utilization
0.33 to 0.46 (`runs/literature-review/b1/dev-20261101-kp/transcripts/`). The documentation types this event only in
the Agent SDK, not in the `claude -p` pages, so we treat it as working but undocumented. I did not see a seven-day
(weekly) event; it may appear only near that limit.

**The rule.** The runner reads utilization after every call:
1. **At 90 % or more of the five-hour window**, it starts no new item and lets the items already running finish.
   Each item (one stage for one paper) is already all-or-nothing and resumable, so an item boundary is the natural
   place to pause.
2. It then writes a **pause report**:
   - what ran since the last pause;
   - where every paper stands;
   - what is left;
   - the reset time.

   It prints the same report and stops. You decide whether to continue.
3. **"Keep running" mode** (`run --continuous`). Instead of stopping, the runner waits until the reset time plus a
   few minutes, then resumes by itself. It writes the pause report each time.
4. **The weekly limit always stops the line and waits for you**, in both modes. Its reset can be days away, and the
   last 10 % is yours for other work.
5. **Fallbacks.** If a call reports no utilization, the runner applies a per-run call cap that you set. If a call
   fails with "You've hit your session limit," the runner records the item as not done and pauses; it does not
   count the item as failed. The documentation gives no specific exit code for this failure, so the runner matches
   the message.

The 90 % threshold is a setting, recorded in the run manifest. With several items running in parallel, the
utilization can end a few percent past it, because running items are allowed to finish.

---

## 2. Cutting agents down to what they need: the build

### 2.1 How a decision step runs now, and how it would run

Today the station script runs `claude -p` for each judgment, whether it reads claims, checks them, drafts the memo
or reviews it (`scripts/agentlib.py:159`). That starts a full Claude Code session, which brings:
- Claude Code's own system prompt and tool definitions;
- a tool loop in which the model decides to open page images with the Read tool, one call at a time;
- a sandbox profile, a permission file and hooks, all needed **because** the model has tools;
- a transcript, which we then parse.

The proposed way is **one request to the Messages API**. The script itself builds the whole input:
- the instructions;
- the page images, or the numbered paragraphs, it has chosen;
- a JSON schema the answer must follow.

One response comes back, and the API guarantees it matches the schema. The parameter is `output_config.format`
with `type: "json_schema"` (from the bundled Claude API reference, checked today). There is no tool loop, no file
access, no transcript to parse, and the model cannot choose its own inputs.

This is what **"move the single-call decision steps to the raw Messages API"** means. "Raw" means we call the API
directly with the `anthropic` Python package, without Claude Code in between.

### 2.2 What gets built: one call module

A single file, `scripts/lr/llm.py`, is the only place a model is called for decision steps. In outline (pseudocode;
the parameter names were checked against the API reference):

```python
def decide(step, model, effort, instructions, inputs, schema, run):
    request = {"model": model, "max_tokens": ..., "system": instructions,       # cached prefix
               "messages": [{"role": "user", "content": inputs}],             # page images + numbered text
               "output_config": {"format": {"type": "json_schema", "schema": schema}, "effort": effort}}
    key = sha256(canonical_json(request))         # same input -> same key; a rerun reuses a stored answer
    response = client.messages.create(**request)
    # refusal / max_tokens -> recorded as a failure, never parsed as an answer
    answer = json.loads(first_text_block(response))
    validate(answer, schema)                      # a second, local check
    ledger.append(step, model, effort, key, response.id, response.usage, stop_reason, seconds, sdk_version)
    store(run / "untrusted/calls" / f"{key}.json", request, response)
    return answer
```

Every station that makes a judgment calls `decide()` with its own schema. The ledger row replaces today's transcript
as the record of what was sent and what came back.

### 2.3 Each step, converted

| Step (old code) | Input the script builds | Output (schema) | Then checked by |
|---|---|---|---|
| Write search queries (new) | discipline brief | `queries[]` | you see them in the brief approval |
| Relevance screen (new) | title, abstract, discipline brief | `{decision, discipline, reason}` | coverage ledger; you see a random sample of rejects |
| Choose citations to follow (new) | one paper's reference list | `follow[]` with reasons | stopping rule (§19, section 1) |
| Extract facts (E1) | page images + numbered paragraphs of a group of pages | `facts[]`: statement, type, quote, paragraph ID, `conditions[]` (each with quote and paragraph ID) | quote check, now per paragraph |
| Assumption checklist (E4, part 0) | the whole paper, no facts | `assumptions[]` with paragraph IDs | used as the completeness reference |
| Check facts (E4, part 2) | facts of up to 4 pages, those pages, the checklist | per fact: verdict; per condition: necessary/scope/unnecessary (§4.4); missing conditions with quotes | quote check on the missing-condition quotes |
| Blind reread (E4 part 1, X1) | page images + only "equation 13, p. 4" | transcription | compared with the extracted fact by script |
| Revise (E5) | the failed fact, the checker's report, the pages | revised fact | goes through the checks again (at most 2 rounds) |
| Draft memo (Y1) | verified facts as JSON | memo with fact IDs | citation check (Y2) |
| Skeptic review (Y3) | memo + facts | objections[] | you, at memo approval |

The checks in the last column are all scripts already, or small extensions of them.

### 2.4 What this removes

Four things go away:
- **The worker sandbox and permission files**, for every step above. A request with no tools cannot fetch a URL,
  run a command or write a file, so a paper with planted instructions can at most corrupt the answer. The quote
  check and the checker catch a corrupted answer. The planted-instruction test (T-trifecta) shrinks to "the answer
  still passes the checks."
- **The Claude Code system prompt and tool definitions**, which today ride along in every call. That is inferred:
  the synthesizer calls on the trial run used about 23,000 input tokens per call (`runs/literature-review/b1/
  dev-20261101-kp/calls.jsonl`). I have not measured how much of that is Claude Code's own text.
- **Transcript parsing**, and the version churn of the `claude` command line.
- **A source of variation**: the model choosing which pages to open. The script decides the inputs.

### 2.5 What it costs, and one thing that may block it

- **Recalibration.** The checker was calibrated as a Claude Code agent at Claude Code's default effort. A different
  engine is a different checker, so the calibration (C0) and the build tests run again. Through the API, Opus 5.5
  defaults to `medium` effort, so the effort must be set explicitly and recorded.
- **Sign-in: this may block the raw API.** A research agent read the Anthropic docs today. What they say:
  - The Messages API takes an API key from the Claude Console, workload identity federation, or App Attest. Your
    enterprise Claude login is not among them (platform docs, "Authentication").
  - API calls are billed per token to a Console organization. They do not draw on your seat's usage allowance.
    Enterprise plans are not eligible for API credits ("API credits for subscribers").
  - `ant auth login` signs in to the Console, not to claude.ai.

  Checked in this shell: no `ANTHROPIC_API_KEY` is set, `ant` is not installed, and `venv_workflows` has `anthropic`
  0.46.0, far behind the current 2.4.0. So the raw API needs **a Console API key from an organization that pays for
  it.** You'd have to ask whether Caltech or the group provides one. One more catch: once `ANTHROPIC_API_KEY` is set
  in a shell, `claude` itself uses the key instead of your login, so the key must be given only to the call module.

### 2.6 Two ways to run the same single call

Because of the sign-in question, the call module gets **two interchangeable back ends**. Every station calls the
same `decide()`:

| | **A. Messages API** | **B. Claude Code, as a single call** |
|---|---|---|
| How | `client.messages.create(...)` with `output_config.format` | `claude -p` with **no tools** (`--tools ""`), one turn, `--json-schema` for the answer, and the inputs given in the call instead of read by the model |
| Sign-in | Console API key | your enterprise login (as now) |
| Draws on | the Console organization's spend limits | your seat's five-hour and weekly allowance; the 90 % rule of section 1.6 applies |
| Extra input per call | none | Claude Code's system prompt |
| Isolation | complete: nothing of ours is loaded | needs the right flags (section 3.6) |

Both remove the tool loop, the model's choice of inputs and the need for a worker sandbox.
- **B works today** without anything new from Caltech.
- **A is the cleaner engine** if you get a key.

The manifest records which back end ran, so a switch between them shows up as a different engine and triggers a
recalibration.

One thing for B must be tested before we rely on it: sending page images inside the call. Claude Code's
`--input-format stream-json` accepts messages as input; I have not verified that it accepts image blocks in them.
If it doesn't, B gives the model the Read tool limited to that call's page images.

---

## 3. The agent frameworks

### 3.1 What our Python stage runner is

The stage runner is `scripts/literature_review/run_b1.py` (566 lines). It is a loop over a fixed graph. The graph is
data, not code: each stage in `workflows/literature-review/flowchart.json` has a `run` block, for example the stage
that publishes your decisions (L1b):

```json
{"per": "paper", "requires": ["U3a"], "station": "curator.py", "gate": "none",
 "steps": [["curator.py", "publish", "--run", "{run}", "--paper-id", "{paper}", "--after-u3"]],
 "done": [{"path": "summary_publish_u3_{paper}.json", "absent": ["failed"]}]}
```

That block says four things:
- the stage runs once per paper;
- it may start only when your claim review (U3a) is done;
- it runs exactly this command;
- it counts as done when that summary file exists and doesn't record a failure.

`run_b1.py run` repeats these steps (`run_b1.py:400–471`):
1. **Find the next ready item**: a stage and paper whose required stages are done and whose own done files are
   missing.
2. **Your gate?** A stage done by you (approve brief, claim review, approve memo) writes a "waiting" file under
   `gates/`, prints the `run_b1.py approve` command, and is held. The line goes on with other papers.
3. **Otherwise run the step's command** through `wf_launch.py`. That runs the station script under a macOS sandbox
   profile, and only if the script is on an allow-list that only you edit. The output goes to a log file.
4. **Check the done files.** The runner trusts the files on disk, not the script's own report.
5. **Record the result** in the run's ledger, `runs.db`: stage, paper, start, seconds, ok or failed, log path.
6. **Stop on a failure**, naming the log. Otherwise go back to step 1.

Other commands:
- `init` makes a run folder and its manifest;
- `status` prints the stage table;
- `approve` records your decision on a gate;
- `freeze` and `thaw` lock and unlock the definitions;
- `verify` checks that nothing frozen has changed.

In the execute phase, `run` refuses to start if any frozen file changed.

This gives the same things LangGraph offers: a fixed graph, resumable state (the files and ledger), human gates
and a full record. It does so in about 600 lines we can read.

### 3.2 Raw Messages API

See section 2.1. In one sentence: the stage runner stays as it is, but a station that needs a judgment makes one
API request with a schema instead of starting a Claude Code agent.

### 3.3 How the Claude Agent SDK works, and where we'd use it

It is **not** a script that calls a model again and again. The loop belongs to the model:
1. you give it a goal and a set of tools;
2. it decides which tool to call;
3. the harness runs the tool and returns the result;
4. the model decides the next call, and so on until it says it is done.

That is exactly what `claude -p` does. The Agent SDK is that same harness as a Python library.

Under the hood, the Python package starts the bundled `claude` command line as a subprocess and talks to it in JSON
messages over standard input and output (Agent SDK docs, "Hosting," "The subprocess model"). So an SDK agent is a
`claude -p` session that your Python program drives and listens to, message by message.

Sign-in: the SDK quickstart says to set an `ANTHROPIC_API_KEY` from the Console. The docs never say that your
subscription login works for SDK scripts, and they don't forbid it for your own use either. (The rule they state
is about developers *offering* claude.ai login in their products.) It probably works, because the SDK runs the same
`claude` program, but that is inferred and must be tested.

What the SDK adds over our `claude -p` calls:
- **Our own Python functions as tools.** For example `openalex_search`, `crossref_lookup` and `already_have(doi)`.
  The agent then never needs Bash or general web access.
- **Hooks as Python functions** in the same program. An example is a check that refuses a tool call outside the
  discipline brief.
- **Typed message objects** instead of a stream of JSON lines to parse.
- **Turn and budget limits** set in code.

**Where it fits for us: only the scout.** The scout searches the web and the databases for one discipline, decides
what to look at next from what it found, and returns proposed papers with reasons. Every step after that is the
script and single calls. A sketch (structure only; the exact names come from the SDK docs at build time):

```python
tools = [openalex_search, semantic_scholar_search, crossref_lookup, already_have]   # our functions
options = {"model": ..., "tools": tools + ["WebSearch"], "max_turns": 40, "output_schema": PROPOSALS}
for message in agent_sdk.query(prompt=discipline_brief, options=options):
    ledger.log(message)          # every tool call and result is recorded
proposals = last_result(message) # -> the script dedupes, checks identity, fetches
```

### 3.4 LangGraph

Agreed: not now. Learn the API and Agent SDK first. If you want LangGraph later, a small side project will teach it.

### 3.5 Why a script instead of a foreman agent

What the foreman does today (`workflows/literature-review/agents/beth-foreman.md`):
- it runs `run_b1.py status` and `run_b1.py run`;
- it reports what happened;
- it is forbidden to edit, approve, run a station or read paper text.

It assigns no work. `run_b1.py` already decides what runs next. The foreman is an Opus session that reads the
script's output and paraphrases it.

You named two jobs where you'd want an agent's independence. Both are better done by code:

| Job | How code does it | Why that is stronger than an agent |
|---|---|---|
| **No two workers review the same paper** | Every paper has one row per stage in the ledger. The runner hands each worker one paper ID; a worker cannot pick papers. Duplicates across disciplines are merged by DOI, OpenAlex ID and normalized title before reading. | It is impossible, not just unlikely. An agent tracking assignments in its context can lose track. That is a known team failure mode (MAST: agents redoing each other's work, losing the conversation history). |
| **Keep workers on task** | A worker's only output is JSON that must match its schema; the scripts then check every quote and fact. Off-task output fails validation and is retried at most twice, then goes to you. | An off-task answer cannot enter the library, whoever notices it. An agent judging another agent's focus is a second opinion from the same kind of model. Our rule is to prefer a deterministic check to a second opinion. |

The judgment the foreman seemed to offer, "what should happen now that something odd occurred," still exists. It
belongs to this interactive session, with you present. The status script tells me what stopped and why, and we
decide together. When the line runs unattended, the runner follows fixed rules: retry, skip the paper, or pause
(section 1.6).

### 3.6 The sentence that was unclear

It was two separate caveats.
1. **"The Agent Builder wind-down dates."** In the frameworks table I said OpenAI's no-code Agent Builder is
   reportedly being shut down. OpenAI's own page refused our request, so that claim rests only on secondary sources.
   It doesn't affect our build; it only qualifies one cell of the comparison table.
2. **`--setting-sources ""` versus `--bare`.** Every headless call we make passes `--setting-sources ""`
   (`scripts/agentlib.py:160`). It is meant to keep the worker from loading anything of ours: CLAUDE.md, memory,
   hooks, settings. Claude Code's documentation now recommends a flag called `--bare` for scripted calls. The open
   question was whether our flag isolates as completely as `--bare` does. If it doesn't, a worker could be reading
   instructions we never meant it to see.

   What the docs for Claude Code 2.1.317 say:
   - **`--bare` cannot be used with your login.** It skips hooks, skills, plugins, MCP, auto memory and CLAUDE.md,
     but it also never reads your login: it accepts only an API key.
   - **`--setting-sources ""` does less than we assumed.** It drops the settings files and, inferred from the docs,
     the CLAUDE.md files tied to them. It does **not** stop *auto memory*, the memory files Claude Code keeps for
     this project. Our workers don't set `CLAUDE_CODE_DISABLE_AUTO_MEMORY=1`, so a worker may be loading the
     session's memory index. This is unverified on our runs. I parked it in the backlog as something to check before
     the freeze (commit 969bf76f).
   - **`--safe-mode`** is a newer flag that turns off CLAUDE.md, auto memory, hooks, skills and MCP but keeps your
     login. It is the likely right flag for back end B and for the scout. It must be tested with a canary: plant a
     marker in memory and confirm that the worker can't see it.

---

## 4. Decisions

### 4.1 Why nine papers took so long

The time went into **building and testing the line, not reading papers.** Machine time per stage on the one-paper
trial run (`runs/literature-review/b1/dev-20261101-kp/runs.db`, ledger table `line`):

| Stage | Minutes |
|---|---|
| prepare paper (P1: parse text, equations, page images) | 4.4 |
| extract claims (E1) | 1.4 |
| check claims (E4) | 3.2 |
| revise claims (E5) | 2.6 |
| check copy, quote check, publish, backup | under 0.1 each |
| draft memo (Y1), citation check (Y2), skeptic review (Y3), once per run | 1.3, 0.8, 1.8 |

That is about 12 minutes of machine time per paper. Run one after another, nine papers take about 2 hours. Run four
at a time, they take about 30 minutes. The days went into:
- designing the stations;
- the calibration (plants and controls);
- the security and determinism tests;
- fixing what the tests found.

That is a one-time cost, which only pays off when it is spread over many papers.

You're right that for nine papers the line saves no time over reading them yourself. What it gives at nine papers
is checks you would not do by hand: every quote is found on its page by script, every equation and value is reread
blind, and every condition is checked against a whole-paper list of assumptions. Its time advantage only appears at
hundreds of papers.

On scope: agreed, the review is one part of the project. I suggest:
- **Finish the nine-paper methods review on the current engine** as the pilot, with only the two paused items
  before it. It answers the methods question that model development needs.
- **Do not let the large-scale build hold up the model work** (Pillar 2). The build can go on alongside it.
- **Use the pilot as the baseline.** The new engine reruns the same nine papers and the same calibration, and has to
  match or beat the pilot. That is the same-budget comparison we owe anyway.

### 4.2 Disciplines

The union of your list and ours, in your order:

| Order | Discipline | From |
|---|---|---|
| 1 | Numerical methods (phase field, vapor transport, interface methods) | both |
| 2 | Planetary observations and setting (large-scale planetary physics, lunar polar cold traps, observations) | both |
| 3 | Experiments: sublimation; ice–sediment (regolith); ice sintering and metamorphism | yours + our "ice sintering" |
| 4 | Sublimation physics (theory, kinetics) | ours |
| 5 | Material properties (conductivities, diffusivities, densities) | ours |

I placed 4 and 5 last because you didn't order them. Material properties may work best on demand: a property is
looked up when the model needs it, rather than as a sweep of its own.

### 4.3 Extraction

Agreed: two layers (section 1.3).

### 4.4 The condition rule, written so a checker can apply it

Your definitions:
- **necessary:** without the condition, the fact is not true;
- **complete:** every necessary condition is present.

One sharpening is needed. Most papers never say whether a fact would fail without a condition. They only state the
fact under it. A reader cannot know that the fact is false outside the tested range, only that the paper didn't
establish it there. Dropping such a condition is exactly the error our "dropped condition" plants test for. So each
condition gets one of three labels:

| Label | Test | What happens |
|---|---|---|
| **necessary** | the paper states or shows that the fact fails, or changes, without it | kept |
| **scope** | the paper asserts the fact only under it, and says nothing about elsewhere | kept: it is necessary for the fact *as the paper supports it* |
| **unnecessary** | the paper itself states or shows that the fact holds without it | flagged; the checker gives the quote that shows this |

The only change from your rule is that "scope" is treated as necessary, because the evidence for the fact stops at
that condition.

**Complete** means that every necessary or scope condition the paper attaches to the fact is listed. That includes
conditions stated elsewhere in the paper: model assumptions in an earlier section, a figure caption, a table header.
The checker's whole-paper assumption list is what it checks completeness against. Each missing condition must come
with a quote and a paragraph ID, so a script can confirm it exists.

This will also help with the paused control review: all nine flagged gold statements involve a dropped condition.

### 4.5 The memo

Kept.

### 4.6 Stage names, and what the letters meant

The letters were not mine; they come from the first design (WORKFLOW revision 9):
- U = you (user);
- C = calibration;
- A = acquire;
- P = prepare;
- E = extract and check;
- L = library;
- Y = synthesis ("S" was taken);
- X = cross-check;
- BK = backup.

Several of them don't say what they mean, which is your point.

Proposal: **the letter is the phase, and the number is the order within the phase.** You'd write "R3: check facts".
Your stages keep the name "your …" so they are easy to spot.

| Phase | Code: name | Old code |
|---|---|---|
| **F, find** | F1: approve discipline brief (you) · F2: write queries · F3: search and chain citations · F4: scout the web · F5: screen relevance · F6: approve shortlist (you) | U0 (part), new |
| **G, get** | G1: fetch · G2: check copy · G3: prepare paper · G4: index paragraphs | new, A2, P1, new |
| **R, read and check** | R1: extract facts · R2: quote check · R3: check facts · R4: revise facts · R5: blind reread | E1, E3, E4, E5, X1 |
| **L, library** | L1: publish · L2: your fact review · L3: publish your decisions | L1, U3a, L1b |
| **M, memo** | M1: draft memo · M2: citation check · M3: skeptic review · M4: approve memo (you) | Y1, Y2, Y3, U4 |
| **K, keep** | K1: backup | BK |
| **C, calibrate** | C1: calibrate checker | C0 |

The reviews also get names instead of B1–B5: the methods review, the planetary review, and so on.

---

## 5. Agreed so far

- Phase 1 runs as scripts plus single model calls, one discipline brief per discipline, run in parallel. Only the
  scout is an agent.
- Two layers: everything indexed; every fact of core papers extracted and checked. New questions go through Layer 1
  first.
- Our checking, revision rounds and claim review stay as built. The claim review adds paragraph numbers.
- The memo stays.
- Everything Scenario C does that we listed is adopted. The budget cap is based on the usage limit (section 1.6).
- Agents are cut to what each step needs: decision steps become single calls with a schema (back end A or B,
  section 2.6); the foreman becomes the
  status script plus this session.
- The stage runner stays. LangGraph comes later, if at all.
- The disciplines are the union, in your order. The condition rule uses three labels.

## 6. Questions before I write the plan

1. **API key.** Can you get a Claude Console API key, through Caltech or the group? If you can, decision steps use
   back end A and the scout can use the Agent SDK. If you can't, everything runs on your login: back end B, and the
   scout on `claude -p`. Note that API calls are billed per token to whoever owns the key, which is not true of your
   seat.
2. **Stage codes.** Do you accept the phase-letter scheme in section 4.6 (F, G, R, L, M, K, C)?
3. **The condition rule.** Do you accept the three labels (necessary, scope, unnecessary), with scope kept as
   necessary (section 4.4)?
4. **Disciplines 4 and 5.** Do sublimation physics and then material properties come last, and should properties
   be looked up on demand rather than swept (section 4.2)?
5. **Usage cap.** Is 90 % of the five-hour window the right threshold? Should the weekly limit always stop and wait
   for you (section 1.6)?

The plan's order would be:
1. resume the paused methods review: the misjudged claim (c006) and the nine control rulings, then the
   isolation fix and the freeze;
2. run the nine papers as the pilot, on the current engine;
3. build the new pieces (call module, usage cap, paragraph IDs, Phase 1 scripts) and rerun the pilot on them as the
   baseline comparison;
4. scale up one discipline at a time, numerical methods first.
