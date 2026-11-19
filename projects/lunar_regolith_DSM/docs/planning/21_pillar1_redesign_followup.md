# Pillar 1 redesign: follow-up answers (round 3)

2026-11-03. This answers your responses to `20_pillar1_redesign_answers.md`: your two new questions first, then your
points in the order of note 20. Section 9 collects what I need from you.

Noted from your reply and taken out of the design: **there is no weekly limit.** The usage cap watches only the
five-hour window.

---

## 1. What is a tool loop?

A model cannot run anything. It can only write text. A **tool** is a function that the *harness* offers to the
model. The harness is the program around the model: Claude Code, the Agent SDK, or our own script. Each tool is
offered with three things:
- a name;
- a description;
- the shape of its arguments.

When the model wants a tool, its reply says "call `Read` with `file_path = p007.png`." Then:
1. **the harness** runs the function, under its permission rules;
2. the harness sends the result back to the model as part of the conversation;
3. the model reads the result and either asks for another tool or gives its final answer.

That repetition is the **tool loop**: model, then tool, then model, then tool, until the model says it is done.

What can be a tool:
- **Built-in Claude Code tools:** Read, Write, Edit, Bash, Grep, WebSearch, WebFetch, and others.
- **Your own Python function**, offered in one of three ways:
  - through an MCP server, a small separate program that serves your functions as tools (works in Claude Code
    and the Agent SDK);
  - as an in-process function in the Agent SDK;
  - as a JSON tool definition in the raw API.
- **A Python script.** This is the sense you described. It is *not* a tool by itself. The agent can run it only if
  it has the Bash tool, and then the tool is Bash: the model writes a shell command and the harness runs it.

**An example from our own line.** Today the claim reader (E1) gets the paper's page images as files. The model asks
for `Read p001.png` and gets the image, then asks for `Read p002.png`, and so on. Every round trip sends the whole
conversation back to the model, so 20 pages means about 20 turns of a growing conversation. **The model decides
which pages it opens and in what order.**

**The single call that replaces it.** The script puts the 20 page images and the numbered paragraphs into **one**
request, and one answer comes back. There is no loop and nothing for the model to choose. That is what "no tool
loop" meant in note 20.

When the loop is worth having: when the next input depends on what the model found. "Search, read the hits, decide
what to search next" is the scout. You can't write that list of inputs in advance, so there the loop is the point.

## 2. What is a schema? Is it our contract?

A **schema** (a JSON Schema) is a machine-readable description of the *shape* of one answer:
- which fields it has;
- each field's type;
- which fields are required;
- which values are allowed.

A shortened example for one extracted fact:

```json
{"type": "object", "additionalProperties": false,
 "required": ["statement", "type", "quote", "paragraph", "conditions"],
 "properties": {
   "statement":  {"type": "string"},
   "type":       {"enum": ["equation", "value", "relation", "method", "observation"]},
   "quote":      {"type": "string"},
   "paragraph":  {"type": "string", "pattern": "^p[0-9]+\\.¶[0-9]+$"},
   "conditions": {"type": "array", "items": {"type": "object",
                  "required": ["text", "quote", "paragraph"], "additionalProperties": false,
                  "properties": {"text": {"type": "string"}, "quote": {"type": "string"},
                                 "paragraph": {"type": "string"}}}}}}
```

When a schema is attached to the call, the API guarantees that the answer has exactly this shape. A missing
paragraph number, an extra field or an unknown type is impossible, not just unlikely. We already use schemas today:
every worker call passes `--json-schema` (`scripts/agentlib.py:162`).

**A schema is not the contract; it is one part of it.**

| | Schema | Contract (`workflows/literature-review/B1_CONTRACT.md`) |
|---|---|---|
| Covers | the shape of one model answer | what the whole line promises: inputs (papers, questions), what each stage must deliver, checks, gates, progress |
| Read by | the API and our validation code | you, me, the stage runner (its progress table) |
| Guarantees | the format | the process |

**A schema says nothing about truth.** A fact can fit the schema perfectly and quote the wrong page. The scripts
after the call check truth: the quote check, the blind reread and the checker.

---

## 3. When an agent is justified (Lectures 9 and 10)

Agreed: this becomes the plan's rule. Each model step in the plan must name the reason it exists. A step with no
reason is a script.
- **Information:** it must see different information, or must *not* see something. The checker must not see the
  reader's reasoning, and the blind reread must not see the claim.
- **Permissions:** it needs different rights. The scout may search the web but can't write to the library.
- **Judgment:** it needs a different kind of judgment, and a separate context keeps that judgment clean. Relevance
  screening is one case; skeptic review is another.
- **Breadth:** the work divides into independent parts. One example is one discipline brief per discipline.

These four reasons justify a **separate context**, meaning a separate call. They don't by themselves justify a tool
loop. So every step answers three questions in order:
1. **Does it need a model at all?** Only judgment does. Otherwise it is a script.
2. **Does it need its own context?** Only for information, permissions, judgment or breadth. Otherwise it joins
   another call.
3. **Does it need a loop?** Only when its next input depends on its own earlier findings. Otherwise it is a single
   call.

Applied to the plan:

| Step | Model? | Own context because | Loop? | Result |
|---|---|---|---|---|
| Search, chain citations, fetch, check copy, publish | no | — | — | script |
| Write queries, screen relevance, choose citations | yes | breadth (per discipline), judgment | no | single call |
| Scout | yes | permissions (web), breadth | **yes** | agent |
| Extract facts | yes | judgment | no | single call |
| Check facts | yes | **information**: must not see the reader's reasoning | no | single call |
| Blind reread | yes | **information**: must not see the claim | no | single call |
| Revise | yes | judgment | no | single call |
| Draft memo, skeptic review | yes | judgment; the skeptic also for information | no | single call |
| Run the line | no | — | — | script (stage runner) |
| Talk with you about the run | yes | permissions (see section 6) | yes | the foreman, interactive |

---

## 4. Making the sign-in question a non-issue (note 20, section 2.6)

### 4.1 Why "try A" runs into your billing rule

You wrote that we can't bill anything extra and should draw on the usage limit. That rule decides between A and B:
- **An API key always bills per token to a Console organization.** The docs say Console API usage does not draw on
  a Claude seat's allowance, and Enterprise plans are not eligible for API credits.
- So **back end A (raw API with a key) is "extra billing" by definition.** That holds unless someone else pays for
  the key: the course or a Caltech Console organization.
- **Back end B uses your login and your five-hour allowance.** It already follows your rule.

**The way to make this a non-issue:** build the call module so that a step doesn't care which back end runs it (the
interface in note 20, section 2.2). Then:
- **B is the default.** It works today and draws on your allowance.
- **A switches on only if** a key turns up that the course or Caltech pays for. You may hear about one: you
  mentioned classmates using SDK agents, and the course may have handed out keys.

Switching back ends changes nothing in the stations. The manifest records which back end ran, and a switch means
recalibrating.

### 4.2 If you do get a key: how to set it so it doesn't break this session

Before you start, a trap: **if `ANTHROPIC_API_KEY` is exported in your shell, Claude Code itself uses the key
instead of your login.** This session and every `claude -p` worker would then bill to the key. So the key must never
go into `~/.zshrc` or any environment-wide setting.

1. **Ask first** (instructor or TA): does the course provide a Claude Console organization or key, and who is billed
   for it? Don't create a key on a paid personal or Caltech account without knowing who pays.
2. **Create the key** in the Claude Console (platform.claude.com): Settings, API keys, Create key. Name it
   `lunar-lr`.
3. **Store it in the macOS Keychain**, typed into a prompt so it never lands in your shell history:
   ```
   security add-generic-password -a "$USER" -s anthropic-lunar-lr -w
   ```
   (It prompts for the password; paste the key there.)
4. **Tell me it is there.** The call module will read it with
   `security find-generic-password -s anthropic-lunar-lr -w` and pass it only to `anthropic.Anthropic(api_key=...)`.
   It never goes into an environment variable, a file in the project, or a log.
5. **Test.** I run one tiny call: about 20 tokens, billed to the key's owner. It shows whether the key works and
   which organization it bills.

### 4.3 Updating the `anthropic` package: a conflict

`venv_workflows` has `anthropic` 0.46.0 only because the PDF converter Marker pulls it in. Marker 1.10.2 pins
`anthropic >=0.46.0,<0.47.0`, and the newest Marker (2.1.0) pins `anthropic >=0.69.0,<1`. Marker is our fallback
parser and the units workflow's OCR. So the current `anthropic` (2.4.0, the version that matches the current docs)
**cannot go into `venv_workflows` without breaking Marker's requirement.**

Options:
- **(a) A small separate environment, `venv_llm`.** It holds only `anthropic` 2.4.0 and `claude-agent-sdk`,
  pinned in `requirements-llm.lock`. The call module and the scout run in it. There is already a precedent: the
  paper parser has its own environment (`venv_docparse`) for the same kind of reason. **Recommended.**
- (b) Upgrade Marker to 2.1.0 and `anthropic` to 0.125.0, the last 0.x release. This changes a parser used by two
  workflows, which then have to be checked again.
- (c) Install 2.4.0 anyway and ignore Marker's pin. Not recommended: `pip check` would fail, and Marker's own
  Claude code path could break.

I haven't installed anything yet. Choose (a), (b) or (c).

---

## 5. The Agent SDK compared with what we have

You're right that Claude Code agents defined in `.md` files can already have their own tools (through MCP servers)
and hooks (shell commands). Here is how the three options differ:

| | Claude Code agents (`.md` + `claude -p`) | Claude Agent SDK | Raw Anthropic SDK (`anthropic`) |
|---|---|---|---|
| What it is | the Claude Code program, configured by a markdown file | **the same Claude Code program, driven from Python** | a thin Python wrapper over the Messages API; no Claude Code at all |
| Who runs the tool loop | Claude Code | Claude Code (the SDK starts it as a subprocess) | nobody, unless you write the loop |
| Built-in tools (Read, Bash, WebSearch) | yes | yes | no |
| Your own tools | via an MCP server (a separate program) | **Python functions in the same program** | JSON tool definitions you execute yourself |
| Hooks | shell commands, configured in files | **Python functions**, with access to your program's state (for example the run ledger) | you write them into your loop |
| Sign-in | your login | probably your login (to test) | API key only |
| Our use today | every worker, through `agentlib.py` | none | none |

**The honest comparison: `agentlib.py` is a small homemade Agent SDK.** It is 276 lines that start `claude -p`,
parse its stream of JSON messages, apply the agent's `.md` settings and log every call. The real SDK does the same
job with more polish:
- tools and hooks as plain Python, with no separate MCP server process;
- typed messages;
- the ability to interrupt a run or approve each tool call from code;
- sessions that can be resumed or forked.

It doesn't give us a capability we lack. Its value is convenience, and **learning the current standard**, which
your classmates are also using.

### Your picture of it

You described two kinds of agents:
1. loose `.md` agents that respond to changing situations;
2. tight SDK agents for repetitive tasks: a Python file calls an agent with specific tools, and an orchestrating
   agent such as the foreman calls that Python file.

Three corrections:
- **The loop is inside the agent, not in your Python file.** Your Python file calls the agent once per task, for
  example once per discipline. Inside that one call, the agent loops over its tools. A loop in your Python, such as
  "for each discipline", is just a script calling a function.
- **Loose versus tight is a matter of configuration, not of `.md` versus SDK.** Both can be tight: few tools,
  narrow permissions, a schema, turn limits. Both can be loose. Tightness comes from what you give the agent.
- **The orchestrator stays code (the stage runner), not an agent.** A foreman agent may *start* the runner and talk
  with you about it (section 6), but the runner decides what runs next.

### The rule I suggest

- **Interactive helpers you talk to** (writer, design critic, foreman): `.md` agents in Claude Code. They are the
  loosest, as you said.
- **Line workers that need a tool loop** (the scout, and only the scout): the **Agent SDK**, if your login works
  with it, which gives you the learning value. Otherwise `claude -p` through `agentlib.py`, as now.
- **Line steps that are single calls** (everything else): the call module, back end B (or A).

Over time `agentlib.py` could be replaced by the SDK, but that isn't needed for the redesign.

## 6. A foreman agent, or the foreman in the main session?

**Keep a foreman agent, as the interface you talk to while a run executes.** Note 20 was too blunt. The point of
"replace the foreman with a script" was that the *control* belongs to the script: what runs next, which worker
gets which paper, when to pause. Talking with you about the run is a legitimate agent job, and a separate foreman is
justified on two of the four reasons:
- **Permissions.** In the execute phase the definitions are frozen. The foreman can run the line, read status and
  explain, but it cannot edit code, record your approvals or read paper text (its hooks enforce this). The main
  session has broad rights because it builds things. Keeping the two apart is what stops "just a quick fix"
  mid-run.
- **Information.** The foreman starts clean, with the run's status and the plan, not our long design conversation.

So:
- **the main session** (this one, CLAUDE.md) is for designing and building;
- **the foreman session** (`beth-foreman`, which already exists) is for running, with you present;
- **the stage runner** controls the line.

The foreman's job list gets rewritten for the new design: run, read the pause report, explain it to you, give you
the exact commands for your decisions, park ideas in the backlog. It drops anything the runner already does.

## 7. Flags and OpenAI

**OpenAI.** It is only in note 19 because you asked whether scientists at Anthropic, OpenAI or DeepMind prefer a
framework. The project uses Claude only, so nothing about OpenAI affects the build.

**`--safe-mode` instead of `--bare`: yes, for our purpose.** From `claude --help` (2.1.317):
- **`--safe-mode`** turns off all customizations: CLAUDE.md, skills, installed plugins, hooks, MCP servers, custom
  commands and agents, output styles, and more. **Sign-in, model choice, built-in tools and permissions work
  normally.** Admin policy settings still apply.
- **`--bare`** also skips those, plus LSP and plugin syncing, but it accepts only an API key.

So `--safe-mode` gives the isolation we want while keeping your login. Two cautions:
- **Auto memory isn't named in the help text** ("and more"). The research agent read the docs as saying safe mode
  covers it. I'll confirm with a canary: plant a marker in memory and check that a worker can't see it.
- **Safe mode turns off hooks.** Our workers rely on `--tools`, the permission file passed with `--settings`, and
  the macOS sandbox, not on hooks, so nothing is lost. The canary test also checks that `--settings` permissions
  still apply.

This goes with the auto-memory item already in the backlog, before the freeze.

## 8. Speeding up the 12 minutes per paper

Where the 12 minutes go, and what each change would do:

| Stage | Min | Why it takes that long | Faster by |
|---|---|---|---|
| prepare paper | 4.4 | the local parser reads every layout block as LaTeX or text on this Mac's GPU | runs **ahead**: all papers can be prepared before reading starts, for example overnight; it never waits on the model |
| extract claims | 1.4 | tool loop: one Read per page, each round trip resends the conversation | single call: all pages at once |
| check claims | 3.2 | several separate agent sessions (checklist, blind transcription, verdicts per 4 pages), run one after another | single calls, and the groups of 4 pages run **in parallel** |
| revise claims | 2.6 | an agent session per failed claim | single calls, in parallel; only failed claims |

The bigger gain is running **many papers at once**. Twelve minutes per paper is the time for one paper going
through the line. Eight papers in parallel take about the same wall-clock time as one, so about 1.5 minutes each.
The real limit is then the five-hour usage window, which the cap of note 20 (section 1.6) handles, not the clock.

I won't promise numbers before measuring them: the pilot rerun on the new engine gives the before-and-after. One
lever needs your decision. Lower effort on simple steps, such as relevance screening, would be faster. But your
standing rule is that workers run at the default effort, and effort changes only as a diagnostic. I suggest keeping
that rule for now.

---

## 9. What I need from you

1. **`anthropic` package:** (a) separate `venv_llm` (recommended), (b) upgrade Marker, or (c) ignore the pin
   (section 4.3)?
2. **API key:** ask the course whether a Console key is provided and who pays. Until then the plan uses back end B
   on your login, which already follows your no-extra-billing rule (section 4.1).
3. **Foreman:** keep `beth-foreman` as the interactive run interface, with control in the runner (section 6)?
4. **Agent SDK for the scout only,** if your login works with it (section 5)?

Everything else in note 20 is agreed, including the stage codes, the condition rule, the discipline order and
dropping the weekly limit. Once you answer these, I'll write the two plan documents.
