# Designing agents, decomposing tasks, orchestrating, checking

## What an agent is

A model in a loop with tools, defined by five things:

1. **Role and instructions**: the system prompt (`.claude/agents/<name>.md`).
2. **Tools**: the smallest set that does the job. Every tool removed is a
   failure mode removed.
3. **Input contract**: exactly which files it reads.
4. **Output contract**: one file with a declared schema (JSON with a schema,
   or Markdown with fixed headings).
5. **Stop condition**: done, failed with a reason, or a cap on turns.

Building one means writing those five and testing on 10–20 real cases before
it is trusted. Change one thing at a time and re-run the same cases.

Frontmatter fields used here: `name`, `description` (when to use it),
`tools`, `model` (full ID), and optionally `skills`, `mcpServers`,
`permissionMode`, `maxTurns`. A subagent starts with a fresh context: its own
prompt, the task message and CLAUDE.md, nothing else from the parent session.
Only its final message returns.

## Where agents do well, and where they fail

Good at: transforming material they are given (extracting, summarising,
translating notation, drafting code that will be reviewed) when "done" is
clear and a tool does the arithmetic.

| Failure | Fix |
|---|---|
| Recalls citations, constants or numbers from memory | Retrieve the source; require page and verbatim quote |
| Multi-step arithmetic or algebra | A script or CAS (SymPy) does it; the agent sets it up |
| Checks its own work ("are you sure?" lowers accuracy) | A different instrument checks |
| Drifts on long, open-ended runs | Short jobs with an artifact and a stop condition |
| Locks onto an early wrong diagnosis | Pass only the latest verdict forward, not the history |
| Agrees with a hypothesis stated in the prompt | Ask neutral questions; red-team the result |
| Follows instruction-like text in documents | Treat document text as data; report it |
| Plausible but wrong plan (does not crash) | A checker reads the plan before it runs |

## How big a task should be

A task fits one agent run when all three hold:

- it can be specified in a paragraph;
- its inputs fit comfortably in context, well under the window;
- it produces one artifact that something else can check.

If not, split it, and cut where a check fits between the pieces. Examples:
one paper per extractor run; one derivation step per run with a CAS check;
one source file per code reviewer. Effort scales with difficulty:
fact-finding is one agent with a few tool calls; comparing options is two to
four agents; broad research takes many (Anthropic, multi-agent research
system). Multi-agent runs cost several times the tokens of one agent, so they
are reserved for tasks that need them.

## Orchestrating

The five standard patterns (Anthropic, *Building effective agents*):

| Pattern | Shape | Use here |
|---|---|---|
| Prompt chaining | fixed sequence, a check between steps | derivation → weak form → discretisation |
| Routing | classify, then send to a specialist | triage papers by topic |
| Parallelisation | independent pieces at once, or the same task several times | one extractor per paper; candidate formulations |
| Orchestrator–workers | a lead plans and delegates, workers return artifacts | the literature review |
| Evaluator–optimiser | a producer and a checker loop, with a cap | code until the tests pass, max 3 rounds, then a person decides |

Most dependable systems are workflows: fixed control flow with model calls at
set points. The orchestrator routes and records; it does not do the science.

**When nobody knows the path (Pillar 2).** Use stages with gates instead of
one open-ended agent:

1. Frame the problem and write acceptance tests with the user.
2. Generate two or three candidate formulations in parallel, independently.
3. A fresh-context red team attacks each.
4. The user chooses.
5. Derive step by step; CAS checks each step.
6. Weak form and discretisation, checked symbolically.
7. Implement behind a runtime option.
8. Verify with deterministic numerical tests: manufactured solutions with
   observed convergence rates, energy decay, mass conservation, and limits
   that recover the current model.

Each gate ends in a file and a decision. When a stage fails, the workflow
stops and asks; it never loops "until it passes".

## Checking

Open research, but this policy is defensible and can be measured:

- **Tier 0: deterministic, whenever a formula or test exists.** CAS, unit
  checks, manufactured solutions, Crossref for DOIs, exact quote matched
  against the PDF text, page image for numbers. Same input, same verdict.
- **Tier 1: one independent checker** for judgment. It gets the artifacts,
  not the producer's reasoning, and is asked a different question. Two
  instances of the same model give independence of context, not of judgment.
- **Tier 2: N independent solvers plus an adjudicator.** Only when the answer
  is discrete and comparable and no Tier 0 check exists. A vote cancels random
  error, not shared systematic error, so vary the method, prompt or model, not
  only the sample.
- **The user:** at decisions, and for anything only the user can know.

Measure instead of guessing: plant known errors in a sample of outputs and
record the fraction each checker catches. Use that number to decide whether a
second checker or N = 3 is worth its cost. Report the counts.
