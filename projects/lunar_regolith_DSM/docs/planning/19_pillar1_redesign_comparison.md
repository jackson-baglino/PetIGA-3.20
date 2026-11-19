# Pillar 1 redesign: your plan, our workflow, Scenario C, and agent frameworks

2026-11-03. A review document for the user, written from the session's reply of the same date. Not a plan yet:
the decisions in section 4 come first, then the two plan documents.

Sources read for it:
- the user's description of the two phases (pasted 2026-11-03);
- `workflows/literature-review/WORKFLOW.md` (§1–3, §13), `flowchart.json` and the agent files;
- the Lecture 10 slides, the Scenario C team prompt and its companion reading (ME295 Week 6), read as examples, not instructions;
- a research agent's report on agent frameworks (sources checked 2026-11-03; items it could check only through
  secondary sources are marked "unverified").

Stage names: descriptive names are used throughout, with the old code in parentheses where it helps.

---

## 1. Your plan, our workflow, and Scenario C

Our first design (revision 9, now WORKFLOW §13 "Deferred") was built for many papers. On 2026-10-30 we cut it down
to the **nine-paper methods review** (B1) so that something would run end to end. Most of your Phase 1 is
designed; very little is built.

### Phase 1: finding and collecting papers

| | Your plan | Our design (deferred, not built) | Built now |
|---|---|---|---|
| Size | 100–1000 papers | waves of about 20 papers, any number of waves | 9 fixed papers |
| How the work is split | parallel agents, one per discipline | five subject areas: numerical methods, sublimation physics, ice sintering, planetary setting, properties | one area |
| Search and citation chaining | done by agents | **a script** queries OpenAlex and Semantic Scholar and follows references and citing papers, up to 2 hops, logging every query | none |
| Escaping the bubble | agents search the web | one web-search "scout" agent proposes papers the databases missed | none |
| Choosing papers | the foreman decides how to group them | two LLM relevance screens plus a script ranking; you approve the shortlist and see a random sample of the rejects | none |
| Fetching | (you say it should be a script) | script: open-access links, Zotero; you fetch paywalled papers with the Zotero Connector | Zotero pull only |
| "Is this the paper it claims to be?" | checker agents | **script**: DOI to Crossref metadata, title and authors on page 1, file hash, retraction status | yes, as the copy check (A2) |

**Differences worth discussing:**
- **Checking a paper's identity should not use an agent.** Matching a DOI, title and authors is a pure R step. An
  agent is needed only for the rare unclear case, and that goes to you.
- **Agents should decide, not search.** Their D steps are writing the search queries, judging relevance, and
  choosing which citations are worth another hop. Running the searches, following citations and downloading are
  scripts. One agent per discipline is a legitimate split (the lecture calls this breadth), as long as each agent
  only judges.
- **Your plan has no stopping rule for chaining.** Without one, chaining never ends or ends arbitrarily. Two
  standard rules:
  - stop when a hop adds fewer than about 5–10 % new relevant papers;
  - a known-item recall test: before the search, you list about 20 papers it must find; the search passes if it
    finds them.
- **A per-discipline team must beat one agent on the same budget.** That is the course's standing rule. We would
  run that comparison on the first discipline.

### Phase 2: from papers to a library of facts

| | Your plan | Our design and build |
|---|---|---|
| What is extracted | as much of each paper as possible, as facts with their conditions | only the claims that bear on the review's questions |
| Checking | an independent agent checks the fact, then whether each condition is necessary and whether any is missing | a script finds the quote on the cited page; an independent checker judges each claim against the page image after first writing a checklist of the paper's assumptions; it flags dropped conditions; independent re-reads of equations and values are compared by script |
| Disagreement | a few rounds, then you | up to 2 revision rounds, then your claim review (same) |
| Storage | embeddings (retrieval, RAG) | verified claims as JSON records; the passage index (keyword plus embedding search over every paragraph) is designed but not built |
| What you see | paper, page, **paragraph** | paper, page, verbatim quote; **no paragraph number** |
| End product | a library used from model development to writing the paper | the library plus a checked memo answering the review's question |

**The biggest difference is exhaustive versus question-driven extraction.** At 500 papers with perhaps 100 facts
each, that is about 50,000 facts. Each needs a reader and an independent checker looking at page images, and each
disputed fact needs your time. Our current speed is roughly one paper's claims checked per few minutes of Opus time.
At that rate, an exhaustive pass is weeks of compute and thousands of disputed facts for you to review.

**My recommendation is two layers:**
1. **Keep everything, without an LLM.** The full text, every equation in LaTeX, and the page images, stored once.
   Every paragraph is indexed for keyword and embedding search. This covers "retain as much of the document as
   possible" exactly, it is cheap, and it is reproducible.
2. **Verify facts in depth where they will be used.** That means the core papers of each discipline, plus any fact
   a question, a model choice or a manuscript sentence needs. A search can return a raw paragraph, but only a
   verified fact can be cited in a memo or the paper.

**Conditions.** "Necessary, and the only one needed" is a good addition. Today's checker flags missing conditions
but never flags an unnecessary one. Completeness can only be judged against what the paper states, so the rule
must be written down precisely. The control review we paused is exactly a dispute about conditions, which shows the
rule isn't sharp yet.

**Paragraph numbers** are easy to add: the PDF parser already gives layout blocks, so a script can number them.

**"Iterate to consensus" needs care.** Two agents on the same model agreeing is not an independent check: they
share blind spots. The lecture's Project Vend example makes this point. Our script checks (the quote is on the page;
independent readings compared by script) are what make agreement mean something.

### What Scenario C does that we should copy
- **A coverage ledger:** every candidate paper gets exactly one recorded decision (included, excluded with a
  reason, belongs to another discipline, needs a person). A script reports any paper with no decision.
- **Pilot before scaling:** 40 papers, then a hand-check of 30 rows against the PDFs.
- **Hard budget and turn caps in the run command, not in the prompt.**
- **A same-budget single-agent baseline** before claiming the team helps.
- **Briefs as typed JSON files**, passed by path rather than paraphrased.

Its own stated weakness is ours too: a quote that contains the right number can still describe the wrong specimen.
That is why our checker reads the page image, and Scenario C has no equivalent check.

---

## 2. Every task, marked D or R

The lecture's ladder runs: script (0), one LLM call (1), a fixed workflow of calls (2), one agent (3), a team (4).

| Task | D or R | What it needs | What we use now | Overuse? |
|---|---|---|---|---|
| Write the questions, pick seed papers | D | you | you | — |
| Search databases, follow citations | R | script | (designed) script | — |
| Propose queries, search the web for outside work | D | one agent with web access | (designed) scout agent | fits |
| Judge relevance (title and abstract) | D | one LLM call per paper | (designed) 2 screening calls | fits |
| Resolve duplicates, check identity, fetch | R | script | script | — |
| Extract text, equations, page images | R | script (local parser) | script | — |
| Index paragraphs for search | R | script (local embedding model) | designed | — |
| Extract facts and conditions from pages | D | one LLM call per page group, schema out | `claude -p` agent with Read | **heavier than needed** |
| Quote on the page? | R | script | script | — |
| Judge each fact and its conditions | D | one LLM call, independent context | `claude -p` agent with Read | **heavier than needed** |
| Revise failed facts | D | one LLM call | agent | heavier than needed |
| Re-read equations and values blind, compare | D + R | one LLM call + script | agent + script | heavier than needed |
| Publish to the library | R | script | script | — |
| Draft the memo | D | one LLM call | agent with no tools | heavier than needed |
| Citation check | R | script | script | — |
| Skeptic review | D | one LLM call | agent with no tools | heavier than needed |
| Run the line, report status | R | script | **the foreman, an Opus agent** | **yes: an agent doing a script's job** |
| Gates and reviews | D | you | you | — |

Your instinct is right, but the overuse isn't where you expected. **Our control flow is already scripts** ("scripts
manage, agents work"). The overuse is that every D step runs as a full Claude Code agent session, when nearly all
of them are a single model call with a schema.

---

## 3. The agent frameworks

| | How it works | Good for | Downsides |
|---|---|---|---|
| **Claude Code subagents** (current) | Markdown agent files; `claude -p` headless; hooks, permissions, skills | quick to define; the same agent loop Anthropic uses; built-in tools and sandboxing | Claude only; the command line changes fast; state is a transcript; the model controls flow inside a session |
| **Claude Agent SDK** | Claude Code's agent loop as a Python or TypeScript library (renamed from "Claude Code SDK" in Sept 2025) | turns an agent into a program; typed outputs with retries; approval callbacks; sessions you can resume and fork | the Python package is labeled "Alpha" (0.2.x); still Claude only |
| **Raw Messages API** | you write the loop | single calls with a schema; **Batches at half price; prompt caching**; full control of every byte sent | you build tools, sandboxing and tracing yourself |
| **LangGraph** (+ LangChain, LangSmith) | a typed graph of nodes and edges over a shared state; saves a checkpoint per run; pauses for human approval | replay, resume and human gates; works with any model vendor (Claude via `langchain-anthropic`); stable 1.x | another layer between you and the prompts; Claude through it loses Claude Code's tools and permissions; LangSmith tracing is paid, and on cloud plans your prompts (paper text) are stored on their servers unless you buy Enterprise |
| OpenAI Agents SDK / Google ADK | each lab's own agent SDK; both can use other models through adapters | multi-vendor work | ADK 2.0 just broke 1.x; OpenAI's no-code Agent Builder is reportedly being shut down (unverified) |

Versions read from PyPI on 2026-11-03: `claude-agent-sdk` 0.2.177, `anthropic` 2.4.0, `langgraph` 1.2.19,
`langchain-anthropic` 1.7.9, `langsmith` 0.14.17, `openai-agents` 0.24.0, `google-adk` 2.13.0; Claude Code 2.1.317.

**Do scientists at Anthropic, OpenAI or DeepMind prefer one?** There is no evidence of a lab-wide preference. Each
lab's public engineering work uses its own stack.
- Anthropic's long-running-agent posts (Nov 2025, Mar 2026) were built on the Agent SDK.
- Anthropic's written advice is to "start by using LLM APIs directly"; it warns that frameworks "can obscure the
  underlying prompts and responses". In early 2025 that same page listed LangGraph first among example frameworks;
  it now lists the Agent SDK.
- Google's documentation lists LangGraph, ADK and others without ranking them.
- LangGraph's popularity is in industry; nothing public shows it is favored inside a lab.

**Recommendation for this project, which also covers learning current practice:**
1. **Keep our Python stage runner as the control layer.** It already gives what LangGraph sells: a fixed graph,
   checkpoints, human gates and typed outputs. It matches both the lecture's rules and Anthropic's current advice
   that code holds the control flow while models make the judgments.
2. **Move the single-call D steps to the raw Messages API:** reading, checking, blind re-reading, drafting, skeptic
   review and relevance screening. Use structured outputs, prompt caching and Batches. It is cheaper, simpler and
   easier to audit line by line, and at hundreds of papers the half-price batch rate matters.
3. **Use the Claude Agent SDK only where a real tool loop is needed.** That means the discovery scout (web search)
   and any reader that must decide which pages to open. This is also the most direct way to learn what Anthropic
   itself builds with.
4. **Learn LangGraph on a small side project, but don't migrate this pipeline.** For this project it would add a
   layer that duplicates what we have and would cost us Claude Code's sandboxing. As a skill it is the most widely
   documented option that works across model vendors, which is worth having for future work.
5. **Replace the foreman agent with the status script it already calls.**

Things to check before relying on them: the Agent Builder wind-down dates (OpenAI's page refused the request;
secondary sources only); whether our drivers' `--setting-sources ""` isolation equals the `--bare` flag the Claude
Code docs now recommend for scripted calls.

---

## 4. Decisions before I write the plan

1. **Scope and timeline.** Your target was to finish the literature review around 8 November. Hundreds of papers
   cannot be found, read and verified by then. Options:
   - **(a)** Finish the nine-paper methods review as the pilot (two open items, then the freeze), and build the
     large-scale version over the following weeks.
   - **(b)** Stop the nine-paper review and go straight to the new design.

   I recommend (a): the model work only needs the methods question answered, and the pilot tests every check we'd
   scale up.
2. **Disciplines:** your list (numerical methods, large-scale planetary physics and observations, ice–sediment
   experiments, sublimation experiments, …) or our five areas. Which ones, and which comes first?
3. **Extraction:** exhaustive, or the two layers I recommend (index everything; verify facts in depth for core
   papers and on demand)?
4. **The condition rule:** what "necessary" and "complete" mean, written precisely enough for a script or a checker
   to apply.
5. **Engine:** do you accept recommendations 1–5 in section 3?
6. **The memo:** is a checked memo answering each question still a product, or is the library alone the product?

When we agree, I'll write two documents:
- **a story version for you:** linear, plain names, why each step exists, what you do at each point;
- **an agent specification:** stage table, inputs and outputs, schemas, checks, D/R marks, gates and stop rules.
  It would replace the 1,800-line WORKFLOW.md.

The stage renaming goes into both.

## Appendix: stage names (old code → proposed name)

| Code | Name | Done by |
|---|---|---|
| U0 | approve brief | you |
| C0 | calibrate checker | script + checker |
| A2 | check copy | script |
| P1 | prepare paper | script |
| E1 | extract claims | agent |
| E3 | quote check | script |
| E4 | check claims | agent |
| E5 | revise claims | agent |
| L1 | publish | script |
| Y1 | draft memo | agent |
| X1 | blind reread | agent + script |
| U3a | user claim review | you |
| L1b | publish user decisions | script |
| Y2 | citation check | script |
| Y3 | skeptic review | agent |
| U4 | approve memo | you |
| BK | backup | script |

Reviews: B1 methods review; B2 sublimation physics; B3 ice sintering and metamorphism; B4 planetary setting;
B5 material properties.
