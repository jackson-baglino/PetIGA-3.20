# Project pillars

Started 2026-10-21, at the restart of the agentic-workflow project for ME295
(AI Agents for Scientific Discovery). The first attempt (2026-10-12 to
2026-10-21) is treated as a practice round and archived under
`scratch/2026-10-21_practice_round/` and the `archive/2026-10-21/*` git tags.
This document grows as each pillar is designed. The detailed research problem
will be added when the user presents it.

## Purpose

Use a designed, documented agentic workflow to answer a research question from
the user's own work: modeling the interaction of ice, water vapor (air) and
regolith/sediment. The course grades the workflow and its evidence
(reproducibility 20 %, technical depth 30 %, results 25 %), so how each result
was obtained matters as much as the result.

The standard for every agent and every step: we are scientists and engineers.
The aim is to be factually correct in everything we state and produce, not to
reach an answer.

## Two campaigns

**Campaign A: new model development.** This is the agentic course project,
organized as the pillars below.

**Campaign B: production science on the current code.** The solver on `main`
runs well and has no known problems. It is used, in parallel with Campaign A,
to run simulations on related scientific questions. It is also the fallback if
the new method does not work out. Campaign B does not wait on Campaign A.

## Pilot — Synthesize `docs/`

A small first workflow: read the existing design documents, check their
claims against the code, and write one concise, verified description of the
current model. It builds the reading, writing and fact-checking agents that
Pillar 1 and later reviews will reuse.

## Pillar 1 — Literature review (numerical methods)

Survey the numerical-methods literature relevant to improving the model. The
output is a set of verified, cited findings, not a reading list. The pillar
ends by presenting model options for discussion with the user, unless that
step is placed at the start of Pillar 2. Later reviews (Pillar 3, and others)
may reuse the same agents or need new ones. That is decided when they are
designed.

## Pillar 2 — Model reformulation and implementation

The largest and hardest pillar. Steps, in order:

1. Choose among the model options from Pillar 1, in discussion with the user.
2. Derive the governing equations.
3. Derive the weak form and the discretization.
4. Implement in the solver (behind runtime options, on a `feat/*` branch).
5. Verify and test.

The model will borrow ideas from existing work but is a new formulation and a
new implementation, so nobody knows the exact path in advance. The workflow
must allow for that: planning in stages, checks at each stage, and the user
deciding at the branch points.

## Pillar 3 — Physics literature and science runs

Extend the literature review to physics research and observations. Use it to
define and run a simulation study that produces new, useful scientific
information. Further ideas will be considered only if this works.

## How each pillar is designed

Before any pillar runs, the user and Claude design its workflow together:

- which agents do what, and on which model;
- the tools, skills and MCP servers each agent has;
- how results are checked;
- how the run is recorded so it can be reproduced;
- where the user approves or decides.

Each design is drawn as a detailed flowchart in an `.html` file. A workflow is
built only after the user approves its design.

## Course milestones

| Date | Deliverable |
|---|---|
| Sun 18 Oct | Lab Report 1 (submitted) |
| Sun 1 Nov | Lab Report 2: 1-page proposal, memory decision table, PERMISSIONS.md and settings.json, audit output |
| Sun 15 Nov | Lab Report 3: first working piece |
| Sun 29 Nov | Lab Report 4: end-to-end run on the user's problem, failure modes stated |
| Wed 2 Dec | Final presentation |
| Dec 4 | 4-page IEEE paper with AI Use Declaration; public repository with a reproducible environment and a pinned model ID (the project guideline says Dec 9) |

Lab Report dates are from the Lab Guideline.
