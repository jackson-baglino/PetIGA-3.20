# 24. A math agent that derives the diffuse-domain model: design and build prompt

Draft, 2026-11-15. Part 1 is the design: what the agent is and which features it needs. Part 2 is a
prompt to paste into a new session that builds it. Part 3 lists the choices that are yours. Nothing
here is built yet.

## 1. Design

### 1.1 What the agent does

The agent is a mathematician. It takes the current two-phase model (ice, temperature, vapor; the
equations in `src/assembly.c`, `Residual_A1`) and derives the version in which the grains are
described by a smooth indicator w instead of mesh edges. This is the diffuse-domain method (DDM). It
must produce:

1. the strong form of all three equations;
2. the weak form that a finite-element code assembles;
3. the derivatives needed for the Jacobian;
4. checks that the result is right.

Every step must be either cited from a paper in the nine-paper set or derived from earlier steps,
and it must be possible to check it.

A librarian helps it find things in those papers. The math agent can also open any paper itself.

### 1.2 Two agents, no more

| Agent | Model (pinned) | Job | Tools |
|---|---|---|---|
| `rick-deriver` (or `beth-deriver`) | `claude-fable-5-1` (or `claude-opus-5-5`) | Derives; writes the derivation and its check scripts | Read, Grep, Glob, Write and Edit in its own folder, Bash limited to running its check scripts, Agent (librarian only) |
| `morty-librarian` | `claude-sonnet-5-5` | Finds things in the papers and returns them word for word | Read, Grep, Glob on `library/` only |

The librarian finds and quotes; it does not interpret. The math agent reads the quoted page itself
before it relies on a quote, so a step never rests only on the librarian's word.

**How they talk to each other.** As far as I know, a Claude Code subagent cannot start another
subagent. Check this in the current docs before building. If it holds, the math agent must run as
the session's main agent (`claude --agent rick-deriver`), with the librarian as its subagent.

Lesson from the foreman: a session that is parked into a background job loses `--agent`. Check the
banner at the start of each session.

### 1.3 Features, and why each one is needed

1. **A ledger of steps.** Each step of the derivation has a number and one type:
   - **cited:** a paper, an equation number, a page and a verbatim quote;
   - **derived:** from steps i and j, by a named operation (product rule, divergence theorem,
     integration by parts, ...);
   - **assumption (A#):** stated, with what changes if it is wrong;
   - **choice (C#):** a modeling decision that is yours, for example insulating grains or the
     contact angle.

   The derivation document is written from the ledger, so every equation traces to its source (the
   project rule).
2. **Symbolic checks in scripts (sympy), not in prose.** One script per claim:
   - **Reduction:** setting w = 1 gives back today's equations term by term, compared with
     `assembly.c` line by line.
   - **Conservation:** the weighted ice + vapor equations conserve total water in a closed box. This
     is the divergence theorem with the w weighting.
   - **Energy:** the weighted free energy does not increase for the ice equation.
   - **Strong and weak form agree:** integrating the weak form by parts gives back the strong form.
   - **Jacobian:** symbolic derivatives of each residual term, ready to code; later, the C version is
     checked with `-snes_test_jacobian`.
   - **Units:** every term in a row has the same units (`memory/semantic/units.md`).

   **sympy is not installed** in either venv today. It needs a pinned version in a requirements file;
   that's your decision because of the pin rule.
3. **The sharp-interface limit.** This is the reason the method works. As the grain interface width
   eps_s → 0, the weighted equation must recover the boundary condition the grain imposes: no flux
   for vapor and heat with insulating grains, and the contact angle for ice. The agent does this the
   way Li et al. 2009 and Lervåg & Lowengrub 2015 do, by matched asymptotics, and states the order of
   accuracy. A small 1D numerical check in Python (not the solver) backs it up.
4. **Manufactured solutions.** For each equation the agent writes the source term that makes a chosen
   smooth field an exact solution. Later these become code tests; they are the strongest test a new
   solver term can get.
5. **A map to the code.** A table lists each term of the weak form, the `assembly.c` line it changes
   or the new code it needs, and the run option that switches it on (`-ddm`, `-ddm_eps_s`,
   `-ddm_floor`). The implementer (a later session) works from this table.
6. **Stopping at your decisions.** The agent stops and asks at every choice (C#):
   - insulating or conducting grains;
   - 90° first or any contact angle;
   - the floor υ on w;
   - which equations are weighted by w and which by φ_a = w(1 − ψ).

   It doesn't pick a choice for you and carry on.
7. **Deterministic gates (hooks).**
   - A Stop hook refuses to end the session while a check script fails, or while a ledger step has
     neither a source nor a derivation.
   - A PreToolUse hook blocks writes outside the derivation folder and edits to `src/`.

   These are scripts, not the agent's opinion of its own work.
8. **Independence, so the result can be evaluated.** `docs/planning/18_f2_implementation_plan.md`
   already holds a derivation of the same model (the F2 pore-space form), and 16 and 23 repeat parts
   of it. If the agent cannot read those files, its derivation is independent. A script and you can
   then compare it with doc 18, term by term. That gives the course one clean agent result:
   - did the agent reach the same equations;
   - did it find an error in doc 18;
   - or did doc 18 catch an error of the agent's.

   This also answers the TA's point about a quantitative evaluation better than another
   literature-review pass.
9. **One fresh-context review at the end.** `beth-writer`, in review mode, reads the finished
   derivation and the papers it cites, but not the agent's notes. It reports problems with concrete
   failure cases.
10. **Fixed notation.** A notation table the agent must use, and may not add to without saying so:
    ψ, w, w̃, υ, φ_a, eps, eps_s, σ_s, and γ_ia, γ_is, γ_as.
11. **Reproducibility.**
    - The model is pinned in the frontmatter.
    - The agent files are committed before the first run.
    - Each session's transcript and the commit it ran on are recorded in the derivation folder's
      README.
12. **Deliver in stages, each one checked by you:**
    1. strong form with 90° walls and insulating grains;
    2. the eps_s → 0 limit;
    3. weak form and Jacobian;
    4. the contact-angle term for θ ≠ 90°.

### 1.4 Permissions the agent needs

- **Read:**
  - `src/`, `include/`;
  - `library/papers/*/` (content, claims, cards);
  - `memory/semantic/`;
  - `docs/interface_width_conventions.md`.
- **Write and Edit:** `studies/ddm/derivation/` only.
- **Bash:** only `venv_lunar/bin/python studies/ddm/derivation/checks/<file>.py` (or the venv that
  gets sympy).
- **None:** web access, edits to `src/`, simulation runs.

Note: in this session the auto mode classifier blocked me from listing `library/papers/*/content`.
Paper content will need an explicit allow rule for these agents, which goes in `.claude/`; you
approve that.

### 1.5 Fit with the course deadline

Lab Report 3 asks for a demo of the smallest working piece. Stage 1 alone could serve as that demo:
the strong form of one equation (vapor, say), derived with citations and passing its reduction and
conservation scripts. It fits a day and shows the agent working on the model rather than on the
literature-review line. Building the solver code (plan 23) would then move to Lab Report 4.
That's for you to weigh.

## 2. Prompt for the session that builds the agents

Paste this into a new session in this project.

````text
Goal: design with me, then build, two agents that derive the diffuse-domain (DDM) version of our
model. Read first: docs/planning/24_math_agent_prompt.md (the design), CLAUDE.md, memory/USER.md,
project/PILLARS.md, src/assembly.c (Residual_A1 and Jacobian_A1), include/NASA_types.h, and the
`new-agent` skill. Do not read docs/planning/16, 18 or 23 into the agents' instructions: the
derivation must be independent so we can compare it with doc 18 afterward. You (this session) may
read them.

Agents:
1. rick-deriver (model claude-fable-5-1; or beth-deriver on claude-opus-5-5 if I choose so). Runs as
   the main agent of its own session (`claude --agent rick-deriver`). Derives the DDM model in
   stages: (1) strong form, 90-degree walls, insulating grains; (2) the eps_s -> 0 limit by matched
   asymptotics; (3) weak form and Jacobian derivatives; (4) contact angle theta != 90 degrees.
   Starting point: the three equations in Residual_A1 (ice phase field, heat, vapor), the grain
   indicator w = (1 + tanh(d / (2 eps_s))) / 2, and the nine papers in library/papers/.
2. morty-librarian (model claude-sonnet-5-5). Subagent of the deriver. Given a question, it returns
   matching passages from library/papers/ verbatim, each with paper, page, equation number and file
   path. It never interprets or paraphrases. Tools: Read, Grep, Glob on library/ only.

Features the deriver must have (details and reasons in doc 24 §1.3):
- A step ledger (studies/ddm/derivation/ledger.json): each step is cited (paper, eq., page, quote),
  derived (from steps i, j by a named operation), an assumption A# or a choice C#. The derivation
  document (derivation.md, LaTeX math) is written from the ledger.
- Before relying on any librarian quote, it opens that page itself.
- One sympy check script per claim in studies/ddm/derivation/checks/: reduction to Residual_A1 at
  w = 1, term by term; water conservation in a closed box; free-energy decrease for the ice
  equation; weak form <-> strong form; units; symbolic Jacobian entries. Each script prints PASS or
  FAIL and the inputs it used.
- Manufactured-solution source terms for each equation.
- A table mapping each weak-form term to the assembly.c line it changes and the run option that
  turns it on (-ddm, -ddm_eps_s, -ddm_floor).
- It stops and asks me at every choice C# (insulating vs conducting grains, contact angle, floor,
  which equations are weighted by w vs phi_a = w(1 - psi)). It never picks one and continues.
- Fixed notation: psi, w, w_tilde, upsilon, phi_a, eps, eps_s, sigma_s, gamma_ia, gamma_is,
  gamma_as. Any new symbol is declared in the ledger.
- No web, no edits outside studies/ddm/derivation/, no runs of the solver.

Deterministic gates:
- A Stop hook script that fails while any check script fails or any ledger step lacks a source or
  a derivation.
- A PreToolUse hook that blocks writes outside studies/ddm/derivation/ and any edit to src/.
- sympy pinned in a requirements file (ask me which venv; numpy and scipy must not change).
- An allow rule so these agents can read library/papers/*/ (auto mode currently blocks it).

Evaluation, after stage 1: a script compares the deriver's strong form with doc 18 term by term;
then beth-writer in review mode reads derivation.md and the cited pages (not the ledger notes) and
reports problems with concrete failure cases.

Process: propose the agent files, hooks, permissions and the folder layout to me first, with
anything you would change in this design and why. Build only after I approve. Commit the
definitions before the first derivation session. Pin every model ID. Check in the current Claude
Code docs whether a subagent can start another subagent, and whether `--agent` survives a
background job; tell me what you found and how you checked.
````

## 3. Your decisions

1. **Model for the deriver.** Fable 5.1 (`rick-`) or Opus 5.5 (`beth-`). I'd use the most capable
   model here: a derivation is a long chain where one wrong sign carries through.
2. **Effort.** Leave it at the default, or raise it for the deriver? The project rule (default
   effort, raise it only as a diagnostic) was written for the literature-review line, so it's your
   call whether it applies here.
3. **Independence.** Hide docs 16, 18 and 23 from the deriver? I recommend yes; it is what makes the
   result an evaluation.
4. **sympy.** Which venv gets it, and at what pinned version?
5. **Flowchart.** USER.md asks for an HTML flowchart for every workflow. With two agents, a short
   diagram in this file may be enough; your call.
6. **Lab Report 3.** Use stage 1 of the derivation as the demo, or keep the solver tests T0–T3 from
   plan 23?
