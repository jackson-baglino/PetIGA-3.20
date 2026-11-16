# 25. The model line: answers, findings, and the literature-review update

Written 2026-11-16 with the approved model-line plan. The design itself is
`workflows/model-development/WORKFLOW.md` and its flowchart. This note holds
what the design rests on. Tags: **[checked]** means read in a file or run, with
the place given; **[inferred]** means judgment; **[unverified]** means from the
user's dictated notes and not yet checked. Supersedes doc 24.

## 1. Answers to the user's questions

### 1.1 Lean: not for this project [inferred]

- **What Lean proves:** Lean proves a stated theorem inside Mathlib's formal
  analysis. Mathlib's support for Sobolev spaces, weak forms and matched
  asymptotics is thin, so formalizing one weak-form derivation would take longer
  than the course has left.
- **What it misses:** our real error sources are code that does not match the
  math, sign errors, missing weights, Jacobian mismatches and unit slips. A Lean
  proof addresses none of them.
- **What the line uses instead:**
  - independent SymPy checks of every derivation step;
  - manufactured solutions with a measured order of convergence;
  - finite-difference Jacobian tests;
  - conservation and energy tests on the running code.
- **If one claim ever needs a real proof,** the installed `math-proof` plugin
  writes a checked informal one.

### 1.2 Why the B1 foreman is not reused

- **Its job is now the line script's.** It ran a long interactive session that
  started `run_b1.py`, reported and stopped at gates.
- **It is tied to B1.** Its Bash whitelist and guard hook are written for B1,
  and it is frozen in b1-v2.
- **It is the layer that stalled** [checked: memory
  `feedback_b1_line_foreman_only.md`; transcripts `f83f89f4`, `9accda1c`]:
  - it lost `--agent` in a background job and after `/clear`;
  - it could not wait on a job;
  - it advised an approval that the line refuses.
- **The new foreman** (`beth-model-foreman`) does one thing. It triages a
  failure from its evidence, as one headless call with fresh context.

### 1.3 Corrections to doc 24

- **Nesting.** Doc 24 assumed a subagent cannot start another subagent. The
  current Claude Code docs say nesting is allowed, up to three levels
  (code.claude.com/docs/en/sub-agents, read 2026-11-16). The design does not
  depend on nesting either way.
- **Fable.** Fable is not available on our enterprise account [checked
  2026-11-16]. A `claude -p` request for `claude-fable-5-1` returned success
  but was served by `claude-opus-5-5` (`modelUsage`). It was a silent
  fallback, not an error. A daily SessionStart hook
  (`~/.claude/hooks/fable_check.sh`) now records this.

## 2. Why agents stopped before, and what changes [checked]

The headless worker calls were reliable. 222 B1 transcripts show no turn-limit
errors and no permission denials. The failures were all in the interactive
orchestration layer:

| Failure | Evidence | Design answer |
|---|---|---|
| `--agent` lost in a background job or after `/clear` | memory `feedback_b1_line_foreman_only.md`; transcript `f83f89f4` | No long-lived agent session. The line script owns the loop and the state. |
| An agent ends its turn while waiting on a job | transcript `9accda1c` ("the rules for this session don't allow a wait command") | Agents never wait. Solver runs are foreground, under a timeout. |
| Five-hour usage limit hit in long Opus work | `.claude/ACTIVITY_LOG.md:39`; sessions `41462a50`, `78c4a8ed` | One call at a time. On a limit error the script stops cleanly and prints the resume command. |
| Context compacted 5× in long main sessions | sessions `bd0b31b7`, `41462a50` | Fresh context per call. Hand-off is through files only. |
| Silent refusals (WebSearch in `dontAsk`; `Write(...)` rules) | `RETROSPECTIVE.md:341-348`; `project/LOG.md:179-181` | Each agent's permission profile is dry-run before a real run. The `agentlib.py` audit flags any denial. |
| "Done" when not done | Claude Code best-practices page | A Stop hook runs the stage's check script, and the output schema requires evidence. |

## 3. Solver facts that shape the performance work [checked]

- **The convergence test is close to empty.** In SI units the residual blocks
  are about 1e-13 (ice), 1e-14 (T) and 1e-24 (vapor). `-snes_atol` is 1e-6,
  so every field passes its absolute test at Newton iteration 1. Sources:
  `inputs/solver.opts:130-131, 171-180`; `src/snes_convergence.c:140-145`.
  Nondimensionalization fixes this.
- **dt control.**
  - The Newton-count rule is hand-patched into PETSc (`alpha1.c:197-288`,
    marked `//AM`). It grows dt by 1.334 whenever Newton takes 5 iterations
    or fewer, which is always.
  - So in practice dt is limited by three things only:
    - `-dtmax` = 2·tau_sub, set by hand;
    - the interface CFL applied after the fact (Δφ ≤ 0.2 per step,
      `monitoring.c:485-501`);
    - the bounds rollback.
  - No diffusive, thermal or vapor time scale is used a priori.
  - Because the rule lives in a patched PETSc, a stock PETSc would behave
    differently [inferred].
- **Solver.** Newton with basic line search, BiCGStab, ASM with ILU(3)
  (`inputs/solver.opts`). Nothing has been compared.
- **Rank sizing.**
  - `scripts/lib/alloc.sh:37` targets 200k DOFs per rank. That value belongs to
    the run coupled with the homogenization model (`effective_thermal_cond`), so
    it does not apply here (user, 2026-11-16).
  - This line uses 20k–100k DOFs per rank and at most 8 of the 12 local cores,
    leaving cores for post-processing. The machine has 12 cores and 64 GiB.
- **Verification.**
  - No manufactured-solution or convergence tests exist.
  - `-test_wall_jacobian` checks dR/dU only, not the shift·dR/dU_t terms
    (`lunar_main.c:1843-1930`).
- **Thin-interface correction: the code contradicts itself.**
  - The correction is ON by default (`lunar_main.c:690`).
  - The option's help text says "Default 0", and that the terms inflate β by
    1.22× for one-sided vapor (`lunar_main.c:691-696`).
  - Memory says the user wants it ON.
  - The derivation must settle this for the new model. It also touches paper
    2's frozen model; that is only flagged here.
- **Earth-air constitutive laws.** ρ_vs uses ρ_air·0.62·P_vs/(P_atm − P_vs)
  (`material_properties.c:180`), and D_v uses Massman (T/T0)^1.81
  (`material_properties.c:135`). Neither holds in vacuum. M1:Scope decides when
  lunar laws come in.

## 4. The user's dictated notes, read skeptically

The notes were dictated, so typos are expected and ignored.

**Adopted:**
- a spec first: problem, mechanisms in and out, notation, sign conventions and
  units (M1:Scope);
- scaling and nondimensionalization before the equations;
- the generator never checks its own work (the verifier writes the derivation
  checks as well as the tests);
- an assumptions ledger;
- tests locked before the code exists, with tolerances set by the user;
- building in layers;
- regression to the two-phase model where ψ ≡ 1;
- separate refinement in eps and the grain width;
- a ξ sensitivity test;
- smooth floors on products, not hard clamps. Today's vapor weight
  max(1 − φ, air_lim) is a clamp [checked, `assembly.c` vapor row];
- the user's time goes to claims that are both disputed and important.

**Not adopted now, and why:**
- **Embedding RAG, citation-graph expansion, historian and scout agents, Elo
  tournaments and cross-model extraction:** these would rebuild B1. Keyword
  search is enough for this library.
- **The FEniCS port, surface diffusion, and the variable-α thin-interface
  model:** separate projects.

**Claims to check before anything relies on them** [unverified]:
- the synthesis note's paper summaries: F&P Eq. 3.9a/b and App. A–C; K&P
  Eqs. 42–44; M&F "Boyer–Lapuerta", "eps < 1e-8 m but uses larger", Table S3;
- "Moure & Fu SI is equivalent to K&P";
- every reference marked "(verify)";
- K&P Sec. III B lamellae as a regression case.

**Settled by the user (2026-11-16):**
- **Notation.** ψ is the fixed grain indicator, and φ ∈ [0, 1] is the ice
  phase field in the pore. Then φ_i = ψφ, φ_a = ψ(1 − φ), φ_s = 1 − ψ, and
  ψ̃ = ψ + υ is the floored indicator. eps is the ice interface parameter and
  eps_s the grain-edge width. This replaces doc 18's (w, ψ).
- **Wetting term.** The contact-angle term is known and correctly implemented
  in the current model. Only its conversion into the new formulation is
  derived.
- **Mixed conductivity.**
  - The rule is the arithmetic–harmonic (tensorial) blend of Nicoli, Plapp &
    Henry, Phys. Rev. E 84, 046707 (2011), doi:10.1103/PhysRevE.84.046707
    (Crossref-checked 2026-11-16).
  - Regolith k and c_p come either from the literature for the lunar
    composition or from a proxy material with well-known properties. The user
    chooses at G:Scope.
- **Comparison test.** The new model is checked against today's phase-field
  model with no sediment, on the neck width of two sintering grains.
- **AGU is separate.** Its runs are on the HPC (paper 2); this line runs locally
  (paper 3).

## 5. The literature review (B1): an honest update

**State at close** [checked: `run_b1.py status`, run `20261105-pilot2`,
2026-11-16]:
- **Papers and claims.** Nine papers were prepared, read, checked and
  published, with 207 new claims and 16 new versions of Kaempfer & Plapp
  claims. Kaempfer & Plapp 2009 and Moure & Fu 2024 (the published copy) are
  both in the library.
- **X1 (blind reread).** It waits at its gate. All 19 disagreeing reads have a
  decision (13 right, 2 wrong, 4 unsure), but the gate is not approved.
- **Disputed claims.** 12 claims from E5 are undecided. Three of them are on
  the two papers the model leans on most: Kaempfer & Plapp c001 and c012, and
  Moure & Fu c011.
- **The Q1 memo draft fails its own citation check.** 39 of 179 sentences are
  uncited.
- **The checker-versus-baseline evaluation is complete and frozen.** Its
  control rate is optimistic, because the controls were not held out (doc 23
  §2).

**What it delivered:** a checked, page-cited claim library for nine papers,
and one clean agent evaluation.

**What it did not deliver:** a finished, cited answer to its question.

**What it cost:** about 17 working days against a 3–5 day budget (doc 23 §2).

**From now on**, the review is a library used on demand, not a line:
1. **B1 is closed** as it stands (`B1_CONTRACT.md`, section "Closed").
2. **The three open claims the model rests on are decided first,** in the
   user's own terminal (`u3_review.py`). The rest stay open.
3. **The math agent searches the library** with
   `scripts/model_development/library_search.py`, a deterministic keyword
   index. Before relying on a passage, it opens the page itself.
4. **Every cited step of the derivation becomes a claim** for the user to
   confirm at the next gate. The quote must match the page (`library_search.py
   quote`). This puts the doc 24 rule into practice: report back before using
   what a paper says.
5. **New papers enter content-only** through
   `scripts/model_development/add_paper.py`. It runs the same MinerU extraction
   B1 used and makes no claims. The first is Nicoli et al. 2011.
