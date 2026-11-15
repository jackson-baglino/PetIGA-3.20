# Plan for the rest of ME295: stop growing the agent pipeline, build the diffuse-domain model

Written 2026-11-15 (Sun), after the TA's feedback on Lab Report 2. Two views went into it:
the main session's, and an independent reviewer agent's, started with fresh context and
given only the files, your comments and the TA's (§8). Tags: **[checked]** means read in a
file or computed by a script, with the place given; **[inferred]** means judgment.

## 1. The decision

You were right, and the TA says the same thing: the literature-review line (B1) is bigger
than this course needs, and its size costs more than it gains. We stop building it. What
it has already produced is kept and reported as it is:
- nine papers read and checked, with 230 claims published;
- one finished, frozen comparison of the agent against a script-only baseline.

From today the work is the model: the diffuse-domain formulation of ice, vapor and heat in
the pore space between regolith grains. This is the F2 model, derived on 2026-10-18 in
`docs/planning/18_f2_implementation_plan.md` and never implemented. It will be:
- built in this solver behind a runtime option, on a `feat/*` branch;
- tested on a handful of small local cases, each with a known answer;
- shown as the Lab Report 3 demo.

Agents stay at the scale of one helper. This session builds and runs, one fresh-context
reviewer audits the code against the derivation, and deterministic tests decide.

## 2. The TA's three points, answered with what we have

**Point 1: too many agents (the foreman is redundant; the blind reread X1 and the checker E4 could merge).**
Agreed.
- **The foreman** added a failure mode (on 2026-11-04 a session silently lost its agent
  setting) and no check that `run_b1.py` does not already make.
- **The line stops where it is** [checked, `run_b1.py status`, 2026-11-15]:
  - publishing (L1) is done for all nine papers, and the memo draft (Y1) is written;
  - the blind reread (X1) has run and waits at its gate, with 19 equation or value reads
    that disagree;
  - your claim review (U3) has 31 undecided items: those 19 plus the skeptic's 12;
  - the memo fails its own citation check: 39 of 179 sentences uncited
    (`summary_memo_Q1.v1.json`).
- **Nothing further is built.** TODO Phase 3, the v2 rebuild, is dropped.
- **For the report:** present E4 and X1 as one verification step, as the TA suggests.
- **The lesson to state:** six model roles plus a foreman were more orchestration than
  one person's nine-paper review needs. The contract's own build budget was 3 to 5 working
  days; the line took about 17. Doc 20 §4.1 already said "for nine papers the line saves no
  time over reading them yourself."

**Point 2: scope; one pillar with a full frozen evaluation.** This already exists, for the B1
checker against a script-only baseline [checked]:
- **Setup:** the held-out set, never used for tuning, was scored once under the frozen
  prompts b1-v2 on 2026-11-04.
- **Sources:** `runs/literature-review/b1/20261105-pilot2/calibration_heldout.json/calibration.json`
  and `baseline_heldout/baseline.json`.

| On the held-out set | Scripts + checker agent | Script-only baseline |
|---|---|---|
| Quote and page plants, which the script should catch (14) | 14 | 1 |
| Sign and exponent plants (6) | 6 | 5 |
| Plants that reach the checker (31) | 30 (lower bound 85.6 %) | 21 (lower bound 51.5 %) |
| Correct controls flagged (40) | 3 | 14 |

The lower bounds are one-sided 95 % Clopper–Pearson. The plant comparison is clean. The
control row is not held out (Point 3). For 3 of 40, the one-sided 95 % upper bound is
18.3 % (scipy, 2026-11-15), so even taken at face value the data cannot show the rate is
below the 10 % limit.

The other evaluations proposed in Lab Report 2 (Pillar 2's derivation pipeline, Pillar 3's
diagnostician) are not built. A small model-side evaluation is optional for Lab Report 4
(§6).

**Point 3: the 16.7 % false flags without leakage. Half done; the TA is right about the
other half.**
- **The planted errors were split before tuning:** 32 tuning, 45 held out
  (`calibrate.py split`, WORKFLOW.md §5 "Split"). The held-out plants were scored once
  (`attempts_so_far = 1`).
- **The controls were not split** [checked 2026-11-15]. The tuning run `calib-dev8` and the
  held-out run use the same 40 controls, G001–G040; §5 splits only plants.
- **What the tuning did to those controls.** The rounds after Lab Report 2 (whose 5 of 30
  was the last tuning run then):
  - changed the checker prompt because of two controls, G004 and G018;
  - rewrote seven gold statements because of flags on them.
- **So** the held-out "3 of 40 controls flagged" was measured on controls the tuning had
  seen, and is likely optimistic.
- **The fix** is 20–40 fresh controls on pages the tuning never used, scored once with the
  frozen checker: about 1–2 hours of your time and a 40-minute run.
- **Recommendation:** state the limitation in Lab Report 3, and do the fix only if there is
  time before the final paper.

## 3. The model

**What "diffuse domain" means here.**
- **Today:** the solver meshes only the pore space. Regolith grains are the domain
  boundary, and the wall enters as a boundary flux with a contact angle,
  3M·cosθ·h′(φ)/6 (`src/assembly.c`, `Residual_A1`, `FACE_WALL`). Every new grain
  arrangement needs a new body-fitted mesh.
- **The diffuse-domain method** (Li et al. 2009) instead puts the whole box on a plain
  Cartesian mesh and marks the grains with a smooth pore indicator w(x).
  - The equations are multiplied by w.
  - The surface conditions come out of that weighting, or out of a volume term in which
    |∇w| plays the part of the surface.

**Unknowns.** Still 3 per node: ψ (ice fraction of the pore space), T and ρ_v.
- **Given:** w(x) = ½(1 + tanh(d/(2 eps_s))), where d is the signed distance to the grain
  surfaces. It is evaluated analytically at each quadrature point, so there is no fourth
  unknown and the output format is unchanged.
- **Derived:** φ_i = w ψ, φ_a = w (1 − ψ), φ_s = 1 − w.

**Equations, in conservative form** (doc 18 §1, with its revisions in §7–§8):
- **Ice row:** the mass, gradient, double-well and source terms are each multiplied by
  w̃ = w + υ (a floor, `-ddm_floor`, default 1e-6).
  - At 90° that is the whole wall condition: no flux of ψ into the grain (Yu, Chen &
    Thornton 2012, Eq. 20 with cos θ = 0).
  - For other angles add the volume term −|∇w|·3M·cosθ·h′(ψ)/6·N: today's boundary flux,
    with the surface replaced by |∇w| (the step taken by Li et al. 2009 and Yu et al. 2012).
- **Vapor row:** storage and diffusion are weighted by φ_a = w(1 − ψ), so there is no vapor
  flux into grains. The phase-change coupling uses the unfloored w, so water balances
  exactly (doc 18 §7).
- **Heat row:** weighted by w̃. The grains are insulating, which is what today's mesh
  boundary does, so the new model can be compared directly with body-fitted runs. Grain
  conduction is new physics and gets its own option later; it needs regolith properties
  from a source you choose.
- **Reduction check:** with w ≡ 1 the ice row is (1 + υ) times today's row, so the new code
  reduces to the current model.

**Where each piece comes from in the B1 paper set** (to be checked against the page before
it goes in the report; the B1 library holds the checked claims with page numbers, which you
read in review mode):

| Piece | Paper |
|---|---|
| Base sublimation phase-field model | Kaempfer & Plapp 2009 |
| Diffuse-domain method, weighted equations, floor, \|∇w\| surface term | Li, Lowengrub, Rätz & Voigt 2009 |
| Convergence in eps_s (first order expected for no-flux conditions) | Lervåg & Lowengrub 2015 |
| Smoothed-boundary method, contact angle in weighted form (Eq. 20) | Yu, Chen & Thornton 2012 |
| Time scaling of diffusion and source together | Moure & Fu 2024 |
| Why no third phase field is needed | Folch & Plapp 2005; Plapp 2011; Steinbach 2009 |

### Implementation, in order (branch `feat/ddm` from `main`; `-ddm 0` by default, so existing runs are unchanged)

1. **Options:**
   - `-ddm`;
   - grains given either as `-ddm_wall_y` (one flat embedded wall, for tests) or as a grain
     file of circles (x, y, R per line);
   - `-ddm_eps_s` (default eps);
   - `-ddm_floor`.

   `-wall_faces` is refused together with `-ddm 1`.
2. **The indicator:** w and ∇w at each quadrature point, from the grain geometry.
3. **Residual:** a DDM branch of `Residual_A1` for the three rows; 90° first.
4. **Jacobian:** the same scaling. w is fixed in space, so each new term is the old term
   times w or w̃. The θ ≠ 90° wall term adds an h″ term times |∇w|.
   - Checked with `-snes_test_jacobian`.
   - Fallback if this runs late: the finite-difference Jacobian (`-snes_fd_color`), which
     is slower but correct.
5. **Initial condition:** ψ from the existing ice geometries, extended into the grain along
   the wall normal (doc 18 §8). For a flat wall, ψ(x, y_w) is copied below the wall.
6. **The contact-angle term for θ ≠ 90°.** Done today if steps 1–5 pass early; otherwise
   by Lab Report 4.
7. **Output:** the grain geometry is in the staged opts, so post-processing recomputes w and
   masks plots to w ≥ 0.5.

### Tests

Each test is small, 1D or 2D, takes minutes, and is run one at a time.

| # | Test | Reference | Pass |
|---|---|---|---|
| T0 | Regression: `-ddm 0` on a short existing case; then `-ddm 1` with no grains | current `main` output | identical with `-ddm 0`; agrees to solver tolerance with `-ddm 1` |
| T1 | `-snes_test_jacobian` on one grain | finite differences | error ~1e-10 (house standard) |
| T2 | Closed box, one grain, `-xi_v 1`, full coupling | total water ∫ w(ρ_i ψ + (1 − ψ) ρ_v) is conserved | drift at solver tolerance |
| T3 | Flat embedded wall, ice–air interface meeting it, 90° | 90° contact | within 2° as eps_s and eps vary |
| T4 | 1D ice front in a pore whose width varies slowly (w linear in x) | V = −3M eps d(ln w)/dx (doc 18 §2; 0.998 of it in the toy solve) | within a few % |
| T5 | θ = 60° on the flat embedded wall (after step 6) | body-fitted θ∞ = 60.11° on the 142×48 channel case (`studies/contact_angle/verification/dtmax_study.csv`) | within 1–2° |
| Demo | Ice bridge between two embedded circular grains, sublimating | none; qualitative | runs a set physical time; ParaView snapshot saved to `review_scratch` |

- **The minimum for "working model" today:** T0, T1, T2 and T3.
- **Then, as time allows:** T4, T5 and the demo.
- **For Lab Report 4:** the curved-grain orthogonality metric (|n_ψ · n_w| < 0.03, doc 18
  §7), eps_s halving for the convergence order, and a DDM run against a body-fitted run on
  the existing wedge geometry.

## 4. Timeline

| When | What | Who |
|---|---|---|
| Sun 15 Nov, afternoon | Decide §7; branch; steps 1–3 | you decide; I build |
| Sun evening | Steps 4–5; build; T1, then T0, T2, T3, one at a time | I run (needs your yes, D5) |
| Mon 16 Nov, day | T4, demo, T5 if step 6 is done; fresh-context review of the code against doc 18; fixes | I run; you look at the snapshots |
| Mon evening | Lab Report 3 submitted (one day late) | you |
| to Sun 29 Nov (Lab Report 4) | analytic Jacobian if the FD one was used; curved grain; eps_s convergence; DDM vs body-fitted; failure modes stated; optional model-side agent evaluation (§6) | together |
| Wed 2 Dec | presentation | you |
| 4 or 9 Dec | 4-page IEEE paper (the syllabus says 4 Dec, the Project Guideline 9 Dec; ask the TA) | you, with the writer helper |

## 5. Lab Report 3

**What is required** [checked by the reviewer: `Lab-Guideline.docx`, "Lab Report 3" and
"What Every Lab Report Includes"; `Lab_Report_Template_Agent_Project.docx`, sections A–E]:
- the first working piece and a genuine demo, where "the smallest slice of your system that
  can actually run end to end" counts;
- the repository link with an up-to-date README;
- the progress report (done / remaining / bottlenecks) and next steps;
- about 5 pages, with an AI Use Declaration.

**Proposed content, in the template's order:**
1. **The pivot and why:** the TA's three points, your own concern, and the numbers in §2,
   including the shared-controls limitation.
2. **Done:** the B1 line and its frozen held-out comparison (the table in §2); the pilot run
   (nine papers, 230 published claims, a memo draft that fails its citation check).
3. **The working piece (C3, the demo):** the diffuse-domain model, its equations with the
   paper behind each term, the T0–T3 results as numbers, and snapshots or a short recording
   of the demo run.
4. **Remaining and bottlenecks:** tokens and orchestration in B1 (the lesson); the DDM risks
   in §6.
5. **Next two weeks:** the Lab Report 4 row of §4.

## 6. The agents in this plan, and the risks

**Agents, kept small.** The Project Guideline's minimum, as quoted by the reviewer from
`Project-Guideline.docx`, "Minimum Technical Expectations":
- an LLM agent with a domain tool;
- human approval points;
- deterministic checks plus a non-agent baseline;
- reliability across repeated runs.

It also says a smaller system evaluated carefully is stronger than a broad one. In this plan:
- **Builder and runner:** this Claude Code session. Its domain tool is the solver, through
  `run_lunar.sh`. You approve at the branch points (§7), and scripts decide pass or fail.
- **One reviewer agent** with fresh context, which audits the code against doc 18 before the
  tests count.
- **Optional for Lab Report 4: a small model-side evaluation.**
  - Plant about 12 bugs in the DDM code: a sign error in the wall term, a missing w in the
    mass term, a missing |∇w|, a wrong h″ in the Jacobian, eps_s below the mesh resolution,
    and so on.
  - Ask a single headless diagnostician call to name each bug's cause from the test outputs,
    3 times per bug.
  - Compare against a fixed rule table, for example "Jacobian test fails → code".
  - This adds the repeated-runs reliability the course asks for, on the model rather than
    on B1. About 36 calls.

**Risks:**
- **Time** [inferred]. The Jacobian and the initial-condition extension are the slow parts.
  The finite-difference fallback and the 90°-first order keep "working today" reachable.
- **The grain interior** [inferred]. The floor υ scales rows inside grains by about 1e-6.
  That may slow Newton and the linear solver, or upset the per-field convergence tests, so
  watch the step sizes. ψ inside the grain can disturb the contact line unless it is extended
  properly.
- **Resolution** [checked, doc 18 §7 item 8]. Wall accuracy is first order in eps_s.
  - The grain band needs 3 to 5 elements, so meshes are larger than body-fitted ones.
- **Measurement** [inferred]. The contact-angle script assumes the wall is the domain edge;
  it needs a w-mask. An error there would look like a physics failure.
- **Not validated, even today** [checked, `memory/MEMORY.md`, solver section]. Only the
  static contact angle is validated in the body-fitted model.

## 7. Decisions for you

- **D1. Stop B1 now** (recommended).
  - The 31 review items and the X1 gate stay open, and the memo is reported as an
    unreviewed draft.
  - Optional: check by hand only the 10–15 claims the model rests on (the equations of Li
    2009 and Yu 2012), in review mode.
- **D2. Heat:** insulating grains first, matching today's mesh boundary (recommended).
  Grain conduction later, behind its own option, with properties from a source you choose.
- **D3. Contact angle:** 90° first, then the |∇w| term for any angle (recommended). The
  alternative is 90° only this month.
- **D4. Demo geometry:** an ice bridge between two grains (recommended), or three grains.
- **D5. Running the tests:** may I run T0–T5 myself, one after another (each a single small
  local run), or should I hand you the commands? The project limit is one run without your
  specific instruction.
- **D6. The controls:** state the limitation only (recommended for now), or write fresh
  held-out controls later.

## 8. The independent review

Started with fresh context, reading only the course documents, Lab Report 2, the planning
documents, the B1 records and the solver. It was not told the main session's view.

**Where both views agree:**
- B1 is over-built for the course;
- report the frozen checker-versus-baseline result as the completed evaluation;
- the controls were not held out. Both views found this independently: the reviewer from
  `SPEC.md` §6.3 and TODO 1.3–1.5, the main session from the run files;
- the model to build is doc 18's F2;
- the minimum for today is the 90° case with a regression test, the Jacobian test and one
  wall test.

**What the reviewer added, now in the plan:**
- insulating grains, which removes the need for regolith properties today and makes a
  direct body-fitted comparison possible;
- the |∇w| volume term for any contact angle, with the existing body-fitted 60.11° result as
  its reference;
- the planted-bug diagnostician as an optional, bounded model-side evaluation;
- the course minimums from the Project Guideline;
- the 4 vs 9 Dec deadline conflict.

**What the reviewer flagged for you to check directly:**
- the sign of the |∇w| wall term against `src/assembly.c` lines 81–100, before it is coded;
- the deadline, with the TA.
