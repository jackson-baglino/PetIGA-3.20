# Bound-Constrained Newton Solve for the Ice Phase Field

> **Purpose:** Reference document explaining how the ice phase-field DOF is
> kept within physical bounds during the nonlinear solve.
> **Code:** `src/permafrost2.c` (bound setup), `src/snes_convergence.c`
> (custom per-DOF convergence test), `inputs/solver.opts` (SNES/KSP options).
> **Last updated:** 2026-06-23

---

## 1. The problem

The ice phase field φ_i is only physically meaningful on [0, 1] (0 = pure
air/vapor, 1 = pure ice). The discretized Allen–Cahn evolution equation for
φ_i is solved implicitly each time step with Newton's method inside PETSc's
`TS`/`SNES` machinery. An *unconstrained* Newton step has no notion of this
bound: a large trial step (particularly at the larger end of the adaptive
`dtmax` range, or near a sharp ice–sediment interface where the B-spline
solution can overshoot) can push individual DOFs outside [0, 1], producing
unphysical states that then feed into the next residual evaluation and can
destabilize the solve.

Two complementary mechanisms guard against this:

1. **Function-domain rejection** (`SNESSetFunctionDomainError`) — catches a
   bad trial iterate *during* the line search, before it is accepted.
2. **Hard variable bounds on the DOF vector** via a variational-inequality
   (VI) SNES — the subject of this document.

---

## 2. Variational-inequality (VI) bound constraints

PETSc's `SNESVISetVariableBounds(snes, Xl, Xu)` turns the plain Newton solve
into a **bound-constrained nonlinear complementarity problem**: instead of
solving `F(x) = 0` directly, it solves the equivalent VI

$$x_l \le x \le x_u, \qquad (x - x_l)^T F(x) \ge 0, \qquad (x_u - x)^T F(x) \ge 0$$

so that any DOF sitting exactly at its lower (upper) bound is only allowed to
have a non-negative (non-positive) residual component — i.e. the solver
treats a DOF clamped at a bound as a complementarity condition, not as a
violated equation. This is solved with PETSc's semismooth Newton method
(`-snes_type vinewtonssls` in `inputs/solver.opts`), which is the practical
mechanism that keeps φ_i inside its bounds at **every** Newton iteration of
**every** time step, not just at convergence.

### Why bounding the DOFs bounds the whole field

The bound is set on the *control-point* (DOF) vector, not on φ_i pointwise.
This is sufficient because the field value at any point in the domain is

$$\phi_i(\mathbf{x}) = \sum_k N_k(\mathbf{x})\, c_k$$

a **convex combination** of nearby control-point values c_k (B-spline basis
functions are non-negative and sum to 1 — the partition-of-unity property).
A convex combination of numbers in `[a, b]` is itself in `[a, b]`, so bounding
every c_k to `[a, b]` automatically bounds φ_i(**x**) to `[a, b]` everywhere,
including between control points. No extra pointwise constraint is needed.

### Where it's set in code

`src/permafrost2.c`, right after the custom SNES convergence test is wired up:

```c
Vec Xl, Xu;
ierr = IGACreateVec(iga, &Xl); CHKERRQ(ierr);
ierr = IGACreateVec(iga, &Xu); CHKERRQ(ierr);
ierr = VecStrideSet(Xl, 0, -0.01);           CHKERRQ(ierr);   // ice lower bound
ierr = VecStrideSet(Xu, 0,  1.01);           CHKERRQ(ierr);   // ice upper bound
ierr = VecStrideSet(Xl, 1, PETSC_NINFINITY); CHKERRQ(ierr);   // temperature: unbounded
ierr = VecStrideSet(Xu, 1, PETSC_INFINITY);  CHKERRQ(ierr);
ierr = VecStrideSet(Xl, 2, PETSC_NINFINITY); CHKERRQ(ierr);   // vapor density: unbounded
ierr = VecStrideSet(Xu, 2, PETSC_INFINITY);  CHKERRQ(ierr);
ierr = VecStrideSet(Xl, 3, PETSC_NINFINITY); CHKERRQ(ierr);   // sediment phase: unbounded
ierr = VecStrideSet(Xu, 3, PETSC_INFINITY);  CHKERRQ(ierr);
ierr = SNESVISetVariableBounds(nonlin, Xl, Xu); CHKERRQ(ierr);
```

Only DOF 0 (φ_i) carries a finite bound. Temperature and vapor density are
genuinely unbounded physical fields. The sediment phase φ_s is currently left
unbounded too (it is frozen after `t_sed_freeze`, so it never needs the VI
machinery in practice).

---

## 3. Why the bound isn't exactly [0, 1]

The mathematically correct bound is strict `[0, 1]`. In production the bound
is currently set slightly wider, **`[-0.01, 1.01]`**, as an intentional,
documented tradeoff — not an oversight. The reason is a separate, known issue
in the custom convergence test:

### The underlying issue: trivial per-DOF convergence

`src/snes_convergence.c` (`SNESDOFConvergence`) implements a *per-DOF*
convergence test — each of φ_i, T, ρ_v is checked independently against
relative (`rtol`), absolute (`atol`), and step (`stol`) criteria, rather than
using PETSc's default single combined residual norm. This gives much finer
control over which physics is allowed to call a step "converged," but it
interacts badly with strict `[0, 1]` bounds:

- When φ_i is exactly at a VI-active bound (0 or 1), its residual there is
  governed by the complementarity condition, not the PDE residual — so the
  per-DOF residual norm for that DOF block trivially satisfies `atol` almost
  immediately.
- Confirmed by direct A/B testing: `atol = 1e-6` and `atol = 1e-8` give
  bit-for-bit identical results (the criterion that's actually firing is the
  trivial one, not a meaningfully tight tolerance), while `atol = 1e-20`
  makes Newton stagnate and never converge — the true residual floor for some
  DOF blocks sits below float64-meaningful precision.
- Net effect: with a strict `[0, 1]` bound, nearly every step converges in a
  single Newton iteration regardless of step size, and small ice-cap features
  can **pulse** (shrink for several steps, jump back up one step, plateau at
  a small nonzero residual) instead of evolving smoothly toward extinction —
  a DOF pinned exactly at the bound gets "stuck" by the complementarity
  treatment rather than continuing to move the way the (noisy) single-step
  dynamics actually want.

This is a genuine per-DOF-tolerance-design problem (each field's natural
residual scale needs its own calibrated `atol`, not a single shared value) —
not something fixable by changing the bound alone. It is tracked as a
post-conference follow-up.

### The practical workaround

Loosening the hard VI bound slightly beyond `[0, 1]` gives the Newton step
just enough room to move before the trivial-convergence/pinning behavior
kicks in, which is what actually lets coarsening (grain shrinkage, Ostwald
ripening between caps) show up visibly instead of stalling. The exact amount
of slack needed is empirical and depends on `-dtmax` (smaller `dtmax` →
gentler steps → less slack needed):

| Bound | dtmax | Notes |
|---|---|---|
| `[0, 1]` | 1e5 | original, physically correct; pulsing/stalling observed |
| `[-0.05, 1.05]` | 1e5 | adds slack; still some pulsing at this dtmax |
| `[-0.1, 1.1]` | 1e5 | more slack; even 1052/1052 steps still trivially converge in 1 iteration |
| `[-0.05, 1.05]` | 1e4 | smaller dtmax alone reduces the need for slack |
| `[0, 1]` | 1e4 | strict bound retested at the smaller dtmax |
| **`[-0.01, 1.01]`** | **1e4** | **current production setting** — small slack, smaller dtmax |

In the most recent production run (2-week, `dtmax = 1e4`, bound
`[-0.01, 1.01]`), the actual excursion stayed within **±0.003** of `[0, 1]` —
well inside the allowed ±0.01 — with 0 diverged steps, confirming the smaller
`dtmax` does most of the stabilizing work and only a small amount of bound
slack is needed on top of it.

This is technically incorrect in the sense that it allows φ_i a small
unphysical excursion outside `[0, 1]`; it is an intentional, bounded tradeoff
to get visibly-evolving, presentable results. The real fix — per-DOF-block
`atol` matched to each field's natural residual scale — is deferred.

---

## 4. Summary for a slide

- Ice phase field is bound-constrained directly on its control points via a
  **variational-inequality semismooth Newton solve**
  (`-snes_type vinewtonssls`), not by post-hoc clipping.
- Bounding the DOFs is sufficient to bound the field everywhere, because
  B-spline basis functions form a convex partition of unity.
- The production bound is `[-0.01, 1.01]` rather than the physically exact
  `[0, 1]` — a small, deliberate, documented amount of slack that works
  around a known per-DOF convergence-tolerance issue, validated empirically
  to keep the true excursion almost an order of magnitude inside the allowed
  range (±0.003 actual vs. ±0.01 allowed).
