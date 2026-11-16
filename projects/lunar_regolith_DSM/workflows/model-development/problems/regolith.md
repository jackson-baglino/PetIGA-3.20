# Problem brief: ice, water vapor and heat in a regolith pore space

This is the brief for M1:Scope of the model line. It states the problem, what is already decided by the
user, and the constraints. It does not say which model to build: choosing the model type is M2:Selection.

## The physical system

Water ice in the pore space of lunar (or Martian) regolith: grains of rock that do not move, ice that
sublimates and deposits, water vapor that diffuses through the open pore space, and heat that conducts
through ice, grains and pore space, with latent heat released or absorbed at the ice surface. Processes of
interest: sublimation and deposition, cryosuction (ice moving toward cold regions under a temperature
gradient), ice aggregation and the evolution of the ice shape and its contact with the grains. Ice–ice
sintering by itself is not the focus. Adhesion physics at the ice–grain contact is out of scope; the
contact angle is a parameter to be swept.

## What exists today

The solver (PETSc 3.20 / PetIGA, C, generalized-alpha in time) solves a two-phase model with three
unknowns per node: the ice phase field φ ∈ [0, 1], the temperature T and the vapor density ρ_v. The pore
space is the mesh: grains are the domain boundary, and the wall condition is a boundary flux with a contact
angle. Read:
- `src/assembly.c`: `Residual_A1` and `Jacobian_A1` (the equations as implemented, including the wall
  term, which is validated for the static contact angle and is to be kept as the reference);
- `src/material_properties.c`: the constitutive laws (saturation vapor density, vapor diffusivity, thermal
  conductivity, heat capacity);
- `src/lunar_main.c`: parameters, the kinetic coefficients and the time-scale estimate tau_sub;
- `include/NASA_types.h`: the parameter set (`AppCtx`);
- `inputs/solver.opts`: the solver settings.
The model it follows is Kaempfer & Plapp 2009 (Phys. Rev. E 79, 031502) with the time scaling of Moure &
Fu 2024 (Cryst. Growth Des., doi 10.1021/acs.cgd.4c00539); both are in the library.

Known problems with today's approach:
- every new grain arrangement needs a new body-fitted mesh;
- an earlier three-phase phase-field attempt (grains as a third field) grew a spurious air film between ice
  and grain, where the ice field relaxed against the stationary grain field;
- the constitutive laws are for Earth air (saturation vapor density through the dry-air density and
  atmospheric pressure; a binary vapor-in-air diffusivity), which do not hold in vacuum;
- in SI units the residual blocks are about 1e-13 (ice), 1e-14 (T) and 1e-24 (vapor), so the absolute
  Newton tolerance passes at the first iteration and the per-field convergence test means little;
- the time step is controlled by Newton iteration counts, a hand-set maximum step and an after-the-fact
  limit on the change of φ per step; no time scale is used a priori;
- the thin-interface correction is on by default, but the option's own help text says it inflates the
  kinetic coefficient for one-sided vapor diffusion (`src/lunar_main.c`, option `-thin_iface_corr`).

## Decided by the user (2026-11-16)

- **Notation**, if the model uses a fixed grain indicator: ψ is the indicator (1 in the pore, 0 in the
  grain), φ ∈ [0, 1] the ice phase field in the pore, φ_i = ψφ, φ_a = ψ(1 − φ), φ_s = 1 − ψ, ψ̃ = ψ + υ a
  floored indicator, eps the ice interface parameter, eps_s the grain-edge width. Also σ_s for a sediment
  well coefficient and γ_ia, γ_is, γ_as for the interfacial energies.
- **The contact-angle (wetting) term** in today's model is known and correctly implemented; a new model
  starts from it.
- **Mixed thermal conductivity** uses the arithmetic–harmonic (tensorial) blend of Nicoli, Plapp & Henry
  2011 (Phys. Rev. E 84, 046707): arithmetic along an interface, harmonic across it.
- **Regolith k and c_p**: from the literature for the lunar regolith composition, or from a proxy material
  with well-known properties; the user chooses at G:Scope.
- **Comparison test**: where no sediment is present, the new model must agree with today's phase-field
  model; one test compares the neck width of two sintering ice grains over time.

## Constraints

- The implementation is in this solver (PetIGA, C1 splines on structured patches), behind a runtime
  option, default off, so existing runs are unchanged. Another code is a possible answer in M2 only if the
  case for it is strong and stated.
- Tests run on this machine: at most 8 ranks, 20k–100k DOFs per rank (3 DOFs per node today), at most 30
  minutes per run.
- Performance is part of the model: nondimensional form, a time-step rule that uses a priori time scales,
  and solver and preconditioner choices, all justified by measurement.
- Deadlines: a working model with its tests by early December 2026.

## What the line must produce

The scope and scaling (M1), the choice of model type (M2), the derivation (M3), the discretization spec
(M4), the locked tests (M5), the implementation (M6), the mutation check (M7), the performance study (M8)
and a teaching document (M9). See `workflows/model-development/WORKFLOW.md`.
