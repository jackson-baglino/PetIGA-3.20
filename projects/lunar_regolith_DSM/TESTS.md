# Testing Guide — Permafrost Phase-Field Simulation

This document describes how to systematically test the permafrost code. Because the
code solves a nonlinear coupled PDE system with no exact closed-form solution, testing
is organized around three complementary strategies:

1. **Smoke tests** — does the code run at all?
2. **Regression tests** — do scalar outputs stay stable across code changes?
3. **Physics verification tests** — do the results match known analytic or semi-analytic
   behaviors (conservation laws, limiting cases, published results)?

---

## Directory Layout

```
inputs/tests/
├── test1_IceCap.opts            # 2D layered ice cap
├── test2_CapillaryBridge.opts   # 2D capillary neck between two grains
├── test3_EnclosedGrainPair.opts # 2D enclosed grain pair
├── test_1D_IceSlab.opts         # 1D centered ice slab (new)
└── test_1D_FlatInterface.opts   # 1D flat ice-air interface (new)
```

---

## How to Run a Test

```bash
# Serial
mpirun -np 1 ./permafrost -options_file inputs/tests/<test>.opts \
    -output_path outputs/<test_name>

# Parallel (e.g., 4 ranks)
mpirun -np 4 ./permafrost -options_file inputs/tests/<test>.opts \
    -output_path outputs/<test_name>
```

Set the `folder` environment variable for the SSA output file:
```bash
export folder=outputs/<test_name>
mkdir -p $folder
mpirun -np 1 ./permafrost -options_file inputs/tests/<test>.opts \
    -output_path $folder
```

---

## Test Suite

### T1 — Smoke Test (1D Slab)

**File:** `inputs/tests/test_1D_IceSlab.opts`  
**Purpose:** Verify that the code initializes, runs one time step, and produces
output in 1D without crashing.  
**Pass criteria:**
- Program exits with code 0
- At least one `sol_*.dat` file is written
- `TOT_ICE` in the monitor output is approximately `0.3 * Lx` (30% of domain)
- No NaN or Inf in any field (check monitor output)

---

### T2 — Smoke Test (1D Flat Interface)

**File:** `inputs/tests/test_1D_FlatInterface.opts`  
**Purpose:** Verify that a 1D ice-air interface with a temperature gradient runs
and that temperature BCs are applied correctly.  
**Pass criteria:**
- Program exits with code 0
- Temperature at boundaries matches `temp0 ± grad_temp0 * Lx/2`
- No NaN or Inf in any field

---

### T3 — Mass Conservation (Phase Field)

**Setup:** Run `test_1D_IceSlab.opts` with a short time horizon (`t_final = 100`).  
**Check:** In the monitor output:
```
TOT_ICE + TOT_AIR ≈ Lx   (for sediment-free 1D domain)
```
The sum of phase fractions must be conserved to within solver tolerance.  
**Why it matters:** The Allen-Cahn equation is not conservative by itself; the
Lagrange multiplier constraint (`phi_i_t + phi_a_t = 0`) enforces it. If this
check fails, the constraint is not working.

---

### T4 — Symmetry Test (1D Slab, No Gradient)

**Setup:** Run `test_1D_IceSlab.opts` with `grad_temp0 = 0` and `humidity = 1.0`
(no driving force).  
**Check:** The solution at step N should be symmetric about `x = 0.5 Lx`:
```
phi_ice(x) ≈ phi_ice(Lx - x)   for all x
```
Extract from `sol_00010.dat` and check numerically.  
**Why it matters:** Verifies that no asymmetric artifacts are introduced by the
parallelization or boundary conditions.

---

### T5 — Convergence Test (Mesh Refinement, 1D)

**Setup:** Run `test_1D_IceSlab.opts` with increasing `Nx`: 64, 128, 256.  
**Check:** The interface width of the phase field at `t=0` (from the initial
condition) should converge at rate p+1 (expected 3rd order for p=2).  
Measure the maximum gradient `max |d phi/dx|` from the initial solution files
and verify it scales as `O(h^{p+1})`.  
**Why it matters:** Verifies that the IGA B-spline basis is correctly resolving
the diffuse interface and that the IC function is smooth.

---

### T6 — 2D Capillary Bridge Regression

**File:** `inputs/tests/test2_CapillaryBridge.opts`  
**Purpose:** Regression test for the 2D capillary neck geometry.  
**Check:** Record `TOT_ICE`, `I-A INTERF`, and `TRIPL_JUNC` at the final step.
After any code change, re-run and verify these values have not changed by more
than 0.5%.  
**Reference values:** (fill in after a known-good run)
```
Step  TOT_ICE    I-A INTERF   TRIPL_JUNC
  0   <value>    <value>      <value>
100   <value>    <value>      <value>
```

---

### T7 — Parallel Reproducibility

**Setup:** Run `test_1D_IceSlab.opts` with 1, 2, and 4 MPI ranks.  
**Check:** The monitor output `TOT_ICE`, `TOT_AIR`, `TEMP`, `TOT_RHOV` should be
bit-identical (or within floating-point rounding) across all three runs.  
**Why it matters:** Verifies that `MPI_Bcast` and `MPI_Allreduce` calls are
correct and that domain decomposition does not affect results.

---

### T8 — Temperature BC Test (1D Flat Interface)

**File:** `inputs/tests/test_1D_FlatInterface.opts`  
**Purpose:** Verify that fixed-temperature boundary conditions are applied
correctly.  
**Check:** At `t=0`, the temperature field should be:
```
T(x) = temp0 + grad_temp0 * (x - 0.5 * Lx)
T(0) = temp0 - grad_temp0 * Lx / 2
T(Lx) = temp0 + grad_temp0 * Lx / 2
```
These BCs should be maintained throughout the simulation.

---

### T9 — Vapor Saturation at Equilibrium

**Setup:** Run `test_1D_IceSlab.opts` with `humidity = 1.0` and `t_final = 0`
(initialization only, no time stepping).  
**Check:** At every node, `rhov = hum0 * rhoVS(T)`. Extract from the initial
solution file and verify the ratio `rhov / rhoVS(T)` is within 0.1% of `hum0`.  
**Why it matters:** Verifies that `RhoVS_I` and the IC function are consistent.

---

### T10 — SNES Convergence Rate

**Setup:** Add `-snes_monitor` to any test run.  
**Check:** The SNES residual should decrease at least geometrically each Newton
iteration (reduction factor < 0.5 per step). Quadratic convergence is expected
once the iterate is close to the solution.  
**Why it matters:** If convergence is slow or non-monotone, the Jacobian
approximation (currently FD) may be inaccurate or the time step may be too large.

---

## Checking Output Files

Solution files are written as `sol_NNNNN.dat` in the output directory. To inspect
them, use the provided Python plotting scripts:

```bash
python scripts/plotpermafrost.py --dir outputs/<test_name> --step 10
```

To extract scalar time-series data:
```bash
cat outputs/<test_name>/SSA_evo.dat
# Columns: sub_interf/eps  tot_ice  t  step
```

---

## Known Limitations / Not Yet Tested

- No automated pass/fail infrastructure (manual inspection required for now)
- Exact reference values for regression tests (T6) not yet recorded
- 3D test cases not yet created
- Adjoint / sensitivity tests not implemented
- The analytic Jacobian (`Jacobian()` in `assembly.c`) is not registered; FD
  Jacobian is used. A test comparing FD vs. analytic Jacobian would be valuable
  once the analytic version is complete.
