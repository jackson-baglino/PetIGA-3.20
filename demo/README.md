# README: Phase Field Simulation Code with PETSc and PetIGA

## Overview
This code implements a phase field model for simulating ice, vapor, and thermal interactions using PETSc and PetIGA. It utilizes the finite element method for solving the governing equations, incorporating adaptive time stepping, nonlinear solvers, and IGA (Isogeometric Analysis) formulations.

## Features
- Solves phase field equations with temperature and vapor density fields.
- Includes thermal conductivity, heat capacity, and density computations.
- Implements phase transition kinetics based on Gibbs-Thomson relationships.
- Uses adaptive time stepping to ensure numerical stability.
- Supports initialization from predefined conditions or input files.
- Outputs simulation results at specified time intervals.
- Parallelized using MPI for efficient large-scale computations.

## Dependencies
To compile and run this code, ensure you have the following installed:
- [PETSc](https://petsc.org/release/)
- [PetIGA](https://bitbucket.org/dalcinl/petiga)
- MPI (Message Passing Interface)
- C compiler (e.g., `mpicc`)

## Compilation
To compile the code, use:
```sh
make all
```
Or compile manually with:
```sh
mpicc -o phase_field_simulation phase_field_simulation.c -I$(PETSC_DIR)/include -L$(PETSC_DIR)/lib -lpetsc -I$(PETIGA_DIR)/include -L$(PETIGA_DIR)/lib -lpetiga
```

## Running the Simulation
Run the compiled executable with:
```sh
mpirun -np <num_procs> ./phase_field_simulation
```
where `<num_procs>` is the number of MPI processes.

### Environment Variables
The code reads several environment variables for configuration:
- `Nx, Ny, Nz`: Number of grid points in x, y, and z directions.
- `Lx, Ly, Lz`: Domain size in x, y, and z directions.
- `delt_t`: Initial time step.
- `t_final`: Final simulation time.
- `n_out`: Number of output steps.
- `humidity`: Initial humidity.
- `temp`: Initial temperature.
- `grad_temp0X, grad_temp0Y, grad_temp0Z`: Initial temperature gradient.
- `dim`: Simulation dimensionality (2D or 3D).
- `eps`: Interface width parameter.

Example usage:
```sh
export Nx=100 Ny=100 Nz=1 Lx=0.01 Ly=0.01 Lz=0.001 delt_t=1e-4 t_final=1.0 n_out=10 humidity=0.8 temp=-10 grad_temp0X=0 grad_temp0Y=0 grad_temp0Z=0 dim=2 eps=1e-6
mpirun -np 4 ./phase_field_simulation
```

## Input and Output Files
- **Initial Conditions:**
  - The simulation can start from an input file or generate initial conditions internally.
  - If reading from a file, set `readFlag=1` and provide the path via `inputFile`.

- **Output Data:**
  - Time evolution of ice phase, temperature, and vapor density is stored in `sol_#####.dat` files.
  - Global statistics (e.g., total ice, temperature) are saved in `SSA_evo.dat`.

## Main Functions and Modules
### `SNESDOFConvergence`
Checks the convergence of the nonlinear solver (SNES) at each iteration and prints residual norms.

### `ThermalCond`, `HeatCap`, `Density`
Computes the thermal conductivity, heat capacity, and density of the mixture based on phase fractions.

### `Residual`
Defines the system of equations to be solved by PETSc.

### `Jacobian`
Computes the Jacobian matrix for the nonlinear solver.

### `Integration`
Computes domain integrals for output quantities such as total ice volume.

### `Monitor`
Tracks the simulation progress and outputs relevant statistics.

### `OutputMonitor`
Handles file output for post-processing.

### `InitialIceGrains`, `InitialSedGrains`
Generates or reads initial conditions for ice grains and sediments.

### `FormInitialCondition2D`, `FormInitialCondition3D`
Initializes the phase field, temperature, and vapor density fields based on input parameters.

## Adaptive Time Stepping
The code uses an adaptive time stepping strategy:
- Adjusts time step size based on nonlinear solver iterations.
- Ensures stability by enforcing minimum (`dtmin`) and maximum (`dtmax`) time step sizes.
- Uses a scaling factor (`factor=10^(1/8)`) to increase/decrease step size.

## Parallel Execution
The code is designed for distributed computing with MPI:
- Each MPI process handles a portion of the computational domain.
- Communication is managed through PETSc’s parallel data structures.
- Output is collected and synchronized across processes.

## Notes and Best Practices
- Ensure PETSc and PetIGA are compiled with MPI support.
- Choose appropriate grid resolution (`Nx, Ny, Nz`) based on available computational resources.
- Monitor output logs to detect numerical instability.
- Adjust `eps` to control interface sharpness (smaller values lead to thinner interfaces).
- Use `mpirun` with more processes for large 3D simulations.

## Troubleshooting
### Segmentation Faults or PETSc Errors
- Ensure all environment variables are set correctly.
- Check the correctness of input files.
- Try running with fewer MPI processes (`mpirun -np 1`).

### Slow Performance
- Increase `Nx, Ny, Nz` gradually to find an optimal balance between accuracy and performance.
- Use `-log_view` to analyze PETSc performance bottlenecks.

## Acknowledgments
This code is built using PETSc and PetIGA for high-performance finite element simulations. Thanks to the PETSc and PetIGA development teams for providing powerful computational tools.

---
For questions or contributions, please contact the project maintainers.

