# 🧊 Effective Thermal Conductivity Simulation

This project simulates the effective thermal conductivity of a porous ice medium
using a phase field approach and isogeometric analysis. It utilizes **PETSc** 
and **PetIGA** for numerical solving and is intended for performance and 
extensibility.

---

## 📁 Folder Structure

```
.
├── include/        # Header files
├── inputs/         # Initial condition and parameter input files
├── obj/            # Compiled object files (auto-generated)
├── postprocess/    # Post-processing scripts and analysis tools
├── scripts/        # Run scripts and automation tools
├── src/            # Source code for simulation
├── Makefile        # Build instructions
└── README.md       # You're here!
```

---

## 🛠️ Build Instructions

1. Ensure PETSc and PetIGA are installed and environment variables are set:
   ```bash
   export PETSC_DIR=/path/to/petsc
   export PETSC_ARCH=arch-name
   export PETIGA_DIR=/path/to/petiga
   ```

2. Compile the code:
   - Optimized mode:
     ```bash
     make all
     ```
   - Debug mode:
     ```bash
     make debug
     ```

3. Clean build artifacts:
   ```bash
   make clean
   ```

---

## 🚀 Running the Simulation

Use the provided run scripts in `scripts/`:

- Optimized run:
  ```bash
  ./scripts/run_effective_k_ice.sh
  ```

- Debug run (includes LLDB and PETSc diagnostics):
  ```bash
  ./scripts/run_effective_k_ice_DEBUG.sh
  ```

Simulation results are saved in timestamped folders under:

```
~/SimulationResults/ThermalConductivity/
```

---

## 📊 Post-Processing

Post-processing tools and example notebooks can be found in the `postprocess/` 
directory. You may include:

- VTK field visualization
- Binary output parsers
- Plotting thermal gradient results
- Python or MATLAB scripts for analysis

---

## 🧪 Testing (Optional)

Add minimal test cases to the `inputs/` directory. For example:

- Constant initial temperature field
- Analytical solutions for thermal diffusion

These help validate numerical accuracy and performance.

---

## 📌 Development Notes

- Code uses PETSc vectors for temperature and phase fields
- Supports 2D and 3D domains
- Interface width (`eps`) computed automatically
- Run configuration controlled via environment variables in the scripts

---

## ✅ TODO

- [ ] Add validation test case
- [ ] Automate post-processing summary
- [ ] Allow for input files