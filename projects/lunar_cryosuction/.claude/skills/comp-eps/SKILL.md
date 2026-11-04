---
name: comp-eps
description: Compute the phase-field interface parameter eps and the mesh (Nx, Ny, Nz) for given domain/grain sizes, temperature and attachment coefficient, using preprocess/comp_eps.py (Kaempfer and Plapp 2009 bounds). Use whenever eps, mesh size or interface resolution is chosen or checked.
allowed-tools: Bash(venv_lunar/bin/python preprocess/comp_eps.py *)
---

# comp-eps

eps is a decay-length scale, not the visible band width. Never estimate or
rescale it by eye; always run the script.

## Run

From the project root:

```bash
venv_lunar/bin/python preprocess/comp_eps.py --Lx <m> [--Ly <m>] [--Lz <m>] --T0 <degC> --alpha <alpha_c>
```

Add `--patch inputs/geometry/<family>/<name>.opts` only when the user asks
for the geometry file to be updated; it writes `-eps` and `-Nx/-Ny/-Nz`.

## Report

Copy these from the script output verbatim:

1. `eps` and `h`, and which bound is binding.
2. `Nx, Ny[, Nz]`, total elements and DOFs.
3. **Elements across the 5–95 % band.** This is the resolution figure. The
   1–99 % count may be given after it, labeled as such.
4. Every `[NOTE]` and any check that is not `[OK]`.
5. The `solver options` line.

Give the exact command you ran. If the script fails or a check fails, report
that and stop; do not substitute numbers of your own.
