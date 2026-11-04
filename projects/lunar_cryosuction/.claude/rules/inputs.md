---
paths:
  - "inputs/**"
  - "preprocess/**"
  - "postprocess/**"
  - "studies/**"
---

# Input files and quarantine folders

Live inputs are `inputs/solver.opts` plus `inputs/geometry/<family>/` and
`inputs/experiment/<family>/`; older ones are in `inputs/scratch/`. Before
writing a new input file, search there (`find inputs/scratch -name
'*<name>*'`); if it exists, `git mv` it back to `inputs/geometry/<family>/` or
`inputs/experiment/<family>/`. The same holds for `preprocess/scratch/`,
`postprocess/scratch/` and `studies/scratch/`: everything in them is pending
deletion unless pulled back. The root `scratch/` folder is git-ignored local
storage (the archived practice round, a stashed run, and the HPC output
fallback when `$SCRATCH` is unset). Never use it as a source for current work.
