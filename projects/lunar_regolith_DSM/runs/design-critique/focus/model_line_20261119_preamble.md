# The model line of a thesis project: design, plan and record, bundled for review

This file bundles many files of one research repository, because you may read only this file. The project
is paper 3 of a PhD thesis: a phase-field model (C, PETSc/PetIGA) of ice, water vapor and heat in the pore
space of lunar regolith. The model line is an agentic system meant to derive, code and test a new model
(a diffuse-domain model), first by rederiving the current model.

Each part below starts with a header giving its path in the repository, the git ref it was taken from and
its sha256. A part is a copy of the file as committed; nothing in it was edited for this review.

Dates, for orientation:
- the model line's design (WORKFLOW.md, part 4) was approved by the user on 2026-11-16;
- the agent system was rebuilt from 2026-11-17, one agent at a time (docs 28 and 29);
- the blind rederivation tasks T04-T07 ran on 2026-11-18 (log, doc 30, the run record);
- doc 31 (2026-11-19) is a proposal the user has not approved.

Where your instructions speak of "the library" or of what is "recorded as fact", read: the derivation
records, the model equations, the C code and the conclusions that the thesis will rely on.

The file is long. Read it to the end, in parts if needed.
