# `docs/tex/` — the constraint set as LaTeX

Transcriptions of `preprocess/comp_eps.py`, for dropping into a talk or a
write-up. `comp_eps.py` remains the source of truth; if it changes, these
change with it.

| file | contents |
|---|---|
| `constraints_kinetics.tex` | $\rho_{vs}$, $D_v$, $d_0$, $D^*_{ia}$; $\beta_{HK}$ vs $\beta_{sub}$ (the scaled/unscaled trap); where $\alpha_c$ comes from |
| `constraints_eps.tex` | the width convention, the four bounds (B-HEAT, B-VAPOR, B-KINETIC, B-CURV), the mesh rule, and the reported-but-unenforced diagnostics |
| `constraints_derived.tex` | M&F SI Eq. (9): $\lambda_{sub}$, $\tau_{sub}$, $M_{sub}$, $\alpha_{sub}$, the $\alpha_{sub}/M_{sub}$ identity, and the physical constants table |
| `constraints_libbrecht.tex` | the Libbrecht $\sigma_0(T)$ chain and the numbers behind `studies/libbrecht_kinetics/` |
| `constraints_standalone.tex` | compilable wrapper: `pdflatex constraints_standalone.tex` |

The first four are **fragments**: no preamble, nothing above `\paragraph`, so
they `\input` cleanly into a beamer frame. They need `amsmath`; the standalone
wrapper also loads `geometry`, `amssymb`, `textcomp` and `fontenc`.

Cross-references run between fragments (`constraints_libbrecht.tex` cites
`eq:bkinetic` from `constraints_eps.tex`), so include them together, or drop
the `\eqref`s if you split them across slides.
