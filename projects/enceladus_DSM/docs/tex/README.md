# `docs/tex/` — the constraint set as LaTeX

Transcriptions of `preprocess/comp_eps.py`, for dropping into a talk or a
write-up. `comp_eps.py` remains the source of truth; if it changes, these
change with it.

## For PowerPoint (IguanaTeX)

**`constraints_iguanatex.txt`** — 13 numbered snippets, one constraint each,
separated by dashed rules. Copy one block into IguanaTeX's *New LaTeX display*
box and hit Generate. Each block is self-contained, uses only `amsmath` /
`amssymb` (already in IguanaTeX's default preamble), and every one is
regression-tested to compile inside `\[ ... \]`, which is what IguanaTeX wraps
it in. The header lines are `%`-commented, so sweeping one up in a sloppy copy
is harmless.

Blocks: [1] the four ε bounds on one slide · [2] ε as their minimum ·
[3] the mesh rule · [4] the thin-interface parameter δ · [5] β_HK vs β_sub ·
[6] the thermodynamic inputs · [7] Libbrecht's α_c · [8] the whole chain and
why it blows up · [9] the supersaturation gap · [10] the derived phase-field
parameters · [11] the width convention · [12] the unenforced diagnostics ·
[13] the kinetic fraction.

## For a written document

| file | contents |
|---|---|
| `constraints_kinetics.tex` | $\rho_{vs}$, $D_v$, $d_0$, $D^*_{ia}$; $\beta_{HK}$ vs $\beta_{sub}$ (the scaled/unscaled trap); where $\alpha_c$ comes from |
| `constraints_eps.tex` | the width convention, the four bounds (B-HEAT, B-VAPOR, B-KINETIC, B-CURV), the mesh rule, and the reported-but-unenforced diagnostics |
| `constraints_derived.tex` | M&F SI Eq. (9): $\lambda_{sub}$, $\tau_{sub}$, $M_{sub}$, $\alpha_{sub}$, the $\alpha_{sub}/M_{sub}$ identity, and the physical constants table |
| `constraints_libbrecht.tex` | the Libbrecht $\sigma_0(T)$ chain and the numbers behind `studies/libbrecht_kinetics/` |
| `constraints_standalone.tex` | compilable wrapper: `pdflatex constraints_standalone.tex` (run it twice for the cross-references) |

The first four are **fragments**: no preamble, nothing above `\paragraph`, so
they `\input` cleanly into a beamer frame. They need `amsmath`; the standalone
wrapper also loads `geometry`, `amssymb`, `textcomp` and `fontenc`.

Cross-references run between fragments (`constraints_libbrecht.tex` cites
`eq:bkinetic` from `constraints_eps.tex`), so include them together, or drop
the `\eqref`s if you split them across slides. The IguanaTeX snippets have no
cross-references at all, by design.

## Checking the snippets still compile

```bash
python3 - <<'EOF'
import re, pathlib, subprocess, tempfile, os
txt = pathlib.Path("docs/tex/constraints_iguanatex.txt").read_text()
blocks = re.findall(r"^-{80}\n(.*?)^-{80}$", txt, flags=re.S | re.M)
d = tempfile.mkdtemp(); bad = []
for i, b in enumerate(blocks, 1):
    src = ("\\documentclass{article}\n\\usepackage{amsmath,amssymb,amsfonts}\n"
           "\\begin{document}\n\\[" + b.strip() + "\\]\n\\end{document}\n")
    f = os.path.join(d, f"b{i}.tex"); open(f, "w").write(src)
    r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                        f"-output-directory={d}", f], capture_output=True, text=True)
    if r.returncode: bad.append(i)
print(f"{len(blocks)} blocks; " + ("all compile" if not bad else f"failed {bad}"))
EOF
```
