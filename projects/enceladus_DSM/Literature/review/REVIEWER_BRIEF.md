# Reviewer brief — k_eff ice-sintering manuscript fact-check

Three reviewers (A, B, C) work from this same brief, each **independently**. Do not open
another reviewer's folder (`Literature/review/<other letter>/`) or its `manuscript_outline_new_<other>.tex`.
Write only inside your own folder and your own outline_new file.

## Hard rules
- **No web access** (no WebSearch, WebFetch or curl). Every source is local.
- Do not edit `manuscript_outline.tex`, the draft, the .bib, or anything in `papers/` / `text/`. Do not commit.
- Use page evidence for every verdict: `<citekey> p.N` (the PDF page number from the
  `=== page N ===` markers, plus the printed journal page when it differs), with a short quote of ≤ 2 sentences.
- Cite only closely matching work. Do not add loosely analogous citations; a suggested citation must
  directly support the specific sentence.
- Take DOIs only from the PDF itself or the bib. Never invent them.

## Inputs (root: /Users/jacksonbaglino/PetIGA-3.20/projects/enceladus_DSM)
- **Outline (primary target):** `studies/keff_sintering/manuscript_outline.tex`. This is the live content plan:
  for each paragraph it says what the paragraph should include. Citations appear as
  `\texttt{key}` inside `\Cite{...}` boxes, plus author-year names in prose.
- **Draft:** `Literature/manuscript/k_eff-ice-sintering-paper/manuscript_draft_v2.tex`, its bib
  `Eff-Thermal-Cond-Paper.bib`, and `Tables/tables.tex`. The draft's old Discussion/Conclusion are
  notes only. Check **every `\cite`** in the draft; claims without citations are checked in the outline.
- **Library:** `Literature/MANIFEST.md` (citekey → file, sources, known gaps), `Literature/text/<key>.txt`
  (full text with page markers) and `Literature/papers/<key>.pdf`. Use the PDF (Read tool, `pages`)
  whenever equations, symbols, tables or figure values are garbled in the text file. Check
  `Literature/review/CHECK_COLLECTION.md` if it exists, since it lists the papers whose text layer is unreliable.
- Missing papers (no local copy): Hobbs & Mason 1964 (HTML abstract only), Kuczynski 1949, Kingery & Berg 1955,
  Murphy & Koop 2005, Seiferlin 1996, Lin 2023, Pokrifka 2023, Asenath-Smith 2025. Mark claims that rest on these
  as `not-checkable (missing PDF)`. You may note what other library papers say about the same point.

## Project conventions the suggestions must respect
- Notation: no index notation, a single dot for products, `\phi` = phase field, `\varphi` = porosity,
  `\zeta` = indicator, ice/air subscripts.
- k_eff results are steady-state homogenization only. Do not use them to justify claims about how
  metamorphism itself evolves beyond what the manuscript computes.
- The power-law sintering exponents attributed to Molaro et al. are **our fits** to their data, not values they report.
- Our model has no surface diffusion. Do not suggest claims that imply it does.
- Manuscript framing: lab-observable sintering, with no extrapolation to Enceladus temperatures. The Intro surveys
  metamorphism broadly.
- Round-1 review style: correctness, typos, major gaps, flow, concision. No filler.

## Procedure
1. Read the outline and the draft in full first. Build a list of every factual claim and every citation.
2. **Read every paper in the library start to finish**, in the order given in your launch message.
   The Auriault et al. 2009 book (465 pp) comes last; read it in full as well.
   **Immediately after finishing each paper**, write `Literature/review/<X>/notes/<citekey>.md` containing:
   - a 2–4 line summary: what was done, method, conditions (T, density, grain size, geometry);
   - the findings relevant to this manuscript, each with a page locator and a short quote or number;
   - which outline/draft claims it supports, contradicts or qualifies;
   - candidate new citations or claims for the manuscript, if any.
   **Resume rule:** if your context is compacted or you are restarted, list `notes/` and skip papers that
   already have a notes file. Notes on disk are your memory.
3. Fill in `Literature/review/<X>/claims_table.md` with one row per claim or citation:
   `id | location (outline line / draft line) | claim (short) | cited keys | verdict | evidence (key p.N + quote) | suggested fix`.
   The verdict is one of: `supported`, `partially`, `unsupported`, `wrong-citation`, `contradicted`, `not-checkable`.
   Also check numbers (values, exponents, temperatures, densities, k_eff ratios) against the source.
4. Add a section `## Suggested additions` to the claims table: new citations to existing claims, and new
   claims or context the manuscript should include, each with evidence. Add `## Typos / bib problems` too.
5. Write `studies/keff_sintering/manuscript_outline_new_<X>.tex` as a full copy of the outline with your
   suggestions applied inline:
   - add `\definecolor{revX}{RGB}{...}` plus `\newcommand{\add}[1]{\textcolor{revX}{#1}}` and
     `\newcommand{\del}[1]{\textcolor{red}{[del: #1]}}`, so deletions remain visible;
   - put a `% [REVIEW-<X> #id] reason + evidence` comment on the line before each change, where `#id` matches the claims table.
   - Compile it in a scratch directory to check it builds:
     `pdflatex -interaction=nonstopmode -halt-on-error -output-directory <scratchpad>/<X> <file>`.
     Do not leave .aux/.log files in the repo.
6. Final report (returned to the coordinator): counts by verdict; the 10 most important problems;
   the top suggested additions; papers you could not read fully, and why.

## Addendum (after the collection check, see `CHECK_COLLECTION.md`)
- `libbrechtMeasurementsSurfaceAttachment2013` was the wrong paper (arXiv:1302.1231). It has been replaced
  with the correct J. Cryst. Growth 377:1–8 article in both `papers/` and `text/`. If you have already
  written notes for it, redo them.
- **The Read tool's `pages` option does not work** (no pdftoppm). To see a PDF page, render it:
  `mutool draw -q -r 110 -o <your scratch>/<key>_pN.png Literature/papers/<key>.pdf N`, then Read the PNG.
- **The text layer silently drops minus signs, μ and ×** in many Wiley/AGU/Elsevier/AIP PDFs: −10 °C reads "10 °C",
  10 μm reads "10 mm", and 3.19×10⁻¹⁰ reads "3.1910−10". Superscripts are flattened everywhere (10⁴ reads "104").
  OCR'd older papers have digit errors. **Before writing any verdict that hinges on a number, sign, unit, exponent
  or equation, confirm it on the rendered page image.** The list of affected papers is in `CHECK_COLLECTION.md`.
