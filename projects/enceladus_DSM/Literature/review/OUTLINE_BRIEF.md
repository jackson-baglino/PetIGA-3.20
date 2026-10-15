# Phase 2 brief: write manuscript_outline_new_<X>.tex (reviewers A, B, C)

You are reviewer X, working **independently**. Do not open another reviewer's folder or file.
Root: /Users/jacksonbaglino/PetIGA-3.20/projects/enceladus_DSM

## What changed since phase 1
- **The draft has been updated**:
  - **Current draft:** `Literature/manuscript/k_eff-ice-sintering-paper/manuscript_draft_v2.tex` (+ `Eff-Thermal-Cond-Paper.bib`, `Tables/tables.tex`).
  - **Previous draft:** `Literature/manuscript_2026-10-14/manuscript_draft_v2.tex`.
  - **What changed:** the Introduction was restructured. It has a new paragraph 4 on Enceladus k_eff/Choukroun, and paragraphs were renumbered 5–9. Figure 1's caption is now vapor-only. In Methods, T is in kelvin (eqs. for P_vs and D_v), the latent term uses rho_i, and P_atm is renamed. There is a new appendix "Model Scaling Parameters". The bib gained Hobbs_and_Mason, kingery_regelation_1960, thomas_crushing_1994 and noauthor_sintering_nodate.
  - **Before you start:** run `diff` on the two drafts. Your phase-1 claims table uses the old line numbers, so re-map every row to the current draft, and drop or update rows whose text changed.
- **The outline to revise:** `studies/keff_sintering/manuscript_outline.tex`. It is a content plan that says what each paragraph should include.

## Your inputs (keep context small)
- Your phase-1 work: `Literature/review/X/claims_table.md` (A and C have one; B does not) and `Literature/review/X/notes/`.
  **Do not `cat` the whole notes folder.** Use `grep -l`/`grep -n` to find the notes relevant to a claim, then open those.
- `Literature/review/REVIEWER_BRIEF.md`: the conventions, hard rules and Addendum still apply. That means no web access, page evidence
  for every verdict, the notation rules, k_eff as steady-state only, Molaro exponents as our fits, no surface diffusion in our model, and
  cite only closely matching work.
- `Literature/review/TOKEN_AUDIT.md`: many notes for uncited papers are **title/abstract-level stubs**, not full reads.
  Before you base any suggestion on a note, check that the note quotes page-located content. If it does not, open the
  relevant pages of `Literature/text/<key>.txt` (use `grep -n` to locate them, then `sed -n` on a page range; never `cat` a whole paper).
- For any number, sign, unit or exponent, confirm it on the rendered page:
  `mutool draw -q -r 110 -o <scratch>/<key>_pN.png Literature/papers/<key>.pdf N`, then Read the PNG.

## Steps
1. (B only) Build `Literature/review/B/claims_table.md` from your notes against the current draft and outline,
   in the format given in REVIEWER_BRIEF.md.
   (A and C) Update your claims table to the current draft. Add rows for new or changed text, especially the new Intro paragraph 4.
2. Write `studies/keff_sintering/manuscript_outline_new_X.tex`: a full copy of the outline with your changes applied.
   - Keep the outline's structure and macros. Add `\definecolor{revX}{RGB}{<your colour>}`,
     `\newcommand{\add}[1]{\textcolor{revX}{#1}}` and `\newcommand{\del}[1]{\textcolor{red}{[del: #1]}}`.
   - Put a `% [REVIEW-X #id] reason; evidence: key p.N "quote"` comment before every change. The `#id` must match the claims table.
   - The outline also tracks the draft. Where the outline's description of the draft text is stale, update it to the current draft.
   - **Prioritise.** At the top, after `\begin{document}`, add a short section "Reviewer X: top findings" listing the
     ≤ 10 most important corrections or additions, in order of importance. Do not pad: a loosely related citation is worse than none.
3. Compile in your scratch directory:
   `pdflatex -interaction=nonstopmode -halt-on-error -output-directory <scratch> studies/keff_sintering/manuscript_outline_new_X.tex`.
   Fix any errors. Do not leave build files in the repo.
4. Report back in ≤ 300 words: verdict counts, your top 10, and any claims you could not check.

## Addendum: suggested manuscript snippets (user request, 2026-10-15)
Where a finding calls for new or revised manuscript text, include a ready-to-paste snippet of draft prose,
**written in the author's voice**. Match the current draft's style:
- first-person plural ("we");
- plain, direct sentences, with `\citet`/`\citep` citations;
- the same notation and macros;
- no filler or hedging stacks.

Learn the voice from the current draft's Introduction and Methods before you write. In the outline, put each snippet in
`\Snippet{...}`, defined as
`\newcommand{\Snippet}[1]{\par\smallskip{\small\textit{Suggested text:} #1}\par\smallskip}`, directly under the paragraph it belongs to.
Tag it with `% [REVIEW-X #id]`. Every factual statement and citation in a snippet must be backed by evidence in your claims table.
Snippets are for the important points only: corrections, missing context and key citations. Do not rewrite paragraphs wholesale.
