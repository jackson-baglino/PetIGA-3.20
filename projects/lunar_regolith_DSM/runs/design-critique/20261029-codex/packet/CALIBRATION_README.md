# Calibration sets for the literature-review workflow (C0, T-res, T-ret)

The user writes these files; the pipeline never does. They are the answer
keys the checks are measured against (WORKFLOW.md §5 C0, §12 T-res and T-ret).
Write in your own words, with page numbers; no copied paragraphs. Agents never
see these files: each agent runs in an empty folder holding only what a
script hands it.

| File | What you write | Size | Needed for |
|---|---|---|---|
| `planted_errors.yaml` | a true claim and your altered version of it | 10 to 12 | core C0 (first) |
| `tres_pages.yaml` | Moure & Fu page numbers; one or two Kaempfer & Plapp pages transcribed | 2–3 pages + 1–2 pages | core C0 |
| `gold_statements.yaml` | facts any good reading must find | Kaempfer & Plapp 12, Li 10 (B1); Molaro 10 (done), Blackford 8 (B3) | branch C0 |
| `tret_questions.yaml` | questions you would ask of the B1 papers, with answer pages | 25 | the search test for B1 (also Lab 5) |
| `borderline_papers.yaml` | relevant papers a screener might wrongly exclude | 10–20 | not now: first large branch (B4) |
| `papers.yaml` | the calibration papers (filled in) | — | — |

What is not yours to write:
- the plants drafted by `beth-planter` (Opus), which you approve one by one;
- the mechanical plants a script makes from verified claims.

Both are generated into `generated/`, for you to approve:
- `generated/script_plants.yaml`: made by `make_plants.py` from your gold
  statements. There are four kinds: a changed digit, a flipped sign or
  exponent, a swapped unit, and a page moved by one or two. Identifiers such
  as Fig. 1 or SA81 are never changed. For each plant, set `approved: true`
  if it is wrong in exactly the way its note says and isn't accidentally
  still true; otherwise set it to `false`. Rerunning the script overwrites the
  file, so approve plants once your gold statements are final.
- `generated/model_plants.yaml`: drafted by `beth-planter`, later.

When the plants are in, `calibrate.py split --heldout-size N` fixes the
tuning and held-out sets once, before any scoring. All your own plants are
held out. To pass, the held-out set needs at least 19 plants that reach the
checker (E4): your 10 or more, plus approved content plants. Page plants
don't count toward the 19, because a script catches them.

Order that unblocks the most:
1. `planted_errors.yaml`
2. the T-res pages
3. the B1 gold statements (Kaempfer & Plapp, Li)
4. `tret_questions.yaml`, once the B1 seeds are chosen
5. Blackford (B3), which can wait until B3 is seeded

Roles of the papers (user's choice, 2026-10-23; branches 2026-10-27):
- **Kaempfer & Plapp 2009** (B1, B2): gold statements, planted errors, T-res.
- **Li et al. 2009** (B1, diffuse domain; this week's homework): gold
  statements; other roles are the user's choice.
- **Molaro et al. 2019** (B3): gold statements, planted errors.
- **Blackford 2007** (B3, a review): gold statements.
- **Moure & Fu 2024** (B1): T-res, from the authors' LaTeX source
  (arXiv 2407.21213v1), which you only spot-check.

**Editing.** Edit the files in any text editor.
- Lines starting with `#` are comments.
- Keep the indentation.
- Leave unused slots empty, or delete them.
- Wrap text that contains a colon followed by a space in single quotes.

Check a file at any time, as often as you like:

    venv_workflows/bin/python scripts/literature_review/check_calibration.py

It reports, per entry, YAML errors (with the line), missing or invalid fields,
pages that do not exist in the paper (printed page numbers are mapped to the
PDF), and duplicates. It never prints paper text. Empty slots are ignored.
