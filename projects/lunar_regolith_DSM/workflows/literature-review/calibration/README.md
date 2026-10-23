# Calibration sets for the literature-review workflow (C0 and T-res)

The user writes these files; the pipeline never does. They are the answer
keys the checks are measured against (WORKFLOW.md §5, C0, and §12, T-res).
Write in your own words with page numbers; no copied paragraphs. Agents never
see these files: each agent runs in an empty folder holding only the files a
script hands it.

| File | What it holds | Target size |
|---|---|---|
| `papers.yaml` | the calibration papers, with DOI and local file (filled in) | 4 papers |
| `gold_statements.yaml` | statements any good reading must find | about 30 over 3 papers |
| `planted_errors.yaml` | a correct claim and your altered version of it | at least 10 |
| `borderline_papers.yaml` | relevant papers a screener might wrongly exclude | 10–20 |
| `tres_pages.yaml` | equation-dense pages and their equations, for the resolution test | 3–4 pages |

Roles of the papers (user's choice, 2026-10-23):

- **Gold statements:** Kaempfer & Plapp 2009 (journal article, about 12),
  Molaro et al. 2019 (long journal article, about 10), Blackford 2007 (review,
  about 8).
- **Planted errors:** Kaempfer & Plapp 2009 and Molaro et al. 2019, the papers
  the user knows best.
- **T-res:** Moure & Fu 2024, whose arXiv version (2407.21213v1) includes the
  authors' LaTeX source, so its equations come from that source and the user
  only spot-checks them; plus one or two pages of Kaempfer & Plapp 2009
  transcribed by the user (a typeset two-column journal page, the realistic
  case).

Edit the files in any text editor. Lines starting with `#` are comments. Keep
the indentation. When a file is done, say so; a script checks it (schema,
pages exist in the PDF, no duplicates) before it is used.
