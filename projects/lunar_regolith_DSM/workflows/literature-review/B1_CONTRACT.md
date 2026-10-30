# B1 run contract

Draft, 2026-10-29, for the user's approval at U0 (WORKFLOW.md §3). This file
fixes what the first branch reads, what it must produce, which stages run and
who approves what. The design it follows is WORKFLOW.md revision 10. Changes
after approval are made as a new version of this file, noted at the bottom.

## Question

**Q1, model-critical** (WORKFLOW.md §9, item 41). How can the stationary
grains be coupled with the evolving ice–vapor interface without a spurious air
film? The air film forms where the evolving ice field relaxes against the
stationary sediment field, so that φ_i + φ_s ≠ 1 in the band. The memo should
name the methods the nine papers offer for this (diffuse domain, smoothed
boundary, multiphase field with a sum constraint, coupled sharp and diffuse
descriptions), the conditions under which each was derived or tested, and
which of them fit the model.

**Secondary questions** (not model-critical; answered only from what the nine
papers say on the way):
- **Q2.** How do these methods impose boundary conditions for vapor diffusion
  (no flux, or a prescribed flux) on a diffuse stationary wall, and how
  accurate is the result as the interface width goes to zero?
- **Q3.** Do any of them treat surface diffusion along the ice surface, or
  other vapor-transport paths, and how?

## Papers

Nine papers: the seven seeds (WORKFLOW.md §9, item 38), Kaempfer & Plapp 2009
(the model the project builds on, and a calibration paper), and Moure & Fu
2024 in its arXiv v1 copy (the group's own model; its LaTeX source is
available). Hashes are sha256 of the files on disk on 2026-10-29. *Research*
is `~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/Projects/lunar_regolith_DSM/`.

| ID | DOI | File (under *Research*) | Pages | sha256 (first 12) | Status | Role (proposed) |
|---|---|---|---|---|---|---|
| li_2009 | 10.4310/cms.2009.v7.n1.a4 | `Literature/li_2009.pdf` (publisher, pp. 81–107) | 27 | 6c8ca5bbb242 | on disk; P1 run in development | core |
| lervag_lowengrub_2015 | 10.4310/cms.2015.v13.n6.a6 | `Literature/incoming/b1/zotero_CGTLZK54/10.4310_cms.2015.v13.n6.a6.pdf` | 32 | accd72942c5e | pulled from Zotero 2026-10-29; copy kind not checked (A2) | core |
| yu_chen_thornton_2012 | 10.1088/0965-0393/20/7/075008 | `Literature/incoming/b1/zotero_Z4WZUXEG/10.1088_0965-0393_20_7_075008.pdf` | 59 | df13f8d53001 | pulled from Zotero; copy kind not checked | core |
| folch_plapp_2005 | 10.1103/physreve.72.011602 | `Literature/incoming/b1/zotero_FVY8URBB/10.1103_physreve.72.011602.pdf` | 37 | 62caf0861949 | pulled from Zotero; copy kind not checked | core |
| joshi_jaiman_2018 | 10.1002/nme.5961 | `Literature/incoming/b1/zotero_8BJ7Y9VL/10.1002_nme.5961.pdf` | 25 | 01bbc5a8132d | pulled from Zotero; copy kind not checked | core |
| steinbach_2009 | 10.1088/0965-0393/17/7/073001 | `Literature/incoming/b1/zotero_57KIMZAG/10.1088_0965-0393_17_7_073001.pdf` | 32 | 87e43331819a | pulled from Zotero; copy kind not checked | context (review) |
| plapp_2011 | 10.1103/physreve.84.031601 | `Literature/incoming/b1/zotero_XCWNQ9Q7/10.1103_physreve.84.031601.pdf` | 36 | 74b1c5c05d1f | pulled from Zotero; copy kind not checked | context |
| kaempfer_plapp_2009 | 10.1103/PhysRevE.79.031502 | `Literature/kaempfer_plapp_2009.pdf` | 17 | 2141892c2e9c | on disk; P1 run in development | context (the base model; calibration paper) |
| moure_fu_2024_arxiv_v1 | 10.48550/arXiv.2407.21213 | `Literature/moure_fu_2024_arxiv_v1/moure_fu_2024_arxiv_v1.pdf`; source `source_2407.21213v1.tar.gz` (f1e3e9fa94cd) | 48 | e265eb8a28ac | on disk with LaTeX source | context (the group's model) |

*Core* papers are read with the user's review of key pages (U3); *context*
papers are read and checked the same way, but the user reviews only the claims
the memo uses. The Zotero copies' hashes match `zotero_pull.json`.

**Copy kinds.** A2 records whether each Zotero copy is the publisher's
version, an accepted manuscript or a preprint. Several page counts are larger
than the published articles are likely to be (Folch & Plapp 37, Plapp 36, Yu,
Chen & Thornton 59), which suggests manuscripts or preprints. Pages cited from
a copy that is not the version of record are given as PDF pages, and the copy
is named on every claim.

## Outputs

1. P1 content for each paper (WORKFLOW.md §6), written once.
2. Claims for each paper, checked by E3 and E4; the reader's symbol lists and
   the checker's applicability checklists as analysis records.
3. One memo, `memos/b1/Q1.vN.md`: statements citing claims; model options
   with their conditions; what the nine papers do not answer; leads (cited
   works that look essential, not read in B1). Q2 and Q3 are short sections
   of the same memo.
4. The calibration report: the diagnostic during development and the frozen
   final evaluation (WORKFLOW.md §5).
5. The repeat-run comparison on one paper (WORKFLOW.md §8), for the course.
6. Run manifests, and the dated backup with its restore test.

## Stages that run, and exceptions to the full design

Run: U0, A2, P1, E1, E3, E4 (steps 0, 1 and 2), E5, L1, Y1, V1 (X1 and U3),
Y2, Y3, U3, U4 (WORKFLOW.md §3).

Not run in B1 (WORKFLOW.md §13): discovery and screening (D1, D2, D2c, D3,
SC0 to SC2, U1, D2g, R1); acquisition by script and Zotero upload (A1, Z1;
the copies are already on disk); the second reader E2 and merging; LaTeX and
JATS text and source comparison (Moure & Fu's source is used only in T-res and
as a reference in X1); automated tables and figures; the passage index and the
reference desk; the general units operator layer; X1's random 20 % and U3b's
random 60; tiers and waves.

Order: build and test (WORKFLOW.md §12), calibration diagnostic, frozen final
evaluation, then the nine papers, then the memo.

## Stopping rule and budget

B1 stops when all nine papers have been read and checked and the memo names
the supported options, their conditions and the unresolved questions, or when
the build budget of about 3 to 5 working days runs out, whichever is first. In
the second case the remaining automation is cut, not the checks: the papers
left are read with the parts already built, and the user's review covers the
rest. After B1 the user returns to model work (Campaign B, then Pillar 2).

## Who approves what

| Item | Approves |
|---|---|
| This contract, the questions and the paper roles | user (U0) |
| Labels on each B1 gold statement naming the question that needs it | user (U0) |
| A copy that is not the version of record | user (A2) |
| The PDF parser and the page method (whole pages, tiles or zoom crops) | user, from T-parser and T-res results |
| Every equation or value the memo cites, when the reads disagree | user (V1) |
| The four meaning fields of every prose claim a model option rests on | user (U3) |
| Disputed claims; missed statements | user (U3) |
| The symbol and scaling table for each equation a model option uses | user |
| The memo | user (U4) |
| Prompts, schemas and model IDs frozen for the final evaluation | main session proposes, user approves |

## Open before approval

- ~~Confirm or change the core and context roles above.~~ Confirmed by the user, 2026-10-30.
- Label each of the 22 B1 gold statements (Kaempfer & Plapp 12, Li 10) with
  Q1, Q2, Q3 or none; recall is counted only on those labeled. Done 2026-10-30:
  Li 6 [Q1, Q2] and 4 none; Kaempfer & Plapp 1 Q1 and 11 none.
- ~~Confirm the arXiv v1 copy of Moure & Fu as the copy read~~ Confirmed 2026-10-30; the user judges it the same content as the journal version, and any Moure & Fu claim the memo relies on is checked in the journal version before the memo is final (its pages match
  the LaTeX source; the journal PDF has a different layout).
- ~~Confirm Q2 and Q3 as written.~~ Confirmed by the user, 2026-10-30.
- `papers.yaml` still lists Li's calibration roles as unset and the Zotero
  seeds without a file; they are filled when this contract is approved.

## Versions

- v1 draft, 2026-10-29.
- v1 approved by the user, 2026-10-30, including the arXiv v1 copy of Moure & Fu.
