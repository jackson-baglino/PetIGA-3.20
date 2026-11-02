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
2024 in its published copy (the group's own model; the arXiv v1 LaTeX source
stays on disk; X1 does not compare it yet). Hashes are sha256 of the files on disk (v2,
2026-11-01). *Research*
is `~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/Projects/lunar_regolith_DSM/`.

| ID | DOI | File (under *Research*) | Pages | sha256 (first 12) | Status | Role (proposed) |
|---|---|---|---|---|---|---|
| li_2009 | 10.4310/cms.2009.v7.n1.a4 | `Literature/li_2009.pdf` (publisher, pp. 81–107) | 27 | 6c8ca5bbb242 | on disk; A2 2026-11-01: no DOI on pp. 1–2, publisher line found (kind unknown; user decides at A2) | core |
| lervag_lowengrub_2015 | 10.4310/cms.2015.v13.n6.a6 | `Literature/incoming/b1/zotero_CGTLZK54/10.4310_cms.2015.v13.n6.a6.pdf` | 32 | accd72942c5e | pulled from Zotero 2026-10-29; arXiv, the only copy (user decides at A2) | core |
| yu_chen_thornton_2012 | 10.1088/0965-0393/20/7/075008 | `Literature/incoming/b1/zotero_FPFXFI8J/10.1088_0965-0393_20_7_075008.pdf` | 42 | b3e677615fcc | published copy, added to Zotero by the user 2026-11-01 (replaces Z4WZUXEG); A2 2026-11-01 (dev-20261101-a2): publisher, DOI found | core |
| folch_plapp_2005 | 10.1103/physreve.72.011602 | `Literature/incoming/b1/zotero_9QQPQWFH/10.1103_physreve.72.011602.pdf` | 27 | 06d68698a9d2 | published copy, added 2026-11-01 (replaces FVY8URBB); A2 2026-11-01 (dev-20261101-a2): publisher, DOI found | core |
| joshi_jaiman_2018 | 10.1002/nme.5961 | `Literature/incoming/b1/zotero_C66WP9T7/10.1002_nme.5961.pdf` | 25 | f80cfd405907 | published copy, added 2026-11-01 (replaces 8BJ7Y9VL); A2 2026-11-01 (dev-20261101-a2): publisher, DOI found | core |
| steinbach_2009 | 10.1088/0965-0393/17/7/073001 | `Literature/incoming/b1/zotero_57KIMZAG/10.1088_0965-0393_17_7_073001.pdf` | 32 | 87e43331819a | pulled from Zotero; A2 2026-11-01 (dev-20261101-a2): publisher, DOI found | context (review) |
| plapp_2011 | 10.1103/physreve.84.031601 | `Literature/incoming/b1/zotero_8LSUVADQ/10.1103_physreve.84.031601.pdf` | 15 | 0713491c5dd2 | published copy, added 2026-11-01 (replaces XCWNQ9Q7); A2 2026-11-01 (dev-20261101-a2): publisher, DOI found | context |
| kaempfer_plapp_2009 | 10.1103/PhysRevE.79.031502 | `Literature/kaempfer_plapp_2009.pdf` | 17 | 2141892c2e9c | on disk; A2 2026-11-01: publisher, DOI found | context (the base model; calibration paper) |
| moure_fu_2024 | 10.1021/acs.cgd.4c00539 | `Literature/incoming/b1/zotero_7663PFVQ/10.1021_acs.cgd.4c00539.pdf` | 14 | 61268456ef25 | published copy (user decision 2026-11-01; replaces the arXiv v1 copy, whose LaTeX source stays on disk; it was used in T-res and is not yet compared in X1); A2 2026-11-01 (dev-20261101-a2): publisher, DOI found | context (the group's model) |

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
- ~~Confirm the arXiv v1 copy of Moure & Fu as the copy read~~ Confirmed 2026-10-30, then replaced by the published copy on 2026-11-01 (v2); the user judges it the same content as the journal version, and any Moure & Fu claim the memo relies on is checked in the journal version before the memo is final (its pages match
  the LaTeX source; the journal PDF has a different layout).
- ~~Confirm Q2 and Q3 as written.~~ Confirmed by the user, 2026-10-30.
- `papers.yaml` still lists Li's calibration roles as unset and the Zotero
  seeds without a file; they are filled when this contract is approved.

## Status

### Line status (generated)

<!-- BEGIN GENERATED: lr_status.py -->
Run `dev-20261101-a2`. B1 · DEVELOP · run dev-20261101-a2 · next: P1 yu_chen_thornton_2012 · 1 gate(s) waiting · 0 review item(s) undecided

| Stage | Station | Done | Gate | Waiting |
|---|---|---|---|---|
| U0 | built | 1/1 | blocking |  |
| C0 | built | 1/1 | blocking |  |
| A2 | built | 7/8 | blocking | 1 |
| P1 | built | 0/8 | none |  |
| E1 | built | 0/8 | none |  |
| E3 | built | 0/8 | review |  |
| E4 | built | 0/8 | none |  |
| E5 | built | 0/8 | review |  |
| L1 | built | 0/8 | review |  |
| Y1 | built | 0/1 | none |  |
| X1 | built | 0/1 | blocking |  |
| U3a | built | 0/1 | blocking |  |
| L1b | built | 0/8 | none |  |
| Y2 | built | 0/1 | none |  |
| Y3 | built | 0/1 | review |  |
| U4 | built | 0/1 | blocking |  |
| BK | built | 0/1 | none |  |
<!-- END GENERATED -->

The line status above is computed by `scripts/literature_review/lr_status.py`
from the run folder of the newest run made with `run_b1.py`; it is not edited
by hand. The table below is the development history, written by hand up to
2026-11-01, and is kept until the first run under `run_b1.py` replaces it.

### Development history (hand-written, to 2026-11-01)

Where B1 stands, one row per node of `flowchart.json` (every node, so the
progress chart can be drawn from this table) and one row per contract output
that has no node. Read at the start of every session; updated in the same
commit as the work it records. The "Remaining" list in `.claude/ACTIVITY_LOG.md`
is copied from here, not written from memory. Rendered as
`flowchart_b1_progress.html` by
`scripts/flowchart/render_flowchart.py workflows/literature-review/flowchart.json --progress workflows/literature-review/B1_CONTRACT.md`.

Status words: **done**; **dev** (done only in a development run, before the
prompts are frozen, so not a B1 result); **partial**; **not-run** (built, not
yet run); **not-built**; **outside** (not run in B1, see the stages section).

Here: Y3
Next by the contract: C0

"Here" is where the development run `runs/literature-review/b1-dev2` has
reached. "Next by the contract" is where the contract's order puts us: the
frozen final evaluation comes before the nine papers and the memo.

| Node | Stage | Status | Evidence, or what is missing |
|---|---|---|---|
| UB | B0 branches and seeds | done | B1 seeds chosen 2026-10-28; recorded in `calibration/papers.yaml` |
| RBR | branches.yaml | partial | no `branches.yaml`; the B1 seeds live in `calibration/papers.yaml` |
| U0 | Brief | done | this contract, approved 2026-10-30 |
| RB | brief.vN.json | done | this contract stands in for the brief |
| D1 | discover.py | outside | |
| D2 | morty-scout | outside | |
| RS | search snapshots | outside | |
| D3 | resolve.py | outside | |
| RC | candidates.json | outside | |
| SC0 | screen_baseline.py | outside | |
| C0 | calibrate.py | partial | diagnostic done on the tuning set (calib-dev5: plants 14/14, controls flagged 5/30); prompts not frozen, held-out evaluation not run |
| SC1 | jerry-screener | outside | |
| SC2 | compare_screens.py | outside | |
| U1 | approve shortlist | outside | |
| RL | shortlist.json | outside | |
| A1 | acquire.py | outside | copies already on disk |
| RP | Literature/incoming | done | nine PDFs on disk; hashes in the Papers table |
| U2 | supply missing papers | done | Zotero pull 2026-10-29 |
| A2 | check_copy.py | partial | Moure & Fu arXiv v1 approved by the user; script not built; seven Zotero copies "copy kind not checked" |
| P1 | prepare.py | partial | parser (MinerU) and page method (whole pages, 2,000 px) chosen; P1 content only for Li, Kaempfer & Plapp, Molaro (p1-dev); reads use page images |
| RG | content/ | partial | as P1 |
| P2 | build_index.py | outside | |
| RI | index/ | outside | |
| Z1 | zotero_add.py | outside | |
| E1 | morty-reader | dev | 9/9 papers, b1-dev2; E2 second read outside |
| RR | reads | dev | b1-dev2/untrusted/<paper>/read.json |
| E3 | quote check | dev | 9/9, in read_check.py (not a separate check_extraction.py) |
| E3b | check_units.py | outside | |
| E4 | beth-checker | dev | 9/9, steps 0-2 (last four 2026-10-30) |
| RV | verdicts.json | dev | b1-dev2 |
| X1 | blind re-read (V1) | dev | dev-20261101-kp: 6/6 agree; U3a requires X1 |
| E5 | revise failed claims | dev | dev-20261101-kp: one claim disputed (c006), sent to U3 |
| U3a | guided reading (U3) | dev | dev-20261101-kp: u3_review.py, 6 decisions by the user 2026-11-02 (model-option claims and the disputed claim; not the full guided reading) |
| U3b | random audit | outside | B1 runs only U3's cited-claims part, not the random 60 |
| L1 | curator.py | dev | dev-20261101-kp |
| L1b | curator.py --after-u3 | not-run | publishes the user's U3 decisions (2026-11-02); added after the dev run's Y2 |
| LB | library/ | not-built | |
| L2 | audit.py | outside | |
| BK | backup archive | dev | dev-20261101-kp: 491 files, restore test passed |
| R1 | retag reader | outside | |
| Y1 | beth-synthesizer | dev | Q1 v1, v2 (revise step); Q2/Q3 sections and leads missing |
| Y2 | check_citations.py | dev | v2 passes |
| Y3 | skeptic | dev | v2: 88 supported, 12 overstated, 1 unsupported; 12 omissions (1 high) |
| RM | memos | dev | in b1-dev2/untrusted/memos, not `memos/b1/` |
| U4 | user decides | dev | dev-20261101-kp: approved (dev run) 2026-11-02 |
| — | run manifests (output 6) | not-built | no manifest.json in b1-dev2 or calib-dev*; per-call usage is in calls.jsonl |
| — | repeat-run comparison (output 5) | not-run | |
| — | restore test (output 6) | dev | T-restore passed 2026-11-02: 208 files restored, 0 hash mismatches, Y2 unchanged |

## Versions

- v1 draft, 2026-10-29.
- v1 approved by the user, 2026-10-30, including the arXiv v1 copy of Moure & Fu.
- v2, 2026-11-01: the user chose the published copy of Moure & Fu (ID `moure_fu_2024`) and published copies of Yu, Chen & Thornton, Folch & Plapp, Joshi & Jaiman and Plapp, which the user added to Zotero after A2 found preprint copies; Lervåg & Lowengrub stays the arXiv copy (the only one available). Rows recomputed from the Zotero originals; `run_b1.py` checks each path and hash.
