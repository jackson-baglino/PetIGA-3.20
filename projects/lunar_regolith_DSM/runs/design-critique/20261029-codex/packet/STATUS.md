# Status of the literature-review build, 2026-10-29

WORKFLOW.md still says "Nothing is built yet"; that line is out of date. This
note is the current state.

## Who and why
- One user, a PhD student, running a phase-field model of ice sintering and
  sublimation in regolith (PETSc/PetIGA). The literature review exists to
  serve that model: find and check the physics, parameters and numerical
  methods the model and its manuscripts rely on.
- The model is the priority. A conference is about early December 2026.
  Time spent on the review pipeline competes directly with model work.
- The work runs in Claude Code. Readers and checkers are Claude calls; Codex
  is not part of the pipeline (revision 9).

## Built and tested (59 unit tests pass)
- `check_calibration.py`: validates the user's calibration files.
- `make_plants.py`: makes mechanical planted errors (changed digit, flipped
  sign or exponent, wrong unit, wrong page) from the user's gold statements.
- `calibrate.py`: split into tuning and held-out sets (plants from one gold
  statement stay on one side), scoring with Clopper-Pearson bounds.
- `zotero_pull.py`: copies PDFs the user saved in Zotero, via Zotero's local
  API only, with sha256 records. Used live on 2026-10-29 for the six B1 seeds.
- `prepare.py` (P1, first part) and `marker_run.py`: 300 dpi page images,
  tiles, Marker text per page (offline, no LLM), displayed equations with
  their printed numbers, headings, PDF-to-printed page map, manifest with
  hashes and tool versions. Runs under a macOS sandbox profile with no
  network.
- `check_permissions.py` and the permission files.

## Calibration material written by the user
- 40 gold statements on 4 papers (Kaempfer & Plapp 2009, Li et al. 2009,
  Molaro et al. 2019, Blackford 2007).
- 11 hand-planted errors of judgment kinds (scope widened, contradicted,
  dropped condition, wrong attribution, ...).
- 57 script plants, all approved by the user.
- T-res pages (equation-dense pages with known content): in progress.
- Retrieval questions (T-ret): deferred until the library desk exists.

## Lessons from the first real runs (P1)
- The first tiler cut pages using blank rows of ink only. Coverage was 1.0
  and no page failed, yet a human look found a table cut mid-row, a plot
  separated from its axis labels and a caption spread over three tiles. Tiles
  are now cut along Marker's layout blocks (tables and figures with captions
  are never cut).
- P1 reported 47 numbered equations for Li et al.; the user counted 58.
  Marker had kept every number in a form the parser missed. Fixed; 58 found.
- Marker takes roughly 3 to 6 minutes per paper on this laptop.

## Not built yet (planned order)
1. Rest of P1/P2: text from arXiv LaTeX source, tables (two extractions
   compared), figure crops, source comparison, `curator.py` publishing,
   passage index (`build_index.py`).
2. Reader and checker prompts, blind reader, math OCR, unit checks,
   comparison of reads.
3. Calibration run (plants scored against readers and checker).
4. Launcher, sandbox and permission tests.
5. The first branch (B1, numerical methods; 9 seed papers), then the desk
   that answers questions from the library, then the retrieval test.
