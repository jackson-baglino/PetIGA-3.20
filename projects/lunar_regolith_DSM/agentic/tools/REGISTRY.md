# Tool registry

Where tools live: Claude Code has no recognised `.claude/tools/` folder (the `.claude/` folder holds configuration: agents, settings, hooks, skills). Our tools are ordinary tested code under `agentic/` and are listed here with who may call them and what they change. Headless agent workers have no tools at all: only drivers and the human call these.

| tool | kind | what it does | command | writes | tests |
|---|---|---|---|---|---|
| mesh_params | program | Sizes the interface width eps and the mesh from domain, temperature and kinetics (replaces preprocess/comp_eps.py). | `python -m agentic.tools.mesh_params --Lx 6e-4 --Ly 6e-4 --T0 -20 --alpha 1e-3` | none (pure) | agentic/tests/test_mesh_params.py |
| sympy_check | program | Recomputes every claim of a derivation exactly (equalities, derivatives, integrals, ODE residuals, Euler-Lagrange, limits); safe parser; reports which steps are only retyped algebra. | `python -m agentic.tools.sympy_check derivation.json  (on branch agent/math-lane: agentic/tools/sympy_check)` | none | agentic/tests/test_sympy_check.py (9) |
| run_lane | pipeline | Math workflow driver: spec -> derive -> symbolic check -> adversary -> optional blind replicate -> human gate. | `python -m agentic.math.run_lane --task T.json [--replicate]` | derivations/<task>/ in the workspace | agentic/tests/test_math_driver.py (4) |
| seed_test | program | Plants known errors in a correct derivation and records which layer catches each. | `python -m agentic.math.seed_test / seed_c3b` | stdout | - |
| f2 numerics suite | program | Toy and exact numerical checks of the pore-space formulation, plus a regression runner. | `python agentic/math/numerics/run_all_f2_checks.py` | stdout | 15 checks |
| discover / triage / fetch_oa | pipeline | Literature search (OpenAlex), abstract scoring by the smallest model, open-access PDF fetch. | `python -m agentic.lit.discover | triage | fetch_oa` | workspace registry, PDFs in Literature/*_incoming | - |
| prepare_text / run_extract / verify_record | pipeline | Page-marked text, one worker per paper, quote matching against the page. | `python -m agentic.lit.run_wave --manifest M.json` | workspace records | verification files |
| check_claims / check_images / consistency | pipeline | Independent claim check (Opus), page-image check of equations, agreement between two readings. | `python -m agentic.lit.check_claims | check_images | consistency` | verification files | - |
| synthesize / audit_memo | pipeline | Gate memo from verified records; citation audit of a memo. Both read the library through the desk. | `python -m agentic.lit.synthesize | audit_memo` | memos in the workspace | - |
| library.migrate / publish | program (curator) | Publishes records, texts, memos, searches and candidates into the library (hashed, versioned, read-only). | `python -m agentic.library.migrate` | library (the only writer) | agentic/tests/test_library.py |
| library.desk | program (reference desk) | The only way agents read the library; policy table, purpose required, hash-checked, logged. | `python -m agentic.library.desk ask --as <who> --why '...' <kind> k=v` | access log in the library (append only) | agentic/tests/test_library.py |
| library.verify / audit | program | Integrity check of every file; plain-language audit by the auditor agent. | `python -m agentic.library.verify | audit` | report outside the library | agentic/tests/test_library.py |
| writing.evidence / style / checks / render / run_report | pipeline | Fact ledger with sources, style profile, deterministic report checks, Word and PDF renderer, writing driver. | `python -m agentic.writing.run_report ...` | report files in the workspace | agentic/tests/test_writing_checks.py (10) |
