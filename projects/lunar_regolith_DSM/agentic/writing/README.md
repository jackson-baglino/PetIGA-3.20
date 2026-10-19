# agentic/writing: a writing workflow for short course reports

Built 2026-10-19 as the third workflow after literature (`agentic/lit`) and derivation (`agentic/math`). Goal: plain, concise text in the author's own voice, where every number traces to a source.

```
facts (evidence.py, no model) -> outline (Sonnet) -> draft with {F:id} placeholders (Sonnet) -> deterministic checks (checks.py, no model)
   -> two independent critics (Opus): clarity/voice and rubric -> revise (Sonnet), best version kept -> render (md, docx, pdf) + checklist for the human
```
- `evidence.py` collects facts from the repository, the ledgers, the records and the audit files. Each fact has a source; facts only the author can confirm are flagged.
- `style.py` measures the author's own writing (sentence length, pronouns) from a manuscript; `style_profile.json` holds numbers only.
- `checks.py` checks length, long sentences, raw numbers, unknown fact ids, unglossed jargon, the "I"/"we" voice, AI-sounding words, required sections and backticked paths.
- `run_report.py` runs the loop; `render.py` writes Markdown, Word and PDF.
- Agents: `.claude/agents/report-{outliner,drafter,critic,rubric-auditor}.md`.
Run: `python -m agentic.writing.run_report --brief B --guideline G --facts F --profile P --samples S --out-dir D`.
Tests: `python agentic/tests/test_writing_checks.py`.
Lessons from the first runs: revision rounds can make a draft longer and worse, so the driver keeps the best-scoring version; pasting sentence-like facts into placeholders garbles the prose, so those facts are paraphrased and tagged `<!-- F:id -->`; a thin ledger makes the drafter leave many author placeholders, so the ledger must hold what the report needs.
