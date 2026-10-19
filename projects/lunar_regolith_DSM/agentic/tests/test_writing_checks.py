"""Known-case tests of the deterministic report checks (no model calls)."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from agentic.writing import checks

F = {"n": {"id": "n", "value": 50, "text": "papers", "source": "x", "needs_human": False},
     "repo_url": {"id": "repo_url", "value": "https://x.org", "text": "link", "source": "y", "needs_human": True}}
BASE = "We read {F:n} papers. Done: a pipeline. Remaining: tests. Bottlenecks: time. Next steps: more. {F:repo_url} AI Use Declaration: we used Claude."

def kinds(d, **kw): return [p["kind"] for p in checks.check(d, F, **kw)["problems"]]

def test_clean_draft_passes():
    assert kinds(BASE) == [], kinds(BASE)

def test_render_fills_value():
    assert "50" in checks.render("We read {F:n} papers.", F)

def test_unknown_fact_flagged():
    assert "unknown_fact" in kinds(BASE + " {F:zzz}")

def test_raw_number_flagged_small_allowed():
    assert "raw_number" in kinds(BASE + " We read 4382 papers.")
    assert "raw_number" not in kinds(BASE + " We built 3 tools in 2026.")

def test_long_sentence_flagged():
    assert "long_sentence" in kinds(BASE + " " + " ".join(["word"] * 40) + ".")

def test_jargon_gloss():
    assert "jargon_unglossed" in kinds(BASE + " The phase-field model runs.")
    assert "jargon_unglossed" not in kinds(BASE + " The phase-field model (a way of drawing blurry boundaries) runs.")

def test_gloss_before_term():
    assert "jargon_unglossed" not in kinds(BASE + " Ice turns straight to vapor (sublimation) here.")

def test_voice_and_length():
    assert "voice" in kinds(BASE + " I think so.")
    assert "too_long" in kinds(BASE + " " + "word " * 1100 + ".", max_words=1000)

def test_missing_section():
    assert "missing_section" in kinds("We read {F:n} papers.")

def test_ai_sounding():
    assert "ai_sounding" in kinds(BASE + " This is a crucial step.")

if __name__ == "__main__":
    f = 0
    for n, fn in list(globals().items()):
        if n.startswith("test_"):
            try: fn(); print("PASS", n)
            except Exception as e: f += 1; print("FAIL", n, type(e).__name__, e)
    sys.exit(1 if f else 0)
