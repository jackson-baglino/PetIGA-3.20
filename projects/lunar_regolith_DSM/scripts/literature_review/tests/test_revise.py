"""Tests for revise.py (E5) with a fake agent call: revise, withdraw, dispute after two rounds."""
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import revise  # noqa: E402

GOOD_QUOTE = "the diffuse domain method converges at first order in the interface width"
PAGE = f"Some text. {GOOD_QUOTE} as shown in the analysis. Other text."


def claim(cid, quote="not on the page at all but long enough words here", page=1):
    return {"id": cid, "text": "t", "kind": "result", "quote": quote, "page_pdf": page, "location": "s2",
            "questions": ["Q1"], "basis": "derived", "conditions": [], "version": 1}


def verdict(cid, text="supported"):
    return {"claim_id": cid, "text": text, "quote_on_page": text == "supported", "equation": "none",
            "parameter": "none", "dropped_conditions": [], "reason": "r"}


class Rounds(unittest.TestCase):
    def test_revise_withdraw_dispute(self):
        with tempfile.TemporaryDirectory() as tmp:
            textdir = Path(tmp)
            (textdir / "p001.txt").write_text(PAGE)
            state = {}
            for cid in ("B", "C", "D"):
                state[cid] = {"claim": claim(cid), "verdict": verdict(cid, "overstated"),
                              "e3": {"ok": False}, "status": "pending", "history": []}
            state["A"] = {"claim": claim("A", GOOD_QUOTE), "verdict": verdict("A"), "e3": {"ok": True},
                          "status": "confirmed", "history": []}
            calls = []

            def call(show, text, schema, agent, **tag):
                calls.append((tag["step"], tag["round"]))
                if schema == "revise":
                    out = []
                    for cid in ("B", "C", "D"):
                        if f'"id": "{cid}"' not in text:
                            continue
                        if cid == "C":
                            out.append({"claim_id": "C", "action": "withdrawn", "reason": "not supported"})
                        else:
                            out.append({"claim_id": cid, "action": "revised", "reason": "narrowed",
                                        "claim": dict(claim(cid, GOOD_QUOTE), text="narrower")})
                    return {"revisions": out}
                ids = [cid for cid in ("B", "D") if f'"id": "{cid}"' in text]
                return {"verdicts": [verdict(cid, "supported" if cid == "B" else "unsupported") for cid in ids]}

            rounds = revise.run_rounds(state, {1: Path("p001.png")}, textdir, {"items": []}, call)
            self.assertEqual(rounds, 2)
            self.assertEqual({k: s["status"] for k, s in state.items()},
                             {"A": "confirmed", "B": "revised", "C": "withdrawn", "D": "disputed"})
            self.assertEqual(state["D"]["claim"]["version"], 3)
            self.assertEqual(state["B"]["claim"]["version"], 2)
            self.assertEqual(calls, [("E5.read", 1), ("E5.check", 1), ("E5.read", 2), ("E5.check", 2)])

    def test_confirmed_rule(self):
        self.assertTrue(revise.confirmed(verdict("x"), {"ok": True}))
        self.assertFalse(revise.confirmed(verdict("x"), {"ok": False}))
        self.assertFalse(revise.confirmed(dict(verdict("x"), dropped_conditions=[{"text": "c"}]), {"ok": True}))
        self.assertFalse(revise.confirmed(dict(verdict("x"), equation="corrected"), {"ok": True}))
        self.assertFalse(revise.confirmed(None, {"ok": True}))


if __name__ == "__main__":
    unittest.main()
