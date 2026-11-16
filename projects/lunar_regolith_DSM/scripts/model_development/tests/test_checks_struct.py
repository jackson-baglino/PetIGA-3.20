"""Tests for checks_struct.py (ledger typing and sourcing, sections) and mutate.py's mutation operators."""
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import checks_struct as cs  # noqa: E402
import mutate  # noqa: E402


def ledger(steps, symbols=({"symbol": "T", "meaning": "temperature", "units": "K"},)):
    return {"symbols": list(symbols), "steps": steps}


class LedgerTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.run = Path(self.tmp.name)
        (self.run / "derivation").mkdir()

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, led, md="S1 S2 S3"):
        (self.run / "derivation" / "ledger.json").write_text(json.dumps(led))
        (self.run / "derivation" / "derivation.md").write_text(md)

    def test_good_ledger(self):
        self.write(ledger([
            {"id": "S1", "type": "assumption", "statement": "k = const", "label": "A1", "if_wrong": "k(T)"},
            {"id": "S2", "type": "derived", "statement": "T_t = k T_xx", "from": ["S1"], "operation": "Fourier law"},
            {"id": "S3", "type": "choice", "statement": "FD or FEM", "label": "C1", "options": ["FD", "FEM"]}]))
        self.assertEqual(cs.check_ledger(self.run), [])

    def test_forward_reference_and_missing_fields(self):
        self.write(ledger([
            {"id": "S1", "type": "derived", "statement": "x", "from": ["S2"], "operation": "sub"},
            {"id": "S2", "type": "derived", "statement": "y"},
            {"id": "S3", "type": "cited", "statement": "z", "paper_id": "p"}]))
        probs = " | ".join(cs.check_ledger(self.run))
        self.assertIn("S1: cites S2, which is not an earlier step", probs)
        self.assertIn("S2: derived without", probs)
        self.assertIn("S3: cited without", probs)

    def test_step_not_in_markdown(self):
        self.write(ledger([{"id": "S1", "type": "assumption", "statement": "a", "label": "A1", "if_wrong": "b"},
                           {"id": "S12", "type": "assumption", "statement": "a", "label": "A2", "if_wrong": "b"}]),
                   md="only S12 here")
        self.assertIn("derivation.md does not mention steps ['S1']", " ".join(cs.check_ledger(self.run)))

    def test_sections_in_order(self):
        p = self.run / "derivation" / "spec_problem.md"
        p.write_text("\n".join(f"## {h}\ntext" for h in cs.M1_SECTIONS))
        self.assertEqual(cs.sections_in_order(p, cs.M1_SECTIONS), [])
        p.write_text("## Scales\n## Problem\n")
        probs = cs.sections_in_order(p, cs.M1_SECTIONS)
        self.assertTrue(any("out of order" in x for x in probs))
        self.assertTrue(any("Notation" in x for x in probs))

    def test_cli_exit_codes(self):
        with redirect_stdout(io.StringIO()):
            self.assertEqual(cs.main(["M3", "--run", str(self.run)]), 1)


class MutationOperatorTest(unittest.TestCase):
    def test_operators(self):
        kinds = dict(mutate.mutations_for("  R += 3.0*M*eps*grad - w*0.5;  // note + here"))
        self.assertEqual(kinds["minus_to_plus"].strip(), "R += 3.0*M*eps*grad + w*0.5;")
        self.assertIn("6.0", kinds["double_literal"])
        self.assertNotIn("plus_to_minus", kinds)          # the + is only in the comment
        self.assertEqual(mutate.mutations_for("x = y;"), [])


if __name__ == "__main__":
    unittest.main()
