"""compare_reads.py (X1): agreement needs every available read equal and none missing."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import compare_reads as cr  # noqa: E402

EQ = {"id": "p-c001", "page_pdf": 3, "equation": {"latex": r"\partial_t \phi = M \nabla^2 \mu", "number": "(12)"}}
VAL = {"id": "p-c002", "page_pdf": 4, "parameter": {"quantity": "alpha", "value": "1e-3", "unit_printed": "m/s"}}


class Compare(unittest.TestCase):
    def test_equation_formatting_agrees_sign_does_not(self):
        same = {"latex": r"\partial_{t}\phi=M\nabla^{2}\mu"}
        r = cr.compare_item(EQ, same, same, {(3, "12"): r"\partial_t\phi = M\,\nabla^2\mu"})
        self.assertTrue(r["agree"])
        self.assertEqual(r["independent_reads"], 3)
        r = cr.compare_item(EQ, same, {"latex": r"\partial_t \phi = -M \nabla^2 \mu"}, {})
        self.assertFalse(r["agree"])

    def test_value_notation_and_not_found(self):
        r = cr.compare_item(VAL, None, {"parameter": {"quantity": "alpha", "value": "1\\times10^{-3}",
                                                      "unit_printed": "m/s"}}, {})
        self.assertTrue(r["agree"])
        r = cr.compare_item(VAL, None, {"not_found": True}, {})
        self.assertFalse(r["agree"])
        r = cr.compare_item(VAL, None, {"parameter": {"quantity": "alpha", "value": "1e-4", "unit_printed": "m/s"}}, {})
        self.assertFalse(r["agree"])

    def test_number_key(self):
        self.assertEqual(cr.number_key("(12a)"), "12a")
        self.assertIsNone(cr.number_key("(unnumbered)"))


if __name__ == "__main__":
    unittest.main()
