"""Offline tests for add_paper.py (no network, no MinerU, no library writes)."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import add_paper as ap  # noqa: E402


class AddPaperTest(unittest.TestCase):
    def test_default_paper_id(self):
        meta = {"doi": "10.1103/PhysRevE.84.046707", "year": 2011,
                "authors": [{"family": "Nicoli", "given": "M."}]}
        pid = ap.default_paper_id(meta)
        self.assertRegex(pid, r"^nicoli_2011_[0-9a-f]{4}$")
        # case of the DOI does not change the id
        meta2 = dict(meta, doi=meta["doi"].lower())
        self.assertEqual(pid, ap.default_paper_id(meta2))

    def test_default_paper_id_accents(self):
        meta = {"doi": "10.1/x", "year": 2015, "authors": [{"family": "Lervåg", "given": "K."}]}
        self.assertTrue(ap.default_paper_id(meta).startswith("lervag_2015_"))

    def test_pdf_from_incoming(self):
        with tempfile.TemporaryDirectory() as tmp:
            review = Path(tmp) / "model-line"
            (review / "zotero_ABC").mkdir(parents=True)
            (review / "zotero_pull.json").write_text(json.dumps({"copies": [
                {"doi": "10.1103/physreve.84.046707", "path": "zotero_ABC/x.pdf"}]}))
            with mock.patch.object(ap, "INCOMING", Path(tmp)):
                p = ap.pdf_from_incoming("model-line", "10.1103/PhysRevE.84.046707")
                self.assertEqual(p, review / "zotero_ABC" / "x.pdf")
                with self.assertRaises(SystemExit):
                    ap.pdf_from_incoming("model-line", "10.9/none")


if __name__ == "__main__":
    unittest.main()
