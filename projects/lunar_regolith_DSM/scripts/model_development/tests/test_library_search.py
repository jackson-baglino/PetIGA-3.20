"""Tests for library_search.py on a synthetic library (no real paper text).

Run: venv_workflows/bin/python -m unittest discover -s scripts/model_development/tests -t .
"""
import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import library_search as ls  # noqa: E402


def make_library(root):
    p = root / "library" / "papers" / "smith_2001_ab12"
    pages = p / "content" / "copyhash0000aaaa" / "pages"
    pages.mkdir(parents=True)
    (pages / "p001.md").write_text(
        "Title of the paper\n\nThe vapor diffusivity D is 2.2e-5 m^2/s at 273 K.\n\n"
        "A second paragraph about grain boundaries.\n", encoding="utf-8")
    (pages / "p002.md").write_text(
        "The interface width W must be small compared with D beta.\n\n"
        "Kinetic coefficient beta = 1/alpha times a constant.\n", encoding="utf-8")
    (p / "content" / "copyhash0000aaaa" / "equations.jsonl").write_text(
        json.dumps({"latex": r"\partial_t \phi = M \nabla^2 \phi", "number": "3", "page_pdf": 2}) + "\n",
        encoding="utf-8")
    claims = p / "claims"
    claims.mkdir()
    (claims / "smith_2001_ab12.c001.v1.json").write_text(json.dumps(
        {"claim": {"text": "D value", "quote": "The vapor diffusivity D is 2.2e-5 m^2/s at 273 K.", "page_pdf": 1}}))
    # an older version that must be ignored in favor of v2
    (claims / "smith_2001_ab12.c002.v1.json").write_text(json.dumps(
        {"claim": {"quote": "old wording", "page_pdf": 2}}))
    (claims / "smith_2001_ab12.c002.v2.json").write_text(json.dumps(
        {"claim": {"quote": "The interface width W must be small compared with D beta.", "page_pdf": 2}}))
    q = root / "library" / "papers" / "jones_2010_cd34" / "content" / "h2" / "pages"
    q.mkdir(parents=True)
    (q / "p001.md").write_text("Unrelated text on sintering necks.\n", encoding="utf-8")
    return root / "library"


class LibrarySearchTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.library = make_library(root)
        self.index = root / "index" / "lib.sqlite"
        self.log = root / "index" / "q.jsonl"
        with redirect_stdout(io.StringIO()):
            ls.main(["--index", str(self.index), "build", "--library", str(self.library)])

    def tearDown(self):
        self.tmp.cleanup()

    def run_cli(self, argv):
        buf = io.StringIO()
        code = 0
        with redirect_stdout(buf):
            try:
                ls.main(["--index", str(self.index)] + argv)
            except SystemExit as e:
                code = e.code or 0
        return code, buf.getvalue()

    def test_counts(self):
        c = ls.summary_counts(self.index)
        self.assertEqual(c["papers"], 2)
        self.assertEqual(c["by_kind"]["page"]["chunks"], 6)
        self.assertEqual(c["by_kind"]["equation"]["chunks"], 1)
        self.assertEqual(c["by_kind"]["claim"]["chunks"], 2)  # c002 v1 superseded

    def test_search_finds_page_and_logs(self):
        code, out = self.run_cli(["search", "vapor diffusivity", "--kind", "page", "--log", str(self.log)])
        r = json.loads(out)
        self.assertFalse(r["abstain"])
        self.assertEqual(r["hits"][0]["page_pdf"], 1)
        self.assertEqual(r["hits"][0]["paper_id"], "smith_2001_ab12")
        self.assertIn("2.2e-5", r["hits"][0]["text"])  # verbatim
        logged = [json.loads(x) for x in self.log.read_text().splitlines()]
        self.assertEqual(logged[-1]["query"], "vapor diffusivity")

    def test_search_abstains(self):
        code, out = self.run_cli(["search", "zzzqqq", "--log", str(self.log)])
        self.assertTrue(json.loads(out)["abstain"])

    def test_search_syntax_safe(self):
        # FTS operators and quotes in the query must not raise
        code, out = self.run_cli(["search", 'beta" OR NEAR(x) -"', "--log", str(self.log)])
        self.assertEqual(code, 0)

    def test_quote_exact_and_whitespace(self):
        code, _ = self.run_cli(["quote", "The vapor  diffusivity D\nis 2.2e-5 m^2/s", "--paper",
                                "smith_2001_ab12", "--page", "1"])
        self.assertEqual(code, 0)

    def test_quote_wrong_page_or_changed_number_fails(self):
        code, _ = self.run_cli(["quote", "The vapor diffusivity D is 2.2e-5", "--paper",
                                "smith_2001_ab12", "--page", "2"])
        self.assertEqual(code, 1)
        code, _ = self.run_cli(["quote", "The vapor diffusivity D is 2.2e-6", "--paper",
                                "smith_2001_ab12", "--page", "1"])
        self.assertEqual(code, 1)

    def test_quote_spanning_paragraphs(self):
        code, _ = self.run_cli(["quote", "compared with D beta. Kinetic coefficient", "--paper",
                                "smith_2001_ab12", "--page", "2"])
        self.assertEqual(code, 0)

    def test_selfcheck(self):
        code, out = self.run_cli(["selfcheck"])
        t = json.loads(out)["total"]
        self.assertEqual(t, {"claims": 2, "found_on_page": 2, "top1_same_page": 2})

    def test_rebuild_is_deterministic(self):
        import sqlite3
        def dump(path):
            con = sqlite3.connect(path)
            rows = con.execute("SELECT * FROM chunks ORDER BY chunk_id").fetchall()
            files = con.execute("SELECT * FROM files ORDER BY path").fetchall()
            con.close()
            return rows, files
        first = dump(self.index)
        with redirect_stdout(io.StringIO()):
            ls.main(["--index", str(self.index), "build", "--library", str(self.library)])
        self.assertEqual(first, dump(self.index))


if __name__ == "__main__":
    unittest.main()
