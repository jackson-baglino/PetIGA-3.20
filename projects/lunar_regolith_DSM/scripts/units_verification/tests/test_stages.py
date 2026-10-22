"""Tests for the stage scripts (S0, S1, S1b helpers, S3) on synthetic inputs. Run from the project root:
    venv_workflows/bin/python -m unittest discover -s scripts/units_verification/tests -v
"""
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import compare_units  # noqa: E402
import locate_sources  # noqa: E402
import marker_extract  # noqa: E402
import render_pages  # noqa: E402
import uvlib  # noqa: E402

TABLE = """# t

| symbol | meaning | unit | source | page | quote (text layer) | verified |
|---|---|---|---|---|---|---|
| `M0` | mobility | m s^-1 | dimensional-analysis |  |  | code |
| `Ci` | heat capacity per volume of ice (printed as Ci) | J m^-3 K^-1 | src_a | 4 | q | v |
| `beta` | kinetic coefficient (printed as β) | s m^-1 | src_a | 4 | q | v |
| `x_s` | humidity ratio (printed as xs) | 1 (dimensionless) | src_a | 4 | q | v |
| `L` | pore size (printed as L) | m | src_a | 9 | q | v |
"""
MARKER_P4 = r"""| Symbol | Meaning | Value | Ref |
|---|---|---|---|
| $C_i, C_a$ | heat capacity of ice, air per volume | $1.8 \times 10^{6}$ , $1.4 \times 10^{3}$ J m <sup>-3</sup> K <sup>-1</sup> | [30] |
| β | interface kinetic coefficient | in simulations $\beta = 5.5 \times 10^5 \text{ s m}^{-1}$ | |
| $x_s$ | humidity ratio | kg/kg | |
"""
MARKER_P9 = "We take an estimated pore size of $L = 10^{-4}$ m.\n"


def read(symbol, found=True, unit=None, value=None, printed=None):
    return {"symbol": symbol, "found": found, "printed_as": printed, "meaning": "m", "value_as_printed": value,
            "unit_as_printed": unit, "location": "Table I"}


class TableReading(unittest.TestCase):

    def test_real_table(self):
        rows = uvlib.read_table()
        lit = uvlib.literature_rows(rows)
        self.assertEqual(len(rows), 62)
        self.assertEqual(len(lit), 60)
        self.assertTrue(all(r["page"] for r in lit))
        self.assertTrue(all(r["printed"] for r in lit))
        self.assertEqual({r["source"] for r in lit}, set(uvlib.LITERATURE_SOURCES))

    def test_synthetic_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "units.md"
            p.write_text(TABLE)
            rows = uvlib.read_table(p)
            self.assertEqual([r["symbol"] for r in rows], ["M0", "Ci", "beta", "x_s", "L"])
            self.assertEqual(rows[2]["printed"], "β")
            self.assertIsNone(rows[0]["page"])
            p.write_text(TABLE + TABLE.splitlines()[-1] + "\n")
            with self.assertRaises(uvlib.StageError):
                uvlib.read_table(p)


class Compare(unittest.TestCase):
    """S3 verdicts on a synthetic Marker page; includes the four C0 error kinds."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        (self.tmp / "units.md").write_text(TABLE)
        self.mdir = self.tmp / "marker"
        self.mdir.mkdir()
        (self.mdir / "src_a_p4.md").write_text(MARKER_P4)
        (self.mdir / "src_a_p9.md").write_text(MARKER_P9)
        self.rows = {r["symbol"]: r for r in uvlib.read_table(self.tmp / "units.md")}
        self.index = {"Ci": {"pages": [4]}, "beta": {"pages": [4]}, "x_s": {"pages": [4]}, "L": {"pages": [9]}}
        self.good = {("src_a", 4): {"Ci": read("Ci", unit="J m⁻³ K⁻¹", value="1.8 × 10⁶"),
                                    "beta": read("beta", unit="s m⁻¹", value="5.5 × 10⁵"),
                                    "x_s": read("x_s", unit="kg/kg")},
                     ("src_a", 9): {"L": read("L", unit="m", value="10⁻⁴")}}

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def verdict(self, sym, reads=None, row=None):
        return compare_units.compare_row(row or self.rows[sym], self.good if reads is None else reads, self.mdir, self.index)

    def test_all_agree(self):
        for sym in ("Ci", "beta", "x_s", "L"):
            with self.subTest(sym=sym):
                r = self.verdict(sym)
                self.assertEqual(r["verdict"], "AGREE", r["reasons"])

    def test_flipped_exponent_sign_in_table(self):
        r = self.verdict("beta", row=dict(self.rows["beta"], unit="s m^1"))
        self.assertEqual(r["verdict"], "DISAGREE")
        self.assertTrue(any(x.startswith("unit:") for x in r["reasons"]))

    def test_wrong_prefix_in_table(self):
        r = self.verdict("L", row=dict(self.rows["L"], unit="mm"))
        self.assertEqual(r["verdict"], "DISAGREE")

    def test_swapped_unit_in_table(self):
        r = self.verdict("Ci", row=dict(self.rows["Ci"], unit="W m^-3 K^-1"))
        self.assertEqual(r["verdict"], "DISAGREE")

    def test_wrong_page_in_table(self):
        reads = {("src_a", 4): dict(self.good[("src_a", 4)], L=read("L", found=False))}  # p9 not read yet
        r = self.verdict("L", reads=reads, row=dict(self.rows["L"], page=4))
        self.assertEqual(r["verdict"], "NOT_FOUND")
        self.assertEqual(r["reread_pages"], [9])
        reads[("src_a", 9)] = {"L": read("L", unit="m", value="10⁻⁴")}
        r = self.verdict("L", reads=reads, row=dict(self.rows["L"], page=4))
        self.assertEqual(r["verdict"], "DISAGREE")
        self.assertIn("page: found on p9, not on the cited p4", r["reasons"])

    def test_read_and_table_agree_but_marker_does_not(self):
        r = self.verdict("Ci", row=dict(self.rows["Ci"], unit="J m^-3"),
                         reads={("src_a", 4): {"Ci": read("Ci", unit="J m^-3", value="1.8 × 10⁶")}})
        self.assertEqual(r["verdict"], "DISAGREE")
        self.assertTrue(any(x.startswith("marker: unit") for x in r["reasons"]))

    def test_value_not_in_marker(self):
        reads = {("src_a", 4): {"Ci": read("Ci", unit="J m^-3 K^-1", value="1.8 × 10⁻⁶")}}
        r = self.verdict("Ci", reads=reads)
        self.assertEqual(r["verdict"], "DISAGREE")
        self.assertIn("marker: printed value not found near the symbol", r["reasons"])

    def test_missing_read_and_unparsable_read_unit(self):
        self.assertEqual(self.verdict("Ci", reads={})["verdict"], "MISSING")
        reads = {("src_a", 4): {"Ci": read("Ci", unit="joules per something")}}
        r = self.verdict("Ci", reads=reads)
        self.assertEqual(r["verdict"], "DISAGREE")
        self.assertTrue(any(x.startswith("read_unit_unparsed") for x in r["reasons"]))

    def test_marker_only_mode(self):
        r = compare_units.compare_row(self.rows["Ci"], {}, self.mdir, self.index, marker_only=True)
        self.assertEqual(r["verdict"], "MARKER_AGREES", r["reasons"])


class MarkerSplit(unittest.TestCase):

    def test_split_pages(self):
        sep = "-" * 48
        md = f"{{0}}{sep}\n\npage one\n\n{{1}}{sep}\n\n\n{{2}}{sep}\n\npage three\n"
        pages = marker_extract.split_pages(md)
        self.assertEqual(pages[1].strip(), "page one")
        self.assertEqual(pages[2].strip(), "")
        self.assertEqual(pages[3].strip(), "page three")


class RenderPages(unittest.TestCase):

    def test_wanted_pages(self):
        rows = [{"source": "a", "page": 1}, {"source": "a", "page": 5}, {"source": "b", "page": None}]
        self.assertEqual(render_pages.wanted_pages(rows, {"a": 5, "b": 3}), {("a", 1), ("a", 2), ("a", 4), ("a", 5)})


@unittest.skipUnless(shutil.which("mutool"), "mutool not installed")
class LocateSources(unittest.TestCase):
    """S0 on generated PDFs: DOI on page 1, title check, identical-render copies, differing copies."""

    def make_pdf(self, path, text, extra_font=False):
        spec = Path(self.tmp) / "content.txt"
        fonts = "%%Font F1 Helvetica\n" + ("%%Font F2 Times-Roman\n" if extra_font else "")
        spec.write_text(f"%%MediaBox 0 0 300 200\n{fonts}BT /F1 9 Tf 10 150 Td ({text}) Tj ET\n")
        path.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["mutool", "create", "-o", str(path), str(spec)], check=True, capture_output=True)

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.lit = Path(self.tmp) / "lit"
        self.spec = {"p": {"doi": "10.1/abc", "title": "A Paper Title"},
                     "si": {"doi": None, "title": "Supporting Information: Something", "prefer": "sub/si.pdf"}}

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_finds_and_deduplicates(self):
        self.make_pdf(self.lit / "paper.pdf", "A Paper Title doi 10.1/abc")
        self.make_pdf(self.lit / "other.pdf", "Unrelated 10.2/xyz")
        self.make_pdf(self.lit / "sub/si.pdf", "Supporting Information: Something")
        self.make_pdf(self.lit / "si_copy.pdf", "Supporting Information: Something", extra_font=True)
        self.assertNotEqual(uvlib.sha256(self.lit / "sub/si.pdf"), uvlib.sha256(self.lit / "si_copy.pdf"))
        result, problems, n, _ = locate_sources.locate(self.lit, self.spec)
        self.assertEqual(problems, [])
        self.assertEqual(n, 4)
        self.assertTrue(result["p"]["path"].endswith("paper.pdf"))
        self.assertTrue(result["si"]["path"].endswith("sub/si.pdf"))
        self.assertEqual(len(result["si"]["identical_copies"]), 1)

    def test_stops_on_missing_or_ambiguous(self):
        self.make_pdf(self.lit / "sub/si.pdf", "Supporting Information: Something")
        self.make_pdf(self.lit / "si_v2.pdf", "Supporting Information: Something, revised")
        result, problems, _, _ = locate_sources.locate(self.lit, self.spec)
        self.assertNotIn("p", result)
        self.assertNotIn("si", result)
        self.assertTrue(any("render differently" in p for p in problems))
        self.assertTrue(any(p.startswith("p: no PDF") for p in problems))

    def test_doi_without_title_is_a_problem(self):
        self.make_pdf(self.lit / "paper.pdf", "doi 10.1/abc only")
        _, problems, _, _ = locate_sources.locate(self.lit, {"p": self.spec["p"]})
        self.assertTrue(any("does not show the title" in p for p in problems))


if __name__ == "__main__":
    unittest.main()
