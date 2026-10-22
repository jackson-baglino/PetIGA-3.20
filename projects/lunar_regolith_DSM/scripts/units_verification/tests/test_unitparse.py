"""Tests for unitparse. Run from the project root:
    venv_workflows/bin/python -m unittest discover -s scripts/units_verification/tests -v
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import unitparse as up  # noqa: E402
import uvlib  # noqa: E402

J_PER_M3K = up.parse_unit("J m^-3 K^-1")


class ParseUnit(unittest.TestCase):

    def test_equivalent_notations(self):
        for s in ["J m^-3 K^-1", "J m⁻³ K⁻¹", "J m−3 K−1", r"\mathrm{J\,m^{-3}\,K^{-1}}", "J/(m3 K)", "J/(m^3·K)",
                  "J m <sup>-3</sup> K <sup>-1</sup>", "J/m^3/K", "J/m3 K", r"\text{J}\,\text{m}^{-3}\,\text{K}^{-1}"]:
            with self.subTest(s=s):
                self.assertTrue(up.same_unit(up.parse_unit(s), J_PER_M3K))

    def test_celsius_is_a_kelvin_difference(self):
        self.assertTrue(up.same_unit(up.parse_unit("W/(m ◦C)"), up.parse_unit("W m^-1 K^-1")))
        self.assertTrue(up.same_unit(up.parse_unit("J/(Kg ◦C)"), up.parse_unit("J kg^-1 K^-1")))
        self.assertTrue(up.same_unit(up.parse_unit(r"\text{W}/(\text{m}\,{}^{\circ}\text{C})"),
                                     up.parse_unit("W m^-1 K^-1")))
        self.assertTrue(up.same_unit(up.parse_unit("°C"), up.parse_unit("K")))

    def test_dimensionless(self):
        for s in ["1 (dimensionless)", "(-)", "-", "–", "kg/kg", "1", "dimensionless", "m/m"]:
            with self.subTest(s=s):
                self.assertEqual(up.parse_unit(s).dims, up.ONE.dims)

    def test_derived_units(self):
        self.assertTrue(up.same_unit(up.parse_unit("N/m"), up.parse_unit("J m^-2")))
        self.assertTrue(up.same_unit(up.parse_unit("Pa"), up.parse_unit("N m^-2")))
        self.assertTrue(up.same_unit(up.parse_unit("W"), up.parse_unit("J/s")))
        self.assertTrue(up.same_unit(up.parse_unit("1/s"), up.parse_unit("s^-1")))
        self.assertTrue(up.same_unit(up.parse_unit("m2/s"), up.parse_unit("m^2 s^-1")))

    def test_prefixes_change_the_scale_not_the_dimension(self):
        mm, m = up.parse_unit("mm"), up.parse_unit("m")
        self.assertEqual(mm.dims, m.dims)
        self.assertFalse(up.same_unit(mm, m))
        self.assertAlmostEqual(up.parse_unit("μm").scale, 1e-6)
        self.assertAlmostEqual(up.parse_unit("µm").scale, 1e-6)   # micro sign, not Greek mu
        self.assertAlmostEqual(up.parse_unit("kg").scale, 1.0)
        self.assertAlmostEqual(up.parse_unit("g").scale, 1e-3)

    def test_planted_error_kinds_are_detected(self):
        """The C0 error kinds: flipped exponent sign, wrong prefix, swapped unit."""
        pairs = [("s m^-1", "s m^1"), ("m^2 s^-1", "m^2 s^1"), ("J m^-3", "kJ m^-3"), ("m", "mm"),
                 ("W m^-1 K^-1", "J m^-1 K^-1"), ("kg m^-3", "kg m^-2"), ("J K^-1", "J kg^-1")]
        for good, bad in pairs:
            with self.subTest(good=good, bad=bad):
                self.assertFalse(up.same_unit(up.parse_unit(good), up.parse_unit(bad)))

    def test_never_guesses(self):
        for s in ["furlong", "J m^-3 K^-", "m)(", "J/(m", "us", "h", "C"]:
            with self.subTest(s=s):
                with self.assertRaises(up.UnitError):
                    up.parse_unit(s)
        with self.assertRaises(up.UnitError):
            up.parse_unit(None)

    def test_every_unit_in_the_table_parses(self):
        for row in uvlib.read_table():
            with self.subTest(symbol=row["symbol"], unit=row["unit"]):
                up.parse_unit(row["unit"])


class ParseNumbers(unittest.TestCase):

    def test_forms(self):
        cases = {
            "1.8 × 10⁶, 1.4×10^{3} J m^{-3} K^{-1}": [1.8e6, 1.4e3],
            "−10 °C": [-10.0],
            "2.178 × 10−5 m2 s−1": [2.178e-5],          # the text layer's form
            r"$10^{-3}$ to $10^{-2}$ kg m$^{-3}$": [1e-3, 1e-2],
            "286.9, 461.5 J Kg <sup>-1</sup> K <sup>-1</sup>": [286.9, 461.5],
            r"$5.5 \times 10^{5} \text{ s m}^{-1}$": [5.5e5],
            "": [],
            None: [],
        }
        for text, want in cases.items():
            with self.subTest(text=text):
                got = up.parse_numbers(text)
                self.assertEqual(len(got), len(want), got)
                for g, w in zip(got, want):
                    self.assertAlmostEqual(g / w if w else g, 1.0 if w else 0.0, places=12)

    def test_exponents_and_subscripts_are_not_numbers(self):
        self.assertEqual(up.parse_numbers("T0 = 261 K, T7 = 260 K"), [261.0, 260.0])
        self.assertEqual(up.parse_numbers("kg m^-3"), [])
        self.assertEqual(up.parse_numbers("m2 s-1"), [])

    def test_containment(self):
        self.assertTrue(up.numbers_contained([1.8e6], [1.8e6, 1.4e3]))
        self.assertFalse(up.numbers_contained([1.8e-6], [1.8e6, 1.4e3]))
        self.assertTrue(up.numbers_contained([], [1.0]))


class Symbols(unittest.TestCase):

    def test_flatten(self):
        self.assertEqual(up.flatten(r"$\rho_{i}$"), "ρ_i")
        self.assertEqual(up.flatten(r"C_{\rm i}"), "C_i")
        self.assertEqual(up.flatten(r"D_v^{\mathrm{STP}}"), "D_v^STP")
        self.assertEqual(up.flatten(r"x = 1 \tag{3}"), "x = 1")
        self.assertEqual(up.flatten("10<sup>28</sup>"), "10^28")

    def test_symbol_matching(self):
        rx = up.symbol_regex("ρi")
        for s in ["ρ_i", "ρi", "ρ_ i", "(ρ_i)"]:
            self.assertTrue(rx.search(s), s)
        for s in ["ρ_is", "ρ_a", "xρ_i"]:
            self.assertFalse(rx.search(s), s)
        self.assertTrue(up.symbol_regex("vn,j").search("v_n,j"))
        self.assertTrue(up.symbol_regex("vn,j").search("v_n, j"))
        self.assertFalse(up.symbol_regex("β").search("β_j"), "a bare symbol must not match its subscripted form")
        self.assertTrue(up.symbol_regex("β").search("β = 5"))
        self.assertFalse(up.symbol_regex("Dv").search("D_v0"))

    def test_printed_key(self):
        self.assertEqual(up.printed_key("ρ_i"), up.printed_key("ρi"))
        self.assertEqual(up.printed_key(r"\rho_{i}"), "ρi")
        self.assertEqual(up.printed_key("v_{n,j}"), up.printed_key("vn,j"))
        self.assertEqual(up.printed_key("ϵ"), up.printed_key("ε"))
        self.assertNotEqual(up.printed_key("ρ_v"), up.printed_key("ρ_vs"))
        self.assertNotEqual(up.printed_key("D_STP"), up.printed_key("D_v^STP"))
        self.assertEqual(up.printed_key(None), "")

    def test_short_symbols_count_only_in_math_or_the_symbol_column(self):
        flat = up.flatten_markdown("a is the spacing, and $a$ = 3 m\n| a | spacing | 3 m |\n")
        hits = up.find_symbol(flat, "a")
        self.assertEqual(len(hits), 2)          # the $a$ and the table's first cell, not the article "a"

    def test_unit_runs(self):
        runs = up.unit_runs("= 1.8 × 10^6, 1.4 × 10^3 J m ^-3 K ^-1 , e.g. us and x_s")
        self.assertEqual(len(runs), 1)
        self.assertTrue(up.same_unit(runs[0][3], J_PER_M3K))
        runs = up.unit_runs("(m/s) | 2 × 10^-7")
        self.assertEqual(runs[0][2], "(m/s)")
        self.assertEqual(up.unit_runs("in simulations β = 5.5 s m^-1)")[0][2], "s m^-1")


if __name__ == "__main__":
    unittest.main()
