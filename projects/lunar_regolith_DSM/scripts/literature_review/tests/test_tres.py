"""Tests for tres.py: the LaTeX canonicalizer must not merge different math (Codex finding 5)."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tres  # noqa: E402

SAME = [
    (r"$\displaystyle \rho_{\mathrm{eq}}(\phi,T)$", r"\rho_{eq}(\phi, T)"),
    (r"\left( 1 + \epsilon^{-1} \right)", r"(1+\epsilon^{-1})"),
    (r"\tau\frac{\partial\phi}{\partial t}=W^{2}\nabla^2\phi,", r"\tau \dfrac{\partial \phi}{\partial t} = W^2 \nabla^{2} \phi"),
    (r"{\rm d} x \, \quad y", r"dxy"),
    (r"a = b \tag{3}", r"a=b"),
    (r"F^{\rm tri} = {\rm exp}(x)", r"F^{tri}=\exp(x)"),
    (r"x_{i}+{y}", r"x_i+y"),
    ("W \u226a 10\u207b\u00b3 m", r"W \ll 10^{-3} m"),
    ("2 \u00d7 10^{-5}", r"2\times10^{-5}"),
]
DIFFERENT = [
    (r"\epsilon^{-1}", r"\epsilon^{1}"),
    (r"\mathbf{u}", r"u"),
    (r"\nabla\cdot(D\nabla u)", r"\nabla(D\nabla u)"),
    (r"x_{ij}", r"x_{ji}"),
    (r"\frac{1+\phi}{2}", r"\frac{1-\phi}{2}"),
    (r"\partial_t u", r"\partial_x u"),
    (r"(a+b)c", r"a+bc"),
    (r"F^{\rm tri}", r"F^{t}ri"),
    (r"\frac{a}{b}", r"\frac ab c"),
]


class Canon(unittest.TestCase):
    def test_formatting_variants_match(self):
        for a, b in SAME:
            self.assertEqual(tres.canon(a), tres.canon(b), (a, b))

    def test_different_math_differs(self):
        for a, b in DIFFERENT:
            self.assertNotEqual(tres.canon(a), tres.canon(b), (a, b))

    def test_score(self):
        ref = [{"label": "(1)", "latex": "a=b"}, {"label": "(2)", "latex": "c=d"}, {"label": "(3)", "latex": "e"}]
        got = [{"label": "1", "latex": "a = b"}, {"label": "(2)", "latex": "c=-d"}]
        self.assertEqual([r["result"] for r in tres.score(ref, got)], ["exact", "wrong", "missed"])


if __name__ == "__main__":
    unittest.main()
