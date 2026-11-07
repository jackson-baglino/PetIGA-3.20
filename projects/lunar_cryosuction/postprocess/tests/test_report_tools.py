"""Tests for the batch-report helpers: postprocess/run_health.py and
scripts/md2pdf.py (Markdown conversion only; the Chrome print step is not run).

    venv_cryosuction/bin/python -m unittest postprocess.tests.test_report_tools
"""
import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
sys.path.insert(0, HERE)
import fixtures as fx  # noqa: E402
import md2pdf  # noqa: E402
import run_health  # noqa: E402


def write_logs(run, n_steps, bounds=(-1e-19, 1.0), rollbacks=0, cfl=0, ok=True):
    with open(os.path.join(run, "outp.txt"), "w") as fh:
        for _ in range(n_steps):
            fh.write(f"  BOUNDS: phi_ice [{bounds[0]:.4e}, {bounds[1]:.4e}]  phi_air [0, 1]\n")
        for k in range(rollbacks):
            fh.write(f"[WARN] Phase field out of bounds at step {k} — deferring rollback\n")
        for k in range(cfl):
            fh.write(f"[WARN] Interface-CFL violated at step {k}: max |dphi| 1 > 0.5, rolling back, dt 1 -> 0.5\n")
        fh.write("RESOLUTION SUMMARY (after the 1.00e+02 s start-up window; max Pi inside it 9, not judged)\n")
    with open(os.path.join(run, "job.o123"), "w") as fh:
        fh.write("Wall time:               1h 02m 03s\nMPI ranks:               4\n")
        if ok:
            fh.write("  ✅ Workflow completed successfully\n")


class TestRunHealth(unittest.TestCase):
    def make(self, d, t_final=1000.0, t_end=1000.0, **kw):
        run = fx.make_run(d, experiment={"-t_final": str(t_final), "-dtmax": "100"})
        open(os.path.join(run, "igasol.dat"), "w").close()
        t = np.linspace(0, t_end, 11)
        fx.write_ssa(run, t, 1e-8 * np.ones_like(t), rhov=1e-12)
        with open(os.path.join(run, "solver_evo.dat"), "w") as fh:
            fh.write("# step t dt newton_its krylov_its krylov_per_newton rejections\n0 0 1 0 0 0 0\n")
            for k in range(1, 11):
                dt = 100.0 if k > 2 else 10.0
                fh.write(f"{k} {t[k]} {dt} {1 if k != 5 else 3} 20 20.0 {1 if k > 7 else 0}\n")
        fx.write_resolution(fx_run := run, t[1:], 1e-12 * np.r_[1e3, np.ones(9)])
        write_logs(fx_run, 10, **kw)
        return run

    def test_counts(self):
        with tempfile.TemporaryDirectory() as d:
            run = self.make(d, rollbacks=2, cfl=1)
            r = run_health.health(run)
            self.assertTrue(r["completed"])
            self.assertEqual(r["steps"], 10)
            self.assertAlmostEqual(r["newton_mean"], 1.2)
            self.assertEqual(r["newton_max"], 3)
            self.assertEqual(r["rejected"], 1)
            self.assertEqual((r["rb_band"], r["rb_cfl"]), (2, 1))
            self.assertAlmostEqual(r["at_dtmax"], 0.8)
            self.assertEqual(r["ranks"], 4)
            self.assertEqual(r["wall"], "1h 02m 03s")
            # the 1e-9 m/s spike is at t = 100 s, inside the 100 s window: not judged
            self.assertLess(r["Pi_max"], 0.5)

    def test_short_run_is_not_completed(self):
        with tempfile.TemporaryDirectory() as d:
            run = self.make(d, t_final=2000.0)
            self.assertFalse(run_health.health(run)["completed"])

    def test_failed_job_is_not_completed(self):
        with tempfile.TemporaryDirectory() as d:
            run = self.make(d, ok=False)
            self.assertFalse(run_health.health(run)["completed"])


class TestMd2Pdf(unittest.TestCase):
    def test_nested_lists_keep_numbering(self):
        html = md2pdf.convert("1. one\n2. two\n   - a\n   - b\n3. three\n\n4. four\n")
        self.assertEqual(html.count("<ol>"), 1)
        self.assertEqual(html.count("<ul>"), 1)
        self.assertEqual(html.count("<li>"), 6)
        self.assertLess(html.index("<ul>"), html.index("three"))

    def test_table_and_inline(self):
        html = md2pdf.convert("| ID | text |\n|---|---|\n| 3A-1 | **bold** and `code` |\n")
        self.assertIn('<td style="white-space:nowrap">3A-1</td>', html)
        self.assertIn("<strong>bold</strong>", html)
        self.assertIn("<code>code</code>", html)

    def test_image_printed_at_its_dpi(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        with tempfile.TemporaryDirectory() as d:
            fig = plt.figure(figsize=(4, 2))
            fig.savefig(os.path.join(d, "f.png"), dpi=150)
            plt.close(fig)
            html = md2pdf.convert("![cap](f.png)\n", base=d)
            self.assertIn("width:4.000in;height:2.000in", html)


if __name__ == "__main__":
    unittest.main()
