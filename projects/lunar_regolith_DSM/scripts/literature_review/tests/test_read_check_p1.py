"""read_check.py takes its page images from P1's content, and refuses to read without it."""
import hashlib
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import read_check as rc  # noqa: E402


class FromP1(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.pdf = self.root / "x.pdf"
        self.pdf.write_bytes(b"%PDF-fake")
        self.run = self.root / "run"
        h = hashlib.sha256(b"%PDF-fake").hexdigest()[:16]
        self.images = self.run / "untrusted" / "content" / "pa" / h / "images"

    def tearDown(self):
        self.tmp.cleanup()

    def test_refuses_without_p1(self):
        with self.assertRaises(SystemExit) as e:
            rc.main(["--run", str(self.run), "--paper-id", "pa", "--pdf", str(self.pdf), "--stages", "E1"])
        self.assertIn("no P1 content", str(e.exception))

    def test_copies_and_scales_p1_images(self):
        self.images.mkdir(parents=True)
        Image.new("L", (2550, 3300), 255).save(self.images / "p001.png")
        src = rc.p1_images(self.run, "pa", self.pdf)
        self.assertEqual(src, self.images)
        out = self.root / "out"
        with mock.patch.object(rc.subprocess, "run") as run:        # only the text layer is made by mutool
            pages = rc.render(self.pdf, out, src)
        self.assertEqual([p.name for p in pages], ["p001.png"])
        self.assertEqual(max(Image.open(pages[0]).size), rc.EDGE)
        self.assertTrue(all("txt" in c.args[0] for c in run.call_args_list))
        self.assertEqual(Image.open(self.images / "p001.png").size, (2550, 3300))   # P1's image untouched


if __name__ == "__main__":
    unittest.main()
