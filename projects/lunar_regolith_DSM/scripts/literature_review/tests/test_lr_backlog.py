"""Tests for lr_backlog.py."""
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lr_backlog  # noqa: E402


class Backlog(unittest.TestCase):
    def test_add_goes_under_open(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "BACKLOG.md"
            p.write_text("# B\n\nintro\n\n## Open\n\n## Decided\n\n- old\n")
            lr_backlog.add("E4 run x", "Split  the checker.", "Too long.", p)
            lr_backlog.add("session", "Second idea.", "Because.", p)
            items = lr_backlog.open_items(p)
            self.assertEqual(len(items), 2)
            self.assertIn("Split the checker. Why: Too long.", items[0])
            self.assertTrue(p.read_text().rstrip().endswith("- old"))


if __name__ == "__main__":
    unittest.main()
