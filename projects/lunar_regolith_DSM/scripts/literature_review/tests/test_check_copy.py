"""Tests for check_copy.py (A2): copy kind from first-page signals."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import check_copy as cc  # noqa: E402

DOI = "10.1103/PhysRevE.79.031502"


class Kind(unittest.TestCase):
    def test_kinds(self):
        self.assertEqual(cc.signals(f"Phys Rev E ... DOI: {DOI} ©2009 The American Physical Society", DOI)[1],
                         "publisher")
        self.assertEqual(cc.signals("arXiv:1105.1670v1 [cond-mat] 9 May 2011", DOI)[1], "arxiv")
        self.assertEqual(cc.signals(f"Author manuscript; doi {DOI}", DOI)[1], "accepted_manuscript")
        self.assertEqual(cc.signals("A paper with no stamps at all.", DOI)[1], "unknown")
        s, k = cc.signals("Published by International Press", DOI)
        self.assertEqual((k, s["doi_found"], s["publisher_line"]), ("unknown", False, True))

    def test_doi_match_ignores_case_and_spaces(self):
        self.assertTrue(cc.signals("doi:10.1103/physreve.79. 031502", DOI)[0]["doi_found"])


if __name__ == "__main__":
    unittest.main()
