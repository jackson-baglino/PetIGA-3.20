"""Tests for check_calibration.py on planted faults in temporary copies of the files.

    venv_workflows/bin/python -m unittest discover -s scripts/literature_review/tests
"""
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import check_calibration as cc  # noqa: E402


class CheckCalibration(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        for f in os.listdir(cc.CAL):
            if f.endswith(".yaml"):
                shutil.copy(os.path.join(cc.CAL, f), self.tmp)
        self.saved = cc.CAL
        cc.CAL = self.tmp

    def tearDown(self):
        cc.CAL = self.saved
        shutil.rmtree(self.tmp)

    def write(self, name, text):
        with open(os.path.join(self.tmp, name), "w") as f:
            f.write(text)

    def problems(self, check):
        problems, counts = [], {}
        check(cc.papers(problems), problems, counts)
        return problems, counts

    def test_unquoted_colon_is_reported_with_its_line(self):
        self.write("planted_errors.yaml", "errors:\n  - paper: li_2009\n    note: a: b: c\n")
        problems, _ = self.problems(cc.check_planted)
        self.assertTrue(any("not valid YAML at line 3" in p for p in problems))

    def test_dropped_condition_needs_condition_page(self):
        self.write("planted_errors.yaml", "errors:\n  - paper: li_2009\n    page: 85\n"
                   "    kind: dropped_condition\n    correct_claim: x holds for small eps\n"
                   "    planted_claim: x holds\n    note: condition dropped\n")
        problems, counts = self.problems(cc.check_planted)
        self.assertTrue(any("needs condition_page" in p for p in problems))
        self.assertEqual(counts["planted errors"][0], 1)

    def test_printed_page_outside_the_paper(self):
        self.write("planted_errors.yaml", "errors:\n  - paper: li_2009\n    page: 12\n"
                   "    kind: wrong_unit\n    correct_claim: a\n    planted_claim: b\n    note: c\n")
        problems, _ = self.problems(cc.check_planted)
        self.assertTrue(any("printed page 12 is not in the paper" in p for p in problems))

    def test_unknown_paper(self):
        self.write("tret_questions.yaml", "questions:\n  - kind: exact\n    question: q?\n"
                   "    paper: nobody_2020\n    page: 3\n    answer: a\n")
        problems, _ = self.problems(cc.check_tret)
        self.assertTrue(any("'nobody_2020' is not in papers.yaml" in p for p in problems))

    def test_absent_question_must_leave_answer_empty(self):
        self.write("tret_questions.yaml", "questions:\n  - kind: absent\n    question: q?\n"
                   "    answer: something\n")
        problems, _ = self.problems(cc.check_tret)
        self.assertTrue(any("absent question leaves" in p for p in problems))

    def test_basis_must_be_a_listed_one(self):
        self.write("gold_statements.yaml", "statements:\n  - paper: li_2009\n    page: 85\n"
                   "    kind: result\n    statement: x\n    basis: measured\n"
                   "  - paper: li_2009\n    page: 86\n    kind: result\n    statement: y\n"
                   "    basis: cited\n    typo: z\n")
        problems, _ = self.problems(cc.check_gold)
        self.assertEqual(len(problems), 1)
        self.assertIn("basis 'measured'", problems[0])

    def test_value_must_be_one_number(self):
        self.write("gold_statements.yaml", "statements:\n  - paper: li_2009\n    page: 85\n"
                   "    kind: parameter\n    statement: x\n    value: '[1e-3, 1e-1]'\n    unit: '1'\n"
                   "  - paper: li_2009\n    page: 86\n    kind: parameter\n    statement: y\n"
                   "    value: 1.8e6\n    unit: J m^-3 K^-1\n")
        problems, _ = self.problems(cc.check_gold)
        self.assertEqual(len(problems), 1)
        self.assertIn("value must be one number", problems[0])

    def test_empty_slots_are_ignored(self):
        self.write("planted_errors.yaml", "errors:\n  - paper: li_2009\n    page:\n    kind:\n"
                   "    correct_claim:\n    planted_claim:\n    note:\n")
        problems, _ = self.problems(cc.check_planted)
        self.assertEqual(problems, [])


if __name__ == "__main__":
    unittest.main()
