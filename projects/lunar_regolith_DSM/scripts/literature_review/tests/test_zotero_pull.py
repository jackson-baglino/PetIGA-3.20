"""Tests for zotero_pull.py on a fake local API: no network, no real Zotero.

    venv_workflows/bin/python -m unittest discover -s scripts/literature_review/tests
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import zotero_pull as zp  # noqa: E402


def item(key, doi=None, extra=None, kind="journalArticle"):
    d = {"itemType": kind}
    if doi:
        d["DOI"] = doi
    if extra:
        d["extra"] = extra
    return {"key": key, "data": d}


def att(key, link="imported_url", ctype="application/pdf"):
    return {"key": key, "data": {"itemType": "attachment", "linkMode": link, "contentType": ctype}}


class Fake(zp.Local):
    def __init__(self, lists, files):
        self.lists, self.files = lists, files

    def get_all(self, path, params=None):
        return self.lists.get(path, [])

    def file_path(self, a):
        return self.files.get(a["key"])


class Pull(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.store = os.path.join(self.tmp, "zotero_storage")
        self.out = os.path.join(self.tmp, "incoming")
        os.makedirs(self.store)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def pdf(self, name, body=b"%PDF-1.7 fake"):
        p = os.path.join(self.store, name)
        with open(p, "wb") as f:
            f.write(body)
        return p

    def fake(self):
        lists = {
            "/items/top": [item("I1", doi="https://doi.org/10.1103/PhysRevE.72.011602"),
                           item("I2", extra="arXiv: 1\nDOI: 10.48550/arXiv.2407.21213"),
                           item("I3", doi="10.1002/nme.5961"),
                           item("I4", doi="10.9999/unrelated")],
            "/items/I1/children": [att("A1"), {"key": "N1", "data": {"itemType": "note"}}],
            "/items/I2/children": [att("A2", link="imported_file")],
            "/items/I3/children": [att("A3", link="linked_url")],
            "/items/I4/children": [att("A4")],
            "/collections": [{"key": "C1", "data": {"name": "lunar_regolith_DSM"}},
                             {"key": "C2", "data": {"name": "B1", "parentCollection": "C1"}}],
            "/collections/C2/items/top": [item("I1", doi="10.1103/physreve.72.011602")],
        }
        files = {"A1": self.pdf("a1.pdf"), "A2": self.pdf("a2.pdf", b"%PDF-1.5 other"),
                 "A4": self.pdf("a4.pdf")}
        return Fake(lists, files)

    def test_dois_are_found_in_a_fetch_list(self):
        text = "- [Folch](https://doi.org/10.1103/PhysRevE.72.011602)\n- doi:10.1002/nme.5961.\n"
        self.assertEqual(zp.dois_in(text), ["10.1002/nme.5961", "10.1103/physreve.72.011602"])

    def test_matches_by_doi_copies_and_records_hashes(self):
        wanted = {"10.1103/physreve.72.011602", "10.48550/arxiv.2407.21213", "10.1002/nme.5961",
                  "10.4310/cms.2015.v13.n6.a6"}
        s = zp.pull(self.fake(), wanted, "b1", self.out, today="2026-10-29")
        self.assertEqual(sorted(s["copied"]), ["10.1103/physreve.72.011602", "10.48550/arxiv.2407.21213"])
        self.assertEqual(s["no_pdf"], ["10.1002/nme.5961"])                    # a link, not a stored file
        self.assertEqual(s["not_in_zotero"], ["10.4310/cms.2015.v13.n6.a6"])
        with open(os.path.join(self.out, "b1", "zotero_pull.json")) as f:
            man = json.load(f)
        self.assertEqual(len(man["copies"]), 2)
        for c in man["copies"]:
            dest = os.path.join(self.out, "b1", c["path"])
            self.assertEqual(zp.sha256(dest), c["sha256"])
            self.assertTrue(c["path"].startswith("zotero_A"))
        # the unrelated item was not copied
        self.assertFalse(any(c["attachment_key"] == "A4" for c in man["copies"]))

    def test_rerun_copies_nothing(self):
        wanted = {"10.1103/physreve.72.011602"}
        zp.pull(self.fake(), wanted, "b1", self.out)
        s = zp.pull(self.fake(), wanted, "b1", self.out)
        self.assertEqual(s["copied"], [])
        self.assertEqual(s["already"], ["10.1103/physreve.72.011602"])

    def test_same_file_under_a_new_attachment_is_not_copied_again(self):
        z = self.fake()
        zp.pull(z, {"10.1103/physreve.72.011602"}, "b1", self.out)
        z.lists["/items/top"].append(item("I5", doi="10.1103/physreve.72.011602"))
        z.lists["/items/I5/children"] = [att("A5")]
        z.files["A5"] = z.files["A1"]
        s = zp.pull(z, {"10.1103/physreve.72.011602"}, "b1", self.out)
        self.assertEqual(s["copied"], [])

    def test_a_file_that_is_not_a_pdf_is_skipped(self):
        z = self.fake()
        z.files["A1"] = self.pdf("bad.pdf", b"<html>login page</html>")
        s = zp.pull(z, {"10.1103/physreve.72.011602"}, "b1", self.out)
        self.assertEqual(s["not_a_pdf"], ["10.1103/physreve.72.011602"])
        self.assertFalse(os.path.exists(os.path.join(self.out, "b1", "zotero_pull.json")))

    def test_dry_run_writes_nothing(self):
        s = zp.pull(self.fake(), {"10.1103/physreve.72.011602"}, "b1", self.out, dry_run=True)
        self.assertEqual(s["copied"], ["10.1103/physreve.72.011602"])
        self.assertFalse(os.path.exists(self.out))

    def test_collection_by_full_path(self):
        z = self.fake()
        self.assertEqual(zp.find_collection(z, "lunar_regolith_DSM/B1"), "C2")
        s = zp.pull(z, {"10.48550/arxiv.2407.21213"}, "b1", self.out, collection="lunar_regolith_DSM/B1")
        self.assertEqual(s["not_in_zotero"], ["10.48550/arxiv.2407.21213"])

    def test_only_the_local_api_is_allowed(self):
        with self.assertRaises(SystemExit):
            zp.Local("https://api.zotero.org/users/1")


if __name__ == "__main__":
    unittest.main()
