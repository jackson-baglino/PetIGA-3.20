"""Tests for the agent driver, S2/S4 helpers and C0 scoring (no agent is called). Run from the project root:
    venv_workflows/bin/python -m unittest discover -s scripts/units_verification/tests -v
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import agents  # noqa: E402
import calibrate  # noqa: E402
import check_rows  # noqa: E402
import read_pages  # noqa: E402
import unitparse as up  # noqa: E402
import uvlib  # noqa: E402


class AgentDefinitions(unittest.TestCase):

    def test_pinned_models_and_names(self):
        """CLAUDE.md: full model ID in the frontmatter; the name prefix matches the model family."""
        prefix = {"claude-sonnet-": "morty-", "claude-opus-": "beth-", "claude-haiku-": "jerry-", "claude-fable-": "rick-"}
        for name in ("morty-unit-reader", "beth-unit-checker"):
            a = agents.load_agent(name)
            with self.subTest(name=name):
                self.assertRegex(a["model"], r"^claude-[a-z]+-\d+-\d+$")
                fam = next(k for k in prefix if a["model"].startswith(k))
                self.assertTrue(name.startswith(prefix[fam]))
                self.assertEqual(a["tools"], ["Read"])


class Validate(unittest.TestCase):

    def test_schema(self):
        self.assertEqual(agents.validate({"reads": []}, read_pages.SCHEMA), [])
        good = {"symbol": "a", "found": False, "printed_as": None, "meaning": None, "value_as_printed": None,
                "unit_as_printed": None, "location": None}
        self.assertEqual(agents.validate({"reads": [good]}, read_pages.SCHEMA), [])
        self.assertTrue(agents.validate({"reads": [dict(good, found="yes")]}, read_pages.SCHEMA))
        self.assertTrue(agents.validate({"reads": [dict(good, extra=1)]}, read_pages.SCHEMA))
        bad = dict(good)
        del bad["unit_as_printed"]
        self.assertTrue(agents.validate({"reads": [bad]}, read_pages.SCHEMA))
        v = {"symbol": "a", "verdict": "maybe", "printed_as": None, "unit": None, "value": None, "meaning": None,
             "page": 4, "location": None, "reason": ""}
        self.assertTrue(agents.validate(v, check_rows.S4_SCHEMA))
        self.assertEqual(agents.validate(dict(v, verdict="cannot_tell"), check_rows.S4_SCHEMA), [])
        self.assertTrue(agents.validate(dict(v, verdict="cannot_tell", page=True), check_rows.S4_SCHEMA))


class Audit(unittest.TestCase):

    def events(self, tool_uses, model="claude-sonnet-5-5"):
        ev = [{"type": "system", "subtype": "init", "model": model}]
        ev += [{"type": "assistant", "message": {"content": [{"type": "tool_use", "name": n, "input": i}]}}
               for n, i in tool_uses]
        ev.append({"type": "result", "subtype": "success", "is_error": False, "structured_output": {}})
        return ev

    def test_audit(self):
        agent = {"model": "claude-sonnet-5-5", "tools": ["Read"]}
        ok = self.events([("Read", {"file_path": "/tmp/c/page.png"}), ("StructuredOutput", {})])
        self.assertEqual(agents._audit(ok, agent, "/tmp/c", ["page.png"])[0], [])
        rel_ok = self.events([("Read", {"file_path": "page.png"})])
        self.assertEqual(agents._audit(rel_ok, agent, "/tmp/c", ["page.png"])[0], [])
        rel_bad = self.events([("Read", {"file_path": "../units.md"})])
        self.assertTrue(agents._audit(rel_bad, agent, "/tmp/c", ["page.png"])[0])
        probs = agents._audit(self.events([("Read", {"file_path": "/x/units.md"})]), agent, "/tmp/c", ["page.png"])[0]
        self.assertTrue(any("read /x/units.md" in p for p in probs))
        probs = agents._audit(self.events([("Bash", {})]), agent, "/tmp/c", ["page.png"])[0]
        self.assertTrue(any("used tool Bash" in p for p in probs))
        probs = agents._audit(self.events([], model="claude-opus-5-5"), agent, "/tmp/c", ["page.png"])[0]
        self.assertTrue(any("session model" in p for p in probs))


class Blindness(unittest.TestCase):

    def test_reader_task_has_no_units_values_or_meanings(self):
        """The task is a fixed preamble plus [{id, look_for}] per symbol; nothing else from the table."""
        import json
        rows = uvlib.literature_rows(uvlib.read_table())
        for (src, pg), grp in read_pages.cited_groups(rows).items():
            task = read_pages.task_text(grp)
            head, _, body = task.partition("in this order:\n")
            with self.subTest(page=f"{src}_p{pg}"):
                self.assertEqual(head, read_pages.task_text([]).partition("in this order:\n")[0])
                self.assertNotIn(src, task)
                items = json.loads(body)
                self.assertEqual(items, [{"id": r["symbol"], "look_for": r["printed"]} for r in grp])


class S4Consistency(unittest.TestCase):

    def test_judge(self):
        t = up.parse_unit("s m^-1")
        c = {"read": {"unit_as_printed": "s m^-1"}}
        out = {"verdict": "table_right", "unit": "s/m"}
        self.assertEqual(check_rows.judge_s4(out, c, t), "consistent")
        self.assertTrue(check_rows.judge_s4(dict(out, unit="m/s"), c, t).startswith("inconsistent"))
        self.assertTrue(check_rows.judge_s4(dict(out, verdict="read_right", unit="m s^-1"), c, t)
                        .startswith("inconsistent"))
        self.assertTrue(check_rows.judge_s4(dict(out, unit="parsecs"), c, t).startswith("inconsistent"))
        self.assertEqual(check_rows.judge_s4({"verdict": "cannot_tell", "unit": None}, c, t), "consistent")


class DefinitionPages(unittest.TestCase):

    def test_related_first_then_own_earliest_without_cited(self):
        entry = {"pages": [9, 2, 6], "related_pages": {"paper": [5, 3]}}
        self.assertEqual(check_rows.definition_pages("si", 9, entry), [("paper", 3), ("paper", 5), ("si", 2)])
        self.assertEqual(check_rows.definition_pages("si", 9, {"pages": [9]}), [])

    def test_index_has_related_pages_for_si_rows(self):
        import marker_extract
        rel = marker_extract.related_sources()
        self.assertEqual(rel["moure_fu_2024_SI"], ["moure_fu_2024"])
        self.assertEqual(rel["kaempfer_plapp_2009"], [])


class Calibration(unittest.TestCase):

    def test_planted_table_round_trips(self):
        real = uvlib.literature_rows(uvlib.read_table())
        text, plants = calibrate.planted_table(real, "| h |\n|---|")
        self.assertEqual(len(plants), 10)
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "t.md"
            p.write_text(text)
            rows = {r["symbol"]: r for r in uvlib.read_table(p)}
        self.assertEqual(len(rows), 14)
        self.assertEqual(rows["beta"]["unit"], "s m^1")
        self.assertEqual(rows["k"]["page"], 13)
        self.assertEqual(rows["rho_i"]["printed"], "ρi")
        orig = {r["symbol"]: r for r in real}
        for p_ in plants:   # every plant really changes its field, and changes a unit's dimension or scale
            self.assertNotEqual(p_["original"], p_["planted"])
            if p_["field"] == "unit":
                self.assertFalse(up.same_unit(up.parse_unit(p_["original"]), up.parse_unit(p_["planted"])))
        for s in calibrate.CONTROLS:
            self.assertEqual(rows[s]["unit"], orig[s]["unit"])

    def test_score(self):
        real = uvlib.literature_rows(uvlib.read_table())
        _, plants = calibrate.planted_table(real, "")
        comp = {"rows": [{"symbol": p["symbol"], "verdict": "DISAGREE" if p["field"] != "meaning" else "AGREE",
                          "reasons": []} for p in plants]
                + [{"symbol": s, "verdict": "AGREE", "reasons": []} for s in calibrate.CONTROLS]}
        s4 = {p["symbol"]: {"verdict": "read_right", "unit": p["original"] if p["field"] == "unit" else None,
                            "page": p["original"] if p["field"] == "page" else None}
              for p in plants if p["field"] != "meaning"}
        s4m = {p["symbol"]: {"meaning_verdict": "table_wrong"} for p in plants if p["field"] == "meaning"}
        s4m.update({s: {"meaning_verdict": "matches"} for s in calibrate.CONTROLS})
        r = calibrate.score(plants, comp, s4, s4m)
        self.assertTrue(r["pass"], r["criteria"])
        self.assertEqual(r["control_false_alarms"], 0)
        s4["beta"]["unit"] = "s m^1"   # the checker sides with the planted unit: not a rejection
        s4["a"]["verdict"] = "table_right"
        r = calibrate.score(plants, comp, s4, s4m)
        self.assertEqual(r["criteria"]["S4_unit_page"], "6/8")
        self.assertFalse(r["pass"])


if __name__ == "__main__":
    unittest.main()
