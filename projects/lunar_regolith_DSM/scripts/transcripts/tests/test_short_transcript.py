"""Tests for short_transcript. Run from the project root:
    venv_workflows/bin/python -m unittest discover -s scripts/transcripts/tests -v
"""
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import short_transcript as st  # noqa: E402


def user(content, **kw):
    return json.dumps({"type": "user", "message": {"content": content}, **kw})


def asst(*blocks, model="claude-opus-5-5"):
    return json.dumps({"type": "assistant", "message": {"model": model, "content": list(blocks)}})


def text(t):
    return {"type": "text", "text": t}


TOOL = {"type": "tool_use", "id": "t1", "name": "Bash", "input": {}}
RESULT = [{"type": "tool_result", "tool_use_id": "t1", "content": "ok"}]


class ParseTurns(unittest.TestCase):

    def test_keeps_only_text_after_last_tool_call(self):
        turns, models = st.parse_turns([
            user("Read the paper."),
            asst(text("I'll look first."), TOOL), user(RESULT),
            asst({"type": "thinking", "thinking": ""}, text("Here is the answer.")),
        ])
        self.assertEqual(len(turns), 1)
        self.assertEqual(turns[0]["reply"], ["Here is the answer."])
        self.assertEqual(models, ["claude-opus-5-5"])

    def test_cut_off_prompt_is_replaced_by_full_copy(self):
        turns, _ = st.parse_turns([user("Weekly homework this"), user("Weekly homework this week."),
                                   asst(text("Done."))])
        self.assertEqual([t["text"] for t in turns], ["Weekly homework this week."])
        self.assertEqual(turns[0]["reply"], ["Done."])

    def test_skips_meta_commands_and_marks_notices(self):
        turns, _ = st.parse_turns([
            user("A question?"), asst(text("Working on it.")),
            user("<system-reminder>x</system-reminder>", isMeta=True),
            user("<task-notification>done</task-notification>"), asst(text("Finished.")),
            user("<command-name>/export</command-name>"), user("<local-command-stdout>ok</local-command-stdout>"),
            user("Next?"),
        ])
        self.assertEqual([t["kind"] for t in turns], ["prompt", "notice", "prompt"])
        self.assertEqual(turns[1]["reply"], ["Finished."])
        self.assertEqual(turns[2]["reply"], [])

    def test_markdown_numbers_prompts_and_applies_labels(self):
        turns, _ = st.parse_turns([user("one"), asst(text("r1")), user("two\nlines")])
        md, n = st.to_markdown(turns, "T", "intro", st.parse_labels("1=Setup, 2=Round 1"))
        self.assertEqual(n, 2)
        self.assertIn("## Prompt 1 (Setup)", md)
        self.assertIn("> two\n> lines", md)
        self.assertIn("[No reply in the log", md)


if __name__ == "__main__":
    unittest.main()
