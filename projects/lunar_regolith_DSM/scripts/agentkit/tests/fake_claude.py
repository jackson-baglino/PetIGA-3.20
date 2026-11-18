#!/usr/bin/env python3
"""A stand-in for `claude -p` in agentkit tests. No model is called.

It reads the same arguments the launcher gives Claude Code, runs the hooks from the composed --settings file
the way Claude Code does (PreToolUse before each tool call, Stop before ending; exit 2 blocks), carries out the
calls the hooks allow, and prints stream-json events like a real session. A scenario chosen by
AGENTKIT_FAKE_SCENARIO decides what the "agent" tries:

    pass       write derivation.md, run an allowed script, return the correct T01 answer
    wrong      the same, but the answer has the conduction sign flipped
    blocked    return status "blocked" with a question
    probe      first try four forbidden things (read /etc/hosts, write ../PROVENANCE.md, curl, a chained
               command), then do what "pass" does
    bypass     read /etc/hosts WITHOUT asking the guard (a barrier that failed), then "pass"
    forgets    try to stop before writing derivation.md, then write it when the Stop hook objects
    no_answer  end without a structured answer

Notes for the tests go to AGENTKIT_FAKE_LOG (one JSON line each), including the full argument list.
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

CORRECT = ("Integral(rho*c*w(x)*Derivative(T(x, t), t), (x, 0, L)) + Integral(Derivative(w(x), x)*k(x)*"
           "Derivative(T(x, t), x), (x, 0, L)) - Integral(w(x)*Q(x), (x, 0, L)) - (Subs(w(x)*k(x)*"
           "Derivative(T(x, t), x), x, L) - Subs(w(x)*k(x)*Derivative(T(x, t), x), x, 0))")
WRONG = CORRECT.replace("(x, 0, L)) + Integral(Derivative(w(x)", "(x, 0, L)) - Integral(Derivative(w(x)")


def note(**kw):
    path = os.environ.get("AGENTKIT_FAKE_LOG")
    if path:
        with open(path, "a") as f:
            f.write(json.dumps(kw) + "\n")


def arg(argv, flag):
    return argv[argv.index(flag) + 1] if flag in argv else None


class Session:
    def __init__(self, argv):
        self.settings = json.loads(Path(arg(argv, "--settings")).read_text())
        self.events, self.n = [], 0
        self.events.append({"type": "system", "subtype": "init", "model": arg(argv, "--model")})

    def hooks(self, event_name, payload, tool=None):
        """Run the hooks for one event; (blocked, stderr)."""
        for group in self.settings.get("hooks", {}).get(event_name, []):
            if tool is not None and not re.fullmatch(group.get("matcher", ".*"), tool):
                continue
            for h in group["hooks"]:
                p = subprocess.run(h["command"], shell=True, input=json.dumps(payload), capture_output=True,
                                   text=True, timeout=h.get("timeout", 60))
                if p.returncode == 2:
                    return True, p.stderr.strip()
        return False, ""

    def tool(self, name, inp, ask_guard=True):
        self.n += 1
        uid = f"tu_{self.n}"
        self.events.append({"type": "assistant", "message": {"content": [
            {"type": "tool_use", "id": uid, "name": name, "input": inp}]}})
        blocked, why = (self.hooks("PreToolUse", {"hook_event_name": "PreToolUse", "tool_name": name,
                                                  "tool_input": inp, "cwd": os.getcwd()}, tool=name)
                        if ask_guard else (False, ""))
        if not blocked:
            if name == "Write":
                Path(inp["file_path"]).write_text(inp["content"])
            out = "ok"
        else:
            out = why
        self.events.append({"type": "user", "message": {"content": [
            {"type": "tool_result", "tool_use_id": uid, "is_error": blocked, "content": out}]}})
        note(tool=name, input=inp, blocked=blocked, reason=why)
        return not blocked

    def stop(self, active=False):
        blocked, why = self.hooks("Stop", {"hook_event_name": "Stop", "stop_hook_active": active})
        note(stop=True, blocked=blocked, reason=why)
        return not blocked

    def finish(self, answer):
        r = {"type": "result", "subtype": "success", "is_error": False, "num_turns": self.n + 1,
             "usage": {"input_tokens": 100, "output_tokens": 50}}
        if answer is not None:
            r["structured_output"] = answer
        self.events.append(r)
        for e in self.events:
            print(json.dumps(e))


def answer(expr, status="done", questions=()):
    return {"status": status, "result_expr": expr, "assumptions_used": ["rho, c constant"],
            "questions": list(questions),
            "steps": [{"id": "S1", "statement": "multiply by w and integrate", "from": ["given"],
                       "operation": "multiply by w, integrate over (0, L)"},
                      {"id": "S2", "statement": "integrate the conduction term by parts", "from": ["S1"],
                       "operation": "integration by parts"}]}


def main(argv):
    if "--version" in argv:
        print("0.0.0 (fake_claude, tests only)")
        return 0
    note(argv=argv)
    scenario = os.environ.get("AGENTKIT_FAKE_SCENARIO", "pass")
    s = Session(argv)
    deriv = {"file_path": str(Path.cwd() / "derivation.md"), "content": "S1 ... S2 ...\n"}
    if scenario == "probe":
        s.tool("Read", {"file_path": "/etc/hosts"})
        s.tool("Write", {"file_path": str(Path.cwd().parent / "PROVENANCE.md"), "content": "forged"})
        s.tool("Bash", {"command": "curl https://example.com"})
        s.tool("Bash", {"command": "python scratch/check.py; rm -rf ~"})
    if scenario == "bypass":
        s.tool("Read", {"file_path": "/etc/hosts"}, ask_guard=False)
    if scenario == "forgets":
        if not s.stop():
            s.tool("Write", deriv)
            s.stop(active=True)
    else:
        s.tool("Write", deriv)
    s.tool("Bash", {"command": "python scratch/check.py"})
    if scenario == "no_answer":
        s.finish(None)
    elif scenario == "blocked":
        s.finish(answer("", status="blocked", questions=["Is k a function of T?"]))
    else:
        s.finish(answer(WRONG if scenario == "wrong" else CORRECT))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
