#!/usr/bin/env python3
"""Check the literature-review permission files for the lethal trifecta.

Reads (never writes) the settings files in .claude/permissions/literature-review/
and the project's .claude/settings.json and settings.local.json, which are the
main session's normal profile (WORKFLOW.md §4). For every agent and session it
works out the three legs of the lethal trifecta (Willison 2025) from what the
files allow:

- untrusted text: the inputs.json class "untrusted", or a web tool allowed;
- a way out: a web tool allowed, or Bash with network to any host (the model
  provider itself does not count: every call's input goes there by design);
- private data: the inputs.json class "user_private".

It fails if an entry holds all three and is not in accepted_risks.json, and
on structural faults: an agent that may edit or write, a sandbox that can fail
open, a review mode or normal profile missing a store deny or the read hook,
a launcher script without its profile.

    venv_workflows/bin/python scripts/literature_review/check_permissions.py           # check
    venv_workflows/bin/python scripts/literature_review/check_permissions.py --write   # also refresh the tables

The tables go between the GENERATED markers in PERMISSIONS.md and in
workflows/literature-review/WORKFLOW.md §4. Without --write, stale tables fail.
"""
import argparse
import json
import os
import sys

WEB_TOOLS = {"WebFetch", "WebSearch"}
# Untrusted stores (relative to the project) that review mode and the normal
# profile must deny to the Read tool (WORKFLOW.md §4, §6).
STORE_DENIES = [
    "Read(/library/papers/**)", "Read(/library/memos/**)", "Read(/library/history/**)",
    "Read(/library/views/**)", "Read(/library/.staging/**)", "Read(/library/events.jsonl)",
    "Read(/library/catalog.sqlite)", "Read(/library/references.bib)",
    "Read(/library_practice/**)", "Read(/index/**)",
]
# Path-scoped Write(...) rules are ignored by Claude Code; Edit(...) covers every file-editing tool.
NORMAL_EXTRA_DENIES = ["Read(/runs/literature-review/**/untrusted/**)", "Edit(/library/**)"]
REVIEW_DENIES = ["WebFetch", "WebSearch", "Bash(curl *)", "Bash(wget *)", "Bash(git push *)",
                 "Bash(security *)"]
HOOK = "review_read_allowlist.py"
BEGIN = "<!-- BEGIN GENERATED: check_permissions.py (literature-review) -->"
END = "<!-- END GENERATED -->"


def load(path):
    with open(path) as f:
        return json.load(f)


def tool_names(rules):
    return {r.split("(")[0] for r in rules}


def sandbox_faults(sb, where):
    faults = []
    if not sb.get("enabled"):
        faults.append(f"{where}: sandbox not enabled")
    if sb.get("failIfUnavailable") is not True:
        faults.append(f"{where}: sandbox.failIfUnavailable is not true (would fail open)")
    if sb.get("allowUnsandboxedCommands") is not False:
        faults.append(f"{where}: sandbox.allowUnsandboxedCommands is not false")
    return faults


def agent_entry(name, s, classes):
    allow = s.get("permissions", {}).get("allow", [])
    deny = s.get("permissions", {}).get("deny", [])
    tools = tool_names(allow)
    sb = s.get("sandbox", {})
    faults = sandbox_faults(sb, name)
    for t in ("Edit", "Write", "NotebookEdit"):
        if t in tools or t not in deny:
            faults.append(f"{name}: {t} must be denied")
    if "WebFetch" in tools:
        faults.append(f"{name}: no agent may fetch pages")
    net = sb.get("network", {}).get("allowedDomains", [])
    untrusted = "untrusted" in classes or bool(tools & WEB_TOOLS)
    way_out = bool(tools & WEB_TOOLS) or ("Bash" in tools and bool(net))
    private = "user_private" in classes
    return {"name": name, "kind": "agent", "tools": sorted(allow) or ["none"],
            "network": "web search" if "WebSearch" in tools else (", ".join(net) or "none"),
            "inputs": classes, "legs": (untrusted, way_out, private), "faults": faults}


def review_entry(s, classes):
    deny = s.get("permissions", {}).get("deny", [])
    faults = sandbox_faults(s.get("sandbox", {}), "review-mode")
    for r in REVIEW_DENIES + STORE_DENIES:
        if r not in deny:
            faults.append(f"review-mode: missing deny {r}")
    hooks = json.dumps(s.get("hooks", {}).get("PreToolUse", []))
    if HOOK not in hooks or "Read|Grep|Glob" not in hooks:
        faults.append("review-mode: the Read/Grep/Glob allow-list hook is not registered")
    fs = s.get("sandbox", {}).get("filesystem", {})
    if "~/" not in fs.get("denyRead", []):
        faults.append("review-mode: sandbox.filesystem.denyRead must include ~/")
    net = s.get("sandbox", {}).get("network", {}).get("allowedDomains", [])
    # Untrusted text is kept out by the denies, the hook and the sandbox (checked above).
    untrusted = "untrusted" in classes
    way_out = bool(net)
    private = "user_private" in classes
    return {"name": "review-mode", "kind": "session", "tools": ["all but web tools"],
            "network": ", ".join(net) or "none", "inputs": classes,
            "legs": (untrusted, way_out, private), "faults": faults}


def normal_entry(project, local, classes):
    deny = project.get("permissions", {}).get("deny", []) + local.get("permissions", {}).get("deny", [])
    allow = project.get("permissions", {}).get("allow", []) + local.get("permissions", {}).get("allow", [])
    faults = [f"normal: missing deny {r}" for r in STORE_DENIES + NORMAL_EXTRA_DENIES if r not in deny]
    web_allowed = sorted(r for r in allow if tool_names([r]) & WEB_TOOLS)
    web_denied = WEB_TOOLS <= tool_names([r for r in deny if "(" not in r])
    untrusted = "untrusted" in classes
    way_out = not web_denied
    private = "user_private" in classes
    return {"name": "normal", "kind": "session", "tools": ["all"],
            "network": ("; ".join(web_allowed) + "; others ask") if web_allowed else "web tools ask",
            "inputs": classes, "legs": (untrusted, way_out, private), "faults": faults}


def launcher_faults(perm_dir):
    faults = []
    path = os.path.join(perm_dir, "launcher.allow")
    with open(path) as f:
        lines = f.readlines()
    for ln in lines:
        ln = ln.split("#")[0].strip()
        if not ln:
            continue
        script, profile = ln.split()
        sb = os.path.join(perm_dir, "profiles", profile + ".sb")
        if not os.path.exists(sb):
            faults.append(f"launcher.allow: {script} names missing profile {profile}")
            continue
        with open(sb) as f:
            text = f.read()
        if "(deny network*)" not in text or "(deny file-write*)" not in text:
            faults.append(f"profiles/{profile}.sb: must start from no network and no writes")
    return faults


def collect(root):
    perm = os.path.join(root, ".claude", "permissions", "literature-review")
    inputs = load(os.path.join(perm, "inputs.json"))
    accepted = load(os.path.join(perm, "accepted_risks.json"))
    entries = []
    adir = os.path.join(perm, "agents")
    for fn in sorted(os.listdir(adir)):
        if fn.endswith(".settings.json"):
            name = fn[: -len(".settings.json")]
            if name not in inputs["agents"]:
                entries.append({"name": name, "kind": "agent", "tools": [], "network": "", "inputs": [],
                                "legs": (False, False, False),
                                "faults": [f"{name}: no input classes in inputs.json"]})
                continue
            entries.append(agent_entry(name, load(os.path.join(adir, fn)), inputs["agents"][name]))
    for name in inputs["agents"]:
        if not os.path.exists(os.path.join(adir, name + ".settings.json")):
            entries.append({"name": name, "kind": "agent", "tools": [], "network": "", "inputs": [],
                            "legs": (False, False, False),
                            "faults": [f"{name}: listed in inputs.json but has no settings file"]})
    entries.append(review_entry(load(os.path.join(perm, "sessions", "review-mode.json")),
                                inputs["sessions"]["review-mode"]))
    cdir = os.path.join(root, ".claude")
    local = os.path.join(cdir, "settings.local.json")
    entries.append(normal_entry(load(os.path.join(cdir, "settings.json")),
                                load(local) if os.path.exists(local) else {},
                                inputs["sessions"]["normal"]))
    faults = launcher_faults(perm)
    for e in entries:
        faults += e["faults"]
        if all(e["legs"]):
            if e["name"] in accepted:
                e["accepted"] = True
            else:
                faults.append(f"{e['name']}: holds all three legs of the trifecta "
                              "(untrusted text, a way out, private data)")
    return entries, faults


def table(entries):
    yes = lambda b: "yes" if b else "no"
    rows = ["| Agent or session | Allowed tools | Network | Inputs | Untrusted text | Way out | Private data |",
            "|---|---|---|---|---|---|---|"]
    for e in entries:
        u, w, p = e["legs"]
        name = f"`{e['name']}`" + (" (accepted risk)" if e.get("accepted") else "")
        tools = ", ".join(f"`{t}`" if t not in ("none", "all", "all but web tools") else t
                          for t in e["tools"])
        rows.append(f"| {name} | {tools} | {e['network']} | {', '.join(e['inputs'])} | "
                    f"{yes(u)} | {yes(w)} | {yes(p)} |")
    return ("Generated from `.claude/permissions/literature-review/` and the project "
            "settings by `scripts/literature_review/check_permissions.py`; do not edit by hand.\n\n"
            + "\n".join(rows) + "\n")


def refresh(path, body, write):
    with open(path) as f:
        text = f.read()
    if BEGIN not in text or END not in text:
        return f"{path}: GENERATED markers missing"
    head, rest = text.split(BEGIN, 1)
    _, tail = rest.split(END, 1)
    new = head + BEGIN + "\n" + body + END + tail
    if new == text:
        return None
    if write:
        with open(path, "w") as f:
            f.write(new)
        return None
    return f"{path}: permission table is stale (run with --write)"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    ap.add_argument("--write", action="store_true", help="refresh the generated tables")
    args = ap.parse_args(argv)
    entries, faults = collect(args.root)
    body = table(entries)
    for doc in ("PERMISSIONS.md", "workflows/literature-review/WORKFLOW.md"):
        f = refresh(os.path.join(args.root, doc), body, args.write)
        if f:
            faults.append(f)
    print(body)
    if faults:
        print("FAIL")
        for f in faults:
            print(" -", f)
        return 1
    print("OK: no agent or session holds all three legs, except accepted risks:",
          ", ".join(e["name"] for e in entries if e.get("accepted")) or "none")
    return 0


if __name__ == "__main__":
    sys.exit(main())
