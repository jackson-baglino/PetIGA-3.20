#!/usr/bin/env python3
"""T-restore: restore one paper and the memo's citations from the backup, check hashes, rerun Y2 (WORKFLOW.md §12).

    restore_test.py --run <run folder>

From the backup named in <run>/summary_backup.json it copies into an empty folder: the latest memo and its claim
list, every library claim file the memo cites, and every library file of the first paper the memo cites. Each
copied file's sha256 must match the backup manifest, and Y2 (check_citations.check) on the restored memo and
restored claims must give the same result as on the run's own files. Adds restore_ok and the counts to
summary_backup.json; exits 1 if the test fails.
"""
import argparse
import json
import re
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import backup as bk  # noqa: E402
import check_citations as cx  # noqa: E402


def restore_test(run, question="Q1"):
    run = Path(run).resolve()
    s = json.loads((run / "summary_backup.json").read_text())
    dest = Path(s["dest"])
    man = json.loads((dest / "manifest.json").read_text())["files"]
    if bk.sha256(dest / "manifest.json") != s["manifest_sha256"]:
        raise SystemExit("backup manifest changed since the backup")
    k = cx.latest_version(dest / "run" / "untrusted" / "memos", question)
    memo_rel = f"run/untrusted/memos/{question}.v{k}.md"
    claims_rel = f"run/untrusted/memos/{question}.v{k}.claims.json"
    memo = (dest / memo_rel).read_text()
    cited = sorted(set(re.findall(r"\[([a-z0-9_]+-c\d{3})", memo)))
    want = {memo_rel, claims_rel}
    for cid in cited:
        paper = cid.rsplit("-c", 1)[0]
        want |= {r for r in man if r.startswith(f"library/papers/{paper}/claims/{cid}.v")}
    if cited:
        first = cited[0].rsplit("-c", 1)[0]
        want |= {r for r in man if r.startswith(f"library/papers/{first}/")}
    mismatches = []
    with tempfile.TemporaryDirectory() as tmp:
        for r in sorted(want):
            if r not in man:
                mismatches.append(f"{r}: not in the manifest")
                continue
            t = Path(tmp) / r
            t.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(dest / r, t)
            if bk.sha256(t) != man[r]["sha256"]:
                mismatches.append(f"{r}: hash differs")
        restored = cx.check((Path(tmp) / memo_rel).read_text(), json.loads((Path(tmp) / claims_rel).read_text()))
    own = cx.check((run / "untrusted" / "memos" / f"{question}.v{k}.md").read_text(),
                   json.loads((run / "untrusted" / "memos" / f"{question}.v{k}.claims.json").read_text()))
    ok = not mismatches and restored == own
    s.update({"restore_ok": ok, "restored_files": len(want), "restore_mismatches": mismatches[:20],
              "restored_y2_pass": restored["pass"]})
    (run / "summary_backup.json").write_text(json.dumps(s, indent=1))
    return s


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--question", default="Q1")
    a = ap.parse_args(argv)
    s = restore_test(a.run, a.question)
    print(json.dumps({k: s[k] for k in ("restore_ok", "restored_files", "restored_y2_pass")}))
    return 0 if s["restore_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
