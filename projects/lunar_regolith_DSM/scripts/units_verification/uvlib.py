"""Shared helpers for the units-verification workflow (workflows/units-verification/WORKFLOW.md).

- read the units table (memory/semantic/units.md) into rows
- hash files, record git and tool versions
- merge a stage's record into runs/units-verification/<run_id>/manifest.json
- refuse to run a real stage when the workflow's definitions are not committed
"""
import datetime
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]                      # projects/lunar_regolith_DSM
UNITS_MD = ROOT / "memory" / "semantic" / "units.md"
WORKFLOW_DIR = ROOT / "workflows" / "units-verification"
# UV_RUNS_DIR redirects runs for tests and development; real runs use the default.
RUNS_DIR = Path(os.environ.get("UV_RUNS_DIR", ROOT / "runs" / "units-verification"))
LITERATURE = Path.home() / ("Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/Projects/"
                            "lunar_regolith_DSM/Literature")
# Paths whose committed state defines the workflow; a real run refuses if any is modified.
DEFINITION_PATHS = ["workflows/units-verification", "scripts/units_verification", ".claude/agents"]
COLUMNS = ("symbol", "meaning", "unit", "source", "page", "quote", "verified")
LITERATURE_SOURCES = ("kaempfer_plapp_2009", "moure_fu_2024", "moure_fu_2024_SI")


class StageError(RuntimeError):
    """A stage stops: the message says why and what the user should do."""


def read_table(path=UNITS_MD):
    """Rows of the units table as dicts (COLUMNS); `printed` is the form after '(printed as ...)', if any."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith("| `"):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split(" | ")]
        if len(cells) != len(COLUMNS):
            raise StageError(f"units table: expected {len(COLUMNS)} cells, got {len(cells)}: {line[:80]}")
        row = dict(zip(COLUMNS, cells))
        row["symbol"] = row["symbol"].strip("`")
        row["page"] = int(row["page"]) if row["page"] else None
        m = re.search(r"\(printed as ([^)]+)\)", row["meaning"])
        row["printed"] = m.group(1).strip() if m else None
        rows.append(row)
    symbols = [r["symbol"] for r in rows]
    dup = sorted({s for s in symbols if symbols.count(s) > 1})
    if dup:
        raise StageError(f"units table: duplicate symbols {dup}")
    return rows


def literature_rows(rows):
    """Rows in scope for this pass: literature symbols only (user decision, 2026-10-22)."""
    return [r for r in rows if r["source"] in LITERATURE_SOURCES]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def run(cmd, **kw):
    """Run a command; return stdout. Raises StageError with stderr on failure."""
    p = subprocess.run(cmd, capture_output=True, text=kw.pop("text", True), **kw)
    if p.returncode != 0:
        raise StageError(f"command failed ({p.returncode}): {' '.join(map(str, cmd))}\n{p.stderr}")
    return p.stdout


def git_state():
    commit = run(["git", "rev-parse", "HEAD"], cwd=ROOT).strip()
    dirty = run(["git", "status", "--porcelain", "--", *DEFINITION_PATHS], cwd=ROOT).strip()
    return {"commit": commit, "definitions_modified": dirty.splitlines() if dirty else []}


def require_committed(allow_dirty):
    """A workflow run uses only committed definitions (CLAUDE.md, Reproducibility)."""
    g = git_state()
    if g["definitions_modified"] and not allow_dirty:
        raise StageError("workflow definitions have uncommitted changes; commit first:\n  "
                         + "\n  ".join(g["definitions_modified"]))
    return g


def tool_versions():
    v = {}
    for name, cmd in (("claude_code", ["claude", "--version"]), ("mutool", ["mutool", "-v"])):
        try:
            p = subprocess.run(cmd, capture_output=True, text=True)
            v[name] = (p.stdout + p.stderr).strip().splitlines()[0]
        except (OSError, IndexError):
            v[name] = None
    return v


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def new_run_id():
    return datetime.datetime.now().strftime("%Y%m%d-%H%M%S")


def run_dir(run_id):
    d = RUNS_DIR / run_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def load_manifest(rdir):
    p = Path(rdir) / "manifest.json"
    return json.loads(p.read_text()) if p.exists() else {}


def update_manifest(rdir, stage, record, allow_dirty=False):
    """Merge one stage's record into the run manifest; each stage carries its own git commit and time."""
    rdir = Path(rdir)
    m = load_manifest(rdir)
    m.setdefault("workflow", "units-verification")
    m.setdefault("run_id", rdir.name)
    m.setdefault("created", now())
    m.setdefault("stages", {})
    record = dict(record)
    record["git"] = git_state()
    record["allow_dirty"] = bool(allow_dirty)
    record["tools"] = tool_versions()
    record["finished"] = now()
    m["stages"][stage] = record
    tmp = rdir / "manifest.json.tmp"
    tmp.write_text(json.dumps(m, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, rdir / "manifest.json")
    return m


def abspath(shown):
    """Inverse of rel(): a path recorded with ~ or relative to the project root, as an absolute Path."""
    p = Path(shown).expanduser()
    return p if p.is_absolute() else ROOT / p


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, path)


def rel(path):
    """Path relative to the project root when inside it; else with the home folder shown as ~."""
    path = Path(path).resolve()
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        home = str(Path.home())
        s = str(path)
        return "~" + s[len(home):] if s.startswith(home) else s
