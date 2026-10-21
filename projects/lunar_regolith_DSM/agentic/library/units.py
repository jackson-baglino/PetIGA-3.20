"""The units standard: one registry of symbols, SI units and dimensions that every derivation and every agent brief must use.

Same protection model as the library:
  * the curator code path is the only writer; each change is an event in the log;
  * adding a convention for a NEW symbol is the normal operation;
  * an existing convention is frozen: changing it needs the user's own instruction (a quoted sentence) and human_approved=True;
  * every entry carries its origin: a paper, page and verbatim quote, or an explicit user instruction, or a stated dimensional-analysis basis;
  * readers get a read-only view (`render_markdown`, `lookup`); agents never open the files (deny rules in .claude/settings.json).

  python -m agentic.library.units show                      # the standard as markdown
  python -m agentic.library.units lookup eps
  python -m agentic.library.units add --symbol X --meaning ... --unit "kg m^-3" --kind parameter --paper moure_fu_2024 --page 4 --quote "..."   # curator use
"""
from __future__ import annotations

import argparse
import contextlib
import json
import os
import re
import sqlite3
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONV = Path(os.environ.get("CONVENTIONS_DIR", ROOT / "research" / "Conventions"))

SCHEMA = """
CREATE TABLE IF NOT EXISTS units (symbol TEXT PRIMARY KEY, meaning TEXT, unit TEXT, dim TEXT, kind TEXT, origin TEXT, paper_id TEXT, page INTEGER,
    quote TEXT, basis TEXT, verified TEXT, added_ts TEXT, added_by TEXT);
CREATE TABLE IF NOT EXISTS events (ts TEXT, actor TEXT, action TEXT, symbol TEXT, detail TEXT);
"""
BASE = ("m", "kg", "s", "K", "mol")
UNITS = {"m": (1, 0, 0, 0, 0), "kg": (0, 1, 0, 0, 0), "s": (0, 0, 1, 0, 0), "K": (0, 0, 0, 1, 0), "mol": (0, 0, 0, 0, 1), "1": (0,) * 5,
         "N": (1, 1, -2, 0, 0), "Pa": (-1, 1, -2, 0, 0), "J": (2, 1, -2, 0, 0), "W": (2, 1, -3, 0, 0), "g": (0, 1, 0, 0, 0)}   # g is converted by the user-facing note, never stored as g
KINDS = {"variable", "parameter", "constant", "derived", "scale"}


class Frozen(PermissionError):
    """An existing convention cannot be changed without the user's own instruction."""


def dim_of(unit: str) -> tuple:
    """Dimension vector (m, kg, s, K, mol) of a unit string such as 'kg m^-3', 'J/(kg K)', 'm^2/s', '1'."""
    unit = unit.strip()
    if unit in ("", "1", "-", "dimensionless"):
        return (0,) * 5
    num, _, den = unit.partition("/")
    v = [0] * 5

    def add(part: str, sign: int) -> None:
        part = part.replace("(", " ").replace(")", " ")
        for tok in part.split():
            m = re.fullmatch(r"([A-Za-z]+)(?:\^(-?\d+))?", tok)
            if not m or m.group(1) not in UNITS or m.group(1) == "g":
                raise ValueError(f"unknown unit token {tok!r} in {unit!r} (use SI base units, N, Pa, J, W)")
            e = int(m.group(2) or 1) * sign
            v[:] = [a + e * b for a, b in zip(v, UNITS[m.group(1)])]
    add(num, 1); add(den, -1)
    return tuple(v)


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


class Standard:
    """Read side."""

    def __init__(self, root: Path | None = None):
        self.root = Path(root or CONV)

    @property
    def db(self) -> Path:
        return self.root / "units.sqlite"

    def ro(self) -> sqlite3.Connection:
        con = sqlite3.connect(f"file:{self.db}?mode=ro", uri=True); con.row_factory = sqlite3.Row
        return con

    def lookup(self, symbol: str) -> dict | None:
        con = self.ro()
        try:
            r = con.execute("SELECT * FROM units WHERE symbol=?", (symbol,)).fetchone()
            return dict(r) if r else None
        finally:
            con.close()

    def all(self) -> list[dict]:
        con = self.ro()
        try:
            return [dict(r) for r in con.execute("SELECT * FROM units ORDER BY kind, symbol")]
        finally:
            con.close()

    def render_markdown(self) -> str:
        rows = self.all()
        out = ["# Units standard (trusted; frozen entries; changed only on the user's own instruction)", "",
               "| symbol | meaning | SI unit | kind | origin |", "|---|---|---|---|---|"]
        for r in rows:
            org = f"{r['paper_id']} p.{r['page']}" if r["origin"] == "paper" else r["origin"] + (f" ({r['basis']})" if r["basis"] else "")
            out.append(f"| {r['symbol']} | {r['meaning']} | {r['unit']} | {r['kind']} | {org} |")
        return "\n".join(out) + "\n"

    def brief_text(self) -> str:
        """Compact text for agent briefs."""
        return "; ".join(f"{r['symbol']} [{r['unit']}]" for r in self.all())


class Curator(Standard):
    """Write side: the only code path that changes the standard."""

    def __init__(self, root: Path | None = None, actor: str = "curator"):
        super().__init__(root); self.actor = actor

    @contextlib.contextmanager
    def rw(self):
        self.root.mkdir(parents=True, exist_ok=True)
        if self.db.exists():
            self.db.chmod(0o644)
        con = sqlite3.connect(self.db); con.executescript(SCHEMA)
        try:
            yield con; con.commit()
        finally:
            con.close(); self.db.chmod(0o444)

    def add(self, symbol: str, meaning: str, unit: str, kind: str, *, paper_id: str = "", page: int = 0, quote: str = "", basis: str = "",
            origin: str = "paper", user_instruction: str = "", verified: str = "unverified") -> None:
        """Add a convention for a NEW symbol. Refuses to touch an existing one."""
        dim_of(unit)                                   # rejects a unit string that does not parse
        if kind not in KINDS:
            raise ValueError(f"kind must be one of {sorted(KINDS)}")
        if origin == "paper" and not (paper_id and page and quote.strip()):
            raise ValueError("a paper-origin entry needs paper_id, page and a verbatim quote")
        if origin == "user" and not user_instruction.strip():
            raise ValueError("a user-origin entry needs the user's instruction, quoted")
        if origin == "dimensional-analysis" and not basis.strip():
            raise ValueError("a dimensional-analysis entry needs its basis (the equation and the paper page it comes from)")
        if origin not in ("paper", "user", "dimensional-analysis"):
            raise ValueError("origin must be paper, user or dimensional-analysis")
        with self.rw() as con:
            if con.execute("SELECT 1 FROM units WHERE symbol=?", (symbol,)).fetchone():
                raise Frozen(f"{symbol} already has a convention; it can only be changed with change(), on the user's instruction")
            con.execute("INSERT INTO units VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)", (symbol, meaning, unit, json.dumps(dim_of(unit)), kind, origin, paper_id, page, quote,
                                                                              basis or user_instruction, verified, now(), self.actor))
            con.execute("INSERT INTO events VALUES (?,?,?,?,?)", (now(), self.actor, "add", symbol, f"{unit} ({origin})"))

    def change(self, symbol: str, *, unit: str | None = None, meaning: str | None = None, user_instruction: str, human_approved: bool) -> None:
        """Change an existing convention. Only on the user's own words (quoted) with human_approved=True."""
        if not (user_instruction.strip() and human_approved):
            raise Frozen("changing a convention needs the user's own instruction (quoted) and human_approved=True")
        with self.rw() as con:
            r = con.execute("SELECT unit, meaning FROM units WHERE symbol=?", (symbol,)).fetchone()
            if not r:
                raise KeyError(symbol)
            new_unit = unit or r[0]; new_mean = meaning or r[1]
            con.execute("UPDATE units SET unit=?, dim=?, meaning=? WHERE symbol=?", (new_unit, json.dumps(dim_of(new_unit)), new_mean, symbol))
            con.execute("INSERT INTO events VALUES (?,?,?,?,?)", (now(), self.actor, "change", symbol, f"{r[0]} -> {new_unit}; user: {user_instruction[:300]}"))

    def mark_verified(self, symbol: str, how: str) -> None:
        """Record that the entry was confirmed visually against the page image (the only edit allowed besides change())."""
        with self.rw() as con:
            con.execute("UPDATE units SET verified=? WHERE symbol=?", (how, symbol))
            con.execute("INSERT INTO events VALUES (?,?,?,?,?)", (now(), self.actor, "verify", symbol, how))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("show")
    lk = sub.add_parser("lookup"); lk.add_argument("symbol")
    ad = sub.add_parser("add")
    for k in ("symbol", "meaning", "unit", "kind"):
        ad.add_argument("--" + k, required=True)
    ad.add_argument("--paper", default=""); ad.add_argument("--page", type=int, default=0); ad.add_argument("--quote", default="")
    ad.add_argument("--basis", default=""); ad.add_argument("--origin", default="paper"); ad.add_argument("--user-instruction", default="")
    a = ap.parse_args()
    if a.cmd == "show":
        print(Standard().render_markdown())
    elif a.cmd == "lookup":
        print(json.dumps(Standard().lookup(a.symbol), indent=1))
    else:
        Curator().add(a.symbol, a.meaning, a.unit, a.kind, paper_id=a.paper, page=a.page, quote=a.quote, basis=a.basis, origin=a.origin, user_instruction=a.user_instruction)
        print("added", a.symbol)
    return 0


if __name__ == "__main__":
    sys.exit(main())
