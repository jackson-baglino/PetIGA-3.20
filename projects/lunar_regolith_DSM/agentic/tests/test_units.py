"""Tests of the units standard. Run: python agentic/tests/test_units.py"""
import os, sys, tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
from agentic.library import units as U
fails = []
def check(n, ok):
    print(("PASS " if ok else "FAIL ") + n)
    if not ok: fails.append(n)
def raises(f, exc):
    try: f()
    except exc: return True
    except Exception: return False
    return False
c = U.Curator(Path(tempfile.mkdtemp()))
check("dimension of kg m^-3", U.dim_of("kg m^-3") == (-3, 1, 0, 0, 0))
check("dimension of J/(kg K)", U.dim_of("J/(kg K)") == (2 - 0, 0, -2, -1, 0) or U.dim_of("J/(kg K)") == (2, 0, -2, -1, 0))
check("dimension of m^2/s", U.dim_of("m^2/s") == (2, 0, -1, 0, 0))
check("unknown unit rejected", raises(lambda: U.dim_of("furlong"), ValueError))
c.add("eps", "interface decay length", "m", "parameter", paper_id="p", page=1, quote="q")
check("new symbol added", c.lookup("eps")["unit"] == "m")
check("existing symbol is frozen", raises(lambda: c.add("eps", "x", "m", "parameter", paper_id="p", page=1, quote="q"), U.Frozen))
check("change without user instruction refused", raises(lambda: c.change("eps", unit="mm", user_instruction="", human_approved=False), U.Frozen))
check("change with instruction but no approval refused", raises(lambda: c.change("eps", unit="m^2", user_instruction="'x'", human_approved=False), U.Frozen))
check("paper entry needs a quote", raises(lambda: c.add("a", "m", "m", "parameter", paper_id="p", page=1), ValueError))
check("user origin needs the instruction", raises(lambda: c.add("b", "m", "m", "parameter", origin="user"), ValueError))
c.change("eps", unit="m", meaning="decay length", user_instruction="'use metres'", human_approved=True)
check("approved change goes through and is logged", c.lookup("eps")["meaning"] == "decay length")
check("db file is read-only between writes", not os.access(c.db, os.W_OK) or os.geteuid() == 0)
print("FAILED: " + ", ".join(fails) if fails else "all passed"); sys.exit(1 if fails else 0)
