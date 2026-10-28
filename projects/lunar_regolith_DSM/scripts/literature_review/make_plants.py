#!/usr/bin/env python3
"""Make mechanical planted errors from the user's gold statements (C0 core, WORKFLOW.md §5).

Each plant changes a verified statement in exactly one mechanical way:

- changed_digit:              one digit of a number in the text or the value field;
- equation_sign_or_exponent:  a sign or an exponent flipped (-55 -> 55, 10^{90} -> 10^{-90},
                              1.8e6 -> 1.8e-6);
- wrong_unit:                 a unit swapped for another of the same kind (microns -> millimeters,
                              minutes -> hours, J m^-3 -> J m^-2);
- wrong_page:                 the cited page moved by one or two, inside the paper.

Numbers that identify things (Fig. 1, equation (31), Table 2, SA81) are never changed.
Quote plants (fabricated_quote) need a claim's verbatim quote, so they are made later
from verified claims, not from gold statements written in the user's words.

The output is for the user to review: every plant is listed with the statement it
came from. The choice of which plants to keep is random but fixed by --seed.

    venv_workflows/bin/python scripts/literature_review/make_plants.py [--seed 20261028] [--per-statement 2]
"""
import argparse
import hashlib
import os
import random
import re
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_calibration as cc  # noqa: E402

OUT = os.path.join(cc.CAL, "generated", "script_plants.yaml")

# A number inside an identifier or a reference is not a value.
REF_BEFORE = re.compile(r"(Fig\.?|Figure|Figs\.?|[Ee]q\.?|[Ee]quations?|Table|Sec\.?|Section|"
                        r"Ref\.?|p\.|pp\.)\s*\(?$")
NUMBER = re.compile(r"(?<![A-Za-z_\\{\d.,])(-|−)?(\d{1,3}(?:,\d{3})+|\d+(?:\.\d+)?)(?![A-Za-z_\d]|\^\{)")
EXPONENT = re.compile(r"\^\{(-?)(\d+)\}")
UNIT_SWAPS = [
    ("microns", "millimeters"), ("micron", "millimeter"), ("micrometers", "millimeters"),
    ("nanometers", "micrometers"), ("nanometer", "micrometer"),
    ("minutes", "hours"), ("hours", "days"), ("years", "days"), ("seconds", "minutes"),
    (r"^\circ$C", "$ K"), (r"^\circ$ C", "$ K"),
]
UNIT_FIELD_SWAPS = [(r"m\^-3", "m^-2"), (r"m\^-2", "m^-3"), (r"K\^-1", "K"), (r"s\^-1", "s"),
                    (r"\bkg\b", "g"), (r"\bmol\b", "kmol")]


def value_numbers(text):
    """(start, end, sign, digits) of numbers that are values, not references."""
    out = []
    for m in NUMBER.finditer(text):
        before = text[max(0, m.start() - 12):m.start()]
        if REF_BEFORE.search(before):
            continue
        if re.search(r"[A-Za-z]$", before) and not re.search(r"\s$", before):
            continue
        if "^{" in text[max(0, m.start() - 2):m.start()]:
            continue                          # exponents are handled by EXPONENT
        out.append((m.start(), m.end(), m.group(1) or "", m.group(2)))
    return out


def changed_digit(text, rng):
    nums = value_numbers(text)
    if not nums:
        return None
    s, e, sign, digits = rng.choice(nums)
    pos = [i for i, c in enumerate(digits) if c.isdigit()]
    i = rng.choice(pos)
    old = digits[i]
    choices = [d for d in "0123456789" if d != old and not (i == pos[0] and d == "0" and len(pos) > 1)]
    new = digits[:i] + rng.choice(choices) + digits[i + 1:]
    start = s + len(sign)
    return text[:start] + new + text[start + len(digits):], f"digit changed: {digits} -> {new}"


def changed_value_digit(value, rng):
    """A digit of the mantissa of a value field (1.8e6 -> 1.3e6), never of its exponent."""
    m = re.fullmatch(r"(-?)(\d+(?:\.\d+)?)([eE][-+]?\d+)?", value.strip())
    if not m:
        return None
    mant = m.group(2)
    pos = [i for i, c in enumerate(mant) if c.isdigit()]
    i = rng.choice(pos)
    choices = [d for d in "0123456789" if d != mant[i] and not (i == pos[0] and d == "0" and len(pos) > 1)]
    new = mant[:i] + rng.choice(choices) + mant[i + 1:]
    out = m.group(1) + new + (m.group(3) or "")
    return out, f"digit changed: {value} -> {out}"


def sign_or_exponent(text, value, rng):
    opts = []
    for m in EXPONENT.finditer(text):
        flipped = "" if m.group(1) else "-"
        opts.append(("text", m.start(), m.end(), "^{" + flipped + m.group(2) + "}",
                     f"exponent flipped: {m.group(0)} -> ^{{{flipped}{m.group(2)}}}"))
    for s, e, sign, digits in value_numbers(text):
        if sign:
            opts.append(("text", s, e, digits, f"sign dropped: {sign}{digits} -> {digits}"))
    if value is not None and re.fullmatch(r"-?\d+(\.\d+)?[eE][-+]?\d+", str(value)):
        v = str(value)
        mant, exp = re.split(r"[eE]", v)
        new = mant + "e" + (exp[1:] if exp.startswith("-") else "-" + exp.lstrip("+"))
        opts.append(("value", 0, 0, new, f"exponent flipped: {v} -> {new}"))
    if not opts:
        return None
    where, s, e, new, note = rng.choice(opts)
    if where == "value":
        return text, new, note
    return text[:s] + new + text[e:], value, note


def wrong_unit(text, unit, rng):
    opts = []
    for a, b in UNIT_SWAPS:
        if a in text:
            opts.append(("text", a, b))
    for a, b in UNIT_FIELD_SWAPS:
        if unit and re.search(a, str(unit)):
            opts.append(("unit", a, b))
    if not opts:
        return None
    where, a, b = rng.choice(opts)
    if where == "unit":
        new = re.sub(a, b, str(unit), count=1)
        return text, new, f"unit changed: {unit} -> {new}"
    return text.replace(a, b, 1), unit, f"unit changed: {a.strip()} -> {b.strip()}"


def wrong_page(page, pmap, rng):
    main = page[0] if isinstance(page, list) else page
    if not isinstance(main, int):
        return None
    shifts = [d for d in (-2, -1, 1, 2) if pmap is None or pmap.to_pdf(main + d)]
    if not shifts:
        return None
    new = main + rng.choice(shifts)
    return new, f"page moved: {main} -> {new}"


def plants_for(i, e, pp, rng):
    text = str(e.get("statement") or "")
    out = []
    base = {"paper": e["paper"], "from_statement": i, "location": e.get("location"),
            "kind_of_claim": e.get("kind"), "correct_claim": text, "correct_page": e.get("page"),
            "correct_value": e.get("value"), "correct_unit": e.get("unit")}
    r = changed_digit(text, rng)
    if r:
        out.append(dict(base, kind="changed_digit", planted_claim=r[0], planted_page=e.get("page"),
                        note=r[1]))
    if e.get("value") is not None:
        r = changed_value_digit(str(e["value"]), rng)
        if r:
            out.append(dict(base, kind="changed_digit", planted_claim=text, planted_value=r[0],
                            planted_page=e.get("page"), note="value " + r[1]))
    r = sign_or_exponent(text, e.get("value"), rng)
    if r:
        p = dict(base, kind="equation_sign_or_exponent", planted_claim=r[0],
                 planted_page=e.get("page"), note=r[2])
        if r[1] != e.get("value"):
            p["planted_value"] = r[1]
        out.append(p)
    r = wrong_unit(text, e.get("unit"), rng)
    if r:
        p = dict(base, kind="wrong_unit", planted_claim=r[0], planted_page=e.get("page"), note=r[2])
        if r[1] != e.get("unit"):
            p["planted_unit"] = r[1]
        out.append(p)
    info = pp.get(e["paper"])
    pmap = None
    if info and info["file"]:
        info["map"] = info["map"] or cc.PageMap(info["file"])
        pmap = info["map"] if info["map"].offset is not None else None
    r = wrong_page(e.get("page"), pmap, rng)
    if r:
        out.append(dict(base, kind="wrong_page", planted_claim=text, planted_page=r[0], note=r[1]))
    return out


def plant_id(prefix, p):
    """A stable ID from the plant's content, so editing other entries never renumbers it."""
    key = "|".join(str(p.get(k)) for k in ("paper", "kind", "correct_claim", "planted_claim",
                                           "planted_page", "planted_value", "planted_unit"))
    return prefix + "-" + hashlib.sha256(key.encode()).hexdigest()[:8]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--seed", type=int, default=20261028)
    ap.add_argument("--per-statement", type=int, default=2,
                    help="plants kept per statement, chosen at random from those possible")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args(argv)
    problems = []
    pp = cc.papers(problems)
    gold = cc.load_yaml("gold_statements.yaml", problems)
    if problems:
        print("\n".join(problems))
        return 1
    rng = random.Random(args.seed)
    plants = []
    for i, e in enumerate(gold.get("statements") or [], 1):
        if not e.get("statement") or e.get("page") is None:
            continue
        cands = plants_for(i, e, pp, rng)
        rng.shuffle(cands)
        cands.sort(key=lambda p: p["kind"] == "wrong_page")   # content kinds first; pages fill in
        plants += cands[: args.per_statement]
    for p in plants:
        p["id"] = plant_id("S", p)
        p["approved"] = None
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    head = ("# GENERATED by scripts/literature_review/make_plants.py "
            f"(seed {args.seed}, {args.per_statement} per statement) from gold_statements.yaml.\n"
            "# Mechanical planted errors for the core calibration (WORKFLOW.md §5). Review each one:\n"
            "# set approved to true if the planted version is wrong in exactly the way the note says\n"
            "# and is not accidentally still true, false otherwise. Rerunning the script overwrites\n"
            "# this file, so record approvals only once the gold statements are final.\n")
    with open(args.out, "w") as f:
        f.write(head)
        yaml.safe_dump({"plants": plants}, f, sort_keys=False, allow_unicode=True, width=100)
    kinds = {}
    for p in plants:
        kinds[p["kind"]] = kinds.get(p["kind"], 0) + 1
    print(f"{len(plants)} plants -> {os.path.relpath(args.out, cc.ROOT)}: "
          + ", ".join(f"{k} {v}" for k, v in sorted(kinds.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
