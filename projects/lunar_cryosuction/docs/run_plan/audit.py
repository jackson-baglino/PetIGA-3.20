"""List text outside Word equations that still holds symbols or relations."""
import re, sys, zipfile
from lxml import etree
W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
M = "{http://schemas.openxmlformats.org/officeDocument/2006/math}"
x = etree.fromstring(zipfile.ZipFile(sys.argv[1]).read("word/document.xml"))
pat = re.compile(r"[Ͱ-Ͽ⁰¹²³⁴⁵⁶⁷⁸⁹⁻×≈≤≥±∓∫∇·]|(?<![\w/.-])(G|H|h|R|T0|d0|l|v)\s*[=<>≈]|\b(d0|T0|R0|h_open|D_v|w_t|w_c|v_n|α_c)\b")
n = 0; nm = 0
for p in x.iter(W + "p"):
    nm += len(list(p.iter(M + "oMath")))
    txt = "".join(t.text or "" for t in p.iter(W + "t") if not any(a.tag == M + "oMath" for a in t.iterancestors()))
    if "--" in txt or "venv_" in txt:   # Appendix B command lines
        continue
    for m in pat.finditer(txt):
        n += 1
        print(repr(txt[max(0, m.start() - 40): m.end() + 30]))
print("equations:", nm, " leftover hits:", n)
