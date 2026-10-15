"""Flag outline bullets whose quoted draft text is no longer in the draft.
Usage: python3 -I outline_stale_check.py <outline.tex> <draft.tex> <out.json>
GONE = every quoted fragment missing (likely addressed), PARTIAL = some, OPEN = all present.
Math and citations are stripped before matching, so confirm GONE/PARTIAL by eye."""
import re,sys,json
O,D=sys.argv[1],sys.argv[2]
def norm(t):
    t=re.sub(r'(?<!\\)%.*','',t)
    t=re.sub(r'\\(citep|citet|cref|ref|label|eqref)\{[^}]*\}',' ',t)
    t=re.sub(r'\$[^$]*\$',' ',t)
    t=re.sub(r'\\[a-zA-Z]+\*?',' ',t)
    t=re.sub(r'[{}~\\]',' ',t)
    t=re.sub(r'[^a-z0-9 ]',' ',t.lower())
    return re.sub(r'\s+',' ',t).strip()
draft=norm(open(D).read())
lines=open(O).read().split('\n')
start=next(i for i,l in enumerate(lines) if 'end literature check' in l)
items=[];cur=None
for i in range(start,len(lines)):
    l=lines[i]
    if re.match(r'\s*\\item\b',l) or re.match(r'\s*\\Pg\{',l) or re.match(r'\s*\\(sub)*section',l) or re.match(r'\s*\\end\{itemize\}',l):
        if cur: items.append(cur)
        cur={'line':i+1,'text':l} if re.match(r'\s*\\item\b',l) else None
        if re.match(r'\s*\\Pg\{',l): items.append({'line':i+1,'text':l,'pg':True})
    elif cur: cur['text']+='\n'+l
if cur: items.append(cur)
out=[]
for it in items:
    if it.get('pg'): out.append(('PG',it['line'],it['text'][:110])); continue
    qs=re.findall(r"``(.+?)''",it['text'],re.S)
    frags=[]
    for q in qs:
        for f in re.split(r'\\ldots|\[[^\]]*\]',q):
            n=norm(f)
            if len(n)>=12: frags.append(n)
    if not frags: st='NOQUOTE'
    else:
        hit=[f in draft for f in frags]
        st='OPEN' if all(hit) else ('GONE' if not any(hit) else 'PARTIAL')
    kind=re.search(r'\\t[CGTA]|\\(Have|Miss|Stale|Cut|Note)',it['text'])
    out.append((st,it['line'],(kind.group(0) if kind else '')+' '+re.sub(r'\s+',' ',it['text'])[:150]))
from collections import Counter
print(Counter(o[0] for o in out))
json.dump(out,open(sys.argv[3],'w'),indent=0)
