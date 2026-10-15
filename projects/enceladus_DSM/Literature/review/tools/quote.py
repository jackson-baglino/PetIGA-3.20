"""Find a phrase in Literature/text/<key>.txt; print PDF page and context.
Usage: python3 -I quote.py <key> "<phrase>" [<key> "<phrase>" ...]   (case-insensitive, whitespace-normalised)"""
import sys,re,os
T=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','text')
a=sys.argv[1:]
for key,ph in zip(a[::2],a[1::2]):
    s=open(os.path.join(T,key+'.txt'),encoding='utf-8',errors='replace').read()
    pages=re.split(r'\n=== page (\d+) ===\n',s)
    hits=0
    for i in range(1,len(pages),2):
        txt=re.sub(r'\s+',' ',pages[i+1]); 
        for m in re.finditer(re.escape(re.sub(r'\s+',' ',ph)),txt,re.I):
            print(f'[{key} p.{pages[i]}] ...{txt[max(0,m.start()-220):m.end()+220]}...'); hits+=1
            if hits>=2: break
        if hits>=2: break
    if not hits: print(f'[{key}] NOT FOUND: {ph}')
    print()
