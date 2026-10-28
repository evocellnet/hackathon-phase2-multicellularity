#!/usr/bin/env python3
"""Verify references against PubMed (E-utilities) or Crossref.
Input lines (stdin): KEY<TAB>QUERY where QUERY is either
  pmid:NNNN | doi:10.xxx | a PubMed search term (first hit taken, up to 3 shown)
Output: KEY | PMID | first author | year | title | journal vol:pages | DOI
"""
import os, sys, json, time, urllib.request, urllib.parse
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(url):
    for i in range(4):
        try:
            with urllib.request.urlopen(url, timeout=30) as r: return json.loads(r.read())
        except Exception as e:
            time.sleep(1.5*(i+1))
    return None
def summ(ids):
    if not ids: return []
    j=get(E+"esummary.fcgi?db=pubmed&retmode=json&id="+",".join(ids))
    out=[]
    for i in ids:
        d=j["result"].get(i,{}) if j else {}
        doi=next((a["value"] for a in d.get("articleids",[]) if a["idtype"]=="doi"),"")
        au=d.get("authors",[]); fa=au[0]["name"] if au else ""
        out.append(f"PMID:{i} | {fa}{' et al.' if len(au)>1 else ''} | {d.get('pubdate','')[:4]} | {d.get('title','')} | {d.get('source','')} {d.get('volume','')}:{d.get('pages','')} | doi:{doi}")
    return out
def crossref(doi):
    j=get("https://api.crossref.org/works/"+urllib.parse.quote(doi))
    if not j: return ["CROSSREF-NOT-FOUND "+doi]
    m=j["message"]; au=m.get("author",[]); fa=(au[0].get("family","")+" "+au[0].get("given","")[:1]) if au else ""
    yr=(m.get("issued",{}).get("date-parts",[[""]])[0][0])
    return [f"CROSSREF | {fa}{' et al.' if len(au)>1 else ''} | {yr} | {(m.get('title') or [''])[0]} | {(m.get('container-title') or [''])[0]} {m.get('volume','')}:{m.get('page','')} | doi:{doi}"]
for line in sys.stdin:
    line=line.rstrip("\n")
    if not line.strip() or "\t" not in line: continue
    k,q=line.split("\t",1)
    if q.startswith("pmid:"): res=summ([q[5:]])
    elif q.startswith("doi:"): res=crossref(q[4:])
    else:
        j=get(E+"esearch.fcgi?db=pubmed&retmode=json&retmax="+os.environ.get("RETMAX","3")+"&term="+urllib.parse.quote(q))
        ids=j["esearchresult"]["idlist"] if j else []
        res=summ(ids) if ids else ["NOT-FOUND"]
    for r in res: print(f"{k}\t{r}")
    time.sleep(0.4)
