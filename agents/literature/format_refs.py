import sys, json, urllib.request, urllib.parse, time
ids=[l.strip() for l in open(sys.argv[1]) if l.strip()]
out=[]
for i in range(0,len(ids),150):
    chunk=ids[i:i+150]
    data=urllib.parse.urlencode({"db":"pubmed","retmode":"json","id":",".join(chunk)}).encode()
    for _t in range(6):
        try:
            j=json.loads(urllib.request.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",data=data,timeout=60).read()); break
        except Exception: time.sleep(5*(_t+1))
    for p in chunk:
        d=j["result"][p]
        au=[a["name"] for a in d.get("authors",[]) if a.get("authtype")=="Author"] or [a["name"] for a in d.get("authors",[])]
        aus=", ".join(au[:3])+(" et al." if len(au)>3 else "")
        doi=next((a["value"] for a in d.get("articleids",[]) if a["idtype"]=="doi"),"")
        yr=d.get("pubdate","")[:4]
        pages=d.get("pages","")
        vol=d.get("volume","")
        loc=f"{vol}:{pages}" if pages else vol
        out.append((au[0].split()[0] if au else "", f"- {aus} {yr}. {d.get('title','').rstrip('.')}. *{d.get('source','')}* {loc}. doi:[{doi}](https://doi.org/{doi}). PMID:{p}"))
    time.sleep(0.5)
import unicodedata
def key(s): return unicodedata.normalize("NFKD",s[0]).encode("ascii","ignore").decode().lower()
for k,line in sorted(out,key=lambda x:(key(x),x[1])): print(line)
