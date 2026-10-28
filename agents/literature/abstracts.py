import sys, urllib.request, re, time
ids=sys.argv[1].split(",")
for i in range(0,len(ids),10):
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&rettype=abstract&retmode=xml&id="+",".join(ids[i:i+10])
    x=urllib.request.urlopen(url,timeout=60).read().decode()
    for art in x.split("<PubmedArticle>")[1:]:
        pm=re.search(r"<PMID[^>]*>(\d+)</PMID>",art).group(1)
        ab=" ".join(re.sub(r"<[^>]+>","",a) for a in re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>",art,re.S))
        print(f"## {pm}: {ab}\n")
    time.sleep(0.5)
