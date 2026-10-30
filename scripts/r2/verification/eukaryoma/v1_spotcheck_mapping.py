"""Task 1: spot-check 200 random protein mappings (150 MD5 + 50 RBH) by sequence
identity between our proteome sequence and the EukaryOMA API sequence."""
import gzip, json, random, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices

P = "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT = f"{P}/results/r2/verification/eukaryoma"
random.seed(20261030)
md5 = pd.read_csv(f"{P}/runs/eukaryoma_xval/protein_map.tsv.gz", sep="\t")
rbh = pd.read_csv(f"{P}/runs/eukaryoma_xval/protein_map_rbh.tsv.gz", sep="\t")
s = pd.concat([md5.sample(150, random_state=1).assign(src="md5"),
               rbh.sample(50, random_state=2).assign(src="rbh")], ignore_index=True)
smap = pd.read_csv(f"{P}/results/eukaryoma_xval/species_map.tsv", sep="\t")
code2x = dict(zip(smap.code, smap.eukaryoma_code))

# our sequences
need = {}
for c, i in zip(s.code, s.our_id):
    need.setdefault(c, set()).add(i)
ours = {}
for c, ids in need.items():
    cur = None
    with open(f"{P}/data/proteomes/euk/{c}.fa") as fh:
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                if cur in ids:
                    ours[cur] = []
                else:
                    cur = None
            elif cur:
                ours[cur].append(line.strip())
ours = {k: "".join(v).rstrip("*") for k, v in ours.items()}

def fetch(oid):
    for a in range(5):
        try:
            with urllib.request.urlopen(f"https://eukaryoma.unil.ch/api/protein/{oid}/", timeout=60) as r:
                return json.load(r)
        except Exception as e:
            time.sleep(2 + 3 * a)
    return None

with ThreadPoolExecutor(6) as ex:
    js = list(ex.map(fetch, s.oma_id))

al = Align.PairwiseAligner(mode="local", open_gap_score=-11, extend_gap_score=-1)
al.substitution_matrix = substitution_matrices.load("BLOSUM62")
rows = []
for r, j in zip(s.itertuples(), js):
    q = ours.get(r.our_id)
    if j is None or q is None:
        rows.append(dict(src=r.src, code=r.code, our_id=r.our_id, oma_id=r.oma_id, status="fetch_fail" if j is None else "our_seq_missing"))
        continue
    t = j["sequence"].rstrip("*")
    if q == t:
        pid, qcov, tcov = 100.0, 1.0, 1.0
    else:
        a = al.align(q, t)[0]
        ident = sum(1 for (qs, qe), (ts, te) in zip(*a.aligned) for k in range(qe - qs) if q[qs + k] == t[ts + k])
        alen = a.shape[1]
        qa = sum(qe - qs for qs, qe in a.aligned[0]); ta = sum(te - ts for ts, te in a.aligned[1])
        pid, qcov, tcov = 100 * ident / alen, qa / len(q), ta / len(t)
    rows.append(dict(src=r.src, code=r.code, our_id=r.our_id, oma_id=r.oma_id,
                     api_species_code=j["species"]["code"], expected_code=code2x.get(r.code),
                     identical=q == t, pident=round(pid, 2), qcov=round(qcov, 3), tcov=round(tcov, 3),
                     status="ok"))
df = pd.DataFrame(rows)
df["species_ok"] = df.api_species_code == df.expected_code
df["pass"] = df.species_ok & (df.identical | ((df.pident >= 95) & (df.qcov >= 0.5) & (df.tcov >= 0.5)))
df.to_csv(f"{OUT}/v1_mapping_spotcheck.tsv", sep="\t", index=False)
print(df.groupby("src").agg(n=("pass", "size"), passed=("pass", "sum"), identical=("identical", "sum"),
                             species_ok=("species_ok", "sum"), min_pid=("pident", "min")))
print(df[~df["pass"]])
