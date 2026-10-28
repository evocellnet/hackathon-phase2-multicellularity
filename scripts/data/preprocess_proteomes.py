#!/usr/bin/env python3
"""Step 2b: build the final per-species proteome FASTA files for Ancestra.

For each selected species (data/species_final*.tsv, raw files located through
data/raw/download_manifest.tsv or the EukProt/LukProt directories):

 1. one protein per gene (source specific)
    - UniProt canonical FASTA (already gene-centric); additionally, accessions that
      share a gene group in *.gene2acc (column 3, '-' = unknown -> never grouped) are
      collapsed to the first-listed accession (UniProt's representative).
    - NCBI: protein -> gene through the GFF Parent chain (CDS -> mRNA -> gene);
      longest protein per gene.
    - EukProt/LukProt genome sets: gene key parsed from the original source ID when
      it encodes isoforms (gene=..., scgid..., .tN, -RA, -mRNA-N, ...); longest per gene.
    - transcriptomes (and EukProt "genome" sets that are translated transcriptomes,
      see TRANSCRIPT_LIKE): longest ORF per Trinity gene (header field before '~~'),
      then CD-HIT -c 0.95 -G 0 -aS 0.9 -n 5 (local coverage: contained fragments join the cluster
      of their longer copy; representative = longest of each cluster).
 2. cleaning: terminal '*' (or '.') trimmed; internal stops and every letter outside
    the 20 standard amino acids replaced by X; '-' and whitespace removed;
    sequences < 50 aa removed; exact duplicate sequences dropped (first kept).
 3. IDs: <code>_<original id with characters outside [A-Za-z0-9_.-] replaced by '_'>.

Outputs: data/proteomes/{euk,prok}/<code>.fa
         data/proteomes/idmaps/<code>.tsv.gz   (original_id, new_id, gene_key, length, status)
         data/proteomes/stats/<code>.json
Usage:   preprocess_proteomes.py [--jobs 5] [--cdhit-threads 4] [--codes A B ...] [--force]
"""
import argparse, glob, gzip, hashlib, json, os, re, shutil, subprocess, sys, tempfile, time
from collections import OrderedDict, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
RAW = f'{P}/data/raw'
OUT = f'{P}/data/proteomes'
BIN = f'{P}/conda_envs/busco/bin'
MINLEN = 50
STD = set('ACDEFGHIKLMNPQRSTVWY')
TRANSCRIPT_LIKE = {'EP00667'}   # Euglena gracilis: EukProt "genome" set = translated transcriptome (EG_transcript_N)
SCRATCH = os.environ.get('PREPROC_TMP', f'{P}/work/preprocess_tmp')


# ----------------------------------------------------------------------------- io
def read_fasta(path):
    op = gzip.open if path.endswith('.gz') else open
    hdr, seq = None, []
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None:
                    yield hdr, ''.join(seq)
                hdr, seq = line[1:].rstrip('\n'), []
            else:
                seq.append(line.strip())
    if hdr is not None:
        yield hdr, ''.join(seq)


def clean(seq):
    s = re.sub(r'[\s\-]', '', seq).upper()
    s = s.rstrip('*.')
    return ''.join(c if c in STD else 'X' for c in s)


def sanitise(x):
    return re.sub(r'[^A-Za-z0-9_.\-]', '_', x)


# ------------------------------------------------------------------- gene keys
def eukprot_gene_key(hdr):
    """Gene key for EukProt/LukProt genome sets (None = protein is its own gene)."""
    m = re.search(r'\bgene=(\S+)', hdr)
    if m:
        return m.group(1)
    m = re.search(r'\bscgid(\d+)', hdr)
    if m:
        return 'scgid' + m.group(1)
    toks = hdr.split()
    if len(toks) < 2:
        return None
    t = toks[1]
    if t.startswith('jgi|'):
        t = t.split('|')[-1]
    if re.match(r'^[A-Z]{2}_\d+\.\d+$', t):          # RefSeq XP_/NP_ accession: no gene info
        return None
    for pat, rep in ((r'-mRNA-\d+$', ''), (r'\.mRNA\.\d+$', ''), (r'-R[A-Z]{1,2}(\.\d+)?$', ''),
                     (r'\.t\d+(-p\d+)?$', ''), (r'(\d)[Tt]\d*$', r'\1'), (r'(g\d+)[a-z]$', r'\1'),
                     (r'(_\d+)\.\d+$', r'\1')):
        t2 = re.sub(pat, rep, t)
        if t2 != t:
            return t2
    return t


def trinity_gene_key(hdr):
    """Gene key for transcriptome sets (None = no gene information; CD-HIT only).
    TransDecoder/Trinity  '... TRINITY_DN0_c0_g1~~TRINITY_DN0_c0_g1_i1.p1 ...' -> TRINITY_DN0_c0_g1
                          '... GENE.c0_g1_i1~~c0_g1_i1.p1 ...'                 -> c0_g1
    old Trinity           '... comp100027_c0_seq1:141-917(+)'                  -> comp100027_c0
    SOAPdenovo-Trans      'Unigene10000_Nuclearia_a|m.30543 ...'               -> Unigene10000_Nuclearia_a
    MMETSP (CAMPEP_...), '@'-numbered sets, plain contig IDs: no gene information."""
    toks = hdr.split()
    for t in toks[1:3]:
        if '~~' in t:
            g = t.split('~~')[0]
            if g.startswith('GENE.'):
                g = g[5:]
            return re.sub(r'((?:^|_)c\d+_g\d+)_i\d+$', r'\1', g)
    m = re.search(r'\b(comp\d+_c\d+)_seq\d+', hdr)
    if m:
        return m.group(1)
    if len(toks) > 1:
        m = re.match(r'(Unigene\d+[^|\s]*)\|', toks[1])
        if m:
            return m.group(1)
        m = re.match(r'(TRINITY_\S+?_g\d+)_i\d+', toks[1])
        if m:
            return m.group(1)
        m = re.match(r'(\S*?DN\d+_c\d+_g\d+)_i\d+', toks[1])   # e.g. Partr_v1_DN10022_c0_g2_i1_m38271
        if m:
            return m.group(1)
    return None


def uniprot_groups(g2a_path):
    acc2grp = {}
    if not os.path.exists(g2a_path):
        return acc2grp
    with gzip.open(g2a_path, 'rt') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 3 or re.search(r'-\d+$', f[1]):
                continue   # isoform rows
            grp = f[2] if f[2] != '-' else None
            acc2grp.setdefault(f[1], grp)
    return acc2grp


def ncbi_gene_map(gff):
    parent, typ, cds2prot = {}, {}, {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
            i = a.get('ID')
            if i:
                typ[i] = f[2]
                if 'Parent' in a:
                    parent[i] = a['Parent'].split(',')[0]
            if f[2] == 'CDS' and 'protein_id' in a:
                cds2prot[a['protein_id']] = (a.get('Parent'), a.get('locus_tag'))
    prot2gene = {}
    for prot, (par, lt) in cds2prot.items():
        cur, n = par, 0
        while cur is not None and typ.get(cur) != 'gene' and cur in parent and n < 10:
            cur, n = parent[cur], n + 1
        prot2gene[prot] = cur if cur is not None else (lt or prot)
    return prot2gene


# ------------------------------------------------------------------- per species
def run_cdhit(records, threads, workdir):
    """records: list of (orig_id, seq). Returns dict orig_id -> representative orig_id."""
    inp, out = f'{workdir}/in.fa', f'{workdir}/out.fa'
    idx = {}
    with open(inp, 'w') as fh:
        for k, (oid, s) in enumerate(records):
            idx[f's{k}'] = oid
            fh.write(f'>s{k}\n{s}\n')
    cmd = [f'{BIN}/cd-hit', '-i', inp, '-o', out, '-c', '0.95', '-G', '0', '-aS', '0.9', '-n', '5',
           '-T', str(threads), '-M', '0', '-d', '0']
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('cd-hit failed: ' + r.stderr[-2000:])
    rep_of, members, rep = {}, [], None
    with open(out + '.clstr') as fh:
        for line in fh:
            if line.startswith('>Cluster'):
                for m in members:
                    rep_of[m] = rep
                members, rep = [], None
            else:
                sid = line.split('>')[1].split('...')[0]
                members.append(idx[sid])
                if line.rstrip().endswith('*'):
                    rep = idx[sid]
    for m in members:
        rep_of[m] = rep
    return rep_of, ' '.join(cmd[1:]).replace(inp, 'IN').replace(out, 'OUT')


def process(row, kingdom, manifest, cdhit_threads, force):
    code, src, dtype = row['code'], row['source'], row['data_type']
    db, acc = src.split(':', 1)
    fa_out = f'{OUT}/{kingdom}/{code}.fa'
    st_out = f'{OUT}/stats/{code}.json'
    if not force and os.path.exists(fa_out) and os.path.exists(st_out):
        return json.load(open(st_out))
    t0 = time.time()
    method = []
    gene_of = {}
    recs = []   # (orig_id, raw_header, seq)
    if db == 'UniProt':
        fa = [p for p in manifest.get(code, []) if p.endswith('.fasta.gz')][0]
        g2a = fa.replace('.fasta.gz', '.gene2acc.gz')
        acc2grp = uniprot_groups(g2a)
        for h, s in read_fasta(fa):
            oid = h.split()[0].split('|')[1] if '|' in h.split()[0] else h.split()[0]
            recs.append((oid, h, s))
            g = acc2grp.get(oid)
            gene_of[oid] = ('G:' + g) if g else None
        method.append('UniProt canonical FASTA; gene2acc group collapse (first listed accession kept)')
        src_file = fa
    elif db == 'NCBI':
        d = f'{RAW}/ncbi/{acc}/ncbi_dataset/data/{acc}'
        p2g = ncbi_gene_map(f'{d}/genomic.gff')
        for h, s in read_fasta(f'{d}/protein.faa'):
            oid = h.split()[0]
            recs.append((oid, h, s))
            gene_of[oid] = p2g.get(oid)
        method.append('NCBI protein.faa; longest protein per gene (GFF CDS->mRNA->gene)')
        src_file = f'{d}/protein.faa'
    elif db in ('EukProt', 'LukProt'):
        if db == 'EukProt':
            src_file = glob.glob(f'{RAW}/eukprot/proteins/{acc}_*.fasta')[0]
        else:
            src_file = [p for p in manifest.get(code, [])][0]
        transcript_like = dtype != 'genome' or acc in TRANSCRIPT_LIKE
        for h, s in read_fasta(src_file):
            oid = h.split()[0]
            recs.append((oid, h, s))
            gene_of[oid] = trinity_gene_key(h) if transcript_like else eukprot_gene_key(h)
        method.append(f'{db} set ({"transcriptome-like" if transcript_like else "genome"}); '
                      + ('longest ORF per Trinity gene' if transcript_like else 'longest protein per parsed gene ID'))
    else:
        raise RuntimeError(src)

    n_raw = len(recs)
    status = OrderedDict()
    seqs = {}
    for oid, h, s in recs:
        if oid in seqs:
            raise RuntimeError(f'{code}: duplicated original ID {oid}')
        seqs[oid] = clean(s)
    # 1. one protein per gene
    best = {}
    for oid, h, s in recs:
        g = gene_of.get(oid)
        if g is None:
            continue
        L = len(seqs[oid])
        if g not in best or L > len(seqs[best[g]]):
            best[g] = oid
    keep_gene = set(best.values())
    for oid, h, s in recs:
        g = gene_of.get(oid)
        status[oid] = 'kept' if (g is None or oid in keep_gene) else f'dropped_isoform_of:{best[g]}'
    n_genes = sum(1 for v in status.values() if v == 'kept')
    # 2. length + exact duplicates
    seen = {}
    for oid in status:
        if status[oid] != 'kept':
            continue
        s = seqs[oid]
        if len(s) < MINLEN:
            status[oid] = 'dropped_short'
            continue
        h = hashlib.md5(s.encode()).hexdigest()
        if h in seen:
            status[oid] = f'dropped_identical_to:{seen[h]}'
        else:
            seen[h] = oid
    n_short = sum(1 for v in status.values() if v == 'dropped_short')
    n_dup = sum(1 for v in status.values() if v.startswith('dropped_identical'))
    # 3. CD-HIT for transcriptome-like sets
    n_cd = 0
    cd_cmd = ''
    transcript_like = db in ('EukProt', 'LukProt') and (dtype != 'genome' or acc in TRANSCRIPT_LIKE)
    if transcript_like:
        os.makedirs(SCRATCH, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=SCRATCH) as wd:
            kept = [(o, seqs[o]) for o, v in status.items() if v == 'kept']
            kept.sort(key=lambda x: -len(x[1]))
            rep_of, cd_cmd = run_cdhit(kept, cdhit_threads, wd)
        for o, _ in kept:
            r = rep_of.get(o, o)
            if r != o:
                status[o] = f'dropped_cdhit_member_of:{r}'
                n_cd += 1
        method.append('CD-HIT ' + cd_cmd)
    # 4. write
    os.makedirs(os.path.dirname(fa_out), exist_ok=True)
    os.makedirs(f'{OUT}/idmaps', exist_ok=True)
    os.makedirs(f'{OUT}/stats', exist_ok=True)
    newid, used = {}, set()
    for oid in status:
        n = f'{code}_{sanitise(oid)}'
        k = 1
        while n in used:
            k += 1
            n = f'{code}_{sanitise(oid)}_{k}'
        used.add(n)
        newid[oid] = n
    tmp = fa_out + '.tmp'
    n_final = 0
    with open(tmp, 'w') as fh:
        for oid, v in status.items():
            if v == 'kept':
                s = seqs[oid]
                fh.write(f'>{newid[oid]}\n')
                for i in range(0, len(s), 60):
                    fh.write(s[i:i + 60] + '\n')
                n_final += 1
    os.replace(tmp, fa_out)
    with gzip.open(f'{OUT}/idmaps/{code}.tsv.gz', 'wt') as fh:
        fh.write('original_id\tnew_id\tgene_key\tlength\tstatus\n')
        for oid, v in status.items():
            fh.write(f'{oid}\t{newid[oid] if v == "kept" else ""}\t{gene_of.get(oid) or ""}\t{len(seqs[oid])}\t{v}\n')
    h = hashlib.md5(open(fa_out, 'rb').read()).hexdigest()
    st = dict(code=code, kingdom=kingdom, source=src, data_type=dtype, source_file=os.path.relpath(src_file, P),
              n_raw=n_raw, n_after_gene_collapse=n_genes, n_dropped_short=n_short, n_dropped_identical=n_dup,
              n_dropped_cdhit=n_cd, n_final=n_final, md5=h, method='; '.join(method),
              seconds=round(time.time() - t0, 1))
    json.dump(st, open(st_out, 'w'), indent=1)
    return st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=5)
    ap.add_argument('--cdhit-threads', type=int, default=4)
    ap.add_argument('--codes', nargs='*')
    ap.add_argument('--sources', nargs='*', help='restrict to source dbs (UniProt NCBI EukProt LukProt)')
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    man = pd.read_csv(f'{RAW}/download_manifest.tsv', sep='\t', dtype=str)
    manifest = defaultdict(list)
    for c, lp in zip(man['code'], man['local_path']):
        manifest[c].append(f'{P}/{lp}')
    rows = []
    for kingdom, f in (('euk', 'species_final.tsv'), ('prok', 'species_final_prok.tsv')):
        d = pd.read_csv(f'{P}/data/{f}', sep='\t', dtype=str)
        rows += [(r, kingdom) for _, r in d.iterrows()]
    # FASTA files of species that are no longer selected are moved out of the Ancestra input directories
    sel = {(k, r['code']) for r, k in rows}
    for k in ('euk', 'prok'):
        for f in glob.glob(f'{OUT}/{k}/*.fa'):
            c = os.path.basename(f)[:-3]
            if (k, c) not in sel:
                os.makedirs(f'{OUT}/_dropped', exist_ok=True)
                shutil.move(f, f'{OUT}/_dropped/{c}.fa')
                print(f'moved unselected {c}.fa to data/proteomes/_dropped/', flush=True)
    if a.codes:
        rows = [(r, k) for r, k in rows if r['code'] in a.codes]
    if a.sources:
        rows = [(r, k) for r, k in rows if r['source'].split(':')[0] in a.sources]
    errs = 0
    with ProcessPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(process, r, k, manifest, a.cdhit_threads, a.force): r['code'] for r, k in rows}
        for fu in as_completed(futs):
            c = futs[fu]
            try:
                st = fu.result()
                print(f"{c}\traw={st['n_raw']}\tgenes={st['n_after_gene_collapse']}\tshort={st['n_dropped_short']}"
                      f"\tident={st['n_dropped_identical']}\tcdhit={st['n_dropped_cdhit']}\tfinal={st['n_final']}", flush=True)
            except Exception as e:
                errs += 1
                print(f'ERROR\t{c}\t{e}', flush=True)
    sys.exit(1 if errs else 0)


if __name__ == '__main__':
    main()
