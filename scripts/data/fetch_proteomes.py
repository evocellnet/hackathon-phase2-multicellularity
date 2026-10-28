#!/usr/bin/env python3
"""Step 2a: download the raw proteome files for the selected species.

  UniProt (release pinned: 2026_03, checked against current_release/relnotes.txt):
      .../reference_proteomes/<Division>/<UP>/<UP>_<taxid>.fasta.gz   (canonical set)
      .../reference_proteomes/<Division>/<UP>/<UP>_<taxid>.gene2acc.gz (gene -> accession)
  NCBI (Datasets v2 REST): genome/accession/<acc>/download?include_annotation_type=PROT_FASTA,GENOME_GFF
  EukProt v3 / LukProt v1.5.1: local archives (downloaded by the PI), only located here.

Writes data/raw/download_manifest.tsv (code, source, url, local_path, md5, bytes, date).
Re-runnable: files already present (and non-empty) are not downloaded again.
Usage: fetch_proteomes.py [--threads 6]
"""
import argparse, datetime, glob, hashlib, os, subprocess, sys, time, zipfile
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
RAW = f'{P}/data/raw'
UNIPROT_RELEASE = '2026_03'
UP_BASE = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes'
NCBI_DL = ('https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download'
           '?include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_GFF')


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def curl(url, out, tries=5):
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return
    tmp = out + '.part'
    for i in range(tries):
        r = subprocess.run(['curl', '-sSfL', '--retry', '3', '-o', tmp, url], capture_output=True, text=True)
        if r.returncode == 0 and os.path.getsize(tmp) > 0:
            os.replace(tmp, out)
            return
        time.sleep(5 * (i + 1))
    raise RuntimeError(f'download failed: {url}: {r.stderr.strip()}')


def check_release():
    r = subprocess.run(['curl', '-sSfL', f'{UP_BASE}/../../relnotes.txt'], capture_output=True, text=True)
    first = r.stdout.splitlines()[0] if r.stdout else ''
    if UNIPROT_RELEASE not in first:
        sys.exit(f'UniProt current_release is "{first}", expected {UNIPROT_RELEASE}: per-proteome files of the '
                 f'pinned release are no longer on the FTP; stop and decide.')
    print('UniProt release check:', first.strip())


def up_taxid():
    m = {}
    with open(f'{RAW}/uniprot/STATS_{UNIPROT_RELEASE}.txt') as fh:
        for line in fh:
            if line.startswith('UP'):
                f = line.split('\t')
                m[f[0]] = f[1].strip()
    return m


def job(row, kingdom, uptax):
    code, src = row['code'], row['source']
    db, acc = src.split(':', 1)
    out = []
    if db == 'UniProt':
        tx = uptax[acc]
        div = {'euk': 'Eukaryota'}.get(kingdom) or ('Archaea' if row['lineage_group'].startswith('Archaea') else 'Bacteria')
        d = f'{RAW}/uniprot/proteomes'
        os.makedirs(d, exist_ok=True)
        for suf in ('fasta.gz', 'gene2acc.gz'):
            url = f'{UP_BASE}/{div}/{acc}/{acc}_{tx}.{suf}'
            lp = f'{d}/{acc}_{tx}.{suf}'
            curl(url, lp)
            out.append((code, src, url, lp))
    elif db == 'NCBI':
        d = f'{RAW}/ncbi/{acc}'
        os.makedirs(d, exist_ok=True)
        url = NCBI_DL.format(acc=acc)
        lp = f'{d}/{acc}.zip'
        curl(url, lp)
        with zipfile.ZipFile(lp) as z:
            names = z.namelist()
            if not any(n.endswith('protein.faa') for n in names):
                raise RuntimeError(f'{acc}: no protein.faa in package ({names})')
            if not os.path.exists(f'{d}/ncbi_dataset'):
                z.extractall(d)
        out.append((code, src, url, lp))
    elif db == 'EukProt':
        fs = glob.glob(f'{RAW}/eukprot/proteins/{acc}_*.fasta')
        if len(fs) != 1:
            raise RuntimeError(f'{acc}: EukProt file not found ({fs})')
        out.append((code, src, 'https://ndownloader.figshare.com/files/34434377 (EukProt_proteins.v03.2021_11_22.tgz)', fs[0]))
    elif db == 'LukProt':
        fs = [f for f in glob.glob(f'{RAW}/lukprot/**/{acc}*', recursive=True)
              if f.endswith(('.fasta', '.fa', '.faa', '.fas'))]
        if len(fs) != 1:
            raise RuntimeError(f'{acc}: LukProt file not found ({fs}); is the 7z extracted?')
        out.append((code, src, 'https://zenodo.org/records/13829058/files/LukProt_v1.5.1_single_species_FASTA.7z', fs[0]))
    else:
        raise RuntimeError(f'unknown source {src}')
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--threads', type=int, default=6)
    ap.add_argument('--only', choices=['UniProt', 'NCBI', 'EukProt', 'LukProt'], nargs='*')
    a = ap.parse_args()
    check_release()
    uptax = up_taxid()
    rows = []
    for kingdom, f in (('euk', 'species_final.tsv'), ('prok', 'species_final_prok.tsv')):
        d = pd.read_csv(f'{P}/data/{f}', sep='\t', dtype=str)
        rows += [(r, kingdom) for _, r in d.iterrows()]
    if a.only:
        rows = [(r, k) for r, k in rows if r['source'].split(':')[0] in a.only]
    res, errs = [], []
    with ThreadPoolExecutor(a.threads) as ex:
        futs = {ex.submit(job, r, k, uptax): r['code'] for r, k in rows}
        for fu, code in futs.items():
            try:
                res += fu.result()
            except Exception as e:  # report, continue
                errs.append((code, str(e)))
    today = datetime.date.today().isoformat()
    man_path = f'{RAW}/download_manifest.tsv'
    old = pd.read_csv(man_path, sep='\t', dtype=str) if os.path.exists(man_path) else None
    new = pd.DataFrame(res, columns=['code', 'source', 'url', 'local_path'])
    new['md5'] = [md5(p) for p in new['local_path']]
    new['bytes'] = [os.path.getsize(p) for p in new['local_path']]
    new['date'] = today
    new['local_path'] = [os.path.relpath(p, P) for p in new['local_path']]
    if old is not None:
        # keep the original download date for unchanged files
        k = old.set_index(['local_path', 'md5'])['date'].to_dict()
        new['date'] = [k.get((p, m), today) for p, m in zip(new['local_path'], new['md5'])]
        old = old[~old['local_path'].isin(new['local_path'])]
        new = pd.concat([old, new])
    new.sort_values(['code', 'local_path']).to_csv(man_path, sep='\t', index=False)
    print(f'manifest rows: {len(new)}; errors: {len(errs)}')
    for c, e in errs:
        print('ERROR', c, e)
    sys.exit(1 if errs else 0)


if __name__ == '__main__':
    main()
