#!/usr/bin/env python3
"""Step 3: BUSCO (v6.1.0, odb10 lineages, offline, protein mode) on every final proteome.

eukaryotes -> eukaryota_odb10; Bacteria -> bacteria_odb10; Archaea -> archaea_odb10.
Runs in runs/busco/<kingdom>/<code>/ (git-ignored); parsed summary ->
data/qc/busco_summary.tsv (code, lineage, C, S, D, F, M, n, busco_version, input_md5); full outputs
copied to data/busco_runs/<code>/.
Default 5 parallel jobs x 4 CPUs = 20 CPUs (project limit). Re-runnable: a finished run is reused only if the md5 of its input FASTA is unchanged.
Usage: run_busco.py [--jobs 5] [--cpu 4] [--codes ...]
"""
import argparse, glob, hashlib, json, os, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

P = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
BIN = f'{P}/conda_envs/busco/bin'
DL = f'{P}/data/raw/busco/busco_downloads'
ENV = dict(os.environ, PATH=f'{BIN}:' + os.environ.get('PATH', ''))


def lineage_for(kingdom, lineage_group):
    if kingdom == 'euk':
        return 'eukaryota_odb10'
    return 'archaea_odb10' if str(lineage_group).startswith('Archaea') else 'bacteria_odb10'


def parse(outdir):
    js = glob.glob(f'{outdir}/short_summary*.json')
    if not js:
        return None
    j = json.load(open(js[0]))
    r = j['results']
    return dict(C=float(r['Complete percentage']), S=float(r['Single copy percentage']),
                D=float(r['Multi copy percentage']), F=float(r['Fragmented percentage']),
                M=float(r['Missing percentage']), n=int(r['n_markers']),
                busco_version=j.get('parameters', {}).get('busco_version', j.get('versions', {}).get('busco', '')),
                lineage=j['lineage_dataset']['name'])


def run(code, kingdom, lineage, cpu):
    base = f'{P}/runs/busco/{kingdom}'
    os.makedirs(base, exist_ok=True)
    out = f'{base}/{code}'
    fa = f'{P}/data/proteomes/{kingdom}/{code}.fa'
    h = hashlib.md5(open(fa, 'rb').read()).hexdigest()
    stamp = f'{base}/{code}.input_md5'
    res = parse(out)
    if res and os.path.exists(stamp) and open(stamp).read().strip() == h:
        return code, dict(res, input_md5=h), 'cached'
    cmd = [f'{BIN}/busco', '-i', fa, '-m', 'proteins', '-l', f'{DL}/lineages/{lineage}', '--offline',
           '--download_path', DL, '-c', str(cpu), '-o', code, '--out_path', base, '-f', '-q']
    r = subprocess.run(cmd, capture_output=True, text=True, env=ENV)
    with open(f'{base}/{code}.log', 'w') as fh:
        fh.write(' '.join(cmd) + '\n' + r.stdout + r.stderr)
    res = parse(out)
    if res is None:
        raise RuntimeError(f'{code}: BUSCO failed (rc={r.returncode}); see {base}/{code}.log')
    open(stamp, 'w').write(h + '\n')
    return code, dict(res, input_md5=h), 'run'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=5)
    ap.add_argument('--cpu', type=int, default=4)
    ap.add_argument('--codes', nargs='*')
    a = ap.parse_args()
    assert a.jobs * a.cpu <= 20, 'CPU limit is 20'
    todo = []
    for kingdom, f in (('euk', 'species_final.tsv'), ('prok', 'species_final_prok.tsv')):
        d = pd.read_csv(f'{P}/data/{f}', sep='\t', dtype=str)
        todo += [(c, kingdom, lineage_for(kingdom, g)) for c, g in zip(d['code'], d['lineage_group'])]
    if a.codes:
        todo = [t for t in todo if t[0] in a.codes]
    rows, errs = [], 0
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(run, c, k, l, a.cpu) for c, k, l in todo]
        for fu in as_completed(futs):
            try:
                code, res, how = fu.result()
                rows.append(dict(code=code, **res))
                print(f"{code}\t{res['lineage']}\tC={res['C']}\tD={res['D']}\tF={res['F']}\tM={res['M']}\t{how}", flush=True)
            except Exception as e:
                errs += 1
                print('ERROR', e, flush=True)
    # keep the full BUSCO output of every selected species under data/busco_runs/<code>/ (single-copy marker
    # sequences run_*/busco_sequences/ and full_table.tsv are used for the concatenated species tree / RER)
    import shutil
    for c, k, l in todo:
        src, dst = f'{P}/runs/busco/{k}/{c}', f'{P}/data/busco_runs/{c}'
        if os.path.isdir(src):
            shutil.copytree(src, dst, dirs_exist_ok=True)
            shutil.copy(f'{P}/runs/busco/{k}/{c}.input_md5', f'{dst}/input_md5.txt')
    os.makedirs(f'{P}/data/qc', exist_ok=True)
    out = f'{P}/data/qc/busco_summary.tsv'
    new = pd.DataFrame(rows)
    if os.path.exists(out) and a.codes:
        old = pd.read_csv(out, sep='\t')
        new = pd.concat([old[~old['code'].isin(new['code'])], new])
    new.sort_values('code').to_csv(out, sep='\t', index=False)
    print(f'wrote {out} ({len(new)} rows); errors={errs}')
    sys.exit(1 if errs else 0)


if __name__ == '__main__':
    main()
