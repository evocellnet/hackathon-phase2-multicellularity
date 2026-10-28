#!/usr/bin/env python3
"""Summarise the reference AFDB download: data/structure/ref_download_summary.tsv, ref_missing.tsv,
and refresh the auto-generated block in docs/STRUCTURE_DATA.md (between AUTO-BEGIN/AUTO-END markers)."""
import os, csv, re, glob, datetime
P = os.environ.get('P', os.getcwd())
sp = list(csv.DictReader(open(f'{P}/data/structure_species.tsv'), delimiter='\t'))
codes = []; 
for r in sp:
    if r['code'] not in codes: codes.append(r['code'])
last = {}
for r in csv.DictReader(open(f'{P}/runs/structure/afdb_ref/download_log.tsv'), delimiter='\t'):
    last[r['acc']] = r          # last row per accession wins (resumed runs re-log cached files)
os.makedirs(f'{P}/data/structure', exist_ok=True)
rows = []; miss = []
for c in codes:
    accs = [l[1:].strip().split('_', 1)[1] for l in open(f'{P}/data/proteomes/euk/{c}.fa') if l[0] == '>']
    ok = [last[a] for a in accs if a in last and last[a]['status'] in ('ok', 'cached')]
    m = [a for a in accs if a in last and last[a]['status'] == 'missing']
    err = [a for a in accs if a in last and last[a]['status'] == 'error']
    pend = len(accs) - len(ok) - len(m) - len(err)
    pl = [float(r['mean_plddt']) for r in ok if r['mean_plddt']]
    rows.append([c, len(accs), len(ok), len(m), len(err), pend, round(sum(pl)/len(pl), 1) if pl else '',
                 sum(1 for x in pl if x >= 70), round(sum(int(r['bytes']) for r in ok)/1e9, 2)])
    miss += [(c, a, 'missing') for a in m] + [(c, a, 'error') for a in err]
H = ['code', 'n_proteins', 'n_downloaded', 'n_missing_404', 'n_error', 'n_pending', 'mean_plddt', 'n_plddt_ge70', 'GB']
with open(f'{P}/data/structure/ref_download_summary.tsv', 'w') as w:
    w.write('\t'.join(H) + '\n'); [w.write('\t'.join(map(str, r)) + '\n') for r in rows]
with open(f'{P}/data/structure/ref_missing.tsv', 'w') as w:
    w.write('code\tacc\tstatus\n'); [w.write('\t'.join(m) + '\n') for m in miss]
tot = [sum(r[i] for r in rows) for i in (1, 2, 3, 4, 5, 7)]
vers = sorted({r['version'] for r in last.values() if r['version'] not in ('', 'None')})
blk = ['| ' + ' | '.join(H) + ' |', '|' + '---|' * len(H)] + ['| ' + ' | '.join(map(str, r)) + ' |' for r in rows]
blk.append(f'| **total** | {tot[0]} | {tot[1]} | {tot[2]} | {tot[3]} | {tot[4]} | | {tot[5]} | {round(sum(r[8] for r in rows),2)} |')
blk.append(f'\nAFDB model version(s) seen: {", ".join(vers)}. Report generated {datetime.date.today()}.')
p = f'{P}/docs/STRUCTURE_DATA.md'
if os.path.exists(p):
    s = open(p).read()
    s = re.sub(r'(<!-- AUTO-BEGIN -->).*?(<!-- AUTO-END -->)', lambda m: m.group(1) + '\n' + '\n'.join(blk) + '\n' + m.group(2), s, flags=re.S)
    open(p, 'w').write(s)
print('\n'.join(blk))
