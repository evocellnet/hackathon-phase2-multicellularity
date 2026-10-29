#!/bin/bash
# Ab initio gene prediction for unannotated volvocine assemblies (no annotation exists on NCBI/UniProt/EukProt/LukProt).
# AUGUSTUS 3.5.0 (conda_envs/augustus) with the pre-trained Chlamydomonas reinhardtii model 'chlamy2011', no hints,
# one transcript per gene (--alternatives-from-sampling=false), protein output.
# Usage: augustus_predict.sh <GCA accession> <n_chunks>   (input: data/raw/ncbi_genomes/<acc>/.../*.fna)
set -euo pipefail
P=$(cd "$(dirname "$0")/../.." && pwd)
ACC=$1; N=${2:-5}
BIN=$P/conda_envs/augustus/bin
export AUGUSTUS_CONFIG_PATH=$P/conda_envs/augustus/config
FNA=$(ls $P/data/raw/ncbi_genomes/$ACC/ncbi_dataset/data/$ACC/*.fna)
OUT=$P/data/raw/augustus/$ACC; mkdir -p $OUT/chunks
# balanced split by length (greedy)
$P/conda_envs/ancestra/bin/python - "$FNA" "$OUT/chunks" "$N" <<'PY'
import sys
fa, od, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
recs, h, s = [], None, []
for l in open(fa):
    if l[0] == '>':
        if h: recs.append((h, ''.join(s)))
        h, s = l[1:].split()[0], []
    else: s.append(l.strip())
recs.append((h, ''.join(s)))
recs.sort(key=lambda x: -len(x[1]))
bins = [[0, []] for _ in range(n)]
for r in recs:
    b = min(bins, key=lambda x: x[0]); b[0] += len(r[1]); b[1].append(r)
for i, (L, rs) in enumerate(bins):
    with open(f'{od}/chunk{i}.fa', 'w') as fh:
        for h, s in rs: fh.write(f'>{h}\n{s}\n')
PY
for i in $(seq 0 $((N-1))); do
  [ -s $OUT/chunks/chunk$i.gff3.done ] && continue
  ( nice -n 10 $BIN/augustus --species=chlamy2011 --gff3=on --protein=on --codingseq=off --uniqueGeneId=true \
      --alternatives-from-sampling=false $OUT/chunks/chunk$i.fa > $OUT/chunks/chunk$i.gff3 2> $OUT/chunks/chunk$i.log \
      && echo ok > $OUT/chunks/chunk$i.gff3.done ) &
done
wait
# proteins from the GFF3 comments ('# protein sequence = [...]')
$P/conda_envs/ancestra/bin/python - $OUT <<'PY'
import sys, glob, re
od = sys.argv[1]
n = 0
with open(f'{od}/proteins.faa', 'w') as out:
    for g in sorted(glob.glob(f'{od}/chunks/chunk*.gff3')):
        ck = g.split('/')[-1].split('.')[0]
        tid, seq, inp = None, [], False
        for l in open(g):
            if '\ttranscript\t' in l or '\tmRNA\t' in l:
                tid = re.search(r'ID=([^;\s]+)', l).group(1)
            if l.startswith('# protein sequence = ['):
                inp = True; seq = [l.split('[', 1)[1].strip()]
            elif inp and l.startswith('# '):
                seq.append(l[2:].strip())
            if inp and seq and seq[-1].endswith(']'):
                s = ''.join(seq).rstrip(']')
                out.write(f'>{ck}_{tid}\n{s}\n'); n += 1; inp = False; seq = []
print('proteins', n)
PY
