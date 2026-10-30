#!/bin/bash
# exhaustive (no k-mer prefilter) Foldseek 3Di+AA re-search of producer-"structurally new" families, full panel DB, filter outside clade later
set -euo pipefail
FS=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/conda_envs/ancestra/bin/foldseek
S=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/results/r2/verification/w5/work/spot
E=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/results/r2/verification/w5/work/exh
cd $E
mk() { nice -n 10 $FS base:createdb $1 $3 --shuffle 0 -v 1; nice -n 10 $FS base:createdb $2 $3_t3 --shuffle 0 -v 1
  mv $3_t3 $3_ss; mv $3_t3.index $3_ss.index; mv $3_t3.dbtype $3_ss.dbtype; }
for R in "$@"; do
  [ -s $R.m8 ] && continue
  python3 - <<PY
ids=set(open('$R.qids').read().split())
for ext in ('aa','3di'):
    out=open('$R.q.'+ext+'.fa','w'); keep=False
    for l in open('$S/$R.'+ext+'.fa'):
        if l[0]=='>': keep=l[1:].strip() in ids
        if keep: out.write(l)
PY
  [ -e db_$R ] || mk $S/$R.aa.fa $S/$R.3di.fa db_$R
  mk $R.q.aa.fa $R.q.3di.fa q_$R
  t0=$(date +%s)
  nice -n 10 $FS search q_$R db_$R aln_$R tmp_$R --alignment-type 2 -e 1e-3 --max-seqs 5000 --exhaustive-search 1 --threads 12 -v 1
  nice -n 10 $FS convertalis q_$R db_$R aln_$R $R.m8.part --format-output query,target,evalue,bits,qstart,qend,qlen,tstart,tend,tlen --threads 12 -v 1
  mv $R.m8.part $R.m8; echo "$R exhaustive done $(( $(date +%s)-t0 )) s"
  rm -rf tmp_$R aln_$R* q_$R*
done
