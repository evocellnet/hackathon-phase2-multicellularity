#!/bin/bash
set -euo pipefail
FS=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/conda_envs/ancestra/bin/foldseek
W=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/results/r2/verification/w5/work/spot
cd $W
mk() { # $1 aa fasta, $2 3di fasta, $3 db name
  nice -n 10 $FS base:createdb $1 $3 --shuffle 0 -v 1
  nice -n 10 $FS base:createdb $2 $3_tmp3di --shuffle 0 -v 1
  mv $3_tmp3di $3_ss; mv $3_tmp3di.index $3_ss.index; mv $3_tmp3di.dbtype $3_ss.dbtype
  rm -f $3_tmp3di*
}
for R in R1 R2 R3 R4 R5 R6 R7 R8; do
  [ -s $R.m8 ] && continue
  python3 - <<PY
ids=set(open('$R.qids').read().split())
for ext in ('aa','3di'):
    out=open('$R.q.'+ext+'.fa','w'); keep=False
    for l in open('$R.'+ext+'.fa'):
        if l[0]=='>': keep=l[1:].strip() in ids
        if keep: out.write(l)
PY
  mk $R.aa.fa $R.3di.fa db_$R
  mk $R.q.aa.fa $R.q.3di.fa q_$R
  # sanity: same ordering of headers in aa and 3di dbs
  cmp <(cut -f2 db_$R.lookup) <(grep '>' $R.3di.fa | sed 's/>//') && echo "$R lookup order OK"
  t0=$(date +%s)
  nice -n 10 $FS search q_$R db_$R aln_$R tmp_$R --alignment-type 2 -e 1e-3 --max-seqs 10000 --threads 12 -v 1
  nice -n 10 $FS convertalis q_$R db_$R aln_$R $R.m8.part --format-output query,target,evalue,bits,qstart,qend,qlen,tstart,tend,tlen --threads 12 -v 1
  mv $R.m8.part $R.m8
  echo "$R search done $(( $(date +%s)-t0 )) s"
  rm -rf tmp_$R aln_$R* db_$R* q_$R* 
done
