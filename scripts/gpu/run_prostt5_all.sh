#!/bin/bash
# Predict ProstT5 3Di for every proteome, one species at a time (resumable). Order file as $1 (default data/structure3di/order_r2b.txt).
# Proteome path: data/structure3di/panel_sources.tsv (dense-panel species) else data/proteomes/{euk,prok}/<code>.fa
P=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity
cd $P
export HF_HOME=$P/data/raw/prostt5/hf
ORDER=${1:-data/structure3di/order_r2b.txt}
for c in $(cat $ORDER); do
  out=data/structure3di/all/$c.3di.fa
  [ -f $out.done ] && continue
  in=$(awk -v c=$c '$1==c{print $2}' data/structure3di/panel_sources.tsv)
  [ -z "$in" ] && in=data/proteomes/euk/$c.fa
  [ -f $in ] || in=data/proteomes/prok/$c.fa
  conda_envs/gpu/bin/python scripts/gpu/prostt5_3di.py $out $in > logs/prostt5/$c.log 2>&1 && touch $out.done
done
echo ALL_DONE
