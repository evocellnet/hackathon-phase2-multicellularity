#!/bin/bash
# Predict ProstT5 3Di for every eukaryotic proteome, one species at a time (resumable; order in data/structure3di/euk_order.txt).
P=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity
cd $P
export HF_HOME=$P/data/raw/prostt5/hf
for c in $(cat data/structure3di/euk_order.txt); do
  in=data/proteomes/euk/$c.fa; [ -f $in ] || in=data/proteomes/prok/$c.fa
  out=data/structure3di/all/$c.3di.fa
  [ -f $out.done ] && continue
  conda_envs/gpu/bin/python scripts/gpu/prostt5_3di.py $out $in > logs/prostt5/$c.log 2>&1 && touch $out.done
done
echo ALL_DONE
