#!/bin/bash
# Add-on tip ladders (DEVIATIONS D-P6): start after the main ladder queue has no work, 12 threads, one at a time.
P=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity; cd $P
PY=$P/conda_envs/ancestra/bin/python; L=$P/logs/r2panels; O=$P/results/r2/panels
main_done() { for X in R1 R2 R3 R4 R5 R6 R7 R8 R4_curated; do [ -e $O/$X/branch_gain_ladder.tsv ] || [ -e $L/${X}_ladder.failed ] || return 1; done; return 0; }
for R in R2 R3 R5 R8 R1 R7; do
  while pgrep -f "panel.py ladder R" >/dev/null || ! main_done; do sleep 60; done
  [ -e $O/$R/tip_gain_ladder.tsv ] && continue
  echo "[$(date '+%F %T')] $R tip ladder start" >> $L/driver.log
  nice -n 10 $PY scripts/r2/panels/panel.py ladder $R --tips-only --threads 12 >> $L/${R}_ladder_tips.log 2>&1
  echo "[$(date '+%F %T')] $R tip ladder rc=$?" >> $L/driver.log
done
