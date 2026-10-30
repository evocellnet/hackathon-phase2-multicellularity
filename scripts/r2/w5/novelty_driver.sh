#!/usr/bin/env bash
# W5 task B driver: for each panel, once (a) the panel agent's ladder and tables exist and (b) ProstT5 3Di is
# complete for all panel species, run struct_novelty.py prep -> search -> classify (12 threads, nice 10).
# Polls every 5 min; stops at the given deadline (epoch seconds, default: 2026-10-30 14:00 CET) or when all done.
set -uo pipefail
P=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity
cd $P
PY=$P/conda_envs/ancestra/bin/python
DEADLINE=${1:-$(date -d "2026-10-30 14:00" +%s)}
PANELS="R2 R1 R4 R3 R6 R5 R7 R8"
while [ $(date +%s) -lt $DEADLINE ]; do
  left=0
  for R in $PANELS; do
    [ -e results/r2/w5/novelty/$R/DONE.json ] && continue
    [ -e results/r2/w5/novelty/$R/FAILED ] && continue
    left=$((left+1))
    L=results/r2/panels/$R
    [ -e $L/work/ladder/classified.parquet ] && [ -e $L/work/ladder/per_query.parquet ] && [ -e $L/branch_table.tsv ] && [ -e $L/event_table.tsv ] || continue
    miss=0; for f in data/dense/$R/proteomes/*; do s=$(basename $f); s=${s%.*}; [ "$R:$s" = "R4:Eudosp" ] && continue; [ -e data/structure3di/all/$s.3di.fa.done ] || miss=1; done  # Eudosp never queued for 3Di (struct_novelty.ALLOW_MISSING)
    [ $miss -eq 0 ] || continue
    echo "[$(date -Iseconds)] $R: start"
    ( [ -e results/r2/w5/novelty/$R/sample.tsv ] || $PY scripts/r2/w5/struct_novelty.py prep $R ) && \
    ( [ -e results/r2/w5/novelty/$R/hits.tsv.gz ] || $PY scripts/r2/w5/struct_novelty.py search $R ) && \
    $PY scripts/r2/w5/struct_novelty.py classify $R \
      || { echo "[$(date -Iseconds)] $R: FAILED"; mkdir -p results/r2/w5/novelty/$R; touch results/r2/w5/novelty/$R/FAILED; }
    echo "[$(date -Iseconds)] $R: end"
  done
  [ $left -eq 0 ] && { echo "all panels done"; break; }
  sleep 300
done
