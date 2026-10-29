#!/usr/bin/env bash
# W5 task A: wait for ProstT5 3Di of all 37 R1 species, build the R1 3Di FASTA, then run Ancestra `fast`
# with --family-rescue 3di (runs/r2/R1_3di) and with --family-rescue 3di --family-rescue-guard
# (runs/r2/R1_3di_guard). The S0 (families) and S2 (scores) stages are copied from the default run
# runs/r2/R1: their stage keys depend only on input md5s and configs (not on paths or thread counts), so
# Ancestra skips them and runs only rescue + merge + output. The guard run also reuses the rescue stage.
# 12 threads, nice 10. Refuses to overwrite existing run directories.
set -uo pipefail
P=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity
cd $P
PY=$P/conda_envs/ancestra/bin/python
IN=runs/r2/R1_3di_input
SPECIES=$(ls data/dense/R1/proteomes | sed 's/\.fa$//')
while true; do
  n=0; for s in $SPECIES; do [ -e data/structure3di/all/$s.3di.fa.done ] || n=$((n+1)); done
  [ $n -eq 0 ] && break
  echo "[$(date -Iseconds)] waiting for 3Di: $n R1 species missing"; sleep 120
done
mkdir -p $IN
[ -e $IN/R1.3di.fa ] || $PY scripts/r2/w5/make_3di_fasta.py runs/r2/R1 $IN/R1.3di.fa || exit 1
md5sum $IN/R1.3di.fa > $IN/R1.3di.fa.md5
for V in 3di 3di_guard; do
  OUT=runs/r2/R1_$V
  if [ -e $OUT ]; then echo "exists: $OUT (skip)"; continue; fi
  mkdir -p $OUT/stages
  cp -a runs/r2/R1/stages/input runs/r2/R1/stages/families runs/r2/R1/stages/scores $OUT/stages/
  EXTRA="--family-rescue 3di --structure-3di $P/$IN/R1.3di.fa"
  if [ $V = 3di_guard ]; then
    [ -d runs/r2/R1_3di/stages/rescue ] && cp -a runs/r2/R1_3di/stages/rescue $OUT/stages/
    EXTRA="$EXTRA --family-rescue-guard"
  fi
  echo "[$(date -Iseconds)] start $V: $EXTRA"
  ANCESTRA_EXTRA="$EXTRA" nice -n 10 scripts/r2/w5/run_ancestra_patched_extra.sh data/dense/R1/proteomes \
    data/dense/R1/guide_tree.nwk $OUT 12 fast 150 --workers 8 --progress-log $OUT/fastleiden_progress.log
  echo "[$(date -Iseconds)] end $V rc=$?"
done
