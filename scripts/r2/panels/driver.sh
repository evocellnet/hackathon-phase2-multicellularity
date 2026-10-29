#!/bin/bash
# Two queues, <= 16 CPUs total, nice 10, panels in order of Ancestra completion (runs/r2/<R>/DONE.json):
#   ladder queue (12 threads): DIAMOND novelty ladder, one panel at a time;
#   side queue (4 threads): parse, phylo, prereg Pfam (pre) and tables, controls, W3, PROGRESS (post).
# Every step is skipped when its output exists (the scripts refuse to overwrite).
P=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity
cd $P
PY=$P/conda_envs/ancestra/bin/python
S=$P/scripts/r2/panels
L=$P/logs/r2panels
O=$P/results/r2/panels
mkdir -p $L
log() { echo "[$(date '+%F %T')] $*" >> $L/driver.log; }
step() { local R=$1 n=$2 out=$3; shift 3
  while pgrep -f "panel.py $n $R( |$)" >/dev/null; do sleep 10; done
  [ -e "$out" ] && return 0
  [ -e "$L/${R}_$n.failed" ] && return 1
  log "$R $n start"; nice -n 10 "$@" >> $L/${R}_$n.log 2>&1; local rc=$?
  log "$R $n rc=$rc"; [ $rc -ne 0 ] && touch $L/${R}_$n.failed; return $rc; }
arrived() { ls -tr $P/runs/r2/*/DONE.json 2>/dev/null | sed 's#.*/r2/\([^/]*\)/DONE.json#\1#' | grep -xE 'R[1-8]|R4_curated'; }
ladder_queue() {
  while true; do
    while pgrep -f "panel.py ladder R" >/dev/null; do sleep 30; done
    local todo=0
    for R in R1 R2 R3 R4 R5 R6 R7 R8 R4_curated; do [ -e $O/$R/branch_gain_ladder.tsv ] || [ -e $L/${R}_ladder.failed ] || todo=1; done
    [ $todo -eq 0 ] && break
    for R in $(arrived); do
      [ -e $O/$R/branch_events.tsv ] || continue
      if [ ! -e $O/$R/branch_gain_ladder.tsv ] && [ ! -e $L/${R}_ladder.failed ] && ! pgrep -f "panel.py ladder $R" >/dev/null; then
        step $R ladder $O/$R/branch_gain_ladder.tsv $PY $S/panel.py ladder $R --threads 12; continue 2; fi
    done
    sleep 60
  done; log "ladder queue empty"; }
side_queue() {
  while true; do
    local all=1 did=0
    for R in $(arrived); do
      if [ -e $O/$R/branch_gain_ladder.tsv ] && [ -e $O/$R/tree_bl.nwk ] && [ ! -e $O/$R/w3_families.tsv.gz ] && [ ! -e $L/${R}_w3.failed ]; then
        step $R tables $O/$R/event_table.tsv $PY $S/panel.py tables $R && \
        step $R controls $O/$R/positive_controls.tsv $PY $S/panel.py controls $R
        step $R w3 $O/$R/w3_families.tsv.gz $PY $S/w3_panel.py $R --cpus 4
        $PY $S/progress.py > /dev/null 2>> $L/progress.log; did=1; break
      fi
      if [ ! -e $O/$R/work/pfam_prereg/.done ] && [ ! -e $L/${R}_prereg_pfam.failed ]; then
        step $R parse $O/$R/branch_events.tsv $PY $S/panel.py parse $R
        step $R phylo $O/$R/tree_bl.nwk $PY $S/panel.py phylo $R --threads 4
        step $R prereg_pfam $O/$R/work/pfam_prereg/.done bash -c "$PY $S/panel.py prereg_pfam $R --cpus 4 && touch $O/$R/work/pfam_prereg/.done"
        did=1; break
      fi
    done
    for R in R1 R2 R3 R4 R5 R6 R7 R8 R4_curated; do [ -e $O/$R/w3_families.tsv.gz ] || [ -e $L/${R}_w3.failed ] || [ -e $L/${R}_tables.failed ] || all=0; done
    [ $all -eq 1 ] && break
    [ $did -eq 0 ] && sleep 60
  done; log "side queue empty"; }
ladder_queue & side_queue & wait
log "driver finished"
