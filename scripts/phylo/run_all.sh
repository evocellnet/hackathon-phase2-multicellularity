#!/bin/bash
# Steps 1-3 for euk then prok. CPU budget (PI, 13:5x): <= 2 single-thread jobs under nice -n 19 until
# runs/euk_fast/DONE.json exists, then 12. Re-evaluated before every step.
set -euo pipefail
P=$(cd "$(dirname "$0")/../.." && pwd); cd "$P"
export PATH=$P/conda_envs/ancestra/bin:$PATH
nt() { if [ -f runs/euk_fast/DONE.json ]; then echo 12; else echo 2; fi; }
log() { echo "[$(date +%FT%T)] $*"; }
for S in euk prok; do
  if [ ! -f runs/phylo/$S/build_markers.DONE.json ]; then
    log "$S step1 jobs=$(nt)"; nice -n 19 python scripts/phylo/build_markers.py --set $S --jobs $(nt); fi
  if [ ! -f runs/phylo/$S/species_tree_bl.DONE.json ]; then
    log "$S step2 threads=$(nt)"; nice -n 19 python scripts/phylo/species_tree_bl.py --set $S --threads $(nt); fi
  if [ ! -f runs/phylo/$S/gene_trees.DONE.json ]; then
    log "$S step3 jobs=$(nt)"; nice -n 19 python scripts/phylo/gene_trees.py --set $S --jobs $(nt); fi
  log "$S done"
done
nice -n 19 python scripts/phylo/summarize.py > /dev/null
log "ALL DONE"
