#!/usr/bin/env bash
# Run Ancestra (frozen freeze-v1) on one proteome set with a memory watchdog.
# Usage: scripts/run_ancestra.sh <proteome_dir> <guide_tree.nwk> <out_dir> <threads> <preset> [max_rss_gb]
# Writes <out_dir>/DONE.json only if Ancestra exits 0 and the OrthoXML verifier passes.
set -uo pipefail
P=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity
export PATH=$P/conda_envs/ancestra/bin:$PATH
PROT=$1; TREE=$2; OUT=$3; THREADS=$4; PRESET=$5; MAXGB=${6:-145}
mkdir -p "$OUT"
LOG="$OUT/run.log"
echo "[$(date -Iseconds)] start ancestra preset=$PRESET threads=$THREADS prot=$PROT tree=$TREE" | tee -a "$LOG"
ancestra run --proteomes "$PROT" --tree "$TREE" --out "$OUT" --threads "$THREADS" --preset "$PRESET" >> "$LOG" 2>&1 &
PID=$!
PEAK=0
while kill -0 $PID 2>/dev/null; do
  # total RSS (kB) of the process tree rooted at PID
  RSS=$(ps -o rss= --ppid $PID -p $PID 2>/dev/null | awk '{s+=$1} END {print s+0}')
  for c in $(pgrep -P $PID); do RSS=$((RSS + $(ps -o rss= --ppid $c 2>/dev/null | awk '{s+=$1} END {print s+0}'))); done
  [ "$RSS" -gt "$PEAK" ] && PEAK=$RSS
  if [ "$RSS" -gt $((MAXGB * 1024 * 1024)) ]; then
    echo "[$(date -Iseconds)] WATCHDOG: RSS ${RSS} kB > ${MAXGB} GB, killing" | tee -a "$LOG"
    pkill -TERM -P $PID; kill -TERM $PID; sleep 10; pkill -KILL -P $PID; kill -KILL $PID
    echo "{\"status\":\"killed_memory\",\"peak_rss_kb\":$PEAK}" > "$OUT/FAILED.json"; exit 2
  fi
  echo "$(date +%s) $RSS" >> "$OUT/rss_trace.txt"
  sleep 30
done
wait $PID; RC=$?
echo "[$(date -Iseconds)] ancestra exit $RC peak_rss_kb=$PEAK" | tee -a "$LOG"
if [ $RC -ne 0 ]; then echo "{\"status\":\"failed\",\"rc\":$RC}" > "$OUT/FAILED.json"; exit $RC; fi
python $P/external/ancestra-mvp-ai/scripts/verify_orthoxml.py "$OUT" >> "$LOG" 2>&1; VRC=$?
echo "{\"status\":\"done\",\"preset\":\"$PRESET\",\"threads\":$THREADS,\"peak_rss_kb\":$PEAK,\"verify_rc\":$VRC,\"finished\":\"$(date -Iseconds)\"}" > "$OUT/DONE.json"
cat "$OUT/DONE.json"
