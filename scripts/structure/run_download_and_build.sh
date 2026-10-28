#!/bin/bash
# Resumable: download AFDB reference models (rate-limited ~3 req/s => ~12-13 h for ~124k), then build Foldseek DBs + report.
# Usage: nohup scripts/structure/run_download_and_build.sh > logs/structure_download.log 2>&1 &
export P=${P:-/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity}; cd $P
PY=$P/conda_envs/ancestra/bin/python
nice -n 19 $PY scripts/structure/download_ref_afdb.py "$@" && \
nice -n 19 $PY scripts/structure/build_refdb.py --threads 4 && \
nice -n 19 $PY scripts/structure/report_ref_download.py
