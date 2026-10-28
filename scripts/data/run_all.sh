#!/bin/bash
# Rebuild the proteome dataset end to end (re-runnable; finished steps are cached).
# Requires data/raw/{eukprot/proteins,lukprot/extracted,taxdump,busco/busco_downloads} (see DATA.md).
# CPU use is capped at 20 (preprocess 5 jobs x 4 CD-HIT threads; BUSCO 5 jobs x 4 CPUs).
set -euo pipefail
P=$(cd "$(dirname "$0")/../.." && pwd)
PY=$P/conda_envs/ancestra/bin/python
L=$P/logs
mkdir -p $L
cd $P
$PY scripts/data/select_final.py                          2>&1 | tee $L/data_select_final.log
$PY scripts/data/fetch_proteomes.py --threads 6            2>&1 | tee $L/data_fetch.log
$PY scripts/data/preprocess_proteomes.py --jobs 5 --cdhit-threads 4 2>&1 | tee $L/data_preprocess.log
$PY scripts/data/run_busco.py --jobs 5 --cpu 4             2>&1 | tee $L/data_busco.log
$PY scripts/data/build_guide_trees.py                     2>&1 | tee $L/data_guide_trees.log
$PY scripts/data/build_transitions.py                     2>&1 | tee $L/data_transitions.log
$PY scripts/data/finalize_tables.py                       2>&1 | tee $L/data_finalize.log
$PY scripts/data/write_data_md.py                         2>&1 | tee $L/data_write_data_md.log
$PY scripts/data/check_dataset.py                         2>&1 | tee $L/data_check.log
