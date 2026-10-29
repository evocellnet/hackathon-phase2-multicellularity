#!/usr/bin/env bash
# Re-run every round-1 errata computation (writes only to results/round1_errata/; ~3 min, 12 threads).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="$ROOT/conda_envs/ancestra/bin/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$ROOT/logs/errata" "$ROOT/results/round1_errata"
nice -n 10 "$PY" "$ROOT/scripts/errata/smd_provenance.py"            > "$ROOT/logs/errata/smd_provenance.log" 2>&1
nice -n 10 "$PY" "$ROOT/scripts/errata/power_table_real_bl.py" --bl placeholder --jobs 12 > "$ROOT/logs/errata/power_placeholder.log" 2>&1
nice -n 10 "$PY" "$ROOT/scripts/errata/power_table_real_bl.py" --bl real --jobs 12        > "$ROOT/logs/errata/power_real.log" 2>&1
nice -n 10 "$PY" "$ROOT/scripts/errata/compare_power.py"             > "$ROOT/logs/errata/compare_power.log" 2>&1
nice -n 10 "$PY" "$ROOT/scripts/errata/positive_controls_classify.py" > "$ROOT/logs/errata/positive_controls_classify.log" 2>&1
OMP_NUM_THREADS=4 nice -n 10 "$PY" "$ROOT/scripts/errata/h9_all_orientations.py" > "$ROOT/logs/errata/h9_all_orientations.log" 2>&1
echo done
