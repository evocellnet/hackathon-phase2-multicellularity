#!/bin/bash
# Regenerate all report v2 figures into report/figures_v2/ (PNG 300 dpi + PDF).
# Figure 5 (cell-death architectures) is drawn by scripts/r2/w4_celldeath/04_figure.py.
set -euo pipefail
P=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity
PY="$P/conda_envs/ancestra/bin/python"
cd "$P/scripts/r2/figures"
for f in fig1_overview fig2_orthology_agreement fig3_novelty_events fig4_family_sets; do
  "$PY" "$f.py"
done
"$PY" "$P/scripts/r2/w4_celldeath/04_figure.py"
for f in fig6_complex_families fig7_structure; do
  "$PY" "$f.py"
done
