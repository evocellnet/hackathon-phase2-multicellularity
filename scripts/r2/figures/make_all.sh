#!/bin/bash
# Regenerate all report v2 figures into report/figures_v2/ (PNG 300 dpi + PDF).
set -euo pipefail
P=/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity
cd "$P/scripts/r2/figures"
for f in fig1_overview fig2_orthology_agreement fig3_novelty_events fig4_family_sets fig5_complex_families fig6_structure; do
  "$P/conda_envs/ancestra/bin/python" "$f.py"
done
