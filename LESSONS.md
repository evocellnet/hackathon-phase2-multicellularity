# LESSONS — process rules for running this project efficiently

One line per rule, with the evidence that motivated it. Pruned at each milestone.

1. Start long, independent jobs (environment build, large downloads) before planning work, so they finish while thinking. (Env build and Pfam/taxdump downloads launched at t=0.)
2. Never `pip install` extra tools into the pinned Ancestra env: metapredict silently upgraded numpy/scipy and pulled torch (11 test failures risk). Use a separate env (`conda_envs/annot`) per tool family and re-run the Ancestra tests after any env change.
3. Background subshells launched with `( ... ) &` inside a Bash call can die silently (Pfam download produced nothing); use the tool's run_in_background mode for downloads.
