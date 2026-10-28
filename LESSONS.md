# LESSONS — process rules for running this project efficiently

One line per rule, with the evidence that motivated it. Pruned at each milestone.

1. Start long, independent jobs (environment build, large downloads) before planning work, so they finish while thinking. (Env build and Pfam/taxdump downloads launched at t=0.)
2. Never `pip install` extra tools into the pinned Ancestra env: metapredict silently upgraded numpy/scipy and pulled torch (11 test failures risk). Use a separate env (`conda_envs/annot`) per tool family and re-run the Ancestra tests after any env change.
3. Background subshells launched with `( ... ) &` inside a Bash call can die silently (Pfam download produced nothing); use the tool's run_in_background mode for downloads.
4. Conda env CLI scripts with '#!/usr/bin/env python3' shebangs run the base python unless the env's bin is first on PATH; always prepend the env bin to PATH. (BUSCO 'No module named busco' cost two rebuilds.)
5. Check `date` against START.txt before writing any timestamp; never estimate elapsed time from agent durations (M1 was mis-stamped by 1 h).
6. Give every agent an explicit CPU cap and require `nice -n 19` for non-critical work; four agents each launching 'a few' jobs pushed load to 73 on 24 CPUs and stalled the critical-path Ancestra run (t≈2h).
