# LESSONS — process rules for running this project efficiently

One line per rule, with the evidence that motivated it. Pruned at milestone M3 (t≈8 h) from 9 raw rules after a fresh-agent retrospective (agents/retrospective_M3.md).

1. Start long independent jobs (env builds, downloads) at t=0, before planning. (Env/Pfam/taxdump launched at t=0; dataset build before review end, D1.)
2. Run a cheap pilot before any long run on new data, and put the result in the tree: component-size/timing/rate-limit. (No pilot for eukaryote Ancestra: giant 73% component cost t=4-6 h, D8; AFDB 3 req/s turned 1-1.5 h into 12-13 h, forcing a Swiss-Prot fallback.)
3. Run every pipeline end-to-end on the small/cheap dataset with positive controls before the expensive one. (Prok dry run recovered P1/P2; patched Leiden verified byte-identical on prokaryotes.)
4. Calibrate each statistic on simulations/the real tree before using it, and fix interpretation thresholds in advance. (Relabelling null 14.5% FPR; split repair failed; Dev-4 2.2x rule; raw Ancestra gains 2.1-2.6x inflated.)
5. Get independent adversarial review before freezing and again before relying on results; the author does not review their own code. (Three design rounds rewrote stems/nulls/statistic; shadowed T_volvocales_gv bug found by another agent.)
6. Log every post-freeze change as a dated deviation (run `date` in a SEPARATE call before writing the entry — same-command date+edit produced three wrong stamps) stating what results had been seen; use `date` against START.txt, never edit old timestamps, and add a corrective entry instead. (Dev-2 and Dev-4 timestamp-fix commits; M1 mis-stamped by 1 h.)
7. Keep the research tree current at every milestone (summary, status boxes, order), separate primary from exploratory results, and label underpowered nulls as such. (Tree last updated M2 at 2.1 h while at 7.5 h; null results marked [-] despite 13-35% power.)
8. Cap CPUs per agent (explicit thread count, `nice -n 19` for non-critical jobs). (Load 73 on 24 CPUs stalled the critical Ancestra run at t~2 h.)
9. Isolate environments and always call env binaries by full path or env-first PATH; never install into the pinned env. (metapredict upgraded numpy/scipy, pulled torch; BUSCO 'No module named busco' cost two rebuilds.)
10. Use the tool's background mode for long jobs and kill by PID, never `pkill -f` with a pattern in your own command line. (`( ... ) &` Pfam download died silently; pkill -f killed the caller twice, exit 144.)
11. Analysis scripts that rewrite a whole results file must be run with all their tests together (or merge, not overwrite): running only H1/H2 overwrote the H4/H5/H12 rows (t≈12 h); a full deterministic rerun fixed it and doubled as a determinism check.
12. Expect subagents to be blocked from writing report `.md` files in this harness; have them write `.txt` and convert, or extract their returned text programmatically (saved ~4,000 words of re-typing at t≈7 h).
13. Pilot GPU tools on 1-5k proteins and check the device is really used (`Use GPU 0` in Foldseek's log; nvidia-smi at 0%) before planning on them; Foldseek's GPU ProstT5 silently ran on CPU, HF autoregressive ProstT5 was 10x slower than the encoder+CNN path (R2, 30 min).
