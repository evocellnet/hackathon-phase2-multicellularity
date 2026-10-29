# LESSONS — process rules for running this project efficiently

One line per rule, with evidence. Pruned at R2 (2026-10-29 ~22:15) from a fresh-agent retrospective (agents/retrospective_R2.txt); previous version agents/LESSONS_before_R2_prune.md.

1. Start long independent jobs (env builds, downloads, GPU runs) at t=0, before planning. (Env/Pfam/taxdump at t=0; ProstT5 3Di launched at R2 t~0.5 h and finished while panels ran.)
2. Pilot before any long run on new data or new hardware (size/timing/rate limit; check the device is really used): giant 73% component cost t=4-6 h (D8); AFDB 3 req/s gave 12-13 h not 1.5 h; Foldseek GPU silently ran on CPU, 15 vs 48 prot/s (R2 pilot).
3. Run each pipeline end-to-end on a small dataset with positive controls before the expensive one, and smoke-test every analysis script/env on a fixture. (Prok dry run recovered P1/P2; glmmTMB ABI failure forced D-W4-1; filter bug D-P11.)
4. Calibrate each statistic on simulations or control data before freezing it, and fix thresholds in advance, including covariate range and tip/internal comparability. (Relabelling null 14.5% FPR; Dev-4 2.2x rule; tip events scored 21x vs internal controls, D-P1/D-P6.)
5. Compute power before pre-registering; if power <~40%, label the test as estimation, not confirmation. (R2-P1 3 events, 18-41% power, ended p=0.072; round 1 power 13-35%.)
6. Get independent adversarial review before freezing and before relying on results; authors do not verify their own work. (3 design rounds; verification b5f0e7c found P3 size caveat; round-1 shadowed-variable bug.)
7. Log every post-freeze change as a dated deviation with a "Seen:" line stating what results were visible; run `date` in a separate call; correct by adding entries, never editing. (Dev-2/Dev-4 timestamp fixes; panels/DEVIATIONS.txt D-P1..D-P11.)
8. Never auto-launch analyses from wait loops; confirmatory outputs are computed once with a declared tag, everything later is labelled sensitivity. (W3 premature run D-P7 19:03; D-P9 'final_8panels'.)
9. Keep the research tree current at each result commit, with primary vs exploratory separated, nulls labelled with power, and decision rules (go/no-go times) recorded when triggered. (Tree lagged M2 by 5 h in round 1; W5/W6 status stale at R2 t=7 h.)
10. Commit scripts and logs with the run that uses them; keep output dirs guarded (refuse-if-exists, DONE.json md5). (20 untracked logs/scripts in git status after 4395cd0.)
11. Cap CPUs per agent (explicit thread count, nice -n 19 for non-critical jobs) and isolate envs; call binaries by full path, never install into the pinned env. (Load 73 on 24 CPUs stalled Ancestra at t~2 h; metapredict upgraded numpy; BUSCO two rebuilds.)
12. Use the tool's background mode for long jobs; kill by explicit PID from a separate call, never pkill -f or a pattern-kill inside a command containing the pattern. (Exit 144 twice; R2 17:27 killed the caller and orphaned a GPU job.)
13. Scripts that rewrite a whole results file must run all their tests together or merge. (H1/H2-only run overwrote H4/H5/H12 rows at t~12 h.)
14. Cross-check key results with a second method or dataset and report disagreement plainly; subagents write .txt not .md. (EukaryOMA agreed on only 29% of stem gains; round-1 errata SMD 1.77->2.06; ~4,000 words saved by extracting subagent text.)
