#!/usr/bin/env python
"""H9 stability across the pre-registered seeds [11, 23, 37, 41, 53] (DESIGN §7; round-1 S6).

Two sources of randomness: (a) the permutation draws (Monte-Carlo stability) and (b) the random
± orientation of the control r (branchstats.h9_values). Seed 11 is the pre-registered primary for
both; other seeds are stability checks only. Writes results/h9/h9_seed_stability.tsv and appends
one registry row per non-primary combination (tier exploratory, label=seed_stability).
"""
import os
import sys
import warnings

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "phylo"))
from mcell import branchstats as bs  # noqa: E402
import h9_rates as H  # noqa: E402

SEEDS = [11, 23, 37, 41, 53]


def main():
    warnings.simplefilter("ignore")
    H.verify_prereg()
    tr, traits, tree, bt, gts = H.setup("euk")
    stems = bs.stem_branches(bt)
    pool = bs.h9_control_branches(tree, bt)
    m9 = bs.match_controls(bt, stems, pool=pool)
    per = pd.read_csv(os.path.join(H.OUT, "h9_transitions.tsv"), sep="\t").set_index("transition_id")
    ctrl = pd.read_csv(os.path.join(H.OUT, "h9_controls.tsv"), sep="\t").set_index("control")
    r_st = {t: per.loc[t, "r"] for t in stems.index}
    rows = []
    for kind in ("draw_seed", "orientation_seed"):
        for s in SEEDS:
            o, d = (11, s) if kind == "draw_seed" else (s, 11)
            v = bs.h9_values(r_st, dict(stems), ctrl["r"].to_dict(), orient=True, seed=o)
            res = bs.relabel_test(v, m9, alternative="two-sided", n_draws=H.N_DRAWS, seed=d).summary()
            rows.append(dict(varied=kind, orientation_seed=o, draw_seed=d, T_obs=res["T_obs"], p=res["p"],
                             p_mcse=res["p_mcse"], n_controls_flipped_vs_seed11=None))
            if s != 11:
                bs.append_registry(H.REG, dict(
                    hypothesis="H9", statistic="median per-gene log rate ratio; relabelling T (seed stability)",
                    n=res["n"], k=res["k"], estimate="NA", ci_low="NA", ci_high="NA", p=float(res["p"]),
                    p_mcse=float(res["p_mcse"]), holm_p="NA", tier="exploratory", seed=d, n_draws=H.N_DRAWS,
                    scheme="permute", notes=f"label=seed_stability; {kind}={s}; orientation_seed={o}; "
                                            f"draw_seed={d}; T_obs={res['T_obs']:.4f}"))
    out = pd.DataFrame(rows).drop(columns="n_controls_flipped_vs_seed11")
    out.to_csv(os.path.join(H.OUT, "h9_seed_stability.tsv"), sep="\t", index=False, float_format="%.5g")
    print(out.to_string())


if __name__ == "__main__":
    main()
