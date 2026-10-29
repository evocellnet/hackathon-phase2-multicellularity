#!/usr/bin/env python
"""Errata item 2: old (frozen, placeholder bl) vs new (real bl) power at identical effects.
Also checks that the placeholder re-run reproduces the frozen results/prereg tables exactly.
Output: results/round1_errata/power_old_vs_new.tsv, power_mde_old_vs_new.tsv"""
import os
import pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
R = lambda *p: os.path.join(ROOT, *p)
old = pd.read_csv(R("results/prereg/power_curves.tsv"), sep="\t")
rep = pd.read_csv(R("results/round1_errata/power_placeholder/power_curves.tsv"), sep="\t")
new = pd.read_csv(R("results/round1_errata/power_real/power_curves.tsv"), sep="\t")
key = ["test", "param_value", "effect"]
for d in (old, rep, new):
    d["param_value"] = d["param_value"].astype(str)
chk = old.merge(rep, on=key, suffixes=("", "_rep"))
ok = ((chk.power_0_05 == chk.power_0_05_rep) & (chk.power_0_01 == chk.power_0_01_rep)).all()
print(f"placeholder re-run reproduces frozen power_curves.tsv exactly: {ok} ({len(chk)} rows)")
m = old.merge(new, on=key, suffixes=("_old_placeholder", "_new_realbl"))
m = m[key + ["power_0_05_old_placeholder", "power_0_05_new_realbl", "power_0_01_old_placeholder",
             "power_0_01_new_realbl", "n_sim_old_placeholder"]].rename(columns={"n_sim_old_placeholder": "n_sim"})
m["diff_0_05"] = m.power_0_05_new_realbl - m.power_0_05_old_placeholder
m.to_csv(R("results/round1_errata/power_old_vs_new.tsv"), sep="\t", index=False, float_format="%.4g")
ot = pd.read_csv(R("results/prereg/power_table.tsv"), sep="\t")
nt = pd.read_csv(R("results/round1_errata/power_real/power_table.tsv"), sep="\t")
cols = ["mde_alpha_0_05", "mde_alpha_0_01", "power_at_1_5x_alpha_0_05", "size_at_effect_1_alpha_0_05"]
t = ot[["test", "variance_value"] + cols].merge(nt[["test", "variance_value"] + cols], on=["test", "variance_value"],
                                                suffixes=("_old", "_new"))
t.to_csv(R("results/round1_errata/power_mde_old_vs_new.tsv"), sep="\t", index=False, float_format="%.4g")
pd.set_option("display.width", 200)
print(t.to_string()); print(m.to_string())
