#!/usr/bin/env python
"""EXPLORATORY within-protein constraint proxy: radical/conservative replacement ratio (Kr/Kc) per transition.

For each of the 250 BUSCO markers (trimmed alignment with the taxa of its gene tree,
runs/phylo/euk/genetrees/<m>.fa) Fitch-parsimony changes are counted on every edge of the pruned guide
tree (mcell.krkc), classified radical/conservative by charge and by polarity (Zhang 2000). Per unit:
LR = log Kr/Kc(M side) − log Kr/Kc(S side), pooled over genes, 2,000 gene bootstraps (seed 11);
scopes: 'clade' (stem + crown vs sister stem + crown) and 'stem' (stem edges only).
Units: the 7 primary stems, the 3 secondary nested stems, and the 10 H9 control nodes (oriented with
the same random signs as the H9 test, seed 11) as the reference distribution; an exploratory
relabelling test (permute, two-sided, 10,000, seed 11) on LR_charge (clade scope).
Outputs results/krkc/*.tsv and summary.md; registry rows tier exploratory.
"""
from __future__ import annotations

import datetime as dt
import json
import math
import os
import sys
import warnings

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "phylo"))
from mcell import krkc, rates, branchstats as bs  # noqa: E402
from scripts.phylo.build_markers import read_fasta  # noqa: E402
import h9_rates as H  # noqa: E402

OUT = os.path.join(ROOT, "results", "krkc")
SEED = 11


def f(x, d=3):
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{d}f}"


def main():
    warnings.simplefilter("ignore")
    os.makedirs(OUT, exist_ok=True)
    tr, traits, tree, bt, gts = H.setup("euk")
    idx = rates._guide_index(tree)
    gdir = os.path.join(ROOT, "runs", "phylo", "euk", "genetrees")
    counts = {}
    for g, gt in gts.items():
        seqs = read_fasta(os.path.join(gdir, f"{g}.fa"))
        assert set(seqs) == set(gt.tips()), g
        counts[g] = krkc.edge_counts(gt, seqs, seed=SEED)
    stems = bs.stem_branches(bt)
    pool = bs.h9_control_branches(tree, bt)
    m9 = bs.match_controls(bt, stems, pool=pool)
    prim = rates.stem_triples(tr, tree, sets=["primary"], kingdom="euk")
    sec = rates.stem_triples(tr, tree, sets=["secondary"], kingdom="euk")
    ctrl = rates.control_triples(tree, pool)
    hc = pd.read_csv(os.path.join(ROOT, "results", "h9", "h9_controls.tsv"), sep="\t").set_index("control")
    h9 = pd.read_csv(os.path.join(ROOT, "results", "h9", "h9_transitions.tsv"), sep="\t").set_index("transition_id")
    sign = hc["orientation_sign"].to_dict()
    allrows = []
    for scope in ("clade", "stem"):
        tab = krkc.unit_table(counts, gts, idx, prim + sec + ctrl, scope=scope)
        tab.to_csv(os.path.join(OUT, f"krkc_gene_table_{scope}.tsv"), sep="\t", index=False)
        sm = krkc.summarise(tab, n_boot=2000, seed=SEED)
        sm["scope"] = scope
        sm["set"] = ["primary" if u in dict(stems) else ("secondary_nested" if u.startswith("T_") else "h9_control")
                     for u in sm.index]
        for s in krkc.SCHEMES:
            o = sm.index.map(lambda u: sign.get(u, 1.0))
            sm[f"LR_{s}_oriented"] = sm[f"LR_{s}"] * o
        allrows.append(sm)
    res = pd.concat(allrows)
    res["h9_r"] = [h9["r"].get(u, hc["r_oriented"].get(u, np.nan)) for u in res.index]
    res.to_csv(os.path.join(OUT, "krkc_units.tsv"), sep="\t", float_format="%.5g")
    # exploratory relabelling test on the clade-scope charge LR (same machinery and controls as H9)
    cl = res[res.scope == "clade"]
    vals = pd.Series({**{stems[t]: cl.loc[t, "LR_charge"] for t in stems.index if t in cl.index},
                      **{c: cl.loc[c, "LR_charge_oriented"] for c in pool if c in cl.index}})
    tests = {}
    for s in krkc.SCHEMES:
        v = pd.Series({**{stems[t]: cl.loc[t, f"LR_{s}"] for t in stems.index if t in cl.index},
                       **{c: cl.loc[c, f"LR_{s}_oriented"] for c in pool if c in cl.index}})
        r = bs.relabel_test(v, m9, alternative="two-sided", n_draws=10_000, seed=SEED).summary()
        tests[s] = r
        bs.append_registry(H.REG, dict(
            hypothesis="Q5/H9 within-protein constraint (Kr/Kc proxy)",
            statistic=f"log Kr/Kc ratio M vs S side ({s} classes, Fitch, stem+crown); relabelling T",
            n=r["n"], k=r["k"], estimate=float(np.median([cl.loc[t, f'LR_{s}'] for t in stems.index])),
            ci_low="NA", ci_high="NA", p=float(r["p"]), p_mcse=float(r["p_mcse"]), holm_p="NA",
            tier="exploratory", seed=SEED, n_draws=10_000, scheme="permute",
            input_md5s=[os.path.join(ROOT, "runs/phylo/euk/genetrees/branch_lengths.tsv")],
            notes=f"label=exploratory; estimate = median LR over primary stems; T_obs={r['T_obs']:.4f}; "
                  f"controls = H9 pool with H9 orientation signs"))
    json.dump(dict(finished=dt.datetime.now().isoformat(timespec="seconds"), n_genes=len(counts),
                   tests={k: {kk: (float(vv) if isinstance(vv, (int, float, np.floating)) else vv)
                              for kk, vv in v.items()} for k, v in tests.items()}),
              open(os.path.join(OUT, "krkc_run.json"), "w"), indent=1)
    render(res, tests, len(counts))
    print(res[res.scope == "clade"][["set", "n_genes", "LR_charge", "ci_lo_charge", "ci_hi_charge", "LR_charge_oriented",
                                      "LR_polarity", "ci_lo_polarity", "ci_hi_polarity", "fracR_M_charge",
                                      "fracR_S_charge", "h9_r"]].round(3).to_string())
    print(res[res.scope == "stem"][["set", "n_genes", "LR_charge", "ci_lo_charge", "ci_hi_charge",
                                     "n_changes_M_charge", "n_changes_S_charge"]].round(3).to_string())
    print(tests)


def render(res, tests, n_genes):
    L = []
    a = L.append
    cl, st = res[res.scope == "clade"], res[res.scope == "stem"]
    a("# Kr/Kc — radical vs conservative replacements (within-protein constraint proxy) — EXPLORATORY\n")
    a("Generated by `scripts/phylo/krkc_run.py` (module `mcell/krkc.py`, tests `tests/test_krkc.py`).\n")
    a("**Changes.** Fitch-parsimony ancestral states on the pruned guide tree of each of the "
      f"{n_genes} BUSCO markers (trimmed alignments; taxa with < 50 % coverage removed; gaps = missing). Polytomies use the generalised Fitch "
      "rule, and ties are resolved by a seeded random choice.\n")
    a("**Classes.** Radical/conservative replacement classes follow Zhang (2000): charge {RHK}/{DE}/neutral, and polarity.\n")
    a("**Statistic.** LR = log Kr/Kc(M side) − log Kr/Kc(S side), pooled over genes, with a 95 % CI from 2,000 gene bootstraps.")
    a("- **Scope.** The M side is the stem plus the crown; the S side is the sister stem plus its crown. LR > 0 means proportionally more radical replacements on the "
      "multicellular side, which is what relaxed constraint predicts.")
    a("- **Reference.** The reference distribution comes from the 10 H9 control nodes (same random orientation as H9).")
    a("- **Validation.** The tests check hand-countable toy alignments, a polytomy, and a planted simulation (LR ≈ 1.0 for 4× up-weighted radical "
      "targets on the M side, CI > 0; null |LR| ≤ 0.02, CI covers 0, across 6 seeds).\n")
    a("## Clade scope (stem + crown vs sister stem + crown)\n")
    a("| unit | set | genes | H9 r | LR charge [95 % CI] | LR polarity [95 % CI] | radical fraction M / S (charge) |")
    a("|---|---|---|---|---|---|---|")
    for u, r in cl.sort_values(["set", "h9_r"], ascending=[False, False]).iterrows():
        lc = r.LR_charge_oriented if r.set == "h9_control" else r.LR_charge
        lp = r.LR_polarity_oriented if r.set == "h9_control" else r.LR_polarity
        sg = -1 if (r.set == "h9_control" and r.LR_charge_oriented != r.LR_charge) else 1
        lo_c, hi_c = sorted([sg * r.ci_lo_charge, sg * r.ci_hi_charge])
        lo_p, hi_p = sorted([sg * r.ci_lo_polarity, sg * r.ci_hi_polarity])
        a(f"| {u} | {r.set} | {int(r.n_genes)} | {f(r.h9_r)} | {f(lc)} [{f(lo_c)}, {f(hi_c)}] | {f(lp)} [{f(lo_p)}, {f(hi_p)}] | "
          f"{f(r.fracR_M_charge)} / {f(r.fracR_S_charge)} |")
    a("\n(The H9 control rows are oriented as in the H9 test; for them the H9 r column shows the oriented control r.)\n")
    a("**Exploratory relabelling test** (same machinery and controls as H9; permute, two-sided, 10,000 draws, seed 11):")
    for s, r in tests.items():
        a(f"- {s}: n = {r['n']}, k = {r['k']}, T_obs = {f(r['T_obs'])}, p = {f(r['p'])} (MC SE {f(r['p_mcse'], 4)})")
    a("\n## Stem scope (stem edges only)\n")
    a("| unit | set | genes | LR charge [95 % CI] | changes on M / S stem edges |")
    a("|---|---|---|---|---|")
    for u, r in st[st.set != "h9_control"].iterrows():
        a(f"| {u} | {r.set} | {int(r.n_genes)} | {f(r.LR_charge)} [{f(r.ci_lo_charge)}, {f(r.ci_hi_charge)}] | "
          f"{int(r.n_changes_M_charge)} / {int(r.n_changes_S_charge)} |")
    prim = cl[cl.set == "primary"]
    pos = prim[prim.h9_r > 0]
    sup = pos[(pos.ci_lo_charge > 0) & (pos.ci_lo_polarity > 0)]
    a("\n## Reading for the H9 wording rule\n")
    a(f"Of the {len(pos)} primary transitions with H9 r > 0 ({', '.join(pos.index)}), {len(sup)} show LR > 0 with both CIs above 0"
      + (f" ({', '.join(sup.index)})" if len(sup) else "") + ". The across-transition relabelling tests are "
      + ("not significant" if all(t["p"] >= 0.05 for t in tests.values()) else "significant for at least one class")
      + f". Oriented control LR ranges: charge {f(cl[cl.set == 'h9_control'].LR_charge_oriented.min())} to "
        f"{f(cl[cl.set == 'h9_control'].LR_charge_oriented.max())}, polarity "
        f"{f(cl[cl.set == 'h9_control'].LR_polarity_oriented.min())} to {f(cl[cl.set == 'h9_control'].LR_polarity_oriented.max())}"
      + (f"; the supported transitions' LR ({', '.join(f'{u}: {f(r.LR_charge)}/{f(r.LR_polarity)}' for u, r in sup.iterrows())}, "
         "charge/polarity) lie inside them, so they do not stand out from unicellular controls" if len(sup) and all(
             (cl[cl.set == 'h9_control'].LR_charge_oriented.min() <= r.LR_charge <= cl[cl.set == 'h9_control'].LR_charge_oriented.max())
             for _, r in sup.iterrows()) else "")
      + ". \"Relaxed constraint\" can be worded for a transition only where its H9 r > 0 and its LR > 0 both hold and stand out from the control range "
      "above. Otherwise H9 stays a \"rate shift\".\n")
    a("**Caveats.**")
    a("- **Divergence bias.** Parsimony undercounts multiple hits, and in real proteins conservative replacements saturate first, so a faster or longer side can show a "
      "higher radical fraction without any relaxation. The H9 controls share this bias, so compare with their range rather than with 0.")
    a("- **Site opportunities.** Opportunities are not normalised (the same alignment columns are used on both sides).")
    a("- **Polarity scheme.** G and C are counted as polar (Zhang 2000).")
    with open(os.path.join(OUT, "summary.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
