import sys, os, warnings
root, P, out = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, root)
warnings.simplefilter("ignore")
import pandas as pd
from mcell import branchstats as bs
print("using", bs.__file__)
tr = pd.read_csv(f"{P}/data/transitions.tsv", sep="\t"); traits = pd.read_csv(f"{P}/data/species_traits.tsv", sep="\t")
tree = bs.load_tree(f"{P}/data/guide_tree_euk.nwk")
bl = bs.branch_lengths_from_tree(tree, f"{P}/data/species_tree_bl_euk.nwk")
bt = bs.build_branch_table(tree, traits, tr[tr.kingdom == "euk"], branch_lengths=bl)
pool = bs.h9_control_branches(tree, bt)
bt["h9_pool"] = bt.index.isin(pool)
bt[["branch_class","transition_id","role","stem_path","branch_length","in_U","U_exclusion","adjacent_to_stem","h9_pool"]].to_csv(out, sep="\t")
print("U", int(bt.in_U.sum()), "h9 pool", pool)
