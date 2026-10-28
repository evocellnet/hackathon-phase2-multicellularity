"""mcell -- analysis code of the multicellularity comparative-genomics project.

Modules
-------
hogs      Streaming OrthoXML (Ancestra HOGs) -> ancestral gene sets and per-branch
          gain / loss / duplication tables.
simtruth  Independent per-branch event counts from the true gene trees of
          ``ancestra.simulate`` (used to validate :mod:`mcell.hogs`).
gainladder Gain-evidence ladder (DESIGN §3.2): representatives of HOG gains, one
          DIAMOND search, clade-filtered sequence novelty, split repair, Pfam
          architecture novelty.
branchstats  Branch table, control pool U, matching, relabelling test, effect
          sizes, REML pooling, Holm, power/MDE and the analysis registry
          (DESIGN.md §4.6; docs/branchstats.md).
"""

__version__ = "0.1.0"
