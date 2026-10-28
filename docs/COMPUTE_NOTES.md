# Compute timing notes (measured on this node, 24 CPU)
- Ancestra unit tests: 147 passed / 6 skipped in ~45 s (after env fix).
- pyhmmer 0.12.3 hmmsearch, Pfam-A 38.2 (30,134 families, GA cutoffs): 6,067 yeast proteins in 40.7 s on 24 CPUs → ~150 proteins/s → ~5.5 h per 3 M proteins.
- metapredict 3.1.01 (CPU torch, 4 threads): 2,000 random proteins (100-800 aa) in 5.1 s → ~400 proteins/s.
