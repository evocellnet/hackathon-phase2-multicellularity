#!/bin/bash
# quick lookup in the merged source index: q.sh "regex" [db]
grep -i -E "$1" $(dirname $0)/../../data/raw/source_index.tsv | grep -E "^${2:-}" | cut -f1-7 | sort -t$'\t' -k1,1 -k6,6gr | head -${N:-25}
