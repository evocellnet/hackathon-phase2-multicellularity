# Shrink per-event log fold-changes (pseudocount LFC and NB SE from 05_scan.py) with ashr (normal mixture, method "shrink").
suppressPackageStartupMessages(library(ashr))
W <- "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/results/r2/w4/scan"
set.seed(11)
for (f in list.files(W, pattern = "^event_.*\\.tsv$")) {
  o <- file.path(W, sub("^event_", "ashr_", f)); if (file.exists(o)) stop("REFUSE: ", o)
  d <- read.delim(file.path(W, f), stringsAsFactors = FALSE)
  a <- ash(d$lfc, d$se, mixcompdist = "normal", method = "shrink")
  r <- data.frame(pfam = d$pfam, lfc_ash = a$result$PosteriorMean, lfc_ash_sd = a$result$PosteriorSD,
                  lfsr = a$result$lfsr, svalue = a$result$svalue)
  write.table(r, o, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(f, nrow(d), "lfsr<0.05:", sum(r$lfsr < 0.05), "\n")
}
