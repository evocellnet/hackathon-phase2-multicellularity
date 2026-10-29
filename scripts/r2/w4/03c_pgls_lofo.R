# Exploratory (post hoc, added after seeing that NB-ARC has the largest per-family effect in P3):
# leave-one-family-out for the P2/P3/P4 composites, and P3 without the NLR pair (NB-ARC + NACHT).
suppressPackageStartupMessages({library(ape); library(phylolm)})
P <- "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
W <- file.path(P, "results/r2/w4")
f_out <- file.path(W, "pgls/pgls_leave_one_family_out.tsv"); if (file.exists(f_out)) stop("REFUSE")
sp  <- read.delim(file.path(W, "pgls_input/species.tsv"), stringsAsFactors = FALSE); rownames(sp) <- sp$code
cnt <- read.delim(file.path(W, "pgls_input/counts_sets.tsv"), row.names = 1, check.names = FALSE)
fam <- read.delim(file.path(W, "pgls_input/set_families.tsv"), stringsAsFactors = FALSE)
tree <- read.tree(file.path(P, "data/species_tree_bl_euk.nwk"))
sp$log_size <- log(sp$n_proteins)
prim <- sp$code[sp$primary == 1]; phy <- keep.tip(tree, prim)
fit <- function(fams) {
  x <- log(as.matrix(cnt[prim, fams, drop = FALSE]) + 0.5); x <- x[, apply(x, 2, sd) > 0, drop = FALSE]
  d <- sp[prim, ]; d$y <- rowMeans(scale(x))
  s <- summary(phylolm(y ~ multicellular + log_size + disorder_frac, data = d, phy = phy, model = "lambda"))
  c(estimate = s$coefficients["multicellular", "Estimate"], se = s$coefficients["multicellular", "StdErr"],
    p_one = pt(s$coefficients["multicellular", "t.value"], nrow(d) - 4, lower.tail = FALSE), lambda = s$optpar)
}
sets <- list(P2_elm = fam$pfam[fam$set == "P2_elm" & fam$subset == "primary"],
             P3_celldeath = fam$pfam[fam$set == "P3_celldeath" & fam$widespread == 1],
             P4_tf = fam$pfam[fam$set == "P4_tf" & fam$widespread == 1])
rows <- list()
for (S in names(sets)) {
  fs <- sets[[S]]
  for (f in c("none", fs)) {
    r <- fit(setdiff(fs, f))
    rows[[length(rows) + 1]] <- data.frame(set = S, dropped = f, dropped_name = if (f == "none") "" else fam$name[fam$pfam == f][1], t(r))
  }
}
r <- fit(setdiff(sets$P3_celldeath, c("PF00931", "PF05729")))
rows[[length(rows) + 1]] <- data.frame(set = "P3_celldeath", dropped = "PF00931+PF05729", dropped_name = "NB-ARC+NACHT (NLR pair)", t(r))
R <- do.call(rbind, rows)
write.table(R, f_out, sep = "\t", quote = FALSE, row.names = FALSE)
print(R[R$set == "P3_celldeath", ])
for (S in names(sets)) { x <- R[R$set == S & R$dropped != "none", ]; cat(S, "LOFO: all positive =", all(x$estimate > 0), " max p =", max(x$p_one), " min est =", min(x$estimate), "\n") }
