# R2-P2 / R2-P3 confirmatory PGLS and sensitivities; exploratory P4 and per-family PGLS.
# Specification: DESIGN_ROUND2.md §6 and results/r2/w4/P3_rule.txt (written before this script was run).
suppressPackageStartupMessages({library(ape); library(phylolm)})
P <- "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
W <- file.path(P, "results/r2/w4")
outp <- function(rel) { f <- file.path(W, rel); if (file.exists(f)) stop("REFUSE: exists ", f)
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE); f }
set.seed(11)
sp  <- read.delim(file.path(W, "pgls_input/species.tsv"), stringsAsFactors = FALSE)
cnt <- read.delim(file.path(W, "pgls_input/counts_sets.tsv"), row.names = 1, check.names = FALSE)
fam <- read.delim(file.path(W, "pgls_input/set_families.tsv"), stringsAsFactors = FALSE)
tr  <- read.delim(file.path(P, "data/transitions.tsv"), stringsAsFactors = FALSE)
tr  <- tr[tr$kingdom == "euk", ]
po  <- tr[tr$primary_eligible == 1, ]
tree <- read.tree(file.path(P, "data/species_tree_bl_euk.nwk"))
rownames(sp) <- sp$code
sp$log_size <- log(sp$n_proteins)
stopifnot(all(sp$code %in% tree$tip.label))

famsets <- list(
  P2_elm = list(primary = fam$pfam[fam$set == "P2_elm" & fam$subset == "primary"],
                rule    = fam$pfam[fam$set == "P2_elm" & fam$subset == "primary" & fam$widespread == 1],
                repeats = fam$pfam[fam$set == "P2_elm"]),
  P3_celldeath = list(primary = fam$pfam[fam$set == "P3_celldeath" & fam$widespread == 1],
                      all     = fam$pfam[fam$set == "P3_celldeath"]),
  P4_tf = list(primary = fam$pfam[fam$set == "P4_tf" & fam$widespread == 1],
               all     = fam$pfam[fam$set == "P4_tf" & fam$absent_primary == 0]))

# response: equal-weight mean of per-family z-scored log(count+0.5), standardised over `codes`
make_y <- function(codes, fams, raw = FALSE) {
  x <- log(as.matrix(cnt[codes, fams, drop = FALSE]) + 0.5)
  s <- apply(x, 2, sd); x <- x[, s > 0, drop = FALSE]
  if (raw) return(list(y = rowMeans(x), nfam = ncol(x)))
  z <- scale(x)
  list(y = rowMeans(z), nfam = ncol(x))
}

fit1 <- function(d, form, term, boot = 0, model = "lambda") {
  phy <- keep.tip(tree, d$code)
  rownames(d) <- d$code
  f <- phylolm(as.formula(form), data = d, phy = phy, model = model, boot = boot)
  s <- summary(f)$coefficients
  est <- s[term, "Estimate"]; se <- s[term, "StdErr"]; tv <- s[term, "t.value"]
  dfree <- nrow(d) - nrow(s)
  out <- data.frame(term = term, estimate = est, se = se, t = tv, df = dfree,
                    ci_lo = est - qt(0.975, dfree) * se, ci_hi = est + qt(0.975, dfree) * se,
                    p_two = s[term, "p.value"], p_one = pt(tv, dfree, lower.tail = FALSE),
                    lambda = ifelse(is.null(f$optpar), NA, f$optpar), n = nrow(d),
                    n_mc = sum(d[[sub(" .*", "", term)]] == 1), logLik = f$logLik,
                    boot_ci_lo = NA, boot_ci_hi = NA, boot_lambda_lo = NA, boot_lambda_hi = NA)
  if (boot > 0) {
    bt <- f$bootstrap
    out$boot_ci_lo <- quantile(bt[, term], 0.025, na.rm = TRUE); out$boot_ci_hi <- quantile(bt[, term], 0.975, na.rm = TRUE)
    if ("lambda" %in% colnames(bt)) {
      out$boot_lambda_lo <- quantile(bt[, "lambda"], 0.025, na.rm = TRUE); out$boot_lambda_hi <- quantile(bt[, "lambda"], 0.975, na.rm = TRUE) }
  }
  out
}

prim <- sp$code[sp$primary == 1]
allc <- sp$code
FORM <- "y ~ multicellular + log_size + disorder_frac"
res <- list(); resp <- list(); k <- 0
add <- function(df, set, famset, variant, nfam) {
  k <<- k + 1; res[[k]] <<- cbind(data.frame(set = set, famset = famset, variant = variant, n_families = nfam), df) }

for (S in names(famsets)) for (FS in names(famsets[[S]])) {
  fams <- famsets[[S]][[FS]]
  main <- (FS == "primary")
  yy <- make_y(prim, fams); d <- sp[prim, ]; d$y <- yy$y[prim]
  resp[[paste(S, FS, "primary148", sep = "|")]] <- data.frame(set = S, famset = FS, species_set = "primary148", code = prim, y = d$y)
  add(fit1(d, FORM, "multicellular", boot = if (main) 1000 else 0), S, FS, "main", yy$nfam)
  if (!main) next
  # model forms
  add(fit1(d, FORM, "multicellular", model = "BM"), S, FS, "BM (lambda=1)", yy$nfam)
  add(fit1(d, "y ~ multicellular + log_size", "multicellular"), S, FS, "no disorder covariate", yy$nfam)
  # interpretable scale: mean log(count+0.5), not standardised -> exp(coef) = geometric-mean fold change
  yr <- make_y(prim, fams, raw = TRUE); dr <- d; dr$y <- yr$y[prim]
  add(fit1(dr, FORM, "multicellular"), S, FS, "unstandardised mean log count", yr$nfam)
  # leave-one-origin-out (standardisation of the primary set kept; rows dropped)
  for (i in seq_len(nrow(po))) {
    drop <- strsplit(po$mc_species[i], ",")[[1]]
    dd <- d[!(d$code %in% drop), ]
    add(fit1(dd, FORM, "multicellular"), S, FS, paste0("LOO -", po$transition_id[i], " (", sum(d$code %in% drop), " sp dropped)"), yy$nfam)
  }
  # no grey-zone species
  add(fit1(d[d$mc_state != "grey", ], FORM, "multicellular"), S, FS, "grey-zone excluded", yy$nfam)
  # all 187 proteomes, re-standardised over 187
  ya <- make_y(allc, fams); da <- sp; da$y <- ya$y[allc]
  da$transcriptome <- as.integer(da$data_type == "transcriptome")
  resp[[paste(S, FS, "all187", sep = "|")]] <- data.frame(set = S, famset = FS, species_set = "all187", code = allc, y = da$y)
  add(fit1(da, FORM, "multicellular"), S, FS, "all 187 proteomes", ya$nfam)
  add(fit1(da, paste(FORM, "+ transcriptome + busco_C"), "multicellular"), S, FS, "all 187 + transcriptome & BUSCO covariates", ya$nfam)
  # complex vs simple among multicellular species
  dm <- d[d$multicellular == 1, ]
  add(fit1(dm, "y ~ complex_DT + log_size + disorder_frac", "complex_DT"), S, FS, "complex(D or T) vs simple, MC species only", yy$nfam)
  add(fit1(dm, "y ~ complex_T + log_size + disorder_frac", "complex_T"), S, FS, "complex(T only) vs simple, MC species only", yy$nfam)
  # three-level state: simple and complex vs unicellular
  d3 <- d; d3$simple <- as.integer(d3$multicellular == 1 & d3$complex_DT == 0); d3$complex <- d3$complex_DT
  r1 <- fit1(d3, "y ~ simple + complex + log_size + disorder_frac", "simple")
  r2 <- fit1(d3, "y ~ simple + complex + log_size + disorder_frac", "complex")
  add(r1, S, FS, "3-state: simple vs uni", yy$nfam); add(r2, S, FS, "3-state: complex(D or T) vs uni", yy$nfam)
}
R <- do.call(rbind, res)
write.table(R, outp("pgls/pgls_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, resp), outp("pgls/responses.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# sister-clade contrasts per primary origin (descriptive; P3_rule item 7)
sc <- list()
for (S in names(famsets)) {
  fams <- famsets[[S]][["primary"]]
  for (ss in c("primary148", "all187")) {
    codes <- if (ss == "primary148") prim else allc
    y <- make_y(codes, fams)$y
    for (i in seq_len(nrow(po))) {
      m <- intersect(strsplit(po$mc_species[i], ",")[[1]], codes)
      s <- intersect(strsplit(po$sister_species[i], ",")[[1]], codes)
      sc[[length(sc) + 1]] <- data.frame(set = S, species_set = ss, origin = po$transition_id[i], n_mc = length(m), n_sister = length(s),
        mean_mc = if (length(m)) mean(y[m]) else NA, mean_sister = if (length(s)) mean(y[s]) else NA,
        diff = if (length(m) && length(s)) mean(y[m]) - mean(y[s]) else NA)
    }
  }
}
SC <- do.call(rbind, sc)
write.table(SC, outp("pgls/sister_contrasts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# per-family PGLS (exploratory), primary species set, log(count+0.5) response
pf <- list()
for (i in seq_len(nrow(fam))) {
  f <- fam$pfam[i]; x <- log(cnt[prim, f] + 0.5)
  if (sd(x) == 0) next
  d <- sp[prim, ]; d$y <- x
  r <- tryCatch(fit1(d, FORM, "multicellular"), error = function(e) NULL)
  if (is.null(r)) next
  pf[[length(pf) + 1]] <- cbind(fam[i, c("set", "pfam", "name", "subset", "widespread", "frac_present_primary")], r)
}
PF <- do.call(rbind, pf)
PF$q_bh_within_set <- ave(PF$p_one, PF$set, FUN = function(p) p.adjust(p, "BH"))
write.table(PF, outp("pgls/per_family_pgls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("done\n")
print(R[R$variant %in% c("main"), c("set", "famset", "n_families", "estimate", "se", "p_one", "lambda", "n", "n_mc")])
