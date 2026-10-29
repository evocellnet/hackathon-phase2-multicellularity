# Independent verifier for R2-P2 / R2-P3 (not using producer code).
suppressMessages({library(phylolm); library(ape)})
P <- "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
OUT <- file.path(P, "results/r2/verification/w4")
log_ <- file(file.path(OUT, "v2_pgls_output.txt"), "w")
say <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = log_) }

M  <- read.delim(gzfile(file.path(OUT, "count_matrix.tsv.gz")), row.names = 1, check.names = FALSE)
cv <- read.delim(file.path(OUT, "species_cov.tsv"))
tr <- read.delim(file.path(P, "data/species_traits.tsv"))
tx <- read.delim(file.path(P, "data/transitions.tsv"))
tree <- read.tree(file.path(P, "data/species_tree_bl_euk.nwk"))
cd  <- read.delim(file.path(P, "data/r2_sets/cell_death.tsv"))
elm <- read.delim(file.path(P, "data/r2_sets/elm_motif_binding.tsv"))

tr <- tr[tr$kingdom == "euk", ]
S <- tr$code[tr$data_type == "genome" & tr$busco_C >= 50]
say("Primary species set n = ", length(S))
cnt <- function(fams, sp = rownames(M)) {
  out <- matrix(0, length(sp), length(fams), dimnames = list(sp, fams))
  ok <- intersect(fams, colnames(M)); out[, ok] <- as.matrix(M[sp, ok]); out }

orig <- tx[tx$primary_eligible == 1, ]
mcl  <- setNames(strsplit(orig$mc_species, ","), orig$transition_id)
sisl <- setNames(strsplit(orig$sister_species, ","), orig$transition_id)
n_orig <- function(fams, cl) { C <- cnt(fams); sapply(fams, function(f) sum(sapply(cl, function(s) any(C[s, f] > 0)))) }
prevS  <- function(fams) colMeans(cnt(fams, S) > 0)

rule <- function(fams, symmetric = FALSE) {
  a <- n_orig(fams, mcl) >= 3; b <- prevS(fams) >= 0.5
  keep <- a | b
  if (symmetric) keep <- keep | n_orig(fams, sisl) >= 3
  fams[keep & prevS(fams) > 0] }

d0 <- data.frame(code = S)
d0$multicellular <- tr$multicellular[match(S, tr$code)]
d0$log_size <- log(cv$n_proteins[match(S, cv$code)])
d0$disorder <- cv$disorder[match(S, cv$code)]
rownames(d0) <- S

response <- function(fams, sp) {
  X <- log(cnt(fams, sp) + 0.5)
  Z <- scale(X)            # per-family z across the species given
  rowMeans(Z) }

fit <- function(y, sp, form = y ~ multicellular + log_size + disorder) {
  d <- d0[sp, ]; d$y <- y[sp]
  tt <- keep.tip(tree, sp)
  m <- phylolm(form, data = d, phy = tt, model = "lambda")
  s <- summary(m)$coefficients
  df <- nrow(d) - length(coef(m))
  tval <- s["multicellular", "t.value"]
  list(m = m, coef = s, lambda = m$optpar, n = nrow(d),
       b = s["multicellular", "Estimate"], se = s["multicellular", "StdErr"],
       p1 = pt(tval, df, lower.tail = FALSE), df = df) }
rep_fit <- function(lbl, f) {
  ci <- f$b + c(-1, 1) * qt(0.975, f$df) * f$se
  say(sprintf("%-55s n=%d b=%+.4f SE=%.4f 95%%CI(t)=[%.3f,%.3f] p1=%.5f lambda=%.3f", lbl, f$n, f$b, f$se, ci[1], ci[2], f$p1, f$lambda)) }

## ---- Task 1: P3 ----
cdf <- cd$pfam
say("Cell-death families in file: ", length(cdf), "; present in >=1 of 187: ", sum(colSums(cnt(cdf)) > 0),
    "; present in >=1 of S: ", sum(prevS(cdf) > 0))
tab <- data.frame(pfam = cdf, name = cd$name, n_mc_origins = n_orig(cdf, mcl), n_sister_origins = n_orig(cdf, sisl),
                  prev_S = round(prevS(cdf), 3), total_S = colSums(cnt(cdf, S)))
tab$rule_primary <- tab$n_mc_origins >= 3 | tab$prev_S >= 0.5
say(paste(capture.output(print(tab, row.names = FALSE)), collapse = "\n"))
fam12 <- rule(cdf)
say("Widespread (rule) families: ", length(fam12), ": ", paste(fam12, collapse = ","))
y3 <- response(fam12, S); f3 <- fit(y3, S)
rep_fit("P3 primary (rule families)", f3)
say(paste(capture.output(print(f3$coef)), collapse = "\n"))

## ---- Task 1: P2 ----
p2f <- elm$pfam[elm$set == "primary"]
p2p <- p2f[prevS(p2f) > 0]
say("P2 primary families: ", length(p2f), "; present in S: ", length(p2p))
f2 <- fit(response(p2p, S), S); rep_fit("P2 primary (all present primary families)", f2)
say(paste(capture.output(print(f2$coef)), collapse = "\n"))
p2r <- rule(p2f); f2r <- fit(response(p2r, S), S); rep_fit(sprintf("P2 sensitivity: rule families (%d)", length(p2r)), f2r)

## ---- Task 2 ----
famsym <- rule(cdf, symmetric = TRUE)
say("Symmetric rule families: ", length(famsym), ": ", paste(setdiff(famsym, fam12), collapse = ","), " added")
rep_fit(sprintf("P3 (a) symmetric rule (%d fams)", length(famsym)), fit(response(famsym, S), S))
famany <- cdf[prevS(cdf) > 0]
rep_fit(sprintf("P3 (b) all families present in S (%d fams)", length(famany)), fit(response(famany, S), S))
famany187 <- cdf[colSums(cnt(cdf)) > 0]
if (!setequal(famany187, famany)) say("note: families present in 187 but not in S: ", paste(setdiff(famany187, famany), collapse=","))
# per-family contributions (diagnostic)
pf <- t(sapply(famany, function(f) { g <- fit(response(f, S), S); c(b = g$b, p1 = g$p1) }))
say("Per-family PGLS (same covariates), diagnostic:")
say(paste(capture.output(print(round(cbind(pf, in_rule = famany %in% fam12, in_sym = famany %in% famsym), 4))), collapse = "\n"))

# (c) permutation of prevalence-matched random sets through rule + PGLS
set.seed(11)
allf <- colnames(M); prevAll <- colMeans(M[S, ] > 0)
pool <- setdiff(allf, cdf)
tgt  <- prevS(cdf)
draw_set <- function() {
  chosen <- character(0)
  for (f in cdf) {
    w <- 0.025
    repeat { cand <- setdiff(pool[abs(prevAll[pool] - tgt[f]) <= w], chosen); if (length(cand) >= 20) break; w <- w * 2 }
    chosen <- c(chosen, cand[sample.int(length(cand), 1)]) }
  chosen }
perm <- data.frame()
for (i in 1:200) {
  s <- draw_set(); r <- rule(s)
  if (length(r) < 1) { perm <- rbind(perm, data.frame(i = i, n_sel = 0, b = NA, p1 = NA, b_all = NA, p1_all = NA)); next }
  g <- fit(response(r, S), S); ga <- fit(response(s[prevAll[s] > 0], S), S)
  perm <- rbind(perm, data.frame(i = i, n_sel = length(r), b = g$b, p1 = g$p1, b_all = ga$b, p1_all = ga$p1))
}
write.table(perm, file.path(OUT, "perm_random_sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ok <- !is.na(perm$p1)
say(sprintf("Permutation (200 prevalence-matched random 22-family sets, seed 11): rule-selected n_sel median %d (range %d-%d)",
            as.integer(median(perm$n_sel)), min(perm$n_sel), max(perm$n_sel)))
say(sprintf("  rule+PGLS: frac p1<0.05 = %.3f (%d/%d), binomial p vs 5%% = %.3g; mean b = %+.4f; frac b>0 = %.3f; frac b >= observed %.3f = %.3f",
            mean(perm$p1[ok] < 0.05), sum(perm$p1[ok] < 0.05), sum(ok), binom.test(sum(perm$p1[ok] < 0.05), sum(ok), 0.05, "greater")$p.value,
            mean(perm$b[ok]), mean(perm$b[ok] > 0), f3$b, mean(perm$b[ok] >= f3$b)))
say(sprintf("  no-selection (all present) PGLS: frac p1<0.05 = %.3f; mean b = %+.4f; frac b_all >= observed = %.3f",
            mean(perm$p1_all < 0.05), mean(perm$b_all), mean(perm$b_all >= f3$b)))

## ---- Task 3 ----
say("P3 covariates: log_size b=", sprintf("%+.4f (p2=%.3g)", f3$coef["log_size", "Estimate"], f3$coef["log_size", "p.value"]),
    "; disorder b=", sprintf("%+.4f (p2=%.3g)", f3$coef["disorder", "Estimate"], f3$coef["disorder", "p.value"]))
rep_fit("P3 without covariates (y ~ multicellular)", fit(y3, S, y ~ multicellular))
rep_fit("P3 log_size only", fit(y3, S, y ~ multicellular + log_size))
big <- S[order(-d0[S, "log_size"])][1:10]
say("10 largest proteomes: ", paste(sprintf("%s(%d,mc=%d)", big, round(exp(d0[big, "log_size"])), d0[big, "multicellular"]), collapse = ", "))
S2 <- setdiff(S, big)
rep_fit("P3 excl. 10 largest, standardisation kept (full S)", fit(y3, S2))
rep_fit("P3 excl. 10 largest, re-standardised", fit(response(fam12, S2), S2))
rep_fit("P3 using proportion log(count/size) instead (diagnostic)", fit(rowMeans(scale(log((cnt(fam12, S) + 0.5) / exp(d0$log_size)))), S, y ~ multicellular + disorder))
close(log_)
