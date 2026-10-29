# Sensitivity: phylogenetic Poisson-lognormal GLMM (MCMCglmm, family "poisson" = log-normal overdispersion via units)
# on per-family counts. Usage: Rscript 04_glmm.R <set: P2_elm|P3_celldeath> <seed> <nitt>
# count_sf ~ family intercepts + log_size + disorder_frac + multicellular; random: phylogenetic species effect,
# non-phylogenetic species effect, family-specific deviation of the multicellular slope; residual = units (overdispersion).
suppressPackageStartupMessages({library(ape); library(MCMCglmm)})
a <- commandArgs(TRUE); S <- a[1]; seed <- as.integer(a[2]); nitt <- as.integer(a[3])
P <- "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity"
W <- file.path(P, "results/r2/w4")
f_out <- file.path(W, sprintf("glmm/glmm_%s_seed%d.tsv", S, seed))
if (file.exists(f_out)) stop("REFUSE: exists ", f_out)
dir.create(dirname(f_out), showWarnings = FALSE)
set.seed(seed)
sp  <- read.delim(file.path(W, "pgls_input/species.tsv"), stringsAsFactors = FALSE)
cnt <- read.delim(file.path(W, "pgls_input/counts_sets.tsv"), row.names = 1, check.names = FALSE)
fam <- read.delim(file.path(W, "pgls_input/set_families.tsv"), stringsAsFactors = FALSE)
fams <- if (S == "P2_elm") fam$pfam[fam$set == S & fam$subset == "primary"] else fam$pfam[fam$set == S & fam$widespread == 1]
prim <- sp$code[sp$primary == 1]
tree <- keep.tip(read.tree(file.path(P, "data/species_tree_bl_euk.nwk")), prim)
tree$node.label <- NULL
tree$edge.length <- tree$edge.length / max(node.depth.edgelength(tree))  # max root-to-tip = 1 (tree not ultrametric)
Ainv <- inverseA(tree, nodes = "TIPS", scale = FALSE)$Ainv
d <- do.call(rbind, lapply(fams, function(f) data.frame(code = prim, fam = f, count = cnt[prim, f])))
d <- merge(d, sp[, c("code", "multicellular", "n_proteins", "disorder_frac")], by = "code")
d$log_size <- log(d$n_proteins) - mean(log(sp$n_proteins[sp$primary == 1]))
d$phylo <- factor(d$code); d$species <- factor(d$code); d$fam <- factor(d$fam)
d$mc <- as.numeric(d$multicellular)
pr <- list(R = list(V = 1, nu = 0.002),
           G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                    G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000),
                    G3 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)))
t0 <- Sys.time()
m <- MCMCglmm(count ~ fam + log_size + disorder_frac + mc - 1, random = ~ phylo + species + idh(mc):fam,
              ginverse = list(phylo = Ainv), family = "poisson", data = d, prior = pr,
              nitt = nitt, burnin = nitt %/% 6, thin = max(1, (nitt - nitt %/% 6) %/% 1000), verbose = FALSE, pr = FALSE)
b <- as.numeric(m$Sol[, "mc"])
vc <- colMeans(m$VCV)
res <- data.frame(set = S, seed = seed, nitt = nitt, n_obs = nrow(d), n_families = length(fams), n_species = length(prim),
                  post_mean_mc = mean(b), hpd_lo = HPDinterval(m$Sol[, "mc"])[1], hpd_hi = HPDinterval(m$Sol[, "mc"])[2],
                  p_one_post = mean(b <= 0), ess_mc = as.numeric(effectiveSize(m$Sol[, "mc"])),
                  fold_change = exp(mean(b)),
                  vc_phylo = vc["phylo"], vc_species = vc["species"], vc_mc_by_family = vc[grep("mc", names(vc))[1]], vc_units = vc["units"],
                  minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
write.table(res, f_out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res)
