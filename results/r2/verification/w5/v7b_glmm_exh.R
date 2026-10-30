suppressMessages({library(lme4)})
W <- "/scratch/spascare/agentic-ai/workspaces/phase2-multicellularity/results/r2/verification/w5/work"
X <- read.delim(file.path(W,"glmm_input_exh.tsv")); X$new <- X$new_exh
X$br <- paste(X$panel,X$branch_child)
X$llen <- log(X$len_m); X$lng <- log(X$n_genes); X$lprot <- log(X$n_out_prot); X$lsp <- log(X$n_out_sp)
X$nested <- as.integer(X$role=="event" & X$ekind=="nested"); X$origin <- as.integer(X$role=="event" & X$ekind=="origin")
sc <- function(v) (v-mean(v))/sd(v)
for (v in c("llen","dis_m","lc_m","lng","depth","lprot","lsp","frac_tx")) X[[paste0("z_",v)]] <- sc(X[[v]])
ctl <- glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))
rep <- function(m, term="event"){ s <- summary(m)$coefficients; e<-s[term,1]; se<-s[term,2]
  sprintf("%s OR %.2f [%.2f, %.2f] p=%.3g", term, exp(e), exp(e-1.96*se), exp(e+1.96*se), s[term,4]) }
out <- c()
m0 <- glmer(new ~ event + (1|panel) + (1|br), data=X, family=binomial, control=ctl)
out <- c(out, paste("GLMM unadjusted (random panel+branch):", rep(m0)))
m1 <- glmer(new ~ event + z_lprot + (1|panel) + (1|br), data=X, family=binomial, control=ctl)
out <- c(out, paste("GLMM + log outside proteins:", rep(m1)))
m2 <- glmer(new ~ event + z_depth + z_lprot + z_frac_tx + (1|panel) + (1|br), data=X, family=binomial, control=ctl)
out <- c(out, paste("GLMM + branch covariates (depth, log out-prot, frac_tx):", rep(m2)))
m3 <- glmer(new ~ event + z_llen + z_dis_m + z_lc_m + z_lng + (1|panel) + (1|br), data=X, family=binomial, control=ctl)
out <- c(out, paste("GLMM + family covariates (log len, disorder, LC, log size):", rep(m3)))
m4 <- glmer(new ~ event + z_depth + z_lprot + z_frac_tx + z_llen + z_dis_m + z_lc_m + z_lng + (1|panel) + (1|br), data=X, family=binomial, control=ctl)
out <- c(out, paste("GLMM full:", rep(m4)))
out <- c(out, "full model coefficients:", capture.output(print(round(summary(m4)$coefficients,4))))
m4b <- glmer(new ~ event + z_depth + z_lsp + z_frac_tx + z_llen + z_dis_m + z_lc_m + z_lng + (1|panel) + (1|br), data=X, family=binomial, control=ctl)
out <- c(out, paste("GLMM full, log outside SPECIES instead of proteins:", rep(m4b)))
m5 <- glmer(new ~ nested + origin + z_depth + z_lprot + z_frac_tx + z_llen + z_dis_m + z_lc_m + z_lng + (1|panel) + (1|br), data=X, family=binomial, control=ctl)
out <- c(out, paste("GLMM full by kind:", rep(m5,"nested"), "|", rep(m5,"origin")))
m5u <- glmer(new ~ nested + origin + (1|panel) + (1|br), data=X, family=binomial, control=ctl)
out <- c(out, paste("GLMM unadjusted by kind:", rep(m5u,"nested"), "|", rep(m5u,"origin")))
# restrict to families >= 150 aa and disorder < 0.5 ("foldable")
Y <- subset(X, len_m>=150 & dis_m<0.5)
out <- c(out, sprintf("Foldable subset (mean len>=150 aa, disorder<0.5): n=%d; event new %.3f (%d) vs control %.3f (%d)", nrow(Y),
   mean(Y$new[Y$event==1]), sum(Y$event==1), mean(Y$new[Y$event==0]), sum(Y$event==0)))
m6 <- glmer(new ~ event + (1|panel) + (1|br), data=Y, family=binomial, control=ctl)
out <- c(out, paste("  GLMM unadjusted:", rep(m6)))
m7 <- glmer(new ~ event + z_depth + z_lprot + z_frac_tx + z_llen + z_dis_m + z_lc_m + z_lng + (1|panel) + (1|br), data=Y, family=binomial, control=ctl)
out <- c(out, paste("  GLMM full:", rep(m7)))
m7k <- glmer(new ~ nested + origin + z_depth + z_lprot + z_frac_tx + z_llen + z_dis_m + z_lc_m + z_lng + (1|panel) + (1|br), data=Y, family=binomial, control=ctl)
out <- c(out, paste("  GLMM full by kind:", rep(m7k,"nested"), "|", rep(m7k,"origin")))
# branch-level regression: logit of branch fraction on covariate means, panel fixed effects, quasibinomial
B <- aggregate(cbind(new, llen, dis_m, lc_m, lng)~ br+panel+event+nested+origin+depth+lprot+frac_tx, data=X, FUN=mean)
B$n <- aggregate(new~br, data=X, FUN=length)$new[match(B$br, aggregate(new~br, data=X, FUN=length)$br)]
mb0 <- glm(new ~ event + panel, data=B, family=quasibinomial, weights=n)
mb <- glm(new ~ event + panel + depth + lprot + frac_tx + llen + dis_m + lc_m + lng, data=B, family=quasibinomial, weights=n)
rq <- function(m,term="event"){s<-summary(m)$coefficients; sprintf("%s OR %.2f [%.2f, %.2f] p=%.3g",term,exp(s[term,1]),exp(s[term,1]-1.96*s[term,2]),exp(s[term,1]+1.96*s[term,2]),s[term,4])}
out <- c(out, sprintf("Branch-level quasibinomial (n=%d branches), panel FE: unadjusted %s ; adjusted %s", nrow(B), rq(mb0), rq(mb)))
mb2 <- glm(new ~ nested + origin + panel + depth + lprot + frac_tx + llen + dis_m + lc_m + lng, data=B, family=quasibinomial, weights=n)
out <- c(out, paste("  branch-level adjusted by kind:", rq(mb2,"nested"), "|", rq(mb2,"origin")))
# covariate imbalance event vs control at branch level
for (v in c("depth","lprot","frac_tx","llen","dis_m","lc_m","lng")) { a<-B[[v]][B$event==1]; b<-B[[v]][B$event==0]
  out <- c(out, sprintf("  branch covariate %s: event mean %.3f, control mean %.3f, SMD %.2f", v, mean(a), mean(b), (mean(a)-mean(b))/sqrt((var(a)+var(b))/2))) }
writeLines(out); writeLines(out, file.path(W,"glmm_results_exh.txt"))
