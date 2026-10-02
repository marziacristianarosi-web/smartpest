# Script di verifica: GLMM Poisson con lme4 (domanda 1). Il file "dati.xlsx" va sostituito con il percorso dell Excel.
suppressMessages({library(readxl); library(glmmTMB); library(lme4)}); set.seed(11)
d <- as.data.frame(read_excel("dati.xlsx")); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
m <- d[!duplicated(key), c("Season","Host_species","N_pianta")]; rownames(m) <- unique(key); m <- m[rownames(pa),]; m$R <- rowSums(pa)
S <- m[m$Host_species=="Quercus_robur",]; S$g <- factor(S$Season); S$pl <- factor(S$N_pianta)
a <- glmer(R ~ g + (1|pl), data=S, family=poisson); b <- glmmTMB(R ~ g + (1|pl), data=S, family=poisson)
cat("logLik glmer", logLik(a), " glmmTMB", logLik(b), "\n")
cat("var pianta glmer", VarCorr(a)$pl[1], " singular", isSingular(a), "\n")
cat("dispersione Pearson", sum(residuals(a,"pearson")^2)/df.residual(a), "\n")
a0 <- update(a, . ~ . - g); L0 <- 2*(logLik(a)-logLik(a0))
# tasso di errore di I tipo del LRT Poisson glmer sotto la permutazione esatta entro pianta (1024 combinazioni)
E <- as.matrix(expand.grid(rep(list(0:1),10))); lev <- levels(S$g); pls <- levels(S$pl)
L <- apply(E,1,function(e){ dd <- S; for(i in 1:10){ if(e[i]==1){ j <- dd$pl==pls[i]; dd$g[j] <- rev(dd$g[j]) } }
  f1 <- suppressMessages(suppressWarnings(glmer(R ~ g + (1|pl), data=dd, family=poisson))); f0 <- suppressMessages(suppressWarnings(glmer(R ~ 1 + (1|pl), data=dd, family=poisson)))
  as.numeric(2*(logLik(f1)-logLik(f0)))})
cat(sprintf("glmer Poisson: %% permutazioni con p<0.05 = %.1f ; quantile 95%% LRT = %.2f (atteso 3.84)\n", 100*mean(pchisq(L,1,lower.tail=FALSE)<0.05), quantile(L,.95)))
