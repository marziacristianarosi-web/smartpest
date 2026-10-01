# Script di verifica (sessione di analisi Corythucha). Il file dati "dati.xlsx" va sostituito con il percorso dell Excel.
suppressMessages({library(readxl); library(glmmTMB)}); set.seed(11)
d <- as.data.frame(read_excel("dati.xlsx")); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
m <- d[!duplicated(key), c("Season","Host_species","N_pianta")]; rownames(m) <- unique(key); m <- m[rownames(pa),]
m$R <- rowSums(pa)
S <- m[m$Host_species=="Quercus_robur",]; S$Season <- factor(S$Season); S$pl <- factor(S$N_pianta)   # domanda 1
H <- m[m$Season=="Summer",]; H$Host <- factor(H$Host_species)                                          # domanda 4
fams <- list(poisson=poisson(), genpois=genpois(), compois=compois())
cat("== 1. Adattamento: scelta della distribuzione (stagione, GLMM con pianta casuale) ==\n")
for (f in names(fams)) { fit <- try(glmmTMB(R ~ Season + (1|pl), data=S, family=fams[[f]]), silent=TRUE)
  if (inherits(fit,"try-error")) { cat(f, ": errore\n"); next }
  ok <- fit$fit$convergence==0 && !any(is.na(sqrt(diag(vcov(fit)$cond))))
  cat(sprintf("%-8s AIC=%.1f  convergenza=%s  var. effetto pianta=%.4f\n", f, AIC(fit), ok, VarCorr(fit)$cond$pl[1])) }
cat("\n== 2. Adattamento: ospite (GLM, una riga per pianta) ==\n")
for (f in names(fams)) { fit <- glmmTMB(R ~ Host, data=H, family=fams[[f]]); cat(sprintf("%-8s AIC=%.1f\n", f, AIC(fit))) }
# Calibrazione dell'errore di I tipo con dati simulati SENZA effetto
calib <- function(dat, form_full, form_null, fam, nsim) {
  null <- glmmTMB(form_null, data=dat, family=fam); sims <- simulate(null, nsim=nsim)
  res <- t(sapply(sims, function(y){ dd <- dat; dd$R <- y
    f1 <- try(suppressWarnings(glmmTMB(form_full, data=dd, family=fam)), silent=TRUE)
    f0 <- try(suppressWarnings(glmmTMB(form_null, data=dd, family=fam)), silent=TRUE)
    if (inherits(f1,"try-error")||inherits(f0,"try-error")||f1$fit$convergence!=0) return(c(NA,NA))
    pw <- summary(f1)$coefficients$cond[2,4]; lr <- 2*(logLik(f1)-logLik(f0)); c(pw, pchisq(as.numeric(lr),1,lower.tail=FALSE)) }))
  c(falliti=mean(is.na(res[,1])), Wald=mean(res[,1]<0.05,na.rm=TRUE), LRT=mean(res[,2]<0.05,na.rm=TRUE)) }
cat("\n== 3. Falsi positivi attesi 5% (dati simulati senza effetto) ==\n")
for (f in c("poisson","genpois")) {
  r1 <- calib(S, R ~ Season + (1|pl), R ~ 1 + (1|pl), fams[[f]], 500)
  r4 <- calib(H, R ~ Host, R ~ 1, fams[[f]], 500)
  cat(sprintf("%-8s Domanda 1 (GLMM, 10 piante x 2): falliti %.0f%%, Wald %.1f%%, LRT %.1f%% | Domanda 4 (GLM, 10+10): falliti %.0f%%, Wald %.1f%%, LRT %.1f%%\n",
     f, 100*r1[1], 100*r1[2], 100*r1[3], 100*r4[1], 100*r4[2], 100*r4[3])) }
