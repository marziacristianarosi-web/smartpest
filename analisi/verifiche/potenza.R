# Script di verifica (sessione di analisi Corythucha). Il file dati "dati.xlsx" va sostituito con il percorso dell Excel.
suppressMessages({library(readxl); library(glmmTMB)}); set.seed(2026)
d <- as.data.frame(read_excel("dati.xlsx")); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
m <- d[!duplicated(key), c("Season","Host_species","N_pianta")]; rownames(m) <- unique(key); m <- m[rownames(pa),]; m$R <- rowSums(pa)
S <- m[m$Host_species=="Quercus_robur",]; S$g <- factor(S$Season); S$pl <- factor(S$N_pianta)
H <- m[m$Season=="Summer",]; H$g <- factor(H$Host_species)
# Parametri di "disturbo" (variabilita') dal modello completo; media generale come livello di base; effetto fissato a priori
fS <- glmmTMB(R ~ g + (1|pl), data=S, family=genpois()); fH <- glmmTMB(R ~ g, data=H, family=genpois())
par <- list(S=list(dat=S, M=mean(S$R), betad=fixef(fS)$disp, theta=fS$fit$par[names(fS$fit$par)=="theta"], re=TRUE),
            H=list(dat=H, M=mean(H$R), betad=fixef(fH)$disp, theta=NULL, re=FALSE))
lrt <- function(dd, re) { f1 <- if (re) R ~ g + (1|pl) else R ~ g; f0 <- if (re) R ~ 1 + (1|pl) else R ~ 1
  a <- try(suppressWarnings(glmmTMB(f1, data=dd, family=genpois())), silent=TRUE); b <- try(suppressWarnings(glmmTMB(f0, data=dd, family=genpois())), silent=TRUE)
  if (inherits(a,"try-error")||inherits(b,"try-error")) return(NA); as.numeric(2*(logLik(a)-logLik(b))) }
sim <- function(p, delta, nsim) { dd <- p$dat; mu0 <- p$M - delta/2; mu1 <- p$M + delta/2
  beta <- c(log(mu0), log(mu1) - log(mu0))
  form <- if (p$re) ~ g + (1|pl) else ~ g
  np <- list(beta=beta, betad=p$betad); if (p$re) np$theta <- p$theta
  ys <- simulate_new(form, newdata=dd, family=genpois(), newparams=np, nsim=nsim)
  sapply(ys, function(y){ dd$R <- y; lrt(dd, p$re) }) }
deltas <- c(0.5, 1, 1.5, 2, 2.5, 3)
out <- list()
for (k in c("S","H")) { p <- par[[k]]
  L0 <- sim(p, 0, 400); crit <- quantile(L0, 0.95, na.rm=TRUE)
  cat(sprintf("\n%s: media di base %.2f | dispersione phi %.2f | valore critico LRT (bootstrap) %.2f invece di 3.84 (chi2)\n", k, p$M, exp(p$betad), crit))
  for (dl in deltas) { L <- sim(p, dl, 200); pw <- mean(L > crit, na.rm=TRUE)
    cat(sprintf("  differenza %.1f taxa/pianta: potenza %.0f%% (IC95 %.0f-%.0f%%)\n", dl, 100*pw, 100*max(0,pw-1.96*sqrt(pw*(1-pw)/200)), 100*min(1,pw+1.96*sqrt(pw*(1-pw)/200)))) } }
