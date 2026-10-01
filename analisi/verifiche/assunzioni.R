# Script di verifica (sessione di analisi Corythucha). Il file dati "dati.xlsx" va sostituito con il percorso dell Excel.
suppressMessages({library(readxl); library(glmmTMB)}); set.seed(7)
d <- as.data.frame(read_excel("dati.xlsx")); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
m <- d[!duplicated(key), c("Season","Host_species","N_pianta")]; rownames(m) <- unique(key); m <- m[rownames(pa),]; m$R <- rowSums(pa)
S <- m[m$Host_species=="Quercus_robur",]; S$g <- factor(S$Season); S$pl <- factor(S$N_pianta)
H <- m[m$Season=="Summer",]; H$g <- factor(H$Host_species); H$pl <- factor(H$N_pianta)
# Residui quantili simulati (stesso principio di DHARMa, Hartig): 1000 simulazioni dal modello stimato
qres <- function(fit, y, nsim=1000) { sims <- as.matrix(simulate(fit, nsim=nsim))
  r <- sapply(seq_along(y), function(i) { lo <- mean(sims[i,] < y[i]); hi <- mean(sims[i,] <= y[i]); runif(1, lo, hi) })
  mu <- fitted(fit); disp_obs <- sum((y-mu)^2); disp_sim <- apply(sims, 2, function(s) sum((s-mu)^2))
  list(r=r, ks=ks.test(r, "punif")$p.value, disp_ratio=disp_obs/mean(disp_sim),
       disp_p=2*min(mean(disp_sim>=disp_obs), mean(disp_sim<=disp_obs)),
       out=sum(y > apply(sims,1,max) | y < apply(sims,1,min))) }
check <- function(nome, dat, f_full, f_het, re) {
  cat("\n==========", nome, "==========\n")
  fit <- glmmTMB(f_full, data=dat, family=genpois())
  q <- qres(fit, dat$R)
  cat(sprintf("A1 Distribuzione (residui uniformi, KS): p = %.3f\n", q$ks))
  cat(sprintf("A2 Dispersione residua (oss/sim = 1 atteso): %.2f, p = %.3f\n", q$disp_ratio, q$disp_p))
  cat(sprintf("A3 Valori fuori dall'intervallo simulato (outlier): %d su %d\n", q$out, nrow(dat)))
  cat("   Residui per gruppo (media dei residui quantili, attesa 0,5):", paste(names(tapply(q$r,dat$g,mean)), round(tapply(q$r,dat$g,mean),2), collapse=" | "), "\n")
  het <- glmmTMB(f_full, dispformula = ~ g, data=dat, family=genpois())
  lr <- as.numeric(2*(logLik(het)-logLik(fit)))
  cat(sprintf("A4 Dispersione uguale nei due gruppi: phi = %s | LRT p = %.3f | AIC omogeneo %.1f vs eterogeneo %.1f\n",
     paste(names(tapply(dat$R,dat$g,var)), round(exp(fixef(het)$disp[1]+c(0,fixef(het)$disp[2])),2), collapse=" / "),
     pchisq(lr,1,lower.tail=FALSE), AIC(fit), AIC(het)))
  # A5 osservazioni influenti: variazione del coefficiente togliendo una pianta alla volta
  b <- fixef(fit)$cond[2]; ids <- unique(dat$N_pianta)
  inf <- sapply(ids, function(p) { dd <- dat[dat$N_pianta!=p,]; f <- try(glmmTMB(f_full, data=dd, family=genpois()), silent=TRUE); if (inherits(f,"try-error")) NA else fixef(f)$cond[2] })
  cat(sprintf("A5 Influenza (rapporto delle medie, exp(coef)): completo %.2f | escludendo una pianta alla volta: min %.2f (senza %s), max %.2f (senza %s)\n",
     exp(b), exp(min(inf,na.rm=T)), ids[which.min(inf)], exp(max(inf,na.rm=T)), ids[which.max(inf)]))
  if (re) { u <- ranef(fit)$cond$pl[,1]; cat(sprintf("A6 Effetti casuali (10 livelli): varianza %.4f; normalita' Shapiro p = %.3f (potenza molto bassa con 10 valori)\n", VarCorr(fit)$cond$pl[1], shapiro.test(u)$p.value)) }
  cat(sprintf("A7 Limite superiore: valore massimo osservato %d su 14 taxa possibili; prob. simulata di superare 14: %.4f\n", max(dat$R), mean(as.matrix(simulate(fit, nsim=500))>14)))
}
check("Domanda 1 - stagione (GLMM)", S, R ~ g + (1|pl), NULL, TRUE)
check("Domanda 4 - ospite (GLM)", H, R ~ g, NULL, FALSE)
