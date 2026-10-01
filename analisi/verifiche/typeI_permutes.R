# Script di verifica (sessione di analisi Corythucha). Il file dati "dati.xlsx" va sostituito con il percorso dell Excel.
suppressMessages({library(readxl); library(glmmTMB); library(permutes)}); set.seed(3)
d <- as.data.frame(read_excel("dati.xlsx")); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
m <- d[!duplicated(key), c("Season","Host_species","N_pianta")]; rownames(m) <- unique(key); m <- m[rownames(pa),]; m$R <- rowSums(pa)
S <- m[m$Host_species=="Quercus_robur",]; S <- S[order(S$N_pianta, S$Season),]; S$pl <- factor(S$N_pianta)
# Dati SENZA effetto stagione con la variabilita' reale: entro ogni pianta si assegna a caso quale valore e' "estate"
nd <- 60; p <- rep(NA, nd)
for (i in 1:nd) { dd <- S; for (pp in levels(dd$pl)) if (runif(1) < .5) { k <- dd$pl==pp; dd$Season[k] <- rev(dd$Season[k]) }
  dd$g <- factor(dd$Season)
  r <- try(suppressWarnings(perm.glmmTMB(R ~ g + (1|pl), data=dd, family=genpois(), nperm=200, progress=FALSE)), silent=TRUE)
  if (!inherits(r,"try-error")) p[i] <- r$p[r$Factor=="gSummer"]
  if (i %% 10 == 0) cat(i, "dataset; quota p<0,05 finora:", round(mean(p[1:i] < .05, na.rm=TRUE),3), "\n") }
cat(sprintf("FINALE permutes::perm.glmmTMB: dataset validi %d/%d | falsi positivi (p<0,05) %.1f%% (atteso 5%%) | p mediano %.2f\n", sum(!is.na(p)), nd, 100*mean(p<.05,na.rm=TRUE), median(p,na.rm=TRUE)))
