# Script di verifica (sessione di analisi Corythucha). Il file dati "dati.xlsx" va sostituito con il percorso dell Excel.
suppressMessages({library(readxl); library(glmmTMB)})
d <- as.data.frame(read_excel("dati.xlsx")); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
m <- d[!duplicated(key), c("Season","Host_species","N_pianta")]; rownames(m) <- unique(key); m <- m[rownames(pa),]; m$R <- rowSums(pa)
S <- m[m$Host_species=="Quercus_robur",]; S <- S[order(S$N_pianta, S$Season),]; S$pl <- factor(S$N_pianta)
lrt <- function(dd, fam) { a <- try(suppressWarnings(glmmTMB(R ~ g + (1|pl), data=dd, family=fam)), silent=TRUE)
  b <- try(suppressWarnings(glmmTMB(R ~ 1 + (1|pl), data=dd, family=fam)), silent=TRUE)
  if (inherits(a,"try-error")||inherits(b,"try-error")) return(NA); as.numeric(2*(logLik(a)-logLik(b))) }
piante <- levels(S$pl); combo <- as.matrix(expand.grid(rep(list(0:1), 10)))
for (fam in list(genpois(), compois())) {
  L <- apply(combo, 1, function(sw) { dd <- S; for (i in which(sw==1)) { k <- dd$pl==piante[i]; dd$Season[k] <- rev(dd$Season[k]) }
    dd$g <- factor(dd$Season); lrt(dd, fam) })
  cat(sprintf("%-8s permutazioni valide %d/1024 | quota con p<0,05 (chi2): %.1f%% | 95o percentile LRT %.2f (chi2 teorico 3,84)\n",
      fam$family, sum(!is.na(L)), 100*mean(pchisq(L,1,lower.tail=FALSE) < 0.05, na.rm=TRUE), quantile(L, 0.95, na.rm=TRUE))) }
