# Script di verifica (sessione di analisi Corythucha). Il file dati "dati.xlsx" va sostituito con il percorso dell Excel.
suppressMessages(library(readxl))
d <- as.data.frame(read_excel("dati.xlsx")); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
m <- d[!duplicated(key), c("Season","Host_species","N_pianta")]; rownames(m) <- unique(key); m <- m[rownames(pa),]
m$R <- rowSums(pa); m$G <- paste(m$Host_species, m$Season)
cat("== Ricchezza per pianta ==\n")
for (g in unique(m$G)) { r <- m$R[m$G==g]; n <- length(r)
  D <- (n-1)*var(r)/mean(r); p_low <- pchisq(D, n-1); 
  cat(sprintf("%-22s valori: %s | media %.1f | var %.2f | var/media %.2f | min-max %d-%d | test dispersione Poisson: chi2=%.2f, gl=%d, P(sottodisp.)=%.4f\n",
    g, paste(sort(r), collapse=" "), mean(r), var(r), var(r)/mean(r), min(r), max(r), D, n-1, p_low)) }
cat("\n== Differenze appaiate estate - autunno (farnia) ==\n")
a <- m[m$G=="Quercus_robur Fall",]; e <- m[m$G=="Quercus_robur Summer",]
dd <- e$R[match(a$N_pianta, e$N_pianta)] - a$R; names(dd) <- a$N_pianta; print(dd)
cat("media", mean(dd), "sd", round(sd(dd),2), "valori distinti", length(unique(dd)), "\n")
cat("\n== Variabili binarie per taxon: struttura delle tabelle 2x2 ==\n")
classif <- function(x, y){ # x, y = n piante positive su 10 nei due gruppi
  if (x+y==0) "assente in entrambi" else if ((x==0 & y==10)|(x==10 & y==0)) "separazione completa" else if (x==0|y==0|x==10|y==10) "separazione quasi-completa" else "nessuna separazione" }
cmp <- list(Stagione=c("Quercus_robur Fall","Quercus_robur Summer"), Ospite=c("Quercus_cerris Summer","Quercus_robur Summer"))
for (k in names(cmp)) { A <- colSums(pa[m$G==cmp[[k]][1],]); B <- colSums(pa[m$G==cmp[[k]][2],])
  cl <- mapply(classif, A, B); cat("\n", k, ":\n"); print(table(cl))
  cat("  taxa con <3 piante positive in totale:", sum(A+B>0 & A+B<3), "\n") }
