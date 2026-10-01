# Script di verifica (sessione di analisi Corythucha). Il file dati "dati.xlsx" va sostituito con il percorso dell Excel.
suppressMessages({library(readxl); library(glmmTMB); library(lme4); library(permutes)}); set.seed(1)
d <- as.data.frame(read_excel("dati.xlsx")); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
m <- d[!duplicated(key), c("Season","Host_species","N_pianta")]; rownames(m) <- unique(key); m <- m[rownames(pa),]; m$R <- rowSums(pa)
S <- m[m$Host_species=="Quercus_robur",]; S$g <- factor(S$Season); S$pl <- factor(S$N_pianta)
# Ricostruisco la risposta "perm_y" che permutes costruisce internamente, per vederne la natura
for (k in c("glmer_poisson","glmmTMB_genpois")) {
  mod <- if (k=="glmer_poisson") glmer(R ~ g + (1|pl), data=S, family=poisson) else glmmTMB(R ~ g + (1|pl), data=S, family=genpois())
  X <- getME(mod,"X"); B <- if (inherits(mod,"glmmTMB")) fixef(mod)$cond else fixef(mod)
  keep <- colnames(X)=="gSummer"; X[,!keep] <- 0
  e <- as.vector(resid(mod) + X %*% B); e <- exp(e)
  cat(sprintf("\n%s: tipo di residuo usato da resid() = %s\n", k, if (inherits(mod,"glmmTMB")) "response (y - media)" else "deviance"))
  cat("  risposta originale (prime 6):", head(S$R), "\n  risposta trasformata che permutes analizza (prime 6):", round(head(e),3), "\n")
}
w <- character(0)
res <- withCallingHandlers(perm.glmmTMB(R ~ g + (1|pl), data=S, family=genpois(), nperm=200, progress=FALSE),
        warning=function(x){ w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")})
cat("\nAvvisi emessi da perm.glmmTMB (unici):\n"); print(unique(substr(w,1,120)))
cat("Colonne restituite:", names(res), "\n")
