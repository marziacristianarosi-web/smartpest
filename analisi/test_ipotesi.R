# Test delle sei ipotesi (Corythucha arcuata - comunità fungina).
# Uso: Rscript test_ipotesi.R <file_excel> <cartella_output> [domande, es. 1,4,2,5,3,6]
# I dati e i risultati NON vanno caricati nel repository.
#
# Quadro inferenziale unico: tutti i p-value sono di permutazione (o test esatti condizionati,
# che ne sono il caso particolare per dati binari). Lo schema di permutazione segue il disegno:
#   stagione (stesse 10 farnie)  -> scambio delle etichette solo entro pianta (2^10 = 1024, enumerazione completa)
#   ospite (piante diverse)      -> riassegnazione libera delle etichette fra piante
suppressMessages({library(readxl); library(glmmTMB); library(vegan)})
args <- commandArgs(trailingOnly = TRUE)
f_dati <- if (length(args) >= 1) args[1] else "dati.xlsx"
out    <- if (length(args) >= 2) args[2] else "risultati_test"
domande <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else as.character(1:6)
ncore  <- max(1, parallel::detectCores() - 0)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
set.seed(20261002)
NPERM_LIBERE <- 9999
alpha <- 0.05

# ---- Dati a livello di pianta (unità campionaria; 12 insetti per pianta in entrambe le stagioni) ----
d  <- as.data.frame(read_excel(f_dati))
tx <- setdiff(names(d), c("Season", "Host_species", "N_pianta", "Replicate", "N_insetti"))
tx <- tx[nzchar(tx) & !grepl("^\\.\\.\\.", tx)]
key <- paste(d$Season, d$Host_species, d$N_pianta, sep = "|")
pa  <- (rowsum(as.matrix(d[tx]), key, reorder = FALSE) > 0) * 1
m   <- d[!duplicated(key), c("Season", "Host_species", "N_pianta")]
rownames(m) <- unique(key); m <- m[rownames(pa), ]
m$R <- rowSums(pa)
stopifnot(all(m$R > 0))   # Jaccard non definito per unità vuote

# Domande 1-3: Q. robur, autunno vs estate, appaiato per pianta
iS <- which(m$Host_species == "Quercus_robur")
iS <- iS[order(m$N_pianta[iS], m$Season[iS])]
S  <- m[iS, ]; S$pl <- factor(S$N_pianta); S$g <- factor(S$Season, levels = c("Summer", "Fall"))
paS <- pa[iS, ]
piante <- levels(S$pl)
stopifnot(all(table(S$pl) == 2))
# Tutte le 1024 permutazioni entro pianta (scambio o no delle due etichette di ciascuna pianta)
combo <- as.matrix(expand.grid(rep(list(0:1), length(piante))))
perm_entro <- t(apply(combo, 1, function(sw) {
  idx <- seq_len(nrow(S))
  for (i in which(sw == 1)) { k <- which(S$pl == piante[i]); idx[k] <- rev(idx[k]) }
  idx }))

# Domande 4-6: estate, Q. robur vs Q. cerris, piante indipendenti
iH <- which(m$Season == "Summer")
H  <- m[iH, ]; H$g <- factor(H$Host_species, levels = c("Quercus_robur", "Quercus_cerris"))
paH <- pa[iH, ]
perm_libere <- t(replicate(NPERM_LIBERE, sample(nrow(H))))

lrt <- function(f1, f0, dd) {
  a <- try(suppressWarnings(glmmTMB(f1, data = dd, family = genpois())), silent = TRUE)
  b <- try(suppressWarnings(glmmTMB(f0, data = dd, family = genpois())), silent = TRUE)
  if (inherits(a, "try-error") || inherits(b, "try-error")) return(NA_real_)
  as.numeric(2 * (logLik(a) - logLik(b)))
}
rapporto <- function(fit) {          # rapporto tra medie (exp del coefficiente) con IC 95% da profilo
  ci <- try(confint(fit, parm = 2, method = "profile"), silent = TRUE)
  if (inherits(ci, "try-error")) ci <- confint(fit, parm = 2, method = "wald")
  exp(c(stima = unname(fixef(fit)$cond[2]), inf = ci[1, 1], sup = ci[1, 2]))
}

# Matrici di dissimilarità di Jaccard (domande 2 e 5)
dS <- vegdist(paS, method = "jaccard", binary = TRUE)
dH <- vegdist(paH, method = "jaccard", binary = TRUE)

if ("1" %in% domande) {
sink(file.path(out, "D1.txt"), split = TRUE)
# ======================= DOMANDA 1: ricchezza ~ stagione ==========================
cat("\n=== D1. Ricchezza, stagione (Q. robur, n = 10 coppie) ===\n")
print(aggregate(R ~ g, S, function(x) c(media = mean(x), ds = sd(x), min = min(x), max = max(x))))
f1 <- glmmTMB(R ~ g + (1 | pl), data = S, family = genpois())
print(summary(f1))
L_obs <- lrt(R ~ g + (1 | pl), R ~ 1 + (1 | pl), S)
L_perm <- unlist(parallel::mclapply(seq_len(nrow(perm_entro)), mc.cores = ncore, function(i) { idx <- perm_entro[i, ]; dd <- S; dd$g <- S$g[idx]; lrt(R ~ g + (1 | pl), R ~ 1 + (1 | pl), dd) }))
ok <- !is.na(L_perm)
cat(sprintf("LRT osservato = %.3f | permutazioni valide = %d/1024 | p esatto = %.4f\n",
            L_obs, sum(ok), mean(L_perm[ok] >= L_obs - 1e-8)))
cat("Rapporto tra medie (autunno/estate) e IC 95% (profilo):\n"); print(round(rapporto(f1), 3))
saveRDS(list(L = L_perm, L_obs = L_obs), file.path(out, "D1_permutazioni.rds"))
sink()
}

if ("4" %in% domande) {
sink(file.path(out, "D4.txt"), split = TRUE)
# ======================= DOMANDA 4: ricchezza ~ ospite ============================
cat("\n=== D4. Ricchezza, ospite (estate, 10 + 10 piante) ===\n")
print(aggregate(R ~ g, H, function(x) c(media = mean(x), ds = sd(x), min = min(x), max = max(x))))
f4 <- glmmTMB(R ~ g, data = H, family = genpois())
print(summary(f4))
L_obs4 <- lrt(R ~ g, R ~ 1, H)
L_perm4 <- unlist(parallel::mclapply(seq_len(nrow(perm_libere)), mc.cores = ncore, function(i) { idx <- perm_libere[i, ]; dd <- H; dd$g <- H$g[idx]; lrt(R ~ g, R ~ 1, dd) }))
ok4 <- !is.na(L_perm4)
cat(sprintf("LRT osservato = %.3f | permutazioni valide = %d/%d | p = %.4f\n", L_obs4, sum(ok4), NPERM_LIBERE,
            (sum(L_perm4[ok4] >= L_obs4 - 1e-8) + 1) / (sum(ok4) + 1)))
cat("Rapporto tra medie (cerro/farnia) e IC 95% (profilo):\n"); print(round(rapporto(f4), 3))
saveRDS(list(L = L_perm4, L_obs = L_obs4), file.path(out, "D4_permutazioni.rds"))
sink()
}

if ("2" %in% domande) {
sink(file.path(out, "D2.txt"), split = TRUE)
cat("\n=== D2. Composizione, stagione: PERMANOVA (Jaccard), permutazioni entro pianta ===\n")
P2 <- perm_entro[-1, ]                       # esclusa l'identità: vegan aggiunge l'osservato
print(adonis2(dS ~ pl + g, data = S, permutations = P2, by = "terms"))
bd2 <- betadisper(dS, S$g)
cat("PERMDISP (verifica dell'omogeneità della dispersione), permutazioni entro pianta:\n")
print(tapply(bd2$distances, S$g, mean)); print(permutest(bd2, permutations = P2))
sink()
}

if ("5" %in% domande) {
sink(file.path(out, "D5.txt"), split = TRUE)
cat("\n=== D5. Composizione, ospite: PERMANOVA (Jaccard), permutazioni libere ===\n")
print(adonis2(dH ~ g, data = H, permutations = NPERM_LIBERE))
bd5 <- betadisper(dH, H$g)
cat("PERMDISP:\n"); print(tapply(bd5$distances, H$g, mean)); print(permutest(bd5, permutations = NPERM_LIBERE))
sink()
}

# ======================= DOMANDE 3 e 6: singoli taxa ==============================
# Esclusione delle ipotesi non verificabili (Tarone 1990) e FDR di Benjamini-Hochberg sulle restanti (Gilbert 2005)
tarone_bh <- function(p, pmin, alpha = 0.05) {
  m <- length(p); K <- 1
  while (K <= m && sum(pmin <= alpha / K) > K) K <- K + 1
  R <- pmin <= alpha / K
  q <- rep(NA_real_, m); q[R] <- p.adjust(p[R], "BH")
  list(K = K, testabili = R, q = q)
}

if ("3" %in% domande) {
sink(file.path(out, "D3.txt"), split = TRUE)
cat("\n=== D3. Singoli taxa, stagione: McNemar esatto (binomiale sulle coppie discordanti) ===\n")
fa <- paS[S$g == "Fall", ]; su <- paS[S$g == "Summer", ]          # righe allineate per pianta
stopifnot(identical(S$pl[S$g == "Fall"], S$pl[S$g == "Summer"]))
t3 <- data.frame(taxon = tx, prev_estate = colSums(su), prev_autunno = colSums(fa),
                 solo_autunno = colSums(fa == 1 & su == 0), solo_estate = colSums(fa == 0 & su == 1))
nd <- t3$solo_autunno + t3$solo_estate
t3$p <- ifelse(nd == 0, 1, mapply(function(b, n) binom.test(b, n)$p.value, t3$solo_autunno, pmax(nd, 1)))
t3$p_min <- pmin(1, 2 * 0.5^nd)
ci <- t(mapply(function(b, n) if (n == 0) c(NA, NA) else binom.test(b, n)$conf.int, t3$solo_autunno, pmax(nd, 1)))
t3$OR_cond <- t3$solo_autunno / t3$solo_estate
t3$OR_inf <- ci[, 1] / (1 - ci[, 1]); t3$OR_sup <- ci[, 2] / (1 - ci[, 2])
tb3 <- tarone_bh(t3$p, t3$p_min); t3$testabile <- tb3$testabili; t3$q_BH <- tb3$q
cat(sprintf("Ipotesi verificabili (Tarone, K = %d): %d su %d\n", tb3$K, sum(tb3$testabili), length(tx)))
print(t3, digits = 3, row.names = FALSE)
write.csv(t3, file.path(out, "D3_taxa_stagione.csv"), row.names = FALSE)
sink()
}

if ("6" %in% domande) {
sink(file.path(out, "D6.txt"), split = TRUE)
cat("\n=== D6. Singoli taxa, ospite: test esatto di Fisher ===\n")
ro <- paH[H$g == "Quercus_robur", ]; ce <- paH[H$g == "Quercus_cerris", ]
n1 <- nrow(ro); n2 <- nrow(ce)
t6 <- data.frame(taxon = tx, prev_farnia = colSums(ro), prev_cerro = colSums(ce))
fis <- lapply(seq_along(tx), function(j) fisher.test(matrix(c(t6$prev_cerro[j], n2 - t6$prev_cerro[j],
                                                              t6$prev_farnia[j], n1 - t6$prev_farnia[j]), 2)))
t6$p <- sapply(fis, `[[`, "p.value")
t6$OR_cerro_farnia <- sapply(fis, function(x) unname(x$estimate))
t6$OR_inf <- sapply(fis, function(x) x$conf.int[1]); t6$OR_sup <- sapply(fis, function(x) x$conf.int[2])
t6$p_min <- sapply(t6$prev_farnia + t6$prev_cerro, function(k) {   # tabella più estrema con lo stesso margine
  if (k == 0 || k == n1 + n2) return(1)
  a <- min(k, n2); fisher.test(matrix(c(a, n2 - a, k - a, n1 - (k - a)), 2))$p.value })
nd6 <- (t6$prev_farnia + t6$prev_cerro) %in% c(0, n1 + n2)   # tabella degenere: OR non stimabile
t6[nd6, c("OR_cerro_farnia", "OR_inf", "OR_sup")] <- NA
tb6 <- tarone_bh(t6$p, t6$p_min); t6$testabile <- tb6$testabili; t6$q_BH <- tb6$q
cat(sprintf("Ipotesi verificabili (Tarone, K = %d): %d su %d\n", tb6$K, sum(tb6$testabili), length(tx)))
print(t6, digits = 3, row.names = FALSE)
write.csv(t6, file.path(out, "D6_taxa_ospite.csv"), row.names = FALSE)
sink()
}

