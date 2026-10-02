# =============================================================================
# 10_test_ipotesi.R - Passaggio 10: test delle sei ipotesi (un test per domanda, scelto a priori)
# Quadro unico: inferenza per permutazione; lo schema di permutazione segue il disegno.
#   Stagione (stesse 10 farnie): scambio delle etichette solo entro pianta, 2^10 = 1024 permutazioni (p esatto)
#   Ospite (piante diverse):     riassegnazione libera delle etichette, 9999 permutazioni (p = (b+1)/(m+1))
# D1/D4 ricchezza:      LRT del modello di Poisson generalizzata (glmmTMB), p di permutazione
# D2/D5 composizione:   PERMANOVA su distanze di Jaccard (vegan::adonis2) + PERMDISP (betadisper, permutest)
# D3/D6 singoli taxa:   McNemar esatto (binom.test sulle coppie discordanti) / Fisher esatto (fisher.test),
#                       esclusione delle ipotesi non verificabili (Tarone 1990) e FDR di Benjamini-Hochberg (Gilbert 2005)
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
set.seed(20261002)
NPERM <- 9999
P_LIBERE <- t(replicate(NPERM, sample(nrow(H))))
rapporto <- function(fit) {   # rapporto tra le medie con IC 95% dal profilo di verosimiglianza
  ci <- try(confint(fit, parm = 2, method = "profile"), silent = TRUE)
  if (inherits(ci, "try-error")) ci <- confint(fit, parm = 2, method = "wald")
  exp(c(stima = unname(fixef(fit)$cond[2]), inf = ci[1, 1], sup = ci[1, 2]))
}

# ---- D1: ricchezza ~ stagione (GLMM, permutazione entro pianta) ----
f1 <- glmmTMB(R ~ g + (1 | pl), data = S, family = genpois())
L1_oss <- lrt_gp(R ~ g + (1 | pl), R ~ 1 + (1 | pl), S)
fileL1 <- file.path(DIR_TAB, "D1_distribuzione_permutazione.rds")
L1 <- if (file.exists(fileL1)) readRDS(fileL1) else unlist(in_parallelo(seq_len(nrow(PERM_ENTRO)), function(i, S, P) {
  dd <- S; dd$g <- S$g[P[i, ]]; lrt_gp(R ~ g + (1 | pl), R ~ 1 + (1 | pl), dd) }, S = S, P = PERM_ENTRO))
r1 <- rapporto(f1)
D1 <- data.frame(Domanda = "D1 ricchezza ~ stagione", Media_gruppo1 = mean(S$R[S$g == "Summer"]), Media_gruppo2 = mean(S$R[S$g == "Fall"]),
                 Rapporto_medie = r1[1], IC95_inf = r1[2], IC95_sup = r1[3], LRT = L1_oss, Permutazioni = sum(!is.na(L1)),
                 p = mean(L1[!is.na(L1)] >= L1_oss - 1e-8))

# ---- D4: ricchezza ~ ospite (GLM, permutazione libera) ----
f4 <- glmmTMB(R ~ g, data = H, family = genpois())
L4_oss <- lrt_gp(R ~ g, R ~ 1, H)
L4 <- unlist(in_parallelo(seq_len(NPERM), function(i, H, P) { dd <- H; dd$g <- H$g[P[i, ]]; lrt_gp(R ~ g, R ~ 1, dd) },
                          H = H, P = P_LIBERE))
saveRDS(L4, file.path(DIR_TAB, "D4_distribuzione_permutazione.rds"))
r4 <- rapporto(f4)
D4 <- data.frame(Domanda = "D4 ricchezza ~ ospite", Media_gruppo1 = mean(H$R[H$g == "Quercus_robur"]), Media_gruppo2 = mean(H$R[H$g == "Quercus_cerris"]),
                 Rapporto_medie = r4[1], IC95_inf = r4[2], IC95_sup = r4[3], LRT = L4_oss, Permutazioni = sum(!is.na(L4)),
                 p = (sum(L4[!is.na(L4)] >= L4_oss - 1e-8) + 1) / (sum(!is.na(L4)) + 1))
ric <- rbind(D1, D4); ric[, -1] <- round(ric[, -1], 4)
print(ric, row.names = FALSE); salva_tab(ric, "T16_ricchezza_D1_D4")
capture.output(summary(f1), summary(f4), file = file.path(DIR_LOG, "modelli_ricchezza.txt"))

# ---- D2 e D5: composizione (PERMANOVA + PERMDISP) ----
dS <- vegdist(PA_S, method = "jaccard", binary = TRUE); dH <- vegdist(PA_H, method = "jaccard", binary = TRUE)
P2 <- PERM_ENTRO[-1, ]                                   # senza l'identità: vegan aggiunge l'osservato
a2 <- adonis2(dS ~ pl + g, data = S, permutations = P2, by = "terms")
b2 <- betadisper(dS, S$g); pd2 <- permutest(b2, permutations = P2)
a5 <- adonis2(dH ~ g, data = H, permutations = NPERM)
b5 <- betadisper(dH, H$g); pd5 <- permutest(b5, permutations = NPERM)
comp <- data.frame(Domanda = c("D2 composizione ~ stagione", "D5 composizione ~ ospite"),
  Pseudo_F = c(a2["g", "F"], a5["g", "F"]), R2 = c(a2["g", "R2"], a5["g", "R2"]), p_PERMANOVA = c(a2["g", "Pr(>F)"], a5["g", "Pr(>F)"]),
  Dispersione_gruppo1 = c(mean(b2$distances[S$g == "Summer"]), mean(b5$distances[H$g == "Quercus_robur"])),
  Dispersione_gruppo2 = c(mean(b2$distances[S$g == "Fall"]), mean(b5$distances[H$g == "Quercus_cerris"])),
  F_PERMDISP = c(pd2$tab[1, "F"], pd5$tab[1, "F"]), p_PERMDISP = c(pd2$tab[1, "Pr(>F)"], pd5$tab[1, "Pr(>F)"]))
comp[, -1] <- round(comp[, -1], 4); print(comp, row.names = FALSE); salva_tab(comp, "T17_composizione_D2_D5")
capture.output(a2, pd2, a5, pd5, file = file.path(DIR_LOG, "permanova_permdisp.txt"))

# ---- D3 e D6: singoli taxa ----
tarone_bh <- function(p, pmin, alpha = 0.05) {
  m <- length(p); K <- 1
  while (K <= m && sum(pmin <= alpha / K) > K) K <- K + 1
  R <- pmin <= alpha / K; q <- rep(NA_real_, m); q[R] <- p.adjust(p[R], "BH")
  list(K = K, testabili = R, q = q)
}
fa <- PA_S[S$g == "Fall", ]; su <- PA_S[S$g == "Summer", ]
stopifnot(identical(S$pl[S$g == "Fall"], S$pl[S$g == "Summer"]))
t3 <- data.frame(Taxon = TAXA, Estate = colSums(su), Autunno = colSums(fa),
                 b_solo_autunno = colSums(fa == 1 & su == 0), c_solo_estate = colSums(fa == 0 & su == 1))
nd <- t3$b_solo_autunno + t3$c_solo_estate
t3$p <- ifelse(nd == 0, 1, mapply(function(b, n) binom.test(b, n)$p.value, t3$b_solo_autunno, pmax(nd, 1)))
t3$p_minimo <- pmin(1, 2 * 0.5^nd)
ci <- t(mapply(function(b, n) if (n == 0) c(NA, NA) else binom.test(b, n)$conf.int, t3$b_solo_autunno, pmax(nd, 1)))
t3$OR_condizionato <- ifelse(nd == 0, NA, t3$b_solo_autunno / t3$c_solo_estate)
t3$OR_inf <- ci[, 1] / (1 - ci[, 1]); t3$OR_sup <- ci[, 2] / (1 - ci[, 2])
tb3 <- tarone_bh(t3$p, t3$p_minimo); t3$Saggiabile <- tb3$testabili; t3$q_BH <- tb3$q
cat(sprintf("D3: ipotesi verificabili (Tarone, K = %d): %d su %d\n", tb3$K, sum(tb3$testabili), length(TAXA)))
print(t3, digits = 3, row.names = FALSE); salva_tab(t3, "T18_taxa_stagione_D3")

ro <- PA_H[H$g == "Quercus_robur", ]; ce <- PA_H[H$g == "Quercus_cerris", ]; n1 <- nrow(ro); n2 <- nrow(ce)
t6 <- data.frame(Taxon = TAXA, Q_robur = colSums(ro), Q_cerris = colSums(ce))
fis <- lapply(seq_along(TAXA), function(j) fisher.test(matrix(c(t6$Q_cerris[j], n2 - t6$Q_cerris[j], t6$Q_robur[j], n1 - t6$Q_robur[j]), 2)))
t6$p <- sapply(fis, `[[`, "p.value")
t6$OR_cerro_su_farnia <- sapply(fis, function(x) unname(x$estimate))
t6$OR_inf <- sapply(fis, function(x) x$conf.int[1]); t6$OR_sup <- sapply(fis, function(x) x$conf.int[2])
t6$p_minimo <- sapply(t6$Q_robur + t6$Q_cerris, function(k) {   # tabella più estrema con lo stesso margine
  if (k == 0 || k == n1 + n2) return(1)
  a <- min(k, n2); fisher.test(matrix(c(a, n2 - a, k - a, n1 - (k - a)), 2))$p.value })
deg <- (t6$Q_robur + t6$Q_cerris) %in% c(0, n1 + n2); t6[deg, c("OR_cerro_su_farnia", "OR_inf", "OR_sup")] <- NA
tb6 <- tarone_bh(t6$p, t6$p_minimo); t6$Saggiabile <- tb6$testabili; t6$q_BH <- tb6$q
cat(sprintf("D6: ipotesi verificabili (Tarone, K = %d): %d su %d\n", tb6$K, sum(tb6$testabili), length(TAXA)))
print(t6, digits = 3, row.names = FALSE); salva_tab(t6, "T19_taxa_ospite_D6")

# ---- Figura: distribuzioni nulle di permutazione e statistica osservata (D1, D4) ----
nulla <- function(L, oss, tit, sotto) {
  ggplot(data.frame(L = L[!is.na(L)]), aes(L)) + geom_histogram(bins = 40, fill = "grey75", colour = "white") +
    geom_vline(xintercept = oss, colour = "#c0392b", linewidth = 0.8) +
    labs(x = "Statistica LRT", y = "Permutazioni", title = tit, subtitle = sotto) + TEMA }
pN <- nulla(L1, L1_oss, "Domanda 1: 1024 permutazioni entro pianta", sprintf("Rosso: valore osservato (LRT = %s; p = %s)", virgola(L1_oss, 1), virgola(D1$p, 3))) +
  nulla(L4, L4_oss, "Domanda 4: 9999 permutazioni libere", sprintf("Rosso: valore osservato (LRT = %s; p = %s)", virgola(L4_oss, 3), virgola(D4$p, 2))) +
  plot_annotation(tag_levels = "A")
salva_fig(pN, "F08_distribuzioni_permutazione", 170, 75)
