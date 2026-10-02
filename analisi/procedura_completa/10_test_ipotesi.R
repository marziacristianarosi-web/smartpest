# =============================================================================
# 10_test_ipotesi.R - Passaggio 10: test delle sei ipotesi (un test per domanda, scelto a priori)
# Quadro unico: inferenza per permutazione; lo schema di permutazione segue il disegno.
#   Stagione (stesse 10 farnie): scambio delle etichette solo entro pianta, tutte le 2^10 = 1024 inversioni (p senza errore Monte Carlo,
#   esatto rispetto allo schema di permutazione sotto scambiabilità entro pianta)
#   Ospite (piante diverse):     riassegnazione libera delle etichette, 9999 permutazioni (p = (b+1)/(m+1))
# D1/D4 ricchezza:      LRT del modello di Poisson generalizzata (glmmTMB), p di permutazione
# D2/D5 composizione:   PERMANOVA su distanze di Jaccard (vegan::adonis2) + PERMDISP (betadisper, permutest)
# D3/D6 singoli taxa:   McNemar esatto (binom.test sulle coppie discordanti) / Fisher esatto (fisher.test),
#                       correzione di Benjamini-Hochberg (1995) sui taxa osservati nel confronto
# Versione per piastra delle domande 4-6 (estate, 4 piastre per pianta): le unità restano le piante, quindi si
# permutano piante intere (con tutte le loro piastre): D4 GLMM con (1|pianta); D5 PERMANOVA sulle piastre;
# D6 numero di piastre positive per pianta (0-4).
# Analisi di sensibilità per la composizione: test W*d robusto a dispersioni diverse (Hamidi et al. 2019;
# pacchetto WdStar, installazione: remotes::install_github("alekseyenko/WdStar")).
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
set.seed(20261002)
NPERM <- 9999
RICALCOLA <- as.logical(Sys.getenv("CORY_RICALCOLA", "TRUE"))   # FALSE = riusa le permutazioni salvate
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
# distribuzione già calcolata nel passaggio 8 (se presente viene riusata)
L1 <- if (file.exists(fileL1)) readRDS(fileL1) else unlist(in_parallelo(seq_len(nrow(PERM_ENTRO)), function(i, S, P) {
  dd <- S; dd$g <- S$g[P[i, ]]; lrt_gp(R ~ g + (1 | pl), R ~ 1 + (1 | pl), dd) }, S = S, P = PERM_ENTRO))
r1 <- rapporto(f1)
D1 <- data.frame(Domanda = "D1 ricchezza ~ stagione", Media_gruppo1 = mean(S$R[S$g == "Summer"]), Media_gruppo2 = mean(S$R[S$g == "Fall"]),
                 Rapporto_medie = r1[1], IC95_inf = r1[2], IC95_sup = r1[3], LRT = L1_oss, Permutazioni = sum(!is.na(L1)),
                 p = mean(L1[!is.na(L1)] >= L1_oss - 1e-8))

# ---- D4: ricchezza ~ ospite (GLM, permutazione libera) ----
f4 <- glmmTMB(R ~ g, data = H, family = genpois())
L4_oss <- lrt_gp(R ~ g, R ~ 1, H)
fileL4 <- file.path(DIR_TAB, "D4_distribuzione_permutazione.rds")
L4 <- if (!RICALCOLA && file.exists(fileL4)) readRDS(fileL4) else unlist(in_parallelo(seq_len(NPERM), function(i, H, P) { dd <- H; dd$g <- H$g[P[i, ]]; lrt_gp(R ~ g, R ~ 1, dd) },
                          H = H, P = P_LIBERE))
saveRDS(L4, fileL4)
r4 <- rapporto(f4)
D4 <- data.frame(Domanda = "D4 ricchezza ~ ospite", Media_gruppo1 = mean(H$R[H$g == "Quercus_robur"]), Media_gruppo2 = mean(H$R[H$g == "Quercus_cerris"]),
                 Rapporto_medie = r4[1], IC95_inf = r4[2], IC95_sup = r4[3], LRT = L4_oss, Permutazioni = sum(!is.na(L4)),
                 p = (sum(L4[!is.na(L4)] >= L4_oss - 1e-8) + 1) / (sum(!is.na(L4)) + 1))
# ---- D4 per piastra: GLMM R ~ ospite + (1|pianta), dispersione diversa per ospite (A.7: A4 non rispettata) ----
f4p <- glmmTMB(R ~ g + (1 | pl), dispformula = ~ g, data = E, family = genpois())
lrt_p <- function(dd) { a <- try(suppressWarnings(glmmTMB(R ~ g + (1 | pl), dispformula = ~ g, data = dd, family = genpois())), silent = TRUE)
  b <- try(suppressWarnings(glmmTMB(R ~ 1 + (1 | pl), dispformula = ~ g, data = dd, family = genpois())), silent = TRUE)
  if (inherits(a, "try-error") || inherits(b, "try-error")) NA else as.numeric(2 * (logLik(a) - logLik(b))) }
L4p_oss <- lrt_p(E)
fileL4p <- file.path(DIR_TAB, "D4p_distribuzione_permutazione.rds")
L4p <- if (!RICALCOLA && file.exists(fileL4p)) readRDS(fileL4p) else unlist(in_parallelo(seq_len(NPERM), function(i, E, P, perm_piante, lrt_p) {
  dd <- E; dd$g <- perm_piante(P[i, ]); lrt_p(dd) }, E = E, P = P_LIBERE, perm_piante = perm_piante, lrt_p = lrt_p, esporta = c("H", "E")))
saveRDS(L4p, fileL4p)
r4p <- rapporto(f4p)
D4pr <- data.frame(Domanda = "D4 ricchezza per piastra ~ ospite", Media_gruppo1 = mean(E$R[E$g == "Quercus_robur"]),
                   Media_gruppo2 = mean(E$R[E$g == "Quercus_cerris"]), Rapporto_medie = r4p[1], IC95_inf = r4p[2], IC95_sup = r4p[3],
                   LRT = L4p_oss, Permutazioni = sum(!is.na(L4p)), p = (sum(L4p[!is.na(L4p)] >= L4p_oss - 1e-8) + 1) / (sum(!is.na(L4p)) + 1))
ric <- rbind(D1, D4, D4pr); ric[, -1] <- round(ric[, -1], 4)
print(ric, row.names = FALSE); salva_tab(ric, "T16_ricchezza_D1_D4")
capture.output(summary(f1), summary(f4), summary(f4p), file = file.path(DIR_LOG, "modelli_ricchezza.txt"))

# ---- Verifica di robustezza (dichiarata, non un secondo test): differenza tra le medie come statistica ----
# Stesse permutazioni dell'LRT; la validità del test di permutazione non dipende dalla statistica (Good 2005).
dm <- function(y, g, a, b) mean(y[g == a]) - mean(y[g == b])
Dd1 <- apply(PERM_ENTRO, 1, function(idx) dm(S$R, S$g[idx], "Summer", "Fall")); D1o <- dm(S$R, S$g, "Summer", "Fall")
Dd4 <- apply(P_LIBERE, 1, function(idx) dm(H$R, H$g[idx], "Quercus_cerris", "Quercus_robur")); D4o <- dm(H$R, H$g, "Quercus_cerris", "Quercus_robur")
rob <- data.frame(Domanda = c("D1", "D4"), Differenza_medie = c(D1o, D4o),
  p_differenza_medie = c(mean(abs(Dd1) >= abs(D1o) - 1e-9), (sum(abs(Dd4) >= abs(D4o) - 1e-9) + 1) / (NPERM + 1)),
  p_LRT = c(D1$p, D4$p),
  Correlazione_rango_LRT_vs_diff = c(cor(L1, abs(Dd1), method = "spearman", use = "complete.obs"),
                                     cor(L4, abs(Dd4), method = "spearman", use = "complete.obs")))
rob[, -1] <- round(rob[, -1], 4); print(rob, row.names = FALSE); salva_tab(rob, "T20_robustezza_statistica")

# ---- D2 e D5: composizione (PERMANOVA + PERMDISP) ----
dS <- vegdist(PA_S, method = "jaccard", binary = TRUE); dH <- vegdist(PA_H, method = "jaccard", binary = TRUE)
P2 <- PERM_ENTRO[-1, ]                                   # senza l'identità: vegan aggiunge l'osservato
a2 <- adonis2(dS ~ pl + g, data = S, permutations = P2, by = "terms")
b2 <- betadisper(dS, S$g); pd2 <- permutest(b2, permutations = P2)
a5 <- adonis2(dH ~ g, data = H, permutations = NPERM)
b5 <- betadisper(dH, H$g); pd5 <- permutest(b5, permutations = NPERM)
# D5 per piastra: 80 piastre, permutazione di piante intere (ogni pianta porta con sé le sue 4 piastre)
dE <- vegdist(PA_E, method = "jaccard", binary = TRUE)
CTRL_P <- how(plots = Plots(strata = E$pl, type = "free"), within = Within(type = "none"), nperm = NPERM)
a5p <- adonis2(dE ~ g, data = E, permutations = CTRL_P)
b5p <- betadisper(dE, E$g); pd5p <- permutest(b5p, permutations = CTRL_P)
comp <- data.frame(Domanda = c("D2 composizione ~ stagione", "D5 composizione ~ ospite", "D5 composizione per piastra ~ ospite"),
  Pseudo_F = c(a2["g", "F"], a5["g", "F"], a5p["g", "F"]), R2 = c(a2["g", "R2"], a5["g", "R2"], a5p["g", "R2"]),
  p_PERMANOVA = c(a2["g", "Pr(>F)"], a5["g", "Pr(>F)"], a5p["g", "Pr(>F)"]),
  Dispersione_gruppo1 = c(mean(b2$distances[S$g == "Summer"]), mean(b5$distances[H$g == "Quercus_robur"]), mean(b5p$distances[E$g == "Quercus_robur"])),
  Dispersione_gruppo2 = c(mean(b2$distances[S$g == "Fall"]), mean(b5$distances[H$g == "Quercus_cerris"]), mean(b5p$distances[E$g == "Quercus_cerris"])),
  F_PERMDISP = c(pd2$tab[1, "F"], pd5$tab[1, "F"], pd5p$tab[1, "F"]),
  p_PERMDISP = c(pd2$tab[1, "Pr(>F)"], pd5$tab[1, "Pr(>F)"], pd5p$tab[1, "Pr(>F)"]))
# Sensibilità: test W*d (statistica Tw2 di Welch per due gruppi, robusta a dispersioni diverse; Hamidi et al. 2019)
# con lo stesso schema di permutazione: tutte le 1024 permutazioni entro pianta (D2), 9999 libere (D5)
Tw2 <- WdStar::Tw2
T2_oss <- Tw2(dS, S$g); T2 <- apply(PERM_ENTRO, 1, function(idx) Tw2(dS, S$g[idx]))
T5_oss <- Tw2(dH, H$g); T5 <- apply(P_LIBERE, 1, function(idx) Tw2(dH, H$g[idx]))
# piastre separate: 80 piastre, permutazioni di piante intere (stesse 9999 riassegnazioni delle piante)
T5p_oss <- Tw2(dE, E$g); T5p <- apply(P_LIBERE, 1, function(idx) Tw2(dE, perm_piante(idx)))
wd <- data.frame(Domanda = c("D2 stagione (1024 permutazioni entro pianta)", "D5 ospite (9999 permutazioni libere)",
                             "D5 ospite per piastra (9999 permutazioni di piante intere)"),
                 Tw2 = round(c(T2_oss, T5_oss, T5p_oss), 3),
                 p = round(c(mean(T2 >= T2_oss - 1e-9), (sum(T5 >= T5_oss - 1e-9) + 1) / (NPERM + 1),
                             (sum(T5p >= T5p_oss - 1e-9) + 1) / (NPERM + 1)), 4))
print(wd, row.names = FALSE); salva_tab(wd, "T24_Wd_sensibilita")
comp[, -1] <- round(comp[, -1], 4); print(comp, row.names = FALSE); salva_tab(comp, "T17_composizione_D2_D5")
capture.output(a2, pd2, a5, pd5, a5p, pd5p, file = file.path(DIR_LOG, "permanova_permdisp.txt"))
# Tabelle complete (formato ANOVA). Nel modello stagionale il termine pianta serve solo a rappresentare
# l'appaiamento: con permutazioni entro pianta non è saggiabile, quindi il suo p-value non è riportato.
tab_pm <- function(a, dom, fattore) { x <- as.data.frame(a)
  fonte <- c(pl = "Pianta", g = fattore, Residual = "Residuo", Total = "Totale")[rownames(x)]
  data.frame(Domanda = dom, Fonte = unname(fonte), gl = x$Df, Somma_quadrati = round(x$SumOfSqs, 4),
             Media_quadrati = round(x$SumOfSqs / x$Df, 4), Pseudo_F = round(x$F, 3), R2 = round(x$R2, 3),
             p = round(x$`Pr(>F)`, 4)) }
pm <- rbind(tab_pm(a2, "D2 stagione (permutazioni entro pianta, 1024)", "Stagione"),
            tab_pm(a5, "D5 ospite (permutazioni libere, 9999)", "Ospite"),
            tab_pm(a5p, "D5 ospite per piastra (permutazioni di piante intere, 9999)", "Ospite"))
pm$p[pm$Fonte == "Pianta"] <- NA; pm$Media_quadrati[pm$Fonte == "Totale"] <- NA
print(pm, row.names = FALSE); salva_tab(pm, "T21_PERMANOVA")
tab_pd <- function(pt, b, g, dom) { x <- as.data.frame(pt$tab); m <- tapply(b$distances, g, mean)
  data.frame(Domanda = dom, Fonte = c("Gruppi", "Residuo"), gl = x$Df, Somma_quadrati = round(x$`Sum Sq`, 4),
             Media_quadrati = round(x$`Mean Sq`, 4), F = round(x$F, 3), Permutazioni = x$N.Perm, p = round(x$`Pr(>F)`, 4),
             Distanza_media_dal_centroide = c(paste(names(m), sprintf("%.3f", m), collapse = "; "), NA)) }
pdd <- rbind(tab_pd(pd2, b2, S$g, "D2 stagione"), tab_pd(pd5, b5, H$g, "D5 ospite"), tab_pd(pd5p, b5p, E$g, "D5 ospite per piastra"))
print(pdd, row.names = FALSE); salva_tab(pdd, "T22_PERMDISP")

# ---- D3 e D6: singoli taxa ----
# Correzione di Benjamini-Hochberg sui taxa osservati nel confronto (presenti in almeno una pianta)
bh <- function(p, osservato) { q <- rep(NA_real_, length(p)); q[osservato] <- p.adjust(p[osservato], "BH"); q }
fa <- PA_S[S$g == "Fall", ]; su <- PA_S[S$g == "Summer", ]
stopifnot(identical(S$pl[S$g == "Fall"], S$pl[S$g == "Summer"]))
t3 <- data.frame(Taxon = TAXA, Estate = colSums(su), Autunno = colSums(fa),
                 b_solo_autunno = colSums(fa == 1 & su == 0), c_solo_estate = colSums(fa == 0 & su == 1))
nd <- t3$b_solo_autunno + t3$c_solo_estate
t3$p <- ifelse(nd == 0, 1, mapply(function(b, n) binom.test(b, n)$p.value, t3$b_solo_autunno, pmax(nd, 1)))
ci <- t(mapply(function(b, n) if (n == 0) c(NA, NA) else binom.test(b, n)$conf.int, t3$b_solo_autunno, pmax(nd, 1)))
t3$OR_condizionato <- ifelse(nd == 0, NA, t3$b_solo_autunno / t3$c_solo_estate)
t3$OR_inf <- ci[, 1] / (1 - ci[, 1]); t3$OR_sup <- ci[, 2] / (1 - ci[, 2])
t3$p_BH <- bh(t3$p, t3$Estate + t3$Autunno > 0)
print(t3, digits = 3, row.names = FALSE); salva_tab(t3, "T18_taxa_stagione_D3")

ro <- PA_H[H$g == "Quercus_robur", ]; ce <- PA_H[H$g == "Quercus_cerris", ]; n1 <- nrow(ro); n2 <- nrow(ce)
t6 <- data.frame(Taxon = TAXA, Q_robur = colSums(ro), Q_cerris = colSums(ce))
fis <- lapply(seq_along(TAXA), function(j) fisher.test(matrix(c(t6$Q_cerris[j], n2 - t6$Q_cerris[j], t6$Q_robur[j], n1 - t6$Q_robur[j]), 2)))
t6$p <- sapply(fis, `[[`, "p.value")
t6$OR_cerro_su_farnia <- sapply(fis, function(x) unname(x$estimate))
t6$OR_inf <- sapply(fis, function(x) x$conf.int[1]); t6$OR_sup <- sapply(fis, function(x) x$conf.int[2])
deg <- (t6$Q_robur + t6$Q_cerris) %in% c(0, n1 + n2); t6[deg, c("OR_cerro_su_farnia", "OR_inf", "OR_sup")] <- NA
t6$p_BH <- bh(t6$p, t6$Q_robur + t6$Q_cerris > 0); t6$p[t6$Q_robur + t6$Q_cerris == 0] <- NA
print(t6, digits = 3, row.names = FALSE); salva_tab(t6, "T19_taxa_ospite_D6")

# D6 per piastra: numero di piastre positive per pianta (0-4); statistica = differenza tra le medie dei due ospiti;
# p da 9999 permutazioni libere delle piante (le stesse della domanda 4); correzione di Benjamini-Hochberg
fr <- FREQ_H
d6 <- function(g) colMeans(fr[g == "Quercus_cerris", , drop = FALSE]) - colMeans(fr[g == "Quercus_robur", , drop = FALSE])
D6o <- d6(H$g); D6perm <- apply(P_LIBERE, 1, function(idx) d6(H$g[idx]))
t6p <- data.frame(Taxon = TAXA, Media_piastre_Q_robur = round(colMeans(fr[H$g == "Quercus_robur", ]), 2),
                  Media_piastre_Q_cerris = round(colMeans(fr[H$g == "Quercus_cerris", ]), 2), Differenza_cerro_meno_farnia = round(D6o, 2))
t6p$p <- sapply(seq_along(TAXA), function(j) (sum(abs(D6perm[j, ]) >= abs(D6o[j]) - 1e-9) + 1) / (NPERM + 1))
oss6 <- colSums(fr) > 0; t6p$p[!oss6] <- NA; t6p$p_BH <- bh(t6p$p, oss6)
print(t6p, row.names = FALSE); salva_tab(t6p, "T25_taxa_ospite_per_piastra")

# Misura dell'effetto complementare per D6 a piastre separate: modello binomiale misto sulle 80 piastre
# (presenza nella piastra ~ ospite + (1 | pianta)); odds ratio cerro/farnia con IC dal profilo di verosimiglianza.
# Solo stima descrittiva: il p-value resta quello di permutazione (scelto a priori). IC non corretti per la molteplicità.
or_bin <- do.call(rbind, lapply(TAXA[oss6], function(t) {
  d <- data.frame(y = PA_E[, t], g = E$g, pl = E$pl)
  m <- try(suppressWarnings(glmmTMB(y ~ g + (1 | pl), family = binomial, data = d)), silent = TRUE)
  ok <- !inherits(m, "try-error") && m$fit$convergence == 0 && isTRUE(m$sdr$pdHess)
  ci <- if (ok) try(suppressWarnings(confint(m, parm = "gQuercus_cerris", method = "profile")), silent = TRUE) else NULL
  se <- if (ok) summary(m)$coefficients$cond[2, 2] else NA
  stimabile <- ok && se < 10 && !inherits(ci, "try-error")   # se enorme = separazione (taxon in una sola piastra)
  data.frame(Taxon = t, Piastre_pos_robur = sum(d$y[d$g == "Quercus_robur"]), Piastre_pos_cerris = sum(d$y[d$g == "Quercus_cerris"]),
             OR_cerro_su_farnia = if (stimabile) round(exp(fixef(m)$cond[2]), 3) else NA,
             IC95_inf = if (stimabile) round(exp(ci[1, 1]), 3) else NA, IC95_sup = if (stimabile) round(exp(ci[1, 2]), 3) else NA,
             Varianza_pianta = if (ok) signif(VarCorr(m)$cond$pl[1], 3) else NA,
             Nota = if (stimabile) "" else "non stimabile (separazione)")
}))
print(or_bin, row.names = FALSE); salva_tab(or_bin, "T27_OR_binomiale_misto_piastre")

# ---- Figura: distribuzioni nulle di permutazione e statistica osservata (D1, D4) ----
nulla <- function(L, oss, tit, sotto) {
  ggplot(data.frame(L = L[!is.na(L)]), aes(L)) + geom_histogram(bins = 40, fill = "grey75", colour = "white") +
    geom_vline(xintercept = oss, colour = "#c0392b", linewidth = 0.8) +
    labs(x = "Statistica LRT", y = "Permutazioni", title = tit, subtitle = sotto) + TEMA }
pN <- nulla(L1, L1_oss, "D1: 1024 permutazioni entro pianta", sprintf("Rosso: valore osservato\n(LRT = %s; p = %s)", virgola(L1_oss, 1), virgola(D1$p, 3))) +
  nulla(L4, L4_oss, "D4: 9999 permutazioni libere", sprintf("Rosso: valore osservato\n(LRT = %s; p = %s)", virgola(L4_oss, 3), virgola(D4$p, 2))) +
  plot_annotation(tag_levels = "A")
salva_fig(pN, "F08_distribuzioni_permutazione", 170, 75)
# pannelli singoli, uno per domanda (usati nella Parte B del report)
salva_fig(nulla(L1, L1_oss, "D1: 1024 permutazioni entro pianta", sprintf("Rosso: valore osservato (LRT = %s; p = %s)", virgola(L1_oss, 1), virgola(D1$p, 3))),
          "F08A_permutazioni_D1", 100, 75)
salva_fig(nulla(L4, L4_oss, "D4: 9999 permutazioni libere", sprintf("Rosso: valore osservato (LRT = %s; p = %s)", virgola(L4_oss, 3), virgola(D4$p, 2))),
          "F08B_permutazioni_D4", 100, 75)
fp4 <- function(p) if (p < 0.001) "p < 0,001" else paste("p =", virgola(p, 3))
salva_fig(nulla(L4p, L4p_oss, "D4 per piastra: 9999 permutazioni di piante", sprintf("Rosso: valore osservato (LRT = %s; %s)", virgola(L4p_oss, 2), fp4(D4pr$p))),
          "F08C_permutazioni_D4_piastra", 100, 75)
