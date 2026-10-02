# =============================================================================
# 11_figure_risultati.R - Passaggio 11: figure dei risultati
# F09 ricchezza per pianta; F10 ordinamento PCoA (Jaccard, vegan::wcmdscale); F11 incidenza dei taxa.
# I valori riportati nelle figure sono letti dalle tabelle prodotte da 10_test_ipotesi.R.
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
ric  <- as.data.frame(read_excel(file.path(DIR_TAB, "T16_ricchezza_D1_D4.xlsx")))   # righe: D1, D4 per pianta, D4 per piastra
comp <- as.data.frame(read_excel(file.path(DIR_TAB, "T17_composizione_D2_D5.xlsx")))
t3   <- as.data.frame(read_excel(file.path(DIR_TAB, "T18_taxa_stagione_D3.xlsx")))
fp <- function(p) if (p < 0.001) "p < 0,001" else paste("p =", virgola(p, if (p < 0.01) 3 else 2))

# ---- F09: ricchezza ----
media_ic <- function(x) { t <- t.test(x); data.frame(y = mean(x), ymin = t$conf.int[1], ymax = t$conf.int[2]) }
Sf <- S; Sf$x <- as.numeric(Sf$g)                       # 1 = estate, 2 = autunno
set.seed(1); sp <- setNames(runif(nlevels(Sf$pl), -0.07, 0.07), levels(Sf$pl)); Sf$xj <- Sf$x + sp[as.character(Sf$pl)]
yl <- c(0, max(PIANTE$R) + 2)
p1a <- ggplot(Sf, aes(x, R)) +
  geom_boxplot(aes(group = g), width = 0.32, fill = NA, colour = "grey55", outlier.shape = NA, linewidth = 0.35) +
  geom_line(aes(xj, group = pl), colour = "grey75", linewidth = 0.35) +
  geom_point(aes(xj, colour = gruppo, shape = gruppo), size = 1.8) +
  stat_summary(aes(group = g), fun.data = media_ic, geom = "errorbar", width = 0, linewidth = 0.6, position = position_nudge(x = 0.25)) +
  stat_summary(aes(group = g), fun = mean, geom = "point", shape = 23, size = 2.6, fill = "white", position = position_nudge(x = 0.25)) +
  annotate("text", x = 1.5, y = yl[2] - 0.5, label = fp(ric$p[1]), size = 2.9, family = FONT) +
  scale_x_continuous(breaks = 1:2, labels = c("Estate", "Autunno"), limits = c(0.6, 2.4)) +
  scale_colour_manual(values = COL) + scale_shape_manual(values = SHP) +
  scale_y_continuous(limits = yl, breaks = seq(0, 12, 2), expand = c(0, 0)) +
  labs(x = NULL, y = "Ricchezza in taxa per pianta", title = expression(italic("Q. robur")*" – stagione")) + TEMA
p1b <- ggplot(H, aes(gruppo, R)) +
  geom_boxplot(width = 0.32, fill = NA, colour = "grey55", outlier.shape = NA, linewidth = 0.35) +
  geom_point(aes(colour = gruppo, shape = gruppo), size = 1.8, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  stat_summary(fun.data = media_ic, geom = "errorbar", width = 0, linewidth = 0.6, position = position_nudge(x = 0.25)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6, fill = "white", position = position_nudge(x = 0.25)) +
  annotate("text", x = 1.5, y = yl[2] - 0.5, label = fp(ric$p[2]), size = 2.9, family = FONT) +
  scale_x_discrete(labels = c(farnia_estate = expression(italic("Q. robur")), cerro_estate = expression(italic("Q. cerris")))) +
  scale_colour_manual(values = COL) + scale_shape_manual(values = SHP) +
  scale_y_continuous(limits = yl, breaks = seq(0, 12, 2), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Estate – specie ospite") + TEMA
salva_fig(p1a + p1b + plot_annotation(tag_levels = "A"), "F09_ricchezza_risultati", 140, 75)
salva_fig(p1a, "F09A_ricchezza_stagione", 80, 75); salva_fig(p1b + labs(y = "Ricchezza in taxa per pianta"), "F09B_ricchezza_ospite", 80, 75)

# ---- F10: PCoA su distanze di Jaccard ----
pcoa <- function(dd, M) {
  o <- wcmdscale(vegdist(M, "jaccard", binary = TRUE), eig = TRUE)
  list(sc = cbind(dd, A1 = o$points[, 1], A2 = o$points[, 2]), pv = 100 * o$eig[1:2] / sum(o$eig[o$eig > 0]))
}
grafico_pcoa <- function(o, titolo, testo, coppie = FALSE, assi = NULL) {
  sc <- o$sc; cen <- aggregate(cbind(A1, A2) ~ gruppo, sc, mean); p <- ggplot(sc, aes(A1, A2))
  if (coppie) {
    w <- reshape(sc[, c("N_pianta", "gruppo", "A1", "A2")], idvar = "N_pianta", timevar = "gruppo", direction = "wide")
    p <- p + geom_segment(data = w, aes(x = A1.farnia_estate, y = A2.farnia_estate, xend = A1.farnia_autunno, yend = A2.farnia_autunno),
                          colour = "grey80", linewidth = 0.3) }
  p + geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.3) +
    # ellisse al 95% per gruppo (distribuzione t multivariata, ggplot2::stat_ellipse): solo descrittiva
    stat_ellipse(aes(colour = gruppo), type = "t", level = 0.95, linewidth = 0.5, show.legend = FALSE) +
    geom_point(aes(colour = gruppo, shape = gruppo), size = 2, alpha = 0.9, position = position_jitter(width = 0.008, height = 0.008, seed = 2)) +
    geom_point(data = cen, aes(fill = gruppo), shape = 21, size = 4, colour = "white", stroke = 0.8) +
    scale_colour_manual(values = COL, labels = ETI_IT, name = NULL, limits = names(COL), drop = FALSE) +
    scale_shape_manual(values = SHP, labels = ETI_IT, name = NULL, limits = names(COL), drop = FALSE) +
    scale_fill_manual(values = COL, guide = "none", limits = names(COL)) + coord_equal() +
    labs(x = if (is.null(assi)) sprintf("PCoA 1 (%.0f%%)", o$pv[1]) else assi[1],
         y = if (is.null(assi)) sprintf("PCoA 2 (%.0f%%)", o$pv[2]) else assi[2], title = titolo, subtitle = testo) +
    TEMA + theme(legend.position = "bottom", legend.text = element_text(size = 8))
}
txt <- function(i) sprintf("PERMANOVA: R² = %s; %s\nPERMDISP: %s", virgola(comp$R2[i]), fp(comp$p_PERMANOVA[i]), fp(comp$p_PERMDISP[i]))
p2a <- grafico_pcoa(pcoa(S, PA_S), expression(italic("Q. robur")*" – stagione"), txt(1))
p2b <- grafico_pcoa(pcoa(H, PA_H), "Estate – specie ospite", txt(2))
salva_fig((p2a + p2b + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + plot_annotation(tag_levels = "A"),
          "F10_PCoA_composizione", 174, 100)
# nei pannelli singoli la legenda mostra solo i gruppi presenti
solo <- function(p, liv) suppressMessages(p + scale_colour_manual(values = COL, labels = ETI_IT[liv], limits = liv, name = NULL) +
  scale_shape_manual(values = SHP, labels = ETI_IT[liv], limits = liv, name = NULL))
LS <- c("farnia_estate", "farnia_autunno"); LH <- c("farnia_estate", "cerro_estate")
salva_fig(solo(p2a, LS), "F10A_PCoA_stagione", 100, 105); salva_fig(solo(p2b, LH), "F10B_PCoA_ospite", 100, 105)

# ---- F12: NMDS su distanze di Jaccard (vegan::metaMDS), rappresentazione complementare alla PCoA ----
# L'NMDS conserva solo l'ordine delle dissimilarità; lo stress indica la qualità della rappresentazione
# (valori bassi = distanze rappresentate fedelmente). È una visualizzazione: il test resta la PERMANOVA.
nmds <- function(dd, M) {
  set.seed(2026)
  o <- metaMDS(vegdist(M, "jaccard", binary = TRUE), k = 2, trymax = 200, trace = FALSE, autotransform = FALSE)
  list(sc = cbind(dd, A1 = o$points[, 1], A2 = o$points[, 2]), stress = o$stress, conv = o$converged)
}
nS <- nmds(S, PA_S); nH <- nmds(H, PA_H)
st <- data.frame(Confronto = c("Stagione (Q. robur)", "Ospite (estate)"), Stress = round(c(nS$stress, nH$stress), 3),
                 Convergenza = c(nS$conv, nH$conv) > 0)   # in vegan 2.6 un codice > 0 indica convergenza
print(st); salva_tab(st, "T23_NMDS_stress")
txtN <- function(o, i) sprintf("Stress = %s\nPERMANOVA: R² = %s; %s", virgola(o$stress, 3), virgola(comp$R2[i]), fp(comp$p_PERMANOVA[i]))
p4a <- grafico_pcoa(nS, expression(italic("Q. robur")*" \u2013 stagione"), txtN(nS, 1), assi = c("NMDS 1", "NMDS 2"))
p4b <- grafico_pcoa(nH, "Estate \u2013 specie ospite", txtN(nH, 2), assi = c("NMDS 1", "NMDS 2"))
salva_fig((p4a + p4b + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + plot_annotation(tag_levels = "A"),
          "F12_NMDS_composizione", 174, 100)
salva_fig(solo(p4a, LS), "F12A_NMDS_stagione", 100, 105); salva_fig(solo(p4b, LH), "F12B_NMDS_ospite", 100, 105)

# ---- Piastre separate (estate, domande 4-6): ricchezza per piastra, ordinamenti delle 80 piastre ----
Ep <- E; Ep$gruppo <- factor(ifelse(E$g == "Quercus_robur", "farnia_estate", "cerro_estate"), levels = levels(PIANTE$gruppo))
p1c <- ggplot(Ep, aes(gruppo, R)) +
  geom_boxplot(width = 0.32, fill = NA, colour = "grey55", outlier.shape = NA, linewidth = 0.35) +
  geom_point(aes(colour = gruppo, shape = gruppo), size = 1.5, alpha = 0.8, position = position_jitter(width = 0.1, height = 0.12, seed = 1)) +
  stat_summary(fun.data = media_ic, geom = "errorbar", width = 0, linewidth = 0.6, position = position_nudge(x = 0.3)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6, fill = "white", position = position_nudge(x = 0.3)) +
  annotate("text", x = 1.5, y = 10.5, label = fp(ric$p[3]), size = 2.9, family = FONT) +
  scale_x_discrete(labels = c(farnia_estate = expression(italic("Q. robur")), cerro_estate = expression(italic("Q. cerris")))) +
  scale_colour_manual(values = COL) + scale_shape_manual(values = SHP) +
  scale_y_continuous(limits = c(0, 11), breaks = seq(0, 10, 2), expand = c(0, 0)) +
  labs(x = NULL, y = "Ricchezza in taxa per piastra (3 insetti)", title = "Estate \u2013 ospite, piastre separate",
       subtitle = "Punti = 80 piastre (4 per pianta)\nrombo = media con IC 95% (descrittivo)") + TEMA
salva_fig(p1c, "F09C_ricchezza_piastre", 90, 80)
p2c <- solo(grafico_pcoa(pcoa(Ep, PA_E), "Estate \u2013 ospite, piastre separate", txt(3)), LH)
salva_fig(p2c, "F10C_PCoA_piastre", 100, 105)
nE <- nmds(Ep, PA_E)
p4c <- solo(grafico_pcoa(nE, "Estate \u2013 ospite, piastre separate", sprintf("Stress = %s\nPERMANOVA: R² = %s; %s", virgola(nE$stress, 3),
            virgola(comp$R2[3]), fp(comp$p_PERMANOVA[3])), assi = c("NMDS 1", "NMDS 2")), LH)
salva_fig(p4c, "F12C_NMDS_piastre", 100, 105)
st <- rbind(st, data.frame(Confronto = "Ospite, piastre separate", Stress = round(nE$stress, 3), Convergenza = nE$conv > 0))
salva_tab(st, "T23_NMDS_stress")

# ---- F11: incidenza dei taxa ----
inc <- aggregate(PA, list(gr = PIANTE$gruppo), sum)
lg <- reshape(inc, direction = "long", varying = TAXA, v.names = "n", timevar = "taxon", times = TAXA, idvar = "gr")
nm <- function(x) { x <- sub("Thrichoderma", "Trichoderma", gsub("_", " ", x)); x <- sub(" [Ss]pp(\\d)$", " sp. \\1", x); sub(" [Ss]pp$", " sp.", x) }
dif <- unlist(inc[inc$gr == "farnia_estate", TAXA]) - unlist(inc[inc$gr == "farnia_autunno", TAXA])
lg$taxon <- factor(nm(lg$taxon), levels = nm(TAXA[order(dif, colSums(PA[, TAXA]))]))
sig <- t3$Taxon[!is.na(t3$q_BH) & t3$q_BH < 0.05]       # taxa con q < 0,05 (domanda 3)
lg$etichetta <- ifelse(lg$n == 0, "", lg$n)
p3 <- ggplot(lg, aes(gr, taxon)) + geom_tile(aes(fill = n), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = etichetta, colour = n >= 6), size = 2.7, family = FONT) +
  scale_fill_gradient(low = "#f1f4f9", high = "#1d4e8f", limits = c(0, 10), breaks = c(0, 5, 10), name = "Piante con\nil taxon (su 10)") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_x_discrete(labels = c(farnia_estate = expression(atop(italic("Q. robur"), "estate")),
                              farnia_autunno = expression(atop(italic("Q. robur"), "autunno")),
                              cerro_estate = expression(atop(italic("Q. cerris"), "estate"))), position = "top") +
  scale_y_discrete(labels = function(x) parse(text = sapply(x, function(t) {
    e <- if (grepl(" sp\\.", t)) sprintf('italic("%s")*" %s"', sub(" sp\\..*", "", t), sub("^\\S+ ", "", t)) else sprintf('italic("%s")', t)
    if (t %in% nm(sig)) paste0(e, '*" *"') else e }))) +
  labs(x = NULL, y = NULL) + TEMA +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title = element_text(size = 8))
salva_fig(p3, "F11_incidenza_taxa", 120, 110)

# ---- F11B: frequenza dei taxa con piastre separate (numero medio di piastre positive per pianta, 0-4) ----
t6p <- as.data.frame(read_excel(file.path(DIR_TAB, "T25_taxa_ospite_per_piastra.xlsx")))
fq <- rbind(data.frame(taxon = t6p$Taxon, gr = "farnia_estate", n = t6p$Media_piastre_Q_robur),
            data.frame(taxon = t6p$Taxon, gr = "cerro_estate", n = t6p$Media_piastre_Q_cerris))
sig6p <- t6p$Taxon[!is.na(t6p$q_BH) & t6p$q_BH < 0.05]   # taxa con q < 0,05 (domanda 6, piastre separate)
fq$gr <- factor(fq$gr, levels = c("farnia_estate", "cerro_estate")); fq$taxon <- factor(nm(fq$taxon), levels = levels(lg$taxon))
p3b <- ggplot(fq, aes(gr, taxon)) + geom_tile(aes(fill = n), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(n == 0, "", formatC(n, format = "f", digits = 1, decimal.mark = ",")), colour = n >= 2.4), size = 2.7, family = FONT) +
  scale_fill_gradient(low = "#f1f4f9", high = "#1d4e8f", limits = c(0, 4), breaks = 0:4, name = "Piastre positive\nper pianta (0-4)") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_x_discrete(labels = c(farnia_estate = expression(atop(italic("Q. robur"), "estate")),
                              cerro_estate = expression(atop(italic("Q. cerris"), "estate"))), position = "top") +
  scale_y_discrete(labels = function(x) parse(text = sapply(x, function(t) {
    e <- if (grepl(" sp\\.", t)) sprintf('italic("%s")*" %s"', sub(" sp\\..*", "", t), sub("^\\S+ ", "", t)) else sprintf('italic("%s")', t)
    if (t %in% nm(sig6p)) paste0(e, '*" *"') else e }))) +
  labs(x = NULL, y = NULL) + TEMA +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title = element_text(size = 8))
salva_fig(p3b, "F11B_frequenza_taxa_piastre", 100, 110)
