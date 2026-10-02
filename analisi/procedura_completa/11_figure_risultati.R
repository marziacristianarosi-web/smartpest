# =============================================================================
# 11_figure_risultati.R - Passaggio 11: figure dei risultati
# F09 ricchezza per pianta; F10 ordinamento PCoA (Jaccard, vegan::wcmdscale); F11 incidenza dei taxa.
# I valori riportati nelle figure sono letti dalle tabelle prodotte da 10_test_ipotesi.R.
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
ric  <- as.data.frame(read_excel(file.path(DIR_TAB, "T16_ricchezza_D1_D4.xlsx")))
comp <- as.data.frame(read_excel(file.path(DIR_TAB, "T17_composizione_D2_D5.xlsx")))
t3   <- as.data.frame(read_excel(file.path(DIR_TAB, "T18_taxa_stagione_D3.xlsx")))
fp <- function(p) if (p < 0.001) "p < 0,001" else paste("p =", virgola(p, if (p < 0.01) 3 else 2))

# ---- F09: ricchezza ----
media_ic <- function(x) { t <- t.test(x); data.frame(y = mean(x), ymin = t$conf.int[1], ymax = t$conf.int[2]) }
Sf <- S; Sf$x <- as.numeric(Sf$g)                       # 1 = estate, 2 = autunno
set.seed(1); sp <- setNames(runif(nlevels(Sf$pl), -0.07, 0.07), levels(Sf$pl)); Sf$xj <- Sf$x + sp[as.character(Sf$pl)]
yl <- c(0, max(PIANTE$R) + 2)
p1a <- ggplot(Sf, aes(x, R)) +
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
  geom_point(aes(colour = gruppo, shape = gruppo), size = 1.8, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  stat_summary(fun.data = media_ic, geom = "errorbar", width = 0, linewidth = 0.6, position = position_nudge(x = 0.25)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6, fill = "white", position = position_nudge(x = 0.25)) +
  annotate("text", x = 1.5, y = yl[2] - 0.5, label = fp(ric$p[2]), size = 2.9, family = FONT) +
  scale_x_discrete(labels = c(farnia_estate = expression(italic("Q. robur")), cerro_estate = expression(italic("Q. cerris")))) +
  scale_colour_manual(values = COL) + scale_shape_manual(values = SHP) +
  scale_y_continuous(limits = yl, breaks = seq(0, 12, 2), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Estate – specie ospite") + TEMA
salva_fig(p1a + p1b + plot_annotation(tag_levels = "A"), "F09_ricchezza_risultati", 140, 75)

# ---- F10: PCoA su distanze di Jaccard ----
pcoa <- function(dd, M) {
  o <- wcmdscale(vegdist(M, "jaccard", binary = TRUE), eig = TRUE)
  list(sc = cbind(dd, A1 = o$points[, 1], A2 = o$points[, 2]), pv = 100 * o$eig[1:2] / sum(o$eig[o$eig > 0]))
}
grafico_pcoa <- function(o, titolo, testo, coppie = FALSE) {
  sc <- o$sc; cen <- aggregate(cbind(A1, A2) ~ gruppo, sc, mean); p <- ggplot(sc, aes(A1, A2))
  if (coppie) {
    w <- reshape(sc[, c("N_pianta", "gruppo", "A1", "A2")], idvar = "N_pianta", timevar = "gruppo", direction = "wide")
    p <- p + geom_segment(data = w, aes(x = A1.farnia_estate, y = A2.farnia_estate, xend = A1.farnia_autunno, yend = A2.farnia_autunno),
                          colour = "grey80", linewidth = 0.3) }
  p + geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.3) +
    geom_point(aes(colour = gruppo, shape = gruppo), size = 2, alpha = 0.9, position = position_jitter(width = 0.008, height = 0.008, seed = 2)) +
    geom_point(data = cen, aes(fill = gruppo), shape = 21, size = 4, colour = "white", stroke = 0.8) +
    scale_colour_manual(values = COL, labels = ETI_IT, name = NULL, limits = names(COL), drop = FALSE) +
    scale_shape_manual(values = SHP, labels = ETI_IT, name = NULL, limits = names(COL), drop = FALSE) +
    scale_fill_manual(values = COL, guide = "none", limits = names(COL)) + coord_equal() +
    labs(x = sprintf("PCoA 1 (%.0f%%)", o$pv[1]), y = sprintf("PCoA 2 (%.0f%%)", o$pv[2]), title = titolo, subtitle = testo) +
    TEMA + theme(legend.position = "bottom", legend.text = element_text(size = 8))
}
txt <- function(i) sprintf("PERMANOVA: R² = %s; %s\nPERMDISP: %s", virgola(comp$R2[i]), fp(comp$p_PERMANOVA[i]), fp(comp$p_PERMDISP[i]))
p2a <- grafico_pcoa(pcoa(S, PA_S), expression(italic("Q. robur")*" – stagione"), txt(1), coppie = TRUE)
p2b <- grafico_pcoa(pcoa(H, PA_H), "Estate – specie ospite", txt(2))
salva_fig((p2a + p2b + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + plot_annotation(tag_levels = "A"),
          "F10_PCoA_composizione", 174, 100)

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
