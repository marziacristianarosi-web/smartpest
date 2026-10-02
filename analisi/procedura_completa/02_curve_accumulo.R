# =============================================================================
# 02_curve_accumulo.R - Passaggio 2: curve di accumulo dei taxa
# (a) tra piante: vegan::specaccum(method = "exact") = rarefazione per campioni (Mao Tau)
# (b) entro pianta: accumulo sulle 4 piastre delle piante estive (sforzo entro pianta)
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")

# (a) tra piante -------------------------------------------------------------
acc <- do.call(rbind, lapply(levels(PIANTE$gruppo), function(g) {
  sa <- specaccum(PA[PIANTE$gruppo == g, ], method = "exact")
  data.frame(gruppo = g, piante = sa$sites, taxa = sa$richness, ds = sa$sd)
}))
acc$gruppo <- factor(acc$gruppo, levels = levels(PIANTE$gruppo))
acc$inf <- pmax(0, acc$taxa - 1.96 * acc$ds); acc$sup <- acc$taxa + 1.96 * acc$ds
riep <- do.call(rbind, lapply(split(acc, acc$gruppo), function(x) data.frame(
  Gruppo = ETI[as.character(x$gruppo[1])], Taxa_10_piante = x$taxa[10], Taxa_5_piante = round(x$taxa[5], 1),
  Percentuale_a_5 = round(100 * x$taxa[5] / x$taxa[10]), Incremento_ultima_pianta = round(x$taxa[10] - x$taxa[9], 2))))
print(riep); salva_tab(riep, "T03_curve_accumulo")

pA <- ggplot(acc, aes(piante, taxa, colour = gruppo, fill = gruppo)) +
  geom_ribbon(aes(ymin = inf, ymax = sup), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(aes(shape = gruppo), size = 1.6) +
  scale_colour_manual(values = COL, labels = ETI_IT) + scale_fill_manual(values = COL, labels = ETI_IT) +
  scale_shape_manual(values = SHP, labels = ETI_IT) +
  scale_x_continuous(breaks = 1:10) + scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2), expand = c(0, 0)) +
  labs(x = "Numero di piante", y = "Numero cumulato di taxa", colour = NULL, fill = NULL, shape = NULL,
       title = "Accumulo tra piante") + TEMA + theme(legend.position = c(0.68, 0.2), legend.text = element_text(size = 7.5))

# (b) entro pianta: 4 piastre da 3 insetti (solo estate; in autunno le piastre sono registrate insieme)
est <- grezzi[grezzi$Season == "Summer", ]
entro <- do.call(rbind, lapply(split(est, paste(est$Host_species, est$N_pianta)), function(x) {
  m <- as.matrix(x[TAXA]); tot <- sum(colSums(m) > 0)
  sa <- specaccum(m, method = "exact")
  data.frame(ospite = x$Host_species[1], pianta = x$N_pianta[1], piastre = sa$sites, frazione = sa$richness / tot)
}))
entro$gruppo <- ifelse(entro$ospite == "Quercus_robur", "farnia_estate", "cerro_estate")
me <- aggregate(frazione ~ gruppo + piastre, entro, function(v) c(media = mean(v), min = min(v), max = max(v)))
me <- data.frame(me[1:2], me$frazione)
print(me); salva_tab(me, "T04_completezza_entro_pianta")
pB <- ggplot(me, aes(piastre, media, colour = gruppo)) +
  geom_linerange(aes(ymin = min, ymax = max), position = position_dodge(0.25), linewidth = 0.5) +
  geom_line(position = position_dodge(0.25), linewidth = 0.7) +
  geom_point(aes(shape = gruppo), position = position_dodge(0.25), size = 1.8) +
  scale_colour_manual(values = COL, labels = ETI_IT) + scale_shape_manual(values = SHP, labels = ETI_IT) +
  scale_y_continuous(limits = c(0, 1.02), labels = function(x) paste0(100 * x, "%"), expand = c(0, 0)) +
  labs(x = "Numero di piastre (3 insetti ciascuna)", y = "Taxa della pianta rilevati",
       title = "Accumulo entro pianta (estate)", colour = NULL, shape = NULL) +
  TEMA + theme(legend.position = c(0.65, 0.2), legend.text = element_text(size = 7.5))
salva_fig(pA + pB + plot_annotation(tag_levels = "A"), "F01_curve_accumulo", 170, 80)
