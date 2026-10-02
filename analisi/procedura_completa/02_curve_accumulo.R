# =============================================================================
# 02_curve_accumulo.R - A.2: curve di accumulo dei taxa (rarefazione per campioni)
# Funzione: iNEXT::iNEXT(q = 0, datatype = "incidence_freq") tra piante; entro pianta la stessa formula
# analitica (iNEXT non stima le curve su 4 sole piastre quando molti taxa sono in tutte). La parte interpolata
# coincide con la curva di accumulo per campioni (Mao Tau; Colwell et al. 2012), con IC 95% bootstrap.
# (a) tra piante: 1-10 piante per gruppo; (b) entro pianta: 1-4 piastre delle piante estive.
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
set.seed(2026)

# (a) tra piante -------------------------------------------------------------
inc <- lapply(levels(PIANTE$gruppo), function(g) { m <- PA[PIANTE$gruppo == g, ]; c(nrow(m), colSums(m)) })
names(inc) <- levels(PIANTE$gruppo)
acc <- iNEXT(inc, q = 0, datatype = "incidence_freq", size = 1:10, nboot = 200)$iNextEst$size_based
acc <- acc[acc$t <= 10, ]; acc$gruppo <- factor(acc$Assemblage, levels = levels(PIANTE$gruppo))
riep <- do.call(rbind, lapply(split(acc, acc$gruppo), function(x) data.frame(
  Gruppo = ETI[as.character(x$gruppo[1])], Taxa_10_piante = x$qD[x$t == 10],
  IC95_10 = sprintf("%.1f-%.1f", x$qD.LCL[x$t == 10], x$qD.UCL[x$t == 10]),
  Taxa_5_piante = round(x$qD[x$t == 5], 1), Percentuale_a_5 = round(100 * x$qD[x$t == 5] / x$qD[x$t == 10]),
  Incremento_ultima_pianta = round(x$qD[x$t == 10] - x$qD[x$t == 9], 2))))
print(riep); salva_tab(riep, "T03_curve_accumulo")

pA <- ggplot(acc, aes(t, qD, colour = gruppo, fill = gruppo)) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.15, colour = NA) +
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
  # stessa formula di rarefazione per incidenza usata da iNEXT (Chao et al. 2014):
  # E[S(t)] = somma sui taxa di 1 - C(T - Y_k, t) / C(T, t), con T = 4 piastre e Y_k = piastre con il taxon k
  y <- colSums(m); y <- y[y > 0]; T <- nrow(m)
  St <- sapply(1:T, function(t) sum(1 - choose(T - y, t) / choose(T, t)))
  data.frame(ospite = x$Host_species[1], pianta = x$N_pianta[1], piastre = 1:T, frazione = St / tot)
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
