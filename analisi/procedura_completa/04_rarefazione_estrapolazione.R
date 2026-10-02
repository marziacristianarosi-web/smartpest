# =============================================================================
# 04_rarefazione_estrapolazione.R - Passaggio 4: rarefazione ed estrapolazione (numeri di Hill)
# Funzioni: iNEXT::iNEXT(q = 0, 1, 2; datatype = "incidence_freq"; endpoint = 20 = 2T)
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
set.seed(2026)
inc <- lapply(levels(PIANTE$gruppo), function(g) { m <- PA[PIANTE$gruppo == g, ]; c(nrow(m), colSums(m)) })
names(inc) <- levels(PIANTE$gruppo)
out <- iNEXT(inc, q = c(0, 1, 2), datatype = "incidence_freq", endpoint = 20, nboot = 200)
asy <- out$AsyEst; asy$Assemblage <- ETI[asy$Assemblage]
print(asy); salva_tab(asy, "T06_hill_asintotici")
s <- out$iNextEst$size_based
tab <- s[s$t %in% c(5, 10, 20), c("Assemblage", "t", "Method", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC")]
tab$Assemblage <- ETI[tab$Assemblage]; salva_tab(tab, "T07_rarefazione_estrapolazione")

s$gruppo <- factor(s$Assemblage, levels = names(COL)); s0 <- s[s$Order.q == 0, ]
oss <- s0[s0$Method == "Observed", ]
zona <- annotate("rect", xmin = 10, xmax = 20, ymin = -Inf, ymax = Inf, fill = "grey93")
p1 <- ggplot(s0, aes(t, qD, colour = gruppo, fill = gruppo)) + zona +
  annotate("text", x = c(5.5, 15), y = 11.4, label = c("rarefazione (dati)", "estrapolazione"), size = 2.6,
           colour = "grey35", family = FONT) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.15, colour = NA) +
  geom_line(data = s0[s0$t <= 10, ], linewidth = 0.7) + geom_line(data = s0[s0$t >= 10, ], linewidth = 0.7, linetype = "22") +
  geom_point(data = oss, aes(shape = gruppo), size = 2) +
  scale_colour_manual(values = COL, labels = ETI_IT) + scale_fill_manual(values = COL, labels = ETI_IT) +
  scale_shape_manual(values = SHP, labels = ETI_IT) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) + scale_y_continuous(limits = c(0, 12), expand = c(0, 0)) +
  labs(x = "Numero di piante", y = "Ricchezza in taxa (q = 0)", title = "Rarefazione ed estrapolazione",
       colour = NULL, fill = NULL, shape = NULL) + TEMA + theme(legend.position = c(0.7, 0.22), legend.text = element_text(size = 7.5))
p2 <- ggplot(s0, aes(t, SC, colour = gruppo, fill = gruppo)) + zona +
  geom_ribbon(aes(ymin = SC.LCL, ymax = SC.UCL), alpha = 0.15, colour = NA) +
  geom_line(data = s0[s0$t <= 10, ], linewidth = 0.7) + geom_line(data = s0[s0$t >= 10, ], linewidth = 0.7, linetype = "22") +
  geom_point(data = oss, aes(shape = gruppo), size = 2) +
  scale_colour_manual(values = COL) + scale_fill_manual(values = COL) + scale_shape_manual(values = SHP) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20)) +
  scale_y_continuous(limits = c(0.3, 1.01), labels = function(x) virgola(x, 1), expand = c(0, 0)) +
  labs(x = "Numero di piante", y = "Copertura del campione", title = "Completezza del campionamento") + TEMA
salva_fig(p1 + p2 + plot_annotation(tag_levels = "A"), "F02_rarefazione_estrapolazione", 170, 80)
