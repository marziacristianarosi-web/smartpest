# =============================================================================
# 05_distribuzione.R - Passaggio 5: distribuzione delle variabili risposta
# (a) ricchezza per pianta: rapporto varianza/media e test di dispersione di Poisson
#     D = (n-1)·s²/media ~ chi² con n-1 gl; P(D <= osservato) = evidenza di sottodispersione
# (b) incidenza dei singoli taxa: struttura delle tabelle 2x2 e separazione
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")

disp <- do.call(rbind, lapply(levels(PIANTE$gruppo), function(g) {
  r <- PIANTE$R[PIANTE$gruppo == g]; n <- length(r); D <- (n - 1) * var(r) / mean(r)
  data.frame(Gruppo = ETI[g], Valori = paste(sort(r), collapse = " "), Media = mean(r), Varianza = round(var(r), 2),
             Var_su_media = round(var(r) / mean(r), 2), Chi2 = round(D, 2), gl = n - 1,
             P_sottodispersione = round(pchisq(D, n - 1), 4))
}))
print(disp, row.names = FALSE); salva_tab(disp, "T08_dispersione_ricchezza")

# differenze appaiate estate - autunno
dif <- S$R[S$g == "Summer"] - S$R[S$g == "Fall"]
cat("Differenze appaiate estate - autunno:", dif, "| media", mean(dif), "\n")

# struttura delle tabelle 2x2 per taxon (separazione)
classe <- function(a, b, n = 10) {
  if (a + b == 0) "assente in entrambi" else if (a + b == 2 * n) "presente in tutte le piante" else if ((a == 0 & b == n) | (a == n & b == 0)) "separazione completa"
  else if (a %in% c(0, n) | b %in% c(0, n)) "separazione quasi completa" else "nessuna separazione" }
sep <- rbind(
  data.frame(Confronto = "Stagione", Taxon = TAXA, Gruppo1 = colSums(PA_S[S$g == "Summer", ]), Gruppo2 = colSums(PA_S[S$g == "Fall", ])),
  data.frame(Confronto = "Ospite", Taxon = TAXA, Gruppo1 = colSums(PA_H[H$g == "Quercus_robur", ]), Gruppo2 = colSums(PA_H[H$g == "Quercus_cerris", ])))
sep$Struttura <- mapply(classe, sep$Gruppo1, sep$Gruppo2)
print(table(sep$Confronto, sep$Struttura)); salva_tab(sep, "T09_separazione_taxa")

# Figura: distribuzione della ricchezza con la varianza attesa sotto Poisson
pd <- PIANTE; pd$gr <- pd$gruppo
lim <- data.frame(gr = levels(pd$gr), m = tapply(pd$R, pd$gr, mean)[levels(pd$gr)])
lim$inf <- lim$m - sqrt(lim$m); lim$sup <- lim$m + sqrt(lim$m)            # media ± DS attesa sotto Poisson
lim$sinf <- lim$m - tapply(pd$R, pd$gr, sd)[lim$gr]; lim$ssup <- lim$m + tapply(pd$R, pd$gr, sd)[lim$gr]
lim$x <- as.numeric(factor(lim$gr, levels = levels(pd$gr)))
p1 <- ggplot(pd, aes(gr, R)) +
  geom_rect(data = lim, aes(xmin = x + 0.22, xmax = x + 0.36, ymin = inf, ymax = sup), inherit.aes = FALSE, fill = "grey85") +
  geom_linerange(data = lim, aes(x = x + 0.29, ymin = sinf, ymax = ssup), inherit.aes = FALSE, linewidth = 0.8) +
  geom_point(data = lim, aes(x = x + 0.29, y = m), inherit.aes = FALSE, shape = 23, fill = "white", size = 2.2) +
  geom_boxplot(width = 0.42, fill = NA, colour = "grey55", outlier.shape = NA, linewidth = 0.35) +
  geom_dotplot(aes(fill = gr), binaxis = "y", stackdir = "center", binwidth = 0.35, dotsize = 0.9, colour = NA) +
  scale_fill_manual(values = COL) +
  scale_x_discrete(labels = c(farnia_estate = "Q. robur\nestate", farnia_autunno = "Q. robur\nautunno", cerro_estate = "Q. cerris\nestate")) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 2)) +
  labs(x = NULL, y = "Ricchezza in taxa per pianta", title = "Ricchezza per pianta",
       subtitle = "Boxplot e punti: piante; rombo e barra: media ± DS\nFascia grigia: media ± DS attesa con Poisson") + TEMA
disp$Gruppo <- factor(disp$Gruppo, levels = ETI[levels(PIANTE$gruppo)])
p2 <- ggplot(disp, aes(Gruppo, Var_su_media)) +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey40") +
  geom_col(aes(fill = levels(PIANTE$gruppo)), width = 0.55) +
  annotate("text", x = 0.6, y = 1.06, label = "Poisson: varianza = media", hjust = 0, size = 2.6, family = FONT) +
  scale_fill_manual(values = COL) + scale_y_continuous(limits = c(0, 1.15), expand = c(0, 0), labels = function(x) virgola(x, 1)) +
  labs(x = NULL, y = "Varianza / media", title = "Indice di dispersione", subtitle = "Valori < 1 = sottodispersione\n ") + TEMA +
  scale_x_discrete(labels = function(x) sub(", ", "\n", x))
salva_fig(p1 + p2 + plot_layout(widths = c(1.3, 1)) + plot_annotation(tag_levels = "A"), "F03_distribuzione_ricchezza", 170, 80)
