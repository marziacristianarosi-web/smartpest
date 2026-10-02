# =============================================================================
# 09_potenza.R - Passaggio 9: potenza per la ricchezza con 10 piante per gruppo
# Simulazione dal modello di Poisson generalizzata (glmmTMB::simulate_new) con la variabilità
# stimata dai dati e differenze fissate a priori (0,5-3 taxa per pianta) attorno alla media generale.
# Valore critico dell'LRT ricavato dalle simulazioni senza effetto (non dal chi², che è anticonservativo).
# La potenza è una proprietà del disegno: non usa l'effetto osservato.
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
set.seed(2026)
fS <- glmmTMB(R ~ g + (1 | pl), data = S, family = genpois()); fH <- glmmTMB(R ~ g, data = H, family = genpois())
par <- list(
  `Domanda 1 (stagione, appaiato)` = list(dat = S, M = mean(S$R), betad = fixef(fS)$disp,
                                         theta = fS$fit$par[names(fS$fit$par) == "theta"], re = TRUE),
  `Domanda 4 (ospite, indipendente)` = list(dat = H, M = mean(H$R), betad = fixef(fH)$disp, theta = NULL, re = FALSE))
simula_lrt <- function(p, delta, nsim) {
  mu0 <- p$M - delta / 2; mu1 <- p$M + delta / 2
  np <- list(beta = c(log(mu0), log(mu1) - log(mu0)), betad = p$betad); if (p$re) np$theta <- p$theta
  ys <- simulate_new(if (p$re) ~ g + (1 | pl) else ~ g, newdata = p$dat, family = genpois(), newparams = np, nsim = nsim)
  unlist(in_parallelo(ys, function(y, dat, re) { dd <- dat; dd$R <- y
    if (re) lrt_gp(R ~ g + (1 | pl), R ~ 1 + (1 | pl), dd) else lrt_gp(R ~ g, R ~ 1, dd) }, dat = p$dat, re = p$re))
}
DELTA <- c(0.5, 1, 1.5, 2, 2.5, 3)
pot <- do.call(rbind, lapply(names(par), function(k) { p <- par[[k]]
  crit <- quantile(simula_lrt(p, 0, 400), 0.95, na.rm = TRUE)
  do.call(rbind, lapply(DELTA, function(dl) { L <- simula_lrt(p, dl, 200); pw <- mean(L > crit, na.rm = TRUE)
    data.frame(Domanda = k, Differenza_taxa = dl, Valore_critico = round(crit, 2), Potenza = pw,
               IC_inf = max(0, pw - 1.96 * sqrt(pw * (1 - pw) / 200)), IC_sup = min(1, pw + 1.96 * sqrt(pw * (1 - pw) / 200))) }))
}))
print(pot, row.names = FALSE); salva_tab(pot, "T15_potenza")
pP <- ggplot(pot, aes(Differenza_taxa, Potenza, linetype = Domanda, shape = Domanda)) +
  geom_hline(yintercept = 0.8, colour = "grey60", linetype = "22") +
  geom_errorbar(aes(ymin = IC_inf, ymax = IC_sup), width = 0.06, linetype = 1, colour = "grey50") +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(100 * x, "%")) + scale_x_continuous(breaks = DELTA) +
  labs(x = "Differenza vera di ricchezza (taxa per pianta)", y = "Potenza (α = 0,05)", linetype = NULL, shape = NULL,
       title = "Potenza del test sulla ricchezza con 10 piante per gruppo",
       subtitle = "200 simulazioni per punto; linea grigia: potenza 80%") + TEMA + theme(legend.position = c(0.7, 0.25))
salva_fig(pP, "F07_potenza", 130, 85)
