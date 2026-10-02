# =============================================================================
# 07_assunzioni.R - Passaggio 7: verifica delle assunzioni dei modelli per la ricchezza
# Residui quantilici randomizzati (Dunn & Smyth 1996; stesso principio di DHARMa, Hartig 2022):
# per ogni osservazione r = U(P(Y < y), P(Y <= y)) calcolato su 1000 simulazioni dal modello stimato.
# Se il modello è corretto i residui sono uniformi su (0, 1).
# A1 forma della distribuzione (KS sull'uniforme)   A2 dispersione residua (oss/sim)
# A3 valori anomali                                A4 dispersione uguale nei gruppi (dispformula)
# A5 osservazioni influenti (leave-one-out)        A6 normalità degli effetti casuali
# A7 limite superiore (14 taxa)
# Funzioni: glmmTMB(), simulate(), ks.test(), shapiro.test()
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
set.seed(7)
NSIM <- 1000

residui <- function(fit, y) {
  sims <- as.matrix(simulate(fit, nsim = NSIM))
  r <- sapply(seq_along(y), function(i) runif(1, mean(sims[i, ] < y[i]), mean(sims[i, ] <= y[i])))
  mu <- fitted(fit); d_oss <- sum((y - mu)^2); d_sim <- apply(sims, 2, function(s) sum((s - mu)^2))
  list(r = r, mu = mu, sims = sims, d_oss = d_oss, d_sim = d_sim)
}
verifica <- function(nome, dat, form, re) {
  fit <- glmmTMB(form, data = dat, family = genpois()); q <- residui(fit, dat$R)
  het <- glmmTMB(form, dispformula = ~ g, data = dat, family = genpois())
  lr <- as.numeric(2 * (logLik(het) - logLik(fit)))
  b <- fixef(fit)$cond[2]; ids <- unique(as.character(dat$N_pianta))
  loo <- sapply(ids, function(p) { f <- try(suppressWarnings(glmmTMB(form, data = dat[dat$N_pianta != p, ], family = genpois())), silent = TRUE)
                                    if (inherits(f, "try-error")) NA else fixef(f)$cond[2] })
  ris <- data.frame(Modello = nome,
    A1_KS_p = round(ks.test(q$r, "punif")$p.value, 3),
    A2_dispersione_oss_su_sim = round(q$d_oss / mean(q$d_sim), 2),
    A2_p = round(2 * min(mean(q$d_sim >= q$d_oss), mean(q$d_sim <= q$d_oss)), 3),
    A3_outlier = sum(dat$R > apply(q$sims, 1, max) | dat$R < apply(q$sims, 1, min)),
    A4_LRT_dispersione_p = round(pchisq(lr, 1, lower.tail = FALSE), 3),
    A5_rapporto_medie = round(exp(b), 2), A5_min_loo = round(exp(min(loo, na.rm = TRUE)), 2),
    A5_max_loo = round(exp(max(loo, na.rm = TRUE)), 2),
    A6_Shapiro_p = if (re) round(shapiro.test(ranef(fit)$cond$pl[, 1])$p.value, 3) else NA,
    A7_P_oltre_14 = mean(q$sims > 14))
  list(ris = ris, q = q, dat = dat, nome = nome)
}
v1 <- verifica("Domanda 1: GLMM, R ~ stagione + (1|pianta)", S, R ~ g + (1 | pl), TRUE)
v4 <- verifica("Domanda 4: GLM, R ~ ospite", H, R ~ g, FALSE)
tab <- rbind(v1$ris, v4$ris); print(t(tab)); salva_tab(tab, "T12_assunzioni")

# Figura diagnostica: 3 pannelli per modello
pannelli <- function(v, titolo) {
  n <- length(v$q$r); d <- data.frame(att = (seq_len(n) - 0.5) / n, oss = sort(v$q$r))
  a <- ggplot(d, aes(att, oss)) + geom_abline(colour = "grey50", linetype = "22") + geom_point(size = 1.4) +
    coord_equal() + scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Quantili attesi (uniforme)", y = "Residui quantilici", title = titolo,
         subtitle = sprintf("A1 KS p = %s", virgola(v$ris$A1_KS_p, 3))) + TEMA
  e <- data.frame(pred = rank(v$q$mu, ties.method = "average") / n, r = v$q$r, g = v$dat$g)
  b <- ggplot(e, aes(pred, r)) + geom_hline(yintercept = c(0.25, 0.5, 0.75), colour = "grey70", linetype = "22") +
    geom_jitter(aes(shape = g), width = 0.01, height = 0, size = 1.6) + scale_shape_manual(values = c(16, 2), labels = c(Summer = "Estate", Fall = "Autunno", Quercus_robur = "Q. robur", Quercus_cerris = "Q. cerris")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Valori predetti (rango)", y = "Residui quantilici", shape = NULL, title = "Residui vs predetti") +
    TEMA + theme(legend.position = "bottom")
  c <- ggplot(data.frame(s = v$q$d_sim), aes(s)) + geom_histogram(bins = 30, fill = "grey75", colour = "white") +
    geom_vline(xintercept = v$q$d_oss, colour = "#c0392b", linewidth = 0.8) +
    labs(x = "Somma dei quadrati dei residui", y = "Simulazioni", title = "Dispersione residua",
         subtitle = sprintf("A2 oss/sim = %s; p = %s", virgola(v$ris$A2_dispersione_oss_su_sim), virgola(v$ris$A2_p, 3))) + TEMA
  a + b + c
}
salva_fig(pannelli(v1, "Domanda 1 (GLMM)") / pannelli(v4, "Domanda 4 (GLM)") + plot_annotation(tag_levels = "A"),
          "F05_diagnostica_assunzioni", 180, 140)
