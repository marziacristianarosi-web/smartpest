# =============================================================================
# 07_assunzioni.R - A.7: verifica delle assunzioni dei modelli per la ricchezza
# Residui quantilici randomizzati (Dunn & Smyth 1996), come in DHARMa (Hartig): per ogni osservazione
# r = U(P(Y < y), P(Y <= y)) su 1000 simulazioni dal modello stimato. Se il modello è corretto, r ~ Uniforme(0, 1).
# Nei modelli misti le simulazioni sono CONDIZIONATE agli effetti di pianta stimati (impostazione predefinita di
# DHARMa dalla versione 0.5.0): si simula solo il livello della distribuzione, con le medie stimate per ciascuna pianta.
# A1 forma della distribuzione (KS sull'uniforme)     A2 dispersione residua (osservata / simulata)
# A3 valori anomali                                  A4 dispersione uguale nei gruppi (dispformula)
# A5 osservazioni influenti (esclusione di una pianta) A6 normalità degli effetti casuali
# A7 limite superiore (14 taxa)
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
set.seed(7)
NSIM <- 1000
dgp <- function(y, mu, phi) {   # densità della Poisson generalizzata (parametrizzazione di glmmTMB: var = mu*phi)
  th <- mu / sqrt(phi); la <- 1 - 1 / sqrt(phi)
  v <- suppressWarnings(exp(log(th) + (y - 1) * log(th + la * y) - th - la * y - lgamma(y + 1)))
  v[!is.finite(v) | th + la * y <= 0] <- 0; v }
sim_cond <- function(fit, nsim) {   # simulazione condizionata: medie stimate (con gli effetti di pianta) e dispersione stimata
  mu <- fitted(fit); phi <- sigma(fit); sup <- 0:40
  P <- sapply(mu, function(m) { p <- dgp(sup, m, phi); p / sum(p) })
  sapply(seq_len(nsim), function(k) apply(P, 2, function(p) sample(sup, 1, prob = p)))
}
verifica <- function(nome, dat, form, re) {
  fit <- glmmTMB(form, data = dat, family = genpois()); y <- dat$R
  sims <- sim_cond(fit, NSIM); mu <- fitted(fit)
  r <- sapply(seq_along(y), function(i) runif(1, mean(sims[i, ] < y[i]), mean(sims[i, ] <= y[i])))
  d_oss <- sum((y - mu)^2); d_sim <- apply(sims, 2, function(s) sum((s - mu)^2))
  het <- glmmTMB(form, dispformula = ~ g, data = dat, family = genpois())
  lr <- as.numeric(2 * (logLik(het) - logLik(fit)))
  b <- fixef(fit)$cond[2]; ids <- unique(as.character(dat$N_pianta))
  loo <- sapply(ids, function(p) { f <- try(suppressWarnings(glmmTMB(form, data = dat[dat$N_pianta != p, ], family = genpois())), silent = TRUE)
                                    if (inherits(f, "try-error")) NA else fixef(f)$cond[2] })
  ris <- data.frame(Modello = nome, N_oss = length(y),
    A1_KS_p = round(ks.test(r, "punif")$p.value, 3),
    A2_dispersione_oss_su_sim = round(d_oss / mean(d_sim), 2),
    A2_p = round(2 * min(mean(d_sim >= d_oss), mean(d_sim <= d_oss)), 3),
    A3_outlier = sum(y > apply(sims, 1, max) | y < apply(sims, 1, min)),
    A4_LRT_dispersione_p = round(pchisq(lr, 1, lower.tail = FALSE), 3),
    A5_rapporto_medie = round(exp(b), 2), A5_min_loo = round(exp(min(loo, na.rm = TRUE)), 2),
    A5_max_loo = round(exp(max(loo, na.rm = TRUE)), 2),
    A6_Shapiro_p = if (re) round(shapiro.test(ranef(fit)$cond$pl[, 1])$p.value, 3) else NA,
    A7_P_oltre_14 = mean(sims > 14))
  list(ris = ris, r = r, mu = mu, d_oss = d_oss, d_sim = d_sim, dat = dat)
}
v1 <- verifica("D1: GLMM, R ~ stagione + (1|pianta)", S, R ~ g + (1 | pl), TRUE)
v4 <- verifica("D4 per pianta: GLM, R ~ ospite", H, R ~ g, FALSE)
v4p <- verifica("D4 per piastra: GLMM, R ~ ospite + (1|pianta)", E, R ~ g + (1 | pl), TRUE)
tab <- rbind(v1$ris, v4$ris, v4p$ris); print(t(tab)); salva_tab(tab, "T12_assunzioni")

ETI_G <- c(Summer = "Estate", Fall = "Autunno", Quercus_robur = "Q. robur", Quercus_cerris = "Q. cerris")
pannelli <- function(v, titolo) {
  n <- length(v$r); d <- data.frame(att = (seq_len(n) - 0.5) / n, oss = sort(v$r))
  a <- ggplot(d, aes(att, oss)) + geom_abline(colour = "grey50", linetype = "22") + geom_point(size = 1.2) +
    coord_equal() + scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Quantili attesi (uniforme)", y = "Residui quantilici", title = titolo,
         subtitle = sprintf("A1 KS p = %s", virgola(v$ris$A1_KS_p, 3))) + TEMA
  e <- data.frame(pred = rank(v$mu, ties.method = "average") / n, r = v$r, g = v$dat$g)
  b <- ggplot(e, aes(pred, r)) + geom_hline(yintercept = c(0.25, 0.5, 0.75), colour = "grey70", linetype = "22") +
    geom_jitter(aes(shape = g), width = 0.01, height = 0, size = 1.4) + scale_shape_manual(values = c(16, 2), labels = ETI_G) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Valori predetti (rango)", y = "Residui quantilici", shape = NULL, title = "Residui vs predetti") +
    TEMA + theme(legend.position = "bottom")
  c <- ggplot(data.frame(s = v$d_sim), aes(s)) + geom_histogram(bins = 30, fill = "grey75", colour = "white") +
    geom_vline(xintercept = v$d_oss, colour = "#c0392b", linewidth = 0.8) +
    labs(x = "Somma dei quadrati dei residui", y = "Simulazioni", title = "Dispersione residua",
         subtitle = sprintf("A2 oss/sim = %s; p = %s", virgola(v$ris$A2_dispersione_oss_su_sim), virgola(v$ris$A2_p, 3))) + TEMA
  a + b + c
}
salva_fig(pannelli(v1, "D1 (GLMM)") / pannelli(v4, "D4 per pianta (GLM)") / pannelli(v4p, "D4 per piastra (GLMM)") +
          plot_annotation(tag_levels = "A"), "F05_diagnostica_assunzioni", 180, 200)
