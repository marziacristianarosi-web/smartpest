# =============================================================================
# 08_errore_primo_tipo.R - A.8: quali metodi di test mantengono il 5% di falsi positivi con questi campioni?
# Confronto, su dati simulati SENZA effetto dal modello stimato (Poisson generalizzata), di:
#   1 Wald z           coefficiente / errore standard confrontato con la normale
#   2 Wald t           stessa statistica confrontata con t di Student, gl = piante - parametri a livello di pianta
#   3 LRT chi²         rapporto di verosimiglianza confrontato con chi² a 1 gl
#   4 LRT bootstrap    distribuzione dell'LRT ricavata simulando dal modello nullo stimato (bootstrap parametrico)
#   5 LRT permutazione distribuzione dell'LRT ricavata permutando le etichette secondo il disegno: livello alfa
#                      garantito solo se le osservazioni sono scambiabili sotto H0 (non richiede simulazione)
# Modelli: D1 (GLMM appaiato), D4 per pianta (GLM), D4 per piastra (GLMM, dispersione per ospite: vedi A.7, A4).
# Più la distribuzione di permutazione completa dell'LRT della domanda 1 (tutte le 1024 inversioni entro pianta).
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
set.seed(11)
NSIM <- 500; NSIM_B <- 200; B <- 39   # insiemi simulati; insiemi per il bootstrap; repliche bootstrap per insieme
RICALCOLA <- as.logical(Sys.getenv("CORY_RICALCOLA", "TRUE"))   # FALSE = rifà solo le figure dai file salvati

MOD <- list(
  D1 = list(nome = "D1 stagione (GLMM, 10 coppie)", dat = S, f1 = R ~ g + (1 | pl), f0 = R ~ 1 + (1 | pl), disp = ~ 1, gl = 9),
  D4 = list(nome = "D4 ospite per pianta (GLM, 10 + 10)", dat = H, f1 = R ~ g, f0 = R ~ 1, disp = ~ 1, gl = 18),
  D4p = list(nome = "D4 ospite per piastra (GLMM, 80 piastre in 20 piante)", dat = E, f1 = R ~ g + (1 | pl),
             f0 = R ~ 1 + (1 | pl), disp = ~ g, gl = 18))

un_test <- function(y, m, boot, B) {
  dd <- m$dat; dd$R <- y
  a <- try(suppressWarnings(glmmTMB(m$f1, dispformula = m$disp, data = dd, family = genpois())), silent = TRUE)
  b <- try(suppressWarnings(glmmTMB(m$f0, dispformula = m$disp, data = dd, family = genpois())), silent = TRUE)
  if (inherits(a, "try-error") || inherits(b, "try-error") || a$fit$convergence != 0) return(rep(NA, 5))
  z <- summary(a)$coefficients$cond[2, 3]; L <- max(0, as.numeric(2 * (logLik(a) - logLik(b))))
  pb <- NA
  if (boot) {   # bootstrap parametrico: simula dal modello nullo stimato su questo insieme di dati
    ys <- simulate(b, nsim = B)
    Lb <- sapply(ys, function(yb) { db <- dd; db$R <- yb
      a2 <- try(suppressWarnings(glmmTMB(m$f1, dispformula = m$disp, data = db, family = genpois())), silent = TRUE)
      b2 <- try(suppressWarnings(glmmTMB(m$f0, dispformula = m$disp, data = db, family = genpois())), silent = TRUE)
      if (inherits(a2, "try-error") || inherits(b2, "try-error")) NA else as.numeric(2 * (logLik(a2) - logLik(b2))) })
    pb <- (sum(Lb >= L - 1e-8, na.rm = TRUE) + 1) / (sum(!is.na(Lb)) + 1)
  }
  c(Wald_z = 2 * pnorm(-abs(z)), Wald_t = 2 * pt(-abs(z), m$gl), LRT_chi2 = pchisq(L, 1, lower.tail = FALSE), LRT_boot = pb, L = L)
}
calibra <- function(m, nsim, boot) {
  nullo <- glmmTMB(m$f0, dispformula = m$disp, data = m$dat, family = genpois())
  sims <- simulate(nullo, nsim = nsim)
  r <- in_parallelo(sims, function(y, m, boot, un_test, B) un_test(y, m, boot, B), m = m, boot = boot, un_test = un_test, B = B)
  do.call(rbind, r)
}
file_cal <- file.path(DIR_TAB, "T13_confronto_metodi_test.xlsx")
if (RICALCOLA || !file.exists(file_cal)) {
  cal <- do.call(rbind, lapply(names(MOD), function(k) {
    m <- MOD[[k]]
    r1 <- calibra(m, NSIM, FALSE); rb <- calibra(m, NSIM_B, TRUE)
    tasso <- function(p) round(100 * mean(p < 0.05, na.rm = TRUE), 1)
    se <- function(p) round(100 * sqrt(mean(p < 0.05, na.rm = TRUE) * (1 - mean(p < 0.05, na.rm = TRUE)) / sum(!is.na(p))), 1)
    data.frame(Modello = m$nome,
               Metodo = c("Wald z", sprintf("Wald t (%d gl)", m$gl), "LRT chi²", sprintf("LRT bootstrap parametrico (B = %d)", B), "LRT permutazione"),
               Falsi_positivi_perc = c(tasso(r1[, 1]), tasso(r1[, 2]), tasso(r1[, 3]), tasso(rb[, 4]), 5),
               ES_MonteCarlo = c(se(r1[, 1]), se(r1[, 2]), se(r1[, 3]), se(rb[, 4]), NA),
               Insiemi_simulati = c(rep(sum(!is.na(r1[, 1])), 3), sum(!is.na(rb[, 4])), NA),
               Nota = c("", "", "", "", "livello alfa se le unità sono scambiabili sotto H0 (non stimato per simulazione)"))
  }))
  print(cal, row.names = FALSE); salva_tab(cal, "T13_confronto_metodi_test")
} else cal <- as.data.frame(read_excel(file_cal))
# Livello nominale: con B = 39 repliche il p bootstrap vale (k + 1)/40; il criterio p < 0,05 è soddisfatto solo con k = 0,
# quindi il tasso atteso sotto H0 è 1/40 = 2,5%. Per gli altri metodi è il 5%.
cal$Livello_nominale_perc <- ifelse(grepl("bootstrap", cal$Metodo), 100 / (B + 1), 5)
cal$Nota[cal$Metodo == "LRT permutazione"] <- "livello alfa se le unità sono scambiabili sotto H0 (non stimato per simulazione)"
salva_tab(cal, "T13_confronto_metodi_test")

# distribuzione di permutazione completa entro pianta dell'LRT della domanda 1 (usata anche dal test, passaggio 10)
fileL1 <- file.path(DIR_TAB, "D1_distribuzione_permutazione.rds")
L1 <- if (!RICALCOLA && file.exists(fileL1)) readRDS(fileL1) else unlist(in_parallelo(seq_len(nrow(PERM_ENTRO)), function(i, S, P) {
  dd <- S; dd$g <- S$g[P[i, ]]; lrt_gp(R ~ g + (1 | pl), R ~ 1 + (1 | pl), dd) }, S = S, P = PERM_ENTRO))
saveRDS(L1, fileL1)
perm_tab <- data.frame(Permutazioni_valide = sum(!is.na(L1)),
                       Rigetto_condizionale_chi2_perc = round(100 * mean(pchisq(L1, 1, lower.tail = FALSE) < 0.05, na.rm = TRUE), 1),
                       Percentile95_permutazione = round(quantile(L1, 0.95, na.rm = TRUE), 2), Percentile95_chi2 = 3.84)
print(perm_tab); salva_tab(perm_tab, "T14_permutazione_vs_chi2_D1")

d <- data.frame(L = L1[!is.na(L1) & L1 < 15])
p1 <- ggplot(d, aes(L)) + geom_histogram(aes(y = after_stat(density)), binwidth = 0.25, boundary = 0, fill = "grey75", colour = "white") +
  stat_function(fun = function(x) dchisq(x, 1), colour = "black", linewidth = 0.6, n = 400, xlim = c(0.05, 15)) +
  geom_vline(xintercept = 3.84, linetype = "22") + geom_vline(xintercept = perm_tab$Percentile95_permutazione, colour = "#c0392b", linewidth = 0.7) +
  coord_cartesian(ylim = c(0, 1.2)) +
  labs(x = "Statistica LRT", y = "Densità", title = "D1: distribuzione nulla dell'LRT",
       subtitle = sprintf("Barre: 1024 permutazioni entro pianta\nCurva: chi² (1 gl); tratteggio: 3,84\nRosso: 95° percentile di permutazione (%s)",
                          virgola(perm_tab$Percentile95_permutazione))) + TEMA
cl <- cal; cl$Metodo <- sub(" \\(.*\\)", "", cl$Metodo); cl$Metodo <- factor(cl$Metodo, levels = unique(cl$Metodo))
cl$Modello <- factor(sub(" \\(.*", "", cl$Modello), levels = unique(sub(" \\(.*", "", cl$Modello)))
p2 <- ggplot(cl, aes(Metodo, Falsi_positivi_perc)) +
  geom_col(aes(fill = Metodo == "LRT permutazione"), width = 0.65, show.legend = FALSE) +
  geom_errorbar(aes(ymin = Falsi_positivi_perc - 1.96 * ES_MonteCarlo, ymax = Falsi_positivi_perc + 1.96 * ES_MonteCarlo), width = 0.2, na.rm = TRUE) +
  geom_errorbar(aes(ymin = Livello_nominale_perc, ymax = Livello_nominale_perc), width = 0.8, linetype = "22", colour = "#c0392b") +
  facet_wrap(~ Modello, ncol = 1) +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "grey85")) +
  labs(x = NULL, y = "Falsi positivi (%)", title = "Falsi positivi per metodo di test",
       subtitle = "Dati simulati senza effetto; barre = IC 95% Monte Carlo\nTratteggio rosso: livello nominale (5%; 2,5% per il bootstrap con B = 39)\npermutazione: livello alfa se le piante sono scambiabili sotto H0") +
  TEMA + theme(strip.background = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1))
salva_fig(p1 + p2 + plot_layout(widths = c(1, 1.3)) + plot_annotation(tag_levels = "A"), "F06_errore_primo_tipo", 190, 150)
