# =============================================================================
# 08_errore_primo_tipo.R - Passaggio 8: i test asintotici mantengono il 5% di falsi positivi con n = 10?
# (a) Calibrazione per simulazione: 500 insiemi di dati simulati dal modello SENZA effetto
#     (glmmTMB::simulate), test di Wald e del rapporto di verosimiglianza (LRT) con riferimento chi².
# (b) Distribuzione di permutazione esatta entro pianta (1024 permutazioni) dell'LRT della domanda 1:
#     quota di permutazioni in cui il test chi² sarebbe "significativo" (atteso 5%).
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
set.seed(11)
NSIM <- 500
RICALCOLA <- as.logical(Sys.getenv("CORY_RICALCOLA", "TRUE"))   # FALSE = rifà solo la figura dai file salvati

calibra <- function(dat, f1, f0, fam_nome) {
  fam <- switch(fam_nome, poisson = poisson(), genpois = genpois())
  nullo <- glmmTMB(f0, data = dat, family = fam); sims <- simulate(nullo, nsim = NSIM)
  r <- in_parallelo(sims, function(y, dat, f1, f0, fam_nome) {
    fam <- switch(fam_nome, poisson = poisson(), genpois = genpois())
    dd <- dat; dd$R <- y
    a <- try(suppressWarnings(glmmTMB(f1, data = dd, family = fam)), silent = TRUE)
    b <- try(suppressWarnings(glmmTMB(f0, data = dd, family = fam)), silent = TRUE)
    if (inherits(a, "try-error") || inherits(b, "try-error") || a$fit$convergence != 0) return(c(NA, NA))
    c(summary(a)$coefficients$cond[2, 4], pchisq(as.numeric(2 * (logLik(a) - logLik(b))), 1, lower.tail = FALSE))
  }, dat = dat, f1 = f1, f0 = f0, fam_nome = fam_nome)
  r <- do.call(rbind, r)
  c(Falliti = mean(is.na(r[, 1])), Wald = mean(r[, 1] < 0.05, na.rm = TRUE), LRT = mean(r[, 2] < 0.05, na.rm = TRUE))
}
if (RICALCOLA || !file.exists(file.path(DIR_TAB, "T13_calibrazione_errore_I_tipo.xlsx"))) {
cal <- do.call(rbind, lapply(c("poisson", "genpois"), function(fm) rbind(
  data.frame(Domanda = "1 (GLMM, 10 coppie)", Famiglia = fm, t(calibra(S, R ~ g + (1 | pl), R ~ 1 + (1 | pl), fm)) * 100),
  data.frame(Domanda = "4 (GLM, 10 + 10)", Famiglia = fm, t(calibra(H, R ~ g, R ~ 1, fm)) * 100))))
cal[, 3:5] <- round(cal[, 3:5], 1); names(cal)[3:5] <- c("Falliti_%", "Wald_%", "LRT_%")
print(cal, row.names = FALSE); salva_tab(cal, "T13_calibrazione_errore_I_tipo")
} else cal <- as.data.frame(read_excel(file.path(DIR_TAB, "T13_calibrazione_errore_I_tipo.xlsx")), check.names = FALSE)

# (b) distribuzione di permutazione entro pianta dell'LRT (usata anche dal test della domanda 1)
fileL1 <- file.path(DIR_TAB, "D1_distribuzione_permutazione.rds")
L1 <- if (!RICALCOLA && file.exists(fileL1)) readRDS(fileL1) else unlist(in_parallelo(seq_len(nrow(PERM_ENTRO)), function(i, S, P) {
  dd <- S; dd$g <- S$g[P[i, ]]; lrt_gp(R ~ g + (1 | pl), R ~ 1 + (1 | pl), dd) }, S = S, P = PERM_ENTRO))
saveRDS(L1, fileL1)
perm_tab <- data.frame(Permutazioni_valide = sum(!is.na(L1)),
                       Quota_p_chi2_sotto_0.05 = round(100 * mean(pchisq(L1, 1, lower.tail = FALSE) < 0.05, na.rm = TRUE), 1),
                       Percentile95_permutazione = round(quantile(L1, 0.95, na.rm = TRUE), 2), Percentile95_chi2 = 3.84)
print(perm_tab); salva_tab(perm_tab, "T14_permutazione_vs_chi2_D1")

d <- data.frame(L = L1[!is.na(L1) & L1 < 15])
p1 <- ggplot(d, aes(L)) + geom_histogram(aes(y = after_stat(density)), binwidth = 0.25, boundary = 0, fill = "grey75", colour = "white") +
  stat_function(fun = function(x) dchisq(x, 1), colour = "black", linewidth = 0.6, n = 400, xlim = c(0.05, 15)) +
  geom_vline(xintercept = 3.84, linetype = "22") + geom_vline(xintercept = perm_tab$Percentile95_permutazione, colour = "#c0392b", linewidth = 0.7) +
  coord_cartesian(ylim = c(0, 1.2)) +
  labs(x = "Statistica LRT", y = "Densità", title = "Domanda 1: distribuzione nulla dell'LRT",
       subtitle = sprintf("Barre: 1024 permutazioni entro pianta\nCurva: chi² (1 gl); tratteggio: 3,84\nRosso: 95° percentile di permutazione (%s)",
                          virgola(perm_tab$Percentile95_permutazione))) + TEMA
cal$Famiglia <- c(poisson = "Poisson", genpois = "Poisson gen.")[cal$Famiglia]
cl <- reshape(cal[, c("Domanda", "Famiglia", "Wald_%", "LRT_%")], direction = "long", varying = c("Wald_%", "LRT_%"),
              v.names = "tasso", timevar = "Test", times = c("Wald", "LRT"))
p2 <- ggplot(cl, aes(Famiglia, tasso, fill = Test)) + geom_col(position = position_dodge(0.7), width = 0.65) +
  scale_fill_manual(values = c(LRT = "grey35", Wald = "grey75"), name = NULL) +
  geom_hline(yintercept = 5, linetype = "22") + facet_wrap(~ Domanda) +
  labs(x = NULL, y = "Falsi positivi (%)", title = "Calibrazione per simulazione",
       subtitle = "500 insiemi di dati senza effetto\nTratteggio: livello nominale 5%\n ") + TEMA + theme(strip.background = element_blank(), legend.position = "top")
salva_fig(p1 + p2 + plot_layout(widths = c(1, 1.2)) + plot_annotation(tag_levels = "A"), "F06_errore_primo_tipo", 180, 85)
