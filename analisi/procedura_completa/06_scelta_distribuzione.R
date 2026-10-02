# =============================================================================
# 06_scelta_distribuzione.R - Passaggio 6: scelta della famiglia di distribuzione per la ricchezza
# Confronto per AIC tra Poisson, Poisson generalizzata (genpois) e Conway-Maxwell-Poisson (compois)
# Funzioni: glmmTMB::glmmTMB(); confronto con lme4::glmer(family = poisson)
# La binomiale negativa non è considerata: ammette solo varianza > media (sovradispersione).
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
fam <- list(Poisson = poisson(), Poisson_generalizzata = genpois(), Conway_Maxwell_Poisson = compois())
riga <- function(nome, f, dd, form) {
  fit <- glmmTMB(form, data = dd, family = f)
  data.frame(Famiglia = nome, AIC = round(AIC(fit), 1), logLik = round(as.numeric(logLik(fit)), 2),
             Convergenza = fit$fit$convergence == 0,
             Var_pianta = if (length(VarCorr(fit)$cond)) round(VarCorr(fit)$cond$pl[1], 4) else NA)
}
t1 <- do.call(rbind, lapply(names(fam), function(n) riga(n, fam[[n]], S, R ~ g + (1 | pl))))
t4 <- do.call(rbind, lapply(names(fam), function(n) riga(n, fam[[n]], H, R ~ g)))
tab <- rbind(cbind(Modello = "Domanda 1: R ~ stagione + (1|pianta)", t1), cbind(Modello = "Domanda 4: R ~ ospite", t4))
tab$dAIC <- ave(tab$AIC, tab$Modello, FUN = function(a) a - min(a))
print(tab, row.names = FALSE); salva_tab(tab, "T10_scelta_famiglia")

# Confronto con lme4 (solo Poisson disponibile per i conteggi; binomiale negativa solo per sovradispersione)
g1 <- lme4::glmer(R ~ g + (1 | pl), data = S, family = poisson)
lme <- data.frame(logLik_glmer = as.numeric(logLik(g1)), Var_pianta = lme4::VarCorr(g1)$pl[1],
                  Adattamento_singolare = lme4::isSingular(g1),
                  Dispersione_Pearson = sum(residuals(g1, "pearson")^2) / df.residual(g1))
print(lme); salva_tab(lme, "T11_confronto_lme4")

# Figura: frequenze osservate e attese (Poisson vs Poisson generalizzata), domanda 4 come esempio a una riga per pianta
fP <- glmmTMB(R ~ g, data = H, family = poisson()); fG <- glmmTMB(R ~ g, data = H, family = genpois())
dgp <- function(y, mu, phi) {   # densità della Poisson generalizzata nella parametrizzazione di glmmTMB (var = mu*phi)
  th <- mu / sqrt(phi); la <- 1 - 1 / sqrt(phi)
  v <- suppressWarnings(exp(log(th) + (y - 1) * log(th + la * y) - th - la * y - lgamma(y + 1)))
  v[!is.finite(v) | th + la * y <= 0] <- 0; v }        # verificata: riproduce logLik() di glmmTMB
yy <- 0:14
attP <- sapply(yy, function(y) sum(dpois(y, fitted(fP))))
attG <- sapply(yy, function(y) sum(dgp(y, fitted(fG), sigma(fG))))
oss <- as.numeric(table(factor(H$R, levels = yy)))
fr <- rbind(data.frame(y = yy, n = attP, mod = "Poisson"), data.frame(y = yy, n = attG, mod = "Poisson generalizzata"))
pF <- ggplot() + geom_col(data = data.frame(y = yy, n = oss), aes(y, n), fill = "grey80", width = 0.8) +
  geom_line(data = fr, aes(y, n, linetype = mod), linewidth = 0.6) + geom_point(data = fr, aes(y, n, shape = mod), size = 1.6) +
  scale_shape_manual(values = c(1, 16)) + scale_x_continuous(breaks = yy) +
  labs(x = "Ricchezza in taxa per pianta", y = "Numero di piante", linetype = NULL, shape = NULL,
       title = "Frequenze osservate (barre) e attese dai modelli",
       subtitle = "Estate, 20 piante (domanda 4)") + TEMA + theme(legend.position = c(0.2, 0.85))
salva_fig(pF, "F04_scelta_distribuzione", 120, 80)
