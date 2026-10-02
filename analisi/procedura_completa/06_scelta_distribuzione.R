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
# distribuzioni discrete: barre affiancate per ciascun valore intero (nessuna linea di raccordo)
fr <- rbind(data.frame(y = yy, n = oss, serie = "Osservate"),
            data.frame(y = yy, n = attP, serie = "Attese: Poisson"),
            data.frame(y = yy, n = attG, serie = "Attese: Poisson generalizzata"))
fr$serie <- factor(fr$serie, levels = c("Osservate", "Attese: Poisson", "Attese: Poisson generalizzata"))
pF <- ggplot(fr, aes(y, n, fill = serie)) +
  geom_col(position = position_dodge(0.85), width = 0.8, colour = "grey20", linewidth = 0.25) +
  scale_fill_manual(values = c("grey80", "white", "grey25"), name = NULL) +
  scale_x_continuous(breaks = yy, expand = c(0.01, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 7.5)) +
  labs(x = "Ricchezza in taxa per pianta (valori interi)", y = "Numero di piante",
       title = "Frequenze osservate e attese dai due modelli",
       subtitle = "Estate, 20 piante (domanda 4)\nFrequenze attese = somma delle probabilità previste per ciascuna pianta") +
  TEMA + theme(legend.position = c(0.2, 0.82))
salva_fig(pF, "F04_scelta_distribuzione", 120, 80)
