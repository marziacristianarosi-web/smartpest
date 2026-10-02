# =============================================================================
# 06_scelta_distribuzione.R - A.6: scelta della distribuzione per la ricchezza (confronto per AIC)
# Candidati: Poisson (varianza = media), binomiale negativa (varianza > media), binomiale su 14 taxa
# (conteggio limitato), Poisson generalizzata e Conway-Maxwell-Poisson (ammettono varianza < media).
# Modelli: D1 R ~ stagione + (1|pianta); D4 per pianta R ~ ospite; D4 per piastra R ~ ospite + (1|pianta).
# Funzioni: glmmTMB::glmmTMB(); confronto con lme4::glmer(family = poisson).
# Figure: frequenze osservate e attese (curve; barre) e rootogramma sospeso (Kleiber & Zeileis 2016).
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")
NTAXA <- length(TAXA)
adatta <- function(fam, dd, rhs) {
  dd$NT <- NTAXA
  f <- if (fam == "Binomiale (14 taxa)") as.formula(paste("cbind(R, NT - R) ~", rhs)) else as.formula(paste("R ~", rhs))
  fm <- switch(fam, "Poisson" = poisson(), "Binomiale negativa" = nbinom2(), "Binomiale (14 taxa)" = binomial(),
               "Poisson generalizzata" = genpois(), "Conway-Maxwell-Poisson" = compois())
  fit <- try(suppressWarnings(glmmTMB(f, data = dd, family = fm)), silent = TRUE)
  if (inherits(fit, "try-error")) return(data.frame(Famiglia = fam, AIC = NA, Convergenza = FALSE))
  ok <- fit$fit$convergence == 0 && isTRUE(fit$sdr$pdHess)
  data.frame(Famiglia = fam, AIC = round(AIC(fit), 1), Convergenza = ok)
}
FAM <- c("Poisson", "Binomiale negativa", "Binomiale (14 taxa)", "Poisson generalizzata", "Conway-Maxwell-Poisson")
mod <- list(`D1: R ~ stagione + (1|pianta)` = list(S, "g + (1 | pl)"),
            `D4 per pianta: R ~ ospite` = list(H, "g"),
            `D4 per piastra: R ~ ospite + (1|pianta)` = list(E, "g + (1 | pl)"))
tab <- do.call(rbind, lapply(names(mod), function(m) cbind(Modello = m, do.call(rbind, lapply(FAM, adatta, dd = mod[[m]][[1]], rhs = mod[[m]][[2]])))))
tab$dAIC <- ave(tab$AIC, tab$Modello, FUN = function(a) a - min(a, na.rm = TRUE))
print(tab, row.names = FALSE); salva_tab(tab, "T10_scelta_famiglia")

# Confronto con lme4 (per i conteggi offre Poisson e binomiale negativa)
g1 <- lme4::glmer(R ~ g + (1 | pl), data = S, family = poisson)
lme <- data.frame(logLik_glmer = as.numeric(logLik(g1)), Var_pianta = lme4::VarCorr(g1)$pl[1],
                  Adattamento_singolare = lme4::isSingular(g1),
                  Dispersione_Pearson = sum(residuals(g1, "pearson")^2) / df.residual(g1))
print(lme); salva_tab(lme, "T11_confronto_lme4")

# ---- Frequenze osservate e attese (domanda 4, una riga per pianta) ----
fP <- glmmTMB(R ~ g, data = H, family = poisson()); fG <- glmmTMB(R ~ g, data = H, family = genpois())
H$NT <- NTAXA; fB <- glmmTMB(cbind(R, NT - R) ~ g, data = H, family = binomial())
dgp <- function(y, mu, phi) {   # densità della Poisson generalizzata nella parametrizzazione di glmmTMB (var = mu*phi)
  th <- mu / sqrt(phi); la <- 1 - 1 / sqrt(phi)
  v <- suppressWarnings(exp(log(th) + (y - 1) * log(th + la * y) - th - la * y - lgamma(y + 1)))
  v[!is.finite(v) | th + la * y <= 0] <- 0; v }        # verificata: riproduce logLik() di glmmTMB
yy <- 0:14
att <- rbind(data.frame(y = yy, n = sapply(yy, function(y) sum(dpois(y, fitted(fP)))), mod = "Poisson"),
             data.frame(y = yy, n = sapply(yy, function(y) sum(dbinom(y, NTAXA, fitted(fB)))), mod = "Binomiale (14 taxa)"),
             data.frame(y = yy, n = sapply(yy, function(y) sum(dgp(y, fitted(fG), sigma(fG)))), mod = "Poisson generalizzata"))
att$mod <- factor(att$mod, levels = c("Poisson", "Binomiale (14 taxa)", "Poisson generalizzata"))
oss <- data.frame(y = yy, n = as.numeric(table(factor(H$R, levels = yy))))
sub <- "Estate, 20 piante (domanda 4); frequenze attese = somma delle probabilità previste per ciascuna pianta"

# (A) curve: punti delle frequenze attese uniti da linee
pC <- ggplot() + geom_col(data = oss, aes(y, n), fill = "grey80", width = 0.8) +
  geom_line(data = att, aes(y, n, linetype = mod), linewidth = 0.6) + geom_point(data = att, aes(y, n, shape = mod), size = 1.6) +
  scale_shape_manual(values = c(1, 2, 16)) + scale_linetype_manual(values = c("solid", "dotted", "22")) +
  scale_x_continuous(breaks = yy) +
  labs(x = "Ricchezza in taxa per pianta", y = "Numero di piante", linetype = NULL, shape = NULL,
       title = "Frequenze osservate (barre) e attese (curve)", subtitle = sub) + TEMA + theme(legend.position = c(0.2, 0.8))
salva_fig(pC, "F04A_frequenze_curve", 140, 85)

# (B) barre affiancate per ciascun valore intero
fr <- rbind(cbind(oss, serie = "Osservate"), setNames(att, c("y", "n", "serie")))
fr$serie <- factor(fr$serie, levels = c("Osservate", levels(att$mod)))
pBar <- ggplot(fr, aes(y, n, fill = serie)) +
  geom_col(position = position_dodge(0.85), width = 0.8, colour = "grey20", linewidth = 0.25) +
  scale_fill_manual(values = c("grey80", "white", "grey60", "grey20"), name = NULL) +
  scale_x_continuous(breaks = yy, expand = c(0.01, 0)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 7.5)) +
  labs(x = "Ricchezza in taxa per pianta (valori interi)", y = "Numero di piante",
       title = "Frequenze osservate e attese (barre)", subtitle = sub) + TEMA + theme(legend.position = c(0.2, 0.8))
salva_fig(pBar, "F04B_frequenze_barre", 140, 85)

# (C) rootogramma sospeso: barre di altezza sqrt(osservate) appese alla curva sqrt(attese);
# se il modello è adeguato il fondo delle barre si trova vicino alla linea dello zero
rt <- merge(att, oss, by = "y", suffixes = c("_att", "_oss"))
rt$top <- sqrt(rt$n_att); rt$bottom <- sqrt(rt$n_att) - sqrt(rt$n_oss)
pR <- ggplot(rt) + geom_hline(yintercept = 0, colour = "grey40") +
  geom_rect(aes(xmin = y - 0.4, xmax = y + 0.4, ymin = bottom, ymax = top), fill = "grey80", colour = "grey40", linewidth = 0.25) +
  geom_line(aes(y, top), colour = "#c0392b", linewidth = 0.6) + geom_point(aes(y, top), colour = "#c0392b", size = 1.2) +
  facet_wrap(~ mod, nrow = 1) + scale_x_continuous(breaks = seq(0, 14, 2)) +
  labs(x = "Ricchezza in taxa per pianta", y = "√ frequenza",
       title = "Rootogramma sospeso",
       subtitle = "Curva rossa: √ frequenze attese; barre appese: √ frequenze osservate; adattamento buono = fondo delle barre vicino a 0") +
  TEMA + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
salva_fig(pR, "F04C_rootogramma", 180, 80)

# ---- Piastre separate: GLM (piastre indipendenti) o GLMM (piastre raggruppate nella pianta)? ----
# Le 4 piastre di una pianta condividono albero e popolazione di insetti: se si somigliano più del caso, un GLM
# tratta 80 piastre come repliche indipendenti (pseudoreplicazione) e sottostima l'errore standard.
gA <- glmmTMB(R ~ g, dispformula = ~ g, data = E, family = genpois())
gB <- glmmTMB(R ~ g + (1 | pl), dispformula = ~ g, data = E, family = genpois())
rp <- residuals(gA, "pearson"); mp <- tapply(rp, E$pl, mean)
cA <- summary(gA)$coefficients$cond[2, ]; cB <- summary(gB)$coefficients$cond[2, ]
glm_glmm <- data.frame(Modello = c("GLM: R ~ ospite (80 piastre indipendenti)", "GLMM: R ~ ospite + (1|pianta)"),
  Rapporto_medie = round(exp(c(cA[1], cB[1])), 3), Errore_standard = round(c(cA[2], cB[2]), 3),
  p_Wald_illustrativo = round(c(cA[4], cB[4]), 4), AIC = round(c(AIC(gA), AIC(gB)), 1),
  Varianza_tra_piante = c(NA, round(VarCorr(gB)$cond$pl[1], 4)),
  Var_medie_residui_per_pianta = c(round(var(mp), 3), NA), Var_attesa_se_indipendenti = c(round(var(rp) / 4, 3), NA))
print(glm_glmm, row.names = FALSE); salva_tab(glm_glmm, "T26_GLM_vs_GLMM_piastre")
