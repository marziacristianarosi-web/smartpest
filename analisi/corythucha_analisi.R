# =============================================================================
#  Comunità fungina trasportata da Corythucha arcuata su Quercus robur e Q. cerris
#  Script di analisi: copertura del campionamento, dispersione, ricchezza,
#  composizione e taxa indicatori (stagione su farnia; ospite in estate)
# =============================================================================
#
#  UNITÀ SPERIMENTALE = PIANTA. Le piastre sono sotto-campioni (pseudo-repliche):
#  le righe estive (4 piastre/pianta) vengono aggregate a livello di pianta, così
#  da renderle confrontabili con l'autunno (già aggregato in una riga/pianta).
#
#  Formato atteso del CSV (formato "largo"):
#    - una riga per piastra (estate) o per pianta (autunno);
#    - colonne descrittive (vedi CONFIGURAZIONE);
#    - una colonna per ciascun taxon fungino, con 0/1 (presenza) oppure conteggi.
# =============================================================================

## ---- 0. CONFIGURAZIONE (adattare ai nomi reali delle colonne) ---------------
file_dati   <- "G:/Il mio Drive/UNIVERSITA'/LAVORO/DATI_GIULIA_TIZIANA/AnalisiDati_Corythucha.xlsx"  # accetta anche .csv
# Nomi delle colonne del file AnalisiDati_Corythucha.csv
col_pianta  <- "N_pianta"   # identificativo della pianta (F1..F10, C1..C10)
col_ospite  <- "Host_species"
col_stag    <- "Season"
col_piastra <- "Plate"      # PL_1..PL_4 in estate; PL1_PL4 in autunno (4 piastre aggregate)
# Colonne numeriche che NON sono taxa. La colonna finale senza nome del CSV è
# una somma di Excel (Replicate + N_insetti + taxa) e viene esclusa a parte.
col_escludi <- c("N_insetti")
lev_robur   <- "Quercus_robur"
lev_cerris  <- "Quercus_cerris"
lev_estate  <- "Summer"
lev_autunno <- "Fall"
# Le stesse 10 farnie sono state campionate sia in estate sia in autunno?
# TRUE = disegno appaiato (misure ripetute); FALSE = piante diverse.
farnie_appaiate <- TRUE    # confermato: stesse 10 farnie in estate e in autunno
n_perm   <- 9999   # permutazioni
min_occ  <- 3      # n. minimo di piante in cui un taxon deve comparire per i test univariati
set.seed(2026)

## ---- 1. PACCHETTI ------------------------------------------------------------
# install.packages(c("readxl", "vegan","iNEXT","mvabund","indicspecies","logistf","glmmTMB"))
library(vegan)
has <- function(p) requireNamespace(p, quietly = TRUE)

## ---- 2. LETTURA E PREPARAZIONE -----------------------------------------------
leggi_csv <- function(f) {
  if (grepl("\\.xlsx?$", f, ignore.case = TRUE))
    return(as.data.frame(readxl::read_excel(f, sheet = 1)))
  prima <- readLines(f, n = 1, warn = FALSE)
  if (lengths(regmatches(prima, gregexpr(";", prima))) >
      lengths(regmatches(prima, gregexpr(",", prima))))
    read.csv2(f, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  else
    read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
}
if (!exists("dati_grezzi")) dati_grezzi <- leggi_csv(file_dati)

meta_cols <- intersect(c(col_pianta, col_ospite, col_stag, col_piastra), names(dati_grezzi))
# Rimuove colonne senza nome (es. ";" finale nel CSV) e quelle da escludere
dati_grezzi <- dati_grezzi[, nzchar(names(dati_grezzi)) & !grepl("^(V|X|\\.\\.\\.)[0-9]+$", names(dati_grezzi))]
taxa_cols <- setdiff(names(dati_grezzi), c(meta_cols, col_escludi))
taxa_cols <- taxa_cols[vapply(dati_grezzi[taxa_cols], is.numeric, logical(1))]
stopifnot(length(taxa_cols) > 0)
non_binari <- taxa_cols[!vapply(dati_grezzi[taxa_cols], function(x) all(x %in% c(0, 1)), logical(1))]
if (length(non_binari)) warning("Colonne non 0/1 trattate come taxa: ", paste(non_binari, collapse = ", "))
cat("Taxa analizzati (", length(taxa_cols), "):", paste(taxa_cols, collapse = ", "), "\n")
dati_grezzi[taxa_cols][is.na(dati_grezzi[taxa_cols])] <- 0

# Gruppo = combinazione ospite × stagione (RE = robur estate, RA = robur autunno, CE = cerris estate)
dati_grezzi$Gruppo <- interaction(dati_grezzi[[col_ospite]], dati_grezzi[[col_stag]],
                                  sep = " | ", drop = TRUE)

# Aggregazione a livello di pianta: un taxon è presente nella pianta se compare
# in almeno una delle piastre (equivale a quanto avvenuto di fatto in autunno).
chiave <- paste(dati_grezzi[[col_ospite]], dati_grezzi[[col_stag]], dati_grezzi[[col_pianta]], sep = "__")
abb_pianta <- rowsum(as.matrix(dati_grezzi[taxa_cols]), chiave, reorder = FALSE)
pa_pianta  <- (abb_pianta > 0) * 1
meta_pianta <- dati_grezzi[!duplicated(chiave), c(col_pianta, col_ospite, col_stag, "Gruppo")]
rownames(meta_pianta) <- unique(chiave)
meta_pianta <- meta_pianta[rownames(pa_pianta), ]
meta_pianta$Ricchezza <- rowSums(pa_pianta)

cat("\nPiante per gruppo:\n");      print(table(meta_pianta$Gruppo))
cat("\nRighe (piastre) per gruppo:\n"); print(table(dati_grezzi$Gruppo))

## ---- 3. DESCRIZIONE: QUALI FUNGHI TRASPORTA L'INSETTO ------------------------
# Frequenza di occorrenza (% di piante in cui il taxon è presente) per gruppo
freq_occ <- t(sapply(split(as.data.frame(pa_pianta), meta_pianta$Gruppo), colMeans)) * 100
freq_occ <- t(round(freq_occ, 1))
freq_occ <- freq_occ[order(-rowMeans(freq_occ)), , drop = FALSE]
cat("\nFrequenza di occorrenza (% piante):\n"); print(freq_occ)
write.csv2(freq_occ, "tab_frequenza_occorrenza.csv")

## ---- 4. SUFFICIENZA DEL CAMPIONAMENTO -----------------------------------------
# 4a. Statistiche di incidenza per gruppo: Sobs, singletons Q1 (taxa in 1 sola pianta),
#     doubletons Q2, Chao2 e copertura del campione (Chao & Jost 2012, dati di incidenza).
copertura_incidenza <- function(pa) {
  T  <- nrow(pa); y <- colSums(pa); y <- y[y > 0]
  U  <- sum(y); Q1 <- sum(y == 1); Q2 <- sum(y == 2)
  cov <- if (Q2 > 0) 1 - Q1 / U * ((T - 1) * Q1 / ((T - 1) * Q1 + 2 * Q2))
         else        1 - Q1 / U * ((T - 1) * (Q1 - 1) / ((T - 1) * (Q1 - 1) + 2))
  chao2 <- if (Q2 > 0) length(y) + (T - 1) / T * Q1^2 / (2 * Q2)
           else        length(y) + (T - 1) / T * Q1 * (Q1 - 1) / 2
  c(Unita = T, Sobs = length(y), Q1 = Q1, Q2 = Q2,
    Chao2 = round(chao2, 1), Completezza = round(length(y) / chao2, 2),
    Copertura = round(cov, 3))
}
suff <- t(sapply(split(as.data.frame(pa_pianta), meta_pianta$Gruppo), copertura_incidenza))
cat("\nSufficienza del campionamento (unità = piante):\n"); print(suff)
write.csv2(suff, "tab_sufficienza_campionamento.csv")

# Solo per l'estate: copertura con unità = piastra (dettaglio intra-pianta)
estate <- dati_grezzi[[col_stag]] == lev_estate
pa_piastra <- (as.matrix(dati_grezzi[estate, taxa_cols]) > 0) * 1
suff_piastra <- t(sapply(split(as.data.frame(pa_piastra), droplevels(dati_grezzi$Gruppo[estate])),
                         copertura_incidenza))
cat("\nEstate, unità = piastra (solo descrittivo):\n"); print(suff_piastra)

# Estate: le 4 piastre bastano a descrivere la comunità di una pianta?
# Per ogni pianta: taxa osservati, Chao2 con 4 piastre e taxa attesi con 1, 2, 3 piastre
suff_pianta <- do.call(rbind, lapply(split(as.data.frame(pa_piastra),
                       paste(dati_grezzi[[col_ospite]], dati_grezzi[[col_pianta]])[estate]), function(m) {
  y <- colSums(m); T <- nrow(m)
  rk <- sapply(1:(T - 1), function(k) sum(1 - choose(T - y[y > 0], k) / choose(T, k)))
  c(copertura_incidenza(m)[c("Sobs", "Q1", "Q2", "Chao2", "Completezza")],
    setNames(round(rk / sum(y > 0) * 100), paste0("%taxa_con_", 1:(T - 1), "_piastre")))
}))
cat("\nEstate, completezza entro pianta (unità = piastra):\n"); print(suff_pianta)
cat("Media:\n"); print(round(colMeans(suff_pianta), 2))
write.csv2(suff_pianta, "tab_sufficienza_piastre_per_pianta.csv")

# 4b. Curve di accumulo (vegan) per gruppo
pdf("fig_curve_accumulo.pdf", width = 7, height = 5)
gr <- levels(meta_pianta$Gruppo); col_gr <- c("#1b9e77", "#d95f02", "#7570b3")[seq_along(gr)]
for (i in seq_along(gr)) {
  sa <- specaccum(pa_pianta[meta_pianta$Gruppo == gr[i], , drop = FALSE], method = "exact")
  plot(sa, add = i > 1, col = col_gr[i], ci.type = "polygon",
       ci.col = adjustcolor(col_gr[i], 0.2), ci.lty = 0, lwd = 2,
       ylim = c(0, max(suff[, "Chao2"])), xlab = "Piante", ylab = "Taxa cumulati")
}
legend("bottomright", gr, col = col_gr, lwd = 2, bty = "n")
dev.off()

# 4c. Rarefazione/estrapolazione basata sulla copertura (iNEXT, se installato)
if (has("iNEXT")) {
  inc_freq <- lapply(split(as.data.frame(pa_pianta), meta_pianta$Gruppo),
                     function(m) c(nrow(m), colSums(m)))
  out_inext <- iNEXT::iNEXT(inc_freq, q = c(0, 1, 2), datatype = "incidence_freq", endpoint = 20)  # estrapolazione affidabile fino a 2 × n (Chao et al. 2014)
  print(out_inext$DataInfo)
  # Numeri di Hill (q=0 ricchezza, q=1 Shannon esponenziale, q=2 Simpson inverso)
  # confrontati a copertura comune (standardizzazione corretta fra gruppi)
  print(iNEXT::estimateD(inc_freq, q = c(0, 1, 2), datatype = "incidence_freq", base = "coverage"))
  pdf("fig_iNEXT.pdf", width = 8, height = 5)
  print(iNEXT::ggiNEXT(out_inext, type = 1)); print(iNEXT::ggiNEXT(out_inext, type = 3))
  dev.off()
}

## ---- 5. DISPERSIONE DEI DATI ---------------------------------------------------
# 5a. Ricchezza per pianta: rapporto varianza/media (Poisson atteso ≈ 1)
disp_ric <- aggregate(Ricchezza ~ Gruppo, meta_pianta,
                      function(x) c(n = length(x), media = mean(x), varianza = var(x),
                                    var_su_media = var(x) / mean(x)))
cat("\nDispersione della ricchezza:\n"); print(do.call(data.frame, disp_ric))
# < 1 sottodispersione (tipica di conteggi limitati superiormente), > 1 sovradispersione.

# 5b. Dispersione multivariata (omogeneità delle dispersioni, PERMDISP)
d_jac <- vegdist(pa_pianta, method = "jaccard", binary = TRUE)
bd <- betadisper(d_jac, meta_pianta$Gruppo, type = "centroid")
cat("\nPERMDISP (distanza dal centroide, Jaccard):\n")
print(permutest(bd, permutations = n_perm, pairwise = TRUE))

# 5c. Relazione media-varianza dei taxa (motiva l'uso di GLM multivariati)
if (has("mvabund")) {
  pdf("fig_media_varianza.pdf"); mvabund::meanvar.plot(mvabund::mvabund(abb_pianta) ~ meta_pianta$Gruppo); dev.off()
}

## ---- Funzioni di supporto per i confronti a due gruppi -------------------------
confronto_ricchezza <- function(md, fattore, appaiato = FALSE) {
  cat("\n--- Ricchezza ~", fattore, "---\n")
  print(tapply(md$Ricchezza, md[[fattore]], summary))
  r  <- md$Ricchezza; g <- droplevels(factor(md[[fattore]]))
  fm <- glm(r ~ g, family = poisson)
  phi <- sum(residuals(fm, type = "pearson")^2) / df.residual(fm)
  cat("Indice di dispersione (Pearson chi2/df) =", round(phi, 2), "\n")
  if (appaiato) {
    md <- md[order(md[[col_pianta]]), ]
    w <- split(md$Ricchezza, droplevels(factor(md[[fattore]])))
    print(wilcox.test(w[[1]], w[[2]], paired = TRUE, exact = FALSE))
    if (has("glmmTMB")) {
      md$g <- droplevels(factor(md[[fattore]])); md$pl <- factor(md[[col_pianta]])
      print(summary(glmmTMB::glmmTMB(Ricchezza ~ g + (1 | pl), md,
                                     family = if (phi < 1) glmmTMB::genpois() else poisson)))
    }
  } else {
    if (phi < 1 && has("glmmTMB")) {
      # Sottodispersione: Poisson generalizzata (gestisce varianza < media)
      print(summary(glmmTMB::glmmTMB(r ~ g, data = data.frame(r, g), family = glmmTMB::genpois())))
    } else if (phi > 1.5) {
      print(summary(glm(r ~ g, family = quasipoisson)))
    } else print(summary(fm))
    print(wilcox.test(r ~ g, exact = FALSE))       # alternativa non parametrica
  }
}

confronto_composizione <- function(pa, md, fattore, appaiato = FALSE) {
  cat("\n--- Composizione ~", fattore, "---\n")
  keep <- rowSums(pa) > 0
  pa <- pa[keep, colSums(pa) > 0, drop = FALSE]; md <- md[keep, ]
  g  <- droplevels(factor(md[[fattore]]))
  d  <- vegdist(pa, method = "jaccard", binary = TRUE)
  if (appaiato) {
    # Disegno a misure ripetute: la pianta entra come blocco e le permutazioni
    # avvengono solo entro pianta (si scambiano le etichette di stagione).
    pl <- factor(md[[col_pianta]])
    print(adonis2(d ~ pl + g, by = "terms",
                  permutations = how(nperm = n_perm, blocks = pl)))
  } else print(adonis2(d ~ g, permutations = how(nperm = n_perm)))
  print(permutest(betadisper(d, g), permutations = n_perm))   # PERMANOVA valida se p non significativo
  nm <- metaMDS(pa, distance = "jaccard", binary = TRUE, k = 2, trymax = 100, trace = 0)
  pdf(paste0("fig_NMDS_", fattore, ".pdf"))
  plot(nm, type = "n", main = paste("NMDS Jaccard – stress =", round(nm$stress, 3)))
  points(nm, display = "sites", pch = 19, col = as.integer(g) + 1)
  ordiellipse(nm, g, kind = "sd", label = TRUE)
  dev.off()
  if (has("mvabund")) {
    Y  <- mvabund::mvabund(pa)
    mg <- mvabund::manyglm(Y ~ g, family = "binomial")
    # Test score con PIT-trap: robusti anche in caso di separazione.
    # Se appaiato, si ricampionano intere piante (block) per rispettare la dipendenza.
    blk <- if (appaiato) factor(md[[col_pianta]]) else NULL
    print(mvabund::anova.manyglm(mg, p.uni = "adjusted", test = "score",
                                 resamp = "pit.trap", nBoot = 999, block = blk))
  }
}

# Taxa associati a un livello: Fisher esatto (o McNemar se appaiato),
# regressione logistica di Firth (stima finita anche con separazione), IndVal.
taxa_indicatori <- function(pa, md, fattore, appaiato = FALSE) {
  cat("\n--- Taxa indicatori ~", fattore, "---\n")
  g  <- droplevels(factor(md[[fattore]]))
  pa <- pa[, colSums(pa) >= min_occ, drop = FALSE]
  res <- do.call(rbind, lapply(colnames(pa), function(tx) {
    y   <- pa[, tx]
    tab <- table(factor(y, 0:1), g)
    sep <- if (any(tab[2, ] == 0) || any(tab[1, ] == 0)) {
      if (all(diag(tab) == 0) || all(diag(tab[2:1, ]) == 0)) "completa" else "quasi-completa"
    } else "no"
    discord <- NA_character_
    p_test <- if (appaiato) {
      # McNemar esatto: test binomiale sulle sole piante discordanti
      o <- order(md[[col_pianta]]); yy <- split(y[o], g[o])
      n10 <- sum(yy[[1]] == 1 & yy[[2]] == 0); n01 <- sum(yy[[1]] == 0 & yy[[2]] == 1)
      discord <- sprintf("solo %s: %d | solo %s: %d", levels(g)[1], n10, levels(g)[2], n01)
      if (n10 + n01 == 0) 1 else binom.test(n10, n10 + n01)$p.value
    } else fisher.test(tab)$p.value
    or_firth <- p_firth <- NA
    if (has("logistf") && !appaiato) {
      ff <- logistf::logistf(y ~ g, data = data.frame(y = y, g = g))
      or_firth <- exp(coef(ff)[2]); p_firth <- ff$prob[2]
    }
    data.frame(Taxon = tx,
               t(setNames(round(100 * tab[2, ] / colSums(tab), 0), paste0("Freq%_", levels(g)))),
               Separazione = sep, Piante_discordanti = discord, p_esatto = p_test,
               OR_Firth = round(or_firth, 2), p_Firth = p_firth, check.names = FALSE)
  }))
  res$q_esatto_FDR <- p.adjust(res$p_esatto, "BH")
  if (!all(is.na(res$p_Firth))) res$q_Firth_FDR <- p.adjust(res$p_Firth, "BH")
  res <- res[order(res$p_esatto), ]
  print(res, row.names = FALSE)
  if (has("indicspecies")) {
    # Se appaiato, le permutazioni avvengono solo entro pianta
    ctrl_iv <- if (appaiato) how(nperm = n_perm, blocks = factor(md[[col_pianta]])) else how(nperm = n_perm)
    iv <- indicspecies::multipatt(as.data.frame(pa), g, func = "IndVal.g", control = ctrl_iv)
    summary(iv, indvalcomp = TRUE)
  }
  invisible(res)
}

## ---- 6. EFFETTO DELLA STAGIONE (solo Q. robur) ----------------------------------
sel_s <- meta_pianta[[col_ospite]] == lev_robur
md_s  <- meta_pianta[sel_s, ]; pa_s <- pa_pianta[sel_s, , drop = FALSE]
if (farnie_appaiate) {
  # Controllo: ogni farnia deve comparire una volta per stagione con lo stesso codice
  tp <- table(md_s[[col_pianta]], md_s[[col_stag]])
  if (any(tp != 1)) { print(tp); stop("Codici pianta non appaiati fra estate e autunno") }
}
confronto_ricchezza(md_s, col_stag, appaiato = farnie_appaiate)
confronto_composizione(pa_s, md_s, col_stag, appaiato = farnie_appaiate)
ind_stag <- taxa_indicatori(pa_s, md_s, col_stag, appaiato = farnie_appaiate)
write.csv2(ind_stag, "tab_taxa_stagione.csv", row.names = FALSE)

# 6b. Scomposizione della beta-diversità per ciascuna farnia (Baselga 2012):
#     dissimilarità di Jaccard = sostituzione di taxa (turnover) + perdita/
#     acquisizione di taxa (nestedness). Indica se la comunità autunnale è un
#     sottoinsieme/ampliamento di quella estiva o se i taxa vengono sostituiti.
if (farnie_appaiate) {
  o  <- order(md_s[[col_pianta]])
  st <- factor(md_s[[col_stag]][o]); X <- pa_s[o, , drop = FALSE]
  E  <- X[st == lev_estate, , drop = FALSE]; A <- X[st == lev_autunno, , drop = FALSE]
  a  <- rowSums(E & A); b <- rowSums(E & !A); c <- rowSums(!E & A)
  jac <- (b + c) / (a + b + c)
  jtu <- 2 * pmin(b, c) / (a + 2 * pmin(b, c))
  beta_pl <- data.frame(Pianta = md_s[[col_pianta]][o][st == lev_estate],
                        Condivisi = a, Solo_estate = b, Solo_autunno = c,
                        Jaccard = round(jac, 2), Turnover = round(jtu, 2),
                        Nestedness = round(jac - jtu, 2))
  cat("\nScomposizione beta-diversità estate vs autunno, per farnia:\n")
  print(beta_pl, row.names = FALSE)
  cat("Medie: Jaccard =", round(mean(jac), 2), "| Turnover =", round(mean(jtu), 2),
      "| Nestedness =", round(mean(jac - jtu), 2), "\n")
  write.csv2(beta_pl, "tab_beta_stagione_per_pianta.csv", row.names = FALSE)
}

## ---- 7. EFFETTO DELL'OSPITE (solo estate: robur vs cerris) ---------------------
sel_o <- meta_pianta[[col_stag]] == lev_estate
md_o  <- meta_pianta[sel_o, ]; pa_o <- pa_pianta[sel_o, , drop = FALSE]
confronto_ricchezza(md_o, col_ospite)
confronto_composizione(pa_o, md_o, col_ospite)
ind_osp <- taxa_indicatori(pa_o, md_o, col_ospite)
write.csv2(ind_osp, "tab_taxa_ospite.csv", row.names = FALSE)

# 7b. Analisi di sensibilità che sfrutta le piastre estive: per ciascun taxon,
#     n. di piastre positive su 4 per pianta (binomiale con effetto casuale pianta).
if (has("glmmTMB")) {
  de <- dati_grezzi[estate, ]
  de$pl <- factor(paste(de[[col_ospite]], de[[col_pianta]]))
  de$osp <- droplevels(factor(de[[col_ospite]]))
  cat("\n--- Sensibilità: livello piastra, effetto casuale pianta ---\n")
  for (tx in colnames(pa_o)[colSums(pa_o) >= min_occ]) {
    de$y <- as.integer(de[[tx]] > 0)
    m <- try(glmmTMB::glmmTMB(y ~ osp + (1 | pl), de, family = binomial), silent = TRUE)
    if (!inherits(m, "try-error")) {
      cf <- summary(m)$coefficients$cond
      cat(sprintf("%-30s logOR = %6.2f  SE = %6.2f  p = %.4f\n", tx, cf[2, 1], cf[2, 2], cf[2, 4]))
    }
  }
}

## ---- 8. MODELLO UNICO A TRE GRUPPI CON CONTRASTI PIANIFICATI (alternativa) ----
# Gruppo con 3 livelli (robur-estate, robur-autunno, cerris-estate): l'interazione
# ospite × stagione NON è stimabile (cella cerris-autunno vuota), ma i due contrasti
# di interesse lo sono, con un'unica stima della varianza residua.
cat("\n--- PERMANOVA a tre gruppi ---\n")
print(adonis2(d_jac ~ Gruppo, data = meta_pianta, permutations = n_perm))
cat("\nAnalisi completata. Tabelle e figure salvate in:", getwd(), "\n")
