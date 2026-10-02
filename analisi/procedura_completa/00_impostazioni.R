# =============================================================================
# 00_impostazioni.R
# Comunità fungina associata a Corythucha arcuata su Quercus robur e Q. cerris
# Impostazioni comuni: percorsi, pacchetti, lettura dei dati, aggregazione per
# pianta, palette e funzioni di servizio. Viene richiamato da tutti gli altri script.
# =============================================================================

# ---- 1. Percorsi -------------------------------------------------------------
# Cartella che contiene il file dei dati (il file originale non viene mai modificato)
# e cartella dei risultati accanto ad esso. Modificare solo CARTELLA se cambia posizione.
# (Le variabili d'ambiente CORY_DATI e CORY_BASE, se impostate, hanno la precedenza.)
CARTELLA    <- Sys.getenv("CORY_DATI", "G:/Il mio Drive/UNIVERSITA'/LAVORO/DATI_GIULIA_TIZIANA")
FILE_DATI   <- file.path(CARTELLA, "AnalisiDati_Corythucha.xlsx")
FOGLIO_DATI <- "2_10_2026_AnalisiDati_Corythuch"      # foglio con i dati corretti
BASE    <- Sys.getenv("CORY_BASE", file.path(CARTELLA, "selezione_test_AnalisiDati_Corythucha"))
DIR_FIG <- file.path(BASE, "grafici")
DIR_TAB <- file.path(BASE, "tabelle")
DIR_LOG <- file.path(BASE, "script", "log")
for (d in c(DIR_FIG, DIR_TAB, DIR_LOG)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- 2. Pacchetti ------------------------------------------------------------
# Installazione (una sola volta):
# install.packages(c("readxl", "writexl", "vegan", "permute", "iNEXT", "glmmTMB",
#                    "lme4", "ggplot2", "patchwork", "ragg"))
suppressMessages({
  library(readxl)     # lettura del file Excel
  library(writexl)    # scrittura delle tabelle in Excel
  library(vegan)      # specaccum, vegdist, adonis2, betadisper, permutest, wcmdscale
  library(permute)    # permutazioni ristrette
  library(iNEXT)      # rarefazione/estrapolazione, Chao2, copertura (numeri di Hill)
  library(glmmTMB)    # GLM/GLMM con distribuzione di Poisson generalizzata
  library(ggplot2)    # grafici
  library(patchwork)  # composizione dei pannelli
})
NCORE <- max(1, parallel::detectCores() - 1)   # processori usati per le permutazioni

# ---- 3. Lettura dei dati -----------------------------------------------------
grezzi <- as.data.frame(read_excel(FILE_DATI, sheet = FOGLIO_DATI), check.names = FALSE)
META <- c("Season", "Host_species", "N_pianta ORIG", "N_pianta", "Piastra", "Replicate", "N_insetti")
TAXA <- setdiff(names(grezzi), META)                     # 14 colonne di presenza/assenza
stopifnot(all(unlist(grezzi[TAXA]) %in% c(0, 1)))

# ---- 4. Aggregazione per pianta (unità campionaria) --------------------------
# Un taxon è presente nella pianta se isolato in almeno una piastra (incidenza).
chiave <- paste(grezzi$Season, grezzi$Host_species, grezzi$N_pianta, sep = "|")
PA <- (rowsum(as.matrix(grezzi[TAXA]), chiave, reorder = FALSE) > 0) * 1
PIANTE <- grezzi[!duplicated(chiave), c("Season", "Host_species", "N_pianta")]
rownames(PIANTE) <- unique(chiave); PIANTE <- PIANTE[rownames(PA), ]
PIANTE$R <- rowSums(PA)                                   # ricchezza in taxa per pianta
PIANTE$gruppo <- factor(ifelse(PIANTE$Host_species == "Quercus_cerris", "cerro_estate",
                        ifelse(PIANTE$Season == "Summer", "farnia_estate", "farnia_autunno")),
                        levels = c("farnia_estate", "farnia_autunno", "cerro_estate"))

# Domande 1-3: stagione, stesse 10 farnie (disegno appaiato)
iS <- which(PIANTE$Host_species == "Quercus_robur")
iS <- iS[order(PIANTE$N_pianta[iS], PIANTE$Season[iS])]
S <- PIANTE[iS, ]; S$pl <- factor(S$N_pianta); S$g <- factor(S$Season, levels = c("Summer", "Fall"))
PA_S <- PA[iS, ]
stopifnot(all(table(S$pl) == 2))

# Domande 4-6: specie ospite in estate (piante indipendenti)
iH <- which(PIANTE$Season == "Summer")
H <- PIANTE[iH, ]; H$g <- factor(H$Host_species, levels = c("Quercus_robur", "Quercus_cerris"))
PA_H <- PA[iH, ]

# Domande 4-6, versione per piastra: estate, 4 piastre (3 insetti ciascuna) per pianta.
# L'ospite è assegnato alla pianta: le unità indipendenti restano 10 + 10 piante; le piastre sono sottocampioni.
E <- grezzi[grezzi$Season == "Summer", c("Host_species", "N_pianta", "Replicate")]
E$g <- factor(E$Host_species, levels = c("Quercus_robur", "Quercus_cerris")); E$pl <- factor(E$N_pianta)
PA_E <- as.matrix(grezzi[grezzi$Season == "Summer", TAXA]); E$R <- rowSums(PA_E)   # ricchezza per piastra
FREQ_H <- rowsum(PA_E, as.character(E$pl))[as.character(H$N_pianta), ]               # piastre positive per pianta (0-4)
# permutazione di piante intere (l'etichetta dell'ospite segue la pianta con tutte le sue piastre)
perm_piante <- function(idx_piante) {   # idx_piante: permutazione delle 20 piante di H
  nuovo <- setNames(as.character(H$g[idx_piante]), H$N_pianta); factor(nuovo[as.character(E$pl)], levels = levels(E$g)) }

# ---- 5. Aspetto grafico ------------------------------------------------------
COL <- c(farnia_estate = "#2a78d6", farnia_autunno = "#eb6834", cerro_estate = "#1baf7a")  # verificata per daltonismo
SHP <- c(farnia_estate = 16, farnia_autunno = 17, cerro_estate = 15)                       # codifica ridondante
ETI <- c(farnia_estate = "Q. robur, estate", farnia_autunno = "Q. robur, autunno", cerro_estate = "Q. cerris, estate")
ETI_IT <- c(farnia_estate = expression(italic("Q. robur")*", estate"), farnia_autunno = expression(italic("Q. robur")*", autunno"),
            cerro_estate = expression(italic("Q. cerris")*", estate"))   # etichette con il nome di specie in corsivo
FONT <- if ("Liberation Sans" %in% systemfonts::system_fonts()$family) "Liberation Sans" else "sans"
TEMA <- theme_classic(base_size = 9, base_family = FONT) +
  theme(axis.text = element_text(colour = "grey20"), axis.line = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3), legend.position = "none",
        plot.tag = element_text(face = "bold", size = 10), plot.title = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        legend.background = element_blank(), legend.key = element_blank(), legend.text.align = 0)

salva_fig <- function(p, nome, w, h) {   # PNG 600 dpi per Word, PDF vettoriale per la rivista
  ggsave(file.path(DIR_FIG, paste0(nome, ".png")), p, width = w, height = h, units = "mm", dpi = 600,
         bg = "white", device = ragg::agg_png)
  ggsave(file.path(DIR_FIG, paste0(nome, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
}
salva_tab <- function(x, nome) write_xlsx(as.data.frame(x), file.path(DIR_TAB, paste0(nome, ".xlsx")))
virgola <- function(x, d = 2) formatC(x, digits = d, format = "f", decimal.mark = ",")

# Rapporto di verosimiglianza tra modello completo e ridotto (Poisson generalizzata)
lrt_gp <- function(f1, f0, dd) {
  a <- try(suppressWarnings(glmmTMB(f1, data = dd, family = genpois())), silent = TRUE)
  b <- try(suppressWarnings(glmmTMB(f0, data = dd, family = genpois())), silent = TRUE)
  if (inherits(a, "try-error") || inherits(b, "try-error")) return(NA_real_)
  as.numeric(2 * (logLik(a) - logLik(b)))
}
# Applicazione in parallelo compatibile con Windows (cluster PSOCK)
in_parallelo <- function(X, FUN, esporta = character(0), ...) {
  cl <- parallel::makeCluster(NCORE); on.exit(parallel::stopCluster(cl))
  parallel::clusterEvalQ(cl, suppressMessages(library(glmmTMB)))
  parallel::clusterExport(cl, c("lrt_gp", esporta), envir = globalenv())
  parallel::parLapply(cl, X, FUN, ...)
}
# Tutte le 2^10 = 1024 permutazioni entro pianta (scambio o no delle due etichette stagionali di ciascuna pianta)
PERM_ENTRO <- local({
  combo <- as.matrix(expand.grid(rep(list(0:1), nlevels(S$pl))))
  t(apply(combo, 1, function(sw) { idx <- seq_len(nrow(S))
    for (i in which(sw == 1)) { k <- which(S$pl == levels(S$pl)[i]); idx[k] <- rev(idx[k]) }
    idx }))
})
cat(sprintf("Dati: %d righe, %d taxa, %d piante-stagione\n", nrow(grezzi), length(TAXA), nrow(PIANTE)))
