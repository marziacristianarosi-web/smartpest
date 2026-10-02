# =============================================================================
# 03_chao2_copertura.R - Passaggio 3: ricchezza stimata (Chao2) e copertura del campione
# Funzioni: iNEXT::ChaoRichness() (Chao2 corretto per il bias, IC log-normale),
#           iNEXT::DataInfo() (copertura del campione, Chao & Jost 2012)
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")

inc <- lapply(levels(PIANTE$gruppo), function(g) { m <- PA[PIANTE$gruppo == g, ]; c(nrow(m), colSums(m)) })
names(inc) <- ETI[levels(PIANTE$gruppo)]
info <- DataInfo(inc, datatype = "incidence_freq")
chao <- do.call(rbind, lapply(inc, ChaoRichness, datatype = "incidence_freq"))
tab <- data.frame(Gruppo = names(inc), Piante = info$T, Taxa_osservati = info$S.obs,
                  Q1_unici = info$Q1, Q2_duplicati = info$Q2,
                  Chao2 = round(chao$Estimator, 1), ES = round(chao$Est_s.e., 2),
                  IC95_inf = round(chao$`95% Lower`, 1), IC95_sup = round(chao$`95% Upper`, 1),
                  Completezza = round(info$S.obs / chao$Estimator, 2), Copertura = round(info$SC, 3))
print(tab, row.names = FALSE); salva_tab(tab, "T05_chao2_copertura")
