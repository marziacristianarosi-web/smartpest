# =============================================================================
# esegui_tutto.R - Esegue l'intera procedura nell'ordine dei passaggi del report.
# Uso: aprire R/RStudio, impostare come cartella di lavoro la cartella "script"
#      (Session > Set Working Directory > To Source File Location) ed eseguire:
#      source("esegui_tutto.R")
# Tempi indicativi su 4 processori: passaggi 1-7 pochi minuti; 8 circa 10 min;
# 9 circa 15 min; 10 circa 10-20 min (9999 permutazioni della domanda 4).
# Ogni passaggio scrive un registro in script/log.
# =============================================================================
source("00_impostazioni.R")
passi <- c("01_struttura_dati.R", "02_curve_accumulo.R", "03_chao2_copertura.R", "04_rarefazione_estrapolazione.R",
           "05_distribuzione.R", "06_scelta_distribuzione.R", "07_assunzioni.R", "08_errore_primo_tipo.R",
           "09_potenza.R", "10_test_ipotesi.R", "11_figure_risultati.R")
for (f in passi) {
  cat("\n==========", f, "==========\n")
  log <- file(file.path(DIR_LOG, sub("\\.R$", ".log", f)), open = "wt")
  sink(log, split = TRUE); t0 <- Sys.time()
  tryCatch(source(f, local = new.env(parent = globalenv())), finally = { sink(); close(log) })
  cat(sprintf("  completato in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
writeLines(capture.output(sessionInfo()), file.path(DIR_LOG, "sessionInfo.txt"))
