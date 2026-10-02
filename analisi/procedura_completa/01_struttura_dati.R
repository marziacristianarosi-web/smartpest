# =============================================================================
# 01_struttura_dati.R - Passaggio 1: struttura del disegno e aggregazione per pianta
# Funzioni: table(), rowsum() (base R)
# =============================================================================
if (!exists("PIANTE")) source("00_impostazioni.R")   # cartella di lavoro = cartella "script"

# In estate ogni riga è una piastra (3 insetti); in autunno le 4 piastre della pianta sono registrate in un'unica riga (PL1_PL4)
disegno <- aggregate(cbind(Righe = 1, Insetti = N_insetti) ~ Season + Host_species, grezzi, sum)
disegno$Piante <- aggregate(N_pianta ~ Season + Host_species, grezzi, function(x) length(unique(x)))$N_pianta
disegno$Insetti_per_pianta <- disegno$Insetti / disegno$Piante
print(disegno)
salva_tab(disegno, "T01_disegno")

matrice <- data.frame(Stagione = PIANTE$Season, Ospite = PIANTE$Host_species, Pianta = PIANTE$N_pianta,
                      Ricchezza = PIANTE$R, PA, check.names = FALSE)
salva_tab(matrice, "T02_matrice_piante_taxa")
