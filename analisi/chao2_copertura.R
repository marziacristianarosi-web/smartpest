# Chao2 (con IC 95% log-normale), completezza, copertura del campione e sforzo aggiuntivo
# Unità = pianta (piastre estive aggregate per pianta); gruppi analizzati separatamente
file_dati <- "G:/Il mio Drive/UNIVERSITA'/LAVORO/DATI_GIULIA_TIZIANA/AnalisiDati_Corythucha.xlsx"
dir.create("out", showWarnings = FALSE)
suppressMessages({library(readxl)})
d <- as.data.frame(read_excel(file_dati)); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
grp <- sub(" [A-Z][0-9]+$","", rownames(pa))
lab <- c("Fall Quercus_robur"="Farnia - autunno","Summer Quercus_robur"="Farnia - estate","Summer Quercus_cerris"="Cerro - estate")
stima <- function(m){
  T <- nrow(m); y <- colSums(m); y <- y[y>0]; S <- length(y); U <- sum(y); Q1 <- sum(y==1); Q2 <- sum(y==2)
  # Chao2 corretto per il bias (Chao 1987; Chao et al. 2009)
  f0 <- if (Q2>0) (T-1)/T*Q1^2/(2*Q2) else (T-1)/T*Q1*(Q1-1)/2
  ch <- S + f0
  # varianza (formula per Q2 = 0 / Q2 > 0, come in iNEXT/SpadeR) e IC log-normale
  if (Q1==0) { se <- 0; lo <- hi <- S } else {
    v <- if (Q2>0) { k <- (T-1)/T; Q2*(k*0.5*(Q1/Q2)^2 + k^2*(Q1/Q2)^3 + k^2*0.25*(Q1/Q2)^4) } else {
      k <- (T-1)/T; k*Q1*(Q1-1)/2 + k^2*Q1*(2*Q1-1)^2/4 - k^2*Q1^4/(4*ch) }
    se <- sqrt(v); K <- exp(1.96*sqrt(log(1+v/f0^2))); lo <- S + f0/K; hi <- S + f0*K }
  # copertura del campione (Chao & Jost 2012; Chao et al. 2014)
  A <- if (Q2>0) (T-1)*Q1/((T-1)*Q1+2*Q2) else if (Q1>0) (T-1)*(Q1-1)/((T-1)*(Q1-1)+2) else 1
  cov <- 1 - Q1/U*A
  # estrapolazione: taxa attesi e copertura con m piante in più (Chao et al. 2014)
  Sm <- function(mm) if (f0==0) S else S + f0*(1-(1-Q1/(T*f0+Q1))^mm)
  Cm <- function(mm) 1 - Q1/U*A^(mm+1)
  piu <- function(fun, target) { for (mm in 0:500) if (fun(mm) >= target) return(mm); NA }
  data.frame(Piante=T, Taxa=S, Q1=Q1, Q2=Q2, Chao2=round(ch,1), ES=round(se,1),
             IC95=sprintf("%.1f - %.1f", lo, hi), Completezza=round(S/ch,2),
             Copertura=round(cov,3), Deficit_copertura=round(1-cov,3),
             Piante_extra_90pc_Chao2=piu(Sm, 0.9*ch), Piante_extra_cop_0.99=piu(Cm, 0.99),
             Taxa_attesi_con_20_piante=round(Sm(10),1))
}
res <- do.call(rbind, lapply(names(lab), function(g) cbind(Gruppo=lab[g], stima(pa[grp==g,]))))
print(res, row.names=FALSE); write.csv2(res, "out/tab_chao2_copertura.csv", row.names=FALSE)
y <- colSums(pa[grp=="Summer Quercus_cerris",]); cat("\nCerro, taxa presenti in 1 sola pianta:", names(y)[y==1], "\n")
