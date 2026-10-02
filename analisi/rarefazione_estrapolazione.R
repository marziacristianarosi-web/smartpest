# Rarefazione/estrapolazione con numeri di Hill (iNEXT), confronto a copertura comune e figura
# Richiede: install.packages(c("readxl","iNEXT","ggplot2","patchwork"))
file_dati <- "G:/Il mio Drive/UNIVERSITA'/LAVORO/DATI_GIULIA_TIZIANA/AnalisiDati_Corythucha.xlsx"
dir.create("out", showWarnings = FALSE)
suppressMessages({library(readxl); library(iNEXT)}); set.seed(2026)
d <- as.data.frame(read_excel(file_dati)); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
grp <- sub(" [^ ]+$","", rownames(pa))
lab <- c("Fall Quercus_robur"="Farnia - autunno","Summer Quercus_robur"="Farnia - estate","Summer Quercus_cerris"="Cerro - estate")
inc <- lapply(names(lab), function(g){ m <- pa[grp==g,]; c(nrow(m), colSums(m)) }); names(inc) <- lab
out <- iNEXT(inc, q=c(0,1,2), datatype="incidence_freq", endpoint=20, nboot=200)
cat("== DataInfo ==\n"); print(out$DataInfo)
cat("\n== Stime asintotiche ==\n"); print(out$AsyEst)
cat("\n== Valori a T=10 e T=20 (q=0,1,2) ==\n")
s <- out$iNextEst$size_based; print(s[s$t %in% c(5,10,20), c("Assemblage","t","Method","Order.q","qD","qD.LCL","qD.UCL","SC")], row.names=FALSE)
Cmin <- min(out$DataInfo$SC); Cmax <- min(sapply(names(inc), function(n) max(s$SC[s$Assemblage==n])))
cat("\nCopertura comune (min osservata):", Cmin, " | max confrontabile (estrap. 2T):", round(Cmax,3), "\n")
e1 <- estimateD(inc, q=c(0,1,2), datatype="incidence_freq", base="coverage", level=round(Cmin,3), nboot=200)
print(e1, row.names=FALSE)
e2 <- estimateD(inc, q=c(0,1,2), datatype="incidence_freq", base="coverage", level=round(Cmax,3), nboot=200)
print(e2, row.names=FALSE)
saveRDS(out, "out/inext.rds"); write.csv2(rbind(cbind(Livello="Cmin", e1), cbind(Livello="Cmax", e2)), "out/tab_hill_copertura_comune.csv", row.names=FALSE)

## ---- Figura ----
suppressMessages({library(ggplot2); library(patchwork)})
s <- out$iNextEst$size_based
col <- c("Farnia - autunno"="#eb6834", "Farnia - estate"="#2a78d6", "Cerro - estate"="#1baf7a")
s$Gruppo <- factor(s$Assemblage, levels=names(col)); s0 <- s[s$Order.q==0,]
oss <- s0[s0$Method=="Observed",]; fin <- s0[s0$t==20,]
tema <- theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
  panel.grid.major.y=element_line(colour="#e6e5e0", linewidth=0.3), axis.line=element_line(colour="#8a8984", linewidth=0.3),
  axis.text=element_text(colour="#52514e"), text=element_text(colour="#0b0b0b"), legend.position="none",
  plot.title=element_text(face="bold", size=12), plot.subtitle=element_text(colour="#52514e", size=10),
  plot.background=element_rect(fill="#fcfcfb", colour=NA))
zone <- function(ytop) list(
  annotate("rect", xmin=10, xmax=20.5, ymin=-Inf, ymax=Inf, fill="#efeee9", alpha=0.6),
  annotate("text", x=5.5, y=ytop, label="DATI REALI\n(da 1 a 10 piante campionate)", size=3.3, colour="#52514e", fontface="bold"),
  annotate("text", x=15.25, y=ytop, label="PREVISIONE\n(se si campionassero 11-20 piante)", size=3.3, colour="#52514e", fontface="bold"))
lin <- function(y) list(
  geom_ribbon(aes(ymin=.data[[paste0(y,".LCL")]], ymax=.data[[paste0(y,".UCL")]]), alpha=0.13, colour=NA),
  geom_line(data=s0[s0$t<=10,], linewidth=1), geom_line(data=s0[s0$t>=10,], linewidth=1, linetype="22"),
  geom_point(data=oss, size=3.2, shape=21, stroke=0.8, colour="#fcfcfb"))
pa <- ggplot(s0, aes(t, qD, colour=Gruppo, fill=Gruppo)) + zone(21) + lin("qD") +
  geom_text(data=fin, aes(label=Gruppo), hjust=0, nudge_x=0.3, size=3.4, colour="#0b0b0b") +
  scale_colour_manual(values=col) + scale_fill_manual(values=col) +
  scale_x_continuous(breaks=c(1,5,10,15,20), limits=c(1,25.5), expand=c(0.01,0)) +
  scale_y_continuous(breaks=seq(0,22,2), expand=c(0,0)) + coord_cartesian(ylim=c(0,22.5)) +
  labs(title="A. Quanti taxa si trovano campionando un certo numero di piante?",
       subtitle="Punto = valore reale con 10 piante; banda = intervallo di confidenza al 95%",
       x="Numero di piante campionate", y="Numero di taxa trovati") + tema
finB <- data.frame(t=20, SC=c(1, fin$SC[fin$Gruppo=="Cerro - estate"]), lab=c("Farnia (autunno ed estate)", "Cerro - estate"))
pb <- ggplot(s0, aes(t, SC, colour=Gruppo, fill=Gruppo)) + zone(0.36) + lin("SC") +
  geom_text(data=finB, aes(t, SC, label=lab), inherit.aes=FALSE, hjust=0, vjust=c(-0.4, 1.2), nudge_x=0.3, size=3.4, colour="#0b0b0b", lineheight=0.9) +
  scale_colour_manual(values=col) + scale_fill_manual(values=col) +
  scale_x_continuous(breaks=c(1,5,10,15,20), limits=c(1,25.5), expand=c(0.01,0)) +
  scale_y_continuous(breaks=seq(0.3,1,0.1), labels=function(x) format(x, decimal.mark=","), expand=c(0,0)) +
  coord_cartesian(ylim=c(0.3,1.035)) +
  labs(title="B. Quanto è completo il campionamento?",
       subtitle="Copertura = 1 significa che un'altra pianta non porterebbe taxa nuovi",
       x="Numero di piante campionate", y="Copertura del campione") + tema
p <- pa / pb
ggsave("out/fig_rarefazione_estrapolazione.png", p, width=9, height=10, dpi=300)
ggsave("out/fig_rarefazione_estrapolazione.pdf", p, width=9, height=10)
