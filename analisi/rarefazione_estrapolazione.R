# Rarefazione/estrapolazione con numeri di Hill (iNEXT), confronto a copertura comune e figura
# Richiede: install.packages(c("readxl","iNEXT","ggplot2","patchwork"))
file_dati <- "G:/Il mio Drive/UNIVERSITA'/LAVORO/DATI_GIULIA_TIZIANA/AnalisiDati_Corythucha.xlsx"
dir.create("out", showWarnings = FALSE)
suppressMessages({library(readxl); library(iNEXT)}); set.seed(2026)
d <- as.data.frame(read_excel(file_dati)); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
grp <- sub(" [A-Z][0-9]+$","", rownames(pa))
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
col <- c("Farnia - autunno"="#eb6834", "Farnia - estate"="#2a78d6", "Cerro - estate"="#1baf7a")
s$Gruppo <- factor(s$Assemblage, levels=names(col))
s$Tratto <- ifelse(s$Method=="Extrapolation", "Estrapolazione", "Rarefazione")
s0 <- s[s$Order.q==0,]
oss <- s0[s0$Method=="Observed",]
tema <- theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
  panel.grid.major.y=element_line(colour="#e6e5e0", linewidth=0.3), axis.line=element_line(colour="#8a8984", linewidth=0.3),
  axis.ticks=element_line(colour="#8a8984"), axis.text=element_text(colour="#52514e"), text=element_text(colour="#0b0b0b"),
  legend.position="top", legend.justification="left", plot.tag=element_text(face="bold"),
  plot.background=element_rect(fill="#fcfcfb", colour=NA))
pa <- ggplot(s0, aes(t, qD, colour=Gruppo, fill=Gruppo)) +
  geom_vline(xintercept=10, colour="#c3c2b7", linewidth=0.4) +
  geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL), alpha=0.12, colour=NA) +
  geom_line(data=s0[s0$t<=10,], linewidth=0.9) +
  geom_line(data=s0[s0$t>=10,], linewidth=0.9, linetype="22") +
  geom_point(data=oss, size=3, shape=21, stroke=0.8, colour="#fcfcfb") +
  scale_colour_manual(values=col) + scale_fill_manual(values=col) +
  scale_x_continuous(breaks=seq(0,20,2), expand=c(0.01,0)) + scale_y_continuous(breaks=seq(0,22,2), expand=c(0,0)) +
  coord_cartesian(ylim=c(0,22)) +
  labs(tag="A", x="Numero di piante", y="Numero di taxa (q = 0)", colour=NULL, fill=NULL) + tema
pb <- ggplot(s0, aes(t, SC, colour=Gruppo, fill=Gruppo)) +
  geom_vline(xintercept=10, colour="#c3c2b7", linewidth=0.4) +
  geom_ribbon(aes(ymin=SC.LCL, ymax=SC.UCL), alpha=0.12, colour=NA) +
  geom_line(data=s0[s0$t<=10,], linewidth=0.9) +
  geom_line(data=s0[s0$t>=10,], linewidth=0.9, linetype="22") +
  geom_point(data=oss, size=3, shape=21, stroke=0.8, colour="#fcfcfb") +
  scale_colour_manual(values=col) + scale_fill_manual(values=col) +
  scale_x_continuous(breaks=seq(0,20,2), expand=c(0.01,0)) + scale_y_continuous(labels=function(x) format(x, decimal.mark=","), expand=c(0,0)) +
  coord_cartesian(ylim=c(0.5,1.005)) +
  labs(tag="B", x="Numero di piante", y="Copertura del campione", colour=NULL, fill=NULL) + tema
p <- (pa | pb) + plot_layout(guides="collect") & theme(legend.position="top", legend.justification="left")
ggsave("out/fig_rarefazione_estrapolazione.png", p, width=10, height=4.8, dpi=300)
ggsave("out/fig_rarefazione_estrapolazione.pdf", p, width=10, height=4.8)
