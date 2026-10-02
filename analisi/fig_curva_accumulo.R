# Figura: curve di accumulo dei taxa per gruppo (unità = pianta), con IC 95%
# Eseguire con: Rscript fig_curva_accumulo.R (impostare il percorso del file Excel qui sotto)
file_dati <- "G:/Il mio Drive/UNIVERSITA'/LAVORO/DATI_GIULIA_TIZIANA/AnalisiDati_Corythucha.xlsx"
dir.create("out", showWarnings = FALSE)
suppressMessages({library(readxl); library(vegan); library(ggplot2)})
d <- as.data.frame(read_excel(file_dati)); tx <- names(d)[6:19]
key <- paste(d$Season, d$Host_species, d$N_pianta)
pa <- (rowsum(as.matrix(d[tx]), key, reorder=FALSE) > 0)*1
grp <- sub(" [^ ]+$","", rownames(pa))
lab <- c("Fall Quercus_robur"="Farnia - autunno", "Summer Quercus_robur"="Farnia - estate", "Summer Quercus_cerris"="Cerro - estate")
col <- c("Farnia - autunno"="#eb6834", "Farnia - estate"="#2a78d6", "Cerro - estate"="#1baf7a")
cur <- do.call(rbind, lapply(names(lab), function(g){
  sa <- specaccum(pa[grp==g,], method="exact")
  data.frame(Gruppo=lab[g], Piante=sa$sites, Taxa=sa$richness, sd=sa$sd)
}))
cur$Gruppo <- factor(cur$Gruppo, levels=names(col))
cur$lo <- pmax(0, cur$Taxa - 1.96*cur$sd); cur$hi <- cur$Taxa + 1.96*cur$sd
# riepilogo numerico
tab <- do.call(rbind, lapply(split(cur, cur$Gruppo), function(x) data.frame(
  Gruppo=x$Gruppo[1], Taxa_10=round(x$Taxa[10],1), Taxa_5=round(x$Taxa[5],1),
  Pct_a_5=round(100*x$Taxa[5]/x$Taxa[10]), Ultima_pianta=round(x$Taxa[10]-x$Taxa[9],2),
  Ultime_3=round(x$Taxa[10]-x$Taxa[7],2))))
print(tab, row.names=FALSE); write.csv2(tab, "out/tab_curva_accumulo.csv", row.names=FALSE)
fine <- cur[cur$Piante==10,]
p <- ggplot(cur, aes(Piante, Taxa, colour=Gruppo, fill=Gruppo)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.12, colour=NA) +
  geom_line(linewidth=0.9) +
  geom_point(size=2.4, stroke=0.6, shape=21, colour="#fcfcfb") +
  geom_text(data=fine, aes(label=paste0(Gruppo, ": ", round(Taxa), " taxa")), hjust=0, nudge_x=0.25,
            colour="#0b0b0b", size=3.6, show.legend=FALSE) +
  scale_colour_manual(values=col) + scale_fill_manual(values=col) +
  scale_x_continuous(breaks=1:10, limits=c(1, 14.2), expand=c(0.01,0)) +
  scale_y_continuous(breaks=seq(0,16,2), limits=c(0,16.5), expand=c(0,0)) +
  labs(x="Numero di piante campionate", y="Numero cumulato di taxa fungini",
       colour=NULL, fill=NULL) +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="#e6e5e0", linewidth=0.3),
        axis.line=element_line(colour="#8a8984", linewidth=0.3), axis.ticks=element_line(colour="#8a8984"),
        text=element_text(colour="#0b0b0b"), axis.text=element_text(colour="#52514e"),
        plot.subtitle=element_text(colour="#52514e", size=10), plot.title=element_text(size=13),
        legend.position="top", legend.justification="left",
        plot.background=element_rect(fill="#fcfcfb", colour=NA))
ggsave("out/fig_curve_accumulo.png", p, width=8, height=5.2, dpi=300)
ggsave("out/fig_curve_accumulo.pdf", p, width=8, height=5.2)
