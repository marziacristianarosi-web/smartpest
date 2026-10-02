# Figure per i test delle sei ipotesi (Corythucha arcuata - comunità fungina).
# Uso: Rscript fig_test_ipotesi.R <file_excel> <cartella_output>
# Fig. 1 ricchezza per pianta; Fig. 2 ordinamento PCoA (Jaccard); Fig. 3 incidenza dei taxa.
# I p-value riportati nelle figure provengono da test_ipotesi.R e vanno aggiornati se i dati cambiano.
suppressMessages({library(readxl); library(vegan); library(ggplot2); library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
f_dati <- if (length(args) >= 1) args[1] else "dati.xlsx"
out    <- if (length(args) >= 2) args[2] else "figure_test"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# Stessa palette delle curve di accumulo (verificata per daltonismo); la forma del punto è una codifica ridondante
col <- c(estate_farnia = "#2a78d6", autunno_farnia = "#eb6834", estate_cerro = "#1baf7a")
shp <- c(estate_farnia = 16, autunno_farnia = 17, estate_cerro = 15)
lab <- c(estate_farnia = expression(italic("Q. robur")*", estate"), autunno_farnia = expression(italic("Q. robur")*", autunno"),
         estate_cerro = expression(italic("Q. cerris")*", estate"))
tema <- theme_classic(base_size = 9, base_family = "Liberation Sans") +   # font metricamente equivalente ad Arial, con corsivo
  theme(axis.text = element_text(colour = "grey20"), axis.line = element_line(colour = "grey40", linewidth = 0.3),
        axis.ticks = element_line(colour = "grey40", linewidth = 0.3), legend.position = "none",
        plot.tag = element_text(face = "bold", size = 10), plot.title = element_text(size = 9, face = "plain"))
salva <- function(p, nome, w, h) {
  ggsave(file.path(out, paste0(nome, ".pdf")), p, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(out, paste0(nome, ".png")), p, width = w, height = h, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
}

# ---- Dati a livello di pianta ----
d  <- as.data.frame(read_excel(f_dati))
tx <- setdiff(names(d), c("Season", "Host_species", "N_pianta", "Replicate", "N_insetti"))
key <- paste(d$Season, d$Host_species, d$N_pianta, sep = "|")
pa  <- (rowsum(as.matrix(d[tx]), key, reorder = FALSE) > 0) * 1
m   <- d[!duplicated(key), c("Season", "Host_species", "N_pianta")]; rownames(m) <- unique(key); m <- m[rownames(pa), ]
m$R <- rowSums(pa)
m$gr <- factor(ifelse(m$Host_species == "Quercus_cerris", "estate_cerro",
                      ifelse(m$Season == "Summer", "estate_farnia", "autunno_farnia")),
               levels = c("estate_farnia", "autunno_farnia", "estate_cerro"))

# ================= Fig. 1: ricchezza per pianta =================
media_ic <- function(x) { t <- t.test(x); data.frame(y = mean(x), ymin = t$conf.int[1], ymax = t$conf.int[2]) }
S <- droplevels(m[m$Host_species == "Quercus_robur", ])
S$x <- as.numeric(factor(S$gr, levels = c("estate_farnia", "autunno_farnia")))
set.seed(1); sp <- setNames(runif(nlevels(factor(S$N_pianta)), -0.07, 0.07), levels(factor(S$N_pianta)))
S$xj <- S$x + sp[S$N_pianta]
yl <- c(0, max(m$R) + 1.5)
p1a <- ggplot(S, aes(x, R)) +
  geom_line(aes(xj, group = N_pianta), colour = "grey75", linewidth = 0.35) +
  geom_point(aes(xj, colour = gr, shape = gr), size = 1.8) +
  stat_summary(aes(group = gr), fun.data = media_ic, geom = "errorbar", width = 0, linewidth = 0.6,
               position = position_nudge(x = 0.25)) +
  stat_summary(aes(group = gr), fun = mean, geom = "point", shape = 23, size = 2.6, fill = "white",
               colour = "grey15", position = position_nudge(x = 0.25)) +
  annotate("text", x = 1.5, y = yl[2] - 0.3, label = "p = 0,002", size = 2.9, family = "Liberation Sans") +
  scale_x_continuous(breaks = 1:2, labels = c("Estate", "Autunno"), limits = c(0.6, 2.4)) +
  scale_colour_manual(values = col) + scale_shape_manual(values = shp) +
  scale_y_continuous(limits = yl, breaks = seq(0, 12, 2), expand = c(0, 0)) +
  labs(x = NULL, y = "Ricchezza in taxa per pianta", title = expression(italic("Q. robur")*" \u2013 stagione")) + tema
H <- droplevels(m[m$Season == "Summer", ])
p1b <- ggplot(H, aes(gr, R)) +
  geom_point(aes(colour = gr, shape = gr), size = 1.8, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  stat_summary(fun.data = media_ic, geom = "errorbar", width = 0, linewidth = 0.6, position = position_nudge(x = 0.25)) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.6, fill = "white", colour = "grey15",
               position = position_nudge(x = 0.25)) +
  annotate("text", x = 1.5, y = yl[2] - 0.3, label = "p = 0,72", size = 2.9, family = "Liberation Sans") +
  scale_x_discrete(labels = c(estate_farnia = expression(italic("Q. robur")), estate_cerro = expression(italic("Q. cerris")))) +
  scale_colour_manual(values = col) + scale_shape_manual(values = shp) +
  scale_y_continuous(limits = yl, breaks = seq(0, 12, 2), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Estate \u2013 specie ospite") + tema
salva(p1a + p1b + plot_annotation(tag_levels = "A"), "Fig1_ricchezza", 140, 75)

# ================= Fig. 2: PCoA su distanze di Jaccard =================
pcoa <- function(dd) {
  o <- wcmdscale(vegdist(pa[rownames(dd), ], "jaccard", binary = TRUE), eig = TRUE)
  pv <- 100 * o$eig[1:2] / sum(o$eig[o$eig > 0])
  list(sc = cbind(dd, A1 = o$points[, 1], A2 = o$points[, 2]), pv = pv)
}
grafico_pcoa <- function(o, titolo, testo, coppie = FALSE) {
  sc <- o$sc
  cen <- aggregate(cbind(A1, A2) ~ gr, sc, mean)
  p <- ggplot(sc, aes(A1, A2))
  if (coppie) {
    w <- reshape(sc[, c("N_pianta", "gr", "A1", "A2")], idvar = "N_pianta", timevar = "gr", direction = "wide")
    p <- p + geom_segment(data = w, aes(x = A1.estate_farnia, y = A2.estate_farnia,
                                        xend = A1.autunno_farnia, yend = A2.autunno_farnia),
                          colour = "grey80", linewidth = 0.3)
  }
  p + geom_hline(yintercept = 0, colour = "grey90", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey90", linewidth = 0.3) +
    geom_point(aes(colour = gr, shape = gr), size = 2, alpha = 0.9,
               position = position_jitter(width = 0.008, height = 0.008, seed = 2)) +
    geom_point(data = cen, aes(fill = gr), shape = 21, size = 4, colour = "white", stroke = 0.8) +
    scale_colour_manual(values = col, labels = lab, name = NULL, limits = names(col), drop = FALSE) +
    scale_shape_manual(values = shp, labels = lab, name = NULL, limits = names(col), drop = FALSE) +
    scale_fill_manual(values = col, guide = "none", limits = names(col)) +
    coord_equal() +
    labs(x = sprintf("PCoA 1 (%.0f%%)", o$pv[1]), y = sprintf("PCoA 2 (%.0f%%)", o$pv[2]), title = titolo, subtitle = testo) +
    tema + theme(legend.position = "bottom", legend.text = element_text(size = 8), plot.subtitle = element_text(size = 7.5, colour = "grey25"))
}
oS <- pcoa(S); oH <- pcoa(H)
p2a <- grafico_pcoa(oS, expression(italic("Q. robur")*" \u2013 stagione"),
                    "PERMANOVA: R² = 0,50; p = 0,002\nPERMDISP: p = 0,001", coppie = TRUE)
p2b <- grafico_pcoa(oH, "Estate \u2013 specie ospite", "PERMANOVA: R² = 0,02; p = 0,78\nPERMDISP: p = 0,62")
salva((p2a + p2b + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + plot_annotation(tag_levels = "A"),
      "Fig2_PCoA", 174, 100)

# ================= Fig. 3: incidenza dei taxa =================
inc <- aggregate(pa, list(gr = m$gr), sum)
lg <- reshape(inc, direction = "long", varying = tx, v.names = "n", timevar = "taxon", times = tx, idvar = "gr")
nm <- function(x) { x <- sub("Thrichoderma", "Trichoderma", gsub("_", " ", x)); x <- sub(" [Ss]pp(\\d)$", " sp. \\1", x); sub(" [Ss]pp$", " sp.", x) }
dif <- unlist(inc[inc$gr == "estate_farnia", tx]) - unlist(inc[inc$gr == "autunno_farnia", tx])
ordine <- tx[order(dif, colSums(pa[, tx]))]
lg$taxon <- factor(nm(lg$taxon), levels = nm(ordine))
sig <- c("Nigrospora_Spp", "Fusarium_spp2", "Fusarium_spp3", "Botryosphaeria_dothidea")   # q < 0,05 (McNemar, domanda 3)
lg$etichetta <- ifelse(lg$n == 0, "", lg$n)
p3 <- ggplot(lg, aes(gr, taxon)) +
  geom_tile(aes(fill = n), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = etichetta, colour = n >= 6), size = 2.7, family = "Liberation Sans") +
  scale_fill_gradient(low = "#f1f4f9", high = "#1d4e8f", limits = c(0, 10), breaks = c(0, 5, 10),
                      name = "Piante con\nil taxon (su 10)") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_x_discrete(labels = c(estate_farnia = expression(atop(italic("Q. robur"), "estate")),
                              autunno_farnia = expression(atop(italic("Q. robur"), "autunno")),
                              estate_cerro = expression(atop(italic("Q. cerris"), "estate"))), position = "top") +
  scale_y_discrete(labels = function(x) parse(text = sapply(x, function(t) {
    e <- if (grepl(" sp\\.", t)) sprintf('italic("%s")*" %s"', sub(" sp\\..*", "", t), sub("^\\S+ ", "", t))
         else sprintf('italic("%s")', t)
    if (t %in% nm(sig)) paste0(e, '*" *"') else e }))) +
  labs(x = NULL, y = NULL) + tema +
  theme(
        axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right",
        legend.title = element_text(size = 8))
salva(p3, "Fig3_incidenza_taxa", 120, 110)
