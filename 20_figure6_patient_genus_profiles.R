# =============================================================================
# 20 - Figure 6: genus-level profiles of individual patients
# Genera with prevalence >= 10% (as in the differential abundance analysis); CLR-transformed relative
# abundances (zeros replaced by half the genus's minimum non-zero value), standardised per genus (z-score).
# Patients ordered within each subtype, and genera ordered, by hierarchical clustering
# (Ward's method; Euclidean distance on CLR values = Aitchison distance for patients).
# Annotation bars: Western food score and abdominal pain intensity (IBS-SSS item 1).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })
Sys.setlocale("LC_ALL", "C.UTF-8")
x <- readRDS("output/analysis_data.rds"); d <- as.data.frame(x$d); genus <- x$genus
d$Subtype <- factor(d$Subtype, levels = c("IBS-C", "IBS-D"))
clr_mat <- function(m) {                       # identical to 11_main_analyses.R
  m <- m / rowSums(m)
  for (j in seq_len(ncol(m))) { z <- m[, j] == 0; if (any(z)) m[z, j] <- min(m[!z, j]) / 2 }
  m <- m / rowSums(m); l <- log(m); l - rowMeans(l)
}
M <- as.matrix(d[, genus]); rownames(M) <- d$Subject
M <- M[, colMeans(M > 0) >= 0.10, drop = FALSE]
B1 <- read.csv("output/B1_DA_subtype.csv"); stopifnot(setequal(colnames(M), B1$Taxon[B1$Level == "Genus"]))   # same 22 genera as Figure 3B
C <- clr_mat(M)
Z <- scale(C)                                  # z-score per genus across all 64 patients

# ordering by hierarchical clustering
row_ord <- colnames(Z)[hclust(dist(t(Z)), method = "ward.D2")$order]
pat_ord <- unlist(lapply(levels(d$Subtype), function(s) {
  ids <- d$Subject[d$Subtype == s]; ids[hclust(dist(C[ids, ]), method = "ward.D2")$order] }))

long <- data.frame(Subject = rep(rownames(Z), ncol(Z)), Genus = rep(colnames(Z), each = nrow(Z)), z = as.vector(Z)) %>%
  left_join(d[, c("Subject", "Subtype")], by = "Subject") %>%
  mutate(Subject = factor(Subject, levels = pat_ord), Genus = factor(Genus, levels = rev(row_ord)))
nC <- sum(d$Subtype == "IBS-C"); nD <- sum(d$Subtype == "IBS-D")
strip_lab <- c("IBS-C" = sprintf("IBS-C (n=%d)", nC), "IBS-D" = sprintf("IBS-D (n=%d)", nD))
FONT <- getOption("fig.font", "Liberation Sans")
gen_lab <- setNames(lapply(levels(long$Genus), function(g) { z <- sub(" spp\\.$", "", g)
  if (z == "Desulfovibrio") bquote(bolditalic(.(z))) else bquote(italic(.(z))) }), levels(long$Genus))
base_th <- theme_minimal(base_size = 8, base_family = FONT) +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.spacing.x = unit(2.5, "mm"), legend.title = element_text(size = 7.5), legend.text = element_text(size = 7),
        legend.key.height = unit(3.2, "mm"), legend.key.width = unit(3, "mm"), plot.margin = margin(0, 2, 0, 2))

ann <- function(var, lab, cols, name, brks, lims, strip = FALSE) {
  a <- d[, c("Subject", "Subtype", var)]; names(a)[3] <- "v"; a$Subject <- factor(a$Subject, levels = pat_ord)
  ggplot(a, aes(Subject, 1, fill = v)) + geom_tile(colour = "white", linewidth = 0.25) +
    facet_grid(. ~ Subtype, scales = "free_x", space = "free_x", labeller = as_labeller(strip_lab)) +
    scale_fill_gradient(low = cols[1], high = cols[2], name = name, breaks = brks, limits = lims,
                        guide = guide_colourbar(direction = "horizontal", title.position = "top", order = 2)) +
    scale_y_continuous(breaks = 1, labels = lab, expand = c(0, 0)) + labs(x = NULL, y = NULL) + base_th +
    theme(axis.text.y = element_text(size = 7.5, colour = "grey15"),
          strip.text = if (strip) element_text(face = "bold", size = 8.5) else element_blank(),
          legend.key.width = unit(6, "mm"), legend.key.height = unit(2.4, "mm"))
}
aW <- ann("Western_score", "Western food score", c("#f2f0f7", "#54278f"), "Western food score", c(0, 5, 10, 15, 20), c(0, 21), strip = TRUE)
aP <- ann("IBSSSS_Q1_pain_intensity", "Pain intensity", c("#f5f5f5", "#3a3a3a"), "Pain intensity (0–100)", c(0, 25, 50, 75, 100), c(0, 100))

hm <- ggplot(long, aes(Subject, Genus, fill = pmax(pmin(z, 2.5), -2.5))) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(. ~ Subtype, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166ac", mid = "#f0efec", high = "#b2182b", midpoint = 0, limits = c(-2.5, 2.5),
                       breaks = c(-2, -1, 0, 1, 2), labels = c("−2", "−1", "0", "1", "2"),
                       name = "CLR abundance\n(z-score)", guide = guide_colourbar(order = 1)) +
  scale_y_discrete(labels = gen_lab, expand = c(0, 0)) + labs(x = "Patients", y = NULL) + base_th +
  theme(axis.text.y = element_text(size = 7.5, colour = "grey15"), strip.text = element_blank(),
        axis.title.x = element_text(size = 8, margin = margin(t = 3)))

# Legends are extracted and stacked in a separate column (patchwork guide collection fails for these scales)
leg <- function(p) { g <- ggplotGrob(p + theme(legend.position = "right", legend.justification = "left"))
  g$grobs[[which(g$layout$name == "guide-box-right")]] }
legs <- wrap_elements(full = grid::grobTree(
  gtable::gtable_col("legs", list(leg(hm), leg(aW), leg(aP)), width = grid::unit(1, "null"),
                     heights = grid::unit(c(30, 13, 13), "mm"))))
noleg <- theme(legend.position = "none")
fig <- ((aW + noleg) / (aP + noleg) / (hm + noleg)) + plot_layout(heights = c(0.045, 0.045, 1)) & theme(plot.margin = margin(1, 2, 1, 2))
fig <- wrap_elements(full = fig) + legs + plot_layout(widths = c(1, 0.23))
source("fig_save.R"); dir.create("figures", showWarnings = FALSE)
save_fig(fig, "figures/Figure6_patient_genus_profiles", H_mm = 120)
cat("Figure 6 written\n")
