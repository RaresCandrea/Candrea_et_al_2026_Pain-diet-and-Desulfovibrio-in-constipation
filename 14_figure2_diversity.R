# =============================================================================
# 14 - Figure 2: alpha and beta diversity, IBS-C vs IBS-D
# A  Shannon index (genus level, targeted panel)   B  Genera detected (richness)
# C  PCoA on Aitchison distance (genus panel) with PERMANOVA
# Statistics are read from the analysis outputs (11_main_analyses.R), not recomputed here.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })
Sys.setlocale("LC_ALL", "C.UTF-8")
res <- readRDS("output/res_11.rds"); d <- res$d
A1 <- read.csv("output/A1_alpha_diversity.csv"); A2 <- read.csv("output/A2_PERMANOVA.csv")
dir.create("figures", showWarnings = FALSE)

FONT <- getOption("fig.font", "Liberation Sans")                            # metric-compatible with Arial
col  <- c("IBS-C" = "#2a78d6", "IBS-D" = "#eb6834")   # colour-blind-safe pair (validated)
n    <- table(d$Subtype)
xlab <- c("IBS-C" = sprintf("IBS-C\n(n=%d)", n[["IBS-C"]]), "IBS-D" = sprintf("IBS-D\n(n=%d)", n[["IBS-D"]]))
fmtp <- function(p) if (p < 0.001) "italic(P)<0.001" else sprintf('italic(P)=="%.3f"', p)

th <- theme_classic(base_size = 8.5, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.35, colour = "grey20"), axis.ticks = element_line(linewidth = 0.35, colour = "grey20"),
        axis.text = element_text(colour = "grey15", size = 8), axis.title = element_text(size = 8.5),
        plot.title = element_text(size = 8.5, hjust = 0.5, face = "plain"), plot.tag = element_text(face = "bold", size = 11),
        legend.position = "none", plot.margin = margin(4, 6, 4, 4))

violin <- function(v, ylab, p, padj, ybreaks = waiver()) {
  dd <- data.frame(Subtype = d$Subtype, y = d[[v]])
  ggplot(dd, aes(Subtype, y)) +
    geom_violin(aes(colour = Subtype), fill = NA, trim = TRUE, linewidth = 0.45, scale = "width", width = 0.75) +
    geom_point(aes(colour = Subtype), position = position_jitter(width = 0.07, height = 0, seed = 1), size = 0.9, alpha = 0.55, stroke = 0) +
    stat_summary(fun = median, fun.min = function(z) quantile(z, .25), fun.max = function(z) quantile(z, .75),
                 geom = "pointrange", colour = "black", size = 0.35, linewidth = 0.6) +
    scale_colour_manual(values = col) + scale_x_discrete(labels = xlab) + scale_y_continuous(breaks = ybreaks) +
    labs(x = NULL, y = ylab, title = parse(text = paste0(fmtp(p), '*";"~', sub('italic\\(P\\)', 'italic(P)[adj]', fmtp(padj))))[[1]]) + th
}
pA <- violin("Shannon_genus", "Shannon index (genus level)", A1$p_MannWhitney[grepl("genus level", A1$Measure)],
             A1$p_adjusted[grepl("genus level", A1$Measure)])
pB <- violin("Richness_genus", "Genera detected (n)", A1$p_MannWhitney[grepl("richness", A1$Measure)],
             A1$p_adjusted[grepl("richness", A1$Measure)], ybreaks = seq(10, 24, 2))

pm <- A2 %>% filter(Analysis == "Aitchison, genus panel, adjusted (marginal terms)", Term == "Subtype")
o <- res$ord
pC <- ggplot(o, aes(PC1, PC2, colour = Subtype)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
  stat_ellipse(type = "norm", level = 0.95, linetype = "dashed", linewidth = 0.4) +
  geom_point(size = 1.3, alpha = 0.85, stroke = 0) +
  scale_colour_manual(values = col) +
  labs(x = sprintf("PCoA1 (%.1f%%)", res$ord_var[1]), y = sprintf("PCoA2 (%.1f%%)", res$ord_var[2]),
       title = parse(text = sprintf('"PERMANOVA,"~italic(R)^2=="%.3f"*","~%s', pm$R2, sub('italic\\(P\\)', 'italic(P)[adj]', fmtp(pm$p))))[[1]]) +
  th + theme(panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.35), axis.line = element_blank(),
             legend.position = "bottom", legend.title = element_blank(), legend.margin = margin(0, 0, 0, 0),
             legend.box.spacing = unit(1, "mm"), legend.text = element_text(size = 8), legend.key.size = unit(3, "mm")) +
  guides(colour = guide_legend(override.aes = list(linetype = 0, size = 1.8, alpha = 1)))

fig <- (pA | pB | pC) + plot_layout(widths = c(1, 1, 1.45)) + plot_annotation(tag_levels = "A")
W <- 174 / 25.4; H <- 66 / 25.4                       # double-column width, 174 mm
source("fig_save.R"); save_fig(fig, "figures/Figure2_diversity", H_mm = H * 25.4)
cat("Figure 2 written\n")
