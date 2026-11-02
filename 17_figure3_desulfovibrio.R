# =============================================================================
# 17 - Figure 3: Desulfovibrio and genus-level differential abundance, IBS-D vs IBS-C
# A  Desulfovibrio relative abundance by subtype (log10 scale)
# B  LinDA log2 fold changes (95% CI) for all genera with prevalence >= 10%
# Statistics are read from the analysis outputs (11_main_analyses.R).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })
Sys.setlocale("LC_ALL", "C.UTF-8")
res <- readRDS("output/res_11.rds"); d <- res$d
B1 <- read.csv("output/B1_DA_subtype.csv")
FONT <- getOption("fig.font", "Liberation Sans"); col <- c("IBS-C" = "#2a78d6", "IBS-D" = "#eb6834")
fmtp <- function(p) if (p < 0.001) "italic(P)<0.001" else sprintf('italic(P)=="%.3f"', p)
th <- theme_classic(base_size = 8.5, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.35, colour = "grey20"), axis.ticks = element_line(linewidth = 0.35, colour = "grey20"),
        axis.text = element_text(colour = "grey15", size = 8), axis.title = element_text(size = 8.5),
        plot.title = element_text(size = 8.5, hjust = 0.5), plot.tag = element_text(face = "bold", size = 11),
        legend.position = "none", plot.margin = margin(4, 6, 4, 4))

# ---- A: Desulfovibrio abundance -----------------------------------------------------------
ds <- B1 %>% filter(Taxon == "Desulfovibrio spp.")
v <- d$`Desulfovibrio spp.`; floor_v <- min(v[v > 0]) / 2          # same pseudo-count as the analysis variable
dd <- data.frame(Subtype = d$Subtype, y = log10(v + floor_v))
n <- table(d$Subtype)
brk <- c(log10(floor_v), log10(c(0.01, 0.1, 1))); lab <- c("n.d.", "0.01", "0.1", "1")   # n.d. = not detected (plotted at the pseudo-count)
pA <- ggplot(dd, aes(Subtype, y)) +
  geom_boxplot(aes(colour = Subtype), fill = NA, width = 0.5, outlier.shape = NA, linewidth = 0.45) +
  geom_point(aes(colour = Subtype), position = position_jitter(width = 0.12, height = 0, seed = 2), size = 1.1, alpha = 0.7, stroke = 0) +
  scale_colour_manual(values = col) +
  scale_x_discrete(labels = c("IBS-C" = sprintf("IBS-C\n(n=%d)", n[["IBS-C"]]), "IBS-D" = sprintf("IBS-D\n(n=%d)", n[["IBS-D"]]))) +
  scale_y_continuous(breaks = brk, labels = lab) +
  labs(x = NULL, y = expression(italic("Desulfovibrio")~"relative abundance (%, log"[10]*" scale)"),
       title = parse(text = paste0(fmtp(ds$p_Wilcoxon), '*";"~', sub("italic\\(P\\)", "italic(P)[adj]", fmtp(ds$p_LinDA))))[[1]]) + th

# ---- B: genus-level LinDA effects ------------------------------------------------------------
g <- B1 %>% filter(Level == "Genus") %>%
  mutate(lo = log2FC - qnorm(.975) * SE, hi = log2FC + qnorm(.975) * SE,
         genus = sub(" spp\\.$", "", Taxon), key = Taxon == "Desulfovibrio spp.", sig = p_LinDA < 0.05) %>%
  arrange(log2FC) %>% mutate(lab = factor(genus, levels = genus))
lab_expr <- setNames(lapply(g$genus, function(z) if (z == "Desulfovibrio") bquote(bolditalic(.(z))) else bquote(italic(.(z)))), g$genus)
pB <- ggplot(g, aes(log2FC, lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = key), width = 0, linewidth = 0.45, orientation = "y") +
  geom_point(aes(colour = key, fill = interaction(key, sig)), shape = 21, size = 1.8, stroke = 0.6) +
  scale_colour_manual(values = c(`FALSE` = "grey45", `TRUE` = "#c0392b")) +
  scale_fill_manual(values = c(`FALSE.FALSE` = "white", `FALSE.TRUE` = "grey45", `TRUE.TRUE` = "#c0392b", `TRUE.FALSE` = "white")) +
  scale_y_discrete(labels = lab_expr, expand = expansion(add = c(0.7, 1.6))) +
  annotate("text", x = -0.15, y = nrow(g) + 1, label = "\u2190 Lower in IBS-D", hjust = 1, size = 2.7, family = FONT, colour = "grey30") +
  annotate("text", x = 0.15, y = nrow(g) + 1, label = "Higher in IBS-D \u2192", hjust = 0, size = 2.7, family = FONT, colour = "grey30") +
  coord_cartesian(clip = "off") +
  labs(x = expression("log"[2]*" fold change, IBS-D vs. IBS-C (95% CI)"), y = NULL) +
  th + theme(axis.text.y = element_text(size = 7.5), panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3))

fig <- (pA | pB) + plot_layout(widths = c(0.75, 1.6)) + plot_annotation(tag_levels = "A")
dir.create("figures", showWarnings = FALSE)
W <- 174 / 25.4; H <- 100 / 25.4
source("fig_save.R"); save_fig(fig, "figures/Figure3_Desulfovibrio", H_mm = H * 25.4)
cat("Figure 3 written\n")
