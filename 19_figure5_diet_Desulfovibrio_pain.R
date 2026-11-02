# =============================================================================
# 19 - Figure 5: Diet, Desulfovibrio and abdominal pain by IBS subtype
# A-C  Raw data for the three key pairs (IBS-C and IBS-D overlaid); annotations give the
#      adjusted Spearman rho and P within each subtype (age, sex, BMI, PPI) and P for interaction.
#      Lines are least-squares fits shown for visual guidance only.
# D    Adjusted Spearman rho (95% CI) for every diet / Desulfovibrio / symptom pair, by subtype;
#      filled = P < 0.05; dagger = P for interaction with subtype < 0.05 (rank-based).
# All statistics are read from the analysis outputs (11_main_analyses.R).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })
Sys.setlocale("LC_ALL", "C.UTF-8")
res <- readRDS("output/res_11.rds"); d <- res$d
d$Subtype <- factor(d$Subtype, levels = c("IBS-C", "IBS-D"))
chain <- bind_rows(read.csv("output/E1_path_A_diet_to_Desulfovibrio.csv"),
                   read.csv("output/E2_path_B_Desulfovibrio_to_symptoms.csv"),
                   read.csv("output/E3_path_C_diet_to_symptoms.csv")) %>% filter(Group != "All")
intx <- read.csv("output/E4_interaction_by_subtype.csv")
FONT <- getOption("fig.font", "Liberation Sans"); col <- c("IBS-C" = "#2a78d6", "IBS-D" = "#eb6834")
n <- table(d$Subtype)
fp <- function(p) if (p < 0.001) "italic(P)<0.001" else sprintf('italic(P)=="%.3f"', p)
fr <- function(r) sprintf('"%s"', gsub("-", "−", sprintf("%.2f", r)))
th <- theme_classic(base_size = 8.5, base_family = FONT) +
  theme(axis.line = element_line(linewidth = 0.35, colour = "grey20"), axis.ticks = element_line(linewidth = 0.35, colour = "grey20"),
        axis.text = element_text(colour = "grey15", size = 7.8), axis.title = element_text(size = 8.2),
        plot.tag = element_text(face = "bold", size = 11), legend.position = "none", plot.margin = margin(4, 6, 2, 4))

# ---- A-C: scatter plots ---------------------------------------------------------------------
v <- d$`Desulfovibrio spp.`; floor_v <- min(v[v > 0]) / 2
d$Desulfo_plot <- log10(v + floor_v)               # identical to the analysis variable Desulfo_log
ybrk <- c(log10(floor_v), log10(c(0.01, 0.1, 1))); ylab_d <- c("n.d.", "0.01", "0.1", "1")
lab_tax <- expression(italic("Desulfovibrio")~"(%, log"[10]*")")
lab_pain <- "Abdominal pain intensity (0–100)"

sc <- function(x, y, pred, outc, xlab, ylab, ylim, jit_x = 0, ybreaks = waiver(), ylabels = waiver()) {
  s <- chain %>% filter(Predictor == pred, Outcome == outc)
  pint <- intx %>% filter(Predictor == pred, Outcome == outc) %>% pull(p_interaction_ranks)
  rc <- s[s$Group == "IBS-C", ]; rd <- s[s$Group == "IBS-D", ]
  txt <- c(paste0('"IBS-C: "*rho==', fr(rc$partial_rho), '*","~', fp(rc$p)),
           paste0('"IBS-D: "*rho==', fr(rd$partial_rho), '*","~', fp(rd$p)),
           paste0('"Interaction: "*', fp(pint)))
  ytxt <- ylim[2] + diff(ylim) * c(0.31, 0.20, 0.09)     # above the panel (clip off)
  ggplot(d, aes(.data[[x]], .data[[y]], colour = Subtype)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
    geom_point(position = position_jitter(width = jit_x, height = 0, seed = 3), size = 1.2, alpha = 0.75, stroke = 0) +
    annotate("text", x = -Inf, y = ytxt, label = txt, parse = TRUE, hjust = 0, size = 2.4, family = FONT,
             colour = c(col, "grey20")) +
    scale_colour_manual(values = col) +
    scale_y_continuous(breaks = ybreaks, labels = ylabels) +
    coord_cartesian(ylim = ylim, clip = "off") +
    labs(x = xlab, y = ylab) + th + theme(plot.margin = margin(26, 6, 2, 4))
}
yl_d <- range(d$Desulfo_plot) + c(-0.1, 0.1)
pA <- sc("Western_score", "IBSSSS_Q1_pain_intensity", "Western food score", "Abdominal pain intensity",
         "Western food score", lab_pain, ylim = c(0, 100), jit_x = 0.15, ybreaks = seq(0, 100, 25)) + labs(tag = "A")
pB <- sc("Western_score", "Desulfo_plot", "Western food score", "Desulfovibrio",
         "Western food score", lab_tax, ylim = yl_d, jit_x = 0.15, ybreaks = ybrk, ylabels = ylab_d) + labs(tag = "B")
pC <- sc("Desulfo_plot", "IBSSSS_Q1_pain_intensity", "Desulfovibrio", "Abdominal pain intensity",
         lab_tax, lab_pain, ylim = c(0, 100), ybreaks = seq(0, 100, 25)) +
  scale_x_continuous(breaks = ybrk, labels = ylab_d) + labs(tag = "C")

# ---- D: forest matrix -----------------------------------------------------------------------
preds <- c("Western food score", "Processed meat", "Red meat", "Red + processed meat", "Desulfovibrio")
outs <- c("Desulfovibrio", "Abdominal pain intensity", "Days with pain", "Abdominal distension", "IBS-SSS total")
out_lab <- c("Desulfovibrio" = 'italic("Desulfovibrio")', "Abdominal pain intensity" = '"Pain intensity"',
             "Days with pain" = '"Days with pain"', "Abdominal distension" = '"Distension"', "IBS-SSS total" = '"IBS-SSS total"')
fd <- chain %>% mutate(Predictor = factor(Predictor, levels = rev(preds)),
                       Outcome = factor(Outcome, levels = outs, labels = out_lab[outs]),
                       Group = factor(Group, levels = c("IBS-C", "IBS-D")), sig = p < 0.05)
dag <- intx %>% filter(p_interaction_ranks < 0.05) %>%
  mutate(Predictor = factor(Predictor, levels = rev(preds)), Outcome = factor(Outcome, levels = outs, labels = out_lab[outs]))
pred_lab_e <- setNames(lapply(preds, function(z) if (z == "Desulfovibrio") quote(italic("Desulfovibrio")) else z), preds)
dodge <- position_dodge(width = 0.55, reverse = TRUE)   # IBS-C above IBS-D
pD <- ggplot(fd, aes(partial_rho, Predictor, colour = Group, group = Group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0, linewidth = 0.45, position = dodge, orientation = "y") +
  geom_point(aes(fill = interaction(Group, sig)), shape = 21, size = 1.7, stroke = 0.6, position = dodge) +
  geom_text(data = dag, aes(x = 0.97, y = Predictor, label = "†"), inherit.aes = FALSE, size = 3, family = FONT,
            colour = "grey15", hjust = 1, vjust = 0.35) +
  scale_colour_manual(values = col, name = NULL,
                      labels = c(sprintf("IBS-C (n=%d)", n[["IBS-C"]]), sprintf("IBS-D (n=%d)", n[["IBS-D"]]))) +
  scale_fill_manual(values = c(`IBS-C.TRUE` = col[["IBS-C"]], `IBS-D.TRUE` = col[["IBS-D"]],
                               `IBS-C.FALSE` = "white", `IBS-D.FALSE` = "white"), guide = "none") +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-0.5, 0, 0.5), labels = c("−0.5", "0", "0.5")) +
  scale_y_discrete(labels = pred_lab_e, drop = FALSE) +
  facet_grid(. ~ Outcome, labeller = label_parsed) +
  labs(x = expression("Adjusted Spearman "*rho*" (95% CI)"), y = NULL) +
  labs(tag = "D") + th + theme(legend.position = "bottom", legend.text = element_text(size = 8), legend.margin = margin(0, 0, 0, 0),
             strip.background = element_blank(), strip.text = element_text(size = 8.2, face = "bold"),
             panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
             panel.spacing.x = unit(3.5, "mm"), axis.text.y = element_text(size = 8))

fig <- wrap_elements(full = (pA | pB | pC)) / wrap_elements(full = pD) + plot_layout(heights = c(1, 1.1))
source("fig_save.R"); dir.create("figures", showWarnings = FALSE)
save_fig(fig, "figures/Figure5_diet_Desulfovibrio_pain", H_mm = 160)
cat("Figure 5 written\n")
