# =============================================================================
# 21 - Figure 7: clinical, lifestyle and dietary factors and gut microbiota composition
# PERMANOVA (Aitchison distance, genera with prevalence >= 10%, marginal terms, adjusted for age and sex;
# + subtype in all patients) from 11_main_analyses.R (output/D1_lifestyle_PERMANOVA.csv).
# Point size = R2; fill = P category; black ring = q < 0.10 (BH FDR across factors within each group).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
Sys.setlocale("LC_ALL", "C.UTF-8")
D1 <- read.csv("output/D1_lifestyle_PERMANOVA.csv")
d <- readRDS("output/res_11.rds")$d
nC <- sum(d$Subtype == "IBS-C"); nD <- sum(d$Subtype == "IBS-D")
cats <- list(Clinical = c("Sex", "Age", "BMI", "PPI use", "Antibiotics (12 months)", "Anxiety disorder"),
             Lifestyle = c("Current smoking", "Physical activity >=150 min/week", "PSQI global score", "Alcohol", "Coffee",
                           "Water intake", "Healthy lifestyle score"),
             Diet = c("Fruits", "Vegetables", "Legumes", "Whole grains", "Red meat", "Processed meat", "White meat", "Fish",
                      "Dairy products", "Eggs", "Processed sweets", "Plant food score", "Western food score"))
stopifnot(setequal(unlist(cats), unique(D1$Factor)))
lab <- c("Physical activity >=150 min/week" = "Physical activity ≥150 min/week", "BMI" = "Body mass index",
         "PSQI global score" = "Sleep quality (PSQI score)", "Antibiotics (12 months)" = "Antibiotics (12 months)")
grp_lab <- c("IBS-C" = sprintf("IBS-C\n(n=%d)", nC), "IBS-D" = sprintf("IBS-D\n(n=%d)", nD), "All" = sprintf("All\n(n=%d)", nC + nD))
pd <- D1 %>%
  mutate(Category = factor(sapply(Factor, function(f) names(cats)[sapply(cats, function(z) f %in% z)]), levels = names(cats)),
         Factor = factor(Factor, levels = rev(unlist(cats))),
         Group = factor(Group, levels = c("IBS-C", "IBS-D", "All")),
         Pcat = cut(p, c(-Inf, 0.01, 0.05, 0.10, Inf), labels = c("< 0.01", "0.01–0.049", "0.05–0.099", "≥ 0.10"), right = FALSE),
         fdr = q_within_group < 0.10)
FONT <- getOption("fig.font", "Liberation Sans")
p <- ggplot(pd, aes(Group, Factor)) +
  geom_point(aes(size = R2, fill = Pcat), shape = 21, colour = "grey55", stroke = 0.3) +
  geom_point(data = filter(pd, fdr), aes(size = R2), shape = 21, fill = NA, colour = "black", stroke = 0.9) +
  scale_fill_manual(values = c("< 0.01" = "#08306b", "0.01–0.049" = "#2171b5", "0.05–0.099" = "#9ecae1", "≥ 0.10" = "#f0f0f0"),
                    name = expression(italic(P)), drop = FALSE, guide = guide_legend(override.aes = list(size = 3.2), order = 1)) +
  scale_size_area(max_size = 5.5, limits = c(0, 0.09), breaks = c(0.02, 0.04, 0.06, 0.08), name = expression(italic(R)^2),
                  guide = guide_legend(order = 2, override.aes = list(fill = "#bdbdbd", colour = "grey55", stroke = 0.3))) +
  scale_x_discrete(labels = grp_lab, position = "top") +
  scale_y_discrete(labels = function(z) ifelse(z %in% names(lab), lab[z], z)) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8.5, base_family = FONT) +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3), axis.text = element_text(colour = "grey15"),
        axis.text.x.top = element_text(size = 8, face = "bold"), axis.text.y = element_text(size = 8),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 90, face = "bold", size = 8.5),
        panel.spacing.y = unit(2.5, "mm"), legend.title = element_text(size = 8), legend.text = element_text(size = 7.5),
        plot.margin = margin(4, 4, 4, 4))
source("fig_save.R"); dir.create("figures", showWarnings = FALSE)
save_fig(p, "figures/Figure7_lifestyle_composition", W_mm = 130, H_mm = 150)
cat("Figure 7 written\n")
