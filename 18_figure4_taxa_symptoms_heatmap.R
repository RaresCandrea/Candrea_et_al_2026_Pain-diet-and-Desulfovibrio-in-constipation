# =============================================================================
# 18 - Figure 4: associations between key taxa and symptoms (adjusted Spearman), heatmap
# Partial Spearman correlations adjusted for age, sex, BMI and PPI use (+ subtype in 'All'),
# Fisher-z 95% CI; FDR (Benjamini-Hochberg) across the 8 taxa x 7 symptoms within each group.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr) })
Sys.setlocale("LC_ALL", "C.UTF-8")
res <- readRDS("output/res_11.rds"); d <- res$d
d$Subtype <- factor(d$Subtype, levels = c("IBS-C", "IBS-D"))

taxa <- c("Desulfovibrio spp.", "Bilophila wadsworthia", "Methanobrevibacter spp.", "Bifidobacterium spp.",
          "Faecalibacterium prausnitzii", "Akkermansia muciniphila", "Prevotella spp.", "Bacteroides spp.")
symp <- c(IBSSSS_Q1_pain_intensity = "Pain intensity", IBSSSS_Q2_pain_days = "Days with pain",
          IBSSSS_Q3_distension = "Distension", IBSSSS_Q4_bowel_dissatisfaction = "Bowel habit dissatisfaction",
          IBSSSS_Q5_interference = "Interference with life", IBSSSS_TOTAL = "IBS-SSS total", BM_week = "Bowel movements/week")
pspear <- function(a, b, dd, cv) {
  ra <- rank(dd[[a]]); rb <- rank(dd[[b]])
  Z <- model.matrix(as.formula(paste("~", paste(cv, collapse = "+"))), dd)
  r <- cor(lm.fit(Z, ra)$residuals, lm.fit(Z, rb)$residuals); k <- ncol(Z) - 1; n <- nrow(dd)
  df <- n - 2 - k; se <- 1 / sqrt(n - 3 - k)
  c(rho = r, lo = tanh(atanh(r) - qnorm(.975) * se), hi = tanh(atanh(r) + qnorm(.975) * se),
    p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df))
}
grp <- list(All = d, `IBS-C` = d[d$Subtype == "IBS-C", ], `IBS-D` = d[d$Subtype == "IBS-D", ])
C3 <- bind_rows(lapply(names(grp), function(g) {
  dd <- grp[[g]]; cv <- c("Age", "Sex", "BMI", "PPI", if (g == "All") "Subtype")
  bind_rows(lapply(taxa, function(t) bind_rows(lapply(names(symp), function(s) {
    r <- pspear(t, s, dd, cv)
    data.frame(Group = g, n = nrow(dd), Taxon = t, Symptom = symp[[s]], partial_rho = r[["rho"]],
               CI_low = r[["lo"]], CI_high = r[["hi"]], p = r[["p"]]) })))) %>%
    mutate(q_within_group = p.adjust(p, "BH")) }))
write.csv(C3, "output/C3_taxa_x_symptoms.csv", row.names = FALSE)

# ---- Figure 4 -------------------------------------------------------------------------------
FONT <- getOption("fig.font", "Liberation Sans")
hm <- C3 %>% filter(Group != "All") %>%
  mutate(Symptom = factor(Symptom, levels = symp),
         Taxon = factor(Taxon, levels = rev(taxa)),
         star = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""),
         lab = sprintf("%.2f%s", partial_rho, star), lab = sub("^-0\\.00", "0.00", lab), lab = gsub("-", "−", lab),
         fdr = q_within_group < 0.10,
         Group = factor(Group, levels = c("IBS-C", "IBS-D"),
                        labels = c(sprintf("IBS-C (n=%d)", sum(d$Subtype == "IBS-C")), sprintf("IBS-D (n=%d)", sum(d$Subtype == "IBS-D")))))
tax_lab <- setNames(lapply(taxa, function(z) if (grepl(" spp\\.$", z)) bquote(italic(.(sub(" spp\\.$", "", z)))~"spp.") else bquote(italic(.(z)))), taxa)
p <- ggplot(hm, aes(Symptom, Taxon, fill = partial_rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(data = filter(hm, fdr), fill = NA, colour = "black", linewidth = 0.7) +
  geom_text(aes(label = lab, colour = abs(partial_rho) > 0.42), size = 2.35, family = FONT) +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_fill_gradient2(low = "#2166ac", mid = "#f0efec", high = "#b2182b", midpoint = 0, limits = c(-0.6, 0.6),
                       oob = scales::squish, breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
                       labels = c("−0.6", "−0.3", "0", "0.3", "0.6"), name = expression(rho*" (adjusted)")) +
  scale_y_discrete(labels = tax_lab) + scale_x_discrete(position = "bottom") +
  facet_wrap(~ Group, nrow = 1) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8.5, base_family = FONT) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, colour = "grey15", size = 7.8),
        axis.text.y = element_text(colour = "grey15", size = 8), strip.text = element_text(face = "bold", size = 9),
        legend.position = "right", legend.key.height = unit(9, "mm"), legend.key.width = unit(2.6, "mm"),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7.5), plot.margin = margin(4, 4, 4, 4))
dir.create("figures", showWarnings = FALSE)
W <- 174 / 25.4; H <- 82 / 25.4
options(width = 200)
print(C3 %>% filter(p < 0.05) %>% mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)
source("fig_save.R"); save_fig(p, "figures/Figure4_taxa_symptoms_heatmap", H_mm = H * 25.4)
cat("Figure 4 written\n")
