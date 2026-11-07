# =============================================================================
# 15 - Supplementary tables: display values (JSON specs) built from the analysis outputs.
# Each spec is rendered to Word by 16_build_supp_docx.js in the style of Tables 1-2.
# Markup in strings: *text* = italic, _{x} = subscript, ^{x} = superscript.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(jsonlite) })
res <- readRDS("output/res_11.rds"); d <- res$d
dir.create("tables/supp", recursive = TRUE, showWarnings = FALSE)
fp <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
mq <- function(z, dec) { q <- quantile(z, c(.5, .25, .75), type = 7)
  sprintf(paste0("%.", dec, "f (%.", dec, "f, %.", dec, "f)"), q[1], q[2], q[3]) }
ci <- function(b, l, h, dec) gsub("-", "\u2212", sprintf(paste0("%.", dec, "f (%.", dec, "f to %.", dec, "f)"), b, l, h))   # typographic minus
nC <- sum(d$Subtype == "IBS-C"); nD <- sum(d$Subtype == "IBS-D")
row <- function(...) list(type = "data", cells = c(...))
sec <- function(label) list(type = "section", cells = label)

# ---- Table S1: alpha and beta diversity --------------------------------------------------
A1 <- read.csv("output/A1_alpha_diversity.csv"); A2 <- read.csv("output/A2_PERMANOVA.csv")
alpha <- list(c("Shannon_genus", "Shannon index, genus level", 2, "genus level"),
              c("Shannon_phylum", "Shannon index, phylum level", 2, "phylum level"),
              c("Richness_genus", "Genera detected, n", 1, "richness"))
rowsA <- lapply(alpha, function(a) { r <- A1[grepl(a[4], A1$Measure), ]; dec <- as.integer(a[3]); v <- a[1]
  row(a[2], mq(d[[v]][d$Subtype == "IBS-C"], dec), mq(d[[v]][d$Subtype == "IBS-D"], dec), fp(r$p_MannWhitney),
      ci(r$Adj_diff_D_vs_C, r$Adj_CI_low, r$Adj_CI_high, dec), fp(r$p_adjusted)) })
lab <- c("Aitchison, genus panel, unadjusted" = "Aitchison distance, genera, unadjusted",
         "Aitchison, genus panel, adjusted (marginal terms)" = "Aitchison distance, genera, adjusted",
         "Aitchison, species panel, adjusted (marginal terms)" = "Aitchison distance, species, adjusted",
         "Bray-Curtis, genus panel, adjusted (sensitivity)" = "Bray–Curtis dissimilarity, genera, adjusted")
rowsB <- lapply(names(lab), function(a) { r <- A2 %>% filter(Analysis == a, Term == "Subtype")
  row(lab[[a]], sprintf("%.3f", r$R2), sprintf("%.2f", r$F), fp(r$p), fp(r$Dispersion_p_subtype)) })
C1 <- read.csv("output/C1_functional_by_subtype.csv")
rowsC <- lapply(list(c("pH", "Faecal pH", 1), c("Firmicutes / bacteroidetes ratio", "Firmicutes/Bacteroidetes ratio", 1)), function(a) {
  r <- C1[C1$Indicator == a[1], ]; v <- d[[a[1]]]
  row(a[2], mq(v[d$Subtype == "IBS-C"], 1), mq(v[d$Subtype == "IBS-D"], 1), fp(r$p_MannWhitney), fp(r$p_adjusted)) })
S1 <- list(
  file = "Supplementary_Table_S1_diversity", label = "Supplementary Table S1.",
  title = "Alpha and beta diversity, faecal pH and Firmicutes/Bacteroidetes ratio by IBS subtype.",
  blocks = list(
    list(caption = "A. Alpha diversity",
         headers = c("Measure", sprintf("IBS-C (n=%d)", nC), sprintf("IBS-D (n=%d)", nD), "*P*",
                     "Adjusted difference (95% CI)", "*P*_{adj}"),
         widths = c(2900, 1600, 1600, 750, 2150, 750), align = c("left", rep("center", 5)), rows = rowsA),
    list(caption = "B. Beta diversity (PERMANOVA, IBS subtype term)",
         headers = c("Distance and model", "*R*^{2}", "Pseudo-*F*", "*P*", "*P*, dispersion"),
         widths = c(4350, 1200, 1200, 1200, 1800), align = c("left", rep("center", 4)), rows = rowsB),
    list(caption = "C. Faecal pH and Firmicutes/Bacteroidetes ratio",
         headers = c("Measure", sprintf("IBS-C (n=%d)", nC), sprintf("IBS-D (n=%d)", nD), "*P*", "*P*_{adj}"),
         widths = c(2900, 2100, 2100, 1325, 1325), align = c("left", rep("center", 4)), rows = rowsC)),
  footnote = paste0(
    "Diversity indices were computed from the relative abundances of the taxa in the laboratory 16S rRNA-based targeted panel. ",
    "(A) Values are median (first quartile, third quartile). *P*, Mann–Whitney U test; adjusted difference (IBS-D minus IBS-C), with 95% confidence interval (CI), ",
    "and *P*_{adj}, linear regression adjusted for age, sex, body mass index and proton pump inhibitor use, with heteroskedasticity-consistent (HC3) standard errors. ",
    "(B) PERMANOVA (9,999 permutations; marginal effects); adjusted models include age, sex, body mass index and proton pump inhibitor use. ",
    "Aitchison distance, Euclidean distance between centred log-ratio-transformed abundances of taxa with prevalence ≥ 10% ",
    "(zeros replaced by half of the smallest non-zero value of each taxon); Bray–Curtis dissimilarity computed on relative abundances. ",
    "*R*^{2}, proportion of variance explained by IBS subtype; *P*, dispersion, permutation test for homogeneity of multivariate dispersion between subtypes. ",
    "(C) Faecal pH and Firmicutes/Bacteroidetes ratio as reported by the laboratory; median (first quartile, third quartile); *P* and *P*_{adj} as in (A) (ratio modelled on the log(1 + x) scale). ",
    "IBS-C, irritable bowel syndrome with predominant constipation; IBS-D, irritable bowel syndrome with predominant diarrhea; ",
    "PERMANOVA, permutational multivariate analysis of variance."))

write_json(S1, "tables/supp/S1.json", auto_unbox = TRUE, pretty = TRUE)
cat("S1 spec written\n")

# ---- Table S2: differential abundance (all taxa) + Desulfovibrio sensitivity analyses ----------
B1 <- read.csv("output/B1_DA_subtype.csv"); B2 <- read.csv("output/B2_Desulfovibrio_sensitivity.csv")
tx <- function(t, level) {                                # italics for genus and species names
  if (level == "Phylum") return(t)
  if (grepl(" spp\\.$", t)) return(paste0("*", sub(" spp\\.$", "", t), "* spp."))
  paste0("*", t, "*")
}
pct <- function(x) sprintf("%.0f", 100 * x)
f2 <- function(x) gsub("-", "−", sprintf("%.2f", x))
rowsDA <- list()
for (lv in c("Phylum", "Genus", "Species")) {
  rowsDA[[length(rowsDA) + 1]] <- sec(c(Phylum = "Phyla", Genus = "Genera", Species = "Species")[[lv]])
  b <- B1 %>% filter(Level == lv) %>% arrange(p_LinDA)
  for (i in seq_len(nrow(b))) { r <- b[i, ]
    rowsDA[[length(rowsDA) + 1]] <- row(tx(r$Taxon, lv), paste0(pct(r$Prev_IBS_C), " / ", pct(r$Prev_IBS_D)),
      sprintf("%.3f", r$Median_IBS_C), sprintf("%.3f", r$Median_IBS_D), fp(r$p_Wilcoxon), fp(r$q_Wilcoxon),
      sprintf("%s (%s to %s)", f2(r$log2FC), f2(r$log2FC - qnorm(.975) * r$SE), f2(r$log2FC + qnorm(.975) * r$SE)),
      fp(r$p_LinDA), fp(r$q_LinDA), fp(r$p_CLR), fp(r$q_CLR)) }
}
lin <- B2 %>% filter(grepl("^LinDA", Analysis))
sens_lab <- c("LinDA, adjusted for age, sex, BMI, PPI (main analysis)" = "LinDA, adjusted for age, sex, BMI and PPI use (main analysis)",
              "LinDA, unadjusted" = "LinDA, unadjusted",
              "LinDA, adjusted for age, sex, BMI (PPI not in model)" = "LinDA, adjusted for age, sex and BMI (PPI use not in the model)",
              "LinDA, additionally adjusted for western food score (rank)" = "LinDA, additionally adjusted for the western food score")
rowsS <- lapply(seq_len(nrow(lin)), function(i) { r <- lin[i, ]
  row(if (grepl("PPI users excluded", r$Analysis)) sub("n = ", "n=", paste0(r$Analysis, ", adjusted for age, sex and BMI")) else sens_lab[[r$Analysis]], sprintf("log_{2}FC %s (%s to %s)", f2(r$log2FC), f2(r$log2FC - qnorm(.975) * r$SE), f2(r$log2FC + qnorm(.975) * r$SE)), fp(r$p)) })
lm_r <- B2 %>% filter(grepl("^Linear", Analysis)); fi_r <- B2 %>% filter(grepl("^Detected", Analysis))
dif <- regmatches(lm_r$Note, regexpr("-?[0-9.]+ \\(-?[0-9.]+ to -?[0-9.]+\\)", lm_r$Note))
rowsS <- c(rowsS, list(
  row("Linear regression on log_{10} relative abundance, adjusted for age, sex, BMI and PPI use", paste("Difference", gsub("-", "−", dif)), fp(lm_r$p)),
  row("Detection (relative abundance > 0), Fisher's exact test",
      sprintf("%d/%d (%.1f%%) vs. %d/%d (%.1f%%)", sum(d$`Desulfovibrio spp.`[d$Subtype == "IBS-C"] > 0), nC, 100 * mean(d$`Desulfovibrio spp.`[d$Subtype == "IBS-C"] > 0),
              sum(d$`Desulfovibrio spp.`[d$Subtype == "IBS-D"] > 0), nD, 100 * mean(d$`Desulfovibrio spp.`[d$Subtype == "IBS-D"] > 0)), fp(fi_r$p))))
S2 <- list(
  file = "Supplementary_Table_S2_differential_abundance", label = "Supplementary Table S2.", landscape = TRUE,
  title = "Differential abundance of gut microbial taxa between IBS-D and IBS-C, and sensitivity analyses for *Desulfovibrio*.",
  blocks = list(
    list(caption = "A. All taxa with prevalence ≥ 10%",
         headers = c("Taxon", "Prevalence, % (IBS-C / IBS-D)", "Median, IBS-C (%)", "Median, IBS-D (%)", "*P*, Wilcoxon", "*q*, Wilcoxon",
                     "LinDA log_{2}FC (95% CI)", "*P*, LinDA", "*q*, LinDA", "*P*, CLR", "*q*, CLR"),
         widths = c(2700, 1500, 1150, 1150, 1050, 1050, 2100, 1000, 1000, 1000, 1000), align = c("left", rep("center", 10)), rows = rowsDA),
    list(caption = "B. *Desulfovibrio* spp.: sensitivity analyses",
         headers = c("Analysis", "Estimate (95% CI)", "*P*"),
         widths = c(8200, 4700, 1800), align = c("left", "center", "center"), rows = rowsS)),
  footnote = paste0(
    sprintf("IBS-C, n=%d; IBS-D, n=%d. Relative abundances (%%) from the laboratory 16S rRNA-based targeted panel; taxa detected in at least 10%% of patients were analysed. ", nC, nD),
    "*P*, Wilcoxon, unadjusted Mann–Whitney U test. LinDA, linear model for differential abundance analysis of compositional data, adjusted for age, sex, body mass index (BMI) and proton pump inhibitor (PPI) use; ",
    "log_{2}FC, log_{2} fold change of IBS-D relative to IBS-C (negative values = lower in IBS-D). CLR, centred log-ratio transformation followed by linear regression with the same covariates and heteroskedasticity-consistent (HC3) standard errors (second method). ",
    "*q*, Benjamini–Hochberg false discovery rate within each taxonomic level. Taxa are ordered by the LinDA *P* value within each level; taxon names as reported by the laboratory."))
write_json(S2, "tables/supp/S2.json", auto_unbox = TRUE, pretty = TRUE)
cat("S2 spec written\n")

# ---- Table S3: key taxa x symptoms (partial Spearman) ------------------------------------------
C3 <- read.csv("output/C3_taxa_x_symptoms.csv")
r2 <- function(x) gsub("-", "−", sub("^-0\\.00$", "0.00", sprintf("%.2f", x)))
cell3 <- function(z) sprintf("%s (%s to %s)", r2(z$partial_rho), r2(z$CI_low), r2(z$CI_high))
rowsT <- list()
for (t in unique(C3$Taxon)) {
  rowsT[[length(rowsT) + 1]] <- sec(tx(t, "Genus"))
  for (s in unique(C3$Symptom)) {
    g <- function(G) C3[C3$Taxon == t & C3$Symptom == s & C3$Group == G, ]
    a <- g("IBS-C"); b <- g("IBS-D"); o <- g("All")
    rowsT[[length(rowsT) + 1]] <- row(s, cell3(a), fp(a$p), fp(a$q_within_group), cell3(b), fp(b$p), fp(b$q_within_group),
                                      cell3(o), fp(o$p), fp(o$q_within_group)) }
}
S3 <- list(
  file = "Supplementary_Table_S3_taxa_symptoms", label = "Supplementary Table S3.", landscape = TRUE,
  title = "Adjusted Spearman correlations between key gut microbial taxa and IBS symptoms, by IBS subtype.",
  blocks = list(list(caption = NULL,
    headers = c("Taxon and symptom", sprintf("IBS-C (n=%d), ρ (95%% CI)", nC), "*P*", "*q*",
                sprintf("IBS-D (n=%d), ρ (95%% CI)", nD), "*P*", "*q*", sprintf("All (n=%d), ρ (95%% CI)", nC + nD), "*P*", "*q*"),
    widths = c(2900, 2050, 800, 800, 2050, 800, 800, 2050, 800, 800), align = c("left", rep("center", 9)), rows = rowsT)),
  footnote = paste0(
    "ρ, Spearman correlation coefficient adjusted for age, sex, body mass index and proton pump inhibitor use (and IBS subtype in all patients), ",
    "computed as the correlation of rank residuals, with Fisher *z* 95% confidence interval (CI). *q*, Benjamini–Hochberg false discovery rate across the 56 taxon–symptom pairs within each group. ",
    "Symptoms: items and total score of the IBS Severity Scoring System (IBS-SSS) and weekly bowel movement frequency. ",
    "Relative abundances from the laboratory 16S rRNA-based targeted panel; taxon names as reported by the laboratory. ",
    "IBS-C, irritable bowel syndrome with predominant constipation; IBS-D, irritable bowel syndrome with predominant diarrhea."))
write_json(S3, "tables/supp/S3.json", auto_unbox = TRUE, pretty = TRUE)
cat("S3 spec written\n")

# ---- Table S4: diet - Desulfovibrio - symptoms (adjusted Spearman), by subtype ----------------
ch4 <- bind_rows(read.csv("output/E1_path_A_diet_to_Desulfovibrio.csv"), read.csv("output/E2_path_B_Desulfovibrio_to_symptoms.csv"),
                 read.csv("output/E3_path_C_diet_to_symptoms.csv"))
ix4 <- read.csv("output/E4_interaction_by_subtype.csv")
cell4 <- function(z) sprintf("%s (%s to %s)", r2(z$partial_rho), r2(z$CI_low), r2(z$CI_high))
ex4 <- c("Western food score", "Processed meat", "Red meat", "Red + processed meat")
oc4 <- c("Abdominal pain intensity", "Days with pain", "Abdominal distension", "IBS-SSS total")
row4 <- function(lbl, pr, oc) {
  g <- function(G) ch4[ch4$Predictor == pr & ch4$Outcome == oc & ch4$Group == G, ]
  a <- g("IBS-C"); b <- g("IBS-D"); o <- g("All"); pi <- ix4$p_interaction_ranks[ix4$Predictor == pr & ix4$Outcome == oc]
  stopifnot(nrow(a) == 1, nrow(b) == 1, nrow(o) == 1, length(pi) == 1)
  row(lbl, cell4(a), fp(a$p), fp(a$q_chain_within_group), cell4(b), fp(b$p), fp(b$q_chain_within_group),
      cell4(o), fp(o$p), fp(pi)) }
rows4 <- c(list(sec("Dietary exposure and *Desulfovibrio* spp.")), lapply(ex4, function(x) row4(x, x, "Desulfovibrio")),
           list(sec("*Desulfovibrio* spp. and symptoms")), lapply(oc4, function(y) row4(sub("Abdominal p", "P", y), "Desulfovibrio", y)))
sec4 <- c("abdominal pain intensity", "days with pain", "abdominal distension", "IBS-SSS total score")
for (k in seq_along(oc4)) { y <- oc4[k]; rows4 <- c(rows4, list(sec(paste("Dietary exposure and", sec4[k]))),
                          lapply(ex4, function(x) row4(x, x, y))) }
S4 <- list(
  file = "Supplementary_Table_S4_diet_Desulfovibrio_symptoms", label = "Supplementary Table S4.", landscape = TRUE,
  title = "Associations between dietary exposures, *Desulfovibrio* spp. and IBS symptoms, by IBS subtype.",
  blocks = list(list(caption = NULL,
    headers = c("Association", sprintf("IBS-C (n=%d), ρ (95%% CI)", nC), "*P*", "*q*",
                sprintf("IBS-D (n=%d), ρ (95%% CI)", nD), "*P*", "*q*", sprintf("All (n=%d), ρ (95%% CI)", nC + nD), "*P*",
                "*P* for interaction"),
    widths = c(2700, 2050, 750, 750, 2050, 750, 750, 2050, 750, 1250), align = c("left", rep("center", 9)), rows = rows4)),
  footnote = paste0(
    "ρ, Spearman correlation coefficient adjusted for age, sex, body mass index and proton pump inhibitor use (and IBS subtype in all patients), ",
    "computed as the correlation of rank residuals, with Fisher *z* 95% confidence interval (CI). ",
    "*q*, Benjamini–Hochberg false discovery rate across the 24 associations shown, within each subtype. ",
    "*P* for interaction, exposure × subtype term in a linear model on ranks with the same covariates and heteroscedasticity-consistent (HC3) standard errors. ",
    "Western food score, sum of the red meat, processed meat and processed sweets scores; red + processed meat, sum of the two meat scores (food frequency questionnaire). ",
    "*Desulfovibrio* spp. relative abundance from the laboratory 16S rRNA-based targeted panel, log_{10}-transformed after adding half of the smallest non-zero value. ",
    "Symptoms: items of the IBS Severity Scoring System (IBS-SSS) and its total score. ",
    "IBS-C, irritable bowel syndrome with predominant constipation; IBS-D, irritable bowel syndrome with predominant diarrhea."))
write_json(S4, "tables/supp/S4.json", auto_unbox = TRUE, pretty = TRUE)
cat("S4 spec written\n")

# ---- Table S5: sensitivity analyses for the diet - Desulfovibrio - symptom associations ----------
lin <- read.csv("output/E1s_linear_sensitivity.csv")
mr <- read.csv("output/E7_mediation.csv"); ml <- read.csv("output/E7s_mediation_linear.csv")
evs <- read.csv("output/E8_E_values.csv")
d2 <- function(x, dec = 2) gsub("-", "−", sprintf(paste0("%.", dec, "f"), x))
cil <- function(b, l, h, dec = 2) sprintf("%s (%s to %s)", d2(b, dec), d2(l, dec), d2(h, dec))
# A. linear models (HC3) on the original scales
rowl <- function(lbl, pr, oc) {
  g <- function(G) lin[lin$Predictor == pr & lin$Outcome == oc & lin$Group == G, ]
  a <- g("IBS-C"); b <- g("IBS-D"); stopifnot(nrow(a) == 1, nrow(b) == 1)
  row(lbl, cil(a$Beta, a$CI_low, a$CI_high), fp(a$p), cil(b$Beta, b$CI_low, b$CI_high), fp(b$p)) }
rowsL <- c(list(sec("Dietary exposure and *Desulfovibrio* spp. (log_{10} %)")), lapply(ex4, function(x) rowl(x, x, "Desulfovibrio")),
           list(sec("*Desulfovibrio* spp. (per 1 log_{10} unit) and symptoms")),
           lapply(oc4, function(y) rowl(sub("Abdominal p", "P", y), "Desulfovibrio", y)))
for (k in seq_along(oc4)) rowsL <- c(rowsL, list(sec(paste("Dietary exposure and", sec4[k]))),
                                     lapply(ex4, function(x) rowl(x, x, oc4[k])))
lcs <- function(o) c("Abdominal pain intensity" = "pain intensity", "Days with pain" = "days with pain", "IBS-SSS total" = "IBS-SSS total",
                      "Abdominal distension" = "abdominal distension", "Desulfovibrio" = "Desulfovibrio")[[o]]
# B. mediation (IBS-C and all patients; not run in IBS-D, where no exposure-outcome association was observed)
pct <- function(x) sub("^\u22120\\.0$", "0.0", d2(100 * x, 1))
rowm <- function(G, ex, oc) {
  a <- mr[mr$Group == G & mr$Exposure == ex & mr$Outcome == oc, ]; l <- ml[ml$Group == G & ml$Exposure == ex & ml$Outcome == oc, ]
  stopifnot(nrow(a) == 1, nrow(l) == 1)
  row(sprintf("%s → %s", ex, lcs(oc)),
      sprintf("%s (%s to %s)", pct(a$ACME), pct(a$ACME_CI_low), pct(a$ACME_CI_high)), fp(a$ACME_p),
      sprintf("%s (%s to %s)", pct(a$ADE), pct(a$ADE_CI_low), pct(a$ADE_CI_high)), fp(a$ADE_p),
      sprintf("%s (%s to %s)", pct(a$Total), pct(a$Total_CI_low), pct(a$Total_CI_high)), fp(a$Total_p),
      if (a$Total_p < 0.05) sprintf("%.0f%%", 100 * a$Prop_mediated) else "–", fp(l$ACME_p)) }
exm <- c("Western food score", "Processed meat", "Red + processed meat"); ocm <- c("Abdominal pain intensity", "Days with pain", "IBS-SSS total")
rowsM <- list()
for (G in c("IBS-C", "All")) {
  rowsM <- c(rowsM, list(sec(if (G == "IBS-C") sprintf("IBS-C (n=%d)", nC) else sprintf("All patients (n=%d)", nC + nD))))
  for (ex in exm) for (oc in ocm) rowsM <- c(rowsM, list(rowm(G, ex, oc))) }
# C. E-values for associations with P < 0.05
rowsE <- list()
for (G in c("IBS-C", "All")) {
  e <- evs[evs$Group == G, ]
  rowsE <- c(rowsE, list(sec(if (G == "IBS-C") sprintf("IBS-C (n=%d)", nC) else sprintf("All patients (n=%d)", nC + nD))))
  for (i in seq_len(nrow(e))) rowsE <- c(rowsE, list(row(
    sprintf("%s → %s", sub("^Desulfovibrio$", "*Desulfovibrio* spp.", e$Predictor[i]),
            sub("^Desulfovibrio$", "*Desulfovibrio* spp.", lcs(e$Outcome[i]))),
    r2(e$partial_rho[i]), fp(e$p[i]), sprintf("%.1f", e$E_value_point[i]), sprintf("%.1f", e$E_value_CI[i])))) }
S5 <- list(
  file = "Supplementary_Table_S5_sensitivity_mediation", label = "Supplementary Table S5.", landscape = TRUE,
  title = "Sensitivity analyses for the associations between dietary exposures, *Desulfovibrio* spp. and IBS symptoms: linear models, mediation analysis and E-values.",
  blocks = list(
    list(caption = "A. Linear regression on the original scales",
         headers = c("Association", sprintf("IBS-C (n=%d), β (95%% CI)", nC), "*P*", sprintf("IBS-D (n=%d), β (95%% CI)", nD), "*P*"),
         widths = c(5200, 3300, 1200, 3300, 1200), align = c("left", rep("center", 4)), rows = rowsL),
    list(caption = "B. Mediation by *Desulfovibrio* spp. (rank-based; effects in percentile points of the outcome)",
         headers = c("Exposure → outcome", "Indirect effect (95% CI)", "*P*", "Direct effect (95% CI)", "*P*",
                     "Total effect (95% CI)", "*P*", "Proportion mediated", "Linear model, *P* for indirect effect"),
         widths = c(3300, 1850, 700, 1850, 700, 1850, 700, 1300, 1700), align = c("left", rep("center", 8)), rows = rowsM),
    list(caption = "C. E-values for associations with *P* < 0.05",
         headers = c("Association", "Adjusted ρ", "*P*", "E-value, point estimate", "E-value, confidence limit"),
         widths = c(5200, 1800, 1400, 2900, 2900), align = c("left", rep("center", 4)), rows = rowsE)),
  footnote = paste0(
    "All models adjusted for age, sex, body mass index and proton pump inhibitor use (and IBS subtype in all patients). ",
    "(A) β, regression coefficient per 1-point increase in the dietary score or per 1 log_{10} unit of *Desulfovibrio* spp. relative abundance, ",
    "with heteroscedasticity-consistent (HC3) 95% confidence interval (CI); outcomes on their original scales (*Desulfovibrio* spp., log_{10} %; symptoms, IBS Severity Scoring System [IBS-SSS] points). ",
    "(B) Causal mediation analysis (Imai et al.) with *Desulfovibrio* spp. as mediator; exposure, mediator and outcome expressed as percentile ranks within the group; ",
    "effects compare the 75th with the 25th percentile of the exposure and are expressed in percentile points of the outcome; nonparametric bootstrap 95% CI (5000 resamples). ",
    "Proportion mediated is reported only when the total effect was significant (*P* < 0.05). The last column gives *P* for the indirect effect in the same model on the original scales. ",
    "Mediation was not assessed in IBS-D, in which no association between dietary exposures and symptoms was observed. ",
    "(C) E-value, minimum strength of association (risk ratio scale) that an unmeasured confounder would need to have with both the exposure and the outcome to fully explain the observed association ",
    "(point estimate) or to shift its confidence interval to include the null (confidence limit), computed from the adjusted Spearman ρ (VanderWeele and Ding). ",
    "IBS-C, irritable bowel syndrome with predominant constipation; IBS-D, irritable bowel syndrome with predominant diarrhea."))
write_json(S5, "tables/supp/S5.json", auto_unbox = TRUE, pretty = TRUE)
cat("S5 spec written\n")

# ---- Table S6: clinical, lifestyle and dietary factors and gut microbiota composition (PERMANOVA) ---------
D1 <- read.csv("output/D1_lifestyle_PERMANOVA.csv")
cats6 <- list(Clinical = c("Sex", "Age", "BMI", "PPI use", "Antibiotics (12 months)", "Anxiety disorder"),
              Lifestyle = c("Current smoking", "Physical activity >=150 min/week", "PSQI global score", "Alcohol", "Coffee",
                            "Water intake", "Healthy lifestyle score"),
              Diet = c("Fruits", "Vegetables", "Legumes", "Whole grains", "Red meat", "Processed meat", "White meat", "Fish",
                       "Dairy products", "Eggs", "Processed sweets", "Plant food score", "Western food score"))
stopifnot(setequal(unlist(cats6), unique(D1$Factor)))
lab6 <- c("Physical activity >=150 min/week" = "Physical activity ≥150 min/week", "BMI" = "Body mass index",
          "PSQI global score" = "Sleep quality (PSQI global score)", "PPI use" = "PPI use^{a}",
          "Antibiotics (12 months)" = "Antibiotics (12 months)")
c6 <- function(f, G) { z <- D1[D1$Factor == f & D1$Group == G, ]
  if (nrow(z) == 0) return(c("–", "–", "–")); c(sprintf("%.3f", z$R2), fp(z$p), fp(z$q_within_group)) }
rows6 <- list()
for (k in names(cats6)) { rows6 <- c(rows6, list(sec(k)))
  for (f in cats6[[k]]) rows6 <- c(rows6, list(row(if (f %in% names(lab6)) lab6[[f]] else f, c6(f, "IBS-C"), c6(f, "IBS-D"), c6(f, "All")))) }
S6 <- list(
  file = "Supplementary_Table_S6_factors_composition", label = "Supplementary Table S6.", landscape = FALSE,
  title = "Clinical, lifestyle and dietary factors and gut microbiota composition, by IBS subtype (PERMANOVA).",
  blocks = list(list(caption = NULL,
    headers = c("Factor", sprintf("IBS-C (n=%d), *R*^{2}", nC), "*P*", "*q*", sprintf("IBS-D (n=%d), *R*^{2}", nD), "*P*", "*q*",
                sprintf("All (n=%d), *R*^{2}", nC + nD), "*P*", "*q*"),
    widths = c(2750, 1100, 700, 700, 1100, 700, 700, 1100, 700, 700), align = c("left", rep("center", 9)), rows = rows6)),
  footnote = paste0(
    "*R*^{2}, proportion of variance in community composition explained by each factor, from permutational multivariate analysis of variance (PERMANOVA; ",
    "Aitchison distance on centred log-ratio-transformed relative abundances of the genera detected in at least 10% of samples; marginal terms; 4999 permutations), ",
    "with each factor tested separately and adjusted for age and sex (and IBS subtype in all patients). ",
    "Dietary frequency scores, water intake and composite scores were entered as ranks within each group. ",
    "*q*, Benjamini–Hochberg false discovery rate across the factors tested within each group. ",
    "^{a}Not tested in IBS-D (2 users). Antibiotics, use in the 12 months preceding the clinical assessment. ",
    "PPI, proton pump inhibitor; PSQI, Pittsburgh Sleep Quality Index. ",
    "IBS-C, irritable bowel syndrome with predominant constipation; IBS-D, irritable bowel syndrome with predominant diarrhea."))
write_json(S6, "tables/supp/S6.json", auto_unbox = TRUE, pretty = TRUE)
cat("S6 spec written\n")

# ---- Table S7: exploratory analyses - Methanobrevibacter and anxiety disorder -------------------------
F1 <- read.csv("output/F1_methanogens.csv"); F2 <- read.csv("output/F2_anxiety_DA.csv"); F2b <- read.csv("output/F2b_anxiety_community.csv")
g7 <- function(pat) { z <- F1[grepl(pat, F1$Analysis, fixed = TRUE), ]; stopifnot(nrow(z) == 1); z }
mbC <- sum(d$Methano_pos[d$Subtype == "IBS-C"]); mbD <- sum(d$Methano_pos[d$Subtype == "IBS-D"])
anx_n <- sum(d$Anxiety_dx == 1)
e1 <- g7("detected: IBS-D vs IBS-C"); e2 <- g7("Desulfovibrio (log10) by"); e3 <- g7("Bowel movements per week by")
s_all <- g7("Spearman (All)"); s_c <- g7("Spearman (IBS-C)"); s_d <- g7("Spearman (IBS-D)")
ci7 <- function(z) sprintf("%s to %s", d2(z$CI_low), d2(z$CI_high))
rowsA7 <- list(
  row(sprintf("Detection, IBS-D vs IBS-C (%d/%d, %.1f%% vs %d/%d, %.1f%%), OR^{a}", mbD, nD, 100 * mbD / nD, mbC, nC, 100 * mbC / nC),
      d2(e1$Estimate), ci7(e1), fp(e1$p)),
  row("*Desulfovibrio* spp. (log_{10} %), detected vs not detected, β^{b}", d2(e2$Estimate), ci7(e2), fp(e2$p)),
  row("Bowel movements per week, detected vs not detected, β^{b}", d2(e3$Estimate), ci7(e3), fp(e3$p)),
  row(sprintf("Correlation with *Desulfovibrio* spp., all patients (n=%d), ρ^{c}", nC + nD), r2(s_all$Estimate), "–", fp(s_all$p)),
  row(sprintf("Correlation with *Desulfovibrio* spp., IBS-C (n=%d), ρ^{c}", nC), r2(s_c$Estimate), "–", fp(s_c$p)),
  row(sprintf("Correlation with *Desulfovibrio* spp., IBS-D (n=%d), ρ^{c}", nD), r2(s_d$Estimate), "–", fp(s_d$p)))
sh <- F2b[grepl("^Shannon", F2b$Analysis), ]; pm <- F2b[grepl("PERMANOVA", F2b$Analysis), ]
d1a <- read.csv("output/D1_lifestyle_PERMANOVA.csv") %>% filter(Group == "All", Factor == "Anxiety disorder")
stopifnot(nrow(sh) == 1, nrow(pm) == 1, nrow(d1a) == 1)
rowsB7 <- list(
  row("Shannon index (genus level), β^{b}", d2(sh$Beta), sprintf("%s to %s", d2(sh$CI_low), d2(sh$CI_high)), fp(sh$p)),
  row("Community composition, PERMANOVA *R*^{2}, adjusted for age, sex, BMI, PPI use and IBS subtype", sprintf("%.3f", pm$Beta), "–", fp(pm$p)),
  row("Community composition, PERMANOVA *R*^{2}, adjusted for age, sex and IBS subtype (Figure 7)", sprintf("%.3f", d1a$R2), "–", fp(d1a$p)))
rowsC7 <- list()
for (lv in c("Genus", "Species")) {
  rowsC7 <- c(rowsC7, list(sec(c(Genus = "Genera", Species = "Species")[[lv]])))
  b <- F2 %>% filter(Level == lv) %>% arrange(p)
  for (i in seq_len(nrow(b))) { r <- b[i, ]
    rowsC7 <- c(rowsC7, list(row(tx(r$Taxon, lv), sprintf("%s (%s to %s)", d2(r$log2FC), d2(r$log2FC - qnorm(.975) * r$SE), d2(r$log2FC + qnorm(.975) * r$SE)),
                                 fp(r$p), fp(r$q)))) } }
S7 <- list(
  file = "Supplementary_Table_S7_exploratory", label = "Supplementary Table S7.", landscape = FALSE,
  title = "Exploratory analyses: *Methanobrevibacter* spp. and anxiety disorder.",
  blocks = list(
    list(caption = sprintf("A. *Methanobrevibacter* spp. (detected in %d of %d patients)", mbC + mbD, nC + nD),
         headers = c("Analysis", "Estimate", "95% CI", "*P*"), widths = c(6000, 1300, 1900, 1000), align = c("left", rep("center", 3)), rows = rowsA7),
    list(caption = sprintf("B. Anxiety disorder (n=%d) and gut microbiota diversity and composition", anx_n),
         headers = c("Analysis", "Estimate", "95% CI", "*P*"), widths = c(6000, 1300, 1900, 1000), align = c("left", rep("center", 3)), rows = rowsB7),
    list(caption = "C. Anxiety disorder and differential abundance (LinDA)",
         headers = c("Taxon", "log_{2} fold change, anxiety vs no anxiety (95% CI)", "*P*", "*q*"), widths = c(4200, 3800, 1100, 1100),
         align = c("left", rep("center", 3)), rows = rowsC7)),
  footnote = paste0(
    "Unless stated otherwise, models were adjusted for age, sex, body mass index (BMI), proton pump inhibitor (PPI) use and IBS subtype. ",
    "^{a}Odds ratio from Firth's penalised logistic regression. ^{b}Linear regression with heteroscedasticity-consistent (HC3) standard errors. ",
    "^{c}Unadjusted Spearman correlation of relative abundances. ",
    "(B) PERMANOVA, Aitchison distance on genera detected in at least 10% of samples; marginal terms. ",
    "(C) LinDA on taxa detected in at least 10% of samples; *q*, Benjamini–Hochberg false discovery rate within each taxonomic level. ",
    "Taxon names as reported by the laboratory. CI, confidence interval."))
write_json(S7, "tables/supp/S7.json", auto_unbox = TRUE, pretty = TRUE)
cat("S7 spec written\n")
