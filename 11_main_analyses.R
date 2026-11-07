# =============================================================================
# 11 - Main analyses on the final database (all results -> one results workbook)
#
# Covariates in every adjusted model: age, sex, BMI, PPI use in the previous 12 months.
# Regression SEs: heteroskedasticity-consistent (HC3). FDR: Benjamini-Hochberg within
# each analysis family (stated in the q column name / README).
# Microbial data: relative abundances (%) from the laboratory 16S-based targeted panel.
#
# A  Microbial community (alpha diversity of the panel, Aitchison PERMANOVA, composition)
# B  Differential abundance IBS-D vs IBS-C (main result: Desulfovibrio) - LinDA + CLR + crude
# C  Functional indicators (by subtype; vs stool frequency); taxa/indicators x symptoms
# D  Lifestyle and diet -> community composition, within each subtype (PERMANOVA)
# E  Diet -> Desulfovibrio -> abdominal pain (paths A/B/C, LinDA, causal mediation, E-values)
# F  Exploratory: methanogens, anxiety, PPI context; sensitivity power
#
# Dietary frequency scores (FFQ items, plant/western scores, healthy lifestyle score) are ORDINAL:
# in every model they enter as ranks (rank-based / partial Spearman approach, ranks computed within
# the group analysed). Linear models on the original scale are reported as sensitivity analyses.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(vegan); library(MicrobiomeStat)
  library(sandwich); library(lmtest); library(brglm2); library(mediation)
})
set.seed(20261031)
x <- readRDS("output/analysis_data.rds"); d <- as.data.frame(x$d)
d$Subtype <- factor(d$Subtype, levels = c("IBS-C", "IBS-D"))
rownames(d) <- d$Subject
genus <- x$genus; species <- x$species; phyla <- x$phyla
COV <- "Age + Sex + BMI + PPI"
R <- list()                                   # every result table goes here

# ---- helpers ------------------------------------------------------------------------
hc <- function(m, term) {
  ct <- coeftest(m, vcov = vcovHC(m, "HC3")); ci <- coefci(m, vcov = vcovHC(m, "HC3"))
  data.frame(Beta = ct[term, 1], CI_low = ci[term, 1], CI_high = ci[term, 2], p = ct[term, 4])
}
q3 <- function(z) sprintf("%.3g (%.3g, %.3g)", median(z), quantile(z, .25), quantile(z, .75))
prev <- function(m) colMeans(m > 0)
clr_mat <- function(m) {                      # multiplicative replacement: half the feature's minimum non-zero
  m <- m / rowSums(m)
  for (j in seq_len(ncol(m))) { z <- m[, j] == 0; if (any(z)) m[z, j] <- min(m[!z, j]) / 2 }
  m <- m / rowSums(m); l <- log(m); l - rowMeans(l)
}
mat <- function(v, rows = seq_len(nrow(d))) { m <- as.matrix(d[rows, v, drop = FALSE]); rownames(m) <- d$Subject[rows]; m }
keep10 <- function(m) m[, prev(m) >= 0.10, drop = FALSE]
pspear <- function(a, b, data, cv = c("Age", "Sex", "BMI", "PPI")) {   # partial Spearman (rank residuals)
  ra <- rank(data[[a]]); rb <- rank(data[[b]])
  if (sd(ra) == 0 || sd(rb) == 0) return(c(rho = NA, p = NA))
  cv <- cv[sapply(cv, function(v) length(unique(data[[v]])) > 1)]
  Z <- model.matrix(as.formula(paste("~", paste(cv, collapse = "+"))), data)
  r <- cor(lm.fit(Z, ra)$residuals, lm.fit(Z, rb)$residuals); df <- nrow(data) - 2 - (ncol(Z) - 1)
  c(rho = r, p = 2 * pt(-abs(r * sqrt(df / (1 - r^2))), df))
}
d$Male <- as.integer(d$Sex == "M"); d$IBS_D <- as.integer(d$Subtype == "IBS-D")
ordinal <- c(grep("^FFQ_.*_score$", names(d), value = TRUE), "FFQ_water_L_day", "Plant_score", "Western_score", "Meat_RP", "HLS")
rank_ord <- function(dd, vars = ordinal) { for (v in intersect(vars, names(dd))) dd[[v]] <- rank(dd[[v]]); dd }
grp <- list(All = d, `IBS-C` = d[d$Subtype == "IBS-C", ], `IBS-D` = d[d$Subtype == "IBS-D", ])
fcov <- function(g) if (g == "All") paste(COV, "+ Subtype") else COV

# =========================================================================================
# A. MICROBIAL COMMUNITY
# =========================================================================================
alpha_v <- c(Shannon_genus = "Shannon index, genus level (panel)", Shannon_phylum = "Shannon index, phylum level (panel)",
             Richness_genus = "Genera detected (panel richness)")
R$A1_alpha_diversity <- bind_rows(lapply(names(alpha_v), function(v) {
  y <- d[[v]]
  adj <- hc(lm(as.formula(paste(v, "~ Subtype +", COV)), d), "SubtypeIBS-D")
  noppi <- hc(lm(as.formula(paste(v, "~ Subtype + Age + Sex + BMI")), d), "SubtypeIBS-D")
  data.frame(Measure = alpha_v[[v]], IBS_C_median_IQR = q3(y[d$Subtype == "IBS-C"]), IBS_D_median_IQR = q3(y[d$Subtype == "IBS-D"]),
             p_MannWhitney = wilcox.test(y ~ d$Subtype, exact = FALSE)$p.value,
             Adj_diff_D_vs_C = adj$Beta, Adj_CI_low = adj$CI_low, Adj_CI_high = adj$CI_high, p_adjusted = adj$p,
             p_adjusted_without_PPI = noppi$p) }))

Gk <- keep10(mat(genus)); Sk <- keep10(mat(species))
Dg <- dist(clr_mat(Gk)); Ds <- dist(clr_mat(Sk)); Dbc <- vegdist(Gk, "bray")
perm <- function(D, f, label, data = d, nperm = 9999) {
  a <- adonis2(as.formula(paste("D ~", f)), data = data, by = "margin", permutations = nperm)
  data.frame(Analysis = label, Term = rownames(a)[seq_len(nrow(a) - 2)], Df = a$Df[seq_len(nrow(a) - 2)],
             R2 = a$R2[seq_len(nrow(a) - 2)], F = a$F[seq_len(nrow(a) - 2)], p = a$`Pr(>F)`[seq_len(nrow(a) - 2)])
}
R$A2_PERMANOVA <- bind_rows(
  perm(Dg, "Subtype", "Aitchison, genus panel, unadjusted"),
  perm(Dg, paste("Subtype +", COV), "Aitchison, genus panel, adjusted (marginal terms)"),
  perm(Ds, paste("Subtype +", COV), "Aitchison, species panel, adjusted (marginal terms)"),
  perm(Dbc, paste("Subtype +", COV), "Bray-Curtis, genus panel, adjusted (sensitivity)"))
R$A2_PERMANOVA$Dispersion_p_subtype <- NA
R$A2_PERMANOVA$Dispersion_p_subtype[R$A2_PERMANOVA$Term == "Subtype"] <- c(
  permutest(betadisper(Dg, d$Subtype), permutations = 9999)$tab$`Pr(>F)`[1],
  permutest(betadisper(Dg, d$Subtype), permutations = 9999)$tab$`Pr(>F)`[1],
  permutest(betadisper(Ds, d$Subtype), permutations = 9999)$tab$`Pr(>F)`[1],
  permutest(betadisper(Dbc, d$Subtype), permutations = 9999)$tab$`Pr(>F)`[1])

comp <- function(v, level) bind_rows(lapply(v, function(t) data.frame(Level = level, Taxon = t,
  Mean_all = mean(d[[t]]), Mean_IBS_C = mean(d[[t]][d$Subtype == "IBS-C"]), Mean_IBS_D = mean(d[[t]][d$Subtype == "IBS-D"]),
  Prevalence_IBS_C = mean(d[[t]][d$Subtype == "IBS-C"] > 0), Prevalence_IBS_D = mean(d[[t]][d$Subtype == "IBS-D"] > 0))))
R$A3_composition <- bind_rows(comp(phyla, "Phylum"), comp(genus, "Genus"), comp(species, "Species")) %>%
  group_by(Level) %>% arrange(desc(Mean_all), .by_group = TRUE) %>% ungroup()
pc <- prcomp(clr_mat(Gk))
ord <- data.frame(Subject = d$Subject, Subtype = d$Subtype, PC1 = pc$x[, 1], PC2 = pc$x[, 2])
ord_var <- 100 * summary(pc)$importance[2, 1:2]

# =========================================================================================
# B. DIFFERENTIAL ABUNDANCE (main result: Desulfovibrio spp.)
# =========================================================================================
meta <- d[, c("Subtype", "Age", "Sex", "BMI", "PPI", "Western_score", "FFQ_processed_meat_score", "Meat_RP",
              "FFQ_red_meat_score", "Anxiety_dx")]
run_linda <- function(M, fml, term, rows = seq_len(nrow(d)), md = meta) {
  out <- linda(feature.dat = t(M[rows, , drop = FALSE] / 100), meta.dat = md[rows, , drop = FALSE], formula = fml,
               feature.dat.type = "proportion", prev.filter = 0.10, is.winsor = TRUE, outlier.pct = 0.03,
               adaptive = TRUE, p.adj.method = "BH", verbose = FALSE)$output[[term]]
  data.frame(Taxon = rownames(out), log2FC = out$log2FoldChange, SE = out$lfcSE, p = out$pvalue, q = out$padj)
}
run_clr <- function(M, rhs, term) {
  C <- clr_mat(keep10(M))
  bind_rows(lapply(colnames(C), function(t) { dd <- cbind(meta, y = C[, t])
    h <- hc(lm(as.formula(paste("y", rhs)), dd), term); data.frame(Taxon = t, CLR_beta = h$Beta, p_CLR = h$p) })) %>%
    mutate(q_CLR = p.adjust(p_CLR, "BH"))
}
crude <- function(v) bind_rows(lapply(v, function(t) { y <- d[[t]]
  data.frame(Taxon = t, Median_IBS_C = median(y[d$Subtype == "IBS-C"]), Median_IBS_D = median(y[d$Subtype == "IBS-D"]),
             Prev_IBS_C = mean(y[d$Subtype == "IBS-C"] > 0), Prev_IBS_D = mean(y[d$Subtype == "IBS-D"] > 0),
             p_Wilcoxon = wilcox.test(y ~ d$Subtype, exact = FALSE)$p.value) }))
da_level <- function(v, level) {
  M <- mat(v); keep <- colnames(M)[prev(M) >= 0.10]
  crude(keep) %>% mutate(q_Wilcoxon = p.adjust(p_Wilcoxon, "BH")) %>%
    left_join(run_linda(M, paste("~ Subtype +", COV), "SubtypeIBS-D") %>% rename(p_LinDA = p, q_LinDA = q), by = "Taxon") %>%
    left_join(run_clr(M, paste("~ Subtype +", COV), "SubtypeIBS-D"), by = "Taxon") %>%
    mutate(Level = level, Direction = ifelse(log2FC < 0, "lower in IBS-D", "higher in IBS-D"),
           Concordant_q10 = q_LinDA < 0.10 & q_CLR < 0.10 & sign(log2FC) == sign(CLR_beta)) %>%
    relocate(Level, Taxon) %>% arrange(p_LinDA)
}
R$B1_DA_subtype <- bind_rows(da_level(phyla, "Phylum"), da_level(genus, "Genus"), da_level(species, "Species"))

# Desulfovibrio: sensitivity analyses
DS <- "Desulfovibrio spp."; G <- mat(genus)
noppi <- which(d$PPI == 0)
sens <- list(
  c("LinDA, adjusted for age, sex, BMI, PPI (main analysis)", run_linda(G, paste("~ Subtype +", COV), "SubtypeIBS-D") %>% filter(Taxon == DS)),
  c("LinDA, unadjusted", run_linda(G, "~ Subtype", "SubtypeIBS-D") %>% filter(Taxon == DS)),
  c("LinDA, adjusted for age, sex, BMI (PPI not in model)", run_linda(G, "~ Subtype + Age + Sex + BMI", "SubtypeIBS-D") %>% filter(Taxon == DS)),
  c(sprintf("LinDA, PPI users excluded (n = %d)", length(noppi)), run_linda(G, "~ Subtype + Age + Sex + BMI", "SubtypeIBS-D", noppi) %>% filter(Taxon == DS)),
  c("LinDA, additionally adjusted for western food score (rank)", run_linda(G, paste("~ Subtype + Western_score +", COV), "SubtypeIBS-D", md = rank_ord(meta, "Western_score")) %>% filter(Taxon == DS)))
R$B2_Desulfovibrio_sensitivity <- bind_rows(lapply(sens, function(s) data.frame(Analysis = s[[1]], log2FC = s$log2FC, SE = s$SE, p = s$p)))
h <- hc(lm(as.formula(paste("Desulfo_log ~ Subtype +", COV)), d), "SubtypeIBS-D")
R$B2_Desulfovibrio_sensitivity <- bind_rows(R$B2_Desulfovibrio_sensitivity,
  data.frame(Analysis = "Linear model on log10 abundance (HC3), adjusted", log2FC = NA, SE = NA, p = h$p,
             Note = sprintf("difference in log10 abundance %.2f (%.2f to %.2f)", h$Beta, h$CI_low, h$CI_high)),
  data.frame(Analysis = "Detected (>0), Fisher's exact test", log2FC = NA, SE = NA,
             p = fisher.test(table(d$Subtype, d[[DS]] > 0))$p.value,
             Note = sprintf("detected in %d/%d IBS-C vs %d/%d IBS-D", sum(d[[DS]][d$Subtype == "IBS-C"] > 0), sum(d$Subtype == "IBS-C"),
                            sum(d[[DS]][d$Subtype == "IBS-D"] > 0), sum(d$Subtype == "IBS-D"))))

# =========================================================================================
# C. FUNCTIONAL INDICATORS AND SYMPTOM CORRELATIONS
# =========================================================================================
bio <- x$bioind; logv <- c("Firmicutes / bacteroidetes ratio", "Lactate production", "Mucin degradation (%)", "LPS-positive bacteria (%)")
R$C1_functional_by_subtype <- bind_rows(lapply(bio, function(v) {
  dd <- d; dd$y <- if (v %in% logv) log1p(d[[v]]) else d[[v]]
  a <- hc(lm(as.formula(paste("y ~ Subtype +", COV)), dd), "SubtypeIBS-D")
  data.frame(Indicator = v, IBS_C_median_IQR = q3(d[[v]][d$Subtype == "IBS-C"]), IBS_D_median_IQR = q3(d[[v]][d$Subtype == "IBS-D"]),
             p_MannWhitney = wilcox.test(d[[v]] ~ d$Subtype, exact = FALSE)$p.value, Model_scale = ifelse(v %in% logv, "log(1+x)", "original"),
             Adj_diff_D_vs_C = a$Beta, Adj_CI_low = a$CI_low, Adj_CI_high = a$CI_high, p_adjusted = a$p) })) %>%
  mutate(q_adjusted = p.adjust(p_adjusted, "BH"))

key_taxa <- c("Desulfovibrio spp.", "Bilophila wadsworthia", "Methanobrevibacter spp.", "Bifidobacterium spp.",
              "Faecalibacterium prausnitzii", "Akkermansia muciniphila", "Prevotella spp.", "Bacteroides spp.")
symptoms <- c(IBSSSS_Q1_pain_intensity = "Abdominal pain intensity", IBSSSS_Q2_pain_days = "Days with pain",
              IBSSSS_Q3_distension = "Abdominal distension", IBSSSS_Q4_bowel_dissatisfaction = "Dissatisfaction with bowel habit",
              IBSSSS_Q5_interference = "Interference with daily life", IBSSSS_TOTAL = "IBS-SSS total", BM_week = "Bowel movements per week")
R$C2_taxa_indicators_x_symptoms <- bind_rows(lapply(names(grp), function(g) {
  dd <- grp[[g]]; cv <- if (g == "All") c("Age", "Sex", "BMI", "PPI", "Subtype") else c("Age", "Sex", "BMI", "PPI")
  bind_rows(lapply(c(key_taxa, bio), function(a) bind_rows(lapply(names(symptoms), function(s) {
    r <- pspear(a, s, dd, cv); data.frame(Group = g, Variable = a, Symptom = symptoms[[s]], partial_rho = r[["rho"]], p = r[["p"]]) })))) %>%
    mutate(q_within_group = p.adjust(p, "BH")) }))

# =========================================================================================
# D. LIFESTYLE AND DIET -> COMMUNITY COMPOSITION (within each subtype; Su et al. 2023 style)
# =========================================================================================
factors <- c(Sex = "Sex", Age = "Age", BMI = "BMI", PPI = "PPI use", Smoking = "Current smoking", PA_active = "Physical activity >=150 min/week",
             PSQI_TOTAL = "PSQI global score", Anxiety_dx = "Anxiety disorder", Antibiotics_12mo = "Antibiotics (12 months)",
             FFQ_coffee_YES_NO = "Coffee", FFQ_alcohol_score = "Alcohol", FFQ_water_L_day = "Water intake",
             FFQ_fruits_score = "Fruits", FFQ_vegetables_score = "Vegetables", FFQ_legumes_score = "Legumes",
             FFQ_whole_grains_score = "Whole grains", FFQ_red_meat_score = "Red meat", FFQ_processed_meat_score = "Processed meat",
             FFQ_white_meat_score = "White meat", FFQ_fish_score = "Fish", FFQ_dairy_score = "Dairy products", FFQ_eggs_score = "Eggs",
             FFQ_processed_sweets_score = "Processed sweets", Plant_score = "Plant food score", Western_score = "Western food score",
             HLS = "Healthy lifestyle score")
R$D1_lifestyle_PERMANOVA <- bind_rows(lapply(names(grp), function(g) {
  dd <- rank_ord(grp[[g]]); Dx <- dist(clr_mat(keep10(mat(genus, which(d$Subject %in% dd$Subject)))))
  bind_rows(lapply(names(factors), function(f) {
    z <- dd[[f]]; if (length(unique(z)) < 2 || (length(unique(z)) == 2 && min(table(z)) < 3)) return(NULL)
    adj <- setdiff(c("Age", "Sex"), f); if (g == "All") adj <- c(adj, "Subtype")
    a <- adonis2(as.formula(paste("Dx ~", f, "+", paste(adj, collapse = " + "))), data = dd, by = "margin", permutations = 4999)
    data.frame(Group = g, Factor = factors[[f]], n = nrow(dd), R2 = a[f, "R2"], F = a[f, "F"], p = a[f, "Pr(>F)"]) })) %>%
    mutate(q_within_group = p.adjust(p, "BH")) }))

# =========================================================================================
# E. DIET -> DESULFOVIBRIO -> ABDOMINAL PAIN
#    Main analysis: partial Spearman correlations (ranks; adjusted for age, sex, BMI, PPI; + subtype in 'All').
#    FDR: BH across the whole chain (paths A + B + C) within each group.
#    Sensitivity: linear models (HC3) on the original scales.
# =========================================================================================
expo <- c(FFQ_processed_meat_score = "Processed meat", FFQ_red_meat_score = "Red meat",
          Meat_RP = "Red + processed meat", Western_score = "Western food score")
outc <- c(IBSSSS_Q1_pain_intensity = "Abdominal pain intensity", IBSSSS_Q2_pain_days = "Days with pain",
          IBSSSS_Q3_distension = "Abdominal distension", IBSSSS_TOTAL = "IBS-SSS total")
MED <- c(Desulfo_log = "Desulfovibrio")
cvs <- function(g) c("Age", "Sex", "BMI", "PPI", if (g == "All") "Subtype")
psp_ci <- function(a, b, dd, cv) {             # partial Spearman + Fisher-z 95% CI (df corrected for covariates)
  r <- pspear(a, b, dd, cv); k <- ncol(model.matrix(as.formula(paste("~", paste(cv, collapse = "+"))), dd)) - 1
  se <- 1 / sqrt(nrow(dd) - 3 - k); z <- atanh(r[["rho"]])
  data.frame(partial_rho = r[["rho"]], CI_low = tanh(z - qnorm(.975) * se), CI_high = tanh(z + qnorm(.975) * se), p = r[["p"]])
}
chain_rank <- bind_rows(lapply(names(grp), function(g) { dd <- grp[[g]]
  one <- function(pth, xs, ys) bind_rows(lapply(names(ys), function(y) bind_rows(lapply(names(xs), function(xv)
    cbind(Path = pth, Group = g, n = nrow(dd), Predictor = xs[[xv]], Outcome = ys[[y]], psp_ci(xv, y, dd, cvs(g)))))))
  bind_rows(one("A: diet -> Desulfovibrio", expo, MED), one("B: Desulfovibrio -> symptoms", MED, outc),
            one("C: diet -> symptoms", expo, outc)) %>% mutate(q_chain_within_group = p.adjust(p, "BH")) }))
R$E1_path_A_diet_to_Desulfovibrio <- chain_rank %>% filter(startsWith(Path, "A"))
R$E2_path_B_Desulfovibrio_to_symptoms <- chain_rank %>% filter(startsWith(Path, "B"))
R$E3_path_C_diet_to_symptoms <- chain_rank %>% filter(startsWith(Path, "C"))

# Sensitivity: linear models on the original scale (HC3)
chain_lin <- bind_rows(lapply(names(grp), function(g) { dd <- grp[[g]]
  one <- function(pth, xs, ys) bind_rows(lapply(names(ys), function(y) bind_rows(lapply(names(xs), function(xv)
    cbind(Path = pth, Group = g, n = nrow(dd), Predictor = xs[[xv]], Outcome = ys[[y]],
          hc(lm(as.formula(paste0("`", y, "` ~ ", xv, " + ", fcov(g))), dd), xv))))))
  bind_rows(one("A: diet -> Desulfovibrio", expo, MED), one("B: Desulfovibrio -> symptoms", MED, outc),
            one("C: diet -> symptoms", expo, outc)) %>% mutate(q_chain_within_group = p.adjust(p, "BH")) }))
R$E1s_linear_sensitivity <- chain_lin

# Does the association differ by subtype? Rank-based (main) and linear (sensitivity) interaction tests
intx <- function(y, xv) {
  dr <- d; dr[[y]] <- rank(d[[y]]); dr[[xv]] <- rank(d[[xv]])
  pr <- hc(lm(as.formula(paste0("`", y, "` ~ ", xv, " * Subtype + ", COV)), dr), paste0(xv, ":SubtypeIBS-D"))$p
  pl <- hc(lm(as.formula(paste0("`", y, "` ~ ", xv, " * Subtype + ", COV)), d), paste0(xv, ":SubtypeIBS-D"))$p
  data.frame(Outcome = c(MED, outc)[[y]], Predictor = c(MED, expo)[[xv]], p_interaction_ranks = pr, p_interaction_linear = pl)
}
R$E4_interaction_by_subtype <- bind_rows(
  bind_rows(lapply(names(expo), function(xv) intx("Desulfo_log", xv))),
  bind_rows(lapply(names(outc), function(y) intx(y, "Desulfo_log"))),
  bind_rows(lapply(names(outc), function(y) bind_rows(lapply(names(expo), function(xv) intx(y, xv))))))
R$E5_Spearman_unadjusted <- bind_rows(lapply(names(grp), function(g) { dd <- grp[[g]]
  pairs <- c(lapply(names(expo), function(xv) c(xv, "Desulfo_log")), lapply(names(outc), function(y) c("Desulfo_log", y)),
             unlist(lapply(names(expo), function(xv) lapply(names(outc), function(y) c(xv, y))), recursive = FALSE))
  bind_rows(lapply(pairs, function(pr) { ct <- suppressWarnings(cor.test(dd[[pr[1]]], dd[[pr[2]]], method = "spearman", exact = FALSE))
    data.frame(Group = g, n = nrow(dd), Variable_1 = c(expo, MED)[[pr[1]]], Variable_2 = c(outc, MED)[[pr[2]]],
               rho = unname(ct$estimate), p = ct$p.value) })) }))

# Compositional check of path A: LinDA, dietary exposure entered as rank (within group)
srb <- c("Desulfovibrio spp.", "Desulfobacter spp.", "Bilophila wadsworthia")
GS <- mat(c(genus, species))
rows_C <- which(d$Subtype == "IBS-C")
md_all <- rank_ord(meta); md_C <- meta; md_C[rows_C, ] <- rank_ord(meta[rows_C, ])
R$E6_LinDA_diet_to_sulfidogenic <- bind_rows(
  lapply(names(expo), function(xv) run_linda(GS, paste("~", xv, "+ Subtype +", COV), xv, md = md_all) %>%
    filter(Taxon %in% srb) %>% mutate(Exposure = expo[[xv]], Group = "All")),
  lapply(names(expo), function(xv) run_linda(GS, paste("~", xv, "+", COV), xv, rows_C, md = md_C) %>%
    filter(Taxon %in% srb) %>% mutate(Exposure = expo[[xv]], Group = "IBS-C"))) %>%
  relocate(Group, Exposure, Taxon) %>% rename(log2FC_per_rank_unit = log2FC)

# Causal mediation (Imai et al.; nonparametric bootstrap, 5000 resamples).
# Main analysis: rank-based - exposure, mediator and outcome as percentile ranks (rank/n) within the group;
#          contrast = 25th vs 75th percentile of the exposure; effects in percentile points of the outcome.
# Sensitivity: original scales (exposure 25th vs 75th percentile value).
# Covariates as numeric 0/1 indicators (Male, IBS-D) so that no bootstrap resample fails on a single-level factor.
COVn <- "Age + Male + BMI + PPI"
cache_f <- "output/mediation_cache.rds"
cache <- if (file.exists(cache_f)) readRDS(cache_f) else list()
med_one <- function(g, xv, y, scale = c("rank", "original")) {
  scale <- match.arg(scale); key <- paste(if (scale == "rank") "rank" else NULL, g, xv, y, sep = "|")
  key <- sub("^\\|", "", key)
  if (!is.null(cache[[key]])) return(cache[[key]])
  dd <- grp[[g]]; ex <- if (g == "All") " + IBS_D" else ""
  if (scale == "rank") { n <- nrow(dd); dd$X <- rank(dd[[xv]]) / n; dd$M <- rank(dd$Desulfo_log) / n; dd$Y <- rank(dd[[y]]) / n }
  else { dd$X <- dd[[xv]]; dd$M <- dd$Desulfo_log; dd$Y <- dd[[y]] }
  # do.call puts the formula and data into the model call, so mediate() can refit it in its bootstrap
  mf <- do.call("lm", list(formula = as.formula(paste0("M ~ X + ", COVn, ex)), data = dd))
  of <- do.call("lm", list(formula = as.formula(paste0("Y ~ M + X + ", COVn, ex)), data = dd))
  lo <- unname(quantile(dd$X, .25)); hi <- unname(quantile(dd$X, .75))
  if (hi == lo) { lo <- min(dd$X); hi <- max(dd$X) }
  mo <- tryCatch(mediate(mf, of, treat = "X", mediator = "M", control.value = lo, treat.value = hi, boot = TRUE, sims = 5000),
                 error = function(e) e)
  if (inherits(mo, "error")) return(data.frame(Scale = scale, Group = g, Exposure = expo[[xv]], Outcome = outc[[y]], Note = conditionMessage(mo)))
  sens <- tryCatch(medsens(mediate(mf, of, treat = "X", mediator = "M", control.value = lo, treat.value = hi, sims = 1000),
                           rho.by = 0.05, sims = 500), error = function(e) NULL)
  out <- data.frame(Scale = scale, Group = g, Exposure = expo[[xv]], Outcome = outc[[y]], n = nrow(dd),
             Contrast = if (scale == "rank") "75th vs 25th percentile" else sprintf("%.1f vs %.1f", hi, lo),
             ACME = mo$d0, ACME_CI_low = mo$d0.ci[1], ACME_CI_high = mo$d0.ci[2], ACME_p = mo$d0.p,
             ADE = mo$z0, ADE_CI_low = mo$z0.ci[1], ADE_CI_high = mo$z0.ci[2], ADE_p = mo$z0.p,
             Total = mo$tau.coef, Total_CI_low = mo$tau.ci[1], Total_CI_high = mo$tau.ci[2], Total_p = mo$tau.p,
             Prop_mediated = mo$n0, Prop_med_CI_low = mo$n0.ci[1], Prop_med_CI_high = mo$n0.ci[2],
             Sensitivity_rho_at_ACME0 = if (is.null(sens)) NA else sens$err.cr.d[1])
  cache[[key]] <<- out; saveRDS(cache, cache_f); out
}
med_all <- function(scale) {
  out <- bind_rows(lapply(c("All", "IBS-C"), function(g) bind_rows(lapply(c("FFQ_processed_meat_score", "Meat_RP", "Western_score"),
    function(xv) bind_rows(lapply(c("IBSSSS_Q1_pain_intensity", "IBSSSS_Q2_pain_days", "IBSSSS_TOTAL"), function(y) med_one(g, xv, y, scale)))))))
  out$Scale <- scale; rownames(out) <- NULL; out
}
# the bootstrap uses the random number stream; the stream is restored afterwards so that later permutation tests
# give the same P values whether the mediation results are computed afresh or read from the cache
rng_before_mediation <- .Random.seed
R$E7_mediation <- med_all("rank")
R$E7s_mediation_linear <- med_all("original")
assign(".Random.seed", rng_before_mediation, envir = .GlobalEnv)

# E-values (VanderWeele & Ding 2017) for the significant main (rank-based) associations:
# partial rho -> d = 2r / sqrt(1 - r^2) -> RR ~ exp(0.91 d); bound from the CI limit closest to the null.
ev <- function(r) { if (is.na(r) || r == 0) return(1); dd <- 2 * abs(r) / sqrt(1 - r^2); rr <- exp(0.91 * dd); rr + sqrt(rr * (rr - 1)) }
R$E8_E_values <- chain_rank %>% filter(p < 0.05) %>% rowwise() %>%
  mutate(lim = if (partial_rho > 0) CI_low else CI_high,
         E_value_point = ev(partial_rho), E_value_CI = if (sign(lim) == sign(partial_rho)) ev(lim) else 1) %>%
  ungroup() %>% dplyr::select(Path, Group, Predictor, Outcome, partial_rho, p, q_chain_within_group, E_value_point, E_value_CI)

# =========================================================================================
# F. EXPLORATORY, CONTEXT, POWER
# =========================================================================================
firth <- function(f, data, term) { m <- glm(f, family = binomial, data = data, method = "brglmFit", type = "AS_mean")
  ci <- confint.default(m); data.frame(OR = exp(coef(m)[[term]]), CI_low = exp(ci[term, 1]), CI_high = exp(ci[term, 2]),
                                       p = summary(m)$coefficients[term, 4]) }
mb <- "Methanobrevibacter spp."
R$F1_methanogens <- bind_rows(
  cbind(Analysis = "Methanobrevibacter detected: IBS-D vs IBS-C (Firth, adjusted)",
        Detail = sprintf("%d/%d IBS-C vs %d/%d IBS-D; Fisher p = %.3f", sum(d$Methano_pos[d$Subtype == "IBS-C"]), sum(d$Subtype == "IBS-C"),
                         sum(d$Methano_pos[d$Subtype == "IBS-D"]), sum(d$Subtype == "IBS-D"), fisher.test(table(d$Subtype, d$Methano_pos))$p.value),
        firth(as.formula(paste("Methano_pos ~ Subtype +", COV)), d, "SubtypeIBS-D") %>% rename(Estimate = OR)),
  cbind(Analysis = "Bowel movements per week by Methanobrevibacter detection (adjusted, + subtype)", Detail = "linear model, HC3",
        hc(lm(as.formula(paste("BM_week ~ Methano_pos + Subtype +", COV)), d), "Methano_pos") %>% rename(Estimate = Beta)),
  cbind(Analysis = "Desulfovibrio (log10) by Methanobrevibacter detection (adjusted, + subtype)", Detail = "linear model, HC3",
        hc(lm(as.formula(paste("Desulfo_log ~ Methano_pos + Subtype +", COV)), d), "Methano_pos") %>% rename(Estimate = Beta)))
for (g in names(grp)) { ct <- cor.test(grp[[g]]$`Desulfovibrio spp.`, grp[[g]][[mb]], method = "spearman", exact = FALSE)
  R$F1_methanogens <- bind_rows(R$F1_methanogens, data.frame(Analysis = paste0("Desulfovibrio vs Methanobrevibacter, Spearman (", g, ")"),
    Detail = "unadjusted", Estimate = unname(ct$estimate), p = ct$p.value)) }

R$F2_anxiety_DA <- bind_rows(
  run_linda(mat(genus), paste("~ Anxiety_dx + Subtype +", COV), "Anxiety_dx") %>% mutate(Level = "Genus"),
  run_linda(mat(species), paste("~ Anxiety_dx + Subtype +", COV), "Anxiety_dx") %>% mutate(Level = "Species")) %>%
  relocate(Level) %>% arrange(p)
R$F2b_anxiety_community <- bind_rows(
  cbind(Analysis = "Shannon (genus) by anxiety, adjusted + subtype", hc(lm(as.formula(paste("Shannon_genus ~ Anxiety_dx + Subtype +", COV)), d), "Anxiety_dx")),
  perm(Dg, paste("Anxiety_dx + Subtype +", COV), "Aitchison PERMANOVA, anxiety (marginal)") %>% filter(Term == "Anxiety_dx") %>%
    transmute(Analysis = Analysis, Beta = R2, p = p))

R$F3_PPI_context <- bind_rows(
  data.frame(Analysis = "PPI use by subtype", Detail = sprintf("%d/%d IBS-C vs %d/%d IBS-D", sum(d$PPI[d$Subtype == "IBS-C"]), sum(d$Subtype == "IBS-C"),
             sum(d$PPI[d$Subtype == "IBS-D"]), sum(d$Subtype == "IBS-D")), Estimate = NA, p = fisher.test(table(d$Subtype, d$PPI))$p.value),
  data.frame(Analysis = "Upper-GI comorbidity by subtype", Detail = sprintf("%d/%d IBS-C vs %d/%d IBS-D", sum(d$UpperGI[d$Subtype == "IBS-C"]), sum(d$Subtype == "IBS-C"),
             sum(d$UpperGI[d$Subtype == "IBS-D"]), sum(d$Subtype == "IBS-D")), Estimate = NA, p = fisher.test(table(d$Subtype, d$UpperGI))$p.value),
  cbind(Analysis = "PPI use: IBS-D vs IBS-C, adjusted for upper-GI comorbidity, age, sex, BMI (Firth OR)", Detail = "",
        firth(PPI ~ Subtype + UpperGI + Age + Sex + BMI, d, "SubtypeIBS-D") %>% rename(Estimate = OR)))

pw_d <- uniroot(function(dd) { n1 <- 35; n2 <- 29; df <- n1 + n2 - 2; ncp <- dd * sqrt(n1 * n2 / (n1 + n2))
  1 - pt(qt(.975, df), df, ncp) + pt(-qt(.975, df), df, ncp) - 0.8 }, c(0.1, 2))$root
r_detect <- function(n) tanh((qnorm(.975) + qnorm(.8)) / sqrt(n - 3))
R$F4_sensitivity_power <- data.frame(
  Comparison = c("IBS-C (n = 35) vs IBS-D (n = 29), two-sided alpha 0.05, 80% power",
                 "Correlation, all patients (n = 64)", "Correlation, IBS-C (n = 35)", "Correlation, IBS-D (n = 29)"),
  Minimum_detectable_effect = c(pw_d, r_detect(64), r_detect(35), r_detect(29)),
  Metric = c("standardised mean difference (Cohen's d)", "|r|", "|r|", "|r|"))

saveRDS(list(R = R, ord = ord, ord_var = ord_var, d = d), "output/res_11.rds")
for (n in names(R)) write.csv(R[[n]], file.path("output", paste0(n, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
cat("Done:", length(R), "result tables\n")
