# =============================================================================
# 08 - Values for Table 1 (clinical characteristics) and Table 2 (lifestyle/diet)
# Continuous: median (Q1, Q3), Mann-Whitney U (normal approximation, continuity
# correction). Binary: Fisher's exact test. >2 categories: Fisher's exact test with
# Monte Carlo p value (100 000 replicates). p values descriptive, not adjusted.
# =============================================================================
suppressPackageStartupMessages(library(dplyr))
set.seed(20261031)
x <- readRDS("output/analysis_data.rds"); d <- x$d
g <- d$Subtype; iC <- g == "IBS-C"; iD <- g == "IBS-D"

# ---- Derived variables as defined in the manuscript ----------------------------------
d$Plant_score4 <- d$FFQ_fruits_score + d$FFQ_vegetables_score + d$FFQ_legumes_score + d$FFQ_whole_grains_score
d$Western_score <- d$FFQ_red_meat_score + d$FFQ_processed_meat_score + d$FFQ_processed_sweets_score
diet_quality <- d$Plant_score4 - d$Western_score
d$HLS_ms <- (d$Smoking == 0) + (d$PSQI_TOTAL <= 5) + (d$PA_active == 1) +
            (diet_quality > median(diet_quality)) + (d$FFQ_alcohol_score <= 1.5)
d$Female <- as.integer(d$Sex == "F")
d$BMI_cat <- cut(d$BMI, c(-Inf, 18.5, 25, 30, Inf), right = FALSE,
                 labels = c("Underweight", "Normal", "Overweight", "Obesity"))
d$VitD_cat <- cut(d$VitD, c(-Inf, 20, 30, Inf), right = FALSE, labels = c("Deficiency", "Insufficiency", "Sufficiency"))
d$SSS_cat <- droplevels(d$IBSSSS_cat)

# ---- Formatting helpers --------------------------------------------------------------
fp <- function(p) if (is.na(p)) "" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
mq <- function(z, dec) { q <- quantile(z, c(.5, .25, .75), type = 7)
  sprintf(paste0("%.", dec, "f (%.", dec, "f, %.", dec, "f)"), q[1], q[2], q[3]) }
np <- function(k, n) sprintf("%d (%.1f%%)", k, 100 * k / n)

rows <- list()
add <- function(tab, label, all, c, dd, p = NA, ptxt = NULL)
  rows[[length(rows) + 1]] <<- data.frame(Table = tab, Label = label, All = all, IBS_C = c, IBS_D = dd,
                                          p = if (is.null(ptxt)) fp(p) else ptxt, p_num = p)
cont <- function(tab, label, v, dec) { z <- d[[v]]
  add(tab, label, mq(z, dec), mq(z[iC], dec), mq(z[iD], dec), wilcox.test(z[iD], z[iC], exact = FALSE)$p.value) }
bin <- function(tab, label, v) { z <- d[[v]]
  add(tab, label, np(sum(z), length(z)), np(sum(z[iC]), sum(iC)), np(sum(z[iD]), sum(iD)),
      fisher.test(table(g, factor(z, 0:1)))$p.value) }
catv <- function(tab, header, v, labels, test = TRUE) { z <- d[[v]]
  p <- if (test) fisher.test(table(g, z), simulate.p.value = TRUE, B = 1e5)$p.value else NA
  add(tab, header, "", "", "", p, if (!test) "—" else NULL)
  for (k in seq_along(levels(z))) { lv <- levels(z)[k]
    add(tab, labels[k], np(sum(z == lv), length(z)), np(sum(z[iC] == lv), sum(iC)), np(sum(z[iD] == lv), sum(iD)), NA) } }

# ---- Table 1 ------------------------------------------------------------------------------
T <- "T1"
bin(T, "Female sex, n (%)", "Female")
cont(T, "Age, years", "Age", 0)
cont(T, "BMI, kg/m²", "BMI", 1)
catv(T, "BMI category, n (%)", "BMI_cat", c("Underweight (<18.5)", "Normal weight (18.5–24.9)", "Overweight (25.0–29.9)", "Obesity (≥30.0)"))
d$Bristol_f <- factor(d$Bristol, levels = c(1, 2, 6, 7))
catv(T, "Bristol Stool Form Scale, n (%)", "Bristol_f", paste("Type", c(1, 2, 6, 7)), test = FALSE)
cont(T, "Bowel movements per week", "BM_week", 0)
cont(T, "IBS-SSS total score", "IBSSSS_TOTAL", 0)
cont(T, "Abdominal pain intensity", "IBSSSS_Q1_pain_intensity", 0)
cont(T, "Days with pain (x10)", "IBSSSS_Q2_pain_days", 0)
cont(T, "Abdominal distension", "IBSSSS_Q3_distension", 0)
cont(T, "Dissatisfaction with bowel habit", "IBSSSS_Q4_bowel_dissatisfaction", 0)
cont(T, "Interference with daily life", "IBSSSS_Q5_interference", 0)
catv(T, "IBS-SSS severity, n (%)", "SSS_cat", levels(d$SSS_cat))
bin(T, "Anxiety disorder", "Anxiety_dx")
bin(T, "Depression", "Depression_dx")
bin(T, "Autoimmune disease", "Autoimmune_disease_dx")
bin(T, "Other comorbidity", "Other_comorbidity_dx")
bin(T, "Upper gastrointestinal comorbidity", "UpperGI")
bin(T, "Antibiotics", "Antibiotics_12mo")
bin(T, "Proton pump inhibitors", "PPI_12mo_prior")
bin(T, "Laxatives", "Laxatives_12mo_prior")
bin(T, "Antidiarrheals", "Antidiarrheals_12mo_prior")
bin(T, "Antidepressants", "Antidepressants_12mo_prior")
cont(T, "Fasting glucose, mg/dL", "Glucose", 0)
cont(T, "Total cholesterol, mg/dL", "TC", 0)
cont(T, "HDL cholesterol, mg/dL", "HDL", 0)
cont(T, "LDL cholesterol, mg/dL", "LDL", 0)
cont(T, "Triglycerides, mg/dL", "TG", 0)
cont(T, "25(OH) vitamin D, ng/mL", "VitD", 1)
catv(T, "Vitamin D status, n (%)", "VitD_cat", c("Deficiency (<20)", "Insufficiency (20–29.9)", "Sufficiency (≥30)"))

# ---- Table 2 ------------------------------------------------------------------------------
T <- "T2"
bin(T, "Current smoking, n (%)", "Smoking")
catv(T, "Physical activity, min/week, n (%)", "PA_cat", c("<150", "150–299", "≥300"))
cont(T, "PSQI global score", "PSQI_TOTAL", 0)
bin(T, "Poor sleep quality (PSQI > 5), n (%)", "PSQI_poor")
cont(T, "Alcohol intake (frequency score)", "FFQ_alcohol_score", 1)
bin(T, "Coffee consumption, n (%)", "FFQ_coffee_YES_NO")
cont(T, "Water intake, L/day", "FFQ_water_L_day", 1)
cont(T, "Healthy lifestyle score (0–5)", "HLS_ms", 0)
foods <- c(Fruits = "FFQ_fruits_score", Vegetables = "FFQ_vegetables_score", Legumes = "FFQ_legumes_score",
           `Whole grains` = "FFQ_whole_grains_score", `Red meat` = "FFQ_red_meat_score",
           `Processed meat` = "FFQ_processed_meat_score", `White meat` = "FFQ_white_meat_score",
           Fish = "FFQ_fish_score", `Dairy products` = "FFQ_dairy_score", Eggs = "FFQ_eggs_score",
           `Processed sweets` = "FFQ_processed_sweets_score")
for (f in names(foods)) cont(T, f, foods[[f]], 1)
cont(T, "Plant food score (0–28)", "Plant_score4", 1)
cont(T, "Western food score (0–21)", "Western_score", 1)

out <- bind_rows(rows)
dir.create("tables", showWarnings = FALSE)
write.csv(out, "tables/T1_T2_values.csv", row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(d, "output/analysis_data_tables.rds")
options(width = 200); print(out[, 1:6], row.names = FALSE)
