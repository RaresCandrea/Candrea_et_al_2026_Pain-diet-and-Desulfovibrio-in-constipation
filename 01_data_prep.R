# =============================================================================
# IBS subtypes (IBS-C vs IBS-D), lifestyle and gut microbiota
# 01 - Data import, cleaning, derived variables (column-name based)
# =============================================================================
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
})
set.seed(20261031)

args <- commandArgs(trailingOnly = TRUE)
xlsx <- if (length(args)) args[1] else "data/PUBLIC_FINAL_Baza_de_date_microbioame_EN.xlsx"
dir.create("output", showWarnings = FALSE)

raw <- suppressMessages(read_excel(xlsx, col_names = FALSE, .name_repair = "minimal"))
h1  <- as.character(unlist(raw[1, ])); h2 <- as.character(unlist(raw[2, ]))
# Row 1 holds block labels for the laboratory section and variable names for the clinical section;
# row 2 holds variable names for the laboratory section and units for the clinical section.
clin_start <- which(h1 == "IBS_subtype")
nm <- ifelse(seq_along(h1) >= clin_start, h1, h2); nm[1] <- "Subject"
nm <- trimws(gsub("\\*+$", "", nm))
dat <- raw[-(1:2), ]; names(dat) <- make.unique(nm)
dat <- dat[!is.na(dat$Subject), ]

# Block boundaries from row-1 labels
blk <- function(lab) which(h1 == lab)
phyla   <- nm[blk("PHYLUM LEVEL"):(blk("GENUS LEVEL") - 1)]
taxa_gs <- nm[blk("GENUS LEVEL"):(which(nm == "Moulds") - 1)]
taxa_gs <- setdiff(taxa_gs, "Methanobacteria")      # class-level duplicate of Methanobrevibacter
bioind  <- nm[blk("BIOINDICATORS"):(blk("PHYLUM LEVEL") - 1)]

# Numeric conversion: every column whose non-missing values all parse as numbers
for (v in names(dat)) {
  x <- dat[[v]]; xn <- suppressWarnings(as.numeric(as.character(x)))
  if (sum(!is.na(x)) > 0 && sum(is.na(xn)) == sum(is.na(x))) dat[[v]] <- xn
}

dat <- dat %>%
  rename(Exclusion = `Rasons of exclusion`,
         PA = `Physical activity level (mins/week)`,
         Bristol = `Bristol Stool Scale`,
         BM_week = `Weekly Bowel Movement Frequency`,
         Age = `Age (years)`, BMI = `BMI (kg/m2)`, Weight = `Weight (kg)`,
         Glucose = `Serum glucose (mg/dL)`, HDL = `HDL cholesterol (mg/dL)`,
         TC = `Total cholesterol (mg/dL)`, TG = `Triglycerides (mg/dL)`,
         LDL = `LDL cholesterol (mg/dL)`, VitD = `Vitamin D (ng/mL)`) %>%
  mutate(
    Subtype  = factor(IBS_subtype, levels = c("C", "D"), labels = c("IBS-C", "IBS-D")),
    Sex      = factor(Sex, levels = c("F", "M")),
    PA_cat   = factor(case_when(grepl("<\\s*150", PA) ~ "<150",
                                grepl("150", PA) ~ "150-299",
                                grepl("300", PA) ~ ">=300"),
                      levels = c("<150", "150-299", ">=300")),
    PA_active = as.integer(PA_cat != "<150"),            # meets WHO >=150 min/week
    Antibiotics_12mo = as.integer(grepl("antibiotic", Other_medications_12mo_prior, ignore.case = TRUE)),
    Moulds_pos = as.integer(Moulds == "pozitiv")
  )
yn <- c("PPI_12mo_prior", "Laxatives_12mo_prior", "Antidiarrheals_12mo_prior",
        "Antidepressants_12mo_prior", "Smoking", "Anxiety_dx", "Depression_dx",
        "Autoimmune_disease_dx", "Other_comorbidity_dx", "FFQ_coffee_YES_NO")
dat[yn] <- lapply(dat[yn], function(x) as.integer(toupper(x) == "YES"))

# Recompute questionnaire totals (workbook formulas sum partially missing items)
dat <- dat %>% mutate(
  IBSSSS_TOTAL = IBSSSS_Q1_pain_intensity + IBSSSS_Q2_pain_days + IBSSSS_Q3_distension +
                 IBSSSS_Q4_bowel_dissatisfaction + IBSSSS_Q5_interference,
  IBSSSS_cat = cut(IBSSSS_TOTAL, c(-Inf, 74, 174, 300, Inf), labels = c("Remission", "Mild", "Moderate", "Severe")),
  PSQI_TOTAL = PSQI_C1_subjective_quality + PSQI_C2_latency + PSQI_C3_duration +
               PSQI_C4_efficiency + PSQI_C5_disturbances + PSQI_C6_sleep_medication +
               PSQI_C7_daytime_dysfunction,
  PSQI_poor = as.integer(PSQI_TOTAL > 5)
)

# ---- Analytic sample ----------------------------------------------------------------
excluded <- dat %>% filter(!is.na(Exclusion)) %>% select(Subject, Subtype, Exclusion)
d <- dat %>% filter(is.na(Exclusion))
d$D   <- as.integer(d$Subtype == "IBS-D")
d$PPI <- d$PPI_12mo_prior

# ---- Pre-defined derived variables -------------------------------------------------------
genus   <- grep("spp\\.$", taxa_gs, value = TRUE)
species <- setdiff(taxa_gs, genus)
shannon <- function(m) { p <- m / rowSums(m); -rowSums(ifelse(p > 0, p * log(p), 0)) }
logp    <- function(v) log10(v + min(v[v > 0]) / 2)   # log10 with half-minimum pseudo-count

d <- d %>% mutate(
  # Diversity of the targeted panel (computed from relative abundances)
  Shannon_genus  = shannon(as.matrix(d[, genus])),
  Shannon_phylum = shannon(as.matrix(d[, phyla])),
  Richness_genus = rowSums(d[, genus] > 0),
  # Dietary scores (FFQ frequency scores)
  Plant_score   = FFQ_fruits_score + FFQ_vegetables_score + FFQ_legumes_score + FFQ_whole_grains_score,
  Western_score = FFQ_red_meat_score + FFQ_processed_meat_score + FFQ_processed_sweets_score,
  Meat_RP       = FFQ_red_meat_score + FFQ_processed_meat_score,
  # Hydrogenotrophic guild
  Desulfo_log   = logp(`Desulfovibrio spp.`),
  Bilophila_log = logp(`Bilophila wadsworthia`),
  Methano_pos   = as.integer(`Methanobrevibacter spp.` > 0),
  PB_ratio      = log10((`Prevotella spp.` + 0.001) / (`Bacteroides spp.` + 0.001)),
  # Upper-GI comorbidity (free-text)
  UpperGI = as.integer(grepl("gastritis|reflux|GERD|hiatal|Barrett|dyspepsia|duoden|ulcer|Helicobacter",
                             Comorbidity_details, ignore.case = TRUE))
)
# Healthy lifestyle score adapted from Ho et al. (Gut 2024): one point per factor
# (non-smoking, PSQI <= 5, physical activity >= 150 min/week, diet quality [plant - western]
#  above the sample median, alcohol frequency score <= 1.5) - as in Table 2
d <- d %>% mutate(Diet_quality = Plant_score - Western_score)
d <- d %>% mutate(
  HLS = (Smoking == 0) + (PSQI_TOTAL <= 5) + (PA_active == 1) +
        (Diet_quality > median(Diet_quality)) + (FFQ_alcohol_score <= 1.5)
)

ffq     <- grep("^FFQ_.*score$", names(d), value = TRUE)
ffq_all <- c(ffq, "FFQ_water_L_day")
serum   <- c("Glucose", "HDL", "TC", "TG", "LDL", "VitD")

saveRDS(list(d = d, all = dat, excluded = excluded, bioind = bioind, phyla = phyla,
             taxa = taxa_gs, genus = genus, species = species,
             ffq = ffq, ffq_all = ffq_all, serum = serum),
        "output/analysis_data.rds")
cat("Screened:", nrow(dat), " Excluded:", nrow(excluded), " Analysed:", nrow(d), "\n")
print(table(d$Subtype))
cat("Phyla:", length(phyla), " Genus:", length(genus), " Species:", length(species), " Bioindicators:", length(bioind), "\n")
