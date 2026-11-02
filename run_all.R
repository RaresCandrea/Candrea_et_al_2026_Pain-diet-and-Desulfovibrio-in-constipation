# Run the full pipeline: Rscript run_all.R [path/to/database.xlsx]
# Requires R >= 4.4 with the packages listed in output/sessionInfo.txt, and Python 3 with
# openpyxl, python-docx, pandas, scipy and statsmodels (tables, workbook, independent checks).
Sys.setlocale("LC_ALL", "C.UTF-8")
args <- commandArgs(trailingOnly = TRUE)
xlsx <- if (length(args)) args[1] else "data/PUBLIC_FINAL_Baza_de_date_microbioame_EN.xlsx"
run <- function(cmd, f, a = character()) { cat("\n====", f, "====\n"); stopifnot(system2(cmd, c(f, a)) == 0) }
run("Rscript", "01_data_prep.R", shQuote(xlsx))          # import, exclusions, derived variables
dir.create("output", showWarnings = FALSE); dir.create("tables/supp", recursive = TRUE, showWarnings = FALSE)
run("Rscript", "08_tables_1_2.R")                         # Table 1 and Table 2 values
run("Rscript", "11_main_analyses.R")                      # all analyses -> output/*.csv
# Word versions of Tables 1-2 require the journal table templates (not distributed); verify them only if present
if (file.exists("tables/Table1_clinical_characteristics.docx")) {
  run("python3", "10_verify_tables_independent.py", c(shQuote(xlsx), "tables"))
} else cat("\n(Tables 1-2 Word files not present: independent check of Word tables skipped; values are in tables/T1_T2_values.csv)\n")
run("python3", "13_verify_results_independent.py", shQuote(xlsx))
suppressPackageStartupMessages({ library(vegan); library(MicrobiomeStat); library(brglm2); library(mediation)
  library(sandwich); library(lmtest); library(dplyr); library(readxl) })
writeLines(capture.output(sessionInfo()), "output/sessionInfo.txt")
run("python3", "12_build_results_workbook.py", c(shQuote(basename(xlsx)), "IBS_subtypes_results.xlsx"))
# Figures (main text) and supplementary tables
run("Rscript", "14_figure2_diversity.R")                  # Figure 2: alpha and beta diversity
run("Rscript", "17_figure3_desulfovibrio.R")              # Figure 3: Desulfovibrio and genus-level differential abundance
run("Rscript", "18_figure4_taxa_symptoms_heatmap.R")      # Figure 4: key taxa x symptoms (also writes output/C3_taxa_x_symptoms.csv)
run("Rscript", "19_figure5_diet_Desulfovibrio_pain.R")    # Figure 5: diet, Desulfovibrio and abdominal pain
run("Rscript", "20_figure6_patient_genus_profiles.R")     # Figure 6: genus-level profiles of individual patients
run("Rscript", "21_figure7_lifestyle_composition.R")      # Figure 7: clinical, lifestyle and dietary factors and community composition
run("Rscript", "15_supp_tables.R")                        # Supplementary Tables 1-7 (JSON specs)
for (s in sprintf("tables/supp/S%d.json", 1:7)) run("node", "16_build_supp_docx.js", c(s, "tables/supp"))
