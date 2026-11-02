# Abdominal pain, Western diet and *Desulfovibrio* in IBS-C and IBS-D — analysis code

Analysis code for the manuscript:

> **Abdominal Pain Is Associated With Western Diet and *Desulfovibrio* in Constipation- but Not Diarrhea-predominant Irritable Bowel Syndrome**
> Candrea AR, Gavrilaș LI, Mocan A, Lázaro Á, Bogi B, Crișan G. *Journal of Neurogastroenterology and Motility* (submitted).

Retrospective cross-sectional study of adults with irritable bowel syndrome (IBS) who underwent gut microbiota
analysis at a nutrition and dietetics clinic in Cluj-Napoca, Romania (May 2021 – October 2026):
64 patients analysed, 35 with constipation-predominant IBS (IBS-C) and 29 with diarrhea-predominant IBS (IBS-D).
Gut microbiota data are relative abundances from a laboratory 16S rRNA-based targeted panel (raw reads not available).

## Data availability
The patient-level database is **not** included in this repository. The anonymised dataset is available from the
corresponding author on reasonable request. Place it at `data/PUBLIC_FINAL_Baza_de_date_microbioame_EN.xlsx`
(or pass its path to `run_all.R`) to reproduce all results.

## How to run
```
Rscript run_all.R data/<database>.xlsx
python3 figure1/make_fig1_pptx.py          # Figure 1 (flow diagram, editable PowerPoint)
```
Outputs: `output/*.csv` (every result table), `tables/` (Table 1–2 values, Supplementary Tables 1–7),
`figures/` (Figures 2–7 as TIFF 600 dpi, PNG and vector PDF).
Editable PowerPoint versions of Figures 2–7: run the figure scripts in R with
`options(fig.pptx = TRUE, fig.font = "Arial")`.

## Scripts
| Script | Content |
|---|---|
| `01_data_prep.R` | import (column-name based), exclusions, questionnaire totals, derived variables (panel Shannon index, dietary scores, healthy lifestyle score) |
| `08_tables_1_2.R` | values for Table 1 (clinical characteristics) and Table 2 (lifestyle and diet) |
| `09_fill_tables_docx.py` | writes Table 1–2 values into the journal Word templates (templates not distributed) |
| `10_verify_tables_independent.py` | independent re-computation of every Table 1–2 cell (Python) |
| `11_main_analyses.R` | all analyses: diversity (Shannon, Aitchison PERMANOVA), differential abundance (LinDA, CLR), taxa–symptom and diet–*Desulfovibrio*–symptom associations (adjusted Spearman, interaction tests, causal mediation, E-values), lifestyle factors and community composition, exploratory analyses, minimum detectable effects |
| `12_build_results_workbook.py` | collects all result tables into one Excel workbook |
| `13_verify_results_independent.py` | independent re-computation of the key results (Python: SciPy, statsmodels) |
| `14_figure2_diversity.R` | Figure 2: alpha diversity and Aitchison PCoA |
| `17_figure3_desulfovibrio.R` | Figure 3: *Desulfovibrio* abundance and genus-level LinDA effects |
| `18_figure4_taxa_symptoms_heatmap.R` | Figure 4: key taxa × symptoms (adjusted Spearman) |
| `19_figure5_diet_Desulfovibrio_pain.R` | Figure 5: diet, *Desulfovibrio* and abdominal pain by subtype |
| `20_figure6_patient_genus_profiles.R` | Figure 6: genus-level profiles of individual patients |
| `21_figure7_lifestyle_composition.R` | Figure 7: clinical, lifestyle and dietary factors and community composition |
| `15_supp_tables.R`, `16_build_supp_docx.js` | Supplementary Tables 1–7 (JSON specifications rendered to Word) |
| `fig_save.R` | common figure export |
| `figure1/make_fig1_pptx.py` | Figure 1 (patient selection flow diagram) |
| `run_all.R` | runs the full pipeline |

## Statistical methods (summary)
Adjusted models include age, sex, body mass index and proton pump inhibitor use (12 months preceding the clinical
assessment); heteroscedasticity-consistent (HC3) standard errors; Benjamini–Hochberg false discovery rate within
analysis families (*q* < 0.10). Ordinal dietary scores are analysed as ranks. Random seed 20261031.
Full details are given in the manuscript (Materials and Methods).

**Mediation results.** The causal mediation estimates (bootstrap, 5000 resamples) reported in the manuscript are stored
in `output/mediation_cache.rds` (summary estimates only, no patient-level data); `11_main_analyses.R` reads them from
this file so that the published values are reproduced exactly. If the file is deleted, the mediation analyses are
recomputed; bootstrap confidence intervals and *P* values may then differ slightly through Monte Carlo variation.

## Software
R 4.5.3 — vegan 2.7-5, MicrobiomeStat 1.4, sandwich 3.1-3, lmtest 0.9-40, brglm2 1.1.0, mediation 4.5.1,
ggplot2 4.0.3, patchwork, dplyr 1.2.1, readxl 1.6.0, jsonlite, ragg 1.5.2, officer 0.7.7, rvg 0.4.2
(complete list in `output/sessionInfo.txt`).
Python 3 — pandas 3.0.5, SciPy 1.18.1, statsmodels 0.15.0, openpyxl, python-docx, python-pptx.
Node.js — docx (Supplementary Tables in Word).

## License
MIT (see `LICENSE`).
