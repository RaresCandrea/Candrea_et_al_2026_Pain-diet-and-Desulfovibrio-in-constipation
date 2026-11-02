"""Collect every result table (output/*.csv from 08 and 11) into one formatted workbook.

Values are outputs of the R analyses (not spreadsheet formulas). Usage:
    python 12_build_results_workbook.py <database_name> <out.xlsx>
"""
import csv
import sys
from datetime import date

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

db_name, out = sys.argv[1], sys.argv[2]
FONT = "Arial"
HEAD = PatternFill("solid", start_color="D9E2F3")
SIG = PatternFill("solid", start_color="E2EFDA")     # q < 0.05 (or p < 0.05 where no q exists)
TREND = PatternFill("solid", start_color="FFF2CC")   # 0.05 <= q < 0.10
thin = Side(style="thin", color="BFBFBF")

SHEETS = [  # (sheet name, csv file, family, what it is)
    ("T1_clinical", "tables/T1_T2_values.csv", "Descriptive",
     "Table 1 values (clinical characteristics). Median (Q1, Q3) + Mann-Whitney U; n (%) + Fisher's exact test (Monte Carlo for >2 categories). Not adjusted for multiple comparisons."),
    ("T2_lifestyle", "tables/T1_T2_values.csv", "Descriptive",
     "Table 2 values (lifestyle and diet). Same methods as Table 1."),
    ("A1_alpha", "output/A1_alpha_diversity.csv", "Secondary - community",
     "Shannon index and richness computed from the relative abundances of the targeted panel. Mann-Whitney U (crude) and linear model (HC3) adjusted for age, sex, BMI, PPI."),
    ("A2_PERMANOVA", "output/A2_PERMANOVA.csv", "Secondary - community",
     "PERMANOVA (adonis2, 9999 permutations, marginal terms) on Aitchison distance (CLR, taxa with prevalence >= 10%); Bray-Curtis as sensitivity. Dispersion = betadisper permutation test."),
    ("A3_composition", "output/A3_composition.csv", "Descriptive",
     "Mean relative abundance (%) and prevalence (proportion > 0) of every panel taxon by subtype (for the composition figure)."),
    ("B1_DA_subtype", "output/B1_DA_subtype.csv", "Main result (Desulfovibrio) / other taxa",
     "Differential abundance IBS-D vs IBS-C for taxa with prevalence >= 10%: crude Wilcoxon; LinDA (relative abundances, adjusted for age, sex, BMI, PPI; log2FC < 0 = lower in IBS-D); CLR + linear model (HC3) as second method. q = BH within taxonomic level."),
    ("B2_Desulfo_sens", "output/B2_Desulfovibrio_sensitivity.csv", "Main result - sensitivity",
     "Desulfovibrio spp. under alternative models and samples."),
    ("C1_functional", "output/C1_functional_by_subtype.csv", "Secondary - functional",
     "Laboratory functional indicators by subtype: crude Mann-Whitney U and adjusted linear model (HC3; skewed indicators on log(1+x)). q = BH across the 7 indicators."),
    ("C2_taxa_x_symptoms", "output/C2_taxa_indicators_x_symptoms.csv", "Secondary - symptoms",
     "Partial Spearman correlations (adjusted for age, sex, BMI, PPI; + subtype in 'All') between selected key taxa / functional indicators and IBS-SSS items, IBS-SSS total and weekly bowel movements. q = BH within group."),
    ("D1_lifestyle_PERMANOVA", "output/D1_lifestyle_PERMANOVA.csv", "Secondary - lifestyle",
     "Variance in community composition (Aitchison, genus panel) explained by each lifestyle/diet factor, within each subtype (PERMANOVA, 4999 permutations, adjusted for age and sex; + subtype in 'All'; dietary scores as ranks). q = BH within group."),
    ("E1_dietToDesulfo", "output/E1_path_A_diet_to_Desulfovibrio.csv", "Secondary - diet chain (main method)",
     "Path A: dietary exposure -> Desulfovibrio. Partial Spearman correlation (ranks; adjusted for age, sex, BMI, PPI; + subtype in 'All'), Fisher-z 95% CI. FFQ scores are ordinal, hence the rank-based method. q = BH across the whole chain (paths A + B + C) within each group."),
    ("E2_DesulfoToSympt", "output/E2_path_B_Desulfovibrio_to_symptoms.csv", "Secondary - diet chain (main method)",
     "Path B: Desulfovibrio -> symptoms (IBS-SSS items and total). Same method and FDR family as E1."),
    ("E3_dietToSympt", "output/E3_path_C_diet_to_symptoms.csv", "Secondary - diet chain (main method)",
     "Path C: dietary exposure -> symptoms. Same method and FDR family as E1."),
    ("E1s_linear_sens", "output/E1s_linear_sensitivity.csv", "Sensitivity",
     "Paths A, B and C as linear models (HC3) on the original scales, same adjustment. Beta per 1 FFQ point; per 1 log10 unit of Desulfovibrio. Linear models are sensitive to single influential observations (e.g. in IBS-C two patients with no processed meat and high pain)."),
    ("E4_interaction", "output/E4_interaction_by_subtype.csv", "Secondary - diet chain",
     "Does the association differ between IBS-C and IBS-D? p for the exposure x subtype interaction (adjusted, HC3), on ranks (main) and on original scales (sensitivity)."),
    ("E5_Spearman_crude", "output/E5_Spearman_unadjusted.csv", "Descriptive",
     "Unadjusted Spearman correlations for the chain, by group (for comparison with the adjusted estimates)."),
    ("E6_LinDA_diet_SRB", "output/E6_LinDA_diet_to_sulfidogenic.csv", "Secondary - diet chain",
     "Compositional check of path A: LinDA association of each dietary exposure (entered as rank within the group) with sulfidogenic taxa; log2 fold-change per rank unit."),
    ("E7_mediation", "output/E7_mediation.csv", "Secondary - diet chain (main method)",
     "Causal mediation (Imai et al., R 'mediation'; nonparametric bootstrap, 5000 resamples) on percentile ranks: exposure, Desulfovibrio and outcome as rank/n within the group; contrast = 75th vs 25th percentile of the exposure; effects in percentile points of the outcome (0.10 = 10 percentile points). ACME = indirect effect (through Desulfovibrio), ADE = direct effect. Sensitivity rho = correlation of unmeasured mediator-outcome confounding at which ACME = 0. Proportion mediated is only interpretable when the total effect is clearly non-zero."),
    ("E7s_mediation_linear", "output/E7s_mediation_linear.csv", "Sensitivity",
     "Same mediation models on the original scales (effects in IBS-SSS points)."),
    ("E8_E_values", "output/E8_E_values.csv", "Secondary - robustness",
     "E-values (VanderWeele & Ding 2017) for the significant partial Spearman associations (r -> d -> approximate risk ratio): minimum strength of an unmeasured confounder (risk-ratio scale) needed to explain the association away; E_value_CI for the CI limit closest to the null."),
    ("F1_methanogens", "output/F1_methanogens.csv", "Exploratory",
     "Hydrogenotrophic guild: Methanobrevibacter detection by subtype (Firth logistic regression), its relation with stool frequency and with Desulfovibrio."),
    ("F2_anxiety_DA", "output/F2_anxiety_DA.csv", "Exploratory",
     "Anxiety disorder vs taxa (LinDA, adjusted for subtype, age, sex, BMI, PPI)."),
    ("F2b_anxiety_comm", "output/F2b_anxiety_community.csv", "Exploratory",
     "Anxiety disorder vs Shannon (Beta = adjusted difference) and composition (Beta column = PERMANOVA R2)."),
    ("F3_PPI_context", "output/F3_PPI_context.csv", "Context",
     "PPI use and upper-GI comorbidity by subtype (gastritis, reflux, hiatal hernia, Barrett, dyspepsia, duodenitis, ulcer, H. pylori from free-text comorbidities)."),
    ("F4_power", "output/F4_sensitivity_power.csv", "Context",
     "Sensitivity power: smallest effect detectable with 80% power at alpha 0.05 given the sample sizes."),
]


def read_csv(path, table=None):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    head, body = rows[0], rows[1:]
    if table:
        ti = head.index("Table")
        keep = [i for i, h in enumerate(head) if h not in ("Table", "p_num")]
        body = [[r[i] for i in keep] for r in body if r[ti] == table]
        head = [head[i] for i in keep]
    return head, body


def conv(v):
    if v in ("", "NA", "NaN"):
        return None
    if v in ("TRUE", "FALSE"):
        return v == "TRUE"
    try:
        f = float(v)
        return int(f) if f.is_integer() and "." not in v and "e" not in v.lower() else f
    except ValueError:
        return v


wb = Workbook()
ws = wb.active
ws.title = "README"
info = [
    ("IBS-C vs IBS-D: lifestyle, diet and gut microbiota - complete results", True),
    (f"Database: {db_name}", False),
    (f"Generated: {date.today().isoformat()} from the R pipeline (R 4.5.3; package versions in output/sessionInfo.txt)", False),
    ("Sample: 64 patients analysed (IBS-C n = 35, IBS-D n = 29) of 68 in the database (4 exclusions)", False),
    ("Values are outputs of the R analyses (not spreadsheet formulas); every number can be regenerated with run_all.R.", False),
    ("", False),
    ("Conventions", True),
    ("Adjusted models: age, sex, BMI and PPI use in the previous 12 months; linear models with HC3 robust standard errors.", False),
    ("q = Benjamini-Hochberg FDR within the analysis family named in each sheet description.", False),
    ("Marks: p < 0.05 in bold; q < 0.05 green cell; 0.05 <= q < 0.10 yellow cell.", False),
    ("Microbiota = relative abundances (%) from the laboratory's 16S-based targeted panel; diversity computed from the panel.", False),
    ("Plant food score = fruits + vegetables + legumes + whole grains; western food score = red meat + processed meat + processed sweets.", False),
    ("Dietary frequency scores are ordinal: in every model they enter as ranks (partial Spearman / rank-based); linear models on the original scale are given as sensitivity analyses.", False),
    ("Healthy lifestyle score (0-5): non-smoking, PSQI <= 5, physical activity >= 150 min/week, diet quality (plant - western) above the median, alcohol frequency score <= 1.5.", False),
    ("", False),
    ("Sheets", True),
]
for i, (txt, bold) in enumerate(info, start=1):
    c = ws.cell(row=i, column=1, value=txt)
    c.font = Font(name=FONT, bold=bold, size=12 if i == 1 else 10)
r0 = len(info) + 1
for j, h in enumerate(["Sheet", "Analysis family", "Content and methods"], start=1):
    c = ws.cell(row=r0, column=j, value=h); c.font = Font(name=FONT, bold=True, size=10); c.fill = HEAD
for k, (name, _, fam, desc) in enumerate(SHEETS, start=r0 + 1):
    ws.cell(row=k, column=1, value=name).hyperlink = f"#'{name}'!A1"
    ws.cell(row=k, column=1).font = Font(name=FONT, size=10, color="0563C1", underline="single")
    ws.cell(row=k, column=2, value=fam).font = Font(name=FONT, size=10)
    c = ws.cell(row=k, column=3, value=desc); c.font = Font(name=FONT, size=10); c.alignment = Alignment(wrap_text=True, vertical="top")
ws.column_dimensions["A"].width = 24; ws.column_dimensions["B"].width = 30; ws.column_dimensions["C"].width = 120

for name, path, fam, desc in SHEETS:
    table = {"T1_clinical": "T1", "T2_lifestyle": "T2"}.get(name)
    head, body = read_csv(path, table)
    sh = wb.create_sheet(name)
    sh.cell(row=1, column=1, value=desc).font = Font(name=FONT, italic=True, size=9)
    sh.cell(row=1, column=1).alignment = Alignment(wrap_text=False)
    for j, h in enumerate(head, start=1):
        c = sh.cell(row=2, column=j, value=h)
        c.font = Font(name=FONT, bold=True, size=10); c.fill = HEAD
        c.alignment = Alignment(wrap_text=True, vertical="center"); c.border = Border(bottom=thin)
    qcols = [j for j, h in enumerate(head) if h.lower().startswith("q")]
    pcols = [j for j, h in enumerate(head) if h.lower() == "p" or h.lower().startswith("p_") or h.endswith("_p")]
    for i, row in enumerate(body, start=3):
        vals = [v for v in row] if table else [conv(v) for v in row]   # Table 1/2 values stay text, as printed
        for j, v in enumerate(vals, start=1):
            c = sh.cell(row=i, column=j, value=v); c.font = Font(name=FONT, size=10)
            if isinstance(v, float):
                h = head[j - 1].lower()
                c.number_format = "0.0000" if (h == "p" or h.startswith(("p_", "q")) or h.endswith("_p")) else "0.000"
        # cell-level marks: p < 0.05 bold; q < 0.05 green, 0.05 <= q < 0.10 yellow
        if table is None:
            for j in pcols:
                if isinstance(vals[j], float) and vals[j] < 0.05:
                    sh.cell(row=i, column=j + 1).font = Font(name=FONT, size=10, bold=True)
            for j in qcols:
                if isinstance(vals[j], float) and vals[j] < 0.10:
                    sh.cell(row=i, column=j + 1).fill = SIG if vals[j] < 0.05 else TREND
    for j, h in enumerate(head, start=1):
        width = max([len(str(h))] + [len(str(r[j - 1])) for r in body]) + 2
        sh.column_dimensions[get_column_letter(j)].width = min(max(width, 9), 48)
    sh.freeze_panes = "B3"
    sh.row_dimensions[2].height = 30
wb.save(out)
print("saved", out, "sheets:", len(wb.sheetnames))
