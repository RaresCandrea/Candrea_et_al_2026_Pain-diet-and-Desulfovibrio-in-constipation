"""Independent re-computation of the key results (Python: scipy/statsmodels, raw workbook),
compared with the R outputs in output/*.csv.  Usage: python 13_verify_results_independent.py <database.xlsx>
"""
import csv
import sys

import numpy as np
import openpyxl
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

ws = openpyxl.load_workbook(sys.argv[1], data_only=True).active
r1 = [c.value for c in ws[1]]; r2 = [c.value for c in ws[2]]
ci = r1.index("IBS_subtype")
names = [(r1[i] if i >= ci else r2[i]) for i in range(len(r1))]; names[0] = "Subject"
names = [str(n).rstrip("*").strip() if n else f"col{i}" for i, n in enumerate(names)]
df = pd.DataFrame([r for r in ws.iter_rows(min_row=3, values_only=True) if r[0]], columns=names)
df = df[df["Rasons of exclusion"].isna()].reset_index(drop=True)
gs, ge = r2.index("Prevotella spp.**"), r2.index("Saccharomyces cerevisiae**")
taxa = [n for n in names[gs:ge + 1] if n != "Methanobacteria"]
genus = [t for t in taxa if t.endswith("spp.")]
for c in genus + ["Age (years)", "BMI (kg/m2)", "FFQ_processed_meat_score", "FFQ_red_meat_score", "IBSSSS_Q1_pain_intensity",
                  "IBSSSS_Q2_pain_days", "IBSSSS_Q3_distension", "pH", "Butyrate production (%)", "Acetate / propionate"]:
    df[c] = df[c].astype(float)
C, D = df["IBS_subtype"] == "C", df["IBS_subtype"] == "D"

P = pd.DataFrame(df[genus].values, columns=genus)
p = P.div(P.sum(axis=1), axis=0)
df["Shannon"] = -(p * np.log(p.where(p > 0, 1))).sum(axis=1)
ds = df["Desulfovibrio spp."]
df["Dlog"] = np.log10(ds + ds[ds > 0].min() / 2)
m = pd.DataFrame({"y": df["Dlog"], "D": D.astype(int), "Age": df["Age (years)"], "Male": (df["Sex"] == "M").astype(int),
                  "BMI": df["BMI (kg/m2)"], "PPI": (df["PPI_12mo_prior"] == "YES").astype(int),
                  "PM": df["FFQ_processed_meat_score"], "Q1": df["IBSSSS_Q1_pain_intensity"]})


def r_csv(name):
    return pd.read_csv(f"output/{name}.csv")


def ols(formula, data, term):
    f = smf.ols(formula, data).fit(cov_type="HC3", use_t=True)
    return f.params[term], f.pvalues[term]


checks = []
def add(what, python, r, tol):
    ok = abs(python - r) <= tol
    checks.append((what, python, r, "OK" if ok else "MISMATCH"))


A1 = r_csv("A1_alpha_diversity").set_index("Measure")
add("Shannon genus: Mann-Whitney p", mannwhitneyu(df.loc[D, "Shannon"], df.loc[C, "Shannon"], method="asymptotic").pvalue,
    A1.loc["Shannon index, genus level (panel)", "p_MannWhitney"], 1e-6)
B1 = r_csv("B1_DA_subtype").set_index("Taxon")
add("Desulfovibrio: Wilcoxon p", mannwhitneyu(ds[D], ds[C], method="asymptotic").pvalue, B1.loc["Desulfovibrio spp.", "p_Wilcoxon"], 1e-6)
add("Desulfovibrio: median IBS-C", ds[C].median(), B1.loc["Desulfovibrio spp.", "Median_IBS_C"], 1e-9)
add("Desulfovibrio: median IBS-D", ds[D].median(), B1.loc["Desulfovibrio spp.", "Median_IBS_D"], 1e-9)
B2 = r_csv("B2_Desulfovibrio_sensitivity")
add("Desulfovibrio detected: Fisher p", fisher_exact(pd.crosstab(D, ds > 0).values).pvalue,
    B2.loc[B2.Analysis.str.startswith("Detected"), "p"].iloc[0], 1e-6)
add("Desulfovibrio log10 ~ subtype (adjusted, HC3): p", ols("y ~ D + Age + Male + BMI + PPI", m, "D")[1],
    B2.loc[B2.Analysis.str.startswith("Linear model"), "p"].iloc[0], 1e-6)
E1 = r_csv("E1s_linear_sensitivity"); mc = m[C.values]
row = E1[(E1.Group == "IBS-C") & (E1.Predictor == "Processed meat") & (E1.Outcome == "Desulfovibrio")].iloc[0]
b, pv = ols("y ~ PM + Age + Male + BMI + PPI", mc, "PM")
add("IBS-C processed meat -> Desulfovibrio: beta", b, row.Beta, 1e-6); add("  ... p", pv, row.p, 1e-6)
E2 = E1
row = E2[(E2.Group == "IBS-C") & (E2.Predictor == "Desulfovibrio") & (E2.Outcome == "Abdominal pain intensity")].iloc[0]
b, pv = ols("Q1 ~ y + Age + Male + BMI + PPI", mc, "y")
add("IBS-C Desulfovibrio -> pain intensity: beta", b, row.Beta, 1e-6); add("  ... p", pv, row.p, 1e-6)
E3 = E1
row = E3[(E3.Group == "IBS-C") & (E3.Outcome == "Abdominal pain intensity") & (E3.Predictor == "Processed meat")].iloc[0]
b, pv = ols("Q1 ~ PM + Age + Male + BMI + PPI", mc, "PM")
add("IBS-C processed meat -> pain intensity: beta", b, row.Beta, 1e-6); add("  ... p", pv, row.p, 1e-6)
E5 = r_csv("E5_Spearman_unadjusted")
row = E5[(E5.Group == "IBS-C") & (E5.Variable_1 == "Processed meat") & (E5.Variable_2 == "Desulfovibrio")].iloc[0]
add("IBS-C processed meat vs Desulfovibrio: Spearman rho", spearmanr(mc.PM, mc.y).statistic, row.rho, 1e-9)
C1 = r_csv("C1_functional_by_subtype").set_index("Indicator")
for v in ["pH", "Butyrate production (%)", "Acetate / propionate"]:
    add(f"{v}: Mann-Whitney p", mannwhitneyu(df.loc[D, v], df.loc[C, v], method="asymptotic").pvalue, C1.loc[v, "p_MannWhitney"], 1e-6)
F3 = r_csv("F3_PPI_context")
add("PPI by subtype: Fisher p", fisher_exact(pd.crosstab(D, df["PPI_12mo_prior"] == "YES").values).pvalue,
    F3.loc[F3.Analysis == "PPI use by subtype", "p"].iloc[0], 1e-9)

# partial Spearman (rank residuals on age, sex, BMI, PPI), IBS-C
from scipy.stats import rankdata, t as tdist
def pspear_py(x, y, dfc):
    Z = np.column_stack([np.ones(len(dfc)), dfc.Age, dfc.Male, dfc.BMI, dfc.PPI])
    rx = rankdata(x); ry = rankdata(y)
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]; ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    r = np.corrcoef(ex, ey)[0, 1]; df = len(dfc) - 2 - 4
    return r, 2 * tdist.sf(abs(r * np.sqrt(df / (1 - r ** 2))), df)
mc = mc.assign(W=df.loc[C.values, "FFQ_red_meat_score"].values + df.loc[C.values, "FFQ_processed_meat_score"].values)
for lab, xcol, ycol, sheet, pred, out in [
        ("IBS-C processed meat -> pain (partial Spearman)", "PM", "Q1", "E3_path_C_diet_to_symptoms", "Processed meat", "Abdominal pain intensity"),
        ("IBS-C Desulfovibrio -> pain (partial Spearman)", "y", "Q1", "E2_path_B_Desulfovibrio_to_symptoms", "Desulfovibrio", "Abdominal pain intensity"),
        ("IBS-C processed meat -> Desulfovibrio (partial Spearman)", "PM", "y", "E1_path_A_diet_to_Desulfovibrio", "Processed meat", "Desulfovibrio"),
        ("IBS-C red+processed meat -> Desulfovibrio (partial Spearman)", "W", "y", "E1_path_A_diet_to_Desulfovibrio", "Red + processed meat", "Desulfovibrio")]:
    rr, pp = pspear_py(mc[xcol], mc[ycol], mc)
    row = r_csv(sheet).query("Group == 'IBS-C' and Predictor == @pred and Outcome == @out").iloc[0]
    add(lab + ": rho", rr, row.partial_rho, 1e-6); add("  ... p", pp, row.p, 1e-6)

w = max(len(c[0]) for c in checks)
for what, py, r, st in checks:
    print(f"{what:<{w}}  python={py:.6g}  R={r:.6g}  {st}")
bad = [c for c in checks if c[3] != "OK"]
print(f"\n{len(checks)} checks, {len(bad)} mismatches")
