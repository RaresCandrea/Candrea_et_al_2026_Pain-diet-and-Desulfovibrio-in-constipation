"""Independent verification of Table 1 / Table 2.

Recomputes every cell directly from the raw Excel workbook (no shared code with the R
pipeline) and compares it with the text in the generated Word tables.
Usage: python 10_verify_tables_independent.py <database.xlsx> <tables_dir>
"""
import sys
from pathlib import Path

import docx
import numpy as np
import openpyxl
from scipy.stats import fisher_exact, mannwhitneyu

xlsx, tdir = sys.argv[1], Path(sys.argv[2])
ws = openpyxl.load_workbook(xlsx, data_only=True).active
r1 = [c.value for c in ws[1]]; r2 = [c.value for c in ws[2]]
ci = r1.index("IBS_subtype")
names = [(r1[i] if i >= ci else r2[i]) for i in range(len(r1))]; names[0] = "Subject"
names = [str(n).rstrip("*").strip() if n else n for n in names]
rows = [dict(zip(names, r)) for r in ws.iter_rows(min_row=3, values_only=True) if r[0]]
rows = [r for r in rows if not r["Rasons of exclusion"]]
C = [r for r in rows if r["IBS_subtype"] == "C"]; D = [r for r in rows if r["IBS_subtype"] == "D"]
groups = [rows, C, D]


def num(r, k):
    return float(r[k])


def yes(r, k):
    return 1 if str(r[k]).strip().upper() == "YES" else 0


# ---- derived (as defined in the table footnotes) ----
for r in rows:
    r["_sss"] = sum(num(r, k) for k in ["IBSSSS_Q1_pain_intensity", "IBSSSS_Q2_pain_days", "IBSSSS_Q3_distension",
                                         "IBSSSS_Q4_bowel_dissatisfaction", "IBSSSS_Q5_interference"])
    r["_psqi"] = sum(num(r, k) for k in ["PSQI_C1_subjective_quality", "PSQI_C2_latency", "PSQI_C3_duration",
                                          "PSQI_C4_efficiency", "PSQI_C5_disturbances", "PSQI_C6_sleep_medication",
                                          "PSQI_C7_daytime_dysfunction"])
    r["_plant"] = sum(num(r, k) for k in ["FFQ_fruits_score", "FFQ_vegetables_score", "FFQ_legumes_score", "FFQ_whole_grains_score"])
    r["_west"] = sum(num(r, k) for k in ["FFQ_red_meat_score", "FFQ_processed_meat_score", "FFQ_processed_sweets_score"])
    r["_pa150"] = 0 if str(r["Physical activity level (mins/week)"]).replace(" ", "").startswith("<") else 1
    r["_abx"] = 1 if "antibiotic" in str(r["Other_medications_12mo_prior"]).lower() else 0
dq = np.array([r["_plant"] - r["_west"] for r in rows]); med_dq = np.median(dq)
for r in rows:
    r["_hls"] = (yes(r, "Smoking") == 0) + (r["_psqi"] <= 5) + r["_pa150"] + ((r["_plant"] - r["_west"]) > med_dq) + (num(r, "FFQ_alcohol_score") <= 1.5)


def fmt_p(p):
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def mq(vals, dec):
    q = np.percentile(vals, [50, 25, 75])          # linear interpolation (= R type 7)
    return "{0:.{d}f} ({1:.{d}f}, {2:.{d}f})".format(*q, d=dec)


def cont(getter, dec):
    v = [np.array([getter(r) for r in g]) for g in groups]
    p = mannwhitneyu(v[2], v[1], use_continuity=True, method="asymptotic").pvalue
    return [mq(x, dec) for x in v] + [fmt_p(p)]


def binary(getter):
    k = [sum(getter(r) for r in g) for g in groups]; n = [len(g) for g in groups]
    p = fisher_exact([[k[1], n[1] - k[1]], [k[2], n[2] - k[2]]]).pvalue
    return [f"{k[i]} ({100 * k[i] / n[i]:.1f}%)" for i in range(3)] + [fmt_p(p)]


def level(getter, value):
    return [f"{sum(getter(r) == value for r in g)} ({100 * sum(getter(r) == value for r in g) / len(g):.1f}%)" for g in groups]


exp = {}  # (table, label) -> [All, C, D, p]  (p None = not checked: Monte Carlo)
T = "T1"
exp[T, "Female sex, n (%)"] = binary(lambda r: r["Sex"] == "F")
exp[T, "Age, years"] = cont(lambda r: num(r, "Age (years)"), 0)
exp[T, "BMI, kg/m²"] = cont(lambda r: num(r, "BMI (kg/m2)"), 1)
bmi = lambda r: 0 if num(r, "BMI (kg/m2)") < 18.5 else 1 if num(r, "BMI (kg/m2)") < 25 else 2 if num(r, "BMI (kg/m2)") < 30 else 3
for i, lab in enumerate(["Underweight (<18.5)", "Normal weight (18.5–24.9)", "Overweight (25.0–29.9)", "Obesity (≥30.0)"]):
    exp[T, lab] = level(bmi, i) + [""]
for t in [1, 2, 6, 7]:
    exp[T, f"Type {t}"] = level(lambda r: int(num(r, "Bristol Stool Scale")), t) + [""]
exp[T, "Bowel movements per week"] = cont(lambda r: num(r, "Weekly Bowel Movement Frequency"), 0)
exp[T, "IBS-SSS total score"] = cont(lambda r: r["_sss"], 0)
for lab, k in [("Abdominal pain intensity", "IBSSSS_Q1_pain_intensity"), ("Days with pain (x10)", "IBSSSS_Q2_pain_days"),
               ("Abdominal distension", "IBSSSS_Q3_distension"), ("Dissatisfaction with bowel habit", "IBSSSS_Q4_bowel_dissatisfaction"),
               ("Interference with daily life", "IBSSSS_Q5_interference")]:
    exp[T, lab] = cont(lambda r, k=k: num(r, k), 0)
sev = lambda r: "Mild" if r["_sss"] < 175 else "Moderate" if r["_sss"] <= 300 else "Severe"
for lab in ["Mild", "Moderate", "Severe"]:
    exp[T, lab] = level(sev, lab) + [""]
for lab, k in [("Anxiety disorder", "Anxiety_dx"), ("Depression", "Depression_dx"), ("Autoimmune disease", "Autoimmune_disease_dx"),
               ("Other comorbidity", "Other_comorbidity_dx"), ("Proton pump inhibitors", "PPI_12mo_prior"),
               ("Laxatives", "Laxatives_12mo_prior"), ("Antidiarrheals", "Antidiarrheals_12mo_prior"),
               ("Antidepressants", "Antidepressants_12mo_prior")]:
    exp[T, lab] = binary(lambda r, k=k: yes(r, k))
import re as _re
exp[T, "Upper gastrointestinal comorbidity"] = binary(lambda r: int(bool(_re.search(r"gastritis|reflux|GERD|hiatal|Barrett|dyspepsia|duoden|ulcer|Helicobacter", str(r.get("Comorbidity_details") or ""), _re.I))))
exp[T, "Antibiotics"] = binary(lambda r: r["_abx"])
for lab, k, dec in [("Fasting glucose, mg/dL", "Serum glucose (mg/dL)", 0), ("Total cholesterol, mg/dL", "Total cholesterol (mg/dL)", 0),
                    ("HDL cholesterol, mg/dL", "HDL cholesterol (mg/dL)", 0), ("LDL cholesterol, mg/dL", "LDL cholesterol (mg/dL)", 0),
                    ("Triglycerides, mg/dL", "Triglycerides (mg/dL)", 0), ("25(OH) vitamin D, ng/mL", "Vitamin D (ng/mL)", 1)]:
    exp[T, lab] = cont(lambda r, k=k: num(r, k), dec)
vd = lambda r: 0 if num(r, "Vitamin D (ng/mL)") < 20 else 1 if num(r, "Vitamin D (ng/mL)") < 30 else 2
for i, lab in enumerate(["Deficiency (<20)", "Insufficiency (20–29.9)", "Sufficiency (≥30)"]):
    exp[T, lab] = level(vd, i) + [""]

T = "T2"
exp[T, "Current smoking, n (%)"] = binary(lambda r: yes(r, "Smoking"))
pa = lambda r: str(r["Physical activity level (mins/week)"]).replace(" ", "")
for lab in ["<150", "150–299", "≥300"]:
    exp[T, lab] = level(pa, lab) + [""]
exp[T, "PSQI global score"] = cont(lambda r: r["_psqi"], 0)
exp[T, "Poor sleep quality (PSQI > 5), n (%)"] = binary(lambda r: int(r["_psqi"] > 5))
exp[T, "Alcohol intake (frequency score)"] = cont(lambda r: num(r, "FFQ_alcohol_score"), 1)
exp[T, "Coffee consumption, n (%)"] = binary(lambda r: yes(r, "FFQ_coffee_YES_NO"))
exp[T, "Water intake, L/day"] = cont(lambda r: num(r, "FFQ_water_L_day"), 1)
exp[T, "Healthy lifestyle score (0–5)"] = cont(lambda r: r["_hls"], 0)
for lab, k in [("Fruits", "fruits"), ("Vegetables", "vegetables"), ("Legumes", "legumes"), ("Whole grains", "whole_grains"),
               ("Red meat", "red_meat"), ("Processed meat", "processed_meat"), ("White meat", "white_meat"), ("Fish", "fish"),
               ("Dairy products", "dairy"), ("Eggs", "eggs"), ("Processed sweets", "processed_sweets")]:
    exp[T, lab] = cont(lambda r, k=k: num(r, f"FFQ_{k}_score"), 1)
exp[T, "Plant food score (0–28)"] = cont(lambda r: r["_plant"], 1)
exp[T, "Western food score (0–21)"] = cont(lambda r: r["_west"], 1)

# ---- Compare with Word tables ----
mism, checked, mid = [], 0, []
for tab, fname in [("T1", "Table1_clinical_characteristics.docx"), ("T2", "Table2_lifestyle_diet.docx")]:
    t = docx.Document(tdir / fname).tables[0]
    seen = set()
    for row in t.rows[1:]:
        cells = row.cells
        if len({c._tc for c in cells}) < 5:
            continue
        lab = cells[0].text.strip(); got = [c.text.strip() for c in cells[1:5]]
        if (tab, lab) not in exp:
            if got[3] not in ("", "—"):
                # category header rows carry a Monte Carlo p value: checked separately
                pass
            continue
        seen.add((tab, lab))
        want = exp[tab, lab]
        for j in range(4):
            if want[j] == "" and got[j] == "":
                continue
            checked += 1
            if want[j] != got[j]:
                mism.append((tab, lab, ["All", "IBS-C", "IBS-D", "p"][j], got[j], want[j]))
    missing = [k for k in exp if k[0] == tab and k not in seen]
    if missing:
        print("Rows expected but not found in Word table:", missing)

print(f"Cells checked: {checked}")
print(f"Mismatches: {len(mism)}")
for m in mism:
    print("  ", " | ".join(m))
