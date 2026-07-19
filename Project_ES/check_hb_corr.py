"""Check hydrophobic metric correlations with risk."""
import csv
import numpy as np

rows = list(csv.DictReader(open("enhanced_pair_analysis.csv", encoding="utf-8-sig")))
risks = np.array([float(r["risk_score"]) for r in rows])

for col in ["delta_gravy", "gravy_asymmetry", "delta_hb_surfscore",
             "delta_hb_sap", "delta_hb_potential"]:
    vals = np.array([float(r[col]) for r in rows if col in r])
    if len(vals) > 1:
        r_val = np.corrcoef(vals, risks)[0, 1]
        print(f"{col:30s}: r={r_val:.3f}")

# Also check per-arm metrics
for col in ["arm1_hb_sap_max", "arm2_hb_sap_max", "arm1_hb_potential_mean",
            "arm2_hb_potential_mean", "arm1_hb_surfscore_sum", "arm2_hb_surfscore_sum"]:
    vals = np.array([float(r[col]) for r in rows if col in r])
    if len(vals) > 1:
        r_val = np.corrcoef(vals, risks)[0, 1]
        print(f"{col:30s}: r={r_val:.3f}")
