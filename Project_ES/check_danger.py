import csv
rows=list(csv.DictReader(open("enhanced_pair_analysis.csv",encoding="utf-8-sig")))
mapping={}
with open("BsAb_mapping.csv",encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if int(r.get("n_arms",0))==2:
            mapping[r["antibody"]]=(r["arm1_pdb"].replace("_fixed.pdb",""),r["arm2_pdb"].replace("_fixed.pdb",""))
es55={}
with open("results/pH5_5/batch_summary_es.csv",encoding="utf-8-sig") as f:
    for r in csv.DictReader(f): es55[r["protein"]]=r

for row in rows:
    ab=row["antibody"]
    if ab not in mapping: continue
    a1,a2=mapping[ab]
    e1=es55.get(a1,{}); e2=es55.get(a2,{})
    hac1=float(e1.get("integral_pos",0)); hac2=float(e2.get("integral_pos",0))
    hic1=float(row.get("arm1_hb_surfscore_sum",0)); hic2=float(row.get("arm2_hb_surfscore_sum",0))
    d1=(hic1>-3.46)and(hac1>28.2); d2=(hic2>-3.46)and(hac2>28.2)
    if d1 and d2:
        print(ab, row["risk_score"], row["pH5.5_delta_pI"], row["pH5.5_pair_pattern"],
              f"hic={hic1:.1f}/{hic2:.1f} hac={hac1:.0f}/{hac2:.0f}")
