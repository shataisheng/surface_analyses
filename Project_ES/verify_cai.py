import csv
rows = list(csv.DictReader(open("enhanced_pair_analysis.csv", encoding="utf-8-sig")))
identical = 0
diff = 0
for r in rows:
    n1 = float(r["pH5.5_arm1_net_charge"])
    n2 = float(r["pH5.5_arm2_net_charge"])
    same = (n1 * n2) >= 0
    if same:
        cai_old = abs(n1 - n2) / max(abs(n1), abs(n2), 1.0)
    else:
        cai_old = abs(abs(n1) - abs(n2)) / max(abs(n1), abs(n2), 1.0)
    cai_new = abs(abs(n1) - abs(n2)) / max(abs(n1), abs(n2), 1.0)
    if abs(cai_old - cai_new) < 1e-10:
        identical += 1
    else:
        diff += 1
        if diff <= 3:
            print(f"DIFF: {r['antibody']} old={cai_old:.6f} new={cai_new:.6f} n1={n1} n2={n2}")

print(f"\nIdentical: {identical}, Different: {diff}")
print("Conclusion: unified formula is mathematically equivalent to the two-case formula.")
