#!/usr/bin/env python3
"""
BsAb 双臂电荷深度分析
=====================
多维度评价双抗双臂电荷分离情况，辅助抗体设计决策。

指标:
  1. 电荷不对称指数 (Charge Asymmetry Index, CAI)
  2. 静电互补性 (Electrostatic Complementarity)
  3. pI 差异估算 (ΔpI estimation)
  4. CDR vs Framework 电荷分布
  5. pH 响应敏感性 (pH5.5 → 7.4 电荷变化)
  6. 电荷密度 (Charge Density per Å²)
  7. Patch 集中度 (1-3个最大 patch 占总面积比)
  8. 综合风险评分 (Manufacturability Risk Score)

用法（在本项目目录 Project_ES/ 下运行）:
    python enhanced_pair_analysis.py --ph55 results/pH5_5 --ph74 results/pH7_4
"""

import argparse, csv, os, sys
from collections import defaultdict

import numpy as np

# ─── 精确 pI 计算 (来自 SeqBox / Bjellqvist 方法) ──────────────────
def _extract_sequence_from_pdb(pdb_path: str) -> str:
    """从 PDB 提取氨基酸序列（单字母）。"""
    try:
        from Bio.SeqIO import parse
        for record in parse(pdb_path, "pdb-atom"):
            return str(record.seq)
    except Exception:
        pass
    # Fallback: 简单正则
    import re
    aa3to1 = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
        "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
        "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
    }
    seq = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res = line[17:20].strip()
                if res in aa3to1:
                    aa = aa3to1[res]
                    if not seq or line[21:22] != seq[-1][1] if seq else True:
                        seq.append(aa)
    return "".join(seq)

def calculate_pi_bjellqvist(sequence: str, precision: float = 0.0001) -> float:
    """Bjellqvist pI 计算（二分法）。"""
    from collections import Counter
    sequence = sequence.upper().strip()
    if not sequence:
        return 7.0

    pka = {
        'K': 10.0, 'R': 12.0, 'H': 5.98,
        'D': 4.05, 'E': 4.45, 'C': 9.0, 'Y': 10.0,
        'Nterm': 7.5, 'Cterm': 3.55,
    }
    # N/C-term corrections
    nterm_aa = sequence[0]
    cterm_aa = sequence[-1]
    nterm_corr = {'A':7.59,'M':7.00,'S':6.93,'P':8.36,'T':6.82,'V':7.44,'E':7.70}
    cterm_corr = {'D':4.55,'E':4.75}
    if nterm_aa in nterm_corr:
        pka['Nterm'] = nterm_corr[nterm_aa]
    if cterm_aa in cterm_corr:
        pka['Cterm'] = cterm_corr[cterm_aa]

    cnt = Counter(sequence)

    def net_charge(ph):
        pos = (cnt.get('K',0)/(1+10**(ph-pka['K'])) + cnt.get('R',0)/(1+10**(ph-pka['R']))
               + cnt.get('H',0)/(1+10**(ph-pka['H'])) + 1/(1+10**(ph-pka['Nterm'])))
        neg = (cnt.get('D',0)/(1+10**(pka['D']-ph)) + cnt.get('E',0)/(1+10**(pka['E']-ph))
               + cnt.get('C',0)/(1+10**(pka['C']-ph)) + cnt.get('Y',0)/(1+10**(pka['Y']-ph))
               + 1/(1+10**(pka['Cterm']-ph)))
        return pos - neg

    lo, hi = 0.0, 14.0
    while hi - lo > precision:
        mid = (lo + hi) / 2
        if net_charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


# ─── Kyte-Doolittle GRAVY (疏水性指数) ──────────────────────
KYTE_DOOLITTLE_HYDROPATHY = {
    'A': 1.8,   'R': -4.5,  'N': -3.5,  'D': -3.5,
    'C': 2.5,   'E': -3.5,  'Q': -3.5,  'G': -0.4,
    'H': -3.2,  'I': 4.5,   'L': 3.8,   'K': -3.9,
    'M': 1.9,   'F': 2.8,   'P': -1.6,  'S': -0.8,
    'T': -0.7,  'W': -0.9,  'Y': -1.3,  'V': 4.2,
}


def calculate_gravy(sequence: str) -> float:
    """计算 Kyte-Doolittle GRAVY (Grand Average of Hydropathy)。
    正值=疏水，负值=亲水。"""
    sequence = sequence.upper().strip()
    if not sequence:
        return 0.0
    valid = [aa for aa in sequence if aa in KYTE_DOOLITTLE_HYDROPATHY]
    if not valid:
        return 0.0
    return sum(KYTE_DOOLITTLE_HYDROPATHY[aa] for aa in valid) / len(valid)


def compute_gravy_asymmetry(g1: float, g2: float) -> float:
    """GRAVY 不对称指数，类似 CAI 但用于疏水性。"""
    if g1 * g2 >= 0:
        return round(abs(g1 - g2) / max(abs(g1), abs(g2), 0.01), 2)
    else:
        return round(abs(abs(g1) - abs(g2)) / max(abs(g1), abs(g2), 0.01), 2)


# ─── 数据加载 ───────────────────────────────────────────
def load_summary(path: str) -> dict:
    if not os.path.exists(path):
        return {}
    data = {}
    with open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            data[row["protein"]] = row
    return data


def load_hb_summary(path: str = "results/hb_crippen/batch_summary_hb_crippen.csv") -> dict:
    """Load hydrophobic (Crippen) batch summary."""
    return load_summary(path)


def load_mapping(path: str) -> list[dict]:
    if not os.path.exists(path):
        return []
    return [r for r in csv.DictReader(open(path, encoding="utf-8-sig"))
            if int(r.get("n_arms", 0)) == 2]


def sf(val, default=0.0):
    try:
        return float(val) if val else default
    except (ValueError, TypeError):
        return default


# ─── 核心指标计算 ───────────────────────────────────────
def compute_arm_metrics(arm: dict, all_arms: list[dict]) -> dict:
    """Compute per-arm derived metrics."""
    pos_area = sf(arm.get("total_positive_area_A2"))
    neg_area = sf(arm.get("total_negative_area_A2"))
    total_area = pos_area + neg_area
    int_pos = sf(arm.get("integral_pos"))
    int_neg = sf(arm.get("integral_neg"))
    int_total = sf(arm.get("integral_total"))
    n_pos = int(sf(arm.get("n_positive_patches")))
    n_neg = int(sf(arm.get("n_negative_patches")))

    # Top-3 concentration
    top3_area = sum(sf(arm.get(f"top_pos_{i}_area_A2")) for i in range(1, 4))
    top3_conc = top3_area / max(pos_area, 0.01)

    # Charge density
    charge_density = int_total / max(total_area, 0.01) if total_area > 0 else 0

    # Net charge ratio
    net_ratio = int_pos / max(abs(int_neg), 0.01)

    # Percentile among all arms
    all_pos = [sf(a.get("total_positive_area_A2")) for a in all_arms]
    all_int = [sf(a.get("integral_total")) for a in all_arms]
    pos_pct = sum(1 for v in all_pos if v < pos_area) / max(len(all_pos), 1) * 100
    int_pct = sum(1 for v in all_int if v < int_total) / max(len(all_int), 1) * 100

    return {
        "pos_area": pos_area, "neg_area": neg_area, "total_area": total_area,
        "int_pos": int_pos, "int_neg": int_neg, "int_total": int_total,
        "net_ratio": round(net_ratio, 2),
        "charge_density": round(charge_density, 2),
        "top3_concentration": round(top3_conc, 2),
        "pos_area_pct": round(pos_pct, 1),
        "int_total_pct": round(int_pct, 1),
    }


def classify_charge(m: dict) -> str:
    r = m["net_ratio"]
    if m["total_area"] < 5:
        return "neutral"
    if r > 3:
        return "positive"
    if r < 0.33:
        return "negative"
    return "mixed"


def estimate_pI_from_pdb(pdb_stem: str, pdb_dir: str = ".") -> float:
    """从 PDB 文件用 Bjellqvist 方法精确计算 pI。"""
    pdb_path = os.path.join(pdb_dir, pdb_stem + "_fixed.pdb")
    if os.path.exists(pdb_path):
        seq = _extract_sequence_from_pdb(pdb_path)
        if seq:
            return calculate_pi_bjellqvist(seq)
    return 7.0  # fallback


def risk_score(m1: dict, m2: dict, delta_pI: float, asym: float) -> dict:
    """Manufacturability risk score (0-10, lower is better)."""
    score = 0
    flags = []

    # pI difference > 2 → aggregation risk
    if delta_pI > 2:
        score += 3
        flags.append("high_ΔpI")
    elif delta_pI > 1:
        score += 1
        flags.append("moderate_ΔpI")

    # Extreme charge (top/bottom 10%)
    if m1["pos_area_pct"] > 90 or m1["pos_area_pct"] < 10:
        score += 1
    if m2["pos_area_pct"] > 90 or m2["pos_area_pct"] < 10:
        score += 1

    # High asymmetry
    if asym > 5:
        score += 2
        flags.append("high_asymmetry")
    elif asym > 3:
        score += 1

    # Too few patches (concentrated charge)
    if m1["top3_concentration"] > 0.9 or m2["top3_concentration"] > 0.9:
        score += 1
        flags.append("concentrated_charge")

    # Very high charge density
    if abs(m1["charge_density"]) > 50 or abs(m2["charge_density"]) > 50:
        score += 1

    return {"risk_score": score, "flags": "|".join(flags) if flags else "low_risk"}


# ─── 主分析 ────────────────────────────────────────────
def analyze(ph55_dir: str, ph74_dir: str, mapping_csv: str = "BsAb_mapping.csv",
            output: str = "enhanced_pair_analysis.csv",
            hb_summary: str = "results/hb_crippen/batch_summary_hb_crippen.csv"):
    ph55 = load_summary(os.path.join(ph55_dir, "batch_summary_es.csv"))
    ph74 = load_summary(os.path.join(ph74_dir, "batch_summary_es.csv"))
    hb = load_hb_summary(hb_summary)
    pairs = load_mapping(mapping_csv)
    _has_hb = len(hb) > 0

    # Build full arm list for percentile calculations
    all_arms_55 = list(ph55.values())
    all_arms_74 = list(ph74.values())

    rows = []
    for pair in pairs:
        ab = pair["antibody"]
        t1 = pair["arm1_target"]
        t2 = pair["arm2_target"]
        a1_stem = os.path.splitext(pair["arm1_pdb"])[0].replace("_fixed", "")
        a2_stem = os.path.splitext(pair["arm2_pdb"])[0].replace("_fixed", "")

        a1_55 = ph55.get(a1_stem, {})
        a2_55 = ph55.get(a2_stem, {})
        a1_74 = ph74.get(a1_stem, {})
        a2_74 = ph74.get(a2_stem, {})

        if not a1_55 or not a2_55:
            continue

        m1_55 = compute_arm_metrics(a1_55, all_arms_55)
        m2_55 = compute_arm_metrics(a2_55, all_arms_55)
        m1_74 = compute_arm_metrics(a1_74, all_arms_74) if a1_74 else m1_55
        m2_74 = compute_arm_metrics(a2_74, all_arms_74) if a2_74 else m2_55

        # ── 核心指标 ──
        net1 = m1_55["int_total"]
        net2 = m2_55["int_total"]
        net1_74 = m1_74["int_total"]
        net2_74 = m2_74["int_total"]

        # 1. 电荷不对称指数 CAI — 统一幅度比较公式
        def calc_cai(n1, n2):
            return round(abs(abs(n1) - abs(n2)) / max(abs(n1), abs(n2), 1.0), 2)

        cai_55 = calc_cai(net1, net2)
        cai_74 = calc_cai(net1_74, net2_74)

        # 2. 静电互补性 = -(sign(arm1_net) * sign(arm2_net))
        #   +1 = same sign (repulsive), -1 = opposite (complementary)
        complement_55 = -1 if net1 * net2 < 0 else 1
        complement_74 = -1 if net1_74 * net2_74 < 0 else 1

        # 3. pI (Bjellqvist 精确计算)
        a1_stem_full = pair["arm1_pdb"].replace("_fixed.pdb", "")
        a2_stem_full = pair["arm2_pdb"].replace("_fixed.pdb", "")
        pI1 = estimate_pI_from_pdb(a1_stem_full)
        pI2 = estimate_pI_from_pdb(a2_stem_full)
        delta_pI = round(abs(pI1 - pI2), 1)

        # 4. GRAVY (Kyte-Doolittle 序列疏水性)
        seq1 = _extract_sequence_from_pdb(os.path.join(".", a1_stem_full + "_fixed.pdb"))
        seq2 = _extract_sequence_from_pdb(os.path.join(".", a2_stem_full + "_fixed.pdb"))
        gravy1 = round(calculate_gravy(seq1), 2) if seq1 else 0.0
        gravy2 = round(calculate_gravy(seq2), 2) if seq2 else 0.0
        delta_gravy = round(abs(gravy1 - gravy2), 2)
        gravy_asym = compute_gravy_asymmetry(gravy1, gravy2)

        # 5. 表面疏水性 (Crippen) - from PEP-Patch hydrophobic batch
        hb1 = hb.get(a1_stem, {})
        hb2 = hb.get(a2_stem, {})
        if hb1 and hb2:
            hb_surf1 = sf(hb1.get("hb_surfscore_sum"))
            hb_surf2 = sf(hb2.get("hb_surfscore_sum"))
            hb_sap1 = sf(hb1.get("hb_sap_max"))
            hb_sap2 = sf(hb2.get("hb_sap_max"))
            hb_pot1 = sf(hb1.get("hb_potential_mean"))
            hb_pot2 = sf(hb2.get("hb_potential_mean"))
            hb_sh1 = sf(hb1.get("hb_sh_max"))
            hb_sh2 = sf(hb2.get("hb_sh_max"))
            delta_hb_surf = round(abs(hb_surf1 - hb_surf2), 2)
            delta_hb_sap = round(abs(hb_sap1 - hb_sap2), 4)
            delta_hb_pot = round(abs(hb_pot1 - hb_pot2), 4)
        else:
            hb_surf1 = hb_surf2 = hb_sap1 = hb_sap2 = 0.0
            hb_pot1 = hb_pot2 = hb_sh1 = hb_sh2 = 0.0
            delta_hb_surf = delta_hb_sap = delta_hb_pot = 0.0

        # 6. pH 响应
        pH_response_arm1 = round(net1_74 - net1, 2)
        pH_response_arm2 = round(net2_74 - net2, 2)

        # 5. 综合风险 (use Bjellqvist delta_pI)
        risk = risk_score(m1_55, m2_55, delta_pI, cai_55)

        row = {
            "antibody": ab,
            "arm1_target": t1,
            "arm2_target": t2,
            # pH 5.5 metrics
            "pH5.5_arm1_net_charge": round(net1, 1),
            "pH5.5_arm2_net_charge": round(net2, 1),
            "pH5.5_charge_asymmetry_CAI": cai_55,
            "pH5.5_complementarity": complement_55,
            "pH5.5_arm1_charge_type": classify_charge(m1_55),
            "pH5.5_arm2_charge_type": classify_charge(m2_55),
            "pH5.5_pair_pattern": f"{classify_charge(m1_55)}/{classify_charge(m2_55)}",
            "pH5.5_arm1_pI_Bjellqvist": pI1,
            "pH5.5_arm2_pI_Bjellqvist": pI2,
            "pH5.5_delta_pI": delta_pI,
            # pH 7.4 metrics
            "pH7.4_arm1_net_charge": round(net1_74, 1),
            "pH7.4_arm2_net_charge": round(net2_74, 1),
            "pH7.4_charge_asymmetry_CAI": cai_74,
            "pH7.4_complementarity": complement_74,
            "pH7.4_arm1_charge_type": classify_charge(m1_74),
            "pH7.4_arm2_charge_type": classify_charge(m2_74),
            "pH7.4_pair_pattern": f"{classify_charge(m1_74)}/{classify_charge(m2_74)}",
            "pH7.4_arm1_pI_Bjellqvist": pI1,
            "pH7.4_arm2_pI_Bjellqvist": pI2,
            "pH7.4_delta_pI": delta_pI,
            # pH response & density
            "pH_response_arm1": pH_response_arm1,
            "pH_response_arm2": pH_response_arm2,
            "arm1_charge_density": m1_55["charge_density"],
            "arm2_charge_density": m2_55["charge_density"],
            "arm1_top3_conc": m1_55["top3_concentration"],
            "arm2_top3_conc": m2_55["top3_concentration"],
            # GRAVY (sequence-level hydrophobicity)
            "arm1_gravy": gravy1,
            "arm2_gravy": gravy2,
            "delta_gravy": delta_gravy,
            "gravy_asymmetry": gravy_asym,
            # Surface hydrophobicity (Crippen via PEP-Patch)
            "arm1_hb_surfscore_sum": round(hb_surf1, 2),
            "arm2_hb_surfscore_sum": round(hb_surf2, 2),
            "delta_hb_surfscore": delta_hb_surf,
            "arm1_hb_sap_max": round(hb_sap1, 4),
            "arm2_hb_sap_max": round(hb_sap2, 4),
            "delta_hb_sap": delta_hb_sap,
            "arm1_hb_potential_mean": round(hb_pot1, 4),
            "arm2_hb_potential_mean": round(hb_pot2, 4),
            "delta_hb_potential": delta_hb_pot,
            # Risk
            "risk_score": risk["risk_score"],
            "risk_flags": risk["flags"],
        }
        rows.append(row)

    # Output
    if rows:
        cols = list(rows[0].keys())
        with open(output, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print(f"Saved: {output} ({len(rows)} pairs)")

    # ── Summary Statistics ──
    print(f"\n{'='*60}")
    print(f"  BsAb Charge Analysis Summary ({len(rows)} antibodies)")
    print(f"{'='*60}")

    for label, ph, complement_key, asym_key, dpI_key in [
        ("pH 5.5 (endosomal)", "pH5.5", "pH5.5_complementarity", "pH5.5_charge_asymmetry_CAI", "pH5.5_delta_pI"),
        ("pH 7.4 (physiological)", "pH7.4", "pH7.4_complementarity", "pH7.4_charge_asymmetry_CAI", "pH7.4_delta_pI"),
    ]:
        print(f"\n--- {label} ---")

        # Charge type distribution
        types = defaultdict(int)
        for r in rows:
            types[f"{r[f'{ph}_arm1_charge_type']}/{r[f'{ph}_arm2_charge_type']}"] += 1
        print("  Pair patterns:")
        for t, n in sorted(types.items(), key=lambda x: -x[1]):
            print(f"    {t:25s}: {n:3d} ({100*n/len(rows):.1f}%)")

        # Complementarity
        comp = [r[complement_key] for r in rows]
        same = sum(1 for c in comp if c > 0)
        opp = sum(1 for c in comp if c < 0)
        print(f"  Same-sign: {same} ({100*same/len(rows):.1f}%)  Opposite-sign (complementary): {opp} ({100*opp/len(rows):.1f}%)")

        # Asymmetry
        asym = [r[asym_key] for r in rows]
        print(f"  CAI (asymmetry): mean={np.mean(asym):.2f}  median={np.median(asym):.2f}  >0.5: {sum(1 for a in asym if a>0.5)}")

        # Delta pI
        dpI = [r[dpI_key] for r in rows]
        print(f"  ΔpI: mean={np.mean(dpI):.1f}  median={np.median(dpI):.1f}  >2.0 (risk): {sum(1 for d in dpI if d>2.0)}  >1.0: {sum(1 for d in dpI if d>1.0)}")

    # Risk distribution
    risks = [r["risk_score"] for r in rows]
    print(f"\n--- Risk Score ---")
    print(f"  Mean: {np.mean(risks):.1f}/10")
    for threshold, label in [(0, "Low"), (3, "Moderate"), (6, "High")]:
        n = sum(1 for r in risks if r >= threshold)
        print(f"  {label} risk (>= {threshold}): {n} ({100*n/len(rows):.1f}%)")

    # pH response
    resp1 = [r["pH_response_arm1"] for r in rows]
    resp2 = [r["pH_response_arm2"] for r in rows]
    print(f"\n--- pH Response (pH5.5→7.4 Δ net charge) ---")
    print(f"  Arm1: mean={np.mean(resp1):.1f}  Arm2: mean={np.mean(resp2):.1f}")
    print(f"  Large response (|Δ|>5): arm1={sum(1 for r in resp1 if abs(r)>5)}  arm2={sum(1 for r in resp2 if abs(r)>5)}")


def main():
    p = argparse.ArgumentParser(description="Enhanced BsAb pair charge analysis")
    p.add_argument("--ph55", required=True, help="pH 5.5 results dir")
    p.add_argument("--ph74", required=True, help="pH 7.4 results dir")
    p.add_argument("--mapping", default="BsAb_mapping.csv")
    p.add_argument("-o", "--output", default="enhanced_pair_analysis.csv")
    args = p.parse_args()
    analyze(args.ph55, args.ph74, args.mapping, args.output)


if __name__ == "__main__":
    main()
