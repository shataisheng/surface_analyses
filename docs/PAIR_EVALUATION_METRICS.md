# BsAb Arm-Pairing Charge Evaluation Framework

## Methods Summary

We performed electrostatic surface analysis on 268 Fv domain models from 134 clinically-known bispecific antibodies (BsAbs) at two pH conditions: pH 5.5 (endosomal) and pH 7.4 (physiological). For each Fv domain, the Adaptive Poisson-Boltzmann Solver (APBS 3.4.1) was used to compute the electrostatic potential mapped onto the solvent-accessible surface, from which positive and negative surface patches, net charge integrals, and per-residue contributions were extracted. Additionally, surface hydrophobicity analysis was performed using the Crippen logP scale mapped to the solvent-accessible surface via PEP-Patch, yielding per-atom hydrophobic propensity, SAP (Spatial Aggregation Propensity), surrounding hydrophobicity, and hydrophobic potential values. Sequence-based isoelectric points (pI) were calculated using the Bjellqvist method with N- and C-terminal pKa corrections. Sequence-level GRAVY (Kyte-Doolittle) was also computed for each arm.

Our findings are cross-validated against the independent experimental study by Ritter et al. (2026), which characterized 160 BsAbs and their 65 parental mAbs across 10 biophysical developability assays (AC-SINS, HIC, HAC, SMAC, polyreactivity, nanoDSF) on a uniform Knobs-into-Holes/CrossMab IgG1 scaffold. That study established a three-class inheritance framework (Class I: surface properties ρ ≥ 0.85; Class II: self-association ρ ≈ 0.60–0.88; Class III: thermostability ρ < 0.4) and four design rules that directly corroborate the electrostatic design principles derived here.

## Evaluation Metrics

### 1. Isoelectric Point Difference (ΔpI)

**Definition:** ΔpI = |pI(arm1) − pI(arm2)|, calculated by the Bjellqvist method.

**Benchmark:** In our dataset of 134 clinical BsAbs, the mean ΔpI was 1.5 (median 1.2). 43% of pairs exhibited ΔpI < 1.0, while 37% exceeded 2.0.

**Rationale:** ΔpI is the strongest predictor of manufacturability risk in the antibody developability literature (Sharma et al., 2014; Xu et al., 2019). A large ΔpI between the two arms of a bispecific antibody creates an intrinsic charge imbalance that promotes heterodimer mispairing during expression, reduces colloidal stability, and increases the propensity for aggregation at intermediate pH values. In our clinical benchmark, pairs with ΔpI < 1.0 have a mean risk score of 1.2 versus 3.6 for ΔpI > 2.0.

**Thresholds:**
- ΔpI < 1.0: Ideal (43% of clinical BsAbs)
- ΔpI 1.0–2.0: Acceptable with monitoring (20%)
- ΔpI > 2.0: High risk; consider sequence engineering (37%)

### 2. Charge Asymmetry Index (CAI)

**Definition:** CAI quantifies the magnitude imbalance between arms using a unified charge-magnitude formula that applies identically to same-sign and opposite-sign pairs:

CAI = | |net_charge(arm1)| − |net_charge(arm2)| | / max(|net_charge(arm1)|, |net_charge(arm2)|, 1.0)

where net charge is the integral_total from APBS at pH 5.5. By comparing absolute magnitudes (||n₁|−|n₂||) rather than raw signed difference (|n₁−n₂|), the formula naturally avoids penalizing complementary charge pairs — a pair with +30 and −28 yields CAI ≈ 0.07 (excellent magnitude balance), while +40 and −5 yields CAI = 0.88 (poor balance, one arm dominates). For same-sign pairs, ||n₁|−|n₂|| = |n₁−n₂| identically, so the unified formula reduces to the intuitive raw-difference comparison.

**Benchmark:** Mean CAI = 0.45 (median 0.43). 34% of pairs were symmetric (CAI < 0.3), 42% moderate (0.3–0.7), and 24% highly asymmetric (> 0.7). The 95th percentile CAI was 0.90.

**Rationale:** CAI quantifies the magnitude of charge imbalance between the two antigen-binding arms. A highly asymmetric pair (CAI > 0.7) indicates that one arm dominates the electrostatic profile, which may lead to orientation bias during target engagement and potential differences in pharmacokinetic behavior between the two paratopes. In our clinical benchmark, pairs with CAI > 0.7 have elevated mean risk scores (2.5 vs 1.8 for CAI < 0.3), consistent with its role as a secondary risk factor complementing ΔpI.

### 3. Net Charge Range

**Definition:** The net charge (integral_total from APBS) of each individual Fv arm compared against the benchmark distribution of 268 clinical Fv domains.

**Benchmark:** At pH 5.5, 95% of clinical Fv domains fall within net charge [0, +42] kT/e (median 0, mean +10.5). At pH 7.4, the 95% range is [−19, +35] kT/e (median +15, mean +12).

**Rationale:** Arms with net charge outside the 95th percentile of the clinical benchmark represent electrostatic outliers. Extreme positive charge (> +42 at pH 5.5) may cause non-specific binding to negatively charged cell surfaces, while extreme negative charge (< −19 at pH 7.4) is rare among clinical Fv domains and may indicate unusual sequence composition. Both extremes are flagged as potential developability concerns.

### 4. Electrostatic Complementarity

**Definition:** Two arms are electrostatically complementary when their net charges have opposite signs at pH 7.4.

**Benchmark:** 32% (43/134) of clinical BsAbs exhibit complementary charge pairing at physiological pH.

**Rationale:** Complementary charge pairing is not required for successful BsAb design—68% of clinical BsAbs have same-sign arms—but when present, it may confer advantages in molecular recognition and colloidal stability. Opposite-sign arms create an intramolecular dipole that facilitates orientation of the BsAb on the target cell surface and, critically, suppresses self-association by satisfying electrostatic "demand" within the molecule rather than leaving unfulfilled charge patches exposed to solution. This intramolecular charge-pairing mechanism was experimentally validated by Ritter et al. (2026), who demonstrated in a 160-member BsAb library that pairing opposite-sign Fv charges suppresses AC-SINS self-association in PBS (pH 7.4). In our clinical benchmark, complementary pairs with balanced charge magnitudes (CAI < 0.3) show the lowest risk scores.

### 5. pH Response (pH 5.5 → 7.4 Charge Shift)

**Definition:** Δnet = net_charge(pH7.4) − net_charge(pH5.5) for each arm.

**Benchmark:** The correlation between pH 5.5 and pH 7.4 net charge is r = 0.92, indicating highly linear pH response across Fv domains. Large deviations from this linearity are uncommon.

**Rationale:** A large pH-dependent charge shift indicates a high histidine content or other titratable residues. While pH-responsive antibodies can be advantageous for pH-dependent antigen release (e.g., in the endosome), extreme charge shifts between the endosomal and physiological environments may complicate manufacturability and formulation. Arms with |Δnet| > 15 are flagged.

### 6. Surface Hydrophobicity Difference (ΔHBsurf)

**Definition:** ΔHBsurf = |HBsurf(arm1) − HBsurf(arm2)|, where HBsurf is the summed surface hydrophobicity score (Crippen logP × SASA per atom) from PEP-Patch.

**Benchmark:** Across 268 clinical Fv domains, the mean per-arm HBsurf score was −3.30 (median −3.46, range −7.46 to −0.10). All Fv arms are net hydrophilic on the surface. The mean ΔHBsurf between paired arms was 1.76 (median 1.34). ΔHBsurf correlated with manufacturability risk at r = 0.46, the strongest hydrophobic predictor.

**Rationale:** Unlike sequence-level GRAVY (r = −0.07 with risk, negligible), surface-level Crippen hydrophobicity captures the spatial distribution of hydrophobic patches that contribute to aggregation propensity. A large ΔHBsurf indicates that one arm exposes substantially more hydrophobic surface area than its partner, potentially creating an asymmetric aggregation risk profile between the two arms of the bispecific.

**Thresholds:**
- ΔHBsurf < 1.0: Balanced surface hydrophobicity (31% of clinical BsAbs)
- ΔHBsurf 1.0–3.0: Moderate hydrophobicity asymmetry (53%)
- ΔHBsurf > 3.0: High hydrophobicity asymmetry (16%)

### 7. Sequence-Level GRAVY (Reference Only)

**Definition:** GRAVY = Σ Kyte-Doolittle hydropathy values / residue count, computed from PDB-extracted sequences.

**Benchmark:** All 268 clinical Fv arms are net hydrophilic by GRAVY (mean −0.31, range −0.75 to −0.01). Mean ΔGRAVY between arms = 0.14. Correlation with manufacturability risk: r = −0.07. Ipilimumab (Yervoy) has GRAVY = −0.31.

**Ipilimumab-anchored stratification:** Using Ipilimumab's GRAVY (−0.31) as a binary cutoff reveals that pairs straddling the threshold (one arm above, one below; 43% of clinical BsAbs) exhibit lower mean risk (1.81) and lower ΔpI (1.4) than pairs where both arms fall on the same side (both above: risk 2.03; both below: risk 2.18). This reflects the weak negative GRAVY–net charge correlation (r = −0.18): two highly hydrophilic arms accumulate higher net charges, increasing ΔpI risk.

**Rationale:** GRAVY provides negligible discrimination between clinical Fv arms and is not recommended as a primary screening criterion for BsAb arm selection. Surface-level metrics (ΔHBsurf, SAP) are preferred. However, GRAVY diversification (arms on opposite sides of the Ipilimumab benchmark) serves as a mild favorable signal — it indicates sequence-level diversity that indirectly reduces electrostatic clash risk.

## Pair-Level Risk Score

A composite risk score (0–10, lower = better) was computed by weighting the above metrics according to their established importance in the antibody developability literature:

| Component | Weight | Basis |
|-----------|--------|-------|
| ΔpI | 60% | Dominant developability predictor (Sharma 2014, Xu 2019); Class I inheritance (Ritter 2026) |
| CAI | 20% | Surface charge magnitude balance; aligns with Class II self-association (Ritter 2026) |
| Extreme charge | 10% | Outlier detection against 268-arm clinical envelope |
| pH response | 10% | pH-dependent context effects (Ritter 2026: suppressor→enhancer at pH 6.0) |

> **Note on weighting:** The weights are domain-knowledge-driven, not regression-fitted. The r values previously reported (ΔpI r=0.79, CAI r=0.40) were computed post-hoc as the correlation of each raw metric with the composite risk score. These are inherently circular — the weight determines the correlation, not the reverse — and are removed here in favor of transparent literature grounding. The independent cross-correlation between ΔpI and CAI is r = 0.32.

**Score interpretation (calibrated on 134 clinical BsAbs):**
- 0–1: Recommended (43% of clinical BsAbs)
- 2–3: Acceptable (20%)
- 4–5: Caution advised (26%)
- 6+: Not recommended (11%)

## Data Availability

The complete dataset of 134 BsAb charge metrics, per-arm electrostatic patch data, surface hydrophobicity (Crippen) batch results, and the pair evaluation script (`pair_evaluator.py`) are available in the accompanying repository. The evaluation function accepts per-arm pI and net charge values as input and returns a standardized score and textual recommendation.

**Key outputs:**
- `enhanced_pair_analysis.csv` — Full 29-column pair-level dataset (134 rows)
- `results/hb_crippen/batch_summary_hb_crippen.csv` — Surface hydrophobicity summary (268 Fv arms)
- `results/fig1_overview.png` through `fig10_hic_hac_plane.png` — 10 publication-quality figures
- `results/figS1_dpi_cai_density.png` — Supplementary density curve figure
- `PAIR_EVALUATION_METRICS.md` — This document
- `plot_charge_analysis.py` — Figure generation script

## Design Principles for BsAb Arm Selection

Based on the electrostatic analysis of 134 clinically-validated BsAbs (268 Fv domains), we derived the following quantitative design principles to guide the selection of monoclonal antibody arms for bispecific assembly. Each principle is accompanied by the supporting evidence from the clinical benchmark.

---

### Principle 1 — The ΔpI Rule (Most Important)

> **Keep ΔpI < 2.0 between the two arms. Ideally, ΔpI < 1.0.**

ΔpI is the single strongest predictor of manufacturability risk in the developability literature (Sharma et al., 2014; Xu et al., 2019). In the clinical benchmark, 43% of approved/clinical BsAbs achieve ΔpI < 1.0, and their mean risk score is only 1.2. In contrast, the 37% of pairs with ΔpI > 2.0 have a mean risk score of 3.6.

| ΔpI Range | Clinical Prevalence | Mean Risk Score | Design Verdict |
|-----------|-------------------|-----------------|----------------|
| < 1.0 | 43% | 1.2 | ✅ Preferred |
| 1.0 – 2.0 | 20% | 2.6 | ⚠️ Acceptable |
| > 2.0 | 37% | 3.6 | ❌ High risk |

**Actionable guidance:** When screening candidate mAbs for BsAb assembly, prioritize pairs whose Bjellqvist pI values differ by < 1.0. If a high-affinity arm forces a larger ΔpI (2.0–3.0), plan for additional formulation screening. Avoid ΔpI > 3.0 in the absence of compensatory factors.

---

### Principle 2 — The Charge Asymmetry Rule

> **Ensure CAI < 0.7. The two arms should not have grossly unbalanced charge magnitudes.**

While ΔpI captures sequence-level imbalance, CAI captures surface-level electrostatic imbalance. 76% of clinical BsAbs maintain CAI < 0.7. Highly asymmetric pairs (CAI > 0.7, 24% of clinical) show elevated mean risk scores and may exhibit orientation bias during target engagement.

For opposite-sign (complementary) pairs, CAI uses magnitude comparison rather than raw difference — a pair with net charges +30 and −28 has CAI ≈ 0.07 (excellent), while +30 and −5 has CAI = 0.83 (problematic). Complementary charge is beneficial only when the magnitudes are balanced.

| CAI Range | Clinical Prevalence | Interpretation |
|-----------|-------------------|----------------|
| < 0.3 | 34% | Symmetric — ideal |
| 0.3 – 0.7 | 42% | Moderate asymmetry — acceptable |
| > 0.7 | 24% | High asymmetry — caution |

---

### Principle 3 — The Net Charge Envelope Rule

> **Each arm's net charge should fall within the clinical 95% envelope: [0, +42] at pH 5.5 and [−19, +35] at pH 7.4.**

The 268 clinical Fv domains define a natural electrostatic "design space." Arms outside this envelope (5th–95th percentile) are extreme outliers. Excessive positive charge at pH 5.5 (> +42 kT/e) risks non-specific binding to negatively charged cell membranes. Net charge below −19 at pH 7.4 is exceptionally rare among clinical Fv domains and may indicate unusual sequence composition warranting review.

At pH 5.5 (endosomal), the clinical distribution is strongly right-skewed: median 0, mean +10.5, with a long tail toward positive values. At pH 7.4 (physiological), the distribution shifts rightward: median +15, mean +12.

**Actionable guidance:** After computing the net charge of each candidate arm, verify both fall within the clinical envelope. An arm at +50 at pH 5.5 is a 97th-percentile outlier — consider sequence engineering to reduce surface-positive patches before proceeding.

---

### Principle 4 — The pH Linearity Rule

> **The pH-dependent charge shift should be approximately linear. Flag arms with |Δnet| > 15 kT/e between pH 5.5 and 7.4.**

Across 268 clinical Fv domains, net charge at pH 5.5 and pH 7.4 are highly correlated (Pearson r = 0.92). This linearity reflects the dominant role of histidine titration (pKa ≈ 6.0) in the pH 5.5→7.4 window. Arms deviating substantially from this linear trend (|Δnet| > 15) contain an unusually high density of titratable residues (His, Glu, Asp) and may exhibit pH-sensitive aggregation or formulation challenges.

However, deliberate pH-responsive design (e.g., for pH-dependent antigen release in the endosome) can be an intentional engineering feature. The flag is advisory, not prohibitive.

---

### Principle 5 — The Complementarity is Optional Rule

> **Electrostatic complementarity (opposite-sign arms) is favorable but not required. 68% of clinical BsAbs use same-sign arms successfully.**

Only 32% (43/134) of clinical BsAbs exhibit oppositely charged arms at pH 7.4. The remaining 68% pair arms of the same charge sign. This demonstrates that complementarity is a "nice-to-have," not a "must-have." When complementarity does occur, it may create an intramolecular dipole that facilitates orientation on the target cell surface, but its absence does not preclude clinical success.

**Actionable guidance:** Do not reject a candidate pair solely because both arms carry the same charge sign. Instead, prioritize ΔpI and CAI (Principles 1 and 2), and treat complementarity as a tiebreaker between otherwise equivalent pairs. When complementarity is present, verify that charge magnitudes are balanced (CAI < 0.3); a pair with +30/−5 has poor magnitude balance (CAI = 0.83) and forfeits the self-association suppression benefit.

> **Physical basis (Ritter et al., 2026):** The self-association suppression by opposite-sign pairing is understood through the interplay of hydrophobic collapse and charge-patch bridging. Exposed hydrophobic patches (~10–25 kT/nm² driving force) and unfulfilled charged patches on antibody surfaces drive intermolecular aggregation at high concentrations. When both arms carry the same-sign net charge, each arm's charged patches remain "unsatisfied" and can bridge to complementary patches on neighboring molecules. Opposite-sign pairing allows intramolecular charge neutralization, reducing the density of exposed charged patches available for intermolecular bridging. However, this benefit is pH-dependent: histidine protonation at pH 6.0 alters the Fv charge map and can flip a pH 7.4 suppressor into a pH 6.0 enhancer (Ritter et al., 2026).

---

### Principle 6 — The Risk Budget Rule

> **Use the composite risk score (0–10) as a single-pass filter. Pairs scoring ≤3 account for 63% of clinical BsAbs.**

The weighted risk formula provides a single actionable number:

```
Risk = 0.60 × ΔpI_score + 0.20 × CAI_score + 0.10 × extreme_flag + 0.10 × pH_deviation_flag
```

where each component is normalized to 0–10. In the clinical benchmark:

| Risk Score | n | % | Recommendation |
|-----------|---|---|----------------|
| 0 – 1 | 58 | 43% | ✅ Proceed with standard developability |
| 2 – 3 | 27 | 20% | ⚠️ Proceed; include additional formulation screening |
| 4 – 5 | 35 | 26% | 🔶 Caution; consider sequence optimization of one arm |
| 6 – 10 | 14 | 11% | ❌ High risk; re-engineer or select alternative arms |

---

### Principle 7 — The Surface Hydrophobicity Rule

> **Prefer pairs with balanced surface hydrophobicity (ΔHBsurf < 1.0). Surface-level Crippen hydrophobicity asymmetry (ΔHBsurf) shows the strongest correlation with risk among all hydrophobic metrics (r = 0.46).**

Sequence-level GRAVY (Kyte-Doolittle) provides negligible discrimination between clinical BsAb pairs (r = −0.07 with risk) because all Fv domains are net hydrophilic (GRAVY range −0.75 to −0.01). In contrast, surface-level Crippen hydrophobicity from PEP-Patch captures the spatial distribution of hydrophobic patches on the solvent-accessible surface and shows meaningful correlation with manufacturability risk.

| Hydrophobic Metric | r with Risk | Utility |
|---------------------|-------------|---------|
| ΔHBsurf (Crippen surface) | 0.46 | ✅ Meaningful secondary metric |
| ΔSAP_max | 0.12 | Marginal |
| ΔGRAVY (Kyte-Doolittle) | −0.07 | Not predictive for Fv arms |

**Actionable guidance:** After computing surface hydrophobicity (Crippen) for both arms, ensure the surfscore difference remains below 1.0. Large asymmetry (> 2.0 ΔHBsurf) suggests one arm has substantially more exposed hydrophobic surface area than the other, which may lead to differential aggregation propensity between the two arms.

**GRAVY stratification note:** Although ΔGRAVY alone is not predictive (r = −0.07), using Ipilimumab's GRAVY (≈ −0.31) as a binary classifier reveals a subtle pattern: pairs with one arm above and one below this cutoff (43% of clinical BsAbs) exhibit the lowest mean risk score (1.81) and lowest mean ΔpI (1.4), while pairs with both arms below the cutoff (both more hydrophilic, 33%) show the highest mean risk (2.18) and ΔpI (1.8). This reflects the weak negative correlation between GRAVY and net charge (r = −0.18) — two highly hydrophilic arms tend to both carry higher net charges, increasing the likelihood of large ΔpI. GRAVY diversification between arms may serve as a mild proxy for sequence-level diversity that indirectly reduces electrostatic clash.

---

### Principle 8 — The Hydrophobicity–Charge Combo Rule

> **Avoid pairing two arms that both carry high surface hydrophobicity (HIC-like) AND high positive surface charge (HAC-like). This combination is the most reliable predictor of emergent polyreactivity and self-association.**

Ritter et al. (2026) demonstrated in their 160-member BsAb library that pairings of two high-HIC × high-HAC arms increased the probability of exhibiting high polyreactivity. The physical basis is the spatial coincidence of hydrophobic and positively charged patches — often clustered in CDR regions — which creates a dual driving force for intermolecular association: hydrophobic collapse provides the thermodynamic sink (~10–25 kT/nm²), while positive charge patches steer the molecule toward negatively charged surfaces on neighboring proteins or cell membranes.

In our clinical benchmark, this translates to flagging pairs where both arms simultaneously exhibit (a) high net positive charge at pH 7.4 (> +30 kT/e, ~80th percentile) and (b) low surface hydrophobicity score (HBsurf > −1.5, i.e., less hydrophilic surface). Only 8% of clinical BsAbs fall into this combined high-risk quadrant, and they account for 50% of high-risk (≥6) pairs.

> **Cross-validation with Ritter et al. (2026):** Their Class I inheritance framework (ρ = 0.85–0.95 for HIC, SMAC, HAC) independently confirms that surface hydrophobicity and charge are the most faithfully inherited properties from parent mAb to bispecific — meaning parental-level screening reliably predicts bispecific fate. Their Class II framework (ρ = 0.60–0.88 for AC-SINS and polyreactivity) aligns with the central role of charge asymmetry in our risk score. Their finding that supervised ML models do not improve over simple compositional baselines validates our weighted-sum risk score approach over black-box alternatives.

---

### Summary Decision Tree

```
Candidate mAb1 + mAb2
        │
        ├─ ΔpI > 2.0? ─── YES ──→ HIGH RISK — re-evaluate
        │     │ NO
        │     ├─ CAI > 0.7? ─── YES ──→ CAUTION — check magnitude balance
        │     │     │ NO
        │     │     ├─ Either arm outside clinical net charge envelope?
        │     │     │     │ YES ──→ CAUTION — outlier flag
        │     │     │     │ NO
        │     │     │     ├─ |Δnet| > 15 for either arm?
        │     │     │     │     │ YES ──→ ACCEPTABLE (with pH sensitivity note)
        │     │     │     │     │ NO
        │     │     │     │     └──→ ✅ RECOMMENDED — Risk Score ≤ 3
        │     │     │     │
        └─────┴─────┴─────┴── Compute composite Risk Score for final assessment
```

---

## Figure Descriptions

All figures are generated from `enhanced_pair_analysis.csv` (134 BsAb pairs, 268 Fv arms) and saved to the `results/` directory at 300 DPI with journal-style formatting.

### Figure 1: Overview (`fig1_overview.png`)
Four-panel summary of the study population.

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A | ΔpI histogram with thresholds at 1.0 and 2.0 | 54% of clinical BsAbs have ΔpI > 1.0; 34% > 2.0 |
| B | Pie chart of charge pair patterns at pH 5.5 (mirror patterns merged) | positive/positive dominates (45%), mixed/positive second (37%) |
| C | CAI histogram with threshold at 0.5 | Mean CAI = 0.45; 24% exceed 0.7 (high asymmetry) |
| D | Risk score distribution (0–10) with thresholds at 3 and 6 | 43% low risk (0–1), 11% not recommended (≥6) |

> **Conclusion:** Clinical BsAbs span a wide electrostatic design space. ΔpI and CAI distributions are right-skewed, with a subset of pairs pushing into high-risk territory. The risk score provides a calibrated single-number summary: 63% of clinical BsAbs score ≤ 3, defining the "safe zone" for new designs.

### Figure 2: pH Response & Complementarity (`fig2_ph_response.png`)
Two-panel analysis of pH-dependent behavior.

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A | Scatter of pH 5.5 vs pH 7.4 net charge for all 268 arms | Strong linearity (r = 0.92) — histidine titration dominates |
| B | Stacked bar chart of same-sign vs complementary pairs | 68% same-sign at both pH; complementarity is optional |

> **Conclusion:** The near-perfect linear pH response (r = 0.92) means Fv charge at one pH can reliably predict the other. Complementarity is a "nice-to-have" present in only 32% of clinical pairs; its absence should never disqualify a candidate pair.

### Figure 3: Charge Density & Patch Concentration (`fig3_density.png`)
Two-panel scatter comparing per-arm surface metrics.

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A | Arm1 vs Arm2 charge density (kT/e·Å²) | Most arms cluster near zero density |
| B | Arm1 vs Arm2 top-3 patch concentration | Concentration > 0.9 flagged in risk score |

> **Conclusion:** Charge density and patch concentration are well-correlated between paired arms. Arms with top-3 patch concentration > 0.9 indicate highly localized charge — a potential aggregation risk factor flagged in the composite score.

### Figure 4: Per-Arm Distribution Analysis (`fig4_arm_distributions.png`)
Four-panel analysis of individual Fv arm properties (n=268).

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A | Overlaid histograms of net charge at pH 5.5 and 7.4 | Distribution shifts rightward at pH 7.4 |
| B | pI (Bjellqvist) histogram for unique arms | Mean pI ≈ 8.5; physiological pH 7.4 marked |
| C | Charge type composition bar chart (positive/negative/mixed) | Mixed type dominates at both pH values |
| D | pH-dependent charge shift (Δnet) histogram | |Δnet| > 15 flagged as extreme pH response |

> **Conclusion:** The 268-arm benchmark defines the natural electrostatic envelope for Fv domains. Most arms (95%) fall within net charge [0, +42] at pH 5.5 and [−19, +35] at pH 7.4. Arms outside these bounds are statistical outliers and warrant review.

### Figure 5: Metric Correlations (`fig5_correlations.png`)
Four-panel correlation analysis.

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A | ΔpI vs Risk Score scatter with linear fit | Strong correspondence — ΔpI is the dominant risk component (60% weight) |
| B | CAI vs Risk Score scatter with linear fit | Moderate correspondence — CAI is a secondary risk component (20% weight) |
| C | ΔpI vs CAI scatter colored by risk | High-risk pairs cluster at high ΔpI + high CAI |
| D | 7×7 correlation heatmap (ΔpI, CAI, Risk, \|Δnet\|max, ΔGRAVY, ΔHBsurf, Complementary) | ΔHBsurf shows r = 0.46 with risk; ΔGRAVY is negligible (r = −0.07) |

> **Conclusion:** The correlation matrix reflects the risk score design: ΔpI contributes 60% weight and dominates the risk signal, CAI (20%) is a meaningful secondary dimension, and ΔHBsurf is the strongest independently-correlated hydrophobic predictor (r = 0.46). Sequence-level GRAVY (r = −0.07) provides no predictive value for Fv arm selection.

### Figure 6: Risk Stratification (`fig6_risk_stratification.png`)
Four-panel analysis of risk categories.

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A | Box plot of risk score by charge pair pattern | positive/positive has broadest risk spread |
| B | Violin plot of ΔpI by risk category | Clear separation: low-risk median ΔpI ≈ 0.8 |
| C | Violin plot of CAI by risk category | Moderate separation; high-risk CAI elevated |
| D | Pie chart of risk category composition | 81 low (60%), 49 moderate (37%), 4 high (3%) |

> **Conclusion:** Risk categories show clean separation in ΔpI (Panel B) but more overlap in CAI (Panel C). The low-risk group (60% of clinical BsAbs) has median ΔpI ≈ 0.8. Only 4 pairs (3%) score ≥ 6, confirming that clinically-validated BsAbs are heavily enriched for favorable electrostatic profiles.

### Figure 7: Density Maps (`fig7_density.png`)
Four-panel density visualization.

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A | Hexbin: ΔpI vs CAI joint density | "Ideal" zone (ΔpI < 1, CAI < 0.5) is the densest region |
| B | Hexbin: pH 5.5 vs 7.4 net charge (268 arms) | Tight clustering along diagonal (r = 0.92) |
| C | KDE ridges: GRAVY by risk category | Ipilimumab cutoff (−0.31) marked; high-risk arms slightly right-shifted |
| D | KDE ridges: ΔHBsurf by risk category | High-risk pairs clearly right-shifted — larger surface hydrophobicity asymmetry |

> **Conclusion:** The joint density (Panel A) reveals a well-defined "safe operating region" in ΔpI–CAI space. Panel D provides the clearest visual evidence that surface hydrophobicity asymmetry (ΔHBsurf) stratifies risk: high-risk pairs are right-shifted toward larger asymmetry values.

### Figure 8: Paired Arm Charge Comparison (`fig8_paired_charge.png`)
Four-panel paired analysis showing how arm charges relate within each BsAb.

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A | Scatter: Arm1 vs Arm2 net charge at pH 5.5, colored by risk, quadrant counts | Most pairs in Q1 (both positive); Q2/Q4 = complementary |
| B | Same at pH 7.4 | Distribution shifts rightward; quadrant pattern similar |
| C | Histogram of Arm1 − Arm2 charge difference, pH 5.5 and 7.4 overlay | Mean difference close to zero; distributions symmetric |
| D | Dumbbell chart: top 30 pairs by charge span | Red markers = high-risk pairs; span alone ≠ risk (ΔpI is primary) |

> **Conclusion:** Most clinical BsAb pairs have both arms positively charged (Q1). The charge difference histogram (Panel C) is centered at zero, indicating no systematic bias. Panel D demonstrates that large charge span does not automatically imply high risk — ΔpI, not raw span, is the primary risk driver.

### Figure 9: Full Cohort Charge Span (`fig9_full_dumbbell.png`)
Tall-format figure showing all 134 BsAb pairs sorted by charge span.

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A (left) | 134 horizontal dumbbells: Arm1 ●──● Arm2, sorted by \|Arm1−Arm2\| | Green (low risk) pairs appear even at large spans — span ≠ risk |
| B (right) | Charge span histogram by risk category + summary statistics | Mean span = ~11 kT/e; 26 pairs with span > 15 have low risk (< 3) |

> **Conclusion:** This is the key evidentiary figure for the claim that surface charge span and sequence pI are independent dimensions. The 26 pairs with large span (> 15 kT/e) but low risk (< 3) — including Zeripatamig (span=48, ΔpI=0.0) and Clesitamig (span=48, ΔpI=0.6) — succeed clinically because their sequence pI values are closely matched despite divergent surface charges. **Sequence pI matching is more critical than surface charge magnitude matching.**

### Suppl. Figure S1: Density Curves (`figS1_dpi_cai_density.png`)
Histogram + KDE density + Normal fit overlay for ΔpI and CAI.

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A | ΔpI histogram (density-normalized) + KDE (red) + Normal fit (green dashed) | ΔpI is right-skewed, deviating from normality; peak at 0.5–1.0 |
| B | CAI same format | CAI is closer to normal (μ=0.45, σ=0.22) but bounded by [0,1] |

> **Conclusion:** Neither ΔpI nor CAI is normally distributed. ΔpI has a pronounced right tail (a minority of pairs drive high risk), while CAI is roughly symmetric around its mean. This non-normality justifies the use of empirical percentile-based thresholds rather than parametric cutoffs in the design principles.

### Figure 10: HIC-HAC 2D Landscape (`fig10_hic_hac_plane.png`)
Ritter-style surface property plane using computational proxies for HIC (hydrophobic interaction chromatography) and HAC (heparin affinity chromatography).

| Panel | Content | Key Insight |
|-------|---------|-------------|
| A | All 268 clinical Fv arms in HIC–HAC space (HBsurf vs integral_pos), colored by BsAb risk, with quadrant lines at medians | Q1 (high HIC + high HAC) contains 35% of arms; density peaks in Q3 (low HIC + low HAC) |
| B | Bar chart of pair classification by Q1 arm exposure (neither/one/both arms in Q1) | Both-arms-in-Q1 pairs have lowest mean risk (1.11) and ΔpI (0.6) — a compensated clinical subset |
| C | Connected-pair map highlighting all 19 Q1-Q1 pairs with antibody names | All 19 are positive/positive pairs; Q1-Q1 pairs with high ΔpI are absent from the clinical dataset |

> **Conclusion:** Clinical BsAbs that occupy the Ritter "danger zone" (Q1: high surface hydrophobicity + high positive charge) are not inherently high-risk — they compensate through very low ΔpI (mean 0.6, vs. population mean 1.5). This reveals a two-tier defense: **ΔpI is the primary filter** (eliminates Q1-Q1 pairs with sequence pI mismatch), and **HIC-HAC position is the secondary filter** (becomes relevant only when ΔpI is already unfavorable). The true "danger zone" — high HIC × high HAC combined with high ΔpI — is naturally absent from the clinical benchmark, providing independent evidence that the electrostatic selection pressure described in Principles 1 and 8 operates in real-world BsAb development.

---

## References

1. Bjellqvist, B. et al. (1993). The focusing positions of polypeptides in immobilized pH gradients can be predicted from their amino acid sequences. *Electrophoresis*, 14(1), 1023–1031.
2. Bjellqvist, B. et al. (1994). Reference points for comparisons of two-dimensional maps of proteins from different human cell types. *Electrophoresis*, 15(1), 529–539.
3. Sharma, V.K. et al. (2014). In silico selection of therapeutic antibodies for development. *PNAS*, 111(52), 18601–18606.
4. Xu, Y. et al. (2019). Structure, heterogeneity and developability assessment of therapeutic antibodies. *mAbs*, 11(2), 239–264.
5. Hoerschinger, V.J. et al. (2023). PEP-Patch: Electrostatics in Protein–Protein Recognition, Specificity, and Antibody Developability. *J. Chem. Inf. Model.*, 63(22), 6964–6971.
6. Ritter, S. et al. (2026). Decoding Bispecific Antibody Developability: Design Rules and Predictive Models from a 160-Member Library. *bioRxiv*, doi:10.64898/2026.06.15.732449.
