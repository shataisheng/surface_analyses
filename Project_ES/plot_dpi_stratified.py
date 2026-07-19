#!/usr/bin/env python3
"""Figure S2: ΔpI-stratified charge analysis of high-ΔpI BsAbs."""

import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import Counter

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9, "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7,
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.linewidth": 0.8,
})

rows = list(csv.DictReader(open("enhanced_pair_analysis.csv", encoding="utf-8-sig")))

# Groups
groups = [
    ("ΔpI < 1.0\n(n=58)", [r for r in rows if float(r["pH5.5_delta_pI"]) < 1.0], "#27ae60"),
    ("ΔpI 1.0–2.0\n(n=27)", [r for r in rows if 1.0 <= float(r["pH5.5_delta_pI"]) < 2.0], "#f39c12"),
    ("ΔpI 2.0–3.0\n(n=24)", [r for r in rows if 2.0 <= float(r["pH5.5_delta_pI"]) < 3.0], "#e67e22"),
    ("ΔpI ≥ 3.0\n(n=25)", [r for r in rows if float(r["pH5.5_delta_pI"]) >= 3.0], "#e74c3c"),
]

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
((ax1, ax2), (ax3, ax4)) = axes

# ── S2a: Charge span by ΔpI group (boxplot) ──
span_data = []
for label, grp, color in groups:
    spans = [abs(float(r["pH5.5_arm1_net_charge"]) - float(r["pH5.5_arm2_net_charge"])) for r in grp]
    span_data.append(spans)

bp = ax1.boxplot(span_data, patch_artist=True, widths=0.5,
                  medianprops={"color": "black", "linewidth": 1.2},
                  flierprops={"marker": "o", "markersize": 3, "alpha": 0.4})
for patch, (_, _, color) in zip(bp["boxes"], groups):
    patch.set_facecolor(color)
    patch.set_alpha(0.4)

ax1.set_xticklabels([g[0] for g in groups], fontsize=7)
ax1.set_ylabel("Charge Span |Arm1 − Arm2| (kT/e)")
ax1.set_title("A  Charge Span by ΔpI Group")
# Add mean labels
for i, data in enumerate(span_data):
    ax1.text(i + 1, np.max(data) + 1.5, f"μ={np.mean(data):.1f}", ha="center", fontsize=7, color="#555")

# ── S2b: Arm1 vs Arm2 net charge scatter, colored by ΔpI ──
for label, grp, color in reversed(groups):  # reversed so high dPI on top
    n1 = [float(r["pH5.5_arm1_net_charge"]) for r in grp]
    n2 = [float(r["pH5.5_arm2_net_charge"]) for r in grp]
    ax2.scatter(n1, n2, c=color, alpha=0.6, s=14, edgecolors="white",
                linewidth=0.3, label=label, zorder=2)

lim = max(abs(np.concatenate([
    np.array([float(r["pH5.5_arm1_net_charge"]) for r in rows]),
    np.array([float(r["pH5.5_arm2_net_charge"]) for r in rows])
]))) + 5
ax2.plot([-lim, lim], [-lim, lim], "k--", linewidth=0.5, alpha=0.3)
ax2.axhline(0, color="gray", linewidth=0.3)
ax2.axvline(0, color="gray", linewidth=0.3)
ax2.set_xlabel("Arm1 Net Charge (kT/e)")
ax2.set_ylabel("Arm2 Net Charge (kT/e)")
ax2.set_title("B  Arm1 vs Arm2 Net Charge by ΔpI")
ax2.legend(frameon=False, fontsize=6.5, loc="upper left")

# ── S2c: Complementary pair fraction by ΔpI group ──
comp_fracs = []
for label, grp, color in groups:
    n_comp = sum(1 for r in grp if float(r["pH7.4_complementarity"]) < 0)
    comp_fracs.append(100 * n_comp / len(grp))

bars = ax3.bar(range(4), comp_fracs, color=[g[2] for g in groups],
               edgecolor="white", linewidth=0.5, alpha=0.85)
for bar, val in zip(bars, comp_fracs):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
             f"{val:.0f}%", ha="center", fontsize=9, fontweight="bold")
ax3.set_xticks(range(4))
ax3.set_xticklabels([g[0] for g in groups], fontsize=7)
ax3.set_ylabel("Complementary Pairs (%)")
ax3.set_title("C  Complementary Fraction by ΔpI")

# ── S2d: Risk score distribution by ΔpI group ──
risk_data = []
for label, grp, color in groups:
    risks = [float(r["risk_score"]) for r in grp]
    risk_data.append(risks)

bp2 = ax4.boxplot(risk_data, patch_artist=True, widths=0.5,
                   medianprops={"color": "black", "linewidth": 1.2},
                   flierprops={"marker": "o", "markersize": 3, "alpha": 0.4})
for patch, (_, _, color) in zip(bp2["boxes"], groups):
    patch.set_facecolor(color)
    patch.set_alpha(0.4)

ax4.set_xticklabels([g[0] for g in groups], fontsize=7)
ax4.set_ylabel("Risk Score (0–10)")
ax4.set_title("D  Risk Score by ΔpI Group")
ax4.axhline(3, color="#888", linestyle=":", linewidth=0.6, alpha=0.5)
ax4.axhline(6, color="#c0392b", linestyle=":", linewidth=0.6, alpha=0.5)
for i, data in enumerate(risk_data):
    ax4.text(i + 1, np.max(data) + 0.3, f"μ={np.mean(data):.1f}", ha="center", fontsize=7, color="#555")

plt.tight_layout(pad=2)
fig.savefig("results/figS2_dpi_stratified.png", dpi=300)
plt.close(fig)
print("Saved: figS2_dpi_stratified.png")

# ── Print detailed stats for ΔpI ≥ 3.0 ──
print("\n=== ΔpI ≥ 3.0: Detailed Profiles ===")
g3 = groups[3][1]
for r in sorted(g3, key=lambda x: -float(x["pH5.5_delta_pI"])):
    ab = r["antibody"][:22]
    dpi = float(r["pH5.5_delta_pI"])
    risk = float(r["risk_score"])
    n1 = float(r["pH5.5_arm1_net_charge"])
    n2 = float(r["pH5.5_arm2_net_charge"])
    cai = float(r["pH5.5_charge_asymmetry_CAI"])
    pat = r["pH5.5_pair_pattern"]
    comp = "C" if float(r["pH7.4_complementarity"]) < 0 else "S"
    print(f"  {ab:22s}  dPI={dpi:.1f}  Risk={risk:.0f}  n1={n1:+.0f}  n2={n2:+.0f}  CAI={cai:.2f}  {pat}  {comp}")
