#!/usr/bin/env python3
"""
BsAb 双臂电荷对比可视化
========================
从 enhanced_pair_analysis.csv 生成科研期刊风格的统计图。
"""

import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from collections import defaultdict

# ── 科研期刊风格设置 ──────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
})

# ── 配色 ──────────────────────────────────────────────
C = {
    "positive": "#e74c3c",
    "negative": "#3498db",
    "mixed": "#f39c12",
    "neutral": "#95a5a6",
    "ph55": "#e67e22",
    "ph74": "#2ecc71",
    "grid": "#e0e0e0",
}


def load_data(path="enhanced_pair_analysis.csv"):
    rows = []
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def sf(r, key, default=0.0):
    try: return float(r.get(key, default))
    except: return default


# ── Figure 4: 单臂分布分析 ──────────────────────────────────
def make_figure4(rows, out_dir):
    """Per-arm distributions: net charge, pI, charge type, pH response."""
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    ((ax1, ax2), (ax3, ax4)) = axes

    # 4a. 所有 268 个 Fv 臂在 pH5.5 & pH7.4 的净电荷分布（叠层直方图）
    arms_55 = [sf(r, "pH5.5_arm1_net_charge") for r in rows] + \
              [sf(r, "pH5.5_arm2_net_charge") for r in rows]
    arms_74 = [sf(r, "pH7.4_arm1_net_charge") for r in rows] + \
              [sf(r, "pH7.4_arm2_net_charge") for r in rows]
    bins = np.linspace(min(min(arms_55), min(arms_74)) - 2, max(max(arms_55), max(arms_74)) + 2, 35)
    ax1.hist(arms_55, bins=bins, color=C["ph55"], alpha=0.6, edgecolor="white", linewidth=0.3, label="pH 5.5")
    ax1.hist(arms_74, bins=bins, color=C["ph74"], alpha=0.6, edgecolor="white", linewidth=0.3, label="pH 7.4")
    ax1.axvline(0, color="gray", linestyle="--", linewidth=0.5)
    ax1.set_xlabel("Net Charge (kT/e)")
    ax1.set_ylabel("Number of Fv Arms")
    ax1.set_title("A  Net Charge Distribution (n=268 arms)")
    ax1.legend(frameon=False, fontsize=7)

    # 4b. pI (Bjellqvist) 分布（去重唯一臂）
    pi_vals = []
    seen = set()
    for r in rows:
        for key in ["pH5.5_arm1_pI_Bjellqvist", "pH5.5_arm2_pI_Bjellqvist"]:
            v = r.get(key, "")
            if v and v not in seen:
                try:
                    pi_vals.append(float(v))
                    seen.add(v)
                except ValueError:
                    pass
    ax2.hist(pi_vals, bins=25, color="#8e44ad", edgecolor="white", alpha=0.85, linewidth=0.3)
    ax2.axvline(7.4, color="#e74c3c", linestyle="--", linewidth=0.8, label="pI = 7.4 (physiological)")
    mean_pi = np.mean(pi_vals)
    ax2.axvline(mean_pi, color="#2ecc71", linestyle=":", linewidth=0.8, label=f"Mean = {mean_pi:.1f}")
    ax2.set_xlabel("Isoelectric Point (pI, Bjellqvist)")
    ax2.set_ylabel("Unique Fv Arms")
    ax2.set_title(f"B  pI Distribution (n={len(pi_vals)} unique arms)")
    ax2.legend(frameon=False, fontsize=7)

    # 4c. 单臂电荷类型组成（pH 5.5 & 7.4 并列条形）
    type_counts_55 = {"positive": 0, "negative": 0, "mixed": 0}
    type_counts_74 = {"positive": 0, "negative": 0, "mixed": 0}
    for r in rows:
        for key_55, key_74 in [("pH5.5_arm1_charge_type", "pH7.4_arm1_charge_type"),
                                ("pH5.5_arm2_charge_type", "pH7.4_arm2_charge_type")]:
            t55 = r.get(key_55, "mixed")
            t74 = r.get(key_74, "mixed")
            if t55 in type_counts_55: type_counts_55[t55] += 1
            if t74 in type_counts_74: type_counts_74[t74] += 1

    x = np.arange(3)
    w = 0.35
    bars1 = ax3.bar(x - w/2, [type_counts_55[t] for t in ["positive", "negative", "mixed"]],
                    w, label="pH 5.5", color=C["ph55"], alpha=0.85, edgecolor="white", linewidth=0.3)
    bars2 = ax3.bar(x + w/2, [type_counts_74[t] for t in ["positive", "negative", "mixed"]],
                    w, label="pH 7.4", color=C["ph74"], alpha=0.85, edgecolor="white", linewidth=0.3)
    ax3.set_xticks(x)
    ax3.set_xticklabels(["Positive", "Negative", "Mixed"])
    ax3.set_ylabel("Number of Fv Arms")
    ax3.set_title("C  Charge Type Composition (n=268 arms)")
    ax3.legend(frameon=False, fontsize=7)
    for bar in bars1:
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                 str(int(bar.get_height())), ha="center", fontsize=7)
    for bar in bars2:
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                 str(int(bar.get_height())), ha="center", fontsize=7)

    # 4d. pH 响应 Δnet = net(pH7.4) − net(pH5.5) 分布
    delta_nets = [sf(r, "pH_response_arm1") for r in rows] + \
                 [sf(r, "pH_response_arm2") for r in rows]
    ax4.hist(delta_nets, bins=30, color="#1abc9c", edgecolor="white", alpha=0.85, linewidth=0.3)
    ax4.axvline(0, color="gray", linestyle="--", linewidth=0.8, label="Δnet = 0")
    ax4.axvline(15, color="#e74c3c", linestyle=":", linewidth=0.8, label="|Δnet| > 15 (flag)")
    ax4.axvline(-15, color="#e74c3c", linestyle=":", linewidth=0.8)
    ax4.set_xlabel("Δnet = net(pH7.4) − net(pH5.5) (kT/e)")
    ax4.set_ylabel("Number of Fv Arms")
    ax4.set_title("D  pH-Dependent Charge Shift (n=268 arms)")
    ax4.legend(frameon=False, fontsize=7)

    plt.tight_layout(pad=2)
    fig.savefig(os.path.join(out_dir, "fig4_arm_distributions.png"), dpi=300)
    plt.close(fig)
    print("Saved: fig4_arm_distributions.png")


# ── Figure 5: 指标相关性 ──────────────────────────────────
def make_figure5(rows, out_dir):
    """Correlation analysis between metrics."""
    fig, axes = plt.subplots(2, 2, figsize=(8, 7.5))
    ((ax1, ax2), (ax3, ax4)) = axes

    delta_pis = np.array([sf(r, "pH5.5_delta_pI") for r in rows])
    cais = np.array([sf(r, "pH5.5_charge_asymmetry_CAI") for r in rows])
    risks = np.array([sf(r, "risk_score") for r in rows])
    ph_resp_max = np.array([max(abs(sf(r, "pH_response_arm1")), abs(sf(r, "pH_response_arm2")))
                            for r in rows])
    delta_hb_surf = np.array([sf(r, "delta_hb_surfscore") for r in rows])
    delta_gravy = np.array([sf(r, "delta_gravy") for r in rows])

    # 5a. ΔpI vs Risk Score
    from numpy.polynomial.polynomial import polyfit
    mask = ~np.isnan(delta_pis) & ~np.isnan(risks)
    x_d, y_d = delta_pis[mask], risks[mask]
    ax1.scatter(x_d, y_d, c="#5b9bd5", alpha=0.4, s=10, edgecolors="none")
    if len(x_d) > 2:
        b, m = polyfit(x_d, y_d, 1)
        xs = np.linspace(x_d.min(), x_d.max(), 50)
        ax1.plot(xs, b + m * xs, color="#e74c3c", linewidth=1.2)
    r_val = np.corrcoef(x_d, y_d)[0, 1]
    ax1.text(0.05, 0.95, f"r = {r_val:.3f}", transform=ax1.transAxes,
             fontsize=8, verticalalignment="top")
    ax1.set_xlabel("ΔpI (Bjellqvist)")
    ax1.set_ylabel("Risk Score")
    ax1.set_title("A  ΔpI vs Risk Score")

    # 5b. CAI vs Risk Score
    mask2 = ~np.isnan(cais) & ~np.isnan(risks)
    x_c, y_c = cais[mask2], risks[mask2]
    ax2.scatter(x_c, y_c, c="#f39c12", alpha=0.4, s=10, edgecolors="none")
    if len(x_c) > 2:
        b2, m2 = polyfit(x_c, y_c, 1)
        xs2 = np.linspace(x_c.min(), x_c.max(), 50)
        ax2.plot(xs2, b2 + m2 * xs2, color="#e74c3c", linewidth=1.2)
    r_val2 = np.corrcoef(x_c, y_c)[0, 1]
    ax2.text(0.05, 0.95, f"r = {r_val2:.3f}", transform=ax2.transAxes,
             fontsize=8, verticalalignment="top")
    ax2.set_xlabel("Charge Asymmetry Index (CAI)")
    ax2.set_ylabel("Risk Score")
    ax2.set_title("B  CAI vs Risk Score")

    # 5c. ΔpI vs CAI，按风险着色
    scatter = ax3.scatter(delta_pis, cais, c=risks, cmap="RdYlGn_r", alpha=0.55, s=12,
                          edgecolors="none", vmin=0, vmax=10)
    ax3.set_xlabel("ΔpI (Bjellqvist)")
    ax3.set_ylabel("Charge Asymmetry Index (CAI)")
    ax3.set_title("C  ΔpI vs CAI (colored by Risk)")
    cbar = plt.colorbar(scatter, ax=ax3, shrink=0.78)
    cbar.set_label("Risk Score", fontsize=7)

    # 5d. 相关性热图
    metrics = {
        "ΔpI": delta_pis,
        "CAI": cais,
        "Risk": risks,
        "|Δnet|max": ph_resp_max,
        "ΔGRAVY": delta_gravy,
        "ΔHBsurf": delta_hb_surf,
    }
    comp_vals = np.array([1 if int(sf(r, "pH5.5_complementarity")) < 0 else 0 for r in rows])
    metrics["Complementary"] = comp_vals

    labels = list(metrics.keys())
    n = len(labels)
    corr_mat = np.zeros((n, n))
    for i, li in enumerate(labels):
        for j, lj in enumerate(labels):
            a = metrics[li]
            b_mat = metrics[lj]
            mask_ij = ~np.isnan(a) & ~np.isnan(b_mat)
            if mask_ij.sum() > 2:
                corr_mat[i, j] = np.corrcoef(a[mask_ij], b_mat[mask_ij])[0, 1]

    im = ax4.imshow(corr_mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax4.set_xticks(range(n))
    ax4.set_yticks(range(n))
    ax4.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax4.set_yticklabels(labels, fontsize=7)
    ax4.set_title("D  Metric Correlation Matrix")
    for i in range(n):
        for j in range(n):
            ax4.text(j, i, f"{corr_mat[i, j]:.2f}", ha="center", va="center",
                     fontsize=6.5, color="white" if abs(corr_mat[i, j]) > 0.5 else "black")
    plt.colorbar(im, ax=ax4, shrink=0.78).set_label("Pearson r", fontsize=7)

    plt.tight_layout(pad=2)
    fig.savefig(os.path.join(out_dir, "fig5_correlations.png"), dpi=300)
    plt.close(fig)
    print("Saved: fig5_correlations.png")


# ── Figure 6: 风险分层 ────────────────────────────────────
def make_figure6(rows, out_dir):
    """Risk stratification analysis."""
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    ((ax1, ax2), (ax3, ax4)) = axes

    delta_pis = np.array([sf(r, "pH5.5_delta_pI") for r in rows])
    cais = np.array([sf(r, "pH5.5_charge_asymmetry_CAI") for r in rows])
    risks = np.array([sf(r, "risk_score") for r in rows])

    low_mask = risks < 3
    mod_mask = (risks >= 3) & (risks < 6)
    high_mask = risks >= 6

    # 6a. 按配对模式的 Risk Score 箱线图
    pattern_risk = defaultdict(list)
    for r in rows:
        raw = r.get("pH5.5_pair_pattern", "?/?")
        parts = raw.split("/")
        key = "/".join(sorted(parts))
        pattern_risk[key].append(sf(r, "risk_score"))
    pat_labels = sorted(pattern_risk.keys(), key=lambda k: -np.mean(pattern_risk[k]))
    pat_data = [pattern_risk[k] for k in pat_labels]
    bp = ax1.boxplot(pat_data, patch_artist=True, widths=0.6,
                     medianprops={"color": "black", "linewidth": 1},
                     flierprops={"marker": "o", "markersize": 3, "alpha": 0.5})
    box_colors = ["#e74c3c", "#3498db", "#f39c12", "#95a5a6", "#9b59b6"]
    for patch, c in zip(bp["boxes"], box_colors[:len(pat_labels)]):
        patch.set_facecolor(c)
        patch.set_alpha(0.35)
    ax1.set_xticklabels([l.replace("/", "/\n") for l in pat_labels], fontsize=6.5)
    ax1.set_ylabel("Risk Score")
    ax1.set_title("A  Risk Score by Charge Pattern")
    ax1.set_ylim(-0.5, 10.5)
    ax1.axhline(3, color="#e67e22", linestyle="--", linewidth=0.7, alpha=0.5)
    ax1.axhline(6, color="#c0392b", linestyle="--", linewidth=0.7, alpha=0.5)

    # 6b. ΔpI 按风险类别的 violin 图
    data_dpi = [delta_pis[low_mask], delta_pis[mod_mask], delta_pis[high_mask]]
    vp = ax2.violinplot(data_dpi, positions=[0, 1, 2], showmeans=True, showmedians=True,
                         widths=0.6)
    for i, body in enumerate(vp["bodies"]):
        body.set_facecolor(["#27ae60", "#e67e22", "#e74c3c"][i])
        body.set_alpha(0.5)
    ax2.set_xticks([0, 1, 2])
    ax2.set_xticklabels(["Low (<3)\nn=" + str(low_mask.sum()),
                         "Moderate (3–5)\nn=" + str(mod_mask.sum()),
                         "High (≥6)\nn=" + str(high_mask.sum())], fontsize=7)
    ax2.set_ylabel("ΔpI (Bjellqvist)")
    ax2.set_title("B  ΔpI by Risk Category")

    # 6c. CAI 按风险类别的 violin 图
    data_cai = [cais[low_mask], cais[mod_mask], cais[high_mask]]
    vp2 = ax3.violinplot(data_cai, positions=[0, 1, 2], showmeans=True, showmedians=True,
                          widths=0.6)
    for i, body in enumerate(vp2["bodies"]):
        body.set_facecolor(["#27ae60", "#e67e22", "#e74c3c"][i])
        body.set_alpha(0.5)
    ax3.set_xticks([0, 1, 2])
    ax3.set_xticklabels(["Low (<3)\nn=" + str(low_mask.sum()),
                         "Moderate (3–5)\nn=" + str(mod_mask.sum()),
                         "High (≥6)\nn=" + str(high_mask.sum())], fontsize=7)
    ax3.set_ylabel("CAI")
    ax3.set_title("C  CAI by Risk Category")

    # 6d. 风险类别组成饼图
    risk_counts = [low_mask.sum(), mod_mask.sum(), high_mask.sum()]
    risk_labels = [f"Low ({risk_counts[0]})", f"Moderate ({risk_counts[1]})", f"High ({risk_counts[2]})"]
    risk_colors_list = ["#27ae60", "#e67e22", "#e74c3c"]
    wedges, texts, autotexts = ax4.pie(risk_counts, labels=None, autopct="%1.1f%%",
        colors=risk_colors_list, pctdistance=0.7, textprops={"fontsize": 8})
    ax4.set_title("D  Risk Category Composition")
    ax4.legend(wedges, risk_labels, loc="center left", bbox_to_anchor=(1, 0.5),
               frameon=False, fontsize=7)

    plt.tight_layout(pad=2)
    fig.savefig(os.path.join(out_dir, "fig6_risk_stratification.png"), dpi=300)
    plt.close(fig)
    print("Saved: fig6_risk_stratification.png")


# ── Figure 7: 密度图 ────────────────────────────────────
def make_figure7(rows, out_dir):
    """Density plots: hexbin, KDE ridges."""
    from scipy.stats import gaussian_kde

    fig, axes = plt.subplots(2, 2, figsize=(8, 7.5))
    ((ax1, ax2), (ax3, ax4)) = axes

    delta_pis = np.array([sf(r, "pH5.5_delta_pI") for r in rows])
    cais = np.array([sf(r, "pH5.5_charge_asymmetry_CAI") for r in rows])
    risks = np.array([sf(r, "risk_score") for r in rows])
    gravy1 = np.array([sf(r, "arm1_gravy") for r in rows])
    gravy2 = np.array([sf(r, "arm2_gravy") for r in rows])
    hb_s1 = np.array([sf(r, "arm1_hb_surfscore_sum") for r in rows])
    hb_s2 = np.array([sf(r, "arm2_hb_surfscore_sum") for r in rows])
    delta_hb = np.array([sf(r, "delta_hb_surfscore") for r in rows])

    # 7a. ΔpI vs CAI 六边形密度图
    mask = ~np.isnan(delta_pis) & ~np.isnan(cais)
    hb = ax1.hexbin(delta_pis[mask], cais[mask], gridsize=25, cmap="YlOrRd",
                     mincnt=1, linewidths=0.1, edgecolors="#cccccc")
    ax1.set_xlabel("ΔpI (Bjellqvist)")
    ax1.set_ylabel("Charge Asymmetry Index (CAI)")
    ax1.set_title("A  Joint Density: ΔpI vs CAI")
    plt.colorbar(hb, ax=ax1, shrink=0.78).set_label("Count", fontsize=7)
    # Annotate quadrants
    ax1.axvline(1.0, color="#2ecc71", linestyle="--", linewidth=0.6, alpha=0.6)
    ax1.axhline(0.5, color="#2ecc71", linestyle="--", linewidth=0.6, alpha=0.6)
    ax1.text(0.3, 0.15, "Ideal", transform=ax1.transAxes, fontsize=7, color="#2ecc71",
             fontweight="bold")

    # 7b. Net charge pH5.5 vs pH7.4 密度 (all 268 arms)
    arms_55_all = np.array([sf(r, "pH5.5_arm1_net_charge") for r in rows] +
                           [sf(r, "pH5.5_arm2_net_charge") for r in rows])
    arms_74_all = np.array([sf(r, "pH7.4_arm1_net_charge") for r in rows] +
                           [sf(r, "pH7.4_arm2_net_charge") for r in rows])
    hb2 = ax2.hexbin(arms_55_all, arms_74_all, gridsize=30, cmap="Blues",
                      mincnt=1, linewidths=0.1, edgecolors="#cccccc")
    ax2.set_xlabel("Net Charge at pH 5.5 (kT/e)")
    ax2.set_ylabel("Net Charge at pH 7.4 (kT/e)")
    ax2.set_title("B  Joint Density: pH 5.5 vs 7.4 Net Charge (n=268 arms)")
    plt.colorbar(hb2, ax=ax2, shrink=0.78).set_label("Count", fontsize=7)
    # Unity line
    lims = [min(arms_55_all.min(), arms_74_all.min()) - 5,
            max(arms_55_all.max(), arms_74_all.max()) + 5]
    ax2.plot(lims, lims, "k--", linewidth=0.5, alpha=0.4)
    # Correlation annotation
    ax2.text(0.05, 0.95, f"r = {np.corrcoef(arms_55_all, arms_74_all)[0,1]:.3f}",
             transform=ax2.transAxes, fontsize=8, verticalalignment="top")

    # 7c. GRAVY density ridges by risk category
    low_mask = risks < 3
    mod_mask = (risks >= 3) & (risks < 6)
    high_mask = risks >= 6

    all_gravy = np.concatenate([gravy1, gravy2])
    x_grid = np.linspace(all_gravy.min() - 0.02, all_gravy.max() + 0.02, 200)
    colors_ridge = ["#27ae60", "#e67e22", "#e74c3c"]
    labels_ridge = [f"Low (n={2*low_mask.sum()})",
                    f"Moderate (n={2*mod_mask.sum()})",
                    f"High (n={2*high_mask.sum()})"]

    for i, (mask, color, label) in enumerate(
        [(low_mask, colors_ridge[0], labels_ridge[0]),
         (mod_mask, colors_ridge[1], labels_ridge[1]),
         (high_mask, colors_ridge[2], labels_ridge[2])]
    ):
        g_vals = np.concatenate([gravy1[mask], gravy2[mask]])
        if len(g_vals) > 2:
            kde = gaussian_kde(g_vals, bw_method=0.15)
            density = kde(x_grid)
            ax3.fill_between(x_grid, density + i * 0.5, i * 0.5,
                             alpha=0.35, color=color, linewidth=0.5, edgecolor=color)
            ax3.plot(x_grid, density + i * 0.5, color=color, linewidth=0.8, label=label)
    ax3.axvline(-0.310, color="#8e44ad", linestyle="--", linewidth=0.8,
                label="Ipilimumab (−0.31)")
    ax3.set_xlabel("GRAVY (Kyte-Doolittle)")
    ax3.set_yticks([])
    ax3.set_title("C  GRAVY Density by Risk Category (n=268 arms)")
    ax3.legend(frameon=False, fontsize=6.5, loc="upper left")

    # 7d. ΔHBsurf density by risk category
    x_grid2 = np.linspace(delta_hb.min() - 0.2, delta_hb.max() + 0.2, 200)
    for i, (mask, color, label) in enumerate(
        [(low_mask, colors_ridge[0], labels_ridge[0].replace("n=", "n=")),
         (mod_mask, colors_ridge[1], labels_ridge[1].replace("n=", "n=")),
         (high_mask, colors_ridge[2], labels_ridge[2].replace("n=", "n="))]
    ):
        vals = delta_hb[mask]
        if len(vals) > 2:
            kde = gaussian_kde(vals, bw_method=0.2)
            density = kde(x_grid2)
            ax4.fill_between(x_grid2, density + i * 0.3, i * 0.3,
                             alpha=0.35, color=color, linewidth=0.5, edgecolor=color)
            ax4.plot(x_grid2, density + i * 0.3, color=color, linewidth=0.8, label=label)
    ax4.set_xlabel("ΔHBsurf (|arm1 − arm2|, Crippen)")
    ax4.set_yticks([])
    ax4.set_title("D  Surface Hydrophobicity Asymmetry by Risk")
    ax4.legend(frameon=False, fontsize=6.5, loc="upper right")

    plt.tight_layout(pad=2)
    fig.savefig(os.path.join(out_dir, "fig7_density.png"), dpi=300)
    plt.close(fig)
    print("Saved: fig7_density.png")


# ── Figure 8: 配对双臂电荷对比 ──────────────────────────
def make_figure8(rows, out_dir):
    """Paired arm charge comparison: scatter, difference distribution, dumbbell."""
    fig, axes = plt.subplots(2, 2, figsize=(9, 8))
    ((ax1, ax2), (ax3, ax4)) = axes

    net1_55 = np.array([sf(r, "pH5.5_arm1_net_charge") for r in rows])
    net2_55 = np.array([sf(r, "pH5.5_arm2_net_charge") for r in rows])
    net1_74 = np.array([sf(r, "pH7.4_arm1_net_charge") for r in rows])
    net2_74 = np.array([sf(r, "pH7.4_arm2_net_charge") for r in rows])
    risks = np.array([sf(r, "risk_score") for r in rows])

    diff_55 = net1_55 - net2_55
    diff_74 = net1_74 - net2_74

    # 8a. pH 5.5: Arm1 vs Arm2 net charge scatter, colored by risk
    scatter_a = ax1.scatter(net1_55, net2_55, c=risks, cmap="RdYlGn_r",
                             alpha=0.55, s=14, edgecolors="none", vmin=0, vmax=10)
    lim55 = max(abs(net1_55).max(), abs(net2_55).max()) + 5
    ax1.plot([-lim55, lim55], [-lim55, lim55], "k--", linewidth=0.5, alpha=0.4)
    ax1.axhline(0, color="gray", linewidth=0.3)
    ax1.axvline(0, color="gray", linewidth=0.3)
    ax1.set_xlabel("Arm1 Net Charge (kT/e)")
    ax1.set_ylabel("Arm2 Net Charge (kT/e)")
    ax1.set_title("A  pH 5.5: Arm1 vs Arm2 Net Charge")
    plt.colorbar(scatter_a, ax=ax1, shrink=0.78).set_label("Risk", fontsize=7)
    # Quadrant counts
    q1 = ((net1_55 > 0) & (net2_55 > 0)).sum()
    q2 = ((net1_55 < 0) & (net2_55 > 0)).sum()
    q3 = ((net1_55 < 0) & (net2_55 < 0)).sum()
    q4 = ((net1_55 > 0) & (net2_55 < 0)).sum()
    ax1.text(0.98, 0.98, f"Q1: {q1}", transform=ax1.transAxes, fontsize=6,
             ha="right", va="top", color="#888888")
    ax1.text(0.02, 0.98, f"Q2: {q2}", transform=ax1.transAxes, fontsize=6,
             ha="left", va="top", color="#888888")
    ax1.text(0.02, 0.02, f"Q3: {q3}", transform=ax1.transAxes, fontsize=6,
             ha="left", va="bottom", color="#888888")
    ax1.text(0.98, 0.02, f"Q4: {q4}", transform=ax1.transAxes, fontsize=6,
             ha="right", va="bottom", color="#888888")

    # 8b. pH 7.4: Arm1 vs Arm2 net charge scatter
    scatter_b = ax2.scatter(net1_74, net2_74, c=risks, cmap="RdYlGn_r",
                             alpha=0.55, s=14, edgecolors="none", vmin=0, vmax=10)
    lim74 = max(abs(net1_74).max(), abs(net2_74).max()) + 5
    ax2.plot([-lim74, lim74], [-lim74, lim74], "k--", linewidth=0.5, alpha=0.4)
    ax2.axhline(0, color="gray", linewidth=0.3)
    ax2.axvline(0, color="gray", linewidth=0.3)
    ax2.set_xlabel("Arm1 Net Charge (kT/e)")
    ax2.set_ylabel("Arm2 Net Charge (kT/e)")
    ax2.set_title("B  pH 7.4: Arm1 vs Arm2 Net Charge")
    plt.colorbar(scatter_b, ax=ax2, shrink=0.78).set_label("Risk", fontsize=7)
    q1_74 = ((net1_74 > 0) & (net2_74 > 0)).sum()
    q2_74 = ((net1_74 < 0) & (net2_74 > 0)).sum()
    q3_74 = ((net1_74 < 0) & (net2_74 < 0)).sum()
    q4_74 = ((net1_74 > 0) & (net2_74 < 0)).sum()
    ax2.text(0.98, 0.98, f"Q1: {q1_74}", transform=ax2.transAxes, fontsize=6,
             ha="right", va="top", color="#888888")
    ax2.text(0.02, 0.98, f"Q2: {q2_74}", transform=ax2.transAxes, fontsize=6,
             ha="left", va="top", color="#888888")
    ax2.text(0.02, 0.02, f"Q3: {q3_74}", transform=ax2.transAxes, fontsize=6,
             ha="left", va="bottom", color="#888888")
    ax2.text(0.98, 0.02, f"Q4: {q4_74}", transform=ax2.transAxes, fontsize=6,
             ha="right", va="bottom", color="#888888")

    # 8c. Charge difference (arm1 − arm2) overlay for pH5.5 and pH7.4
    bins_diff = np.linspace(min(diff_55.min(), diff_74.min()) - 2,
                            max(diff_55.max(), diff_74.max()) + 2, 30)
    ax3.hist(diff_55, bins=bins_diff, color=C["ph55"], alpha=0.6,
             edgecolor="white", linewidth=0.3, label=f"pH 5.5 (μ={diff_55.mean():.1f})")
    ax3.hist(diff_74, bins=bins_diff, color=C["ph74"], alpha=0.6,
             edgecolor="white", linewidth=0.3, label=f"pH 7.4 (μ={diff_74.mean():.1f})")
    ax3.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax3.set_xlabel("Net Charge Difference: Arm1 − Arm2 (kT/e)")
    ax3.set_ylabel("Count")
    ax3.set_title("C  Pairwise Charge Difference Distribution")
    ax3.legend(frameon=False, fontsize=7)

    # 8d. Top 30 pairs by charge span (dumbbell chart)
    charge_span = np.abs(diff_55)
    top_idx = np.argsort(charge_span)[-30:][::-1]  # largest span first
    top_names = [rows[i]["antibody"][:18] for i in top_idx]

    y_pos = range(len(top_idx))
    for i, idx in enumerate(top_idx[::-1]):  # reverse for bottom-to-top
        actual_i = len(top_idx) - 1 - i
        ax4.plot([net1_55[top_idx[actual_i]], net2_55[top_idx[actual_i]]],
                 [i, i], color="#5b9bd5", linewidth=1.2, alpha=0.7, marker="o",
                 markersize=3, markerfacecolor="#e74c3c" if risks[top_idx[actual_i]] >= 4 else "#5b9bd5")
    ax4.axvline(0, color="gray", linewidth=0.5, alpha=0.5)
    ax4.set_yticks(range(len(top_idx)))
    ax4.set_yticklabels([top_names[j] for j in range(len(top_idx)-1, -1, -1)], fontsize=5.5)
    ax4.set_xlabel("Net Charge at pH 5.5 (kT/e)")
    ax4.set_title("D  Top 30 Pairs by Charge Span (Arm1 ●──● Arm2)")
    ax4.tick_params(axis="y", length=0)

    plt.tight_layout(pad=2)
    fig.savefig(os.path.join(out_dir, "fig8_paired_charge.png"), dpi=300)
    plt.close(fig)
    print("Saved: fig8_paired_charge.png")


# ── Figure 9: 全量双臂电荷跨度图 ────────────────────────
def make_figure9(rows, out_dir):
    """All 134 BsAb pairs as a sorted dumbbell chart with risk coloring."""
    net1_55 = np.array([sf(r, "pH5.5_arm1_net_charge") for r in rows])
    net2_55 = np.array([sf(r, "pH5.5_arm2_net_charge") for r in rows])
    risks = np.array([sf(r, "risk_score") for r in rows])
    dpis = np.array([sf(r, "pH5.5_delta_pI") for r in rows])
    names = np.array([r["antibody"][:22] for r in rows])

    charge_span = np.abs(net1_55 - net2_55)
    sort_idx = np.argsort(charge_span)

    fig, axes = plt.subplots(1, 2, figsize=(12, 16),
                              gridspec_kw={"width_ratios": [3, 1]})
    (ax1, ax2) = axes

    # ── 9a. Full 134-pair dumbbell chart ──
    n = len(rows)
    bar_height = 0.7

    for i, idx in enumerate(sort_idx):
        n1, n2 = net1_55[idx], net2_55[idx]
        risk = risks[idx]
        # Color: green→orange→red based on risk
        if risk < 3:
            color = "#27ae60"
        elif risk < 6:
            color = "#e67e22"
        else:
            color = "#e74c3c"

        # Subtle alternating row background
        if i % 2 == 0:
            ax1.axhspan(i - 0.5, i + 0.5, color="#f8f9fa", zorder=0, alpha=0.5)

        # Draw connecting line + endpoints
        ax1.plot([n1, n2], [i, i], color=color, linewidth=max(0.5, 1.8 - risk * 0.15),
                 alpha=0.75, zorder=2, solid_capstyle="round")
        ax1.scatter([n1], [i], s=max(6, 20 - risk * 1.5), color=color, zorder=3,
                    edgecolors="white", linewidth=0.3)
        ax1.scatter([n2], [i], s=max(6, 20 - risk * 1.5), color=color, zorder=3,
                    edgecolors="white", linewidth=0.3)

    ax1.axvline(0, color="gray", linewidth=0.6, alpha=0.4, linestyle="--")
    ax1.set_xlabel("Net Charge at pH 5.5 (kT/e)", fontsize=9)
    ax1.set_title(f"A  All {n} BsAb Pairs: Arm1 ●──● Arm2 Net Charge Span (pH 5.5)",
                  fontsize=10, fontweight="bold")

    # Label only selected pairs: top 5 span, bottom 5 span, some landmarks
    label_indices = set()
    label_indices.update(sort_idx[-5:])   # largest span
    label_indices.update(sort_idx[:5])    # smallest span
    # Plus every high-risk pair
    for i, idx in enumerate(sort_idx):
        if risks[idx] >= 6:
            label_indices.add(idx)

    for idx in label_indices:
        i = list(sort_idx).index(idx)
        n1, n2 = net1_55[idx], net2_55[idx]
        mid = (n1 + n2) / 2
        ax1.text(mid + 0.8, i, names[idx], fontsize=5, va="center",
                 color="#333333", fontweight="bold" if risks[idx] >= 6 else "normal")

    ax1.set_ylim(-1, n)
    ax1.set_yticks([])
    ax1.tick_params(axis="y", length=0)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color="#27ae60", lw=2, label="Low risk (<3)"),
        Line2D([0], [0], color="#e67e22", lw=2, label="Moderate (3–5)"),
        Line2D([0], [0], color="#e74c3c", lw=2, label="High risk (≥6)"),
    ]
    ax1.legend(handles=legend_elements, frameon=False, fontsize=7,
               loc="lower right")

    # ── 9b. Charge span distribution + summary stats ──
    span_bins = np.linspace(0, charge_span.max() + 2, 30)
    for mask, color, label in [
        (risks < 3, "#27ae60", f"Low (n={(risks<3).sum()})"),
        ((risks >= 3) & (risks < 6), "#e67e22", f"Mod (n={((risks>=3)&(risks<6)).sum()})"),
        (risks >= 6, "#e74c3c", f"High (n={(risks>=6).sum()})"),
    ]:
        ax2.hist(charge_span[mask], bins=span_bins, color=color, alpha=0.55,
                 edgecolor="white", linewidth=0.3, label=label)

    ax2.set_xlabel("Charge Span |Arm1 − Arm2| (kT/e)", fontsize=9)
    ax2.set_ylabel("Count", fontsize=9)
    ax2.set_title("B  Charge Span Distribution by Risk", fontsize=10, fontweight="bold")
    ax2.legend(frameon=False, fontsize=7)

    # Summary text box
    summary_text = (
        f"Mean span: {charge_span.mean():.1f} kT/e\n"
        f"Median span: {np.median(charge_span):.1f} kT/e\n"
        f"Span < 5: {(charge_span<5).sum()} ({(charge_span<5).sum()/n*100:.0f}%)\n"
        f"Span > 20: {(charge_span>20).sum()} ({(charge_span>20).sum()/n*100:.0f}%)\n"
        f"Max span: {charge_span.max():.1f}"
    )
    ax2.text(0.95, 0.95, summary_text, transform=ax2.transAxes, fontsize=7,
             va="top", ha="right", bbox=dict(boxstyle="round,pad=0.4",
             facecolor="#f0f0f0", alpha=0.7))

    plt.tight_layout(pad=2)
    fig.savefig(os.path.join(out_dir, "fig9_full_dumbbell.png"), dpi=300)
    plt.close(fig)
    print("Saved: fig9_full_dumbbell.png")


# ── Suppl. Figure S1: ΔpI & CAI 密度曲线 ─────────────────
def make_figure_s1(rows, out_dir):
    """Histogram + KDE + Normal fit for ΔpI and CAI distributions."""
    from scipy.stats import gaussian_kde, norm

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    (ax1, ax2) = axes

    delta_pis = np.array([sf(r, "pH5.5_delta_pI") for r in rows])
    cais = np.array([sf(r, "pH5.5_charge_asymmetry_CAI") for r in rows])

    # ── S1a: ΔpI density ──
    bins_dpi = np.linspace(0, delta_pis.max() + 0.5, 28)
    ax1.hist(delta_pis, bins=bins_dpi, density=True, color="#5b9bd5", alpha=0.5,
             edgecolor="white", linewidth=0.3, label="Histogram (n=134)")

    # KDE
    kde_dpi = gaussian_kde(delta_pis, bw_method=0.25)
    x_dpi = np.linspace(0, delta_pis.max() + 0.5, 200)
    ax1.plot(x_dpi, kde_dpi(x_dpi), color="#e74c3c", linewidth=1.8, label="KDE density")

    # Normal fit
    mu_dpi, sigma_dpi = norm.fit(delta_pis)
    ax1.plot(x_dpi, norm.pdf(x_dpi, mu_dpi, sigma_dpi), color="#2ecc71",
             linewidth=1.2, linestyle="--",
             label=f"Normal fit (μ={mu_dpi:.1f}, σ={sigma_dpi:.1f})")

    ax1.axvline(1.0, color="#e67e22", linestyle=":", linewidth=0.8, alpha=0.7)
    ax1.text(1.05, ax1.get_ylim()[1] * 0.85 if ax1.get_ylim()[1] > 0 else 0.3,
             "ΔpI=1.0", fontsize=7, color="#e67e22")
    ax1.axvline(2.0, color="#c0392b", linestyle=":", linewidth=0.8, alpha=0.7)
    ax1.text(2.05, ax1.get_ylim()[1] * 0.7 if ax1.get_ylim()[1] > 0 else 0.25,
             "ΔpI=2.0", fontsize=7, color="#c0392b")

    ax1.set_xlabel("ΔpI (Bjellqvist)")
    ax1.set_ylabel("Density")
    ax1.set_title("A  ΔpI Distribution with KDE & Normal Fit")
    ax1.legend(frameon=False, fontsize=6.5)

    # ── S1b: CAI density ──
    bins_cai = np.linspace(0, 1.05, 25)
    ax2.hist(cais, bins=bins_cai, density=True, color="#f39c12", alpha=0.5,
             edgecolor="white", linewidth=0.3, label="Histogram (n=134)")

    kde_cai = gaussian_kde(cais, bw_method=0.2)
    x_cai = np.linspace(0, 1.0, 200)
    ax2.plot(x_cai, kde_cai(x_cai), color="#e74c3c", linewidth=1.8, label="KDE density")

    mu_cai, sigma_cai = norm.fit(cais)
    ax2.plot(x_cai, norm.pdf(x_cai, mu_cai, sigma_cai), color="#2ecc71",
             linewidth=1.2, linestyle="--",
             label=f"Normal fit (μ={mu_cai:.2f}, σ={sigma_cai:.2f})")

    ax2.axvline(0.5, color="#e67e22", linestyle=":", linewidth=0.8, alpha=0.7)
    ax2.text(0.51, ax2.get_ylim()[1] * 0.85 if ax2.get_ylim()[1] > 0 else 1.0,
             "CAI=0.5", fontsize=7, color="#e67e22")

    ax2.set_xlabel("Charge Asymmetry Index (CAI)")
    ax2.set_ylabel("Density")
    ax2.set_title("B  CAI Distribution with KDE & Normal Fit")
    ax2.legend(frameon=False, fontsize=6.5)

    plt.tight_layout(pad=2)
    fig.savefig(os.path.join(out_dir, "figS1_dpi_cai_density.png"), dpi=300)
    plt.close(fig)
    print("Saved: figS1_dpi_cai_density.png")


# ── 主图：多面板 ───────────────────────────────────────
def make_figures(rows, out_dir="."):
    os.makedirs(out_dir, exist_ok=True)

    # ─── Figure 1: 综合四面板 ───
    fig, axes = plt.subplots(2, 2, figsize=(8, 7))
    ((ax1, ax2), (ax3, ax4)) = axes

    # 1a. ΔpI 分布直方图
    delta_pis = [sf(r, "pH5.5_delta_pI") for r in rows]
    ax1.hist(delta_pis, bins=30, color="#5b9bd5", edgecolor="white", alpha=0.85, linewidth=0.3)
    ax1.axvline(1.0, color="#e74c3c", linestyle="--", linewidth=0.8, label="ΔpI = 1.0")
    ax1.axvline(2.0, color="#c0392b", linestyle="-", linewidth=0.8, label="ΔpI = 2.0 (risk)")
    ax1.set_xlabel("ΔpI (|arm1 − arm2|)")
    ax1.set_ylabel("Count")
    ax1.set_title("A  Isoelectric Point Difference (Bjellqvist)")
    ax1.legend(frameon=False, fontsize=7)
    ax1.set_xlim(0, max(delta_pis) * 1.05)

    # 1b. 电荷配对模式饼图（合并镜像对 positive/mixed == mixed/positive）
    patterns_55 = defaultdict(int)
    for r in rows:
        raw = r.get("pH5.5_pair_pattern", "?/?")
        parts = raw.split("/")
        key = "/".join(sorted(parts))  # 合并镜像
        patterns_55[key] += 1
    labels = list(patterns_55.keys())
    sizes = list(patterns_55.values())
    colors_pie = ["#e74c3c", "#3498db", "#f39c12", "#95a5a6", "#9b59b6", "#1abc9c"][:len(labels)]
    wedges, texts, autotexts = ax2.pie(sizes, labels=None, autopct="%1.1f%%",
        colors=colors_pie, pctdistance=0.75, textprops={"fontsize": 7})
    ax2.set_title("B  pH 5.5 Pair Charge Pattern")
    ax2.legend(wedges, [l.replace("/", " / ") for l in labels],
               loc="center left", bbox_to_anchor=(1, 0.5), frameon=False, fontsize=6)

    # 1c. CAI 直方图
    cais = [sf(r, "pH5.5_charge_asymmetry_CAI") for r in rows]
    ax3.hist(cais, bins=30, color="#f39c12", edgecolor="white", alpha=0.85, linewidth=0.3)
    ax3.axvline(0.5, color="#c0392b", linestyle="--", linewidth=0.8, label="CAI = 0.5")
    ax3.set_xlabel("Charge Asymmetry Index (CAI)")
    ax3.set_ylabel("Count")
    ax3.set_title("C  Charge Asymmetry (pH 5.5)")
    ax3.legend(frameon=False, fontsize=7)

    # 1d. 风险评分分布
    risks = [int(sf(r, "risk_score")) for r in rows]
    risk_bins = np.arange(-0.5, 10.6, 1)
    ax4.hist(risks, bins=risk_bins, color="#8e44ad", edgecolor="white", alpha=0.85, linewidth=0.3)
    ax4.axvline(3, color="#e67e22", linestyle="--", linewidth=0.8, label="Moderate (≥3)")
    ax4.axvline(6, color="#c0392b", linestyle="-", linewidth=0.8, label="High (≥6)")
    ax4.set_xlabel("Risk Score (0–10)")
    ax4.set_ylabel("Count")
    ax4.set_title("D  Manufacturability Risk Score")
    ax4.legend(frameon=False, fontsize=7)
    ax4.set_xlim(-0.5, 10.5)

    plt.tight_layout(pad=2)
    fig.savefig(os.path.join(out_dir, "fig1_overview.png"), dpi=300)
    plt.close(fig)
    print("Saved: fig1_overview.png")

    # ─── Figure 2: pH 响应与互补性 ───
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.8))
    (ax1, ax2) = axes

    # 2a. pH5.5 vs pH7.4 净电荷散点（每个抗体两个点，arm1+arm2）
    ph55_charges = []
    ph74_charges = []
    for r in rows:
        ph55_charges.append(sf(r, "pH5.5_arm1_net_charge"))
        ph55_charges.append(sf(r, "pH5.5_arm2_net_charge"))
        ph74_charges.append(sf(r, "pH7.4_arm1_net_charge"))
        ph74_charges.append(sf(r, "pH7.4_arm2_net_charge"))

    ax1.scatter(ph55_charges, ph74_charges, c="#5b9bd5", alpha=0.3, s=8, edgecolors="none")
    # Unity line
    lims = [min(min(ph55_charges), min(ph74_charges)) - 5, max(max(ph55_charges), max(ph74_charges)) + 5]
    ax1.plot(lims, lims, "k--", linewidth=0.5, alpha=0.5)
    ax1.set_xlabel("Net Charge at pH 5.5 (kT/e)")
    ax1.set_ylabel("Net Charge at pH 7.4 (kT/e)")
    ax1.set_title("A  pH Response of Fv Arms (n=268)")
    # Annotate correlation
    ax1.text(0.05, 0.95, f"r = {np.corrcoef(ph55_charges, ph74_charges)[0,1]:.2f}",
             transform=ax1.transAxes, fontsize=8, verticalalignment="top")

    # 2b. 互补性条形图
    comp_55 = sum(1 for r in rows if int(sf(r, "pH5.5_complementarity")) < 0)
    same_55 = len(rows) - comp_55
    comp_74 = sum(1 for r in rows if int(sf(r, "pH7.4_complementarity")) < 0)
    same_74 = len(rows) - comp_74

    x = np.arange(2)
    width = 0.35
    ax2.bar(x - width/2, [same_55, comp_55], width, label="pH 5.5",
            color=C["ph55"], alpha=0.85, edgecolor="white", linewidth=0.3)
    ax2.bar(x + width/2, [same_74, comp_74], width, label="pH 7.4",
            color=C["ph74"], alpha=0.85, edgecolor="white", linewidth=0.3)
    ax2.set_xticks(x)
    ax2.set_xticklabels(["Same Sign\n(Repulsive)", "Opposite Sign\n(Complementary)"])
    ax2.set_ylabel("Number of Antibody Pairs")
    ax2.set_title("B  Electrostatic Complementarity")
    ax2.legend(frameon=False, fontsize=7)
    # Add count labels
    for bar, val in [(ax2.patches[0], same_55), (ax2.patches[1], comp_55),
                     (ax2.patches[2], same_74), (ax2.patches[3], comp_74)]:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                 str(val), ha="center", fontsize=7)

    plt.tight_layout(pad=2)
    fig.savefig(os.path.join(out_dir, "fig2_ph_response.png"), dpi=300)
    plt.close(fig)
    print("Saved: fig2_ph_response.png")

    # ─── Figure 3: 电荷密度与 Patch 集中度 ───
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.8))
    (ax1, ax2) = axes

    densities_arm1 = [sf(r, "arm1_charge_density") for r in rows]
    densities_arm2 = [sf(r, "arm2_charge_density") for r in rows]
    top3_1 = [sf(r, "arm1_top3_conc") for r in rows]
    top3_2 = [sf(r, "arm2_top3_conc") for r in rows]

    ax1.scatter(densities_arm1, densities_arm2, c="#8e44ad", alpha=0.35, s=8, edgecolors="none")
    ax1.axhline(0, color="gray", linewidth=0.3)
    ax1.axvline(0, color="gray", linewidth=0.3)
    ax1.set_xlabel("Arm1 Charge Density (kT/e·Å²)")
    ax1.set_ylabel("Arm2 Charge Density (kT/e·Å²)")
    ax1.set_title("A  Charge Density per Arm")

    ax2.scatter(top3_1, top3_2, c="#e67e22", alpha=0.35, s=8, edgecolors="none")
    ax2.axhline(0.9, color="gray", linestyle="--", linewidth=0.5)
    ax2.axvline(0.9, color="gray", linestyle="--", linewidth=0.5)
    ax2.set_xlabel("Arm1 Top-3 Patch Concentration")
    ax2.set_ylabel("Arm2 Top-3 Patch Concentration")
    ax2.set_title("B  Patch Concentration (Top-3 / Total)")
    ax2.set_xlim(0, 1.05); ax2.set_ylim(0, 1.05)

    plt.tight_layout(pad=2)
    fig.savefig(os.path.join(out_dir, "fig3_density.png"), dpi=300)
    plt.close(fig)
    print("Saved: fig3_density.png")

    # ─── Figure 4: 单臂分布分析 ───
    make_figure4(rows, out_dir)

    # ─── Figure 5: 指标相关性 ───
    make_figure5(rows, out_dir)

    # ─── Figure 6: 风险分层 ───
    make_figure6(rows, out_dir)

    # ─── Figure 7: 密度图 ───
    make_figure7(rows, out_dir)

    # ─── Figure 8: 配对双臂电荷对比 ───
    make_figure8(rows, out_dir)

    # ─── Figure 9: 全量双臂电荷跨度 ───
    make_figure9(rows, out_dir)

    # ─── Suppl. Figure S1: ΔpI & CAI 密度曲线 ───
    make_figure_s1(rows, out_dir)

    # ─── Stats summary ───
    print(f"\n{'='*55}")
    print(f"  Figure Summary ({len(rows)} antibody pairs)")
    print(f"{'='*55}")
    print(f"  ΔpI:   mean={np.mean(delta_pis):.1f}  median={np.median(delta_pis):.1f}")
    print(f"          >1.0: {sum(1 for d in delta_pis if d>1.0)} ({100*sum(1 for d in delta_pis if d>1.0)/len(rows):.0f}%)")
    print(f"          >2.0: {sum(1 for d in delta_pis if d>2.0)} ({100*sum(1 for d in delta_pis if d>2.0)/len(rows):.0f}%)")
    print(f"  CAI:   mean={np.mean(cais):.2f}  median={np.median(cais):.2f}")
    print(f"  Risk:  mean={np.mean(risks):.1f}  low={sum(1 for r in risks if r<3)}  mod={sum(1 for r in risks if 3<=r<6)}  high={sum(1 for r in risks if r>=6)}")
    print(f"  pH Response corr: {np.corrcoef(ph55_charges, ph74_charges)[0,1]:.3f}")
    print(f"  Complementary pairs: pH5.5={comp_55}  pH7.4={comp_74}")


if __name__ == "__main__":
    rows = load_data()
    print(f"Loaded {len(rows)} antibody pairs")
    make_figures(rows, out_dir="results")
