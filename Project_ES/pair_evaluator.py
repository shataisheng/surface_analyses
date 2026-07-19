#!/usr/bin/env python3
"""
BsAb 双臂配对评估器
====================
基于 134 个临床阶段双抗的电荷数据，对新抗体对进行配对可行性评分。

用法:
    from pair_evaluator import evaluate_pair
    result = evaluate_pair(
        arm1_net_charge_ph55=25.0, arm1_net_charge_ph74=18.0, arm1_pi=8.5,
        arm2_net_charge_ph55=20.0, arm2_net_charge_ph74=14.0, arm2_pi=8.2,
    )
    print(result.recommendation)
"""

from dataclasses import dataclass, field
import numpy as np

# ── 基准数据（来自 134 对临床双抗 + 268 个 Fv 结构域） ──
BENCHMARK = {
    "n_pairs": 134,
    "n_arms": 268,
    # Arm-level (268 Fv domains)
    "ph55_net_charge": {"mean": 10.5, "median": 0.0, "p5": 0, "p25": 0, "p75": 22, "p95": 42},
    "ph74_net_charge": {"mean": 12.1, "median": 15.1, "p5": -19, "p25": 1, "p75": 26, "p95": 35},
    # Pair-level (134 BsAbs)
    "delta_pi": {"mean": 1.5, "median": 1.2, "safe_max": 1.0, "risk_threshold": 2.0},
    "cai": {"mean": 0.57, "median": 0.46, "symmetric_max": 0.3, "asymmetric_min": 0.7},
    "risk_correlation": {"delta_pi_to_risk": 0.794, "cai_to_risk": 0.404},
}


@dataclass
class PairEvaluation:
    """双臂配对评估结果"""
    # Inputs
    arm1_ph55_charge: float
    arm1_ph74_charge: float
    arm1_pi: float
    arm2_ph55_charge: float
    arm2_ph74_charge: float
    arm2_pi: float

    # Computed
    delta_pi: float = 0.0
    cai_ph55: float = 0.0
    cai_ph74: float = 0.0
    complement_ph74: bool = False
    score: int = 0
    flags: list = field(default_factory=list)
    recommendation: str = ""

    def __post_init__(self):
        self._compute()

    def _compute(self):
        # ΔpI
        self.delta_pi = round(abs(self.arm1_pi - self.arm2_pi), 1)

        # CAI (Charge Asymmetry Index) — sign-aware
        # Same sign: raw difference / max magnitude
        # Opposite sign: magnitude balance (avoid penalizing complementary pairs)
        def compute_cai(n1, n2):
            same_sign = (n1 * n2) >= 0
            if same_sign:
                return abs(n1 - n2) / max(abs(n1), abs(n2), 1.0)
            else:
                return abs(abs(n1) - abs(n2)) / max(abs(n1), abs(n2), 1.0)

        self.cai_ph55 = round(compute_cai(self.arm1_ph55_charge, self.arm2_ph55_charge), 2)
        self.cai_ph74 = round(compute_cai(self.arm1_ph74_charge, self.arm2_ph74_charge), 2)

        # Complementarity (pH 7.4)
        self.complement_ph74 = (self.arm1_ph74_charge * self.arm2_ph74_charge) < 0

        # ── Scoring (0-10, lower is better) ──
        self.score = 0
        self.flags = []

        # ΔpI (strongest predictor, r=0.794 with risk)
        if self.delta_pi < 1.0:
            pass  # ideal
        elif self.delta_pi < 1.5:
            self.score += 1
        elif self.delta_pi < 2.0:
            self.score += 2
            self.flags.append("moderate_delta_pI")
        elif self.delta_pi < 3.0:
            self.score += 4
            self.flags.append("high_delta_pI")
        else:
            self.score += 6
            self.flags.append("extreme_delta_pI")

        # CAI asymmetry
        max_cai = max(self.cai_ph55, self.cai_ph74)
        if max_cai > 0.7:
            self.score += 2
            self.flags.append("high_asymmetry")

        # Extreme charge (outside P5-P95 of benchmark)
        for label, val_55, val_74 in [
            ("arm1", self.arm1_ph55_charge, self.arm1_ph74_charge),
            ("arm2", self.arm2_ph55_charge, self.arm2_ph74_charge),
        ]:
            if val_55 > BENCHMARK["ph55_net_charge"]["p95"]:
                self.score += 1
                self.flags.append(f"{label}_extreme_positive_ph55")
            if val_74 > BENCHMARK["ph74_net_charge"]["p95"]:
                self.score += 1
                self.flags.append(f"{label}_extreme_positive_ph74")
            if val_74 < BENCHMARK["ph74_net_charge"]["p5"]:
                self.score += 1
                self.flags.append(f"{label}_extreme_negative_ph74")

        # Complementary bonus (slightly favorable but not required)
        if self.complement_ph74:
            self.flags.append("complementary_charge")

        # ── Recommendation ──
        if self.score <= 1:
            self.recommendation = (
                f"✅ RECOMMENDED (score={self.score}/10). "
                f"ΔpI={self.delta_pi} is in the safe zone (<1.0). "
                f"This pair is comparable to the safest 40% of clinical BsAbs."
            )
        elif self.score <= 3:
            self.recommendation = (
                f"⚠️ ACCEPTABLE (score={self.score}/10). "
                f"ΔpI={self.delta_pi} is moderate. "
                f"Consider charge engineering if manufacturability issues arise."
            )
        elif self.score <= 5:
            self.recommendation = (
                f"⚠️ CAUTION (score={self.score}/10). "
                f"High ΔpI={self.delta_pi} or asymmetry detected. "
                f"Only 26% of clinical BsAbs score in this range."
            )
        else:
            self.recommendation = (
                f"❌ NOT RECOMMENDED (score={self.score}/10). "
                f"Extreme ΔpI={self.delta_pi} or out-of-range charge. "
                f"High manufacturability risk. Strongly consider redesign."
            )

    def summary(self) -> str:
        lines = [
            f"{'='*55}",
            f"  BsAb Pair Evaluator (benchmark: 134 clinical BsAbs)",
            f"{'='*55}",
            f"  Arm1: net_charge pH5.5={self.arm1_ph55_charge:+.0f}  pH7.4={self.arm1_ph74_charge:+.0f}  pI={self.arm1_pi:.1f}",
            f"  Arm2: net_charge pH5.5={self.arm2_ph55_charge:+.0f}  pH7.4={self.arm2_ph74_charge:+.0f}  pI={self.arm2_pi:.1f}",
            f"  ---",
            f"  ΔpI = {self.delta_pi}  |  CAI(pH5.5) = {self.cai_ph55}  |  CAI(pH7.4) = {self.cai_ph74}",
            f"  Complement(pH7.4) = {self.complement_ph74}  |  Flags: {' '.join(self.flags) if self.flags else 'none'}",
            f"  ---",
            f"  Score: {self.score}/10",
            f"  {self.recommendation}",
            f"{'='*55}",
        ]
        return "\n".join(lines)


def evaluate_pair(
    arm1_net_charge_ph55: float,
    arm1_net_charge_ph74: float,
    arm1_pi: float,
    arm2_net_charge_ph55: float,
    arm2_net_charge_ph74: float,
    arm2_pi: float,
) -> PairEvaluation:
    """
    评估两个抗体的配对可行性。

    Parameters:
        arm1_net_charge_ph55: Arm1 在 pH 5.5 下的净电荷 (integral_total, kT/e)
        arm1_net_charge_ph74: Arm1 在 pH 7.4 下的净电荷
        arm1_pi: Arm1 的等电点 (Bjellqvist)
        arm2_net_charge_ph55, arm2_net_charge_ph74, arm2_pi: Arm2 同理

    Returns:
        PairEvaluation with score (0-10, lower=better) and recommendation.
    """
    return PairEvaluation(
        arm1_ph55_charge=arm1_net_charge_ph55,
        arm1_ph74_charge=arm1_net_charge_ph74,
        arm1_pi=arm1_pi,
        arm2_ph55_charge=arm2_net_charge_ph55,
        arm2_ph74_charge=arm2_net_charge_ph74,
        arm2_pi=arm2_pi,
    )


# ── CLI ───────────────────────────────────────────────
if __name__ == "__main__":
    # Demo: evaluate a hypothetical pair
    print(evaluate_pair(
        arm1_net_charge_ph55=25, arm1_net_charge_ph74=18, arm1_pi=8.5,
        arm2_net_charge_ph55=20, arm2_net_charge_ph74=14, arm2_pi=8.2,
    ).summary())
    print()
    # Demo: risky pair
    print(evaluate_pair(
        arm1_net_charge_ph55=45, arm1_net_charge_ph74=38, arm1_pi=9.5,
        arm2_net_charge_ph55=-5, arm2_net_charge_ph74=-12, arm2_pi=5.0,
    ).summary())
