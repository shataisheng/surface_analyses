#!/usr/bin/env python3
"""
PEP-Patch 批量汇总
==================
扫描目录中的所有 PDB 分析结果，生成汇总 CSV 方便横向对比。

用法:
    python batch_summary.py                              # 扫描当前目录
    python batch_summary.py --dir ./my_results           # 指定目录
    python batch_summary.py --es-only                    # 只汇总静电
    python batch_summary.py --hb-only                    # 只汇总疏水

输出:
    batch_summary_es.csv   — 静电势批量汇总
    batch_summary_hb.csv   — 疏水性批量汇总
"""

import os, sys, re, csv, glob as globmod, argparse
from collections import defaultdict


def find_proteins(directory: str) -> list[str]:
    """Find all unique protein stems from output files."""
    stems = set()
    for f in os.listdir(directory):
        for pattern in ["_es_patches.csv", "_hb_out.npz", "_es-pos.ply", "_hb_residues_detailed.csv"]:
            if pattern in f:
                stem = f.replace(pattern, "")
                stems.add(stem)
                break
    return sorted(stems)


def parse_integrals_from_log(log_path: str) -> dict | None:
    """Parse electrostatic integrals from a run log file."""
    if not os.path.exists(log_path):
        return None
    with open(log_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Look for: Integrals (total, ++, +, -, --):
    #          21.9658 10.4051 27.1930 -5.2272 -1.6493
    match = re.search(r"Integrals.*?\n\s*([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)", content)
    if match:
        return {
            "integral_total": float(match.group(1)),
            "integral_high": float(match.group(2)),
            "integral_pos": float(match.group(3)),
            "integral_neg": float(match.group(4)),
            "integral_low": float(match.group(5)),
        }
    return None


def summarize_electrostatic(stem: str, directory: str) -> dict | None:
    """Summarize electrostatic results for one protein."""
    patches_csv = os.path.join(directory, f"{stem}_es_patches.csv")
    log_path = os.path.join(directory, f"{stem}_es_run.log")

    if not os.path.exists(patches_csv):
        return None

    pos_patches = []
    neg_patches = []
    has_cdr = False
    cdr_count = 0

    with open(patches_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ptype = row.get("type", "")
            area = float(row.get("area", 0)) * 100  # nm² → Å²
            value = float(row.get("value", 0))
            residue = row.get("main_residue", "")
            cdr = row.get("cdr", "")

            patch_info = {"nr": int(row["nr"]), "area_A2": area, "value": value, "main_residue": residue}

            if ptype == "positive":
                pos_patches.append(patch_info)
            elif ptype == "negative":
                neg_patches.append(patch_info)

            if cdr and cdr.strip().lower() == "true":
                has_cdr = True
                cdr_count += 1

    pos_patches.sort(key=lambda x: x["area_A2"], reverse=True)
    neg_patches.sort(key=lambda x: x["area_A2"], reverse=True)

    integrals = parse_integrals_from_log(log_path)

    summary = {
        "protein": stem,
        "n_positive_patches": len(pos_patches),
        "n_negative_patches": len(neg_patches),
        "total_positive_area_A2": round(sum(p["area_A2"] for p in pos_patches), 1),
        "total_negative_area_A2": round(sum(p["area_A2"] for p in neg_patches), 1),
        "top_pos_1_area_A2": round(pos_patches[0]["area_A2"], 1) if pos_patches else "",
        "top_pos_1_residue": pos_patches[0]["main_residue"] if pos_patches else "",
        "top_pos_2_area_A2": round(pos_patches[1]["area_A2"], 1) if len(pos_patches) > 1 else "",
        "top_pos_2_residue": pos_patches[1]["main_residue"] if len(pos_patches) > 1 else "",
        "top_pos_3_area_A2": round(pos_patches[2]["area_A2"], 1) if len(pos_patches) > 2 else "",
        "top_pos_3_residue": pos_patches[2]["main_residue"] if len(pos_patches) > 2 else "",
        "top_neg_1_area_A2": round(neg_patches[0]["area_A2"], 1) if neg_patches else "",
        "top_neg_1_residue": neg_patches[0]["main_residue"] if neg_patches else "",
        "top_neg_2_area_A2": round(neg_patches[1]["area_A2"], 1) if len(neg_patches) > 1 else "",
        "top_neg_2_residue": neg_patches[1]["main_residue"] if len(neg_patches) > 1 else "",
        "top_neg_3_area_A2": round(neg_patches[2]["area_A2"], 1) if len(neg_patches) > 2 else "",
        "top_neg_3_residue": neg_patches[2]["main_residue"] if len(neg_patches) > 2 else "",
    }

    if has_cdr:
        summary["n_cdr_patches"] = cdr_count

    if integrals:
        summary.update(integrals)

    return summary


def summarize_hydrophobic(stem: str, directory: str) -> dict | None:
    """Summarize hydrophobic results for one protein."""
    detailed_csv = os.path.join(directory, f"{stem}_hb_residues_detailed.csv")

    if not os.path.exists(detailed_csv):
        return None

    patches = defaultdict(lambda: {"area_A2": 0, "top_residue": "", "top_frac": 0})

    with open(detailed_csv, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pn = int(row["patch_nr"])
            area = float(row.get("patch_total_area_A2", 0))
            frac = float(row.get("frac_of_patch", 0))
            res = f"{row.get('chain_id','')}:{row.get('res_id','')}"

            patches[pn]["area_A2"] = area
            if frac > patches[pn]["top_frac"]:
                patches[pn]["top_residue"] = res
                patches[pn]["top_frac"] = frac

    sorted_patches = sorted(patches.values(), key=lambda x: x["area_A2"], reverse=True)

    summary = {
        "protein": stem,
        "n_patches": len(sorted_patches),
        "total_patch_area_A2": round(sum(p["area_A2"] for p in sorted_patches), 1),
        "mean_patch_area_A2": round(sum(p["area_A2"] for p in sorted_patches) / len(sorted_patches), 1) if sorted_patches else 0,
    }

    for i in range(min(3, len(sorted_patches))):
        summary[f"top_{i+1}_area_A2"] = round(sorted_patches[i]["area_A2"], 1)
        summary[f"top_{i+1}_residue"] = sorted_patches[i]["top_residue"]

    return summary


def write_csv(rows: list[dict], path: str):
    if not rows:
        print(f"  No data to write to {path}")
        return
    columns = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Saved: {path} ({len(rows)} proteins)")


def main():
    p = argparse.ArgumentParser(description="PEP-Patch batch summary generator")
    p.add_argument("--dir", default=os.getcwd(), help="Directory to scan (default: cwd)")
    p.add_argument("--es-only", action="store_true", help="Only electrostatic")
    p.add_argument("--hb-only", action="store_true", help="Only hydrophobic")
    p.add_argument("--output-dir", help="Output directory for summary CSVs")
    args = p.parse_args()

    directory = os.path.abspath(args.dir)
    out_dir = os.path.abspath(args.output_dir) if args.output_dir else directory

    print(f"Scanning: {directory}")
    stems = find_proteins(directory)
    print(f"Found {len(stems)} proteins\n")

    if not args.hb_only:
        print("=== Electrostatic Summary ===")
        es_rows = []
        for stem in stems:
            s = summarize_electrostatic(stem, directory)
            if s:
                es_rows.append(s)
        if es_rows:
            write_csv(es_rows, os.path.join(out_dir, "batch_summary_es.csv"))

    if not args.es_only:
        print("\n=== Hydrophobic Summary ===")
        hb_rows = []
        for stem in stems:
            s = summarize_hydrophobic(stem, directory)
            if s:
                hb_rows.append(s)
            else:
                # Try older naming (no _hb suffix, from CLI runs)
                s = summarize_hydrophobic_legacy(stem, directory)
                if s:
                    hb_rows.append(s)
        if hb_rows:
            write_csv(hb_rows, os.path.join(out_dir, "batch_summary_hb.csv"))

    total = max(len(set()) if not 'es_rows' in dir() else len(es_rows),
                len(set()) if not 'hb_rows' in dir() else len(hb_rows))
    print(f"\nDone. {total} proteins summarized.")


def summarize_hydrophobic_legacy(stem: str, directory: str) -> dict | None:
    """Try legacy naming (no _hb prefix, from manual CLI runs)."""
    for suffix in ["_residues_detailed.csv", "_out.npz"]:
        if os.path.exists(os.path.join(directory, f"{stem}{suffix}")):
            # Try using unified_analyzer
            try:
                from src.unified_analyzer import analyze
                npz = os.path.join(directory, f"{stem}_out.npz")
                if os.path.exists(npz):
                    result, _ = analyze(npz_path=npz, stem=stem)
                    patches = defaultdict(lambda: {"area_A2": 0, "top_residue": "", "top_frac": 0})
                    for _, row in result.iterrows():
                        pn = int(row["patch_nr"])
                        area = float(row.get("patch_total_area_A2", 0))
                        patches[pn]["area_A2"] = area
                    sorted_p = sorted(patches.values(), key=lambda x: x["area_A2"], reverse=True)
                    return {
                        "protein": stem, "n_patches": len(sorted_p),
                        "total_patch_area_A2": round(sum(p["area_A2"] for p in sorted_p), 1),
                    }
            except Exception:
                pass
    return None


if __name__ == "__main__":
    main()
