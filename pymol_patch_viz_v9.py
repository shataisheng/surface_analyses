#!/usr/bin/env python
"""
PyMOL Charge Patch 可视化脚本 (自运行版)
=========================================
在 PyMOL 中直接 run 此脚本，自动查找 detailed.csv 并着色。

自动发现规则 (按优先级):
  1. 如果传了参数: pymol_patch_viz(csv="...", obj="...")
  2. 在 PDB 所在目录查找 *_detailed.csv
  3. 在当前工作目录查找

用法:
  PyMOL> load Herceptin_fixed.pdb
  PyMOL> run D:/Python/pymol_patch_viz.py

CSV 必须包含列 (支持别名):
  patch_nr, patch_type, patch_total_area_A2 (或 patch_total_area),
  chain_id (或 chain), res_seq (或 resi)
"""

import os
import csv
import re
import sys
import glob
from collections import defaultdict, OrderedDict
from pathlib import Path

try:
    from pymol import cmd, stored
    IN_PYMOL = True
except ImportError:
    IN_PYMOL = False
    print("[pymol_patch_viz] 未在 PyMOL 内运行 (测试模式)")

try:
    import colorsys
    def _hsv_rgb(h, s, v):
        r, g, b = colorsys.hsv_to_rgb(h % 1.0, s, v)
        return [r, g, b]
except ImportError:
    # 简易 HSV->RGB fallback
    def _hsv_rgb(h, s, v):
        if s == 0: return [v, v, v]
        h = (h % 1.0) * 6
        i = int(h)
        f = h - i
        p = v * (1 - s)
        q = v * (1 - s * f)
        t = v * (1 - s * (1 - f))
        if i == 0: return [v, t, p]
        if i == 1: return [q, v, p]
        if i == 2: return [p, v, t]
        if i == 3: return [p, q, v]
        if i == 4: return [t, p, v]
        return [v, p, q]


COLUMN_MAP = {
    "patch_nr":             ["patch_nr"],
    "patch_type":           ["patch_type"],
    "patch_total_area_A2":  ["patch_total_area_A2", "patch_total_area"],
    "chain_id":             ["chain_id", "chain"],
    "res_seq":              ["res_seq", "res_id_num", "resi"],
}


def _resolve_columns(headers):
    resolved = {}
    headers_lower = [h.strip().lower() for h in headers]
    for std_name, aliases in COLUMN_MAP.items():
        for alias in aliases:
            if alias.lower() in headers_lower:
                idx = headers_lower.index(alias.lower())
                resolved[std_name] = headers[idx]
                break
    required = ["patch_nr", "patch_type", "patch_total_area_A2", "chain_id", "res_seq"]
    missing = [r for r in required if r not in resolved]
    if missing:
        raise KeyError(f"CSV 缺少列: {missing}。表头: {headers}")
    return resolved


def parse_detailed_csv(csv_path):
    patches = defaultdict(lambda: {"type": None, "total_area": 0.0, "residues": []})
    # 尝试多种编码
    used_enc = "utf-8-sig"
    for enc in ["utf-8-sig", "utf-8", "gbk", "cp1252", "latin-1"]:
        try:
            with open(csv_path, "r", encoding=enc) as probe:
                probe.read(100)
            used_enc = enc
            break
        except (UnicodeDecodeError, UnicodeError):
            continue
    with open(csv_path, "r", encoding=used_enc) as f:
        reader = csv.DictReader(f)
        cols = _resolve_columns(reader.fieldnames or [])
        for row in reader:
            pn = int(row[cols["patch_nr"]])
            info = patches[pn]
            if info["type"] is None:
                info["type"] = row[cols["patch_type"]].strip().lower()
                info["total_area"] = float(row[cols["patch_total_area_A2"]])
            chain = row[cols["chain_id"]].strip()
            res_seq = row[cols["res_seq"]].strip()
            info["residues"].append((chain, str(res_seq)))
    for pn in patches:
        patches[pn]["residues"] = list(OrderedDict.fromkeys(patches[pn]["residues"]))
    return dict(patches)



# 科学色板预设 (20点采样自 matplotlib)
COLORMAPS = {
    "viridis": [
        [0.267,0.005,0.329],[0.282,0.146,0.462],[0.250,0.274,0.533],
        [0.199,0.388,0.555],[0.156,0.490,0.558],[0.122,0.589,0.546],
        [0.162,0.687,0.499],[0.336,0.777,0.402],[0.586,0.847,0.250],
        [0.866,0.890,0.096],[0.993,0.906,0.144],
    ],
    "plasma": [
        [0.050,0.030,0.528],[0.261,0.013,0.618],[0.430,0.001,0.659],
        [0.584,0.069,0.633],[0.715,0.187,0.546],[0.820,0.307,0.448],
        [0.905,0.430,0.356],[0.967,0.564,0.265],[0.994,0.717,0.177],
        [0.971,0.888,0.146],[0.940,0.975,0.131],
    ],
    "inferno": [
        [0.001,0.000,0.014],[0.093,0.046,0.234],[0.271,0.041,0.412],
        [0.441,0.099,0.432],[0.609,0.159,0.394],[0.770,0.236,0.307],
        [0.898,0.359,0.189],[0.973,0.530,0.053],[0.985,0.728,0.121],
        [0.946,0.931,0.442],[0.988,0.998,0.645],
    ],
    "magma": [
        [0.001,0.000,0.014],[0.083,0.056,0.221],[0.246,0.059,0.448],
        [0.415,0.110,0.505],[0.582,0.172,0.503],[0.755,0.229,0.463],
        [0.909,0.325,0.385],[0.982,0.506,0.372],[0.997,0.705,0.483],
        [0.992,0.899,0.654],[0.987,0.991,0.750],
    ],
    "YlOrRd": [
        [1.000,1.000,0.800],[0.997,0.877,0.518],[0.996,0.773,0.379],
        [0.995,0.646,0.276],[0.991,0.503,0.221],[0.983,0.295,0.162],
        [0.900,0.123,0.115],[0.779,0.026,0.139],[0.600,0.000,0.149],
        [0.502,0.000,0.149],
    ],
}


def _sample_cmap(cmap_name, n_colors):
    cmap = COLORMAPS.get(cmap_name)
    if cmap is None:
        return None
    if n_colors <= len(cmap):
        idx = [int(i*(len(cmap)-1)/max(n_colors-1,1)) for i in range(n_colors)]
        return [cmap[i] for i in idx]
    result = []
    for i in range(n_colors):
        f = i/max(n_colors-1,1)*(len(cmap)-1)
        lo, hi = int(f), min(int(f)+1, len(cmap)-1)
        t = f-lo
        result.append([cmap[lo][j]+(cmap[hi][j]-cmap[lo][j])*t for j in range(3)])
    return result


def _compute_colors(patch_list, hue_start, hue_end,
                    v_range=(0.35, 0.75), s_range=(0.85, 0.25),
                    cmap=None):
    if not patch_list:
        return {}
    if cmap and cmap in COLORMAPS:
        sampled = _sample_cmap(cmap, len(patch_list))
        return {pn: rgb for (pn, _), rgb in zip(patch_list, sampled)}
    areas = [info["total_area"] for _, info in patch_list]
    max_a, min_a = max(areas), min(areas)
    span = max_a - min_a if max_a > min_a else 1.0
    v_min, v_max = v_range
    s_max, s_min = s_range
    colors = {}
    for pn, info in patch_list:
        frac = (info["total_area"] - min_a) / span
        hue = hue_end - (hue_end - hue_start) * frac
        value = v_max - (v_max - v_min) * frac
        saturation = s_min + (s_max - s_min) * frac
        colors[pn] = _hsv_rgb(hue, saturation, value)
    return colors


def _find_detailed_csv(pdb_path=None):
    search_dirs = []
    if pdb_path:
        pdb_dir = os.path.dirname(os.path.abspath(pdb_path))
        search_dirs.append(pdb_dir)
    search_dirs.append(os.getcwd())
    patterns = ["*_detailed.csv", "*detailed*.csv",
                "*patch*residue*detailed*.csv", "*patch_residue*.csv"]
    for sdir in search_dirs:
        for pat in patterns:
            matches = glob.glob(os.path.join(sdir, pat))
            if matches:
                non_sasa = [m for m in matches
                            if "sasa" not in os.path.basename(m).lower()]
                candidates = non_sasa if non_sasa else matches
                candidates.sort(key=lambda x: len(os.path.basename(x)))
                print(f"[pymol_patch_viz] CSV: {candidates[0]}")
                return candidates[0]
    return None


def _resi_sort_key(r):
    s = str(r)
    num = ""
    for c in s:
        if c.isdigit(): num += c
        else: break
    return int(num) if num else 0


def _build_selection(obj, residues, exclude=None):
    chain_resi = defaultdict(set)
    for chain, resi in residues:
        chain_resi[chain].add(str(resi))
    # 去重：排除已在更大 patch 中着色的残基
    if exclude:
        for chain in chain_resi:
            chain_resi[chain] -= exclude.get(chain, set())
    parts = []
    for chain in sorted(chain_resi):
        if not chain_resi[chain]:
            continue
        resi_str = "+".join(sorted(chain_resi[chain], key=_resi_sort_key))
        parts.append(f"({obj} and chain {chain} and resi {resi_str})")
    return " or ".join(parts)


def pymol_patch_viz(csv_path=None, object_name=None,
                    pos_hue=(0.0, 0.16), neg_hue=(0.55, 0.60),
                    cartoon_transparency=0.55, surface_transparency=0.30,
                    show_spheres=True, sphere_scale=0.5,
                    show_surface=True, bg_color="white",
                    pos_cmap=None, neg_cmap=None):
    if not IN_PYMOL:
        print("[pymol_patch_viz] 错误: 需要在 PyMOL 内运行")
        return

    obj_list = cmd.get_object_list("all")
    if not obj_list:
        print("[pymol_patch_viz] 错误: 没有加载任何对象，请先 load PDB")
        return

    if object_name and object_name in obj_list:
        obj = object_name
    else:
        obj = obj_list[-1]
        if object_name and object_name not in obj_list:
            print(f"[pymol_patch_viz] 对象 '{object_name}' 未找到，使用 '{obj}'")

    print(f"[pymol_patch_viz] 对象: {obj}")

    if csv_path is None:
        pdb_path = None
        try:
            pdb_path = cmd.get_object_file(obj)
        except Exception:
            pass
        csv_path = _find_detailed_csv(pdb_path)

    if csv_path is None:
        print("[pymol_patch_viz] 错误: 找不到 detailed CSV")
        print("  请在 PDB 目录放 *_detailed.csv 或手动指定")
        return

    if not os.path.exists(csv_path):
        print(f"[pymol_patch_viz] 错误: {csv_path} 不存在")
        return

    print(f"[pymol_patch_viz] 解析: {csv_path}")
    try:
        patches = parse_detailed_csv(csv_path)
    except Exception as e:
        print(f"[pymol_patch_viz] 解析失败: {e}")
        return

    pos_patches = [(pn, i) for pn, i in patches.items() if i["type"] == "positive"]
    neg_patches = [(pn, i) for pn, i in patches.items() if i["type"] == "negative"]

    if not pos_patches and not neg_patches:
        print("[pymol_patch_viz] 警告: 无 positive/negative patch")
        return

    pos_sorted = sorted(pos_patches, key=lambda x: x[1]["total_area"], reverse=True)
    neg_sorted = sorted(neg_patches, key=lambda x: x[1]["total_area"], reverse=True)
    print(f"[pymol_patch_viz] {len(pos_sorted)} pos + {len(neg_sorted)} neg")

    pos_colors = _compute_colors(pos_sorted, pos_hue[0], pos_hue[1], cmap=pos_cmap)
    neg_colors = _compute_colors(neg_sorted, neg_hue[0], neg_hue[1], cmap=neg_cmap)

    for pn, rgb in pos_colors.items():
        cmd.set_color(f"pos_patch_{pn}", list(rgb))
    for pn, rgb in neg_colors.items():
        cmd.set_color(f"neg_patch_{pn}", list(rgb))

    cmd.hide("everything", obj)
    cmd.show("cartoon", obj)
    cmd.color("grey80", obj)
    cmd.set("cartoon_transparency", cartoon_transparency)

    if show_surface:
        cmd.show("surface", obj)
        cmd.set("surface_quality", 1)
        cmd.set("transparency", surface_transparency, obj)
        cmd.color("grey80", obj)
    else:
        cmd.hide("surface", obj)

    cmd.group("positive_patches", "", "add")
    colored_residues = defaultdict(set)  # chain -> set of resi
    for pn, info in pos_sorted:
        sel_name = f"patch_pos_{pn}"
        sel_expr = _build_selection(obj, info["residues"], exclude=colored_residues)
        try:
            if sel_expr:
                cmd.select(sel_name, sel_expr)
                cmd.color(f"pos_patch_{pn}", sel_name)
                if show_spheres:
                    cmd.show("sphere", sel_name)
                    cmd.set("sphere_scale", sphere_scale, sel_name)
            else:
                # 残基全被排除，仍创建空 selection 供 all_patches 引用
                cmd.select(sel_name, "none")
            cmd.group("positive_patches", sel_name, "add")
        except Exception as e:
            print(f"  [跳过] pos patch {pn}: {e}")
        for chain, resi in info["residues"]:
            colored_residues[chain].add(str(resi))

    cmd.group("negative_patches", "", "add")
    colored_residues = defaultdict(set)
    for pn, info in neg_sorted:
        sel_name = f"patch_neg_{pn}"
        sel_expr = _build_selection(obj, info["residues"], exclude=colored_residues)
        try:
            if sel_expr:
                cmd.select(sel_name, sel_expr)
                cmd.color(f"neg_patch_{pn}", sel_name)
                if show_spheres:
                    cmd.show("sphere", sel_name)
                    cmd.set("sphere_scale", sphere_scale, sel_name)
            else:
                cmd.select(sel_name, "none")
            cmd.group("negative_patches", sel_name, "add")
        except Exception as e:
            print(f"  [跳过] neg patch {pn}: {e}")
        for chain, resi in info["residues"]:
            colored_residues[chain].add(str(resi))

    cmd.zoom(obj)
    cmd.bg_color(bg_color)
    cmd.set("ray_shadow", "off")
    cmd.set("antialias", 2)
    cmd.set("depth_cue", 0)

    # 拼接 all_patches 为原子级 selection (避免 group 嵌套 hide 失效)
    all_names = [f"patch_pos_{pn}" for pn, _ in pos_sorted]
    all_names += [f"patch_neg_{pn}" for pn, _ in neg_sorted]
    if all_names:
        try:
            cmd.select("all_patches", " or ".join(all_names))
        except Exception:
            pass

    print(f"[pymol_patch_viz] 完成: {len(pos_sorted)} pos + {len(neg_sorted)} neg")
    print(f"  分组:")
    print(f"    positive_patches / negative_patches  -- patch 残基")
    print(f"    all_patches     -- 全部 patch (一键操作)")
    if show_spheres:
        print(f"  隐藏小球:  hide spheres, all_patches")
        print(f"  显示小球:  show spheres, all_patches")
    print(f"  隐藏表面:  hide surface")
    print(f"  显示 ribbon:  hide cartoon; show ribbon")


if IN_PYMOL:
    try:
        pymol_patch_viz()
    except Exception as e:
        print(f"[pymol_patch_viz] 自动执行失败: {e}")
        print("  可手动: pymol_patch_viz('detailed.csv', 'object_name')")
