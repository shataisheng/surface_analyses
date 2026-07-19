#!/usr/bin/env python3
"""
PEP-Patch 统一分析器
====================
自动检测输入类型（NPZ 疏水性 / PLY+CSV 静电势），
生成统一的残基级别 patch 详细报告。

用法:
    # 自动检测（给定 PDB，自动查找同目录 NPZ/PLY/CSV）
    python unified_analyzer.py protein.pdb

    # 指定输入文件
    python unified_analyzer.py --npz protein_out.npz
    python unified_analyzer.py --pdb protein.pdb --ply-pos es-pos.ply --ply-neg es-neg.ply --patches es_patches.csv

    # 导入使用
    from unified_analyzer import analyze
    result, stem = analyze(pdb_path="protein.pdb")

输出:
    {stem}_residues_detailed.csv  — 每个残基在每个 patch 中的面积占比
    {stem}_patch_summary.csv      — patch 级别摘要
"""

import argparse, csv, glob, os, sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


# ─── 文件自动检测 ──────────────────────────────────────
def strip_suffixes(name: str) -> str:
    """去除常见输出后缀，还原基础名称。"""
    for sfx in ["_fixed", "_out", "_output", "_patches"]:
        if name.endswith(sfx):
            return name[:-len(sfx)]
    return name


def detect_inputs(pdb_path: str) -> dict:
    """给定 PDB 路径，自动在同目录查找相关输出文件。"""
    pdb_path = os.path.abspath(pdb_path)
    parent = os.path.dirname(pdb_path)
    stem = strip_suffixes(Path(pdb_path).stem)

    # Build candidate stems: all possible prefixes derived from stripping
    candidates = [stem]
    for sfx in ["_out", "_fixed_out", ""]:
        c = stem + sfx if sfx else stem
        if c not in candidates:
            candidates.append(c)

    def find_one(pattern: str) -> str | None:
        # First try matching with candidate stems
        for c in candidates:
            matches = glob.glob(os.path.join(parent, f"{c}{pattern}"))
            if matches:
                return matches[0]
        # Fallback: any matching file in directory
        matches = glob.glob(os.path.join(parent, f"*{pattern}"))
        return matches[0] if matches else None

    return {
        "pdb": pdb_path,
        "stem": stem,
        "npz": find_one("_out.npz") or find_one(".npz"),
        "patches_csv": find_one("_patches.csv"),
        "ply_pos": find_one("-pos.ply"),
        "ply_neg": find_one("-neg.ply"),
        "ply_potential": find_one("-potential.ply"),
    }


# ─── PDB 解析 ──────────────────────────────────────────
def load_pdb_atom_map(pdb_file: str) -> tuple:
    """解析 PDB，返回 (atom_df, atom_coords_Angstrom)。"""
    import mdtraj as md
    from Bio.PDB import PDBParser

    traj = md.load(pdb_file)
    top = traj.topology

    parser = PDBParser(QUIET=True)
    bio_struct = parser.get_structure("p", pdb_file)
    bio_residues = []
    for chain in sorted(bio_struct[0], key=lambda c: c.id):
        for res in chain:
            if res.id[0] != " ":
                continue
            rseq = res.id[1]
            icode = res.id[2].strip()
            res_seq_str = f"{rseq}{icode}" if icode else str(rseq)
            bio_residues.append((chain.id, res.resname, res_seq_str))

    seq_nr_map = {}
    records = []
    for atom in top.atoms:
        res = atom.residue
        ri = res.index
        if ri < len(bio_residues):
            chain_id, res_name, res_seq_str = bio_residues[ri]
        else:
            chain_id = chr(ord("A") + res.chain.index)
            res_name = res.name
            res_seq_str = str(res.resSeq)
        key = (res.chain.index, res.resSeq)
        if key not in seq_nr_map:
            seq_nr_map[key] = len(seq_nr_map) + 1
        seq_nr = seq_nr_map[key]
        records.append({
            "atom_idx": atom.index,
            "atom_name": atom.name,
            "chain_id": chain_id,
            "res_name": res_name,
            "res_seq": res_seq_str,
            "res_id": f"{res_name}{res_seq_str}",
            "seq_nr": seq_nr,
            "seq_res_id": f"{res_name}{seq_nr}",
        })

    atom_df = pd.DataFrame(records)
    atom_coords = traj.xyz[0] * 10.0
    return atom_df, atom_coords


# ─── NPZ 分析（疏水性） ──────────────────────────────
def _build_atom_residue_map(pdb_path: str) -> list:
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("p", pdb_path)
    model = structure[0]
    atom_map = []
    gseq, prev = 0, None
    for chain in sorted(model, key=lambda c: c.id):
        for res in chain:
            if res.id[0] != " ":
                continue
            res_seq_raw = res.id[1]
            insertion = res.id[2].strip()
            res_seq_str = f"{res_seq_raw}{insertion}" if insertion else str(res_seq_raw)
            key = (chain.id, res.resname, res_seq_str)
            if key != prev:
                gseq += 1; prev = key
            rid = f"{res.resname}{res_seq_str}"
            for _ in res:
                atom_map.append((chain.id, res.resname, res_seq_str, gseq, rid))
    return atom_map


def analyze_npz(npz_path: str, pdb_path: str) -> list[dict]:
    npz = np.load(npz_path, allow_pickle=True)
    basename = next(k.replace(":vertices", "") for k in npz.keys() if ":vertices" in k)

    vertices = npz[f"{basename}:vertices"]
    faces = npz[f"{basename}:faces"]
    patch_arr = npz[f"{basename}:data:patch"]
    atom_arr = npz[f"{basename}:data:atom"]
    values_arr = npz.get(f"{basename}:data:values")

    ab = vertices[faces[:, 1]] - vertices[faces[:, 0]]
    ac = vertices[faces[:, 2]] - vertices[faces[:, 0]]
    tri_areas = np.sqrt(np.sum(np.cross(ab, ac) ** 2, axis=1)) / 2
    vert_areas = np.zeros(len(vertices), dtype=np.float64)
    np.add.at(vert_areas, faces.ravel(), np.repeat(tri_areas, 3))
    vert_areas /= 3

    atom_map = _build_atom_residue_map(pdb_path)
    unique_patches = sorted(set(patch_arr[patch_arr >= 0]))

    rows = []
    for pid in unique_patches:
        mask = patch_arr == pid
        per_vert_A2 = vert_areas[mask] * 100
        total_area_A2 = per_vert_A2.sum()

        ptype = "hydrophobic"
        if values_arr is not None:
            ptype = "hydrophobic" if values_arr[mask].mean() > 0 else "hydrophilic"

        centroid = vertices[mask].mean(axis=0)
        mean_dist_A = np.mean(np.linalg.norm(vertices[mask] - centroid, axis=1)) * 10

        res_areas = defaultdict(float)
        res_verts = defaultdict(int)
        res_info = {}

        for idx, i_atom in enumerate(np.flatnonzero(mask)):
            ai = atom_arr[i_atom]
            if 0 <= ai < len(atom_map):
                chain, rname, rseq, gseq, rid = atom_map[ai]
                key = (chain, rseq)
                res_areas[key] += per_vert_A2[idx]
                res_verts[key] += 1
                res_info[key] = (chain, rname, rseq, gseq, rid)

        for (chain, rseq), area in sorted(res_areas.items()):
            info = res_info[(chain, rseq)]
            frac = area / total_area_A2 if total_area_A2 > 0 else 0
            rows.append({
                "patch_nr": int(pid) + 1,
                "patch_type": ptype,
                "patch_total_area_A2": round(total_area_A2, 2),
                "chain_id": info[0],
                "res_name": info[1],
                "res_seq": info[2],
                "res_id": info[4],
                "seq_nr": info[3],
                "seq_res_id": info[4],
                "n_vertices": res_verts[(chain, rseq)],
                "area_A2": round(area, 3),
                "frac_of_patch": round(frac, 4),
                "mean_dist_A": round(mean_dist_A, 2),
            })
    return rows


# ─── PLY 分析（静电势） ──────────────────────────────
def _match_ply_colors(ply_colors, ply_coords_A, type_patches, atom_df, atom_coords) -> dict:
    from scipy.spatial import KDTree

    npoints_list = type_patches.npoints.tolist()
    nr_list = type_patches.nr.tolist()
    unique_colors, counts = np.unique(ply_colors, axis=0, return_counts=True)
    patch_colors = [
        (tuple(c), int(cnt)) for c, cnt in zip(unique_colors, counts)
        if not np.allclose(c, [256, 256, 256])
    ]
    color_to_nr, unresolved = {}, []

    for color, n_verts in patch_colors:
        candidates = [i for i, np_ in enumerate(npoints_list) if np_ == n_verts]
        if not candidates:
            candidates = [min(range(len(npoints_list)),
                             key=lambda i: abs(npoints_list[i] - n_verts))]
        if len(candidates) == 1:
            color_to_nr[color] = nr_list[candidates[0]]
        else:
            unresolved.append((color, n_verts, candidates))

    for color, n_verts, cands in unresolved:
        mask = np.all(ply_colors == np.array(color), axis=1)
        centroid = ply_coords_A[mask].mean(axis=0)
        best_idx, best_dist = None, float("inf")
        for ci in cands:
            patch_row = type_patches[type_patches.nr == nr_list[ci]].iloc[0]
            main_res_id = patch_row["main_residue"]
            matches = atom_df[atom_df["res_id"] == main_res_id]
            if len(matches) == 0:
                continue
            res_centroid = atom_coords[matches["atom_idx"].values].mean(axis=0)
            dist = np.linalg.norm(centroid - res_centroid)
            if dist < best_dist:
                best_dist, best_idx = dist, ci
        color_to_nr[color] = nr_list[best_idx if best_idx is not None else cands[0]]
    return color_to_nr


def analyze_ply(ply_path, patch_type, patches_df, atom_df, atom_coords, ply_scale=100.0) -> pd.DataFrame:
    import plyfile
    from scipy.spatial import KDTree

    ply = plyfile.PlyData.read(ply_path)
    verts, faces_data = ply["vertex"], ply["face"]
    coords = np.stack([verts["x"], verts["y"], verts["z"]], axis=1) * ply_scale
    colors = np.stack([verts["red"], verts["green"], verts["blue"]], axis=1)

    face_arr = np.vstack([faces_data["vertex_indices"][i] for i in range(len(faces_data))])
    e1 = coords[face_arr[:, 1]] - coords[face_arr[:, 0]]
    e2 = coords[face_arr[:, 2]] - coords[face_arr[:, 0]]
    face_areas = 0.5 * np.linalg.norm(np.cross(e1, e2), axis=1)
    vertex_areas = np.zeros(len(coords))
    for i in range(3):
        np.add.at(vertex_areas, face_arr[:, i], face_areas / 3.0)

    type_patches = patches_df[patches_df.type == patch_type].sort_values("nr").reset_index(drop=True)
    if len(type_patches) == 0:
        return pd.DataFrame()

    tree = KDTree(atom_coords)
    color_to_nr = _match_ply_colors(colors, coords, type_patches, atom_df, atom_coords)

    results = []
    for color, patch_nr in color_to_nr.items():
        mask = np.all(colors == np.array(color), axis=1)
        patch_coords, patch_vareas = coords[mask], vertex_areas[mask]
        patch_row = type_patches[type_patches.nr == patch_nr].iloc[0]
        patch_area = float(patch_row["area"]) * 100.0

        dists, atom_idxs = tree.query(patch_coords, k=1)
        res_data = atom_df.iloc[atom_idxs].copy().reset_index(drop=True)
        res_data["vertex_area_A2"] = patch_vareas
        res_data["dist_A"] = dists

        res_summary = (
            res_data
            .groupby(["chain_id", "res_name", "res_seq", "res_id", "seq_nr", "seq_res_id"])
            .agg(
                n_vertices=("atom_idx", "count"),
                area_A2=("vertex_area_A2", "sum"),
                mean_dist_A=("dist_A", "mean"),
            )
            .reset_index()
            .sort_values("area_A2", ascending=False)
        )
        res_summary["patch_nr"] = patch_nr
        res_summary["patch_type"] = patch_type
        res_summary["patch_total_area_A2"] = round(patch_area, 2)
        res_summary["frac_of_patch"] = (res_summary["area_A2"] / res_summary["area_A2"].sum()).round(4)
        results.append(res_summary)

    if not results:
        return pd.DataFrame()
    return pd.concat(results, ignore_index=True)


# ─── 主入口 ────────────────────────────────────────────
OUTPUT_COLUMNS = [
    "patch_nr", "patch_type", "patch_total_area_A2",
    "chain_id", "res_name", "res_seq", "res_id",
    "seq_nr", "seq_res_id",
    "n_vertices", "area_A2", "frac_of_patch", "mean_dist_A",
]


def analyze(
    pdb_path: str = None,
    npz_path: str = None,
    patches_csv: str = None,
    ply_pos: str = None,
    ply_neg: str = None,
    stem: str = None,
) -> tuple:
    """统一分析入口。自动检测输入类型并返回 (DataFrame, stem)。"""
    if pdb_path:
        detected = detect_inputs(pdb_path)
        if not stem:
            stem = detected["stem"]
        npz_path = npz_path or detected["npz"]
        patches_csv = patches_csv or detected["patches_csv"]
        ply_pos = ply_pos or detected["ply_pos"]
        ply_neg = ply_neg or detected["ply_neg"]

    if not stem:
        if npz_path:
            stem = strip_suffixes(Path(npz_path).stem)
        elif pdb_path:
            stem = strip_suffixes(Path(pdb_path).stem)
        else:
            stem = "pep_patch"

    # Prefer PLY+CSV when PDB is provided, NPZ when NPZ is provided
    has_ply = (ply_pos or ply_neg) and patches_csv and os.path.exists(patches_csv)
    has_npz = npz_path and os.path.exists(npz_path)

    if has_ply and not has_npz:
        mode = "electrostatic"
    elif has_npz and not has_ply:
        mode = "hydrophobic"
    elif has_ply and has_npz:
        # Both available: prefer PLY+CSV (electrostatic has richer data)
        mode = "electrostatic"
    else:
        raise FileNotFoundError(
            "Cannot detect input type. Provide --npz (hydrophobic) "
            "or --pdb + --ply-pos/--ply-neg + --patches (electrostatic)."
        )

    if mode == "hydrophobic":
        # Auto-find PDB if not provided
        if not pdb_path:
            pdb_dir = os.path.dirname(npz_path)
            pdb_candidates = glob.glob(os.path.join(pdb_dir, "*.pdb"))
            if pdb_candidates:
                pdb_path = pdb_candidates[0]
        if not pdb_path:
            raise FileNotFoundError("Cannot find PDB for NPZ analysis. Use --pdb.")
        print(f"Mode: Hydrophobic (NPZ)")
        print(f"  NPZ: {npz_path}")
        print(f"  PDB: {pdb_path}")
        rows = analyze_npz(npz_path, pdb_path)
        result = pd.DataFrame(rows)
    else:  # electrostatic
        print(f"Mode: Electrostatic (PLY+CSV)")
        atom_df, atom_coords = load_pdb_atom_map(pdb_path)
        patches_df = pd.read_csv(patches_csv)
        parts = []
        for pp, ptype in [(ply_pos, "positive"), (ply_neg, "negative")]:
            if pp and os.path.exists(pp):
                parts.append(analyze_ply(pp, ptype, patches_df, atom_df, atom_coords))
        if not parts:
            raise FileNotFoundError("No valid PLY files found")
        result = pd.concat(parts, ignore_index=True)

    cols = [c for c in OUTPUT_COLUMNS if c in result.columns]
    result = result[cols].sort_values(["patch_type", "patch_nr", "area_A2"],
                                       ascending=[True, True, False])
    return result, stem


def print_summary(df: pd.DataFrame, stem: str = ""):
    """打印分析摘要。"""
    pos = df[df.patch_type.isin(["positive", "hydrophobic"])]
    neg = df[df.patch_type.isin(["negative", "hydrophilic"])]

    print(f"\n{'='*65}")
    print(f"  {stem}")
    print(f"  {len(pos.patch_nr.unique())} positive/hydrophobic  +  "
          f"{len(neg.patch_nr.unique())} negative/hydrophilic  patches")
    print(f"  {len(df)} residue entries")
    print(f"{'='*65}")

    for label, subset in [("positive/hydrophobic", pos), ("negative/hydrophilic", neg)]:
        if len(subset) == 0:
            continue
        print(f"\n  [{label}]:")
        for pn, grp in list(subset.groupby("patch_nr"))[:5]:
            total = grp.iloc[0]["patch_total_area_A2"]
            top = grp.iloc[0]
            print(f"    Patch {pn:>3d}  {total:>8.1f} A^2  "
                  f"<- {top['chain_id']}:{top['res_id']} ({top['frac_of_patch']*100:.0f}%)")
        if subset.patch_nr.nunique() > 5:
            print(f"    ... and {subset.patch_nr.nunique() - 5} more")


# ─── CLI ───────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description="PEP-Patch unified analyzer")
    p.add_argument("path", nargs="?", help="PDB or NPZ path (auto-detect related files)")
    p.add_argument("--pdb", help="PDB path")
    p.add_argument("--npz", help="NPZ path (hydrophobic)")
    p.add_argument("--patches", help="patches CSV path")
    p.add_argument("--ply-pos", help="positive PLY")
    p.add_argument("--ply-neg", help="negative PLY")
    p.add_argument("-o", "--output", help="output CSV path")
    p.add_argument("--prefix", help="output file prefix")
    args = p.parse_args()

    if args.path:
        if args.path.endswith(".npz"):
            args.npz = args.path
            # Auto-find PDB: try related names, then any PDB in same dir
            base = os.path.splitext(args.path)[0]
            for sfx in ["_out", "_fixed_out", ""]:
                cand = base.replace(sfx, "") + ".pdb"
                if os.path.exists(cand):
                    args.pdb = cand
                    break
            # Fallback: any PDB in same directory
            if not args.pdb:
                pdb_candidates = glob.glob(os.path.join(os.path.dirname(args.path), "*.pdb"))
                if pdb_candidates:
                    args.pdb = pdb_candidates[0]
        elif args.path.endswith(".pdb"):
            args.pdb = args.path
        else:
            args.pdb = args.path

    if not args.pdb and not args.npz:
        p.print_help()
        return

    try:
        result, stem = analyze(
            pdb_path=args.pdb, npz_path=args.npz,
            patches_csv=args.patches, ply_pos=args.ply_pos,
            ply_neg=args.ply_neg, stem=args.prefix,
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return
    except Exception as e:
        print(f"Error: {e}")
        import traceback; traceback.print_exc()
        return

    print_summary(result, stem)

    # 残基级详细 CSV
    out_csv = args.output or f"{stem}_residues_detailed.csv"
    result.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"\nSaved: {out_csv}")

    # Patch 级摘要 CSV
    base = args.output.replace(".csv", "") if args.output else stem
    if "patch_total_area_A2" in result.columns:
        patch_summary = (
            result.groupby(["patch_nr", "patch_type"])
            .agg(
                patch_total_area_A2=("patch_total_area_A2", "first"),
                n_residues=("res_id", "nunique"),
                top_residue=("res_id", "first"),
                top_frac=("frac_of_patch", "first"),
            )
            .reset_index()
            .sort_values(["patch_type", "patch_nr"])
        )
        summary_file = f"{base}_patch_summary.csv"
        patch_summary.to_csv(summary_file, index=False, encoding="utf-8-sig")
        print(f"Saved: {summary_file}")


if __name__ == "__main__":
    main()
