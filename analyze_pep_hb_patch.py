#!/usr/bin/env python3
"""
pep-patch npz 完整分析器
========================
从 npz (含 patch 数据) + PDB 生成详细的 patch-per-residue 报告。

输出格式与 herceptin_patch_residues_detailed.csv 一致。

用法:
  python analyze_pep_patch.py B7H3_out.npz
  python analyze_pep_patch.py --npz out.npz --pdb protein.pdb -o detailed.csv
"""

import csv, os, re, sys, argparse, numpy as np
from collections import defaultdict


def build_atom_residue_map(pdb_path):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("p", pdb_path)
    model = structure[0]
    atom_map = []
    gseq = 0
    prev = None
    for chain in sorted(model, key=lambda c: c.id):
        for res in chain:
            if res.id[0] != " ": continue
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


def triangle_areas(vertices, faces):
    ab = vertices[faces[:, 1]] - vertices[faces[:, 0]]
    ac = vertices[faces[:, 2]] - vertices[faces[:, 0]]
    cross = np.cross(ab, ac)
    return np.sqrt(np.sum(cross ** 2, axis=1)) / 2


def vertex_areas(faces, tri_areas, n_vertices):
    areas = np.zeros(n_vertices, dtype=np.float64)
    np.add.at(areas, faces.ravel(), np.repeat(tri_areas, 3))
    return areas / 3


def analyze_npz(npz_path, pdb_path, pep_compat=False):
    npz = np.load(npz_path, allow_pickle=True)
    basename = next(k.replace(":vertices", "") for k in npz.keys() if ":vertices" in k)

    vertices = npz[f"{basename}:vertices"]
    faces = npz[f"{basename}:faces"]
    patch_arr = npz[f"{basename}:data:patch"]
    atom_arr = npz[f"{basename}:data:atom"]
    values_arr = npz.get(f"{basename}:data:values")

    nv = len(vertices)
    tri_areas = triangle_areas(vertices, faces)      # nm^2
    vert_areas = vertex_areas(faces, tri_areas, nv)   # nm^2

    atom_map = build_atom_residue_map(pdb_path)
    print(f"Surface: {nv} verts, {len(faces)} faces | PDB: {len(atom_map)} atoms")

    unique_patches = sorted(set(patch_arr[patch_arr >= 0]))
    print(f"Patches: {len(unique_patches)}")

    rows = []
    for pid in unique_patches:
        mask = patch_arr == pid
        ptype = "positive" if (values_arr is not None and values_arr[mask].mean() > 0) else "negative"

        if pep_compat:
            # 匹配 pep_patch 报告 (area[vertices].sum())
            vert_idx = np.flatnonzero(mask)
            valid = vert_idx[vert_idx < len(tri_areas)]
            total_area_A2 = tri_areas[valid].sum() * 100
            per_vert_A2 = vert_areas[mask] * 100  # 残基分配仍用顶点面积
        else:
            # 顶点面积法 (总面积 = 各顶点面积之和, 残基占比自洽)
            per_vert_A2 = vert_areas[mask] * 100
            total_area_A2 = per_vert_A2.sum()

        centroid = vertices[mask].mean(axis=0)
        mean_dist_A = np.mean(np.linalg.norm(vertices[mask] - centroid, axis=1)) * 10

        res_areas = defaultdict(float)
        res_verts = defaultdict(int)
        res_info = {}

        for idx, i in enumerate(np.flatnonzero(mask)):
            ai = atom_arr[i]
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


def write_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["patch_nr","patch_type","patch_total_area_A2","chain_id","res_name",
                     "res_seq","res_id","seq_nr","seq_res_id","n_vertices","area_A2",
                     "frac_of_patch","mean_dist_A"])
        for r in rows:
            w.writerow([r[k] for k in ["patch_nr","patch_type","patch_total_area_A2",
                         "chain_id","res_name","res_seq","res_id","seq_nr","seq_res_id",
                         "n_vertices","area_A2","frac_of_patch","mean_dist_A"]])


def print_summary(rows):
    by_patch = defaultdict(list)
    for r in rows:
        by_patch[r["patch_nr"]].append(r)
    pos = [p for p in by_patch.values() if p[0]["patch_type"]=="positive"]
    neg = [p for p in by_patch.values() if p[0]["patch_type"]=="negative"]
    print(f"\n{'='*60}\nPatch Summary: {len(pos)} pos + {len(neg)} neg\n{'='*60}")
    print(f"{'nr':>4} {'type':>10} {'area(A2)':>10} {'res':>5} {'top_residue':>14}")
    print("-"*50)
    for pn in sorted(by_patch):
        rl = by_patch[pn]
        top = max(rl, key=lambda x: x["frac_of_patch"])
        nres = len(set((r["chain_id"],r["res_seq"]) for r in rl))
        print(f"{pn:>4} {rl[0]['patch_type']:>10} {rl[0]['patch_total_area_A2']:>10.2f} {nres:>5} {top['chain_id']}:{top['res_name']}{top['res_seq']:>14}")


def main():
    p = argparse.ArgumentParser(description="pep-patch npz -> detailed CSV")
    p.add_argument("npz", nargs="?", help="npz path")
    p.add_argument("--npz", dest="npz2", help="npz path")
    p.add_argument("--pdb", help="PDB path (auto-detect if omitted)")
    p.add_argument("-o", "--output", help="output CSV path")
    p.add_argument("--pep-compat", action="store_true", help="面积匹配 pep_patch 报告 (area[vertices].sum)")
    args = p.parse_args()
    npz_path = args.npz or args.npz2
    if not npz_path:
        print("Usage: python analyze_pep_patch.py output.npz [--pdb protein.pdb]")
        return
    if not os.path.exists(npz_path):
        print(f"Not found: {npz_path}"); return

    pdb_path = args.pdb
    if not pdb_path:
        base = os.path.splitext(npz_path)[0]
        for sfx in ["_out", "_fixed_out", ""]:
            cand = base.replace(sfx, "") + ".pdb"
            if os.path.exists(cand):
                pdb_path = cand; break
    if not pdb_path or not os.path.exists(pdb_path):
        print("PDB not found, use --pdb"); return

    print(f"NPZ: {npz_path}\nPDB: {pdb_path}")
    rows = analyze_npz(npz_path, pdb_path, pep_compat=args.pep_compat)
    print_summary(rows)
    out = args.output or f"{os.path.splitext(npz_path)[0]}_detailed.csv"
    write_csv(rows, out)
    print(f"\nOutput: {out}")

if __name__ == "__main__":
    main()
