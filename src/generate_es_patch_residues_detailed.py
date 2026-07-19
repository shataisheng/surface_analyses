#!/usr/bin/env python3
"""
generate_patch_residues_detailed.py  (改进版)
==============================================
从 pep_patch_electrostatic 的 PLY 输出提取 patch-per-residue 详细面积。

用法:
  python generate_patch_residues_detailed.py protein.pdb
  python generate_patch_residues_detailed.py protein.pdb -o detailed.csv

输入 (自动查找 PDB 同目录):
  - <stem>-pos.ply       正电荷 patch 表面
  - <stem>-neg.ply       负电荷 patch 表面
  - <stem>_patches.csv   patch 摘要

输出: <stem>_patch_residues_detailed.csv (UTF-8 BOM)

依赖: numpy, pandas, plyfile, mdtraj, scipy, biopython
"""

import numpy as np
import plyfile
import pandas as pd
import mdtraj as md
import argparse, os, glob, sys
from scipy.spatial import KDTree
from pathlib import Path
from Bio.PDB import PDBParser


def find_related_files(pdb_path):
    pdb_path = os.path.abspath(pdb_path)
    parent = os.path.dirname(pdb_path)
    stem = Path(pdb_path).stem
    for sfx in ["_fixed", "_out", "_output"]:
        if stem.endswith(sfx):
            stem = stem[:-len(sfx)]; break

    ply_pos = glob.glob(os.path.join(parent, f"{stem}*-pos.ply"))
    ply_pos = (ply_pos or glob.glob(os.path.join(parent, f"*pos*.ply")))[0] if ply_pos or glob.glob(os.path.join(parent, f"*pos*.ply")) else None
    ply_neg = glob.glob(os.path.join(parent, f"{stem}*-neg.ply"))
    ply_neg = (ply_neg or glob.glob(os.path.join(parent, f"*neg*.ply")))[0] if ply_neg or glob.glob(os.path.join(parent, f"*neg*.ply")) else None
    patches_csv = glob.glob(os.path.join(parent, f"{stem}*patches.csv"))
    patches_csv = patches_csv[0] if patches_csv else None
    return ply_pos, ply_neg, patches_csv


def load_pdb(pdb_file):
    print(f"Loading PDB: {pdb_file}")
    traj = md.load(pdb_file)
    top = traj.topology

    # 用 Bio.PDB 读取含插入码的正确残基标识
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    bio_struct = parser.get_structure("p", pdb_file)
    bio_residues = []
    for chain in sorted(bio_struct[0], key=lambda c: c.id):
        for res in chain:
            if res.id[0] != " ": continue
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
            chain_id = chr(ord("A")+res.chain.index)
            res_name = res.name
            res_seq_str = res.resSeq
        key = (res.chain.index, res.resSeq)
        if key not in seq_nr_map:
            seq_nr_map[key] = len(seq_nr_map) + 1
        seq_nr = seq_nr_map[key]
        records.append({
            "atom_idx": atom.index, "atom_name": atom.name,
            "res_name": res_name, "res_seq": res_seq_str,
            "res_id": f"{res_name}{res_seq_str}",
            "chain_id": chain_id,
            "seq_nr": seq_nr,
            "seq_res_id": f"{res_name}{seq_nr}",
        })
    atom_df = pd.DataFrame(records)
    atom_coords = traj.xyz[0]*10.0
    tree = KDTree(atom_coords)
    print(f"  {len(atom_df)} atoms, {top.n_residues} residues, {top.n_chains} chains")
    return atom_df, atom_coords, tree

def match_colors_to_patches(ply_colors, ply_coords_A, type_patches, atom_df, atom_coords):
    npoints_list = type_patches.npoints.tolist()
    nr_list = type_patches.nr.tolist()
    unique_colors, counts = np.unique(ply_colors, axis=0, return_counts=True)
    patch_colors = [(tuple(c), int(cnt)) for c, cnt in zip(unique_colors, counts) if not np.allclose(c, [256,256,256])]
    color_to_nr, unresolved = {}, []

    for color, n_verts in patch_colors:
        candidates = [i for i, np_ in enumerate(npoints_list) if np_ == n_verts]
        if not candidates:
            candidates = [min(range(len(npoints_list)), key=lambda i: abs(npoints_list[i]-n_verts))]
        if len(candidates) == 1:
            color_to_nr[color] = nr_list[candidates[0]]
        else:
            unresolved.append((color, n_verts, candidates))

    if unresolved:
        res_coord_map = {}
        for seq_res_id, grp in atom_df.groupby("seq_res_id"):
            res_coord_map[seq_res_id] = atom_coords[grp["atom_idx"].values].mean(axis=0)
        for color, n_verts, candidates in unresolved:
            mask = np.all(ply_colors == np.array(color), axis=1)
            centroid = ply_coords_A[mask].mean(axis=0)
            best_idx, best_dist = None, float("inf")
            for ci in candidates:
                patch_row = type_patches[type_patches.nr == nr_list[ci]].iloc[0]
                main_res_id = patch_row["main_residue"]
                matches = atom_df[atom_df["res_id"] == main_res_id]
                if len(matches) == 0:
                    continue
                res_centroid = atom_coords[matches["atom_idx"].values].mean(axis=0)
                dist = np.linalg.norm(centroid - res_centroid)
                if dist < best_dist:
                    best_dist, best_idx = dist, ci
            if best_idx is not None:
                color_to_nr[color] = nr_list[best_idx]
            else:
                color_to_nr[color] = nr_list[candidates[0]]
    return color_to_nr


def extract_patch_residues(ply_path, patch_type, patches_df, atom_df, atom_coords, tree, ply_scale=100.0):
    print(f"\nProcessing {patch_type}: {ply_path}")
    ply = plyfile.PlyData.read(ply_path)
    verts, faces_data = ply["vertex"], ply["face"]
    coords = np.stack([verts["x"], verts["y"], verts["z"]], axis=1)*ply_scale
    colors = np.stack([verts["red"], verts["green"], verts["blue"]], axis=1)

    face_arr = np.vstack([faces_data["vertex_indices"][i] for i in range(len(faces_data))])
    e1, e2 = coords[face_arr[:,1]]-coords[face_arr[:,0]], coords[face_arr[:,2]]-coords[face_arr[:,0]]
    face_areas = 0.5*np.linalg.norm(np.cross(e1,e2), axis=1)
    vertex_areas = np.zeros(len(coords))
    for i in range(3):
        np.add.at(vertex_areas, face_arr[:,i], face_areas/3.0)

    type_patches = patches_df[patches_df.type==patch_type].sort_values("nr").reset_index(drop=True)
    if len(type_patches) == 0:
        return pd.DataFrame()

    color_to_nr = match_colors_to_patches(colors, coords, type_patches, atom_df, atom_coords)
    print(f"  {len(color_to_nr)} color groups matched")

    results = []
    for color, patch_nr in color_to_nr.items():
        mask = np.all(colors==np.array(color), axis=1)
        patch_coords, patch_vareas = coords[mask], vertex_areas[mask]
        patch_row = type_patches[type_patches.nr==patch_nr].iloc[0]
        patch_area = float(patch_row["area"])*100.0

        dists, atom_idxs = tree.query(patch_coords, k=1)
        res_data = atom_df.iloc[atom_idxs].copy().reset_index(drop=True)
        res_data["vertex_area_A2"] = patch_vareas
        res_data["dist_A"] = dists

        res_summary = (res_data
            .groupby(["chain_id","res_name","res_seq","res_id","seq_nr","seq_res_id"])
            .agg(n_vertices=("atom_idx","count"), area_A2=("vertex_area_A2","sum"), mean_dist_A=("dist_A","mean"))
            .reset_index().sort_values("area_A2", ascending=False))
        res_summary["patch_nr"] = patch_nr
        res_summary["patch_type"] = patch_type
        res_summary["patch_total_area_A2"] = round(patch_area, 2)
        res_summary["frac_of_patch"] = (res_summary["area_A2"]/res_summary["area_A2"].sum()).round(4)
        results.append(res_summary)

    if not results:
        return pd.DataFrame()
    out = pd.concat(results, ignore_index=True)
    print(f"  {out.patch_nr.nunique()} patches, {len(out)} residues")
    return out


def main():
    p = argparse.ArgumentParser(description="PLY -> detailed CSV")
    p.add_argument("pdb", help="PDB path")
    p.add_argument("-o", "--output", help="output CSV")
    p.add_argument("--ply-pos", help="pos PLY")
    p.add_argument("--ply-neg", help="neg PLY")
    p.add_argument("--patches", help="patches CSV")
    p.add_argument("--scale", type=float, default=100.0, help="PLY scale")
    args = p.parse_args()

    if not os.path.exists(args.pdb):
        print(f"PDB not found: {args.pdb}"); return

    ply_pos, ply_neg, patches_csv = find_related_files(args.pdb)
    ply_pos = args.ply_pos or ply_pos
    ply_neg = args.ply_neg or ply_neg
    patches_csv = args.patches or patches_csv

    for label, f in [("pos PLY", ply_pos), ("neg PLY", ply_neg), ("patches CSV", patches_csv)]:
        if not f or not os.path.exists(f):
            print(f"Missing: {label} ({f})"); return

    print(f"PDB: {args.pdb}\nPos: {ply_pos}\nNeg: {ply_neg}\nCSV: {patches_csv}")

    atom_df, atom_coords, tree = load_pdb(args.pdb)
    patches_df = pd.read_csv(patches_csv)
    print(f"Patches: {len(patches_df)}")

    pos_res = extract_patch_residues(ply_pos, "positive", patches_df, atom_df, atom_coords, tree, args.scale)
    neg_res = extract_patch_residues(ply_neg, "negative", patches_df, atom_df, atom_coords, tree, args.scale)

    all_res = pd.concat([pos_res, neg_res], ignore_index=True)
    cols = ["patch_nr","patch_type","patch_total_area_A2","chain_id","res_name","res_seq","res_id","seq_nr","seq_res_id","n_vertices","area_A2","frac_of_patch","mean_dist_A"]
    cols = [c for c in cols if c in all_res.columns]
    all_res = all_res[cols].sort_values(["patch_type","patch_nr","area_A2"], ascending=[True,True,False])

    out = args.output or os.path.join(os.path.dirname(os.path.abspath(args.pdb)), f"{Path(args.pdb).stem}_patch_residues_detailed.csv")
    all_res.to_csv(out, index=False, encoding="utf-8-sig")
    print(f"\nSaved: {out}  ({len(all_res)} rows, {all_res.patch_nr.nunique()} patches)")

    for ptype in ["positive","negative"]:
        top = all_res[all_res.patch_type==ptype].sort_values(["patch_nr","area_A2"], ascending=[True,False])
        print(f"\n=== Top 3 {ptype} ===")
        for nr, grp in list(top.groupby("patch_nr"))[:3]:
            total = grp["patch_total_area_A2"].iloc[0]
            res_str = ", ".join(f"{r['res_id']}({r['chain_id']}) {r['frac_of_patch']*100:.0f}%" for _,r in grp.iterrows())
            print(f"  Patch {nr} ({total:.1f} A2): {res_str}")

if __name__ == "__main__":
    main()
