"""Batch run hydrophobic (Crippen) analysis on all BsAb PDBs."""
import sys, os, csv, time, shutil, traceback, re, glob as globmod
from contextlib import redirect_stdout, redirect_stderr

# Allow running from this project folder while importing the upstream library.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from surface_analyses.platform_config import get_config
cfg = get_config()
cfg.setup_path()

from surface_analyses.commandline_hydrophobic import run_hydrophobic
from surface_analyses.structure import load_trajectory_using_commandline_args

# ── Config ──
OUT_DIR = "results/hb_crippen"
PDB_DIR = "."
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load all PDBs ──
pdbs = sorted([f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")])
skip = ["Davutamig", "Ivicentamab", "Opugotamig", "Tovecimig", "Zanidatamab"]
pdbs = [p for p in pdbs if not any(p.startswith(s) for s in skip)]
print(f"Processing {len(pdbs)} PDBs with Crippen scale → {OUT_DIR}")

# ── Process ──
failed = []
for i, pdb_file in enumerate(pdbs):
    pdb_path = os.path.join(PDB_DIR, pdb_file)
    stem = pdb_file.replace("_fixed.pdb", "")
    hb_prefix = f"{stem}_hb"
    
    out_npz = os.path.join(OUT_DIR, f"{hb_prefix}_out.npz")
    if os.path.exists(out_npz):
        print(f"[{i+1}/{len(pdbs)}] SKIP {stem} (already done)")
        continue
    
    print(f"[{i+1}/{len(pdbs)}] {stem} ...", end=" ", flush=True)
    t0 = time.time()
    
    class A: pass
    a = A()
    a.parm = pdb_path; a.trajs = [pdb_path]; a.stride = 1
    a.ref = None; a.protein_ref = None
    a.scale = "crippen"
    a.smiles = None; a.atom_propensities = None
    a.out = out_npz
    a.surftype = "normal"
    a.group_heavy = False
    a.surfscore = True
    a.sap = True
    a.blur_rad = 0.5
    a.sh = True
    a.sh_rad = 0.8
    a.potential = True
    a.rmax = 0.3; a.solv_rad = 0.14
    a.grid_spacing = 0.05; a.rcut = 0.5; a.alpha = 15.0; a.blur_sigma = 0.6
    a.ply_out = None
    a.ply_cmap = None; a.ply_clim = None
    a.patches = False; a.patch_min = 0.12
    a.verbose = False
    
    log_path = os.path.join(OUT_DIR, f"{hb_prefix}_run.log")
    
    try:
        traj = load_trajectory_using_commandline_args(a)
        del a.parm, a.trajs, a.stride, a.ref, a.protein_ref
        
        with open(log_path, "w", encoding="utf-8") as log_f:
            with redirect_stdout(log_f), redirect_stderr(log_f):
                run_hydrophobic(pdb_path, traj, **vars(a))
        
        elapsed = time.time() - t0
        print(f"OK ({elapsed:.0f}s)")
    except Exception as e:
        elapsed = time.time() - t0
        print(f"FAIL ({elapsed:.0f}s): {e}")
        failed.append((stem, str(e)))
        traceback.print_exc()

print(f"\nDone: {len(pdbs)-len(failed)} OK, {len(failed)} failed")
if failed:
    print("Failed:")
    for s, e in failed:
        print(f"  {s}: {e}")

# ── Generate batch summary ──
print("\nGenerating hydrophobic batch summary...")
import numpy as np

rows = []
for pdb_file in pdbs:
    stem = pdb_file.replace("_fixed.pdb", "")
    npz_path = os.path.join(OUT_DIR, f"{stem}_hb_out.npz")
    if not os.path.exists(npz_path):
        continue
    
    try:
        data = np.load(npz_path, allow_pickle=True)
        row = {"protein": stem}
        
        # surfscore (per-atom surface hydrophobicity score = SASA * propensity)
        if "surfscore" in data:
            sf = data["surfscore"]
            row["hb_surfscore_mean"] = round(float(np.mean(sf)), 4)
            row["hb_surfscore_sum"] = round(float(np.sum(sf)), 2)
        
        # SAP (spatial aggregation propensity)
        if "sap" in data:
            sap = data["sap"]
            row["hb_sap_mean"] = round(float(np.mean(sap)), 4)
            row["hb_sap_max"] = round(float(np.max(sap)), 4)
        
        # surrounding hydrophobicity
        sh_key = "surrounding_hydrophobicity"
        if sh_key in data:
            sh = data[sh_key]
            row["hb_sh_mean"] = round(float(np.mean(sh)), 4)
            row["hb_sh_max"] = round(float(np.max(sh)), 4)
        
        # potential (mean value on surface vertices)
        pot_key = "hydrophobic_potential:0:data:values"
        if pot_key in data:
            pot = data[pot_key]
            row["hb_potential_mean"] = round(float(np.mean(pot)), 4)
            row["hb_potential_std"] = round(float(np.std(pot)), 4)
        
        # propensities (per-atom Crippen logP values)
        if "propensities" in data:
            prop = data["propensities"]
            row["hb_propensity_mean"] = round(float(np.mean(prop)), 4)
        
        data.close()
        rows.append(row)
    except Exception as e:
        print(f"  ERROR parsing {stem}: {e}")

if rows:
    cols = list(rows[0].keys())
    summary_path = os.path.join(OUT_DIR, "batch_summary_hb_crippen.csv")
    with open(summary_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"Saved: {summary_path} ({len(rows)} proteins)")
    
    # Quick stats
    for col in cols[1:]:
        vals = [r[col] for r in rows if col in r]
        if vals:
            print(f"  {col}: mean={np.mean(vals):.4f}  median={np.median(vals):.4f}  min={np.min(vals):.4f}  max={np.max(vals):.4f}")
