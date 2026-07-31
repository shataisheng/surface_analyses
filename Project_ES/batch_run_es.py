"""Batch run electrostatic analysis on all BsAb PDBs at given pH."""
import sys, os, csv, time, io, shutil, traceback
from datetime import datetime
from contextlib import redirect_stdout, redirect_stderr

# Allow running from this project folder while importing the upstream library.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from surface_analyses.platform_config import get_config
cfg = get_config()
cfg.setup_path()

from surface_analyses.commandline_electrostatic import run_electrostatics
from surface_analyses.structure import load_trajectory_using_commandline_args

# ── Config ──
PH = float(sys.argv[1]) if len(sys.argv) > 1 else 5.5
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else f"results/pH{str(PH).replace('.','_')}"
PDB_DIR = "."
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load all PDBs ──
pdbs = sorted([f for f in os.listdir(PDB_DIR) if f.endswith(".pdb")])
# Filter out incomplete ones
skip = ["Davutamig", "Ivicentamab", "Opugotamig", "Tovecimig", "Zanidatamab"]
pdbs = [p for p in pdbs if not any(p.startswith(s) for s in skip)]
print(f"Processing {len(pdbs)} PDBs at pH={PH} → {OUT_DIR}")
print(f"Skipped incomplete: {skip}")
print()

# ── Process ──
failed = []
for i, pdb_file in enumerate(pdbs):
    pdb_path = os.path.join(PDB_DIR, pdb_file)
    stem = pdb_file.replace("_fixed.pdb", "")
    es_prefix = f"{stem}_es"
    
    patches_csv = os.path.join(OUT_DIR, f"{es_prefix}_patches.csv")
    res_detailed = os.path.join(OUT_DIR, f"{es_prefix}_residues_detailed.csv")
    if os.path.exists(patches_csv) and os.path.exists(res_detailed):
        print(f"[{i+1}/{len(pdbs)}] SKIP {stem} (already done)")
        continue
    
    print(f"[{i+1}/{len(pdbs)}] {stem} ...", end=" ", flush=True)
    t0 = time.time()
    
    # Build dummy args
    class A: pass
    a = A()
    a.parm = pdb_path; a.trajs = [pdb_path]; a.stride = 1
    a.ref = None; a.protein_ref = None; a.dx = None
    a.apbs_dir = os.path.join(_REPO_ROOT, f"Tools/apbs_{stem}")
    a.probe_radius = 0.14
    a.out = patches_csv
    a.resout = None
    a.patch_cutoff = (2.0, -2.0)
    a.integral_cutoff = (0.3, -0.3)
    a.surface_type = "sas"
    a.ply_out = os.path.join(OUT_DIR, es_prefix)
    a.pos_patch_cmap = "tab20c"; a.neg_patch_cmap = "tab20c"
    a.ply_cmap = "coolwarm_r"; a.ply_clim = None
    a.check_cdrs = False
    a.n_patches = 0; a.size_cutoff = 0.0
    a.gauss_shift = 0.1; a.gauss_scale = 1.0
    a.pH = PH
    a.ion_species = None
    
    log_path = os.path.join(OUT_DIR, f"{es_prefix}_run.log")
    
    try:
        traj = load_trajectory_using_commandline_args(a)
        del a.parm, a.trajs, a.stride, a.ref, a.protein_ref
        
        with open(log_path, "w", encoding="utf-8") as log_f:
            with redirect_stdout(log_f), redirect_stderr(log_f):
                run_electrostatics(traj, **vars(a))
        
        elapsed = time.time() - t0
        print(f"OK ({elapsed:.0f}s)")
        
        # Generate residue-level CSV + patch summary via unified_analyzer
        try:
            from src.unified_analyzer import analyze
            result_df, _ = analyze(
                pdb_path=pdb_path,
                patches_csv=patches_csv,
                ply_pos=os.path.join(OUT_DIR, f"{es_prefix}-pos.ply"),
                ply_neg=os.path.join(OUT_DIR, f"{es_prefix}-neg.ply"),
                stem=es_prefix,
            )
            res_csv = os.path.join(OUT_DIR, f"{es_prefix}_residues_detailed.csv")
            result_df.to_csv(res_csv, index=False, encoding="utf-8-sig")
            print(f"  -> {es_prefix}_residues_detailed.csv")
            
            if "patch_total_area_A2" in result_df.columns:
                patch_summary = (
                    result_df.groupby(["patch_nr", "patch_type"])
                    .agg(
                        patch_total_area_A2=("patch_total_area_A2", "first"),
                        n_residues=("res_id", "nunique"),
                        top_residue=("res_id", "first"),
                        top_frac=("frac_of_patch", "first"),
                    )
                    .reset_index()
                    .sort_values(["patch_type", "patch_nr"])
                )
                summary_file = os.path.join(OUT_DIR, f"{es_prefix}_patch_summary.csv")
                patch_summary.to_csv(summary_file, index=False, encoding="utf-8-sig")
                print(f"  -> {es_prefix}_patch_summary.csv")
        except Exception as e:
            print(f"  (unified_analyzer warning: {e})")
        
        # Clean APBS temp files to save space
        apbs_dir = os.path.join(_REPO_ROOT, f"Tools/apbs_{stem}")
        if os.path.exists(apbs_dir):
            shutil.rmtree(apbs_dir, ignore_errors=True)
            
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
print("\nGenerating batch summary...")
from src.batch_summary import summarize_electrostatic, write_csv

rows = []
for pdb_file in pdbs:
    stem = pdb_file.replace("_fixed.pdb", "")
    s = summarize_electrostatic(stem, OUT_DIR)
    if s:
        rows.append(s)

if rows:
    path = os.path.join(OUT_DIR, "batch_summary_es.csv")
    write_csv(rows, path)
    print(f"Saved: {path} ({len(rows)} proteins)")
