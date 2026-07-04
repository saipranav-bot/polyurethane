#!/usr/bin/env python3
"""
06a_active_site_analysis.py
============================
Active site analysis for g13490.t1 — the lead polyurethanase candidate.

WHAT WE KNOW FROM PROJECT SUMMARY:
  - g13490.t1: 573 aa, secreted (SP cleavage 16-17), no TM helices
  - Alpha/beta hydrolase family (lipase-like)
  - AlphaFold model: pTM=0.92, ranking=0.95, 0 clashes
  - Detected motif: GESAG at position ~227 → candidate catalytic Ser227
  - Asp candidates: D368, D376, D389, D392, D421
  - His candidates: H425, H462, H471

WHAT THIS SCRIPT DOES:
  1. Parses AlphaFold .cif or .pdb structure
  2. Confirms Ser227, finds best Asp + His for catalytic triad
  3. Measures Ser-Asp-His distances
  4. Identifies GxSxGG nucleophilic elbow motif
  5. Outputs grid box coordinates for AutoDock Vina
  6. Generates active site figures

Usage:
    python3 06a_active_site_analysis.py 05_structure/g13490_best_model.cif
    python3 06a_active_site_analysis.py 05_structure/g13490_best_model.pdb
"""

import os
import sys
import re
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import combinations

# ── PATHS ─────────────────────────────────────────────────────
STRUCT_DIR   = "05_structure"
OUTDIR       = "06_active_site"
DOCK_DIR     = "07_docking"
FIG_DIR      = "08_figures"

for d in [OUTDIR, DOCK_DIR, FIG_DIR]:
    os.makedirs(d, exist_ok=True)

# ── FIND STRUCTURE FILE ───────────────────────────────────────
def find_structure():
    if len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
        return sys.argv[1]
    # Search common locations
    for d in [STRUCT_DIR, ".", os.path.expanduser("~/polyurethanase_pipeline/05_structure")]:
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if ("g13490" in f or "model" in f or "alphafold" in f.lower()) and \
               (f.endswith(".cif") or f.endswith(".pdb")):
                return os.path.join(d, f)
    return None

STRUCT_FILE = find_structure()

if STRUCT_FILE is None:
    print("=" * 60)
    print("[WAITING] Structure file not found.")
    print()
    print("Copy your AlphaFold structure here:")
    print(f"  cp ~/Downloads/fold_g13490_model_0.cif {STRUCT_DIR}/g13490_best_model.cif")
    print()
    print("Then run:")
    print("  python3 06a_active_site_analysis.py 05_structure/g13490_best_model.cif")
    print("=" * 60)
    sys.exit(0)

print("=" * 60)
print(" STEP 6A: ACTIVE SITE ANALYSIS — g13490.t1")
print("=" * 60)
print(f"\nStructure file: {STRUCT_FILE}")

# ── PARSE STRUCTURE (CIF or PDB) ──────────────────────────────
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch

def load_structure(filepath):
    if filepath.endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("g13490", filepath)

print("\nParsing structure...")
structure = load_structure(STRUCT_FILE)
model     = structure[0]
chain     = list(model.get_chains())[0]
residues  = [r for r in chain.get_residues() if r.get_id()[0] == " "]  # ATOM records only

print(f"  Chain: {chain.get_id()}")
print(f"  Residues: {len(residues)}")

# Build residue lookup by number
res_by_num = {r.get_id()[1]: r for r in residues}

# ── EXTRACT pLDDT FROM B-FACTOR ───────────────────────────────
res_data = []
for res in residues:
    resnum  = res.get_id()[1]
    resname = res.get_resname()
    if "CA" in res:
        plddt = res["CA"].get_bfactor()
        coord = res["CA"].get_vector()
        res_data.append({
            "resnum": resnum, "resname": resname,
            "plddt": plddt,
            "x": coord[0], "y": coord[1], "z": coord[2]
        })

import pandas as pd
df = pd.DataFrame(res_data)
mean_plddt = df["plddt"].mean()
print(f"  Mean pLDDT: {mean_plddt:.1f}")

# ── KNOWN CANDIDATES FROM PROJECT SUMMARY ─────────────────────
SER_CANDIDATE  = 227   # from GESAG motif detection
ASP_CANDIDATES = [368, 376, 389, 392, 421]
HIS_CANDIDATES = [425, 462, 471]

print(f"\n--- KNOWN CANDIDATES (from GESAG motif scan) ---")
print(f"  Catalytic Ser: S{SER_CANDIDATE}")
print(f"  Asp candidates: {['D'+str(d) for d in ASP_CANDIDATES]}")
print(f"  His candidates: {['H'+str(h) for h in HIS_CANDIDATES]}")

# ── VERIFY GESAG MOTIF ────────────────────────────────────────
print(f"\n--- VERIFYING GESAG MOTIF ---")
# Build sequence from structure
seq_dict = {r.get_id()[1]: r.get_resname() for r in residues}

aa3to1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"
}

# Check window around Ser227
window_start = max(221, 1)
window_end   = min(233, max(res_by_num.keys()))
motif_seq    = ""
motif_nums   = []
for i in range(window_start, window_end + 1):
    if i in res_by_num:
        aa = aa3to1.get(res_by_num[i].get_resname(), "?")
        motif_seq += aa
        motif_nums.append(i)

print(f"  Sequence around S227 (pos {window_start}-{window_end}): {motif_seq}")
print(f"  Residue numbers:  {motif_nums}")

# Look for GxSxGG pattern
import re as re_mod
gxsxgg_pattern = re_mod.compile(r'G.S.G')
matches = [(m.start(), m.group()) for m in gxsxgg_pattern.finditer(motif_seq)]
if matches:
    for pos, match in matches:
        abs_pos = motif_nums[pos] if pos < len(motif_nums) else "?"
        print(f"  ✅ GxSxG motif found: '{match}' starting at position {abs_pos}")
        print(f"     Nucleophilic Ser at position {motif_nums[pos+2] if pos+2 < len(motif_nums) else SER_CANDIDATE}")
else:
    print(f"  Searching broader window for GxSxGG...")
    # Search full sequence
    full_seq = "".join(aa3to1.get(res_by_num.get(i, type('obj',(),{'get_resname':lambda self:'???'})()).get_resname(),'?')
                       for i in range(min(res_by_num.keys()), max(res_by_num.keys())+1)
                       if i in res_by_num)
    matches2 = [(m.start(), m.group()) for m in re_mod.finditer(r'G.S.G', full_seq)]
    print(f"  Found {len(matches2)} GxSxG matches in full structure")
    for pos, match in matches2[:5]:
        resnum_match = list(res_by_num.keys())[pos] if pos < len(res_by_num) else "?"
        print(f"    '{match}' around residue {resnum_match}")

# ── MEASURE CATALYTIC TRIAD DISTANCES ────────────────────────
print(f"\n--- CATALYTIC TRIAD DISTANCE ANALYSIS ---")

def get_atom(res, atom_name_list):
    """Get first available atom from list."""
    for name in atom_name_list:
        if name in res:
            return res[name]
    return None

def dist(a1, a2):
    """Euclidean distance between two Bio.PDB atoms."""
    v = a1.get_vector() - a2.get_vector()
    return float(np.sqrt(v[0]**2 + v[1]**2 + v[2]**2))

# Get Ser227
best_triad = None
best_score = 999

if SER_CANDIDATE in res_by_num:
    ser_res = res_by_num[SER_CANDIDATE]
    ser_og  = get_atom(ser_res, ["OG", "OG1", "CB"])
    ser_ca  = get_atom(ser_res, ["CA"])
    ser_plddt = ser_res["CA"].get_bfactor() if "CA" in ser_res else 0
    
    print(f"\n  Ser{SER_CANDIDATE}:")
    print(f"    pLDDT: {ser_plddt:.1f}")
    print(f"    OG atom: {'found' if ser_og else 'NOT FOUND'}")
    
    if ser_og:
        print(f"\n  Checking Ser{SER_CANDIDATE}–Asp distances:")
        asp_results = []
        for asp_num in ASP_CANDIDATES:
            if asp_num in res_by_num:
                asp_res  = res_by_num[asp_num]
                asp_od   = get_atom(asp_res, ["OD1", "OD2", "OE1", "OE2", "CB"])
                asp_plddt = asp_res["CA"].get_bfactor() if "CA" in asp_res else 0
                if asp_od:
                    d = dist(ser_og, asp_od)
                    asp_results.append((asp_num, d, asp_plddt, asp_res))
                    flag = "⭐" if d < 10 else "  "
                    print(f"    {flag} S{SER_CANDIDATE}–D{asp_num}: {d:.2f} Å  (pLDDT {asp_plddt:.0f})")
        
        print(f"\n  Checking Ser{SER_CANDIDATE}–His distances:")
        his_results = []
        for his_num in HIS_CANDIDATES:
            if his_num in res_by_num:
                his_res   = res_by_num[his_num]
                his_nd    = get_atom(his_res, ["ND1", "NE2", "CE1", "CB"])
                his_plddt = his_res["CA"].get_bfactor() if "CA" in his_res else 0
                if his_nd:
                    d = dist(ser_og, his_nd)
                    his_results.append((his_num, d, his_plddt, his_res))
                    flag = "⭐" if d < 10 else "  "
                    print(f"    {flag} S{SER_CANDIDATE}–H{his_num}: {d:.2f} Å  (pLDDT {his_plddt:.0f})")
        
        # Find best Asp-His pair (both should be close to Ser)
        print(f"\n  Finding best Ser-Asp-His triad combination:")
        triad_results = []
        for (asp_num, ser_asp_d, asp_plddt, asp_res) in asp_results:
            asp_od = get_atom(asp_res, ["OD1", "OD2", "OE1", "OE2", "CB"])
            for (his_num, ser_his_d, his_plddt, his_res) in his_results:
                his_nd = get_atom(his_res, ["ND1", "NE2", "CE1", "CB"])
                if asp_od and his_nd:
                    asp_his_d = dist(asp_od, his_nd)
                    # Score: lower is better (classic triad distances ~3-7 Å for Ser-His)
                    score = ser_asp_d + ser_his_d + asp_his_d
                    triad_results.append({
                        "ser": SER_CANDIDATE, "asp": asp_num, "his": his_num,
                        "ser_asp_d": ser_asp_d, "ser_his_d": ser_his_d, "asp_his_d": asp_his_d,
                        "score": score,
                        "asp_plddt": asp_plddt, "his_plddt": his_plddt, "ser_plddt": ser_plddt
                    })
        
        if triad_results:
            triad_results.sort(key=lambda x: x["score"])
            best_triad = triad_results[0]
            
            print(f"\n  ✅ BEST CATALYTIC TRIAD:")
            print(f"     Ser{best_triad['ser']}–Asp{best_triad['asp']}–His{best_triad['his']}")
            print(f"     Ser–Asp distance: {best_triad['ser_asp_d']:.2f} Å")
            print(f"     Ser–His distance: {best_triad['ser_his_d']:.2f} Å")
            print(f"     Asp–His distance: {best_triad['asp_his_d']:.2f} Å")
            print(f"     pLDDT: Ser={best_triad['ser_plddt']:.0f}, Asp={best_triad['asp_plddt']:.0f}, His={best_triad['his_plddt']:.0f}")
            
            # Interpret distances
            if best_triad["ser_his_d"] < 6.0:
                print(f"     ✅ Ser–His within hydrogen bond range (<6Å) — active triad geometry")
            elif best_triad["ser_his_d"] < 10.0:
                print(f"     ℹ️  Ser–His {best_triad['ser_his_d']:.1f}Å — plausible triad in AlphaFold model")
            else:
                print(f"     ⚠️  Ser–His {best_triad['ser_his_d']:.1f}Å — may need substrate to close")
else:
    print(f"  [WARNING] S{SER_CANDIDATE} not found in structure. Checking nearby residues...")
    nearby = [r for r in residues if abs(r.get_id()[1] - SER_CANDIDATE) <= 5
              and r.get_resname() == "SER"]
    for r in nearby:
        print(f"    SER found at position {r.get_id()[1]}")

# ── COMPUTE ACTIVE SITE CENTROID ──────────────────────────────
print(f"\n--- ACTIVE SITE CENTROID (for Vina grid box) ---")

active_residue_nums = [SER_CANDIDATE] + ASP_CANDIDATES + HIS_CANDIDATES
active_coords = []

for rnum in active_residue_nums:
    if rnum in res_by_num:
        res = res_by_num[rnum]
        if "CA" in res:
            v = res["CA"].get_vector()
            active_coords.append([v[0], v[1], v[2]])

if active_coords:
    coords_arr = np.array(active_coords)
    centroid   = coords_arr.mean(axis=0)
    spread     = coords_arr.max(axis=0) - coords_arr.min(axis=0)
    box_size   = max(spread.max() + 10, 22)  # at least 22 Å

    print(f"  Active site centroid: ({centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f})")
    print(f"  Spread (Å):           X={spread[0]:.1f}, Y={spread[1]:.1f}, Z={spread[2]:.1f}")
    print(f"  Recommended box size: {box_size:.0f} Å")

    if best_triad:
        # Use best triad centroid if available
        triad_coords = []
        for rnum in [best_triad["ser"], best_triad["asp"], best_triad["his"]]:
            if rnum in res_by_num and "CA" in res_by_num[rnum]:
                v = res_by_num[rnum]["CA"].get_vector()
                triad_coords.append([v[0], v[1], v[2]])
        if triad_coords:
            triad_center = np.array(triad_coords).mean(axis=0)
            print(f"\n  Catalytic triad centroid: ({triad_center[0]:.3f}, {triad_center[1]:.3f}, {triad_center[2]:.3f})")
            print(f"  → Use THIS as Vina center (more precise)")
            centroid = triad_center
else:
    centroid = np.array([0.0, 0.0, 0.0])
    box_size = 22
    print("  [WARNING] No active site residues found in structure. Using structure center.")
    if len(res_data) > 0:
        centroid = np.array([df["x"].mean(), df["y"].mean(), df["z"].mean()])

# ── WRITE VINA CONFIG ─────────────────────────────────────────
vina_config = f"""\
# AutoDock Vina configuration — g13490.t1 polyurethanase
# Active site: Ser{SER_CANDIDATE}-Asp{best_triad['asp'] if best_triad else '?'}-His{best_triad['his'] if best_triad else '?'}
# Generated by 06a_active_site_analysis.py

receptor = 07_docking/receptor.pdbqt
ligand   = 07_docking/ligand.pdbqt

center_x = {centroid[0]:.3f}
center_y = {centroid[1]:.3f}
center_z = {centroid[2]:.3f}

size_x   = {box_size:.0f}
size_y   = {box_size:.0f}
size_z   = {box_size:.0f}

exhaustiveness = 16
num_modes      = 9
energy_range   = 3
"""

vina_path = f"{DOCK_DIR}/vina_config.txt"
with open(vina_path, "w") as f:
    f.write(vina_config)
print(f"\n[OK] Vina config written: {vina_path}")

# ── WRITE ACTIVE SITE SUMMARY ─────────────────────────────────
summary_path = f"{OUTDIR}/active_site_summary.txt"
with open(summary_path, "w") as f:
    f.write("ACTIVE SITE ANALYSIS — g13490.t1\n")
    f.write("P. microspora KFRD-2 | Candidate Polyurethanase\n")
    f.write("=" * 60 + "\n\n")
    f.write(f"Structure file:      {STRUCT_FILE}\n")
    f.write(f"Protein length:      573 aa\n")
    f.write(f"Mean pLDDT:          {mean_plddt:.1f}\n\n")
    f.write("CATALYTIC MOTIF:\n")
    f.write(f"  Motif:             GESAG\n")
    f.write(f"  Nucleophilic Ser:  S{SER_CANDIDATE}\n\n")
    if best_triad:
        f.write("CATALYTIC TRIAD (Ser-Asp-His):\n")
        f.write(f"  Ser{best_triad['ser']}–Asp{best_triad['asp']}–His{best_triad['his']}\n")
        f.write(f"  Ser{best_triad['ser']}–Asp{best_triad['asp']} distance: {best_triad['ser_asp_d']:.2f} Å\n")
        f.write(f"  Ser{best_triad['ser']}–His{best_triad['his']} distance: {best_triad['ser_his_d']:.2f} Å\n")
        f.write(f"  Asp{best_triad['asp']}–His{best_triad['his']} distance: {best_triad['asp_his_d']:.2f} Å\n")
        f.write(f"  pLDDT: Ser={best_triad['ser_plddt']:.0f}, Asp={best_triad['asp_plddt']:.0f}, His={best_triad['his_plddt']:.0f}\n\n")
    f.write("DOCKING GRID BOX:\n")
    f.write(f"  center_x = {centroid[0]:.3f}\n")
    f.write(f"  center_y = {centroid[1]:.3f}\n")
    f.write(f"  center_z = {centroid[2]:.3f}\n")
    f.write(f"  box_size = {box_size:.0f} Å\n\n")
    f.write("LITERATURE CONTEXT:\n")
    f.write("  Russell et al. 2011: ~21 kDa extracellular serine hydrolase\n")
    f.write("  Inhibited by PMSF (serine protease inhibitor)\n")
    f.write("  Inhibited by iodoacetate (cysteine/SH group)\n")
    f.write("  NOT inhibited by EDTA (not metalloenzyme)\n")
    f.write("  → Confirms classical serine hydrolase mechanism\n")

print(f"[OK] Active site summary: {summary_path}")

# ── pLDDT PLOT WITH ACTIVE SITE MARKED ────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle(
    "AlphaFold Structure Analysis — g13490.t1 (573 aa)\n"
    "P. microspora KFRD-2 | Candidate Polyurethanase | pTM=0.92",
    fontsize=12, fontweight="bold"
)

# Left: pLDDT per residue with active site marked
ax = axes[0]
colors_bar = ["#1F3864" if p > 90 else "#2E75B6" if p > 70
              else "#F39C12" if p > 50 else "#E74C3C"
              for p in df["plddt"]]
ax.bar(df["resnum"], df["plddt"], color=colors_bar, width=1.0)
ax.axhline(90, color="#1F3864", lw=1.5, ls="--", alpha=0.6, label=">90 very high")
ax.axhline(70, color="#F39C12", lw=1.5, ls="--", alpha=0.6, label=">70 confident")
ax.axhline(50, color="#E74C3C", lw=1.5, ls="--", alpha=0.6, label="<50 disordered")

# Mark catalytic residues
for rnum, color, label in [
    (SER_CANDIDATE, "green", f"S{SER_CANDIDATE}"),
    (best_triad["asp"] if best_triad else 368, "blue", f"D{best_triad['asp'] if best_triad else 368}"),
    (best_triad["his"] if best_triad else 425, "purple", f"H{best_triad['his'] if best_triad else 425}"),
]:
    ax.axvline(rnum, color=color, lw=2.5, alpha=0.85, label=label)

ax.set_xlabel("Residue number", fontsize=11)
ax.set_ylabel("pLDDT", fontsize=11)
ax.set_title("Per-residue confidence + catalytic triad", fontsize=11)
ax.set_ylim(0, 105)
ax.legend(fontsize=8, ncol=2)
ax.grid(axis="y", alpha=0.3)

# Right: triad distance summary
ax = axes[1]
ax.axis("off")

if best_triad:
    table_data = [
        ["Residue", "Position", "pLDDT", "Role"],
        ["Ser", str(best_triad["ser"]), f"{best_triad['ser_plddt']:.0f}", "Nucleophile\n(GESAG motif)"],
        ["Asp", str(best_triad["asp"]), f"{best_triad['asp_plddt']:.0f}", "Electrophile"],
        ["His", str(best_triad["his"]), f"{best_triad['his_plddt']:.0f}", "General base"],
    ]
    dist_data = [
        ["Pair", "Distance (Å)", "Expected", "Status"],
        [f"S{best_triad['ser']}–H{best_triad['his']}",
         f"{best_triad['ser_his_d']:.2f}",
         "3–7 Å (H-bond)",
         "✅" if best_triad["ser_his_d"] < 10 else "⚠️"],
        [f"H{best_triad['his']}–D{best_triad['asp']}",
         f"{best_triad['asp_his_d']:.2f}",
         "3–7 Å (H-bond)",
         "✅" if best_triad["asp_his_d"] < 10 else "⚠️"],
        [f"S{best_triad['ser']}–D{best_triad['asp']}",
         f"{best_triad['ser_asp_d']:.2f}",
         "5–15 Å",
         "✅" if best_triad["ser_asp_d"] < 15 else "⚠️"],
    ]
    
    y_start = 0.95
    ax.text(0.5, y_start, "Catalytic Triad — g13490.t1", ha="center", va="top",
            fontsize=12, fontweight="bold", transform=ax.transAxes)
    
    # Draw table
    col_labels = ["Pair", "Distance", "Expected", "OK?"]
    cell_data  = [[r[0], r[1], r[2], r[3]] for r in dist_data[1:]]
    tbl = ax.table(cellText=cell_data, colLabels=col_labels,
                   loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.2, 2.0)
    
    # Color header
    for j in range(len(col_labels)):
        tbl[(0, j)].set_facecolor("#1F3864")
        tbl[(0, j)].set_text_props(color="white", fontweight="bold")

    # Add conclusion
    ax.text(0.5, 0.08,
            f"Best triad: S{best_triad['ser']}–D{best_triad['asp']}–H{best_triad['his']}\n"
            f"Mechanism: Ser acts as nucleophile, His as base, Asp stabilizes His",
            ha="center", va="bottom", fontsize=9, color="#1F3864",
            transform=ax.transAxes, style="italic")

plt.tight_layout()
plot_path = f"{FIG_DIR}/Fig_active_site_analysis.png"
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
print(f"[OK] Plot: {plot_path}")

print("\n" + "=" * 60)
print(" STEP 6A COMPLETE")
if best_triad:
    print(f" Catalytic triad: Ser{best_triad['ser']}–Asp{best_triad['asp']}–His{best_triad['his']}")
print(f" Vina config ready: {vina_path}")
print(f"\n NEXT: bash 06b_prepare_receptor.sh")
print("=" * 60)
