#!/usr/bin/env python3
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS_DIR = "07_docking/results_corrected"
FIG_DIR     = "08_figures"
OUTDIR      = "07_docking"
os.makedirs(FIG_DIR, exist_ok=True)

LIGANDS = {
    "L1_ethyl_carbamate":  {"label": "Ethyl carbamate\n(urethane unit)",   "color": "#1F3864"},
    "L2_diethyl_adipate":  {"label": "Diethyl adipate\n(ester bond)",      "color": "#2E75B6"},
    "L3_pur_fragment":     {"label": "PUR fragment\n(ester+urethane)",     "color": "#1A5276"},
    "L4_mdi_fragment":     {"label": "MDI fragment\n(aromatic PUR)",       "color": "#6C3483"},
    "L5_dmso_control":     {"label": "DMSO\n(negative ctrl)",              "color": "#7F8C8D"},
}

print("=" * 60)
print(" STEP 7C (CORRECTED): DOCKING RESULTS ANALYSIS")
print(f" Reading from: {RESULTS_DIR}")
print(" Triad: Ser227-His462-Glu226 (corrected)")
print("=" * 60)

if not os.path.isdir(RESULTS_DIR):
    print(f"[ERROR] {RESULTS_DIR} not found.")
    sys.exit(1)

all_results = {}
for ligand_name, meta in LIGANDS.items():
    log_path = os.path.join(RESULTS_DIR, f"{ligand_name}_vina_log.txt")
    if not os.path.exists(log_path):
        print(f"  [SKIP] {ligand_name}: log not found at {log_path}")
        continue
    poses = []
    with open(log_path) as f:
        for line in f:
            m = re.match(r"\s+(\d+)\s+([-\d\.]+)\s+([\d\.]+)\s+([\d\.]+)", line)
            if m:
                poses.append({
                    "pose": int(m.group(1)),
                    "affinity": float(m.group(2)),
                    "rmsd_lb": float(m.group(3)),
                    "rmsd_ub": float(m.group(4)),
                })
    if poses:
        all_results[ligand_name] = {
            **meta,
            "poses": pd.DataFrame(poses),
            "best_energy": min(p["affinity"] for p in poses),
            "n_poses": len(poses),
        }
        print(f"  {ligand_name:<25} best={all_results[ligand_name]['best_energy']:.3f} kcal/mol  ({len(poses)} poses)")

if not all_results:
    print("[ERROR] No docking logs parsed.")
    sys.exit(1)

names  = list(all_results.keys())
bests  = [all_results[n]["best_energy"] for n in names]
labels = [all_results[n]["label"] for n in names]
colors = [all_results[n]["color"] for n in names]

ctrl_energy = all_results.get("L5_dmso_control", {}).get("best_energy", 0)
best_substrate = min(
    [n for n in all_results if n != "L5_dmso_control"],
    key=lambda n: all_results[n]["best_energy"]
)
best_dg = all_results[best_substrate]["best_energy"]
ddg_best = best_dg - ctrl_energy

print(f"\n  Best substrate: {best_substrate} ({best_dg:.3f} kcal/mol)")
print(f"  DMSO control:   {ctrl_energy:.3f} kcal/mol")
print(f"  DDG (best vs DMSO): {ddg_best:.3f} kcal/mol")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
fig.suptitle(
    "Molecular Docking Results — g13490.t1 Candidate Polyurethanase\n"
    "P. microspora KFRD-2 | AutoDock Vina | Active site: Ser227-His462-Glu226",
    fontsize=13, fontweight="bold"
)

ax = axes[0]
x = range(len(names))
bars = ax.bar(x, bests, color=colors, edgecolor="white", width=0.6, zorder=3)
ax.axhline(-4.0, color="#E74C3C", lw=2, ls="--", label="Significant binding (-4 kcal/mol)", zorder=4)
ax.axhline(-6.0, color="#E67E22", lw=1.5, ls=":", label="Good binding (-6 kcal/mol)", zorder=4)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("Best binding affinity (kcal/mol)", fontsize=11)
ax.set_title("Best Docking Affinity per Ligand", fontsize=11)
ax.legend(fontsize=9)
ax.grid(axis="y", alpha=0.3, zorder=0)
ax.set_ylim(min(bests) - 1, 0)
for bar, e in zip(bars, bests):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() - 0.15,
            f"{e:.2f}", ha="center", va="top", fontsize=10, fontweight="bold", color="white")

ax = axes[1]
for name, res in all_results.items():
    poses = res["poses"]
    ax.scatter(poses["pose"], poses["affinity"], c=res["color"], s=60, alpha=0.8,
               label=res["label"].replace("\n", " "), zorder=3)
    ax.plot(poses["pose"], poses["affinity"], color=res["color"], alpha=0.35, lw=1.5)
ax.axhline(-4.0, color="#E74C3C", lw=2, ls="--", alpha=0.7)
ax.set_xlabel("Docking pose rank", fontsize=11)
ax.set_ylabel("Binding affinity (kcal/mol)", fontsize=11)
ax.set_title("Binding Affinity Across All Poses", fontsize=11)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
ax.set_ylim(min(bests) - 1, -1)

plt.tight_layout()
p1 = f"{FIG_DIR}/Fig_docking_energies.png"
plt.savefig(p1, dpi=300, bbox_inches="tight")
plt.close()
print(f"\n[OK] {p1}")

fig, ax = plt.subplots(figsize=(8, 6))
rel_names  = [n for n in names if n != "L5_dmso_control"]
delta_g    = [all_results[n]["best_energy"] - ctrl_energy for n in rel_names]
colors_rel = [all_results[n]["color"] for n in rel_names]
labels_rel = [all_results[n]["label"] for n in rel_names]

bars = ax.barh(range(len(rel_names)), delta_g, color=colors_rel, edgecolor="white", height=0.5)
ax.set_yticks(range(len(rel_names)))
ax.set_yticklabels(labels_rel, fontsize=10)
ax.set_xlabel("ΔG relative to DMSO control (kcal/mol)", fontsize=11)
ax.set_title(
    "Substrate Selectivity — ΔΔG vs Negative Control\n"
    "More negative = better substrate binding over non-substrate",
    fontsize=11, fontweight="bold"
)
ax.axvline(0, color="#566573", lw=1.5, ls="--")
ax.axvline(-2, color="#E74C3C", lw=1.5, ls=":", alpha=0.7, label="ΔΔG = -2 (selective)")
ax.legend(fontsize=9)
ax.grid(axis="x", alpha=0.3)
ax.invert_yaxis()
for bar, delta in zip(bars, delta_g):
    ax.text(bar.get_width() - 0.1, bar.get_y() + bar.get_height()/2,
            f"{delta:.2f}", va="center", ha="right", fontsize=10, fontweight="bold", color="white")

plt.tight_layout()
p2 = f"{FIG_DIR}/Fig_selectivity_vs_control.png"
plt.savefig(p2, dpi=300, bbox_inches="tight")
plt.close()
print(f"[OK] {p2}")

report_path = f"{OUTDIR}/FINAL_CANDIDATE_REPORT_CORRECTED.txt"
with open(report_path, "w") as f:
    f.write("FINAL CANDIDATE REPORT (CORRECTED)\n")
    f.write("Triad: Ser227-His462-Glu226 | Grid center: (-5.194, 0.640, -2.398)\n")
    f.write("=" * 60 + "\n\n")
    for name in sorted(names, key=lambda n: all_results[n]["best_energy"]):
        f.write(f"  {name:<25} {all_results[name]['best_energy']:.3f} kcal/mol\n")
    f.write(f"\nBest substrate: {best_substrate} ({best_dg:.3f} kcal/mol)\n")
    f.write(f"DMSO control:   {ctrl_energy:.3f} kcal/mol\n")
    f.write(f"DDG:            {ddg_best:.3f} kcal/mol\n")
print(f"[OK] {report_path}")

print("\n" + "=" * 60)
print(" CORRECTED ANALYSIS COMPLETE")
print(f" Best: {best_substrate} = {best_dg:.3f} kcal/mol")
print(f" DDG:  {ddg_best:.3f} kcal/mol")
print("=" * 60)
