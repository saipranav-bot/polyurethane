#!/usr/bin/env python3
"""
03_filter_by_size.py
====================
Filters BRAKER predicted proteins to the ~21 kDa range.
Input: 01_extract_proteins/proteins_cleaned.faa
Output: 02_size_filter/candidates_sizefiltered.faa
        → this file gets uploaded to SignalP 6.0

Scientific rationale:
    Russell 2011 SDS-PAGE: polyurethanase = ~21 kDa band
    21 kDa / 110 Da per aa = ~191 aa
    SDS-PAGE error ±15% → 18–24 kDa = 164–218 aa
    We use 135–230 aa (generous) to ensure we don't miss it
    This filter removes ~85-90% of the proteome, making SignalP practical
"""

import os
import sys
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ── CONFIG ────────────────────────────────────────────────────
INPUT_FASTA  = "01_extract_proteins/proteins_cleaned.faa"
OUTPUT_DIR   = "02_size_filter"
OUTPUT_FASTA = f"{OUTPUT_DIR}/candidates_sizefiltered.faa"
OUTPUT_STATS = f"{OUTPUT_DIR}/size_filter_stats.txt"
OUTPUT_PLOT  = "08_figures/01_size_distribution.png"
OUTPUT_INSTRUCTIONS = f"{OUTPUT_DIR}/NEXT_STEP_SIGNALP_INSTRUCTIONS.txt"

MIN_AA = 135   # ~15 kDa lower bound
MAX_AA = 230   # ~25 kDa upper bound
TARGET_AA = 191  # 21 kDa target

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs("08_figures", exist_ok=True)

# ── CHECK INPUT ───────────────────────────────────────────────
if not os.path.exists(INPUT_FASTA):
    print(f"[ERROR] {INPUT_FASTA} not found.")
    print("  Run 02_extract_proteins.py first")
    sys.exit(1)

print("=" * 55)
print(" STEP 3: FILTER BY SIZE (~21 kDa = 135–230 aa)")
print("=" * 55)

# ── LOAD AND FILTER ───────────────────────────────────────────
print(f"\nLoading: {INPUT_FASTA}")
all_records = list(SeqIO.parse(INPUT_FASTA, "fasta"))
all_lengths = [len(r.seq) for r in all_records]
print(f"  Total loaded: {len(all_records):,}")

filtered = [r for r in all_records if MIN_AA <= len(r.seq) <= MAX_AA]
filtered_lengths = [len(r.seq) for r in filtered]
reduction = 100 * (1 - len(filtered) / len(all_records))

print(f"\nSize filter: {MIN_AA}–{MAX_AA} aa ({MIN_AA*110/1000:.0f}–{MAX_AA*110/1000:.0f} kDa)")
print(f"  Before: {len(all_records):,} proteins")
print(f"  After:  {len(filtered):,} proteins")
print(f"  Reduction: {reduction:.1f}%")

# ── WRITE OUTPUT ──────────────────────────────────────────────
SeqIO.write(filtered, OUTPUT_FASTA, "fasta")
print(f"\n[OK] Written: {OUTPUT_FASTA}")

# ── SIZE STATS ────────────────────────────────────────────────
closest = sorted(filtered, key=lambda r: abs(len(r.seq) - TARGET_AA))
print(f"\nTop 5 sequences closest to 21 kDa ({TARGET_AA} aa):")
print(f"  {'Protein ID':<35} {'Length':>7} {'MW (kDa)':>10}")
for rec in closest[:5]:
    print(f"  {rec.id:<35} {len(rec.seq):>7} {len(rec.seq)*110/1000:>9.1f}")

# ── WRITE STATS ───────────────────────────────────────────────
with open(OUTPUT_STATS, "w") as f:
    f.write("SIZE FILTER STATISTICS\n")
    f.write("=" * 50 + "\n\n")
    f.write(f"Input:              {INPUT_FASTA}\n")
    f.write(f"Total proteins:     {len(all_records):,}\n")
    f.write(f"Size range:         {MIN_AA}–{MAX_AA} aa\n")
    f.write(f"After filter:       {len(filtered):,}\n\n")
    f.write("All size-filtered proteins:\n")
    for rec in sorted(filtered, key=lambda r: abs(len(r.seq) - TARGET_AA)):
        f.write(f"  {rec.id:<35} {len(rec.seq):>5} aa  ~{len(rec.seq)*110/1000:.1f} kDa\n")
print(f"[OK] Stats: {OUTPUT_STATS}")

# ── WRITE INSTRUCTIONS FOR SIGNALP ────────────────────────────
with open(OUTPUT_INSTRUCTIONS, "w") as f:
    f.write("=" * 60 + "\n")
    f.write("NEXT STEP: RUN SIGNALP 6.0 (WEB SERVER)\n")
    f.write("=" * 60 + "\n\n")
    f.write("File to upload:\n")
    f.write(f"  {os.path.abspath(OUTPUT_FASTA)}\n\n")
    f.write("Instructions:\n")
    f.write("  1. Go to: https://services.healthtech.dtu.dk/services/SignalP-6.0/\n")
    f.write("  2. Register with your email (free academic)\n")
    f.write("  3. Upload the file above\n")
    f.write("  4. Organism: Eukarya\n")
    f.write("  5. Output format: Long (includes all prediction columns)\n")
    f.write("  6. Mode: Slow (more accurate)\n")
    f.write("  7. Submit and wait (~10-30 min for this file size)\n")
    f.write("  8. Download result CSV\n")
    f.write("  9. Save as: 03_secretome/signalp6_output.csv\n")
    f.write(" 10. Then run: python3 04_filter_secretome.py\n\n")
    f.write("ALSO run TMHMM 2.0:\n")
    f.write("  URL: https://services.healthtech.dtu.dk/services/TMHMM-2.0/\n")
    f.write("  Upload same file, download result, save as: 03_secretome/tmhmm_output.txt\n\n")
    f.write(f"Proteins to submit: {len(filtered)}\n")
print(f"[OK] Instructions: {OUTPUT_INSTRUCTIONS}")

# ── PLOT ──────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle(
    "Protein Size Filter — Targeting ~21 kDa Polyurethanase\n"
    "P. microspora KFRD-2 BRAKER-Predicted Proteome",
    fontsize=12, fontweight='bold'
)

# Left: full proteome (clipped at 3000 aa for readability)
ax = axes[0]
plot_lengths = [l for l in all_lengths if l < 3000]
ax.hist(plot_lengths, bins=100, color='#2E75B6', alpha=0.7, edgecolor='white')
ax.axvspan(MIN_AA, MAX_AA, alpha=0.35, color='#E74C3C',
           label=f'Target zone ({MIN_AA}–{MAX_AA} aa)\nn={len(filtered):,}')
ax.axvline(TARGET_AA, color='#E74C3C', lw=2.5, ls='--', label=f'21 kDa ({TARGET_AA} aa)')
ax.set_xlabel("Protein length (aa)", fontsize=11)
ax.set_ylabel("Number of proteins", fontsize=11)
ax.set_title(f"Full Proteome (n={len(all_records):,}, shown <3000 aa)", fontsize=11)
ax.legend(fontsize=9)
ax.grid(axis='y', alpha=0.3)

# Right: target zone zoom + dual axis
ax = axes[1]
ax.hist(filtered_lengths, bins=30, color='#E74C3C', alpha=0.8, edgecolor='white')
ax.axvline(TARGET_AA, color='#1F3864', lw=2.5, ls='--', label=f'21 kDa ({TARGET_AA} aa)')
ax.set_xlabel("Protein length (aa)", fontsize=11)
ax.set_ylabel("Number of proteins", fontsize=11)
ax.set_title(f"Size-Filtered Candidates (n={len(filtered):,})", fontsize=11)
ax.legend(fontsize=9)
ax.grid(axis='y', alpha=0.3)

ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
kda_ticks = [135, 155, 175, 191, 210, 230]
ax2.set_xticks(kda_ticks)
ax2.set_xticklabels([f'{v*110/1000:.1f}' for v in kda_ticks], fontsize=8)
ax2.set_xlabel("Approximate MW (kDa)", fontsize=9, color='gray')

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"[OK] Plot: {OUTPUT_PLOT}")

print("\n" + "=" * 55)
print(f" STEP 3 COMPLETE")
print(f" Candidates for SignalP: {len(filtered):,}")
print(f"\n ⚠️  NOW DO THIS:")
print(f"   1. Upload {OUTPUT_FASTA}")
print(f"      to https://services.healthtech.dtu.dk/services/SignalP-6.0/")
print(f"   2. Also upload to TMHMM 2.0 on same website")
print(f"   3. See {OUTPUT_INSTRUCTIONS} for full instructions")
print(f"   4. Save results then run: python3 04_filter_secretome.py")
print("=" * 55)
