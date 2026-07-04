#!/usr/bin/env python3
"""
02_extract_proteins.py
======================
Run this AS SOON AS BRAKER finishes.

Takes braker.aa (predicted protein sequences from AUGUSTUS gene models)
and performs quality control before downstream analysis.

Usage:
    python3 02_extract_proteins.py /path/to/braker_out/braker.aa

What BRAKER gives you in braker.aa:
    >g1.t1     <- gene 1, transcript 1
    MKVLSPADKT...
    >g2.t1
    MQLFSRT...
    Each header = one predicted protein from one gene model.
"""

import os
import sys
import re
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ── CONFIG ────────────────────────────────────────────────────
# Accept path as argument or use default
if len(sys.argv) > 1:
    BRAKER_AA = sys.argv[1]
else:
    # Common default locations — try all
    candidates = [
        os.path.expanduser("~/polyurethanase_pipeline/results/01_genome/braker_out/braker.aa"),
        os.path.expanduser("~/braker_out/braker.aa"),
        "braker_out/braker.aa",
        "braker.aa",
    ]
    BRAKER_AA = None
    for c in candidates:
        if os.path.exists(c):
            BRAKER_AA = c
            break
    if BRAKER_AA is None:
        print("[ERROR] braker.aa not found.")
        print("Usage: python3 02_extract_proteins.py /path/to/braker.aa")
        print("Tried:")
        for c in candidates:
            print(f"  {c}")
        sys.exit(1)

OUTPUT_DIR    = "01_extract_proteins"
OUTPUT_FASTA  = f"{OUTPUT_DIR}/all_predicted_proteins.faa"
OUTPUT_CLEAN  = f"{OUTPUT_DIR}/proteins_cleaned.faa"
OUTPUT_STATS  = f"{OUTPUT_DIR}/proteome_qc_stats.txt"
OUTPUT_PLOT   = "08_figures/00_proteome_overview.png"

os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs("08_figures", exist_ok=True)

print("=" * 55)
print(" STEP 2: EXTRACT + QC BRAKER PREDICTED PROTEINS")
print("=" * 55)
print(f"\nInput: {BRAKER_AA}")

# ── LOAD BRAKER OUTPUT ────────────────────────────────────────
print("\nLoading protein sequences from BRAKER output...")
raw_records = list(SeqIO.parse(BRAKER_AA, "fasta"))
print(f"  Raw sequences loaded: {len(raw_records):,}")

# ── QUALITY CONTROL ───────────────────────────────────────────
print("\nRunning quality control filters...")

clean_records = []
qc_stats = {
    "total_raw": len(raw_records),
    "removed_no_start_met": 0,
    "removed_internal_stop": 0,
    "removed_too_short": 0,
    "removed_ambiguous": 0,
    "passed": 0
}

MIN_LENGTH = 50  # Remove very short fragments (<50 aa = likely false predictions)

for rec in raw_records:
    seq = str(rec.seq).upper()
    
    # Remove sequences without start methionine
    if not seq.startswith("M"):
        qc_stats["removed_no_start_met"] += 1
        continue
    
    # Remove sequences with internal stop codons (AUGUSTUS sometimes includes *)
    seq_no_trail = seq.rstrip("*")
    if "*" in seq_no_trail:
        qc_stats["removed_internal_stop"] += 1
        continue
    
    # Clean trailing stop codon
    seq_clean = seq.rstrip("*")
    
    # Remove too-short sequences
    if len(seq_clean) < MIN_LENGTH:
        qc_stats["removed_too_short"] += 1
        continue
    
    # Remove heavily ambiguous sequences (>5% X)
    x_frac = seq_clean.count("X") / len(seq_clean)
    if x_frac > 0.05:
        qc_stats["removed_ambiguous"] += 1
        continue
    
    # Clean up record
    rec.seq = rec.seq.__class__(seq_clean)
    clean_records.append(rec)

qc_stats["passed"] = len(clean_records)

print(f"  QC results:")
print(f"    No start methionine:   {qc_stats['removed_no_start_met']:>6}")
print(f"    Internal stop codons:  {qc_stats['removed_internal_stop']:>6}")
print(f"    Too short (<{MIN_LENGTH} aa):   {qc_stats['removed_too_short']:>6}")
print(f"    >5% ambiguous (X):     {qc_stats['removed_ambiguous']:>6}")
print(f"    ─────────────────────────────")
print(f"    PASSED QC:             {qc_stats['passed']:>6,}")

# ── WRITE OUTPUTS ─────────────────────────────────────────────
SeqIO.write(raw_records, OUTPUT_FASTA, "fasta")
SeqIO.write(clean_records, OUTPUT_CLEAN, "fasta")
print(f"\n[OK] Raw proteins:   {OUTPUT_FASTA}")
print(f"[OK] Clean proteins: {OUTPUT_CLEAN}")

# ── PROTEOME STATISTICS ───────────────────────────────────────
lengths = [len(r.seq) for r in clean_records]
lengths_arr = np.array(lengths)

print(f"\nProteome statistics:")
print(f"  Min length:    {min(lengths)} aa")
print(f"  Max length:    {max(lengths):,} aa")
print(f"  Mean length:   {np.mean(lengths):.1f} aa")
print(f"  Median length: {np.median(lengths):.1f} aa")
print(f"  Proteins <100 aa: {sum(1 for l in lengths if l < 100)}")
print(f"  Proteins 100-500 aa: {sum(1 for l in lengths if 100 <= l <= 500)}")
print(f"  Proteins >1000 aa: {sum(1 for l in lengths if l > 1000)}")

# Target zone for our analysis
in_target = sum(1 for l in lengths if 135 <= l <= 230)
print(f"\n  Proteins in target zone (135–230 aa = 15–25 kDa): {in_target}")

# ── WRITE STATS FILE ──────────────────────────────────────────
with open(OUTPUT_STATS, "w") as f:
    f.write("BRAKER PREDICTED PROTEOME — QC STATISTICS\n")
    f.write("=" * 55 + "\n\n")
    f.write(f"BRAKER output file:  {BRAKER_AA}\n\n")
    f.write("QC FILTER RESULTS:\n")
    for k, v in qc_stats.items():
        f.write(f"  {k:<30}: {v:,}\n")
    f.write("\nSIZE STATISTICS (clean proteome):\n")
    f.write(f"  Min:    {min(lengths)} aa\n")
    f.write(f"  Max:    {max(lengths):,} aa\n")
    f.write(f"  Mean:   {np.mean(lengths):.1f} aa\n")
    f.write(f"  Median: {np.median(lengths):.1f} aa\n")
    f.write(f"\nTARGET ZONE (135-230 aa = 15-25 kDa): {in_target} proteins\n")
    f.write("\nFIRST 10 PREDICTED PROTEINS:\n")
    for rec in clean_records[:10]:
        f.write(f"  {rec.id:<25} {len(rec.seq)} aa\n")

print(f"[OK] Stats written: {OUTPUT_STATS}")

# ── PLOT ──────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
fig.suptitle(
    f"BRAKER2 De Novo Gene Prediction — P. microspora KFRD-2\n"
    f"Predicted Proteome Overview (n={len(clean_records):,} proteins)",
    fontsize=12, fontweight='bold'
)

# Left: full length distribution
ax = axes[0]
ax.hist(lengths_arr[lengths_arr < 2000], bins=80,
        color='#2E75B6', alpha=0.8, edgecolor='white')
ax.axvspan(135, 230, alpha=0.35, color='#E74C3C',
           label=f'Target zone\n(135–230 aa, n={in_target})')
ax.axvline(191, color='#E74C3C', lw=2, ls='--', label='21 kDa (191 aa)')
ax.set_xlabel("Protein length (aa)", fontsize=11)
ax.set_ylabel("Number of proteins", fontsize=11)
ax.set_title("Length Distribution (<2000 aa)", fontsize=11)
ax.legend(fontsize=9)
ax.grid(axis='y', alpha=0.3)

# Middle: zoomed target zone
ax = axes[1]
target_lengths = [l for l in lengths if 135 <= l <= 230]
ax.hist(target_lengths, bins=25, color='#E74C3C', alpha=0.8, edgecolor='white')
ax.axvline(191, color='#1F3864', lw=2.5, ls='--', label='21 kDa (191 aa)')
ax.set_xlabel("Protein length (aa)", fontsize=11)
ax.set_ylabel("Number of proteins", fontsize=11)
ax.set_title(f"Target Zone (135–230 aa)\nn={len(target_lengths)}", fontsize=11)
ax.legend(fontsize=9)
ax.grid(axis='y', alpha=0.3)
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
kda_ticks = [135, 155, 175, 191, 210, 230]
ax2.set_xticks(kda_ticks)
ax2.set_xticklabels([f'{v*110/1000:.1f}' for v in kda_ticks], fontsize=8)
ax2.set_xlabel("Approx. MW (kDa)", fontsize=9, color='gray')

# Right: QC pie chart
ax = axes[2]
qc_labels = ['Passed QC', 'No Met start', 'Internal stop', 'Too short', 'Ambiguous']
qc_vals = [
    qc_stats['passed'],
    qc_stats['removed_no_start_met'],
    qc_stats['removed_internal_stop'],
    qc_stats['removed_too_short'],
    qc_stats['removed_ambiguous']
]
qc_vals_nonzero = [(l, v) for l, v in zip(qc_labels, qc_vals) if v > 0]
labels_nz = [x[0] for x in qc_vals_nonzero]
vals_nz = [x[1] for x in qc_vals_nonzero]
colors_pie = ['#1E8449','#E74C3C','#A04000','#2E75B6','#6C3483'][:len(labels_nz)]
wedges, texts, autotexts = ax.pie(
    vals_nz, labels=labels_nz, colors=colors_pie,
    autopct='%1.1f%%', startangle=90,
    textprops={'fontsize': 9}
)
ax.set_title("QC Filter Results", fontsize=11)

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"[OK] Plot saved: {OUTPUT_PLOT}")

print("\n" + "=" * 55)
print(f" STEP 2 COMPLETE")
print(f" Clean predicted proteins: {len(clean_records):,}")
print(f" Target zone candidates:   {in_target}")
print(f"\n Next: python3 03_filter_by_size.py")
print("=" * 55)
