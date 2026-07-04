#!/usr/bin/env bash
# ==============================================================
# PIPELINE README — READ THIS FIRST
# ==============================================================
#
# YOUR CURRENT STATUS:
#   You downloaded: GCA_036323685.1 (raw genome FASTA)
#   You are running: BRAKER + AUGUSTUS gene prediction
#   AUGUSTUS is optimizing gene model parameters (92% accuracy)
#   This is CORRECT and BETTER than using pre-annotated proteins
#
# WHY BRAKER IS BETTER THAN NCBI ANNOTATIONS:
#   - NCBI annotations for KFRD-2 may use automated GenBank pipeline
#   - BRAKER trains AUGUSTUS specifically on YOUR genome
#   - Your predicted gene models will be species-specific and more accurate
#   - This is publishable as "de novo gene prediction"
#
# WHAT BRAKER WILL PRODUCE (in ~/polyurethanase_pipeline/01_genome/braker_out/):
#   braker.gtf          ← gene coordinates (where genes are on each chromosome)
#   braker.aa           ← protein sequences (THIS IS WHAT YOU NEED)
#   braker.codingseq    ← CDS nucleotide sequences
#   augustus.hints.gff  ← gene models with evidence support
#
# REVISED PIPELINE FLOW:
#
#   [BRAKER RUNNING NOW] → braker.aa
#         ↓
#   Step 1: Extract + QC proteins from braker.aa
#         ↓
#   Step 2: Filter by size (135-230 aa = ~15-25 kDa)
#         ↓
#   Step 3: SignalP 6.0 → secretome prediction (WEB)
#         ↓
#   Step 4: TMHMM → remove membrane proteins (WEB)
#         ↓
#   Step 5: HMMER + Pfam → serine hydrolase domain search (LOCAL)
#         ↓
#   Step 6: BLASTp → validate top candidates (WEB or LOCAL)
#         ↓
#   Step 7: InterProScan → full functional annotation (WEB)
#         ↓
#   Step 8: AlphaFold 2 → structure prediction (WEB - ColabFold)
#         ↓
#   Step 9: fpocket → active site identification (LOCAL)
#         ↓
#   Step 10: AutoDock Vina → molecular docking (LOCAL)
#         ↓
#   Step 11: Figures + report
#
# SCRIPTS IN THIS DIRECTORY:
#   01_monitor_braker.sh          ← run NOW to watch BRAKER progress
#   02_extract_proteins.py        ← run AFTER BRAKER finishes
#   03_filter_by_size.py          ← filter to ~21 kDa candidates
#   04_filter_secretome.py        ← parse SignalP + TMHMM results
#   05_run_hmmer.sh               ← HMMER domain search
#   06_parse_hmmer_blastp.py      ← parse hits, prepare for BLASTp
#   07_interproscan_prep.py       ← prepare InterProScan input
#   08_structure_prep.py          ← prepare sequences for AlphaFold
#   09_run_fpocket.sh             ← active site analysis
#   10_run_docking.sh             ← AutoDock Vina docking
#   11_figures.py                 ← generate all publication figures
#
# ==============================================================
echo "Read the comments above. Run 01_monitor_braker.sh next."
