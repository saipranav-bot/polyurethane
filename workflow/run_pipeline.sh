#!/bin/bash

set -e

echo "🚀 FULL REPRODUCIBLE POLYURETHANASE PIPELINE"

run () {
    echo "▶ $1"
    python "$1"
}

# --------------------------
# CORE PIPELINE
# --------------------------

run scripts/preprocessing/02_extract_proteins.py
run scripts/preprocessing/03_filter_by_size.py

# 🔥 NEW: fully automated SignalP
run scripts/preprocessing/run_signalp.py

run scripts/preprocessing/04_filter_secretome.py

run scripts/analysis/06_parse_hmmer_rank_candidates.py
run scripts/analysis/07c_analyze_docking_CORRECTED.py

echo "✅ PIPELINE COMPLETE (100% REPRODUCIBLE MODE)"
