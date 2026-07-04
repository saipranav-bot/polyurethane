#!/usr/bin/env bash
# ==============================================================
# 05_run_hmmer.sh
# Runs HMMER hmmscan against Pfam-A to find serine hydrolase
# domains in secreted ~21 kDa candidates from BRAKER
# ==============================================================

set -e

echo "========================================"
echo " STEP 5: HMMER SERINE HYDROLASE SEARCH"
echo "========================================"
echo ""

# ── INPUTS ───────────────────────────────────────────────────
QUERY="03_secretome/secreted_candidates.faa"
OUTDIR="04_hmmer"
PFAM_DIR="$OUTDIR/pfam_db"
PFAM_HMM="$PFAM_DIR/Pfam-A.hmm"

mkdir -p "$OUTDIR" "$PFAM_DIR"

# ── CHECK INPUT ───────────────────────────────────────────────
if [[ ! -f "$QUERY" ]]; then
    echo "[ERROR] $QUERY not found"
    echo "  Run 04_filter_secretome.py first"
    exit 1
fi

NSEQ=$(grep -c ">" "$QUERY")
echo "[INFO] Input: $QUERY ($NSEQ sequences)"

# ── CHECK HMMER ───────────────────────────────────────────────
if ! command -v hmmscan &>/dev/null; then
    echo "[ERROR] HMMER not installed."
    echo ""
    echo "Install with ONE of:"
    echo "  conda install -c bioconda hmmer"
    echo "  sudo apt install hmmer"
    echo "  mamba install -c bioconda hmmer"
    exit 1
fi
HMMER_VER=$(hmmscan -h 2>&1 | grep "# HMMER" | head -1)
echo "[OK] HMMER: $HMMER_VER"

# ── DOWNLOAD PFAM-A (if not already downloaded) ───────────────
if [[ -f "${PFAM_HMM}.h3i" ]]; then
    echo "[OK] Pfam-A already downloaded and indexed"
else
    echo ""
    echo "[INFO] Downloading Pfam-A HMM database..."
    echo "  Size: ~500 MB | Time: 5-20 min depending on connection"
    echo ""
    
    wget --show-progress -q \
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz" \
        -O "$PFAM_DIR/Pfam-A.hmm.gz"
    
    echo "[INFO] Decompressing..."
    gunzip -f "$PFAM_DIR/Pfam-A.hmm.gz"
    
    echo "[INFO] Indexing (hmmpress)..."
    hmmpress "$PFAM_HMM"
    echo "[OK] Pfam-A ready"
fi

# ── RUN HMMSCAN (all Pfam domains) ────────────────────────────
echo ""
echo "[INFO] Running hmmscan (all Pfam domains)..."
echo "  This takes 5-30 min for a small secretome"
echo ""

hmmscan \
    --cpu 4 \
    --tblout  "$OUTDIR/hmmer_tblout.txt" \
    --domtblout "$OUTDIR/hmmer_domtblout.txt" \
    --acc \
    -E 0.01 \
    "$PFAM_HMM" \
    "$QUERY" \
    > "$OUTDIR/hmmer_full_output.txt" 2>&1

echo "[OK] HMMSCAN complete"

# ── QUICK PREVIEW ─────────────────────────────────────────────
echo ""
echo "--- PREVIEW: Top hits (E < 1e-5) ---"
grep -v "^#" "$OUTDIR/hmmer_tblout.txt" | \
    awk '$5 < 1e-5 {printf "%-20s %-12s %-10s %s\n", $1, $3, $5, $1}' | \
    sort -k3 -n | head -20

echo ""
echo "--- Total significant hits (E < 1e-5) ---"
NHITS=$(grep -v "^#" "$OUTDIR/hmmer_tblout.txt" | awk '$5 < 1e-5' | wc -l)
echo "  $NHITS hits"

echo ""
echo "========================================"
echo " STEP 5 COMPLETE"
echo " Output: $OUTDIR/hmmer_tblout.txt"
echo " Next:   python3 06_parse_hmmer_rank_candidates.py"
echo "========================================"
