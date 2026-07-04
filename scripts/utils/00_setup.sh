#!/usr/bin/env bash
# ============================================================
# 00_setup.sh
# Run this FIRST on your Ubuntu machine (scrna_t2d conda env)
# Sets up all tools needed for the pipeline
# ============================================================

set -e  # stop on error

echo "========================================"
echo " POLYURETHANASE PIPELINE — SETUP"
echo "========================================"
echo ""

# ── conda env check ─────────────────────────────────────────
if [[ "$CONDA_DEFAULT_ENV" != "scrna_t2d" ]]; then
    echo "[WARNING] You are NOT in scrna_t2d conda env"
    echo "  Run: conda activate scrna_t2d"
    echo "  Then re-run this script"
    echo "  (continuing anyway but tools may install to wrong env)"
fi

echo "[1/7] Installing Python packages..."
pip install biopython requests pandas numpy matplotlib seaborn \
    --break-system-packages -q
echo "  -> biopython, pandas, numpy, matplotlib, seaborn: OK"

echo ""
echo "[2/7] Installing HMMER (for Pfam domain search)..."
conda install -c bioconda hmmer -y -q 2>/dev/null || \
    sudo apt-get install -y hmmer -q 2>/dev/null || \
    echo "  [!] HMMER not auto-installed. Run: conda install -c bioconda hmmer"
hmmscan -h 2>/dev/null | head -1 && echo "  -> HMMER: OK" || echo "  [!] HMMER not found — install manually"

echo ""
echo "[3/7] Installing NCBI Datasets CLI..."
# Try conda first, then pip
conda install -c conda-forge ncbi-datasets-cli -y -q 2>/dev/null || \
    pip install ncbi-datasets-cli --break-system-packages -q 2>/dev/null || \
    echo "  [!] datasets CLI not auto-installed"
datasets --version 2>/dev/null && echo "  -> NCBI datasets CLI: OK" || \
    echo "  [!] datasets CLI not found. Download binary from:"
    echo "      https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets"
    echo "      chmod +x datasets && mv datasets ~/bin/"

echo ""
echo "[4/7] Installing AutoDock Vina..."
conda install -c conda-forge vina -y -q 2>/dev/null || \
    pip install vina --break-system-packages -q 2>/dev/null || \
    echo "  [!] Vina not auto-installed"
vina --version 2>/dev/null && echo "  -> AutoDock Vina: OK" || \
    echo "  [!] Vina not found. Install: conda install -c conda-forge vina"

echo ""
echo "[5/7] Installing Open Babel..."
conda install -c conda-forge openbabel -y -q 2>/dev/null || \
    sudo apt-get install -y openbabel -q 2>/dev/null || \
    echo "  [!] Open Babel not auto-installed"
obabel --version 2>/dev/null | head -1 && echo "  -> Open Babel: OK" || \
    echo "  [!] Open Babel not found. Run: sudo apt install openbabel"

echo ""
echo "[6/7] Installing fpocket..."
conda install -c bioconda fpocket -y -q 2>/dev/null || \
    echo "  [!] fpocket not auto-installed"
fpocket --version 2>/dev/null && echo "  -> fpocket: OK" || \
    echo "  [!] fpocket not found. Run: conda install -c bioconda fpocket"

echo ""
echo "[7/7] Creating pipeline directory structure..."
mkdir -p 01_genome 02_size_filter 03_secretome \
         04_serine_hydrolase 05_structure \
         06_active_site 07_docking 08_figures logs
echo "  -> Directories created: OK"

echo ""
echo "========================================"
echo " SETUP COMPLETE"
echo " Next: bash 01_download_genome.sh"
echo "========================================"
