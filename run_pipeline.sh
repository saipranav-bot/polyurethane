cd ~/polyurethanase_pipeline

# ----------------------------
# CREATE MAIN README (PUBLICATION-GRADE)
# ----------------------------
cat > README.md << 'EOF'
# Polyurethanase Discovery Pipeline

## Overview
This project implements a computational workflow to identify and characterize putative polyurethanase enzymes from genomic datasets using a multi-stage bioinformatics pipeline including genome processing, protein extraction, size filtering, secretome prediction, HMMER-based homology search, and SignalP-based secretion signal detection.

---

## Pipeline Steps

1. Genome preprocessing
2. Protein extraction
3. Size-based filtering
4. Secretome prediction
5. HMMER homology search
6. Signal peptide prediction

---

## Project Structure

- `data/` → raw and processed input datasets  
- `results/` → pipeline outputs organized by stage  
- `scripts/` → modular analysis code  
- `workflow/` → execution pipeline  
- `config/` → reproducible parameter settings  
- `figures/` → visualization outputs  
- `docs/` → documentation

---

## Setup

```bash
git clone <repo_url>
cd polyurethanase_pipeline

conda env create -f environment.yml
conda activate polyurethanase
