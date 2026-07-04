# Polyurethane Enzyme Discovery Pipeline (In Silico)

## Overview

This project presents a reproducible bioinformatics workflow for the identification and prioritization of putative polyurethane-degrading enzymes from a fungal genome. The pipeline integrates genome-derived protein prediction, secretome filtering, protein family annotation, structural analysis, and molecular docking to identify candidates for future experimental validation.

> **Note:** This is an in silico computational discovery pipeline. Computational predictions should be validated experimentally before drawing biological conclusions.

---

## Objectives

* Predict proteins from a fungal genome.
* Identify secreted proteins with enzyme-like characteristics.
* Prioritize candidates based on sequence, structural, and docking analyses.
* Build a reproducible workflow suitable for future extension and publication.

---

## Workflow

```
Genome Assembly
       │
       ▼
BRAKER Protein Prediction
       │
       ▼
Quality Control
       │
       ▼
Protein Size Filtering
(135–230 aa)
       │
       ▼
SignalP 6 Secretome Prediction
       │
       ▼
Secreted Protein Set
       │
       ▼
Pfam/HMMER Functional Annotation
       │
       ▼
Candidate Ranking
       │
       ▼
Structure Prediction
       │
       ▼
Catalytic Site Analysis
       │
       ▼
Molecular Docking
       │
       ▼
Final Candidate Prioritization
```

---

## Current Pipeline Status

Completed:

* Genome-derived protein extraction
* Protein quality control
* Size-based candidate filtering
* SignalP 6 secretome prediction integration
* Secretome FASTA construction
* Pfam annotation using HMMER
* Enzyme candidate ranking
* Structural analysis
* Active-site inspection
* Molecular docking analysis
* Publication-ready figures and reports

Planned Improvements:

* Automated SignalP installation
* TMHMM integration
* Conserved motif analysis
* Phylogenetic analysis
* CAZyme annotation
* Multi-criteria scoring framework
* Wet-lab validation (future work)

---

## Repository Structure

```
scripts/
    preprocessing/
    analysis/
    utils/

workflow/
    run_pipeline.sh

config/
    pipeline_config.yaml

publication_ready/
    sequences/
    structure/
    active_site/
    docking/
    figures/
    reports/

08_figures/
```

---

## Software

* Ubuntu Linux
* Python 3
* BRAKER
* SignalP 6
* HMMER 3
* Pfam
* AutoDock Vina
* Biopython
* pandas
* NumPy

---

## Current Results

Pipeline summary:

* Raw predicted proteins: **19,152**
* Quality-filtered proteins: **19,095**
* Size-filtered proteins (135–230 aa): **2,150**
* SignalP-positive proteins identified
* Secretome FASTA generated
* Pfam/HMMER annotation completed
* Enzyme candidates ranked
* Docking analyses completed
* Publication-ready figures generated

---

## Reproducibility

The project includes:

* Automated preprocessing scripts
* Analysis scripts
* Workflow runner
* Configuration files
* Publication-ready outputs

Large external databases (e.g., Pfam) are intentionally excluded from the repository and should be downloaded separately.

---

## Disclaimer

This repository presents computational predictions generated using established bioinformatics tools. These results represent candidate enzymes and require experimental validation before biological or biochemical conclusions can be made.

---

## Author

**V. Sai Pranav**

GitHub: https://github.com/saipranav-bot

---

## License

This project is released under the MIT License.
