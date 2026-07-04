#!/usr/bin/env bash

set -e

echo "========================================"
echo " STEP 1: DOWNLOAD P. microspora KFRD-2"
echo "========================================"

cd 01_genome

ACCESSION="GCA_041412015.1"
ZIPFILE="pestalotiopsis_kfrd2.zip"
OUTDIR="pestalotiopsis_kfrd2"

mkdir -p $OUTDIR

if command -v datasets &>/dev/null; then

    echo "[INFO] Using NCBI datasets CLI..."

    datasets download genome accession $ACCESSION \
        --include genome,protein,gff3 \
        --filename $ZIPFILE

    unzip -o $ZIPFILE -d $OUTDIR

else
    echo "ERROR: Install NCBI datasets first"
    exit 1
fi


PROTEIN_FASTA=$(find $OUTDIR -name "*.faa" | head -1)

if [[ -z "$PROTEIN_FASTA" ]]; then
    echo "[ERROR] Protein FASTA not found"
    exit 1
fi


cp "$PROTEIN_FASTA" protein.faa


echo ""
echo "--- PROTEIN STATS ---"

TOTAL=$(grep -c ">" protein.faa)

echo "Total proteins: $TOTAL"


if [[ $TOTAL -lt 10000 ]]; then
    echo "[WARNING] Protein count low"
else
    echo "[OK] Protein count looks good"
fi


echo ""
echo "First 3 headers:"
grep ">" protein.faa | head -3


echo ""
echo "========================================"
echo " COMPLETE"
echo " Output: 01_genome/protein.faa"
echo "========================================"
