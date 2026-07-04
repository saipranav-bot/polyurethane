#!/bin/bash

echo "===== POLYURETHANASE PIPELINE CLEANUP START ====="

BASE=~/polyurethanase_pipeline

cd $BASE || exit

echo "1. Removing failed BRAKER runs..."
rm -rf 01_genome/braker_out_failed_run

echo "2. Removing temporary GeneMark training files..."
rm -rf 01_genome/data
rm -rf 01_genome/tmp_opt_P_microspora_KFRD2

echo "3. Removing downloaded archives (already extracted)..."
rm -f 01_genome/*.zip
rm -f 01_genome/*.tar.gz

echo "4. Removing duplicate report..."
rm -f 01_genome/PROJECT_REPORT.md

echo "5. Cleaning unnecessary GeneMark logs..."
rm -f 01_genome/gmes.log
rm -f 01_genome/gmes_linux_64_4.tar.gz
rm -f 01_genome/gms2_linux_64.tar.gz

echo "6. Organizing reports..."

mkdir -p reports

mv PROJECT_REPORT.md reports/PROJECT_REPORT.md 2>/dev/null
mv README.md reports/README.md 2>/dev/null
mv README_1.md reports/README_1.md 2>/dev/null
mv README.txt reports/README.txt 2>/dev/null

echo "7. Creating clean result folders..."

mkdir -p results

mv 08_figures results/figures 2>/dev/null
mv 07_docking results/docking 2>/dev/null
mv 05_structure results/structure 2>/dev/null
mv 06_active_site results/active_site 2>/dev/null

echo "8. Creating analysis folder..."

mkdir -p scripts

mv *.py scripts/ 2>/dev/null
mv *.sh scripts/ 2>/dev/null

echo "9. Removing empty folders..."

find $BASE -type d -empty -delete


echo "===== CLEANUP FINISHED ====="

echo ""
echo "Final structure:"
tree -L 2
