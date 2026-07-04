#!/usr/bin/env bash
# ==============================================================
# 06b_prepare_receptor.sh
# Converts AlphaFold .cif/.pdb to .pdbqt format for AutoDock Vina
#
# INPUT:  05_structure/g13490_best_model.cif  (or .pdb)
# OUTPUT: 07_docking/receptor.pdbqt
# ==============================================================

set -e

echo "========================================"
echo " STEP 6B: PREPARE RECEPTOR FOR DOCKING"
echo "========================================"

STRUCT_DIR="05_structure"
DOCK_DIR="07_docking"
mkdir -p "$DOCK_DIR"

# ── FIND STRUCTURE FILE ───────────────────────────────────────
STRUCT_FILE=""
for f in "$STRUCT_DIR"/g13490*.cif "$STRUCT_DIR"/g13490*.pdb \
          "$STRUCT_DIR"/*model_0*.cif "$STRUCT_DIR"/*model_0*.pdb \
          "$STRUCT_DIR"/fold_*.cif "$STRUCT_DIR"/fold_*.pdb; do
    if [[ -f "$f" ]]; then
        STRUCT_FILE="$f"
        break
    fi
done

if [[ -z "$STRUCT_FILE" ]]; then
    echo "[ERROR] No structure file found in $STRUCT_DIR/"
    echo "  Copy your AlphaFold structure:"
    echo "  cp ~/Downloads/fold_g13490_model_0.cif $STRUCT_DIR/g13490_best_model.cif"
    exit 1
fi

echo "[INFO] Structure file: $STRUCT_FILE"
EXT="${STRUCT_FILE##*.}"

# ── STEP 1: CIF → PDB (if needed) ────────────────────────────
PDB_CLEAN="$DOCK_DIR/receptor_clean.pdb"

if [[ "$EXT" == "cif" ]]; then
    echo ""
    echo "[INFO] Converting CIF → PDB..."
    
    # Try obabel first (best method)
    if command -v obabel &>/dev/null; then
        obabel "$STRUCT_FILE" -O "$DOCK_DIR/receptor_from_cif.pdb" 2>/dev/null
        grep "^ATOM" "$DOCK_DIR/receptor_from_cif.pdb" > "$PDB_CLEAN"
        echo "[OK] CIF converted with Open Babel"
    
    # Try python mmcif parser (fallback)
    else
        echo "[INFO] Open Babel not found — using Python mmCIF parser..."
        python3 - << 'PYEOF'
import sys, os
sys.path.insert(0, '.')
from Bio.PDB import MMCIFParser, PDBIO
import glob

# Find CIF
cif_files = glob.glob("05_structure/*.cif") + glob.glob("05_structure/fold_*.cif")
if not cif_files:
    print("[ERROR] No .cif file found")
    sys.exit(1)

cif_file = cif_files[0]
print(f"[INFO] Converting: {cif_file}")

parser = MMCIFParser(QUIET=True)
struct = parser.get_structure("g13490", cif_file)

io = PDBIO()
io.set_structure(struct)
io.save("07_docking/receptor_from_cif.pdb")
print("[OK] Saved to 07_docking/receptor_from_cif.pdb")
PYEOF
        grep "^ATOM" "$DOCK_DIR/receptor_from_cif.pdb" > "$PDB_CLEAN" 2>/dev/null || \
            cp "$DOCK_DIR/receptor_from_cif.pdb" "$PDB_CLEAN"
        echo "[OK] CIF converted with BioPython"
    fi
else
    # Already PDB — just clean it
    grep "^ATOM" "$STRUCT_FILE" > "$PDB_CLEAN"
    echo "[OK] PDB cleaned (ATOM records only): $PDB_CLEAN"
fi

NATOMS=$(wc -l < "$PDB_CLEAN")
echo "[INFO] Clean receptor: $NATOMS ATOM lines"

# ── STEP 2: ADD HYDROGENS ─────────────────────────────────────
echo ""
echo "[INFO] Adding hydrogens and computing Gasteiger charges..."
PDB_H="$DOCK_DIR/receptor_with_H.pdb"
RECEPTOR_PDBQT="$DOCK_DIR/receptor.pdbqt"

if ! command -v obabel &>/dev/null; then
    echo "[ERROR] Open Babel not installed."
    echo "  Install: conda install -c conda-forge openbabel"
    echo "  OR:      sudo apt install openbabel"
    exit 1
fi

# Add hydrogens at pH 7.4 (physiological)
obabel "$PDB_CLEAN" \
    -O "$PDB_H" \
    -p 7.4 \
    -h 2>/dev/null
echo "[OK] Hydrogens added: $PDB_H"

# Convert to PDBQT (Vina format)
obabel "$PDB_H" \
    -O "$RECEPTOR_PDBQT" \
    --partialcharge gasteiger \
    -xr 2>/dev/null

if [[ -f "$RECEPTOR_PDBQT" ]]; then
    NLINES=$(wc -l < "$RECEPTOR_PDBQT")
    echo "[OK] Receptor PDBQT: $RECEPTOR_PDBQT ($NLINES lines)"
else
    echo "[ERROR] receptor.pdbqt not created. Check Open Babel installation."
    exit 1
fi

# ── STEP 3: VERIFY PDBQT ─────────────────────────────────────
echo ""
echo "--- RECEPTOR PDBQT VERIFICATION ---"
echo "First 5 ATOM lines:"
grep "^ATOM" "$RECEPTOR_PDBQT" | head -5
echo ""
echo "Atom types found:"
grep "^ATOM" "$RECEPTOR_PDBQT" | awk '{print $NF}' | sort | uniq -c | sort -rn | head -10

echo ""
echo "========================================"
echo " STEP 6B COMPLETE"
echo " Receptor: $RECEPTOR_PDBQT"
echo " Next: bash 07a_prepare_ligands.sh"
echo "========================================"
