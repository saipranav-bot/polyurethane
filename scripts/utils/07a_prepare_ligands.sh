#!/usr/bin/env bash
# ==============================================================
# 07a_prepare_ligands.sh
# Generates 3D ligand structures from SMILES and converts
# to PDBQT format for AutoDock Vina docking.
#
# LIGANDS (polyurethane substrate fragments):
#
#   L1: Ethyl carbamate       — simplest urethane unit, benchmarked
#   L2: Diethyl adipate       — polyester PUR ester bond model
#   L3: Impranil DLN fragment — more realistic PUR substrate
#   L4: MDI fragment          — aromatic isocyanate-derived unit
#   L5: DMSO (control)        — should show low/no binding
#
# WHY THESE LIGANDS:
#   Russell 2011: enzyme cleaves ester bonds in polyester PUR
#   Mechanism: Ser nucleophile attacks carbonyl carbon of ester
#   Ligands chosen to probe: (a) urethane bonds (L1)
#                             (b) ester bonds (L2, L3)
#                             (c) aromatic PUR unit (L4)
#                             (d) negative control (L5)
# ==============================================================

set -e

echo "========================================"
echo " STEP 7A: PREPARE DOCKING LIGANDS"
echo "========================================"

DOCK_DIR="07_docking"
mkdir -p "$DOCK_DIR/ligands"

# ── CHECK OBABEL ──────────────────────────────────────────────
if ! command -v obabel &>/dev/null; then
    echo "[ERROR] Open Babel not installed"
    echo "  Install: conda install -c conda-forge openbabel"
    echo "  OR:      sudo apt install openbabel"
    exit 1
fi
echo "[OK] Open Babel: $(obabel --version 2>/dev/null | head -1)"

# ── LIGAND PREPARATION FUNCTION ───────────────────────────────
prepare_ligand() {
    NAME="$1"
    SMILES="$2"
    DESC="$3"
    
    SDF="$DOCK_DIR/ligands/${NAME}.sdf"
    PDBQT="$DOCK_DIR/ligands/${NAME}.pdbqt"
    
    echo ""
    echo "Preparing: $NAME"
    echo "  Desc:   $DESC"
    echo "  SMILES: $SMILES"
    
    # Generate 3D structure from SMILES (best = force field minimized)
    obabel -:"$SMILES" \
        --gen3D best \
        -O "$SDF" 2>/dev/null
    
    if [[ ! -f "$SDF" ]]; then
        echo "  [WARNING] SDF generation failed, trying without --gen3D best..."
        obabel -:"$SMILES" --gen3D -O "$SDF" 2>/dev/null
    fi
    
    if [[ ! -f "$SDF" ]]; then
        echo "  [ERROR] Could not generate 3D structure for $NAME"
        return 1
    fi
    
    # MMFF94 force field minimization for better geometry
    obabel "$SDF" \
        -O "$SDF" \
        --minimize --ff MMFF94 --steps 2000 2>/dev/null || true
    
    # Convert SDF → PDBQT (adds atom types + Gasteiger charges)
    obabel "$SDF" \
        -O "$PDBQT" \
        --partialcharge gasteiger 2>/dev/null
    
    if [[ -f "$PDBQT" ]]; then
        NATOMS=$(grep -c "^ATOM\|^HETATM" "$PDBQT" 2>/dev/null || echo "?")
        echo "  [OK] $PDBQT  ($NATOMS heavy atoms)"
    else
        echo "  [ERROR] PDBQT not created"
    fi
}

# ── PREPARE ALL LIGANDS ───────────────────────────────────────

# L1: Ethyl carbamate — urethane bond -O-C(=O)-NH-
#     SIMPLEST model: tested as polyurethane surrogate in docking studies
prepare_ligand \
    "L1_ethyl_carbamate" \
    "CCOC(N)=O" \
    "Ethyl carbamate (urethane unit) — simplest PUR model"

# L2: Diethyl adipate — polyester ester bond
#     Models the C(=O)-O-C linkage the enzyme hydrolyzes
prepare_ligand \
    "L2_diethyl_adipate" \
    "CCOC(=O)CCCCC(=O)OCC" \
    "Diethyl adipate (polyester ester bond model)"

# L3: Ethylene adipate diurethane — short PUR oligomer
#     Has both ester AND urethane bonds — most realistic
prepare_ligand \
    "L3_pur_fragment" \
    "O=C(NCCNC(=O)OCCCCC(=O)O)OCC" \
    "PUR fragment (ester+urethane — most realistic substrate)"

# L4: 4,4'-Methylenediphenyl urethane fragment — aromatic PUR
#     Models MDI-based polyurethane (most common industrial PUR)
prepare_ligand \
    "L4_mdi_fragment" \
    "O=C(Nc1ccc(Cc2ccc(NC(=O)OCC)cc2)cc1)OCC" \
    "MDI-based PUR fragment (aromatic urethane)"

# L5: DMSO — negative control (should NOT bind well)
prepare_ligand \
    "L5_dmso_control" \
    "CS(C)=O" \
    "DMSO (negative control — no ester/urethane bonds)"

# ── SUMMARY ───────────────────────────────────────────────────
echo ""
echo "--- LIGAND PREPARATION SUMMARY ---"
echo ""
printf "%-25s %-8s %s\n" "Ligand" "Status" "File"
printf "%-25s %-8s %s\n" "------" "------" "----"
for NAME in L1_ethyl_carbamate L2_diethyl_adipate L3_pur_fragment L4_mdi_fragment L5_dmso_control; do
    PDBQT="$DOCK_DIR/ligands/${NAME}.pdbqt"
    if [[ -f "$PDBQT" ]]; then
        printf "%-25s %-8s %s\n" "$NAME" "✅ OK" "$PDBQT"
    else
        printf "%-25s %-8s %s\n" "$NAME" "❌ FAIL" "not created"
    fi
done

echo ""
echo "========================================"
echo " STEP 7A COMPLETE"
echo " Next: bash 07b_run_docking_all.sh"
echo "========================================"
