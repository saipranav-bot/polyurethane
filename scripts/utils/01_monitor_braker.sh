#!/usr/bin/env bash
# ==============================================================
# 01_monitor_braker.sh
# Run this NOW while BRAKER is running.
# It monitors progress and tells you exactly when to move on.
# ==============================================================

BRAKER_DIR="${1:-~/polyurethanase_pipeline/01_genome/braker_out}"
BRAKER_DIR=$(eval echo "$BRAKER_DIR")  # expand ~

echo "========================================"
echo " BRAKER MONITOR"
echo " Watching: $BRAKER_DIR"
echo "========================================"
echo ""

# ── FUNCTION: check BRAKER status ────────────────────────────
check_braker() {
    echo "--- STATUS CHECK: $(date) ---"
    
    # Check if optimize_augustus is still running
    OPTIM_PID=$(pgrep -f "optimize_augustus" 2>/dev/null | head -1)
    BRAKER_PID=$(pgrep -f "braker.pl" 2>/dev/null | head -1)
    AUGUSTUS_PID=$(pgrep -f "augustus" 2>/dev/null | head -1)
    
    if [[ -n "$OPTIM_PID" ]]; then
        echo "[RUNNING] optimize_augustus.pl  (PID: $OPTIM_PID)"
    fi
    if [[ -n "$BRAKER_PID" ]]; then
        echo "[RUNNING] braker.pl             (PID: $BRAKER_PID)"
    fi
    if [[ -n "$AUGUSTUS_PID" ]]; then
        echo "[RUNNING] augustus              (PID: $AUGUSTUS_PID)"
    fi
    if [[ -z "$OPTIM_PID" && -z "$BRAKER_PID" && -z "$AUGUSTUS_PID" ]]; then
        echo "[INFO] No BRAKER/AUGUSTUS processes detected"
    fi
    
    echo ""
    
    # Check for output files
    echo "Output files:"
    for f in braker.aa braker.gtf braker.codingseq augustus.hints.gff; do
        FPATH="$BRAKER_DIR/$f"
        if [[ -f "$FPATH" ]]; then
            SIZE=$(du -sh "$FPATH" | cut -f1)
            LINES=$(wc -l < "$FPATH")
            echo "  ✅ $f  ($SIZE, $LINES lines)"
        else
            echo "  ⏳ $f  (not yet created)"
        fi
    done
    
    echo ""
    
    # Check current accuracy from optimize log
    OPTIM_LOG="$BRAKER_DIR/optimize_augustus.stdout"
    if [[ -f "$OPTIM_LOG" ]]; then
        LAST_TARGET=$(grep "^targets" "$OPTIM_LOG" 2>/dev/null | tail -1)
        if [[ -n "$LAST_TARGET" ]]; then
            BEST=$(echo "$LAST_TARGET" | tr ' ' '\n' | grep -E '^0\.' | sort -n | tail -1)
            echo "Current best accuracy: $BEST (from optimize_augustus log)"
        fi
        LAST_PARAM=$(grep "^improving parameter" "$OPTIM_LOG" 2>/dev/null | tail -1)
        echo "Last optimized: $LAST_PARAM"
    fi
    
    echo ""
}

# ── DETECT BRAKER DIR ─────────────────────────────────────────
if [[ ! -d "$BRAKER_DIR" ]]; then
    echo "[INFO] Directory $BRAKER_DIR not found."
    echo "Searching for braker output directory..."
    FOUND=$(find ~ -name "braker.aa" 2>/dev/null | head -1)
    if [[ -n "$FOUND" ]]; then
        BRAKER_DIR=$(dirname "$FOUND")
        echo "[FOUND] BRAKER output at: $BRAKER_DIR"
    else
        echo "BRAKER output not found yet."
        echo "Common locations to check:"
        echo "  ~/polyurethanase_pipeline/01_genome/braker_out/"
        echo "  ~/braker_out/"
        echo "  ./braker_out/"
    fi
fi

# ── WATCH MODE ───────────────────────────────────────────────
echo "Checking status every 5 minutes. Press Ctrl+C to stop."
echo "When braker.aa appears with protein sequences, run:"
echo "  python3 02_extract_proteins.py"
echo ""

while true; do
    check_braker
    
    # Check if braker.aa is ready and non-empty
    BRAKER_AA="$BRAKER_DIR/braker.aa"
    if [[ -f "$BRAKER_AA" ]]; then
        PROT_COUNT=$(grep -c ">" "$BRAKER_AA" 2>/dev/null || echo 0)
        if [[ $PROT_COUNT -gt 1000 ]]; then
            echo "========================================"
            echo " ✅ BRAKER FINISHED!"
            echo " braker.aa has $PROT_COUNT protein sequences"
            echo ""
            echo " NEXT STEP:"
            echo "   python3 02_extract_proteins.py $BRAKER_AA"
            echo "========================================"
            break
        fi
    fi
    
    # Also check if there are error logs
    ERR_LOG="$BRAKER_DIR/errors"
    if [[ -d "$ERR_LOG" ]]; then
        ERR_COUNT=$(ls "$ERR_LOG"/*.err 2>/dev/null | wc -l)
        if [[ $ERR_COUNT -gt 0 ]]; then
            echo "[WARNING] Error files found in $ERR_LOG/"
            ls "$ERR_LOG"/*.err 2>/dev/null | head -5
        fi
    fi
    
    sleep 300  # check every 5 minutes
done
