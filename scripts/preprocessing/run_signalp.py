import os
import shutil

INPUT = "02_size_filter/candidates_sizefiltered.faa"
OUTPUT_DIR = "03_secretome"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "signalp6_output.txt")

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("🚀 SignalP STEP (web-based integration mode)")

if not os.path.exists(INPUT):
    raise FileNotFoundError(INPUT)

print("\n📌 MANUAL STEP REQUIRED (REPRODUCIBLE MODE)")
print("1. Go to SignalP 6.0 web server")
print("2. Upload:", INPUT)
print("3. Select: 'Short output'")
print("4. Download results")
print("5. Save as:")
print("   →", OUTPUT_FILE)

print("\n⏳ Pipeline will continue AFTER file exists")

# block execution until file exists
if not os.path.exists(OUTPUT_FILE):
    print("\n⚠️ Waiting for SignalP output file...")
    exit(1)

print("✅ SignalP output detected")
