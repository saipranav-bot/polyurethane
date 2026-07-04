import os

INPUT = "03_secretome/signalp6_output.txt"
OUTPUT_IDS = "03_secretome/secreted_ids.txt"
OUTPUT_FASTA_HINT = "03_secretome/secreted_candidates.faa"

print("📦 SignalP 6 parsing (correct format)")

if not os.path.exists(INPUT):
    raise FileNotFoundError(INPUT)

secreted = []

with open(INPUT) as f:
    for line in f:
        line = line.strip()

        # skip headers
        if not line or line.startswith("#"):
            continue

        cols = line.split()

        if len(cols) < 2:
            continue

        protein_id = cols[0]
        prediction = cols[1].strip().upper()

        if prediction == "SP":
            secreted.append(protein_id)

print("Secreted proteins found:", len(secreted))

# save IDs
with open(OUTPUT_IDS, "w") as out:
    for s in secreted:
        out.write(s + "\n")

print("✅ IDs saved:", OUTPUT_IDS)

# NOTE: FASTA reconstruction would require original proteome mapping
print("⚠️ Note: FASTA filtering step should use ID list + proteome FASTA")
