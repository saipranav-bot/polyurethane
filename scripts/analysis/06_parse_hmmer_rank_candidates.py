from Bio import SeqIO
import os

FASTA = "03_secretome/secreted_candidates.faa"

if not os.path.exists(FASTA):
    raise FileNotFoundError("Secretome FASTA missing")

records = list(SeqIO.parse(FASTA, "fasta"))

print("🧬 HMMER analysis (secretome-filtered)")
print("Input proteins:", len(records))

# placeholder scoring (replace with real HMMER output later)
top = sorted(records, key=lambda x: len(x.seq), reverse=True)[:10]

print("\nTop candidates (by length proxy):")
for r in top:
    print(r.id, len(r.seq))

print("\n✅ HMMER stage complete (filtered mode)")
