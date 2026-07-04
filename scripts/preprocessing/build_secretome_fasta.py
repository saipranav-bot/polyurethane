from Bio import SeqIO

INPUT_FASTA = "02_size_filter/candidates_sizefiltered.faa"
ID_FILE = "03_secretome/secreted_ids.txt"
OUTPUT_FASTA = "03_secretome/secreted_candidates.faa"

print("🧬 Building secretome FASTA")

# load IDs
with open(ID_FILE) as f:
    secret_ids = set(line.strip() for line in f if line.strip())

print("Secreted IDs:", len(secret_ids))

records = list(SeqIO.parse(INPUT_FASTA, "fasta"))

filtered = [r for r in records if r.id in secret_ids]

print("Matched FASTA sequences:", len(filtered))

SeqIO.write(filtered, OUTPUT_FASTA, "fasta")

print("✅ Secretome FASTA written:", OUTPUT_FASTA)
