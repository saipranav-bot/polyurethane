import re

hmmer_file = "results_hmmer.txt"

keywords = [
    "hydrolase",
    "esterase",
    "lipase",
    "cutinase",
    "carboxylesterase",
    "alpha/beta hydrolase",
    "serine",
]

hits = []

with open(hmmer_file) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        for k in keywords:
            if k.lower() in line.lower():
                hits.append(line.strip())
                break

print("🧬 Enzyme-related HMMER hits:", len(hits))

print("\nTop matches:")
for h in hits[:20]:
    print(h)
