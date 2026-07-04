from collections import defaultdict

hmmer_file = "results_hmmer.txt"

priority_domains = {
    "Cutinase": 5,
    "Abhydrolase": 4,
    "Lipase": 4,
    "Esterase": 4,
    "Phospholip": 3,
    "Glyco_hydro": 1,
}

scores = defaultdict(float)

current_query = None

with open(hmmer_file) as f:
    for line in f:
        if line.startswith("Query:"):
            current_query = line.split()[1]

        for domain, weight in priority_domains.items():
            if domain.lower() in line.lower():
                if current_query:
                    scores[current_query] += weight

print("\n🧬 TOP ENZYME CANDIDATES:\n")

ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)

for i, (prot, score) in enumerate(ranked[:20], 1):
    print(f"{i:02d}. {prot:10s} score={score:.1f}")
