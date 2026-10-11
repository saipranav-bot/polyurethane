"""Final candidate prioritisation for the in-silico PU-enzyme pipeline.

Reads the ranked candidates produced by ``06_parse_hmmer_rank_candidates.py``
and combines domain priority with HMMER significance into a single score.
"""

import csv
import math
import os

RANKED = "04_hmmer/hmmer_candidates_ranked.csv"
OUT = "04_hmmer/final_candidate_ranking.csv"


def main():
    if not os.path.exists(RANKED):
        raise FileNotFoundError(
            f"{RANKED} not found - run "
            "scripts/analysis/06_parse_hmmer_rank_candidates.py first"
        )

    rows = []
    with open(RANKED) as handle:
        for row in csv.DictReader(handle):
            evalue = float(row["best_evalue"])
            weight = float(row["domain_weight"])
            # Composite score: domain priority boosted by hit significance.
            significance = -math.log10(evalue + 1e-300)
            score = weight * (1.0 + significance / 10.0)
            rows.append(
                {
                    "protein_id": row["protein_id"],
                    "best_domain": row["best_domain"],
                    "best_evalue": evalue,
                    "domain_weight": weight,
                    "score": score,
                }
            )

    rows.sort(key=lambda r: r["score"], reverse=True)

    with open(OUT, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "protein_id",
                "best_domain",
                "best_evalue",
                "domain_weight",
                "score",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print("TOP ENZYME CANDIDATES\n")
    print(f"Written: {OUT}\n")
    for i, row in enumerate(rows[:20], 1):
        print(
            f"{i:02d}. {row['protein_id']:20s} {row['best_domain']:20s} "
            f"E={row['best_evalue']:.2e} score={row['score']:.2f}"
        )


if __name__ == "__main__":
    main()
