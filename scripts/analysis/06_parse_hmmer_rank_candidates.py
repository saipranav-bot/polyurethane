"""Parse hmmscan output and rank secreted candidates by Pfam domain.

Reads the real ``hmmscan --tblout`` file produced by
``scripts/analysis/05_run_hmmer.sh`` and ranks each protein by the best
domain hit it carries. No surrogate/length-based scoring is used.
"""

from collections import defaultdict
import csv
import os

TBL = "04_hmmer/hmmer_tblout.txt"
OUT = "04_hmmer/hmmer_candidates_ranked.csv"

# Pfam domain families indicative of serine-hydrolase / esterase activity,
# with a relative priority weight (higher = more relevant to polyurethane
# (ester/urethane) bond cleavage).
priority_domains = {
    "cutinase": 5.0,
    "abhydrolase": 4.0,
    "lipase": 4.0,
    "esterase": 4.0,
    "phospholip": 3.0,
    "glyco_hydro": 1.0,
}


def domain_weight(target):
    tl = target.lower()
    for domain, weight in priority_domains.items():
        if domain in tl:
            return weight
    return 0.0


def main():
    if not os.path.exists(TBL):
        raise FileNotFoundError(
            f"{TBL} not found - run scripts/analysis/05_run_hmmer.sh first"
        )

    best_evalue = {}
    best_domain = {}

    with open(TBL) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 5:
                continue
            target = cols[0]
            query = cols[2]
            try:
                evalue = float(cols[4])
            except ValueError:
                continue
            if query not in best_evalue or evalue < best_evalue[query]:
                best_evalue[query] = evalue
                best_domain[query] = target

    rows = []
    for query, evalue in best_evalue.items():
        domain = best_domain[query]
        rows.append(
            {
                "protein_id": query,
                "best_domain": domain,
                "best_evalue": evalue,
                "domain_weight": domain_weight(domain),
            }
        )

    rows.sort(key=lambda r: (-r["domain_weight"], r["best_evalue"]))

    with open(OUT, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["protein_id", "best_domain", "best_evalue", "domain_weight"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Parsed {len(rows)} Pfam-annotated proteins")
    print(f"Written: {OUT}\n")
    print("Top candidates:")
    for row in rows[:20]:
        print(
            f"  {row['protein_id']:20s} {row['best_domain']:20s} "
            f"E={row['best_evalue']:.2e} weight={row['domain_weight']}"
        )


if __name__ == "__main__":
    main()
