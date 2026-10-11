"""Extract enzyme-related Pfam hits from the hmmscan tblout file.

Reads the real hmmscan output produced by ``05_run_hmmer.sh`` and keeps
hits whose Pfam domain description matches hydrolytic-enzyme keywords.
"""

import csv
import os

TBL = "04_hmmer/hmmer_tblout.txt"
OUT = "04_hmmer/hmmer_enzyme_hits.csv"

keywords = [
    "hydrolase",
    "esterase",
    "lipase",
    "cutinase",
    "carboxylesterase",
    "alpha/beta hydrolase",
    "serine",
]


def main():
    if not os.path.exists(TBL):
        raise FileNotFoundError(
            f"{TBL} not found - run scripts/analysis/05_run_hmmer.sh first"
        )

    hits = []
    with open(TBL) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) < 5:
                continue
            target, query, evalue = cols[0], cols[2], cols[4]
            if any(k in target.lower() for k in keywords):
                hits.append(
                    {"protein_id": query, "domain": target, "evalue": evalue}
                )

    with open(OUT, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["protein_id", "domain", "evalue"])
        writer.writeheader()
        writer.writerows(hits)

    print(f"Enzyme-related Pfam hits: {len(hits)}")
    print(f"Written: {OUT}\n")
    for hit in hits[:20]:
        print(f"  {hit['protein_id']:20s} {hit['domain']} E={hit['evalue']}")


if __name__ == "__main__":
    main()
