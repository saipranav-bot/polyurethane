#!/usr/bin/env python3

import os
import sys
import pandas as pd
from Bio import SeqIO


HMMER="04_hmmer/hmmer_tblout.txt"
FASTA="03_secretome/secreted_candidates.faa"

OUT="04_hmmer"

os.makedirs(OUT,exist_ok=True)


print("="*60)
print("STEP 6: HMMER CANDIDATE RANKING")
print("="*60)


# load sequences

records={
r.id:r
for r in SeqIO.parse(FASTA,"fasta")
}


# parse hmmer

rows=[]

with open(HMMER) as f:

    for line in f:

        if line.startswith("#"):
            continue

        x=line.split()

        if len(x)<10:
            continue

        rows.append(
        {
        "domain":x[0],
        "accession":x[1],
        "protein":x[2],
        "evalue":float(x[4]),
        "score":float(x[5])
        })


df=pd.DataFrame(rows)


print(
"Total HMMER hits:",
len(df)
)


df=df[df.evalue < 1e-5]


print(
"Significant hits:",
len(df)
)



# keep useful enzyme families

keywords=[
"hydrolase",
"esterase",
"lipase",
"cutinase",
"abhydrolase",
"carboxyl",
"feruloyl",
"COesterase",
"BD-FAE",
"AA9",
"CBM"
]


candidate=df[
df.domain.str.contains(
"|".join(keywords),
case=False,
na=False
)
]


print(
"Enzyme-like hits:",
len(candidate)
)



# best hit per protein

best=(
candidate
.sort_values("evalue")
.groupby("protein")
.first()
.reset_index()
)


out=[]


for _,r in best.iterrows():

    pid=r.protein


    if pid not in records:
        continue


    seq=records[pid]


    out.append(
    {
    "protein_id":pid,
    "length":len(seq.seq),
    "domain":r.domain,
    "evalue":r.evalue,
    "sequence":str(seq.seq)
    })


result=pd.DataFrame(out)


result=result.sort_values("evalue")


result.to_csv(
"04_hmmer/ranked_candidates.csv",
index=False
)


# top 10 fasta

SeqIO.write(
[
records[x]
for x in result.head(10).protein_id
],
"04_hmmer/top10_candidates.faa",
"fasta"
)


# top candidate

SeqIO.write(
[
records[result.iloc[0].protein_id]
],
"04_hmmer/TOP_CANDIDATE_FOR_ALPHAFOLD.faa",
"fasta"
)



print()
print("TOP CANDIDATES")
print(
result.head(10)[
[
"protein_id",
"length",
"domain",
"evalue"
]
]
.to_string(index=False)
)


print()
print("FILES CREATED:")
print("04_hmmer/ranked_candidates.csv")
print("04_hmmer/top10_candidates.faa")
print("04_hmmer/TOP_CANDIDATE_FOR_ALPHAFOLD.faa")

print()
print("STEP 6 COMPLETE")
