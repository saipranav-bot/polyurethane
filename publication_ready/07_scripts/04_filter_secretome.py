#!/usr/bin/env python3

"""
04_filter_secretome.py

SignalP 6 + TMHMM secretome filtering
"""

import os
import sys
import re
import pandas as pd
from Bio import SeqIO


# ---------------- CONFIG ----------------

SIGNALP_FILE = "03_secretome/signalp6_output.txt"
TMHMM_FILE   = "03_secretome/tmhmm_output.txt"

SIZE_FASTA = "02_size_filter/candidates_sizefiltered.faa"

OUT_DIR = "03_secretome"

OUT_FASTA = f"{OUT_DIR}/secreted_candidates.faa"
OUT_TABLE = f"{OUT_DIR}/secreted_candidates_table.csv"
OUT_STATS = f"{OUT_DIR}/secretome_stats.txt"


os.makedirs(OUT_DIR, exist_ok=True)


print("="*60)
print("STEP 4: SECRETOME FILTER")
print("="*60)


# ---------------- CHECK ----------------

for f in [SIGNALP_FILE, TMHMM_FILE, SIZE_FASTA]:

    if not os.path.exists(f):
        print("Missing:", f)
        sys.exit()


# ---------------- LOAD PROTEINS ----------------

records = {
    r.id:r
    for r in SeqIO.parse(SIZE_FASTA,"fasta")
}


print("Input proteins:", len(records))


# ---------------- SIGNALP ----------------

print("\nReading SignalP...")


signalp_ids=set()


with open(SIGNALP_FILE) as f:

    for line in f:

        if line.startswith("#"):
            continue

        parts=line.strip().split()


        if len(parts) < 2:
            continue


        protein_id = parts[0]

        prediction = parts[1]


        if prediction.upper() == "SP":

            signalp_ids.add(protein_id)



print(
    "SignalP positive proteins:",
    len(signalp_ids)
)



# ---------------- TMHMM ----------------

print("\nReading TMHMM...")


tm_ids=set()


with open(TMHMM_FILE) as f:

    for line in f:


        if line.startswith("#"):
            continue


        if line.strip()=="":
            continue


        m=re.search(
            r"PredHel=(\d+)",
            line
        )


        if m:

            if int(m.group(1)) > 0:

                tm_ids.add(
                    line.split()[0]
                )



print(
    "TM proteins removed:",
    len(tm_ids)
)



# ---------------- FILTER ----------------

final=[]

table=[]


for sid in signalp_ids:


    match=None


    if sid in records:

        match=records[sid]


    else:


        for fid,r in records.items():

            if (
                sid in fid
                or fid in sid
                or fid.split(".")[0]
                ==
                sid.split(".")[0]
            ):

                match=r
                break



    if match is None:
        continue



    if match.id in tm_ids:
        continue



    final.append(match)


    table.append(
        {
            "protein_id": match.id,
            "length_aa": len(match.seq),
            "approx_kDa": round(len(match.seq)*110/1000,2),
            "SignalP": "SP",
            "TMHMM": "No_TM"
        }
    )



print(
    "\nSecreted candidates:",
    len(final)
)



if len(final)==0:

    print("No secreted candidates found")
    sys.exit()



# ---------------- OUTPUT ----------------


SeqIO.write(
    final,
    OUT_FASTA,
    "fasta"
)



pd.DataFrame(table).to_csv(
    OUT_TABLE,
    index=False
)



with open(OUT_STATS,"w") as f:

    f.write(
f"""
SECRETOME FILTER RESULTS

Input proteins:
{len(records)}

SignalP positive:
{len(signalp_ids)}

TM removed:
{len(tm_ids)}

Final secreted candidates:
{len(final)}
"""
)



print("\nFILES CREATED:")

print(OUT_FASTA)
print(OUT_TABLE)
print(OUT_STATS)


print("\nSTEP 4 COMPLETE")

print("\nNEXT:")
print("bash 05_run_hmmer.sh")
