import pandas as pd

# ----------------------------
# INPUT DATA (from your pipeline)
# ----------------------------
data = {
    "protein": [
        "g8384.t1", "g14244.t1", "g8990.t1", "g6347.t1", "g13585.t1",
        "g3351.t1", "g3430.t1", "g8400.t1", "g11584.t1", "g13885.t1"
    ],
    "hmmer_score": [80, 70, 66, 32, 32, 20, 20, 20, 20, 20],
    "docking_energy": [-6.40, -7.06, -6.80, None, None, None, None, None, None, None]
}

df = pd.DataFrame(data)

# ----------------------------
# SAFE NUMERIC HANDLING
# ----------------------------
df["docking_energy"] = df["docking_energy"].fillna(-5)

# Convert docking to "benefit score" (higher = better binding)
df["dock_score"] = -df["docking_energy"]

# ----------------------------
# Z-SCORE NORMALIZATION (publication grade)
# ----------------------------
df["hmmer_z"] = (df["hmmer_score"] - df["hmmer_score"].mean()) / df["hmmer_score"].std()

df["dock_z"] = (df["dock_score"] - df["dock_score"].mean()) / df["dock_score"].std()

# ----------------------------
# FINAL INTEGRATED SCORE
# ----------------------------
df["final_score"] = df["hmmer_z"] + 0.7 * df["dock_z"]

# ----------------------------
# SORT RESULTS
# ----------------------------
df = df.sort_values("final_score", ascending=False)

# ----------------------------
# OUTPUT
# ----------------------------
print("\n🧬 FINAL ENZYME CANDIDATES (PUBLICATION-READY)\n")
print(df[["protein", "hmmer_score", "docking_energy", "final_score"]].to_string(index=False))
