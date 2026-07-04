import numpy as np
from Bio.PDB import MMCIFParser

parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("g13490", "results/structure/g13490_best_model.cif")
chain = list(structure[0].get_chains())[0]
res_by_num = {r.get_id()[1]: r for r in chain.get_residues() if r.get_id()[0] == " "}

old_center = np.array([-7.700, 5.399, -3.886])

triad_coords = []
for resnum in [227, 462, 226]:
    res = res_by_num[resnum]
    if "CA" in res:
        v = res["CA"].get_vector()
        triad_coords.append([v[0], v[1], v[2]])

new_center = np.array(triad_coords).mean(axis=0)
shift = np.linalg.norm(new_center - old_center)

print(f"OLD grid center (Ser227-Asp421-His462): {old_center}")
print(f"NEW grid center (Ser227-Glu226-His462): {new_center}")
print(f"Shift distance: {shift:.2f} A")
print()
if shift < 5:
    print("Shift is SMALL — existing docking results remain valid as-is.")
else:
    print("Shift is LARGE — docking should be re-run with corrected center.")
