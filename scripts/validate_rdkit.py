"""Step 2: validate every SMILES in lipid_structures.csv with RDKit.

Run in the lipidomics_data folder (Jupyter: %run validate_rdkit.py, or: python validate_rdkit.py).
Writes lipid_structures_rdkit_validated.csv with the RDKit results added.
"""
import csv
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

IN, OUT = "lipid_structures.csv", "lipid_structures_rdkit_validated.csv"
rows = list(csv.DictReader(open(IN, encoding="utf-8")))

n_ok = 0
for r in rows:
    mol = Chem.MolFromSmiles(r["smiles"])          # parses AND sanitizes (valence, rings, ...)
    if mol is None:
        r.update(rdkit_valid="NO", rdkit_formula="", rdkit_exact_mass="", formula_match="",
                 rdkit_canonical_smiles="", undefined_stereocentres="")
        continue
    f = CalcMolFormula(mol)
    # RDKit writes charged zwitterions as e.g. C42H82NO8P (net charge 0), same as our formula
    undefined = [c for c, lab in Chem.FindMolChiralCenters(mol, includeUnassigned=True) if lab == "?"]
    r.update(rdkit_valid="YES",
             rdkit_formula=f,
             rdkit_exact_mass=round(Descriptors.ExactMolWt(mol), 4),
             formula_match="YES" if f == r["formula"] else "NO",
             rdkit_canonical_smiles=Chem.MolToSmiles(mol),
             undefined_stereocentres=len(undefined))
    n_ok += 1

with open(OUT, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

print(f"RDKit-valid SMILES: {n_ok}/{len(rows)}")
print("Formula mismatches:", sum(r["formula_match"] == "NO" for r in rows))
print("Invalid:", [r["lipid_annotation"] for r in rows if r["rdkit_valid"] == "NO"])
print("Saved", OUT)
