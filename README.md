# From lipid annotations to molecular structures

Assigning a plausible molecular structure (SMILES) to each lipid in a published LC-MS lipidomics dataset, then checking every structure against the reported masses and validating it with **RDKit**.

## Background

Untargeted lipidomics usually reports lipids at the **sum-composition** level, e.g. `PC 34:1`: the class plus the total number of carbons and double bonds in the fatty acyl chains. Cheminformatics tools need explicit structures. This project turns each annotation into one *plausible* molecular species (e.g. `PC 16:0/18:1`) with a full SMILES, and checks it independently.

**Data:** Supplementary Table 1 of Afshinnia *et al.*, *JCI Insight* 2019 (serum lipidomics in 92 American Indians with type 2 diabetes; open access). The table was transcribed from the supplement PDF; see [`data/SOURCE.txt`](data/SOURCE.txt).

## Method

1. **`scripts/build_structures.py`** (rule-based construction)
   - Splits the sum composition into common mammalian fatty acyl chains, using abundance-weighted priors (e.g. 16:0, 18:1, 18:2, 20:4) and very-long-chain priors for sphingomyelins.
   - Places cis double bonds at their usual positions and assembles the SMILES from class-specific headgroup templates. It covers 17 classes: TG, DG, MG, PC, PE, PS, PI, PG, PA, LPC, LPE, SM, CerP, CE, FA, and plasmalogen PE/LPC.
   - Computes formula and monoisotopic mass with a small SMILES parser that doesn't use RDKit, and compares the expected adduct *m/z* with the reported *m/z* (ppm error).
2. **`scripts/validate_rdkit.py`** (independent check)
   - Parses and sanitises every SMILES with RDKit (valence, rings), compares RDKit's formula with the one from step 1, and reports undefined stereocentres.

## Results

| Check | Result |
|---|---|
| Detected features in source table | 398 |
| Unique lipid annotations | 298 (17 classes; most common TG 63, PC 42, SM 38) |
| Chain sum equals annotation | 298 / 298 |
| Valid SMILES in RDKit | 298 / 298 |
| Formula agrees with RDKit | 298 / 298 |
| Mass error ≤ 5 ppm vs reported *m/z* | 296 / 298 (median 0.7 ppm) |

**Flagged:** `LPC 16:0` and `pPC 18:0` have one reported adduct about 13 Da away from the expected *m/z*. This points to an adduct or annotation inconsistency in the source table, not a structure error. Both are kept in the output with their ppm error.

**Limitations:** each structure is *one plausible* species consistent with the annotation, not the uniquely identified one. MS1 data can't resolve chain composition, double-bond positions or sn-position. Glycerol stereocentres in TG/DG are left undefined (213 structures have at least one undefined stereocentre).

## Repository layout

```
data/      source table (transcribed) + citation
scripts/   build_structures.py, validate_rdkit.py
results/   lipid_structures.csv, lipid_structures_rdkit_validated.csv
```

## Reproduce

```bash
pip install rdkit
python scripts/build_structures.py data/Afshinnia2019_SupplTable1_lipid_species.csv results/lipid_structures.csv
cd results && python ../scripts/validate_rdkit.py
```

`build_structures.py` needs only the Python standard library.

## Reference

Afshinnia F. *et al.* Increased lipogenesis and impaired beta-oxidation predict type 2 diabetic kidney disease progression in American Indians. *JCI Insight* 4(21):e130317 (2019). doi:[10.1172/jci.insight.130317](https://doi.org/10.1172/jci.insight.130317)
