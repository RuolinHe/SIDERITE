# Large files archive
COCONUT_r.txt and Tanimoto_COCONUT_SIDERITE.xlsx can be found in [DOI 10.5281/zenodo.10369626](https://zenodo.org/doi/10.5281/zenodo.10369626)

# Prepare clean COCONUT data
The input file COCONUT_r.txt is COCONUT database excluded known siderophores from SIDERITE database.

# Find potential siderophore in COCONUT
Generate COCONUT-01.txt by
```
python isSiderophore1.py
```

# Calculate Tanimoto similarity between potential siderophores in COCONUT and known siderophores from SIDERITE database
The content of COCONUT-01.txt was copied to Summary.xlsx.
```
python SMILES2distance.py
```
The output is Tanimoto_COCONUT_SIDERITE.xlsx

# draw figure
Use matlab script ```Clustering_coconut.m```
