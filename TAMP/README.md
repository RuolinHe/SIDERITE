# Change SMILES into Canonical SMILES
```
python CheckSMILES2.py
```
The content of input file COCONUT4MetFrag.xlsx is the copy of COCONUT_DB.smi.

COCONUT_DB.smi is download from [COCONUT database](https://coconut.naturalproducts.net/download) in 20230303.

The input file Sid_structure_output3.xlsx is generated in this study.

The output is COCONUT4MetFrag_Canonical.xlsx and Sid_structure_unique_Canonical.xlsx.

# Prepare input files for TAMP
These two files were used in ```tmap_code.m``` to generate COCONUT.csv and Sid_tmap20230306.csv.
The later were used to generate TAMP figure.

# TAMP script
Script in [SIDERITE](https://github.com/RuolinHe/SIDERITE/tree/main/TAMP/SIDERITE) folder was used to generate raw Figure 2:
![image](https://github.com/RuolinHe/SIDERITE/assets/76482251/c4aa52f8-b33d-4cfa-9e5f-d383303ff9a5)

Script in [COCONUT](https://github.com/RuolinHe/SIDERITE/tree/main/TAMP/COCONUT) folder was used to generate raw Figure S3:
![image](https://github.com/RuolinHe/SIDERITE/assets/76482251/7b06ae82-c6d9-4706-9e79-8bfcdcb8bb99)
