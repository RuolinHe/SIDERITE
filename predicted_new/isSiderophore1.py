from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

# Load reference SMILES from the reference file
with open('reference.txt', 'r') as file:
    reference_smiles = [line.strip() for line in file]

# Load ligand SMILES from the ligand file
with open('Ligand02.txt', 'r') as file: # 去掉了酚，并且增加3,4-二羟基苯
    ligand_smiles = [line.strip() for line in file]

with open('Ligand_N0.txt', 'r') as file: # 去掉了酚
    ligand_N_smiles = [line.strip() for line in file]

# Load input SMILES from the input file (only the first column)
with open('COCONUT_r.txt', 'r') as file:
    input_lines = file.readlines()
    input_smiles = [line.strip().split()[0] for line in input_lines]

# Set the Tanimoto similarity threshold and create an empty list to store matching SMILES
threshold = 0
output_smiles = []
output_similarities = []

# Convert ligand SMILES to molecules
ligand_mols = [Chem.MolFromSmiles(ligand) for ligand in ligand_smiles]

ligand_N_mols = [Chem.MolFromSmiles(ligand) for ligand in ligand_N_smiles]
# Convert reference SMILES to molecules and generate fingerprints
reference_mols = [Chem.MolFromSmiles(ref) for ref in reference_smiles]
reference_fps = [AllChem.GetMorganFingerprint(mol, 2) for mol in reference_mols]

# Iterate over each input SMILES
for input_smile in input_smiles:
    input_mol = Chem.MolFromSmiles(input_smile)
    # print(input_mol)
    # Check if the input SMILES contains any of the ligand structures
    contains_ligand = any(input_mol.HasSubstructMatch(ligand_mol) for ligand_mol in ligand_mols)
    contains_ligand_N = any(input_mol.HasSubstructMatch(ligand_mol) for ligand_mol in ligand_N_mols)
    # If the input SMILES contains ligand, proceed to similarity comparison
    if contains_ligand and not contains_ligand_N:
        
        # Generate fingerprint for the input molecule
        input_fp = AllChem.GetMorganFingerprint(input_mol, 2)
        similarity_list=[]
        # Iterate over each reference molecule and its fingerprint
        for reference_mol, reference_fp in zip(reference_mols, reference_fps):
            # Calculate the Tanimoto similarity with the reference fingerprint
            similarity_list.append(DataStructs.TanimotoSimilarity(input_fp, reference_fp))
        output_smiles.append(input_smile)
        output_similarities.append(max(similarity_list))
    else:
        # If the input SMILES doesn't contain ligand, skip similarity comparison
        continue

# Write the selected SMILES and similarities to the output file
with open('COCONUT-01.txt', 'w') as file:
    for smiles, similarity in zip(output_smiles, output_similarities):
        file.write(f"{smiles}\t{similarity}\n")
