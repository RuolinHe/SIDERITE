from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import os
from openpyxl import load_workbook
import sys
import getopt
import json5

argv = sys.argv[1:]
outputdic = {}

input_SMILES = Similarity_type = Threshold = ''

try:
    opts, args = getopt.getopt(argv, "h i S T", [
        "i=", "S=", "T="])
except getopt.GetoptError:
    outputdic["Result data"] = False
    outputdic["Reason data"] = "Get opt error"
    outputjson = json5.dumps(outputdic)
    print(outputjson)
    sys.exit()
for opt, arg in opts:
    if opt == '-h':
        print('FindSimilarSiderSorted.py --i <input of SMILES string> --S <Similarity type: Tanimoto, Dice ot Tversky> --T <Threshold of similarity>')
        sys.exit()
    elif opt in ("--i"):
        input_SMILES = arg
    elif opt in ("--S"):
        Similarity_type = arg
    elif opt in ("--T"):
        Threshold = float(arg)

def calculate_similarity_matrix(smiles_list1,smiles_list2, similarity_func, *args):
    """
    计算相似性矩阵
    :param smiles_list: SMILES列表
    :param similarity_func: 相似性计算函数，可以是TanimotoSimilarity、DiceSimilarity或TverskySimilarity
    :param args: 相似性计算函数的额外参数
    :return: 相似性矩阵
    """
    fps1 = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2) for smiles in smiles_list1]
    fps2 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles_list2), 2)
    similarity_matrix = np.zeros(len(smiles_list1))
    for i in range(len(smiles_list1)):
        similarity_matrix[i] = similarity_func(fps1[i], fps2, *args)
    return similarity_matrix

# 从Excel中读取SMILES列表，跳过第一行表头
df1 = pd.read_excel(r'./data/Sider_database_web20230605.xlsx', sheet_name='Unique1', usecols=[21])
reference_smiles_list = df1.iloc[:, 0].tolist()
df1 = pd.read_excel(r'./data/Sider_database_web20230605.xlsx', sheet_name='Unique1', usecols=[0])
SID_list = df1.iloc[:, 0].tolist()
df1 = pd.read_excel(r'./data/Sider_database_web20230605.xlsx', sheet_name='Unique1', usecols=[1])
Name_list = df1.iloc[:, 0].tolist()
with open('./data/Ligand01.txt', 'r') as file: # 去掉了酚
    ligand_smiles = [line.strip() for line in file]
with open('./data/Ligand_N0.txt', 'r') as file: # 去掉了酚
    ligand_N_smiles = [line.strip() for line in file]
ligand_mols = [Chem.MolFromSmiles(ligand) for ligand in ligand_smiles]
ligand_N_mols = [Chem.MolFromSmiles(ligand) for ligand in ligand_N_smiles]
input_mol = Chem.MolFromSmiles(input_SMILES)
if input_mol is None:
    outputdic['Result'] = []
    outputdic["Result data"] = False
    outputdic["Reason data"] = "SMILES input is incorrect"
    outputdic['New'] = 0
else:
    contains_ligand = any(input_mol.HasSubstructMatch(ligand_mol) for ligand_mol in ligand_mols)
    contains_ligand_N = any(input_mol.HasSubstructMatch(ligand_mol) for ligand_mol in ligand_N_mols)
    # If the input SMILES contains ligand, proceed to similarity comparison
    if contains_ligand and not contains_ligand_N:
        Potential = 1
    else:
        Potential = 0

    if Similarity_type == 'Tanimoto':
        Similarity_matrix = calculate_similarity_matrix(reference_smiles_list, input_SMILES, DataStructs.TanimotoSimilarity)
    elif Similarity_type == 'Dice':
        Similarity_matrix = calculate_similarity_matrix(reference_smiles_list, input_SMILES, DataStructs.DiceSimilarity)
    elif Similarity_type == 'Tversky':
        Similarity_matrix = calculate_similarity_matrix(reference_smiles_list, input_SMILES, DataStructs.TverskySimilarity)

    Similarity_matrix_filtered = [x for x in Similarity_matrix if x > Threshold]
    Name_list_filtered = [Name_list[i] for i in range(len(Similarity_matrix)) if Similarity_matrix[i] > Threshold]
    SID_list_filtered = [SID_list[i] for i in range(len(Similarity_matrix)) if Similarity_matrix[i] > Threshold]

    Similarity_matrix_filtered_sorted = sorted(Similarity_matrix_filtered,reverse = True)
    Name_list_filtered_sorted = sorted(Name_list_filtered, key=lambda x: Similarity_matrix_filtered[Name_list_filtered.index(x)], reverse=True)
    SID_list_filtered_sorted = sorted(SID_list_filtered, key=lambda x: Similarity_matrix_filtered[SID_list_filtered.index(x)], reverse=True)
    SID_list_filtered_sorted_padded = ["SID" + str(num).zfill(5) for num in SID_list_filtered_sorted]

    Result = []
    for sid, name, similarity in zip(SID_list_filtered_sorted_padded, Name_list_filtered_sorted, Similarity_matrix_filtered_sorted):
        entry = {"SID": sid, "Name": name, "Similarity": similarity}
        Result.append(entry)

    outputdic['Result'] = Result
    outputdic["Result data"] = True
    outputdic["Reason data"] = ''
    outputdic['New'] = 0
    if Potential == 1:
        if len(Similarity_matrix_filtered_sorted) == 0:
            outputdic['New'] = 1
        elif max(Similarity_matrix_filtered_sorted) != 1:
            outputdic['New'] = 1

outputjson = json5.dumps(outputdic)
print(outputjson)
