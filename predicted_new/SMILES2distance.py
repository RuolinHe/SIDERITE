from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np


def calculate_similarity_matrix(smiles_list, similarity_func, *args):
    """
    计算相似性矩阵
    :param smiles_list: SMILES列表
    :param similarity_func: 相似性计算函数，可以是TanimotoSimilarity、DiceSimilarity或TverskySimilarity
    :param args: 相似性计算函数的额外参数
    :return: 相似性矩阵
    """
    fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), 2) for smiles in smiles_list]
    similarity_matrix = np.zeros((len(smiles_list), len(smiles_list)))
    for i in range(len(smiles_list)):
        for j in range(i, len(smiles_list)):
            similarity_matrix[i][j] = similarity_matrix[j][i] = similarity_func(fps[i], fps[j], *args)
    return similarity_matrix


def save_similarity_matrix_to_excel(smiles_list, similarity_matrix, sheet_name, file_path):
    """
    将相似性矩阵保存到Excel文件中
    :param smiles_list: SMILES列表
    :param similarity_matrix: 相似性矩阵
    :param sheet_name: Sheet名称
    :param file_path: Excel文件路径
    """
    df = pd.DataFrame(similarity_matrix, index=smiles_list, columns=smiles_list)
    df.to_excel(file_path, sheet_name=sheet_name)


# 从Excel中读取SMILES列表，跳过第一行表头
df = pd.read_excel(r'../TAMP/summary.xlsx', sheet_name='Calculate_similarity', usecols=[1])
smiles_list = df.iloc[:, 0].tolist()
tanimoto_matrix = calculate_similarity_matrix(smiles_list, DataStructs.TanimotoSimilarity)
print('Saving data')
save_similarity_matrix_to_excel(smiles_list, tanimoto_matrix, 'Tanimoto', 'Tanimoto_COCONUT_SIDERITE.xlsx')
