from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import Chem
import pandas as pd

data=1
# 从Excel中读取SMILES列表和右侧格子
if data==1:
    df = pd.read_excel(r'D:\SynologyDrive\Project\2-NRP_substrate\data\COCONUT4MetFrag.xlsx', sheet_name='20230303', usecols=[0,1,2,3])
    smiles_list = df.iloc[:, 0].tolist()
else:
    df = pd.read_excel(r'D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_output3.xlsx', sheet_name='Can', usecols=[0,1,2,3])
    smiles_list = df.iloc[:, 1].tolist()



for i, smiles in enumerate(smiles_list):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        df.iloc[i, 2] = 0
    else:
        df.iloc[i, 2] = 1
        df.iloc[i, 3] = Chem.MolToSmiles(mol, isomericSmiles=False)

# 将结果写入Excel文件
if data==1:
    df.to_excel(r'D:\课题组\zhiyuan_Lab\10-Database_resource\Program\COCONUT4MetFrag_Canonical.xlsx', index=False)
else:
    df.to_excel(r'D:\课题组\zhiyuan_Lab\10-Database_resource\Program\Sid_structure_unique_Canonical.xlsx', index=False)