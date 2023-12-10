# import pickle
# import numpy as np
import tmap as tm
import pandas as pd
import scipy.stats as ss
from rdkit.Chem import AllChem
from mhfp.encoder import MHFPEncoder
from faerun import Faerun
from collections import Counter
from matplotlib.colors import ListedColormap
# from matplotlib import pyplot as plt


def main():
    """ Main funciton """
    df = pd.read_csv("Sid_tmap20230306.csv", sep="\t")

    enc = MHFPEncoder(1024)
    lf = tm.LSHForest(1024, 64)

    fps = []
    labels = []
    for i, row in df.iterrows():
        if i != 0 and i % 1000 == 0:
            print(100 * i / len(df))
        mol = AllChem.MolFromSmiles(row["Canonical SMILES"])
        fps.append(tm.VectorUint(enc.encode_mol(mol)))
        label = (
            row["Canonical SMILES"]
            + '__<a target="_blank" href="http://siderite.bdainformatics.org/database?conditions=%5B%5D&query=&sid='
            + row["Siderophore ID"]
            + '">'
            + row["Siderophore ID"]
            + "</a>"
            + "__"
            + row["Siderophore Class"]
            + "__"
            + row["Siderophore name"]
            + "__"
            + row["Siderophore other name"]
            + "__"
            + str(round(row["Molecular Weight"],2))
            + "__"
            + row["Kingdom"]
            + "__"
            + row["Biosynthetic Type"]
            + "__"
            + str(row["Theoretical denticity"])
            + "__"
            + str(row["monomers Number"])
            + "__"
            + row["Precursor"]
            + "__"
            + row["Ligand Type"]
        )
        labels.append(label.replace("'", "´"))

    lf.batch_add(fps)
    lf.index()
    cfg = tm.LayoutConfiguration()
    cfg.node_size = 5 # 数值越大越不容易重叠
    # cfg.mmm_repeats = 2
    # cfg.sl_extra_scaling_steps = 5
    # cfg.k = 20
    # cfg.sl_scaling_type = tm.RelativeToAvgLength
    cfg.k = 50
    cfg.kc = 50
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)

    Siderophore_Small_Class_data=df["Siderophore Small Class"]
    Theoretical_denticity_data=df["Theoretical denticity"]
    Count_data=df["Count"]
    monomers_Number_data=df["monomers Number"]
    coverage_data=df["coverage"]
    Molecular_Weight_data=df["Molecular Weight"]
    Hydroxamate_data=df["Hydroxamate"]
    Catecholate_data=df["Catecholate"]
    Phenolate_data=df["Phenolate"]
    Carboxylate_data=df["Carboxylate"]
    Carboxylate_in_Citrate_data=df["Carboxylate in Citrate"]
    Alpha_Hydroxycarboxylate_data=df["Alpha-Hydroxycarboxylate"]
    Hydroxyphenyloxazoline_data=df["Hydroxyphenyloxazoline"]
    Hydroxyphenylthiazoline_data=df["Hydroxyphenylthiazoline"]
    Alpha_Aminocarboxylate_data=df["Alpha-Aminocarboxylate"]
    Alpha_Hydroxyimidazole_data=df["Alpha-Hydroxyimidazole"]
    Alpha_Hydroxycarboxylate_in_Citrate_data=df["Alpha-Hydroxycarboxylate in Citrate"]
    Diazeniumdiolate_data=df["Diazeniumdiolate"]
    Nitrosophenol_data=df["2-Nitrosophenol"]
    logS_data=df["Predicted logS"]
    diffusion_coefficient_data=df["Predicted diffusion coefficient"]

    Biosynthetic_Type, Biosynthetic_Data = Faerun.create_categories(df["Biosynthetic Type"])
    Sider_Class_Type, Sider_Class_Data = Faerun.create_categories(df["Siderophore Large Class"])
    Phylum_Type, Phylum_Data = Faerun.create_categories(df["Phylum"])
    Kingdom_Type, Kingdom_Data = Faerun.create_categories(df["Kingdom"])
    Precursor_Type, Precursor_Data = Faerun.create_categories(df["Precursor"])

    top_num=17
    top_Precursor = [i for i, _ in Counter(Precursor_Data).most_common(top_num-1)]

    top_Precursor_Type = []
    Precursor_map = [top_num-1] * len(Precursor_Data)
    value = 0
    for i, name in Precursor_Type:
        if i in top_Precursor:
            v = value
            top_Precursor_Type.append((v, name))
            Precursor_map[i] = v
            value += 1
    top_Precursor_Type.append((top_num-1, "Other"))
    Precursor_Data = [Precursor_map[val] for _, val in enumerate(Precursor_Data)]

    # df["Canonical SMILES"] = (
    #     df["Canonical SMILES"]
    #     + '__<a target="_blank" href="https://www.npatlas.org/joomla/index.php/explore/compounds#npaid='
    #     + df["Siderophore ID"]
    #     + '">'
    #     + df["Siderophore ID"]
    #     + "</a>"
    # )

    custom_cmap_class = ListedColormap(
        ["#dc143c", "#00bfff", "#008000","#ffff00", "#ff00ff", "#ffa500", "#66cdaa", '#7fff00',"#0000ff","#ff4500","#db7093",'#f0e68c','#1e90ff','#8b4513','#483d8b','#ee82ee', '#8b008b', '#2f4f4f', '#ff1493','#ffa07a','#00008b','#556b2f','#9acd32','#00ff7f',"#c0c0c0"],
        # ["#ff0000", "#00ff00", "#ffff00", "#ffa500", "#ff00ff", "#00ffff", '#1e90ff',"#008000","#ffffff","#0000ff",'#ee82ee','#00fa9a','#b03060','#bdb76b','#8b4513','#000080','#708090'],
        name="class",
    )
    custom_cmap_Precursor = ListedColormap(
        ["#FF0000", "#FFFF00", "#006400", "#00FF00", "#ADFF2F","#00FF7F","#FF00FF","#EE82EE","#9400D3","#9370DB","#800080","#4B0082","#0000FF","#4169E1","#D3D3D3","#00BFFF","#000000"],
        name="precursor",
    )
    custom_cmap_bio = ListedColormap(
        ["#008080", "#ffa500", "#0000ff", "#ff1493", "#00ff00"],
        name="bio",
    )
    custom_cmap_Phylum = ListedColormap(
        ["#0000FF", "#00BFFF", "#00FF00","#98FB98", "#3CB371", "#006400", "#808000", "#00FA9A","#FFA500","#FFFF00","#F0E68C","#EE82EE","#FF00FF","#9400D3","#2E8B57","#FF0000", "#8B0000","#CD5C5C","#2F4F4F" ,"#000000"],
        # ["#0000ff", "#000080", "#1e90ff", "#00ff00", "#006400", "#00ffff", "#00fa9a","#eee8aa","#ffa500","#ffff00","#000000","#708090","#c71585","#ff0000","#ff00ff","#8b4513","#deb887","#6f00e9","#f97fff","#bd5fff","#d46b00","#f4b900"],
        name="Phylum",
    )
    # f = Faerun(clear_color="#ffffff",view="front", coords=False) # 设置白色背景clear_color，默认为黑色背景
    f = Faerun(clear_color="#ffffff",view="front", coords=False, impress='<font color="#000000">This page is developed by </font> <a href="https://cqb.pku.edu.cn/zyli/info/1006/1045.htm" target="_blank">Ruolin He</a><br /><font color="#000000">from </font> <a href="https://cqb.pku.edu.cn/zyli/index.htm" target="_blank">Zhiyuan Lab of Peking University</a>') # 设置白色背景clear_color，默认为黑色背景 在左上角加上链接，参考https://tmap.gdb.tools/src/dsstox/dsstox.html
    f.add_scatter(
        "Siderophore_Information_Database",
        {
            "x": x,
            "y": y,
            "c": [
                Biosynthetic_Data,
                Kingdom_Data,
                Sider_Class_Data,
                Precursor_Data,
                Phylum_Data,
                Siderophore_Small_Class_data,
                Theoretical_denticity_data,
                Count_data,
                monomers_Number_data,
                coverage_data,
                Molecular_Weight_data,
                Hydroxamate_data,
                Catecholate_data,
                Phenolate_data,
                Carboxylate_data,
                Carboxylate_in_Citrate_data,
                Alpha_Hydroxycarboxylate_data,
                Hydroxyphenyloxazoline_data,
                Hydroxyphenylthiazoline_data,
                Alpha_Aminocarboxylate_data,
                Alpha_Hydroxyimidazole_data,
                Alpha_Hydroxycarboxylate_in_Citrate_data,
                Diazeniumdiolate_data,
                Nitrosophenol_data,
                logS_data,
                diffusion_coefficient_data,
                ],
            "labels": labels,
        },
        shader="smoothCircle",
        point_scale=10.0, # 调整点的大小
        max_point_size=100,
        legend_labels=[Biosynthetic_Type, Kingdom_Type,Sider_Class_Type, top_Precursor_Type, Phylum_Type], # 条目顺序是字母A-Z
        selected_labels=["SMILES", "Siderophore Information Database ID", "Siderophore Class", "Name", "Other name", "Molecular Weight", "Kingdom of Source", "Biosynthetic Type", "Theoretical denticity", "Number of Monomers", "Precursor", "Ligand Type"],
        categorical=[True, True, True, True, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False],
        colormap=[custom_cmap_bio,custom_cmap_bio,custom_cmap_class,custom_cmap_Precursor, custom_cmap_Phylum, "gist_rainbow","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r","winter_r"],
        series_title=[
            "Biosynthetic Type",
            "Kingdom of Source",
            "Siderophore Large Class",
            "Precursor",
            "Phylum",
            "Siderophore Small Class",
            "Theoretical denticity",
            "Count",
            "Number of Monomers",
            "Coverage of Monomers",
            "Molecular Weight",
            "Number of Hydroxamate",
            "Number of Catecholate",
            "Number of Phenolate",
            "Number of Carboxylate",
            "Number of Carboxylate in Citrate",
            "Number of Alpha-Hydroxycarboxylate",
            "Number of Hydroxyphenyloxazoline",
            "Number of Hydroxyphenylthiazoline",
            "Number of Alpha-Aminocarboxylate",
            "Number of Alpha-Hydroxyimidazole",
            "Number of Alpha-Hydroxycarboxylate in Citrate",
            "Number of Diazeniumdiolate",
            "Number of 2-Nitrosophenol",
            "Predicted logS",
            "Predicted diffusion coefficient",
        ],
        has_legend=True,
    )
    f.add_tree("Siderophore_tree", {"from": s, "to": t}, point_helper="Siderophore_Information_Database")
    f.plot(template="smiles")

if __name__ == "__main__":
    main()