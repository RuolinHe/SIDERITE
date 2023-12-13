The script and data used in SIDERITE paper to generate TMAP picture (Figure 2 and Indexing Interface).

Requirement: tmap

You can install tmap by conda.

```
conda create -n tmap python=3.7

conda activate tmap

conda install -c tmap tmap

pip install faerun

pip install matplotlib

conda install scipy

conda install -c rdkit rdkit

conda install -c conda-forge mhfp
```
 

Usage: ```python plot_sid.py```

Then it will use Sid_tmap20230306.csv to generate  index.html and index.js. Open index.html to see result.
