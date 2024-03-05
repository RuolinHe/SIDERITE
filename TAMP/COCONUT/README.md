Because the data is too large, we save it in the [DOI 10.5281/zenodo.10369626](https://zenodo.org/doi/10.5281/zenodo.10369626).

The script and data used in SIDERITE paper to generate TMAP picture (Figure S3 in supplementart material).

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
 

Usage: ```python plot_COCONUT.py```

Then it will use  COCONUT.csv to generate  index.html and index.js. Open index.html to see result.
