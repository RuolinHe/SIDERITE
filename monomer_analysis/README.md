Monomer annotation of siderophores was applied by [rBAN](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0335-x) with custom monomer database based on original Norine database.

The custom monomer database is MonomersDatabase20230227.json

unique.json contains 649 siderophores with unique structure.

rBAN can be downloaded from [here](https://bitbucket.org/sib-pig/rban/downloads/).

After installing rBAN, you can run it by this command:
```
java -jar rBAN-1.0.jar -inputFile unique.json -outputFolder unique -discoveryMode -monomersDB MonomersDatabase20230227.json
```
