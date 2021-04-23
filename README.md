# CASP14_high_accuracy

This repository keeps the scripts used for high-accuracy assessment by the Lupas team in CASP14.

## Content

There are 2 folders:

 - Code: Contains the python scripts used to compile and analyse data from the Prediction Center, as well as the jupyter notebook for data analysis and a script for the calculation of the DipDiff score given two structure of the same protein
 - Data: Contains the data compiled as well the scores calculated for individual targets, for CASP14, CASP13 and CASP12

## The main different scripts

***1. ```analyse_casp.py```***

The script that collects targets and models from the Prediction Center, finds templates from HHsearch results provided by the Prediction Center, and computes addicional scores (e.g. DipDiff). This script generates the files in the "out_tables" folder in the the "Data/CASPx" (where x is a CASP identifier) folder.

External tools required to run it:

 - Muscle (https://bio.tools/MUSCLE)
 - LGA (http://proteinmodel.org/AS2TS/LGA/lga.html)
 - DSSP (https://bio.tools/dssp)
 - pdb-tools (http://www.bonvinlab.org/pdb-tools/)
 - CCP4 with ARP/wARP or:
   - DipCheck (https://bio.tools/dipcheck)
   - MolProbity (https://bio.tools/molprobity)
   - phenix.pdbtools (https://www.phenix-online.org/documentation/reference/pdbtools.html)
   - sculptor (https://bio.tools/sculptor)

Special python modules:
 - pandas
 - seaborn
 - matplotlib
 - BioPython
 - scipy
 - numpy
 - pdbe api

Please edit the paths to these tools in the header of the script.

ATTENTION: this script requires an assessor username and password. Replace USERNAME and PASSWORD by your assessor information in the following line: ```156   password_mgr.add_password(None, top_level_url, 'USERNAME', 'PASSWORD')```

How to run: 

```
ulimit -s unlimited; python3 analyse_casp.py -casp 14 -tmp /tmp/
```

The only input needed to run this script is a casp identifier (e.g. CASP14, casp14 or 14). Only CASP12, 13 or 14 are possible. For parallelization, it can also take the ```n_threads``` parameter. It will then download all targets and the models and all computations will be carried out in the tmp folder as it generates multiple files while runnning. Also, you have to set the memory limit to unlimited due to LGA.


***2. ```assess_casp.py```***

This script shall be run after ```analyse_casp.py``` to compute Z- and ranking scores, and collect method information from abstracts. It will generate a unique table file that can then be processed and analysed with the jupyter notebook. 

External tools required to run it:

 - pdftotext (https://linux.die.net/man/1/pdftotext) (it is assumed to be in the path)

Special python modules:
 - pandas
 - matplotlib
 - numpy

How to run: 

```
python3 assess_casp.py -casp 14 -tmp /tmp/
```

This will take all the data computed for the input casp with ```analyse_casp.py``` and saved in the tmp folder, and compute Z-scores and generate a unique data table.


***3. ```DipDiff.py```***

This script takes 2 pdb files as input and computes the DipDiff. If mode "template" is passed, it will find residue correspondencies with LGA.

External tools required to run it:

 - LGA (http://proteinmodel.org/AS2TS/LGA/lga.html)
 - CCP4 with ARP/wARP or:
   - DipCheck (https://bio.tools/dipcheck)

How to run: 

```
python3 dipdiff.py -target target.pdb -model model.pdb
```

If the model is a template, add ```-mode template``` and run instead:

```
ulimit -s unlimited; python3 dipdiff.py -target target.pdb -model model.pdb -mode template
```

