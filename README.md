[![icon-small.png](https://github.com/ZzakB/MolNotator/blob/main/images/icon_small.png?raw=true)](https://github.com/ZzakB/MolNotator)

MolNotator is a Python package that predicts the actual molecules present in LC-MS/MS data.
The final data is represented in the form of actual molecular networks, representing the predicted molecules as nodes amidst the ions they generated.
The aim of the method is to help users of LCMS to pinpoint the molecules of interest in their data and avoid the effort of sorting through ions to manually find their target compound.
## Features
 
- Predicts molecule nodes in spectral networks
- Dereplicates with spectral and exact mass data
- Adds retention time and adduct filters to dereplication
- Supervised adduct search


## Documentation

**Note:** This README provides instructions for setup and using basic functions of MolNotator.
For more details, see the [paper](https://doi.org/10.1101/2021.12.21.473622).

MolNotator works within a user-defined project folder with a specific file structure. An example is given in the [examples](https://github.com/ZzakB/MolNotator/tree/main/examples) folder:
```
working_directory
|   input.py
|
|___databases
|   |   211005_MIX_LDB.mgf
|   |   211018_COLOTUS_DB.tsv
|   
|___input_files
|   |   200909_LDB_Thermo_NEG.csv
|   |   200909_LDB_Thermo_NEG.mgf
|   |   200912_LDB_Thermo_POS.csv
|   |   200912_LDB_Thermo_POS.mgf
|    
|___params
|   |   fragnotator_table.tsv
|   |   NEG_adduct_table_primary.tsv
|   |   NEG_adduct_table_secondary.tsv
|   |   POS_adduct_table_primary.tsv
|   |   POS_adduct_table_secondary.tsv
|   |   params.yaml
|   |   params_colotus.yaml
|   |   params_ldb_ions.yaml
```

The databases folder contains database files in MGF, TSV or CSV format. Two files are provided in the examples.
The mzmine_out folder contains the input MGF (MS/MS spectra) and CSV (metadata) files for positive, negative or ideally, both ionization modes.
The params folder contains all parameter files for annotation, dereplication and also the folder names to be used in the project:

- The fragnotator table.
- Primary negative mode adduct table.
- Secondary negative mode adduct table.
- Primary positive mode adduct table.
- Secondary positive mode adduct table.
- The main "params.yaml" file for global parameters.
- The secondary params files, one for each dereplication process.

The fragnotator file is a simple two-column table containing the annotation and the corresponding mass difference:

| loss  | mass      |
|-------|-----------|
| CH3   | 15.023475 |
| H2O   | 18.010565 |
| CO    | 27.994915 |

The adduct tables are all formed the same, the primary being used for triangulation and the secondary only being used to annotate the remaining ions once the neutral node is created.

| Adduct\_code     | Adduct        | Charge | Adduct\_mass| Mol\_multiplier| Complexity | Group |
|------------------|---------------|--------|-------------|----------------|------------|-------|
| M1\|m1H\|pC4H11N | [M-H+C4H11N]- | -1     | 72.081324   | 1              | 3          | H     |
| M1\|m1H\|pHCOOH  | [M-H+HCOOH]-  | -1     | 44.997655   | 1              | 3          | H     |
| M1\|m1H\|        | [M-H]-        | -1     | -1.007825   | 1              | 1          | H     |
| M1\|p1Cl\|       | [M+Cl]-       | -1     | 34.968853   | 1              | 2          | Cl    |

The user can add or delete rows in these tables to fit the needs of the experiment, or transfer adduct from the primary to the secondary tables if computing becomes too long. When transferring adduct from the primary to the secondary tables, the less abundant ion species should be prioritized as removing common species such as [M+H]+ would highly impact the triangulation in a negative way. Multiple charge adduct processing is not implemented as of yet, we would suggest only using single charge ions.  

The params yaml file contains all parameters to be used in the project. Each parameter has a short description as a comment and we would suggest using the default values to begin with. 
The other params files are dedicated to the dereplication with parameters specific to the dereplication to be carried out and the database file (in the databases folder) to be used.

Once all parameters are set, use the example MolNotator script provided to start the process. After most steps, CSV files are exported including a node table and an edge table. Networks can thus be visualized after each step using softwares like [Cytoscape](https://cytoscape.org/) by importing the two tables. The final network with molecules, adducts, in-source fragments, with dereplication and a degree of cosine clustering can be opened after the `Cosiner` function. Simplified versions of the network (only neutrals and adducts or neutrals only) can be produced after the `MolNet` function. 

Global networks containing all samples are produced at each step, but they can be divided to contain only the data for each specific sample. To do this, refer to the "export_samples" parameters in the params.yaml file.

## Installation
### Dependencies
Before installing MolNotator, make sure you have the following requirements installed:

- pandas == 1.3.5
- NumPy == 1.20.3
- matchms <= 0.6.2
- tqdm
- PyYaml

These dependencies can be installed using the following command :

```bash
 conda install pandas==1.3.5 numpy==1.20.3 tqdm pyyaml -c conda-forge
 conda install matchms==0.6.2 -c bioconda
```
### Via PyPI
We deploy the MolNotator package to [PyPi](https://pypi.org/project/MolNotator/). You can install MolNotator as a python module with:
```bash
 pip install MolNotator
```
**Note:** This is the recommended way ! 
### From source
If you cannot use the PyPi bundle or want to install MolNotator from source, we suggest these steps:\
Open a terminal and clone this repository using 
```bash
 git clone https://github.com/ZzakB/MolNotator.git
```
Move to the root directory of your MolNotator repository and run the following command in it 
```bash
 pip install .
```
**Note:** Be aware that you still have to install the above mentioned dependencies and link them correctly.

## Usage/Examples
MolNotator depends on a python [input file](https://github.com/ZzakB/MolNotator/blob/main/examples/input_file.py) to be runned. The example here under can be used as a template :
```python
import os 
import yaml
from MolNotator.duplicate_filter import duplicate_filter
from MolNotator.sample_slicer import sample_slicer
from MolNotator.fragnotator import fragnotator
from MolNotator.adnotator import adnotator
from MolNotator.mode_merger import mode_merger
from MolNotator.dereplicator import dereplicator
from MolNotator.cosiner import cosiner
from MolNotator.molnet import molnet
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

wd = 'set/path/'
os.chdir(wd)

for files in os.listdir(os.getcwd()):
    if files not in ['databases','input_files','params']:
        raise Exception('Potential output files already exist! They need to be removed or moved outside the working directory.')
        
with open("./params/params.yaml") as info:
    params = yaml.load(info, Loader=yaml.FullLoader)


# Duplicate filtering on MZmine's MGF and CSV files (NEG):
duplicate_filter(params = params,
                 ion_mode = "NEG")

# Duplicate filtering on MZmine's MGF and CSV files (POS):
duplicate_filter(params = params,
                 ion_mode = "POS")

# Slicing the negative mode MGF file
sample_slicer(params = params,
           ion_mode = "NEG")

# Slicing the positive mode MGF file
sample_slicer(params = params,
           ion_mode = "POS")

# Use fragnotator on the negative mode sliced MGF files
fragnotator(params = params,
            ion_mode = "NEG")

# Use fragnotator on the positive mode sliced MGF files
fragnotator(params = params,
            ion_mode = "POS")

# Use adnotator on the negative mode data
adnotator(params = params,
          ion_mode = "NEG")

# Use adnotator on the positive mode data
adnotator(params = params,
          ion_mode = "POS")

# Use Moder Merger to merge negative and positive mode data :
mode_merger(params = params)


# Dereplicate the data using the database specified in the YAML file
for db_params in params['db_params']:
    print("Dereplicating using the " + db_params + " file...")
    with open("./params/" + db_params) as info:
        db_params = yaml.load(info, Loader=yaml.FullLoader)    
    dereplicator(params = params,
                 db_params = db_params)

# Compute cosine similarity between some nodes.
cosiner(params = params)

# Produce molecular networks, neutral nodes only
molnet(params = params)
```
Then MolNotator can be runned using the input file above mentioned with the following command : 
```python
python input_file.py
```
**Note:** The output is written to the working directory.

## License
MolNotator is published under the [MIT](https://choosealicense.com/licenses/mit/) licence. 
For more information, please read the [LICENSE](https://github.com/ZzakB/MolNotator/blob/main/LICENCE) file.
Using MolNotator in your commercial or non-commercial project is generally possible when giving a proper reference to this project and the related paper.
 

## Citation Information
If you are using MolNotator in your work, please cite :
> [(1) Olivier-Jimenez, D.; Bouchouireb, Z.; Ollivier, S.; Mocquard, J.; Allard, P.-M.; Bernadat, G.; Chollet-Krugler, M.; Rondeau, D.; Boustie, J.; van der Hooft, J. J. J.; Wolfender, J.-L. From Mass Spectral Features to Molecules in Molecular Networks: A Novel Workflow for Untargeted Metabolomics, 2021. https://doi.org/10.1101/2021.12.21.473622.](https://doi.org/10.1101/2021.12.21.473622)
## Contact / Maintainer
Do you have feature requests, found a bug or want to use MolNotator in your project ?  
Please get in touch : `damien.olivier.jimenez@gmail.com`
