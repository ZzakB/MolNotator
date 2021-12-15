import os 
import yaml
from MolNotator.Duplicate_filter import Duplicate_filter
from MolNotator.MGF_sample_slicer import mgf_slicer
from MolNotator.Fragnotator import Fragnotator
from MolNotator.Adnotator import Adnotator
from MolNotator.MGF_updater import MGF_updater
from MolNotator.Mode_merger import Moder_merger
from MolNotator.Dereplicator import Dereplicator
from MolNotator.Cosiner import Cosiner
from MolNotator.MolNet import MolNet

wd = './working_directory'
os.chdir(wd)

for files in os.listdir(os.getcwd()):
    if files not in ['databases','mzmine_out','params']:
        raise Exception('Potential output files already exist! They need to be removed or moved outside the working directory.')
        
with open("./params/params.yaml") as info:
    params = yaml.load(info, Loader=yaml.FullLoader)


# Duplicate filtering on MZmine's MGF and CSV files (NEG):
Duplicate_filter(params = params,
                 ion_mode = "NEG")

# Duplicate filtering on MZmine's MGF and CSV files (POS):
Duplicate_filter(params = params,
                 ion_mode = "POS")

# Slicing the negative mode MGF file
mgf_slicer(params = params,
           ion_mode = "NEG")

# Slicing the positive mode MGF file
mgf_slicer(params = params,
           ion_mode = "POS")

# Use fragnotator on the negative mode sliced MGF files
Fragnotator(params = params,
            ion_mode = "NEG")

# Use fragnotator on the positive mode sliced MGF files
Fragnotator(params = params,
            ion_mode = "POS")

# Use adnotator on the negative mode data
Adnotator(params = params,
          ion_mode = "NEG")

# Use adnotator on the positive mode data
Adnotator(params = params,
          ion_mode = "POS")

# Use Moder Merger to merge negative and positive mode data :
Moder_merger(params = params)

# Update the MGF files and the node tables with SIRIUS formulas and other annotations
MGF_updater(params = params)

# Dereplicate the data using the database specified in the YAML file
for db_params in params['db_params']:
    print("Dereplicating using the " + db_params + " file...")
    with open("./params/" + db_params) as info:
        db_params = yaml.load(info, Loader=yaml.FullLoader)    
    Dereplicator(params = params,
                 db_params = db_params)

# Compute cosine similarity between some nodes.
Cosiner(params = params)

# Produce molecular networks, neutral nodes only
MolNet(params = params)
