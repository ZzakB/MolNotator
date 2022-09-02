import os 
import yaml
from MolNotator import duplicate_filter
from MolNotator.sample_slicer import sample_slicer
from MolNotator.fragnotator import fragnotator
from MolNotator.adnotator import adnotator
from MolNotator.mode_merger import mode_merger
from MolNotator.dereplicator import dereplicator
from MolNotator.cosiner import cosiner
from MolNotator.molnet import molnet

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
