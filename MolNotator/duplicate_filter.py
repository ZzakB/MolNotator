"""duplicate_filter.py - duplicate_filter function for MolNotator"""
import os
import pandas as pd
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters
from matchms.exporting import save_as_mgf
from pandas.core.common import flatten
from MolNotator.others.spectrum_extractor import spectrum_extractor
from MolNotator.others.duplicate_finder import duplicate_finder


def Spectrum_processing(s):
    s = default_filters(s)
    return s

def duplicate_filter(params : dict, ion_mode : str):
    """
    Finds duplicate ions in a metabolomics experiment by loading the CSV and 
    spectrum files and by comparing the RT, m/z and cosine similarity values
    between ions. Ions with close values are deemed duplicates. Along with the 
    other duplicates, they will be deleted, leaving only a representative ion 
    (often the most intense). This deletion is operated on the CSV and the 
    spectrum file, which will both be exported in the duplicate filter folder.
    Parameters
    ----------
    params : dict
        Dictionary containing the global parameters for the process.
    ion_mode : str
        Either "POS" or "NEG", ion mode for the data.
    Returns
    -------
    CSV and MGF files, filtered, in the duplicate filter folder.
    """    
    
    # Load parameters
    index_col = params['index_col']
    rt_field = params['rt_field']
    mz_field = params['mz_field']
    
    if ion_mode == "NEG":
        csv_file = params['neg_csv']
        spectrum_list= params['neg_mgf']
        out_path= params['neg_out_0']
    elif ion_mode == "POS" :
        csv_file = params['pos_csv']
        spectrum_list= params['pos_mgf']
        out_path= params['pos_out_0']
    
    # Save mgf and csv file names for input
    spectrum_file = spectrum_list
    csv_name = csv_file
    
    # create output dir:
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    
    # Load MZmine mgf and csv files
    print("Loading MGF and CSV files...")
    spectrum_list = list(load_from_mgf(f'{params["input_dir"]}{spectrum_list}'))
    csv_file = pd.read_csv(f'{params["input_dir"]}{csv_file}', index_col = index_col)
    
    # Format columns with ion modes:
    new_cols = csv_file.columns.tolist()
    new_cols = [ion_mode + "_" + col if (params['col_suffix'] in col and col[:4] != f"{ion_mode}_") else col for col in new_cols]
    csv_file.columns = new_cols
    
    # Extract data from the MGF file
    print('Extracting MGF metadata...')
    node_table = spectrum_extractor(spectrum_list)
    
    # Merge node_table from the spectrum file, and the CSV table
    node_table['TIC'] = [s.peaks.intensities.sum() for s in spectrum_list] # Add TIC
    node_table[index_col.lower()] = node_table[index_col.lower()].astype(int)
    node_table.set_index(index_col.lower(), inplace = True, drop = True)
    drop_cols = csv_file.columns.intersection(set(node_table.columns))
    node_table.drop(drop_cols, axis = 1, inplace = True)
    node_table = csv_file.merge(node_table, left_index = True, right_index = True)
    node_table[f'{rt_field}'] = node_table[f'{rt_field}'].astype(float)
    node_table[f'{mz_field}'] = node_table[f'{mz_field}'].astype(float)
    
    # Correct rt unit from m to s if relevant
    if params['rt_unit'] == 'm':
        node_table[f'{rt_field}'] = node_table[f'{rt_field}']*60

    # Add spec_id, for the relative position of ions in the spectrum file
    node_table.insert(0, 'spec_id', range(len(node_table)))

    # If this step must be skipped :
    if params['df_skip'] :
        node_table.to_csv(f'{out_path}{csv_name}')
        save_as_mgf(spectrum_list, f'{out_path}{spectrum_file}')
        return
    
    # Filtering the MGF file
    print('Filtering MGF file...')
    spectrum_list = [Spectrum_processing(s) for s in spectrum_list]
    
    # Get duplicates & delete them
    duplicate_table = duplicate_finder(node_table, spectrum_list, params, ion_mode)
    dropped_ions = list(flatten(duplicate_table['dropped']))
    kept_ions = list(set(node_table.index) - set(dropped_ions))
    kept_ions.sort()
    kept_ions = [node_table.loc[i, "spec_id"] for i in kept_ions]
    mgf_file_new = [spectrum_list[i] for i in kept_ions]
    node_table_new = node_table.drop(dropped_ions)
    
    # Reset the spec_id to account for dropped ions
    node_table_new['spec_id'] = range(len(node_table_new))
    
    # Export the data
    print('Exporting MGF and CSV files...')
    save_as_mgf(mgf_file_new, f'{out_path}{spectrum_file}')
    node_table_new.to_csv(f'{out_path}{csv_name}', index_label = index_col)
    perc = round(100*(len(dropped_ions)/len(spectrum_list)),1)
    print('Export finished.')
    print(f'{len(dropped_ions)} ions removed out of {len(spectrum_list)} ({perc}%)')
    return
