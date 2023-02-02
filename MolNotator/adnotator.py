import os
import pandas as pd
import numpy as np
import sys
from tqdm import tqdm
from matchms.importing import load_from_mgf
from matchms.similarity import ModifiedCosine
from MolNotator.others.global_functions import *
from MolNotator.others.rt_slicer import rt_slicer

def adnotator(params : dict, ion_mode : str):
    """
    Finds neutral-adduct pairs.

    Parameters
    ----------
    params : dict
        Dictionary containing the global parameters for the process..
    ion_mode : str
        Either "POS" or "NEG", ion mode for the data.

    Returns
    -------
    None.

    """
    
    # Load parameters:
    rt_error= params['an_rt_error']
    idx_column = params['index_col']
    
    # If the retention time is in minutes
    if params['rt_unit'] == "m":
        rt_error = rt_error/60
    
    if ion_mode == "NEG":
        in_path_full= params['neg_out_0']
        csv_file= params['neg_csv']
        spectrum_file= params['neg_mgf']
        in_path_spec= params['neg_out_1']
        in_path_csv= params['neg_out_2']
        out_path_full= params['neg_out_3_1']
        out_path_samples= params['neg_out_3_2']
        adduct_table_primary= params['an_addtable_primary_neg']
        adduct_table_secondary= params['an_addtable_secondary_neg']
    elif ion_mode == "POS":
        in_path_full= params['pos_out_0']
        csv_file= params['pos_csv']
        spectrum_file= params['pos_mgf']
        in_path_spec= params['pos_out_1']
        in_path_csv= params['pos_out_2']
        out_path_full= params['pos_out_3_1']
        out_path_samples= params['pos_out_3_2']
        adduct_table_primary= params['an_addtable_primary_pos']
        adduct_table_secondary= params['an_addtable_secondary_pos']
    else:
        print('Ion mode must be either "NEG" or "POS"')
        return
    
    # Load files
    csv_file = in_path_full + csv_file
    spectrum_file = in_path_full + spectrum_file
    
    input_files = pd.DataFrame()
    input_files['spectrum_file'] = os.listdir(in_path_spec)
    input_files['base_name'] = input_files['spectrum_file'].copy().str.replace(f'{ion_mode}_', "")
    input_files['base_name'] = input_files['base_name'].str.replace('.mgf', '', regex = False)
    
    adduct_table_primary = pd.read_csv("./params/" + adduct_table_primary, sep = "\t")
    adduct_table_secondary = pd.read_csv("./params/" + adduct_table_secondary, sep = "\t")
    adduct_table_merged= adduct_table_primary.append(adduct_table_secondary, ignore_index = True)
    
    # Create output folder
    if not os.path.isdir(out_path_full) :
        os.mkdir(out_path_full)
    
    # If adnotator was already used, load previous files, otherwise create new
    if os.path.isfile(out_path_full + "cross_sample_annotations.csv"):
        cross_annotations = pd.read_csv(out_path_full + "cross_sample_annotations.csv", index_col = idx_column, dtype = str)
        cross_annotations = cross_annotations.replace({np.nan: None})
        cross_courts = pd.read_csv(out_path_full + "cross_sample_courts.csv", index_col = idx_column)
        cross_courts = cross_courts.replace({np.nan: None})
        cross_houses = pd.read_csv(out_path_full + "cross_sample_houses.csv", index_col = idx_column, dtype = str)
        cross_houses = cross_houses.replace({np.nan: None})
        cross_points = pd.read_csv(out_path_full + "cross_sample_points.csv", index_col = idx_column, dtype = float)
        cross_points = cross_points.replace({np.nan: None})
        cross_rules = pd.read_csv(out_path_full + "cross_sample_rules.csv", index_col = idx_column, dtype = str)
        cross_rules = cross_rules.replace({np.nan: None})
        cross_neutrals = pd.read_csv(out_path_full + "cross_sample_neutrals.csv", index_col = idx_column, dtype = str)
        cross_neutrals = cross_neutrals.replace({np.nan: None})
        ion_ids = get_ion_ids(spectrum_file, params)
    else:
        cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals, ion_ids = cross_sample_tables(spectrum_file, params)
    
    # Start processing sample by sample
    for x in input_files.index :
        sample_base_name = input_files.loc[x, "base_name"]
        if sample_base_name in cross_annotations.columns : continue
        
        file_basename = input_files.at[x, 'base_name']
        print("Processing " + file_basename)
        
        # Load files for the processed sample
        edge_table = pd.read_csv(f"{in_path_csv}{ion_mode}_{file_basename}_edges.csv", 
                                 index_col = 'Index')
        node_table = pd.read_csv(f"{in_path_csv}{ion_mode}_{file_basename}_nodes.csv", 
                                index_col = params['index_col'])
        subspectrum_file = in_path_spec + input_files.loc[x, "spectrum_file"]
        spectrum_list = list(load_from_mgf(subspectrum_file))
        spectrum_list = [Spectrum_processing(s) for s in spectrum_list]
        
        # Create dataframes to store results
        duplicate_table = pd.DataFrame() # Stores duplicate spectra
        merged_table = pd.DataFrame() # Store final results
        
        # Process each spectrum for the sample
        for i in tqdm(node_table.index) :
            
            # Load data for first ion to be analysed
            if abs(int(node_table.loc[i, params['charge_field']])) > 1 : continue
            ion1_rt = node_table.loc[i, params['rt_field']]
            ion1_mz = node_table.loc[i, params['mz_field']]
            ion1_spec_id = node_table.loc[i, 'spec_id']
    
            # Slice the DF to keep only ions coeluting with current ion
            coelution_table = rt_slicer(ion1_rt, rt_error, i, node_table,
                                        params['rt_field'])
            
            # Filter the table to remove multicharged ions
            coelution_table = coelution_table[abs(coelution_table[params['charge_field']]) <= 1]
        
            # Produce neutral table (all neutrals for the considered ion):
            neutral_table = neutral_tabler(ion1_mz, adduct_table_primary)
        
            # Produce ion hypotheses table (all adducts for all neutrals)
            ion_hypotheses_table = ion_hypotheses_tabler(neutral_table)
        
            # Award points for each combination of hypotheses based on the
            # existance of the Ion2 in the coelution table.
            ion_hypotheses_table = point_counter(ion_hypotheses_table,
                                                 coelution_table, params)
            
            # Cosine similarity are integrated to the scores of the hypotheses
            ion_hypotheses_table, tmp_duplicates = cosine_validation(i, node_table, 
                                                                     spectrum_list,
                                                                     ion_hypotheses_table,
                                                                     adduct_table_primary,
                                                                     params)
            
            # Add duplicates found previously to the duplicate table
            duplicate_table = duplicate_table.append(tmp_duplicates, ignore_index = True)
            
            # Species rules points : award points for annotations that can be confirmed
            ion_hypotheses_table = species_rules(ion1_spec_id, ion_hypotheses_table,
                                                 adduct_table_primary, node_table,
                                                 spectrum_list, params, ion_mode)
    
            # Award points for solvent complex confirmation 
            ion_hypotheses_table = complex_points(neutral_table, ion_hypotheses_table)
            
            # Add points if precursor-fragment connexion exists
            ion_hypotheses_table = fragnotator_points(i, ion_hypotheses_table,
                                        edge_table, node_table.loc[i, 'status'])
    
            # Get the adduct code for each ionisation hypothesis
            ion_hypotheses_table = get_adduct_code(ion_hypotheses_table, neutral_table)
            
            # Add final score (weighted points)
            ion_hypotheses_table = weighted_points(ion_hypotheses_table, adduct_table_primary)
        
            # Report all hypotheses for in current ion in a merged adducts table
            merged_table_local = merged_adducts_table(ion_hypotheses_table, i)

            # Add the local results to the global merged table
            if len(merged_table_local.index) == 0 : continue
            merged_table = merged_table.append(merged_table_local)
            merged_table.reset_index(inplace = True, drop = True)

        # Produce the cohort table containing all ion hypotheses in the sample x
        cohort_table = cohort_tabler(neutral_table, merged_table, node_table)

        # Produce the supercohorts table to merge redundant cohorts in the sample x
        supercohorts_table, transition_table = supercohorts_tabler(cohort_table)

        # Produce the court table
        court_table = get_court_table(supercohorts_table)

        # Produce house table
        court_table = house_selection(court_table, supercohorts_table, node_table, transition_table,
                            merged_table, spectrum_list, params)

        # Update cross sample tables with the results
        cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals = cross_sample_report(court_table, cross_annotations, cross_points, cross_courts,
                            cross_houses, cross_rules, cross_neutrals, 
                            sample_base_name, supercohorts_table, 
                            adduct_table_primary, adduct_table_merged, node_table,
                            ion_mode, duplicate_table, spectrum_list, params)
    
        # Add secondary adducts
        cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals = get_secondary_adducts(cross_annotations, cross_points, cross_courts,
                                  cross_houses, cross_rules, cross_neutrals, sample_base_name,
                                  node_table, spectrum_list, adduct_table_primary, 
                                  adduct_table_secondary, ion_mode, params)
    
        # Export data
        cross_annotations.to_csv(out_path_full+"cross_sample_annotations.csv", index_label = idx_column)
        cross_courts.to_csv(out_path_full+"cross_sample_courts.csv", index_label = idx_column)
        cross_houses.to_csv(out_path_full+"cross_sample_houses.csv", index_label = idx_column)
        cross_points.to_csv(out_path_full+"cross_sample_points.csv", index_label = idx_column)
        cross_rules.to_csv(out_path_full+"cross_sample_rules.csv", index_label = idx_column)
        cross_neutrals.to_csv(out_path_full+"cross_sample_neutrals.csv", index_label = idx_column)
    
    # Load and merge fragnoted 
    merged_node_table, merged_edge_table = get_merged_tables(input_files, ion_mode, params)
    merged_edge_table['Adnotation'] = [None]*len(merged_edge_table)
    merged_node_table['Adnotation'] = [None]*len(merged_node_table)
    
    # Resolve the sample-specific unresolved ions by cross examining all samples:
    cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals = annotation_resolver(cross_annotations, cross_points, cross_courts, cross_houses, cross_rules, cross_neutrals)
    
    # Rebuild houses at a cross-sample level to check compatibilities between annotations
    cross_court_table, singletons = get_cross_courts(cross_annotations, cross_courts, cross_houses)    
    
    # Select the most likeley neutrals from the cross_court_table and report
    # the data on the merged node and edge tables.
    spectrum_list = list(load_from_mgf(spectrum_file))
    spectrum_list = [Spectrum_processing(s) for s in spectrum_list]
    merged_node_table, merged_edge_table = cross_neutral_selection(spectrum_list, cross_court_table, cross_annotations,
                                cross_neutrals, merged_node_table, merged_edge_table,
                                adduct_table_merged, ion_mode, params)
    
    
    # Update node table status:
    merged_node_table = update_node_table(merged_node_table, merged_edge_table,
                        cross_annotations, cross_rules, adduct_table_merged, params)
    
    # Update edge table with singletons:
    merged_edge_table = update_edge_table(merged_node_table, merged_edge_table,
                                          adduct_table_merged)
    
    # Export merged tables
    merged_node_table.to_csv(f"{out_path_full}node_table.csv", index_label = idx_column)
    merged_edge_table.to_csv(f"{out_path_full}edge_table.csv", index_label = "Index")
    
    
    # Also export sample wise slices of the tables for viewing at this stage
    if params['an_export_samples'] : 
        samplewise_export(merged_edge_table, merged_node_table, "adnotator", ion_mode, params)
    return
