import os
import pandas as pd
from pandas.core.common import flatten
import sys
from tqdm import tqdm
from matchms.importing import load_from_mgf
from MolNotator.others.species_validators import *
from MolNotator.others.global_functions import *


def mode_merger(params : dict):

    # Load parametes
    skip = params['mm_skip']
    charge_field = params['charge_field']
    mz_field = params['mz_field']
    rt_field = params['rt_field']
    col_suffix = params['col_suffix']
    idx_column = params['index_col']

    out_full= params['mix_out_4_1']
    out_samples= params['mix_out_4_2']
    
    adduct_table_primary_neg= params['mm_addtable_primary_neg']
    adduct_table_secondary_neg= params['mm_addtable_secondary_neg']
    adduct_table_primary_pos= params['mm_addtable_primary_pos']
    adduct_table_secondary_pos= params['mm_addtable_secondary_pos']
    mass_error= params['mm_mass_error']
    prec_mass_error = params['mm_prec_mass_error']
    rt_error= params['mm_rt_error']  
    bnr_list_neg = params['mm_bnr_neg']
    bnr_list_pos = params['mm_bnr_pos']

    # If the retention time is in minutes
    if params['rt_unit'] == "m":
        rt_error = rt_error/60
        
    # Make out dir
    if not os.path.isdir(out_full) :
        os.mkdir(out_full)
    if not os.path.isdir(out_samples) :
        os.mkdir(out_samples)

    if skip == True:
        ion_mode = params['process_mode']
        # Loading files: POS and NEG node tables, edge tables and the global CSV
        if ion_mode == "POS":
            in_path_full_pos= params['pos_out_0']
            pos_spectrum_file= params['pos_mgf']
            pos_csv= params['pos_csv']
            adnotation_pos= params['pos_out_3_1']
            csv_table = pd.read_csv(in_path_full_pos + pos_csv, index_col = idx_column)
            edge_table = pd.read_csv(adnotation_pos + "edge_table.csv", index_col = "Index")
            node_table = pd.read_csv(adnotation_pos + "node_table.csv", index_col = idx_column)
        elif ion_mode == "NEG":
            in_path_full_neg= params['neg_out_0']
            neg_spectrum_file= params['neg_mgf']
            neg_csv= params['neg_csv']
            adnotation_neg= params['neg_out_3_1']
            csv_table = pd.read_csv(in_path_full_neg + neg_csv, index_col = idx_column)
            edge_table = pd.read_csv(adnotation_neg + "edge_table.csv", index_col = "Index")
            node_table = pd.read_csv(adnotation_neg + "node_table.csv", index_col = idx_column)
        else:
            print("Specify ion mode for single ion mode processing in the parameters")
            return

        # Filter csv table
        csv_table = csv_table.loc[:,csv_table.columns.str.contains(col_suffix)]
        csv_table.columns = csv_table.columns.str.replace(col_suffix, "", regex = False)
        csv_table.columns = csv_table.columns.str.replace(f"{ion_mode}_", "", regex = False)
        
        # Add properties 
        node_table.insert(node_table.columns.get_loc('status') + 1 , "status_universal", node_table['status'].copy())
        node_table['status'] = f"{ion_mode.lower()}_"+node_table['status']
        edge_table.insert(edge_table.columns.get_loc('status') + 1 , "status_universal", edge_table['status'].copy())
        edge_table['status'] = f"{ion_mode.lower()}_"+edge_table['status']
        node_table['ion_mode'] = ion_mode
        edge_table['ion_mode'] = ion_mode

        # Add cluster IDs
        node_table = get_cluster_ids(node_table, edge_table)

        # Add samples
        node_table = get_samples(node_table, edge_table, csv_table, ion_mode)

        # Reset IDs in edge table
        node_table.reset_index(inplace = True, drop = False)
        edge_table = refresh_edge_table(node_table, edge_table, idx_column)       
        
        # Export
        node_table.to_csv(out_full + "node_table.csv", index_label = "Index")
        edge_table.to_csv(out_full + "edge_table.csv", index_label = "Index")
        return
    
    # If double mode (and do actual mode merging)
    in_path_full_neg= params['neg_out_0']
    in_path_full_pos= params['pos_out_0']
    neg_spectrum_file= params['neg_mgf']
    pos_spectrum_file= params['pos_mgf']
    pos_csv= params['pos_csv']
    neg_csv= params['neg_csv']
    adnotation_pos= params['pos_out_3_1']
    adnotation_neg= params['neg_out_3_1']
    
    
    # Load and filter the global spectrum files in pos and neg
    print("Loading and filtering NEG spectrum file...")
    neg_spectrum_list = list(load_from_mgf(in_path_full_neg + neg_spectrum_file))
    neg_spectrum_list = [Spectrum_processing(s) for s in neg_spectrum_list]
    print("Loading and filtering POS spectrum file...")
    pos_specturm_list = list(load_from_mgf(in_path_full_pos + pos_spectrum_file))
    pos_specturm_list = [Spectrum_processing(s) for s in pos_specturm_list]
    
    
    #Load adduct tables
    adduct_table_primary_neg = pd.read_csv("./params/" + adduct_table_primary_neg, sep = "\t")
    adduct_table_secondary_neg = pd.read_csv("./params/" + adduct_table_secondary_neg, sep = "\t")
    adduct_table_merged_neg = adduct_table_primary_neg.append(adduct_table_secondary_neg, ignore_index = True)
    adduct_table_primary_pos = pd.read_csv("./params/" + adduct_table_primary_pos, sep = "\t")
    adduct_table_secondary_pos = pd.read_csv("./params/" + adduct_table_secondary_pos, sep = "\t")
    adduct_table_merged_pos = adduct_table_primary_pos.append(adduct_table_secondary_pos, ignore_index = True)
    
    # Produce base neutral adduct tables
    adduct_table_base_neg = [adduct_table_merged_neg.index[adduct_table_merged_neg['Adduct'] == adduct][0] for adduct in bnr_list_neg]
    adduct_table_base_neg = adduct_table_merged_neg.loc[adduct_table_base_neg]
    adduct_table_base_pos = [adduct_table_merged_pos.index[adduct_table_merged_pos['Adduct'] == adduct][0] for adduct in bnr_list_pos]
    adduct_table_base_pos = adduct_table_merged_pos.loc[adduct_table_base_pos]
    
    # Loading files: POS and NEG node tables, edge tables and the global CSV
    neg_table_global = pd.read_csv(in_path_full_neg + neg_csv, index_col = idx_column)
    pos_table_global = pd.read_csv(in_path_full_pos + pos_csv, index_col = idx_column)

    edge_table_all_neg = pd.read_csv(adnotation_neg + "edge_table.csv", index_col = "Index")
    node_table_all_neg = pd.read_csv(adnotation_neg + "node_table.csv", index_col = idx_column)
    edge_table_all_pos = pd.read_csv(adnotation_pos + "edge_table.csv", index_col = "Index")
    node_table_all_pos = pd.read_csv(adnotation_pos + "node_table.csv", index_col = idx_column)
    
    # Setting ion modes in tables
    node_table_all_neg['ion_mode'] = ['NEG']*len(node_table_all_neg)
    edge_table_all_neg['ion_mode'] = ['NEG']*len(edge_table_all_neg)
    node_table_all_pos['ion_mode'] = ['POS']*len(node_table_all_pos)
    edge_table_all_pos['ion_mode'] = ['POS']*len(edge_table_all_pos)
    
    # New indexes will be used when merging POS and NEG tables, old IDs are saved in a dedicated column
    node_table_all_neg.insert(0, idx_column, node_table_all_neg.index)
    node_table_all_pos.insert(0, idx_column, node_table_all_pos.index)
    

    # Import sample to node tables
    
    # filter global tables to keep only samples
    neg_table_global = neg_table_global.loc[:,neg_table_global.columns.str.contains(col_suffix)]
    pos_table_global = pos_table_global.loc[:,pos_table_global.columns.str.contains(col_suffix)]

    # Remove prefixes and suffixes from the columns 
    neg_table_global.columns = neg_table_global.columns.str.replace('NEG_', '', regex = False)
    neg_table_global.columns = neg_table_global.columns.str.replace(col_suffix, '', regex = False)

    pos_table_global.columns = pos_table_global.columns.str.replace('POS_', '', regex = False)
    pos_table_global.columns = pos_table_global.columns.str.replace(col_suffix, '', regex = False)

    # Add a sample column to the node tables
    node_table_all_neg["samples"] = ['']*len(node_table_all_neg)
    node_table_all_pos["samples"] = ['']*len(node_table_all_pos)

    # Get ion IDs (no neutrals)
    neg_ions_idx = node_table_all_neg.index[node_table_all_neg['status'] != "neutral"]
    pos_ions_idx = node_table_all_pos.index[node_table_all_pos['status'] != "neutral"]

    # Add samples to the samples column (for ions)
    print("Adding samples to neg ions...")
    for ion in tqdm(neg_ions_idx):
        samples = list(neg_table_global.columns[neg_table_global.loc[ion] > 0])
        node_table_all_neg.loc[ion, "samples"] = '|'.join(samples)
    print("Adding samples to pos ions...")
    for ion in tqdm(pos_ions_idx):
        samples = list(pos_table_global.columns[pos_table_global.loc[ion] > 0])
        node_table_all_pos.loc[ion, "samples"] = '|'.join(samples)
    
    # Add samples to neutrals
    
    # Get neutral indexes
    neg_neutrals_idx = node_table_all_neg.index[node_table_all_neg['status'] == "neutral"]
    pos_neutrals_idx = node_table_all_pos.index[node_table_all_pos['status'] == "neutral"]

    # Add samples to the samples column (neutrals)
    print("Adding samples to neg neutrals...")
    for neutral in tqdm(neg_neutrals_idx):
        tmp_ions = list(edge_table_all_neg['node_2'][edge_table_all_neg['node_1'] == neutral])
        tmp_samples = []
        for ion in tmp_ions:
            tmp_samples += node_table_all_neg.loc[ion, "samples"].split('|')
        tmp_samples = list(set(tmp_samples))
        tmp_samples.sort()
        tmp_samples = '|'.join(tmp_samples)
        node_table_all_neg.loc[neutral, "samples"] = tmp_samples
    print("Adding samples to pos neutrals...")
    for neutral in tqdm(pos_neutrals_idx):
        tmp_ions = list(edge_table_all_pos['node_2'][edge_table_all_pos['node_1'] == neutral])
        tmp_samples = []
        for ion in tmp_ions:
            tmp_samples += node_table_all_pos.loc[ion, "samples"].split('|')
        tmp_samples = list(set(tmp_samples))
        tmp_samples.sort()
        tmp_samples = '|'.join(tmp_samples)
        node_table_all_pos.loc[neutral, "samples"] = tmp_samples
    
    # Produce the merged mode edge and node tables
    node_table_ions_neg= node_table_all_neg[node_table_all_neg['status'] != "neutral"]
    node_table_ions_pos= node_table_all_pos[node_table_all_pos['status'] != "neutral"]
    merged_node_table = pd.DataFrame()
    merged_node_table = merged_node_table.append(node_table_ions_neg, ignore_index = True)
    merged_node_table = merged_node_table.append(node_table_ions_pos, ignore_index = True)
    merged_edge_table = pd.DataFrame()
    merged_edge_table = merged_edge_table.append(edge_table_all_neg, ignore_index = True)
    merged_edge_table = merged_edge_table.append(edge_table_all_pos, ignore_index = True)
    
    #Check from NEG to POS if there are neutral nodes that might match
    
    # Get only neutrals
    neg_neutrals = node_table_all_neg[node_table_all_neg['status'] == "neutral"].copy()
    pos_neutrals = node_table_all_pos[node_table_all_pos['status'] == "neutral"].copy()

    # Reset the original index column
    neg_neutrals[idx_column] = [None]*len(neg_neutrals)
    pos_neutrals[idx_column] = [None]*len(pos_neutrals)

    # Start merging neutrals, NEG -> POS    
    merged_neutrals = []
    for neutral in tqdm(neg_neutrals_idx):
        
        # Get neutral attributes, and samples
        mz_neg = neg_neutrals.loc[neutral, mz_field]
        rt_neg = neg_neutrals.loc[neutral, rt_field]
        samples_neg = set(neg_neutrals.loc[neutral, "samples"].split('|'))
        
        # Find a neutral in the opposite mode with similar mass and RT
        counter_pos = pos_neutrals[pos_neutrals[mz_field].between(mz_neg - mass_error, mz_neg + mass_error, inclusive = "both")]
        counter_pos = counter_pos[counter_pos[rt_field].between(rt_neg - rt_error, rt_neg + rt_error, inclusive = "both")] 
        if len(counter_pos) == 0 : continue
    
        # Add an RT gap field
        counter_pos['d_rt'] = abs(counter_pos[rt_field]-rt_neg)

        # Find shared samples betwen neutrals (must share as many samples as possible)
        shared_samples = []
        for i in counter_pos.index:
            samples_pos = counter_pos.loc[i, "samples"].split('|')
            shared_samples.append(len(samples_neg.intersection(samples_pos)))

        # add a shared samples column and filter 
        counter_pos['shared_samples'] = shared_samples
        counter_pos = counter_pos[counter_pos['shared_samples'] > 0]
        
        # Select first by most shared samples, then by smallest drt
        counter_pos = counter_pos[counter_pos['shared_samples'] == counter_pos['shared_samples'].max()]
        counter_pos = counter_pos[counter_pos['d_rt'] == counter_pos['d_rt'].min()]
        
        # Stuff that shouldn't happen, but who knows
        if len(counter_pos) > 1 : sys.exit("MORE THAN 1 COUTNER NEUTRAL AT SAME RT")

        # If no shared samples, skip
        elif len(counter_pos) == 0 : continue

        # If match found
        else:
            mz_pos = counter_pos[mz_field].iloc[0]
            rt_pos = counter_pos[rt_field].iloc[0]
            
            # Add to the merged neutrals list
            merged_neutrals.append((neutral, mz_neg, rt_neg, counter_pos.index[0], mz_pos, rt_pos))
    
    # Create the merged neutrals dataframe
    merged_neutrals = pd.DataFrame(merged_neutrals, columns = ['neg_neutral', 'mass_neg',
                                                               'rt_neg', 'pos_neutral',
                                                               'mass_pos', 'rt_pos'])

    # Add a d_rt column
    merged_neutrals['d_rt'] = abs(merged_neutrals['rt_pos'] - merged_neutrals['rt_neg'])

    # Sometimes, several NEG neutrals might match the same POS neutral. 
    # The NEG neutral with the closest RT to the POS neutral is chosen
    for neutral in merged_neutrals['pos_neutral'].unique():
        tmp_table = merged_neutrals[merged_neutrals['pos_neutral'] == neutral].copy()
        if len(tmp_table) > 1 : 
            tmp_table.sort_values('d_rt', inplace = True)
            for i in tmp_table.index[1:]:
                merged_neutrals.drop(i, inplace = True)
    
    # A transitions table is produced containing the neutrals to be merged and their 
    # new ID. Later will also contain non-merged neutrals and their new IDs
    transitions_table = list()
    print("Adding merged neutrals to the merged node table...")
    for i in tqdm(merged_neutrals.index):
        
        # Set the future ID for the neutral node
        new_idx = merged_node_table.index.max() + 1 

        # Get the old neutral IDs
        neutral_neg = merged_neutrals.loc[i, "neg_neutral"]
        neutral_pos = merged_neutrals.loc[i, "pos_neutral"]
        
        # Get neg ions connected to the neg neutral
        neg_ions = list(edge_table_all_neg['node_2'][edge_table_all_neg['node_1'] == neutral_neg])

        # Get neg neutral attributes
        neg_count = len(neg_ions)
        neg_mz = node_table_all_neg.loc[neutral_neg, mz_field]
        neg_rt = node_table_all_neg.loc[neutral_neg, rt_field]
        neg_tic = node_table_all_neg.loc[neutral_neg, "TIC"]
        neg_samples = node_table_all_neg.loc[neutral_neg, "samples"].split('|')

        # Get pos ions connected to the pos neutrals
        pos_ions = list(edge_table_all_pos['node_2'][edge_table_all_pos['node_1'] == neutral_pos])

        # Get pos neutral attributes
        pos_count = len(pos_ions)
        pos_mz = node_table_all_pos.loc[neutral_pos, mz_field]
        pos_rt = node_table_all_pos.loc[neutral_pos, rt_field]
        pos_tic = node_table_all_pos.loc[neutral_pos, "TIC"]
        pos_samples = node_table_all_pos.loc[neutral_pos, "samples"].split('|')
        
        # Produce combined attributes
        mix_mz = round(((neg_mz*neg_count) + (pos_mz*pos_count))/(neg_count + pos_count), 4)
        mix_rt = round(((neg_rt*neg_count) + (pos_rt*pos_count))/(neg_count + pos_count), 3)
        mix_tic = neg_tic + pos_tic
        mix_samples = list(set(neg_samples + pos_samples))
        mix_samples.sort()
        mix_samples = '|'.join(mix_samples)

        # Report data to the merged node table
        merged_node_table.loc[new_idx] = [None]*len(merged_node_table.columns)
        
        merged_node_table.loc[new_idx, [mz_field, rt_field, "TIC", "rule_points",
                                        charge_field, "status", "adduct_count",
                                        'ion_mode', 'samples']] = [mix_mz, 
                                        mix_rt, mix_tic, 0, 0, 'neutral',
                                        len(neg_ions) + len(pos_ions), "MIX", mix_samples]

        # add a record to the transitions table
        transitions_table.append((new_idx, neutral_neg, neutral_pos))
    transitions_table = pd.DataFrame(transitions_table, columns = ['MIX', 'NEG', 'POS'])   
    
    
    # Drop neutrals that were merged and redo the process for the remaining neutrals
    # that might have escaped, this time in POS to NEG.
    neg_neutrals.drop(transitions_table['NEG'], inplace = True)
    pos_neutrals.drop(transitions_table['POS'], inplace = True)
    merged_neutrals = []

    for neutral in tqdm(pos_neutrals.index): # POS -> NEG
    
        # Get neutral attributes
        mz_pos = pos_neutrals.loc[neutral, mz_field]
        rt_pos = pos_neutrals.loc[neutral, rt_field]
        samples_pos = set(pos_neutrals.loc[neutral, "samples"].split('|'))

        # Find neutral in opposite mode with similar mz and RT
        counter_neg = neg_neutrals[neg_neutrals[mz_field].between(mz_pos - mass_error, mz_pos + mass_error, inclusive = "both")]
        counter_neg = counter_neg[counter_neg[rt_field].between(rt_pos - rt_error, rt_pos + rt_error, inclusive = "both")] 
        if len(counter_neg) == 0 : continue
    
        # Add a d_rt column
        counter_neg['d_rt'] = abs(counter_neg[rt_field]-rt_pos)

        # Find samples shared between neutrals
        shared_samples = []
        for i in counter_neg.index:
            samples_neg = counter_neg.loc[i, "samples"].split('|')
            shared_samples.append(len(samples_pos.intersection(samples_neg)))
 
        # Filter by sample count
        counter_neg['shared_samples'] = shared_samples
        counter_neg = counter_neg[counter_neg['shared_samples'] > 0]
        
        # Select by max shared samples, and then by lowest d_rt
        counter_neg = counter_neg[counter_neg['shared_samples'] == counter_neg['shared_samples'].max()]
        counter_neg = counter_neg[counter_neg['d_rt'] == counter_neg['d_rt'].min()]

        # Stuff that shouldn't happen, but just in case
        if len(counter_neg) > 1 : sys.exit("MORE THAN 1 COUTNER NEUTRAL AT SAME RT")

        # If no neutrals survived selection, skip
        elif len(counter_neg) == 0 : continue

        # If match found
        else:
            
            # Get attributes and report to the merged neutrals table
            mz_neg = counter_neg[mz_field].iloc[0]
            rt_neg = counter_neg[rt_field].iloc[0]
            merged_neutrals.append((counter_neg.index[0], mz_neg, rt_neg, neutral, mz_pos, rt_pos))

    # Create a dataframe with merged neutrals
    merged_neutrals = pd.DataFrame(merged_neutrals, columns = ['neg_neutral', 'mass_neg',
                                                               'rt_neg', 'pos_neutral',
                                                               'mass_pos', 'rt_pos'])
    merged_neutrals['d_rt'] = abs(merged_neutrals['rt_pos'] - merged_neutrals['rt_neg'])
    
    # Same as before, if there are several POS neutrals per NEG neutral, choose the
    # closest one by RT
    for neutral in merged_neutrals['pos_neutral'].unique():
        tmp_table = merged_neutrals[merged_neutrals['pos_neutral'] == neutral].copy()
        if len(tmp_table) > 1 : 
            tmp_table.sort_values('d_rt', inplace = True)
            for i in tmp_table.index[1:]:
                merged_neutrals.drop(i, inplace = True)
    
    # Add the new merged neutrals to the transitions table
    print("Adding merged neutrals to the merged node table...")
    for i in tqdm(merged_neutrals.index):
        
        # Get merged neutral new index
        new_idx = merged_node_table.index.max() + 1 
        
        # Get the old indexes
        neutral_neg = merged_neutrals.loc[i, "neg_neutral"]
        neutral_pos = merged_neutrals.loc[i, "pos_neutral"]
        
        # Get the neg ions connected
        neg_ions = list(edge_table_all_neg['node_2'][edge_table_all_neg['node_1'] == neutral_neg])

        # Get neg neutral attributes
        neg_count = len(neg_ions)
        neg_mz = node_table_all_neg.loc[neutral_neg, mz_field]
        neg_rt = node_table_all_neg.loc[neutral_neg, rt_field]
        neg_tic = node_table_all_neg.loc[neutral_neg, "TIC"]
        neg_samples = node_table_all_neg.loc[neutral_neg, "samples"].split('|')
        
        # Get the pos ions connected
        pos_ions = list(edge_table_all_pos['node_2'][edge_table_all_pos['node_1'] == neutral_pos])

        # Get pos neutral attributes
        pos_count = len(pos_ions)
        pos_mz = node_table_all_pos.loc[neutral_pos, "mz"]
        pos_rt = node_table_all_pos.loc[neutral_pos, "rt"]
        pos_tic = node_table_all_pos.loc[neutral_pos, "TIC"]
        pos_samples = node_table_all_pos.loc[neutral_pos, "samples"].split('|')

        # Produce combined neutral attributes
        mix_mz = round(((neg_mz*neg_count) + (pos_mz*pos_count))/(neg_count + pos_count), 4)
        mix_rt = round(((neg_rt*neg_count) + (pos_rt*pos_count))/(neg_count + pos_count), 3)
        mix_tic = neg_tic + pos_tic
        mix_samples = list(set(neg_samples + pos_samples))
        mix_samples.sort()
        mix_samples = '|'.join(mix_samples)

        # Report data to the merged node table
        merged_node_table.loc[new_idx] = [None]*len(merged_node_table.columns)
        
        merged_node_table.loc[new_idx, [mz_field, rt_field, "TIC", "rule_points",
                                        charge_field, "status", "adduct_count",
                                        'ion_mode', 'samples']] = [mix_mz, 
                                        mix_rt, mix_tic, 0, 0, 'neutral',
                                        len(neg_ions) + len(pos_ions), "MIX", mix_samples]
                    
        # Add record to the transitions table
        transitions_table.loc[transitions_table.index.max() + 1] = [new_idx, neutral_neg, neutral_pos]


    # Remove merged neutrals from neg and pos neutral lists
    intersect_neutrals = set(neg_neutrals.index)
    intersect_neutrals = intersect_neutrals.intersection(transitions_table['NEG'])
    neg_neutrals.drop(intersect_neutrals, inplace = True)
    
    intersect_neutrals = set(pos_neutrals.index)
    intersect_neutrals = intersect_neutrals.intersection(transitions_table['POS'])    
    pos_neutrals.drop(intersect_neutrals, inplace = True)
    
    # Add non merged neg neutrals to the merged node table
    print('Adding non-merged neg neutrals to the merged node table...')
    for neutral in tqdm(neg_neutrals.index):
        new_idx = merged_node_table.index.max() + 1 
        merged_node_table.loc[new_idx] = neg_neutrals.loc[neutral].copy()
        new_row = pd.Series([new_idx, neutral, None], index = ["MIX", "NEG", "POS"])
        transitions_table = transitions_table.append(new_row, ignore_index = True)

    # Add non merged pos neutrals to the merged node table
    print('Adding non-merged pos neutrals to the merged node table...')
    for neutral in tqdm(pos_neutrals.index):
        new_idx = merged_node_table.index.max() + 1 
        merged_node_table.loc[new_idx] = pos_neutrals.loc[neutral].copy()
        new_row = pd.Series([new_idx, None, neutral], index = ["MIX", "NEG", "POS"])
        transitions_table = transitions_table.append(new_row, ignore_index = True)
    
    # Replace old IDs by new ones in the merged edge table.
    print('Resetting edge table IDs...')
    for i in tqdm(merged_edge_table.index):
        ion_mode = merged_edge_table.loc[i, "ion_mode"]
        if ion_mode == "NEG" : opposed_mode = "POS"
        else : opposed_mode = "NEG"
        
        # Get node table of same mode as edge table
        merged_node_table_sub = merged_node_table[merged_node_table['ion_mode'] != opposed_mode]

        # If neutral-adduct edge
        if merged_edge_table.loc[i, "status"] == "add_edge":
            # Get old IDs, replace by new
            node_1 = merged_edge_table.loc[i, "node_1"]
            node_2 = merged_edge_table.loc[i, "node_2"]
            
            # Node 1 : get new ID from transitions table
            new_node_1 = int(transitions_table["MIX"][transitions_table[ion_mode] == node_1].iloc[0])
            new_node_2 = merged_node_table_sub.index[merged_node_table_sub[idx_column] == node_2][0]
            merged_edge_table.loc[i, "node_1"] = new_node_1
            merged_edge_table.loc[i, "node_2"] = new_node_2

        # If ion-ion edge
        else:
            # Get old IDs, replace by new
            node_1 = merged_edge_table.loc[i, "node_1"]
            node_2 = merged_edge_table.loc[i, "node_2"]
            new_node_1 = merged_node_table_sub.index[merged_node_table_sub[idx_column] == node_1][0]
            new_node_2 = merged_node_table_sub.index[merged_node_table_sub[idx_column] == node_2][0]
            merged_edge_table.loc[i, "node_1"] = new_node_1
            merged_edge_table.loc[i, "node_2"] = new_node_2
    
    # Search for remaining ion annotations
    
    # Remains table pos : non neutral, non adduct ions
    remains_table_pos = merged_node_table[merged_node_table['status'] != "neutral"]
    remains_table_pos = remains_table_pos[remains_table_pos['status'] != "adduct"] 
    remains_table_pos = remains_table_pos[remains_table_pos['ion_mode'] == "POS"]

    # Remains table neg, same as above
    remains_table_neg = merged_node_table[merged_node_table['status'] != "neutral"]
    remains_table_neg = remains_table_neg[remains_table_neg['status'] != "adduct"] 
    remains_table_neg = remains_table_neg[remains_table_neg['ion_mode'] == "NEG"]

    # Get positive neutral node table to search for neg ions
    node_table_neutrals_pos = merged_node_table[merged_node_table["status"] == "neutral"]
    node_table_neutrals_pos = node_table_neutrals_pos[node_table_neutrals_pos["ion_mode"] == "POS"]

    # Get negative neutral node table to search for pos ions
    node_table_neutrals_neg = merged_node_table[merged_node_table["status"] == "neutral"]
    node_table_neutrals_neg = node_table_neutrals_neg[node_table_neutrals_neg["ion_mode"] == "NEG"]      

    # Link neg neutrals to pos remains
    print('Linking NEG neutrals to single POS ions...')
    candidates_table = list()
    for i in tqdm(node_table_neutrals_neg.index):
        
        # Get neutral attributes and samples
        mol_mass = node_table_neutrals_neg.loc[i, mz_field]
        mol_rt = node_table_neutrals_neg.loc[i, rt_field]
        mol_samples = set(node_table_neutrals_neg.loc[i, "samples"].split('|'))

        # Get ions that could be produced by the neutral
        ion_table = get_ion_table(mol_mass, adduct_table_merged_pos)
        ion_hits = list()
        hit_table_rt = remains_table_pos[remains_table_pos[rt_field].between(mol_rt - rt_error, mol_rt + rt_error, inclusive = "both")].copy()
        
        # Get shared samples
        shares_samples_list = list()
        for j in hit_table_rt.index:
            ion_samples = hit_table_rt.loc[j, "samples"].split('|')
            shared_samples = len(list(mol_samples.intersection(ion_samples)))
            shares_samples_list.append(shared_samples)
 
        # Filter by shared samples
        hit_table_rt['shared_samples'] = shares_samples_list
        hit_table_rt = hit_table_rt[hit_table_rt['shared_samples'] > 0]
        
        # Check if some coeluted ions match with the ion table
        for j in ion_table.index:
            ion_mz = ion_table.loc[j, "ion_mz"]
            hit_table = hit_table_rt[hit_table_rt[mz_field].between(ion_mz - mass_error, ion_mz + mass_error, inclusive = "both")].copy()
            if len(hit_table) >0 :
                ion_hits.append('|'.join(hit_table.index.astype(str)))
            else:
                ion_hits.append(None)

        # Add and keep only hits to the ion table
        ion_table['ion_hits'] = ion_hits
        ion_table = ion_table[~ion_table['ion_hits'].isnull()]

        # Add matches to the candidates tables
        for j in ion_table.index:
            adduct = ion_table.loc[j, "Adduct"]
            for k in ion_table.loc[j, "ion_hits"].split('|'):
                k = int(k)
                d_rt = abs(mol_rt - remains_table_pos.loc[k, rt_field])
                d_mz = abs(mol_mass - remains_table_pos.loc[k, mz_field])
                candidates_table.append((k, i, adduct, d_rt, d_mz))
                
    # Create candidates table dataframe
    candidates_table = pd.DataFrame(candidates_table, columns = ['ion_idx', 'neutral', 'adduct', 'd_rt', 'd_mz'])
    
    # Migration of ions:
    unique_hits = candidates_table['ion_idx'].unique().tolist()
    selected_neutrals = list()
    selected_adducts = list()
    delta_rts = list()
    delta_mzs = list()
    for i in unique_hits:
        tmp_table = candidates_table[candidates_table['ion_idx'] == i]
        selected_neutral = tmp_table['d_rt'].idxmin()
        selected_neutrals.append(tmp_table.loc[selected_neutral, "neutral"])
        selected_adducts.append(tmp_table.loc[selected_neutral, "adduct"])
        delta_rts.append(tmp_table.loc[selected_neutral, "d_rt"])
        delta_mzs.append(tmp_table.loc[selected_neutral, "d_mz"])
    candidates_table= list(zip(unique_hits, selected_neutrals, selected_adducts, delta_rts, delta_mzs))
    candidates_table = pd.DataFrame(candidates_table, columns = ['ion_idx', 'neutral', 'adduct', 'd_rt', 'd_mz'])
    
    # Report results and node and edge tables:
    print('Updating node and edge tables...')
    for i in tqdm(candidates_table.index):
        # Retrieve data for ion and neutral
        ion_idx = candidates_table.loc[i, "ion_idx"]
        neutral_idx = candidates_table.loc[i, "neutral"]
        adduct = candidates_table.loc[i, "adduct"]
        d_rt = candidates_table.loc[i, "d_rt"]
        d_mz = candidates_table.loc[i, "d_mz"]
        spec_id = int(merged_node_table.loc[ion_idx, "spec_id"])
        adduct_code = adduct_table_merged_pos['Adduct_code'][adduct_table_merged_pos['Adduct'] == adduct].iloc[0]
        Species_rule = Validator_choice(adduct_code, "POS")
        rule_points = Species_rule(prec_mass_error, spec_id, pos_specturm_list)
        ion_tic = merged_node_table.loc[ion_idx, 'TIC']
        
        # Update node table with ion data
        merged_node_table.loc[ion_idx, 'status'] = "adduct"
        merged_node_table.loc[ion_idx, 'Adnotation'] = adduct
        merged_node_table.loc[ion_idx, 'rule_points'] = rule_points
        
        # Update node table with neutral data
        merged_node_table.loc[neutral_idx, "TIC"] += ion_tic
        merged_node_table.loc[neutral_idx, "adduct_count"] += 1
        merged_node_table.loc[neutral_idx, "ion_mode"] = "MIX"
        tmp_samples = merged_node_table.loc[neutral_idx, "samples"].split('|')
        tmp_samples += merged_node_table.loc[ion_idx, "samples"].split('|')
        tmp_samples = list(set(tmp_samples))
        tmp_samples.sort()
        tmp_samples = '|'.join(tmp_samples)
        merged_node_table.loc[neutral_idx, "samples"] = tmp_samples

        # Update edge table :
        del_edge = merged_edge_table[merged_edge_table['node_1'] == ion_idx]
        del_edge = del_edge.index[del_edge['node_2'] == ion_idx]
        if len(del_edge) > 0 : merged_edge_table.drop(del_edge[0], inplace = True)
        new_idx = merged_edge_table.index.max() + 1 
        merged_edge_table.loc[new_idx] = [None]*len(merged_edge_table.columns)
        merged_edge_table.loc[new_idx, ["node_1", "node_2", "matched_peaks",
                                        "total_peaks", "matching_score", "rt_gap",
                                        "mz_gap", "status", "Fragnotation", 
                                        "Adnotation", "All_annotations",
                                        "ion_mode"]] = [neutral_idx, ion_idx,
                                                        0, 0, 0, d_rt, d_mz,
                                                        "add_edge", None, adduct,
                                                        adduct, "POS"]
    
    # Link POS neutrals to NEG ions
    print('Linking POS neutrals to single NEG ions...')
    candidates_table = list()
    for i in tqdm(node_table_neutrals_pos.index):
        
        # Get neutral attributes
        mol_mass = node_table_neutrals_pos.loc[i, mz_field]
        mol_rt = node_table_neutrals_pos.loc[i, rt_field]
        mol_samples = set(node_table_neutrals_pos.loc[i, "samples"].split('|'))

        # Get ion table
        ion_table = get_ion_table(mol_mass, adduct_table_merged_neg)

        # Find ion hits
        ion_hits = list()
        hit_table_rt = remains_table_neg[remains_table_neg[rt_field].between(mol_rt - rt_error, mol_rt + rt_error, inclusive = "both")].copy()
        
        # Find shared samples
        shares_samples_list = list()
        for j in hit_table_rt.index:
            ion_samples = hit_table_rt.loc[j, "samples"].split('|')
            shared_samples = len(list(mol_samples.intersection(ion_samples)))
            shares_samples_list.append(shared_samples)

        # Filter by shared samples
        hit_table_rt['shared_samples'] = shares_samples_list
        hit_table_rt = hit_table_rt[hit_table_rt['shared_samples'] > 0]
        
        # Find hits based on ion table possibilities and available ion m/z values
        for j in ion_table.index:
            ion_mz = ion_table.loc[j, "ion_mz"]
            hit_table = hit_table_rt[hit_table_rt[mz_field].between(ion_mz - mass_error, ion_mz + mass_error, inclusive = "both")].copy()
            if len(hit_table) >0 :
                ion_hits.append('|'.join(hit_table.index.astype(str)))
            else:
                ion_hits.append(None)

        # Report hits to the ion table and keep only ions with hits
        ion_table['ion_hits'] = ion_hits
        ion_table = ion_table[~ion_table['ion_hits'].isnull()]

        # Add hits to the candidates tables
        for j in ion_table.index:
            adduct = ion_table.loc[j, "Adduct"]
            for k in ion_table.loc[j, "ion_hits"].split('|'):
                k = int(k)
                d_rt = abs(mol_rt - remains_table_neg.loc[k, rt_field])
                d_mz = abs(mol_mass - remains_table_neg.loc[k, mz_field])
                candidates_table.append((k, i, adduct, d_rt, d_mz))

    # candidates table to dataframe
    candidates_table = pd.DataFrame(candidates_table, columns = ['ion_idx', 'neutral', 'adduct', 'd_rt', 'd_mz'])
    
    # Migration of ions:
    unique_hits = candidates_table['ion_idx'].unique().tolist()
    selected_neutrals = list()
    selected_adducts = list()
    delta_rts = list()
    delta_mzs = list()
    for i in unique_hits:
        tmp_table = candidates_table[candidates_table['ion_idx'] == i]
        selected_neutral = tmp_table['d_rt'].idxmin()
        selected_neutrals.append(tmp_table.loc[selected_neutral, "neutral"])
        selected_adducts.append(tmp_table.loc[selected_neutral, "adduct"])
        delta_rts.append(tmp_table.loc[selected_neutral, "d_rt"])
        delta_mzs.append(tmp_table.loc[selected_neutral, "d_mz"])
    candidates_table= list(zip(unique_hits, selected_neutrals, selected_adducts, delta_rts, delta_mzs))
    candidates_table = pd.DataFrame(candidates_table, columns = ['ion_idx', 'neutral', 'adduct', 'd_rt', 'd_mz'])
    
    # Report results and node and edge tables:
    print('Updating node and edge tables...')
    for i in tqdm(candidates_table.index):
        # Retrieve data for ion and neutral
        ion_idx = candidates_table.loc[i, "ion_idx"]
        neutral_idx = candidates_table.loc[i, "neutral"]
        adduct = candidates_table.loc[i, "adduct"]
        d_rt = candidates_table.loc[i, "d_rt"]
        d_mz = candidates_table.loc[i, "d_mz"]
        spec_id = int(merged_node_table.loc[ion_idx, "spec_id"])
        adduct_code = adduct_table_merged_neg['Adduct_code'][adduct_table_merged_neg['Adduct'] == adduct].iloc[0]
        Species_rule = Validator_choice(adduct_code, "NEG")
        rule_points = Species_rule(prec_mass_error, spec_id, neg_spectrum_list)
        ion_tic = merged_node_table.loc[ion_idx, 'TIC']
        
        # Update node table with ion data
        merged_node_table.loc[ion_idx, 'status'] = "adduct"
        merged_node_table.loc[ion_idx, 'Adnotation'] = adduct
        merged_node_table.loc[ion_idx, 'rule_points'] = rule_points
        
        # Update node table with neutral data
        merged_node_table.loc[neutral_idx, "TIC"] += ion_tic
        merged_node_table.loc[neutral_idx, "adduct_count"] += 1
        merged_node_table.loc[neutral_idx, "ion_mode"] = "MIX"
        tmp_samples = merged_node_table.loc[neutral_idx, "samples"].split('|')
        tmp_samples += merged_node_table.loc[ion_idx, "samples"].split('|')
        tmp_samples = list(set(tmp_samples))
        tmp_samples.sort()
        tmp_samples = '|'.join(tmp_samples)
        merged_node_table.loc[neutral_idx, "samples"] = tmp_samples

        # Update edge table :
        del_edge = merged_edge_table[merged_edge_table['node_1'] == ion_idx]
        del_edge = del_edge.index[del_edge['node_2'] == ion_idx]
        if len(del_edge) > 0 : merged_edge_table.drop(del_edge[0], inplace = True)
        new_idx = merged_edge_table.index.max() + 1 
        merged_edge_table.loc[new_idx] = [None]*len(merged_edge_table.columns)
        merged_edge_table.loc[new_idx, ["node_1", "node_2", "matched_peaks",
                                        "total_peaks", "matching_score", "rt_gap",
                                        "mz_gap", "status", "Fragnotation", 
                                        "Adnotation", "All_annotations",
                                        "ion_mode"]] = [neutral_idx, ion_idx,
                                                        0, 0, 0, d_rt, d_mz,
                                                        "add_edge", None, adduct,
                                                        adduct, "NEG"]
    
    # Produce cluster IDs. Making these is essential before opposed mode singletons
    # search, to filter out singletons already present in molecular clusters.
    node_pool = list(merged_node_table.index)
    singletons = list(merged_edge_table["node_1"][merged_edge_table['status'] == "self_edge"])
    node_pool = list(set(node_pool) - set(singletons))
    cluster_list = []
    cluster_size_list = []
    total_nodes = len(node_pool)
    while len(node_pool) > 0:
        new_cluster = [node_pool[0]]
        cluster_size = 0
        perc = round((1-(len(node_pool)/total_nodes))*100,1)
        sys.stdout.write("\rDefining new clusters : {0}%".format(perc))
        sys.stdout.flush()
        while cluster_size != len(new_cluster):
            cluster_size = len(new_cluster)
            tmp_idx = []
            for i in new_cluster:
                tmp_idx += list(merged_edge_table.index[merged_edge_table['node_1'] == i])
                tmp_idx += list(merged_edge_table.index[merged_edge_table['node_2'] == i])
            new_cluster += list(merged_edge_table.loc[tmp_idx, 'node_1'].astype(int))
            new_cluster += list(merged_edge_table.loc[tmp_idx, 'node_2'].astype(int))
            new_cluster = list(set(new_cluster))
        new_cluster.sort()
        node_pool = list(set(node_pool) - set(new_cluster))
        cluster_size_list.append(len(new_cluster))
        cluster_list.append('|'.join(list(map(str, new_cluster))))

    cluster_table= pd.DataFrame()
    cluster_table['cluster'] = cluster_list
    cluster_table['cluster_size'] = cluster_size_list
    cluster_table.sort_values('cluster_size', ascending = False, inplace = True)
    cluster_table.reset_index(drop = True, inplace = True)
    
    # Identify molecular clusters
    cluster_molecular = list()
    for i in cluster_table.index:
        node_list = cluster_table.loc[i, "cluster"].split('|')
        node_list = list(map(int, node_list))
        tmp_table_1 = merged_node_table.loc[node_list]
        if sum(tmp_table_1['status'] == "neutral") > 0 :
            cluster_molecular.append(True)
        else:
            cluster_molecular.append(False)
    cluster_table["molecular_cluster"] = cluster_molecular
    
    # Adding new cluster indexes to the merged node table
    merged_node_table['cluster_id'] = [-1]*len(merged_node_table)
    print('Assigning new cluster indexes...')
    for i in tqdm(cluster_table.index):
        node_list = list(map(int, cluster_table.loc[i, 'cluster'].split('|')))
        for j in node_list :
            merged_node_table.loc[j, 'cluster_id'] = i
    
    # Connect singletons to singletons (only precursors and nodes from non molecular clusters)
    remains_table = cluster_table.index[~cluster_table["molecular_cluster"]].tolist()
    remains_table = [merged_node_table.index[merged_node_table['cluster_id'] == i].tolist() for i in remains_table]
    remains_table = list(flatten(remains_table)) # Added precursors and fragments from non-molecular clusters
    remains_table += merged_node_table.index[merged_node_table['cluster_id'] == -1].tolist() # Added singleton nodes
    precursor_ions = list(set(merged_node_table.index) - set(remains_table))
    precursor_ions = merged_node_table.loc[precursor_ions]
    precursor_ions = precursor_ions.index[precursor_ions['status'] == "precursor"].tolist()
    remains_table += precursor_ions # Added precursors from molecular clusters
    remains_table.sort()
    remains_table = merged_node_table.loc[remains_table]
    remains_table_pos = remains_table[remains_table['ion_mode'] == "POS"]
    remains_table_neg = remains_table[remains_table['ion_mode'] == "NEG"]

    # Last step, connecting singletons to singletons
    print('Linking NEG singletons to POS singletons...')
    candidates_table = list()
    for i in tqdm(remains_table_neg.index):
        
        # Get ion attributes
        ion_1_mz = remains_table_neg.loc[i, mz_field]
        ion_1_rt = remains_table_neg.loc[i, rt_field]
        ion_1_samples = set(remains_table_neg.loc[i, "samples"].split('|'))
        
        # Produce neutral table and hit table based on these neutrals
        neutral_table = neutral_tabler(ion_1_mz, adduct_table_base_neg)
        hit_table_rt = remains_table_pos[remains_table_pos[rt_field].between(ion_1_rt - rt_error, ion_1_rt + rt_error, inclusive = "both")].copy()

        # Find shared samples
        shares_samples_list = list()
        for j in hit_table_rt.index:
            ion_samples = hit_table_rt.loc[j, "samples"].split('|')
            shared_samples = len(list(ion_1_samples.intersection(ion_samples)))
            shares_samples_list.append(shared_samples)

        # Filter hit table by shared samples
        hit_table_rt['shared_samples'] = shares_samples_list
        hit_table_rt = hit_table_rt[hit_table_rt['shared_samples'] > 0]
        if len(hit_table_rt) == 0 : continue
    
        # Filter by mz
        for j in neutral_table.index:
            ion_table = get_ion_table(neutral_table.loc[j, "neutral_mass"], adduct_table_base_pos)
            for k in ion_table.index:
                hit_table = hit_table_rt[hit_table_rt[mz_field].between(ion_table.loc[k, "ion_mz"] - mass_error, ion_table.loc[k, "ion_mz"] + mass_error, inclusive = "both")]
                for l in hit_table.index:
                    candidates_table.append((i,
                                             neutral_table.loc[j, "Adduct"],
                                             neutral_table.loc[j, "Complexity"],
                                             l,
                                             ion_table.loc[k, "Adduct"],
                                             ion_table.loc[k, "Complexity"]))
    # candidates table to dataframe
    candidates_table = pd.DataFrame(candidates_table, columns = ['neg_ion', 'neg_adduct', 'neg_complexity', 'pos_ion', 'pos_adduct', 'pos_complexity'])

    # Get the best hypotheses
    unique_negs = candidates_table['neg_ion'].unique().tolist()
    tmp_candidates_table = list()
    print('Selecting best NEG annotations...')
    for i in tqdm(unique_negs) :
        tmp_table_1 = candidates_table[candidates_table['neg_ion'] == i]
        unique_adducts = tmp_table_1['neg_adduct'].unique().tolist()
        point_list = list()
        for adduct in unique_adducts:
            tmp_table_2 = tmp_table_1[tmp_table_1['neg_adduct'] == adduct]
            points = 2*len(tmp_table_2) / (tmp_table_2['neg_complexity'].sum() + tmp_table_2['pos_complexity'].sum())
            point_list.append(points)
        selected_adduct = unique_adducts[point_list.index(max(point_list))]
        tmp_table_1 = tmp_table_1[tmp_table_1['neg_adduct'] == selected_adduct]
        for j in tmp_table_1.index:
            tmp_candidates_table.append((tmp_table_1.loc[j, "neg_ion"],
                                     tmp_table_1.loc[j, "neg_adduct"],
                                     tmp_table_1.loc[j, "neg_complexity"],
                                     tmp_table_1.loc[j, "pos_ion"],
                                     tmp_table_1.loc[j, "pos_adduct"],
                                     tmp_table_1.loc[j, "pos_complexity"]))
    candidates_table = pd.DataFrame(tmp_candidates_table, columns = ['neg_ion', 'neg_adduct', 'neg_complexity', 'pos_ion', 'pos_adduct', 'pos_complexity'])

    # Resolve POS ions with multiple annotations
    unique_pos = candidates_table['pos_ion'].unique().tolist()
    tmp_candidates_table = list()
    print('Selecting best POS annotations...')
    for i in tqdm(unique_pos):
        tmp_table_1 = candidates_table[candidates_table['pos_ion'] == i]
        unique_adducts = tmp_table_1['pos_adduct'].unique().tolist()
        point_list = list()
        for adduct in unique_adducts:
            tmp_table_2 = tmp_table_1[tmp_table_1['pos_adduct'] == adduct]
            points = 2*len(tmp_table_2) / (tmp_table_2['neg_complexity'].sum() + tmp_table_2['pos_complexity'].sum())
            point_list.append(points)
        selected_adduct = unique_adducts[point_list.index(max(point_list))]
        tmp_table_1 = tmp_table_1[tmp_table_1['pos_adduct'] == selected_adduct]
        for j in tmp_table_1.index:
            tmp_candidates_table.append((tmp_table_1.loc[j, "neg_ion"], 
                                         tmp_table_1.loc[j, "neg_adduct"],
                                         tmp_table_1.loc[j, "neg_complexity"],
                                         tmp_table_1.loc[j, "pos_ion"],
                                         tmp_table_1.loc[j, "pos_adduct"],
                                         tmp_table_1.loc[j, "pos_complexity"]))      
    tmp_candidates_table = pd.DataFrame(tmp_candidates_table, columns = ['neg_ion', 'neg_adduct', 'neg_complexity', 'pos_ion', 'pos_adduct', 'pos_complexity'])
    candidates_table = candidates_table.append(tmp_candidates_table, ignore_index = True)

    # Filter candidates table to have only unique combinations
    candidates_table['mix_idx'] = candidates_table['neg_ion'].astype(str) + "&" + candidates_table['pos_ion'].astype(str)
    tmp_candidates_table = list()
    for comb in candidates_table['mix_idx'].unique():
        tmp_candidates_table.append((candidates_table[candidates_table['mix_idx'] == comb].iloc[0].tolist()))    
    candidates_table = pd.DataFrame(tmp_candidates_table, columns = ['neg_ion', 'neg_adduct', 'neg_complexity', 'pos_ion', 'pos_adduct', 'pos_complexity', 'mix_idx'])    

    # Make a neutral table:
    neutral_idx = 0
    neutral_table = list()
    while len(candidates_table) > 0 :
        idx = candidates_table.index[0]
        neg_pool = list()
        pos_pool = list()
        idx_pool = list()
        neg_pool.append(candidates_table.loc[idx, "neg_ion"])
        pos_pool.append(candidates_table.loc[idx, "pos_ion"])
        new_neg_pool = neg_pool
        new_pos_pool = pos_pool
        idx_pool.append(idx)
        len_0 = 0
        while len_0 != len(idx_pool):
            len_0 = len(idx_pool)
            for i in new_neg_pool:
                idx_pool += candidates_table.index[candidates_table['neg_ion'] == i].tolist()
            for i in new_pos_pool:
                idx_pool += candidates_table.index[candidates_table['pos_ion'] == i].tolist()
            idx_pool = list(set(idx_pool))
        tmp_table_1 = candidates_table.loc[idx_pool]
        
        for i in tmp_table_1['neg_ion'].unique():
            idx = tmp_table_1.index[tmp_table_1['neg_ion'] == i][0]
            spec_id = int(merged_node_table.loc[i, "spec_id"])
            ion_mz = merged_node_table.loc[i, mz_field]
            ion_rt = merged_node_table.loc[i, rt_field]
            ion_tic = merged_node_table.loc[i, "TIC"]
            ion_adduct = tmp_table_1.loc[idx, "neg_adduct"]
            ion_samples = merged_node_table.loc[i, "samples"]
            adduct_idx = adduct_table_base_neg.index[adduct_table_base_neg["Adduct"] == ion_adduct][0]
            ion_adduct_code = adduct_table_base_neg.loc[adduct_idx, "Adduct_code"]
            mol_mass = neutral_mass_calculator(ion_mz,
                                adduct_table_base_neg.loc[adduct_idx, "Adduct_mass"],
                                adduct_table_base_neg.loc[adduct_idx, "Mol_multiplier"],
                                adduct_table_base_neg.loc[adduct_idx, "Charge"])
            Species_rule = Validator_choice(ion_adduct_code, "NEG")
            ion_rule_points = Species_rule(prec_mass_error, spec_id, neg_spectrum_list)
            neutral_table.append((neutral_idx, mol_mass, ion_mz, ion_rt, i, ion_adduct, ion_tic, ion_rule_points, ion_samples, "NEG"))
            
        for i in tmp_table_1['pos_ion'].unique():
            idx = tmp_table_1.index[tmp_table_1['pos_ion'] == i][0]
            spec_id = int(merged_node_table.loc[i, "spec_id"])
            ion_mz = merged_node_table.loc[i, mz_field]
            ion_rt = merged_node_table.loc[i, rt_field]
            ion_tic = merged_node_table.loc[i, "TIC"]
            ion_adduct = tmp_table_1.loc[idx, "pos_adduct"]
            ion_samples = merged_node_table.loc[i, "samples"]
            adduct_idx = adduct_table_base_pos.index[adduct_table_base_pos["Adduct"] == ion_adduct][0]
            ion_adduct_code = adduct_table_base_pos.loc[adduct_idx, "Adduct_code"]
            mol_mass = neutral_mass_calculator(ion_mz,
                                adduct_table_base_pos.loc[adduct_idx, "Adduct_mass"],
                                adduct_table_base_pos.loc[adduct_idx, "Mol_multiplier"],
                                adduct_table_base_pos.loc[adduct_idx, "Charge"])
            Species_rule = Validator_choice(ion_adduct_code, "POS")
            ion_rule_points = Species_rule(prec_mass_error, spec_id, pos_specturm_list)
            neutral_table.append((neutral_idx, mol_mass, ion_mz,  ion_rt, i, ion_adduct, ion_tic, ion_rule_points, ion_samples, "POS"))            
        
        candidates_table.drop(idx_pool, inplace = True)
        neutral_idx += 1
    
    neutral_table = pd.DataFrame(neutral_table, columns = ['neutral_idx', 'neutral_mass', 'ion_mz', 'ion_rt', 'ion_idx', 'adduct', 'TIC', 'rule_points', 'samples', 'ion_mode'])         

    # Report the results: 
    print('Reporting results for NEG singletons paired to POS singletons...')
    for i in tqdm(neutral_table['neutral_idx'].unique()):
        tmp_table_1 = neutral_table[neutral_table['neutral_idx'] == i].copy()
        mol_mass = tmp_table_1['neutral_mass'].mean()
        mol_rt = tmp_table_1['ion_rt'].mean()
        mol_tic = tmp_table_1['TIC'].sum()
        mol_cluster = merged_node_table['cluster_id'].max() + 1
        adduct_count = len(tmp_table_1)
        tmp_table_1['mz_gap'] = abs(mol_mass - tmp_table_1['ion_mz'])
        tmp_table_1['rt_gap'] = abs(mol_rt - tmp_table_1['ion_rt'])
        mol_samples = '|'.join(tmp_table_1['samples'])
        mol_samples = list(set(mol_samples.split('|')))
        mol_samples.sort()
        mol_samples = '|'.join(mol_samples)
        new_idx = merged_node_table.index.max() + 1
        merged_node_table.loc[new_idx] = [None]*len(merged_node_table.columns)
        merged_node_table.loc[new_idx, [mz_field, rt_field, "TIC", "rule_points",
                                        'status', 'Adnotation', 'adduct_count',
                                        'ion_mode', 'samples', 'cluster_id']] = [mol_mass,
                                        mol_rt, mol_tic, 0, "neutral", None, adduct_count,
                                        "MIX", mol_samples, mol_cluster]
        
        for j in tmp_table_1.index:
            ion_idx = tmp_table_1.loc[j, 'ion_idx']
            adduct = tmp_table_1.loc[j, 'adduct']
            rule_points = tmp_table_1.loc[j, 'rule_points']
            mz_gap = tmp_table_1.loc[j, 'mz_gap']
            rt_gap = tmp_table_1.loc[j, 'rt_gap']
            tmp_mode = tmp_table_1.loc[j, "ion_mode"]
            merged_node_table.loc[ion_idx, 'rule_points'] = rule_points
            merged_node_table.loc[ion_idx, 'status'] = "adduct"
            merged_node_table.loc[ion_idx, 'Adnotation'] = adduct
            merged_node_table.loc[ion_idx, 'cluster_id'] = mol_cluster
            new_edge = merged_edge_table.index.max() + 1
            merged_edge_table.loc[new_edge] = [None]*len(merged_edge_table.columns)
            merged_edge_table.loc[new_edge, ['node_1', 'node_2', 'matched_peaks',
                                            'total_peaks', 'matching_score',
                   'rt_gap', 'mz_gap', 'status', 'Fragnotation', 'Adnotation',
                   'All_annotations', 'ion_mode']] = [new_idx, ion_idx, 0 ,0 ,0 ,
                    rt_gap, mz_gap, "add_edge", None, adduct, adduct, tmp_mode]

            del_edge = merged_edge_table[merged_edge_table['node_1'] == ion_idx]
            del_edge = del_edge[del_edge['node_2'] == ion_idx].index
            if len(del_edge) > 0 :
                merged_edge_table.drop(del_edge[0], inplace = True)

            


    # Connect singletons to singletons (only precursors and nodes from non molecular clusters)
    remains_table = cluster_table.index[~cluster_table["molecular_cluster"]].tolist()
    remains_table = [merged_node_table.index[merged_node_table['cluster_id'] == i].tolist() for i in remains_table]
    remains_table = list(flatten(remains_table)) # Added precursors and fragments from non-molecular clusters
    remains_table += merged_node_table.index[merged_node_table['cluster_id'] == -1].tolist() # Added singleton nodes
    precursor_ions = list(set(merged_node_table.index) - set(remains_table))
    precursor_ions = merged_node_table.loc[precursor_ions]
    precursor_ions = precursor_ions.index[precursor_ions['status'] == "precursor"].tolist()
    remains_table += precursor_ions # Added precursors from molecular clusters
    remains_table.sort()
    remains_table = merged_node_table.loc[remains_table]
    remains_table_pos = remains_table[remains_table['ion_mode'] == "POS"]
    remains_table_neg = remains_table[remains_table['ion_mode'] == "NEG"]

    print('Linking POS singletons to NEG singletons...')
    candidates_table = list()
    for i in tqdm(remains_table_pos.index):
        ion_1_mz = remains_table_pos.loc[i, mz_field]
        ion_1_rt = remains_table_pos.loc[i, rt_field]
        ion_1_samples = set(remains_table_pos.loc[i, "samples"].split('|'))
        neutral_table = neutral_tabler(ion_1_mz, adduct_table_base_pos)
        hit_table_rt = remains_table_neg[remains_table_neg[rt_field].between(ion_1_rt - rt_error, ion_1_rt + rt_error, inclusive = "both")].copy()
        shares_samples_list = list()
        for j in hit_table_rt.index:
            ion_samples = hit_table_rt.loc[j, "samples"].split('|')
            shared_samples = len(list(ion_1_samples.intersection(ion_samples)))
            shares_samples_list.append(shared_samples)
        hit_table_rt['shared_samples'] = shares_samples_list
        hit_table_rt = hit_table_rt[hit_table_rt['shared_samples'] > 0]
        if len(hit_table_rt) == 0 : continue
        for j in neutral_table.index:
            ion_table = get_ion_table(neutral_table.loc[j, "neutral_mass"], adduct_table_base_neg)
            for k in ion_table.index:
                hit_table = hit_table_rt[hit_table_rt[mz_field].between(ion_table.loc[k, "ion_mz"] - mass_error, ion_table.loc[k, "ion_mz"] + mass_error, inclusive = "both")]
                for l in hit_table.index:
                    candidates_table.append((i,
                                             neutral_table.loc[j, "Adduct"],
                                             neutral_table.loc[j, "Complexity"],
                                             l,
                                             ion_table.loc[k, "Adduct"],
                                             ion_table.loc[k, "Complexity"]))
    candidates_table = pd.DataFrame(candidates_table, columns = ['pos_ion', 'pos_adduct', 'pos_complexity', 'neg_ion', 'neg_adduct', 'neg_complexity'])

    # Get the best hypotheses
    unique_pos = candidates_table['pos_ion'].unique().tolist()
    tmp_candidates_table = list()
    print('Selecting best POS annotations...')
    for i in tqdm(unique_pos) :
        tmp_table_1 = candidates_table[candidates_table['pos_ion'] == i]
        unique_adducts = tmp_table_1['pos_adduct'].unique().tolist()
        point_list = list()
        for adduct in unique_adducts:
            tmp_table_2 = tmp_table_1[tmp_table_1['pos_adduct'] == adduct]
            points = 2*len(tmp_table_2) / (tmp_table_2['pos_complexity'].sum() + tmp_table_2['neg_complexity'].sum())
            point_list.append(points)
        selected_adduct = unique_adducts[point_list.index(max(point_list))]
        tmp_table_1 = tmp_table_1[tmp_table_1['pos_adduct'] == selected_adduct]
        for j in tmp_table_1.index:
            tmp_candidates_table.append((tmp_table_1.loc[j, "pos_ion"],
                                         tmp_table_1.loc[j, "pos_adduct"],
                                         tmp_table_1.loc[j, "pos_complexity"],
                                         tmp_table_1.loc[j, "neg_ion"],
                                         tmp_table_1.loc[j, "neg_adduct"],
                                         tmp_table_1.loc[j, "neg_complexity"]))
    candidates_table = pd.DataFrame(tmp_candidates_table, columns = ['pos_ion', 'pos_adduct', 'pos_complexity', 'neg_ion', 'neg_adduct', 'neg_complexity'])

    # Resolve NEG ions with multiple annotations
    unique_neg = candidates_table['neg_ion'].unique().tolist()
    tmp_candidates_table = list()
    print('Selecting best NEG annotations...')
    for i in tqdm(unique_neg):
        tmp_table_1 = candidates_table[candidates_table['neg_ion'] == i]
        unique_adducts = tmp_table_1['neg_adduct'].unique().tolist()
        point_list = list()
        for adduct in unique_adducts:
            tmp_table_2 = tmp_table_1[tmp_table_1['neg_adduct'] == adduct]
            points = 2*len(tmp_table_2) / (tmp_table_2['pos_complexity'].sum() + tmp_table_2['neg_complexity'].sum())
            point_list.append(points)
        selected_adduct = unique_adducts[point_list.index(max(point_list))]
        tmp_table_1 = tmp_table_1[tmp_table_1['neg_adduct'] == selected_adduct]
        for j in tmp_table_1.index:
            tmp_candidates_table.append((tmp_table_1.loc[j, "pos_ion"],
                                         tmp_table_1.loc[j, "pos_adduct"],
                                         tmp_table_1.loc[j, "pos_complexity"],
                                         tmp_table_1.loc[j, "neg_ion"],
                                         tmp_table_1.loc[j, "neg_adduct"],
                                         tmp_table_1.loc[j, "neg_complexity"]))
    tmp_candidates_table = pd.DataFrame(tmp_candidates_table, columns = ['pos_ion', 'pos_adduct', 'pos_complexity', 'neg_ion', 'neg_adduct', 'neg_complexity'])
    candidates_table = candidates_table.append(tmp_candidates_table, ignore_index = True)

    # Keep only uniques
    candidates_table['mix_idx'] = candidates_table['neg_ion'].astype(str) + "&" + candidates_table['pos_ion'].astype(str)
    tmp_candidates_table = list()
    for comb in candidates_table['mix_idx'].unique():
        tmp_candidates_table.append((candidates_table[candidates_table['mix_idx'] == comb].iloc[0].tolist()))    
    candidates_table = pd.DataFrame(tmp_candidates_table, columns = ['pos_ion', 'pos_adduct', 'pos_complexity', 'neg_ion', 'neg_adduct', 'neg_complexity', 'mix_idx'])    

    # Make a neutral table:
    neutral_idx = 0
    neutral_table = list()
    while len(candidates_table) > 0 :
        idx = candidates_table.index[0]
        pos_pool = list()
        neg_pool = list()
        idx_pool = list()
        pos_pool.append(candidates_table.loc[idx, "pos_ion"])
        neg_pool.append(candidates_table.loc[idx, "neg_ion"])
        new_pos_pool = pos_pool
        new_neg_pool = neg_pool
        idx_pool.append(idx)
        len_0 = 0
        while len_0 != len(idx_pool):
            len_0 = len(idx_pool)
            for i in new_pos_pool:
                idx_pool += candidates_table.index[candidates_table['pos_ion'] == i].tolist()
            for i in new_neg_pool:
                idx_pool += candidates_table.index[candidates_table['neg_ion'] == i].tolist()
            idx_pool = list(set(idx_pool))
        tmp_table_1 = candidates_table.loc[idx_pool]
        
        for i in tmp_table_1['pos_ion'].unique():
            idx = tmp_table_1.index[tmp_table_1['pos_ion'] == i][0]
            spec_id = int(merged_node_table.loc[i, "spec_id"])
            ion_mz = merged_node_table.loc[i, mz_field]
            ion_rt = merged_node_table.loc[i, rt_field]
            ion_tic = merged_node_table.loc[i, "TIC"]
            ion_adduct = tmp_table_1.loc[idx, "pos_adduct"]
            ion_samples = merged_node_table.loc[i, "samples"]
            adduct_idx = adduct_table_base_pos.index[adduct_table_base_pos["Adduct"] == ion_adduct][0]
            ion_adduct_code = adduct_table_base_pos.loc[adduct_idx, "Adduct_code"]
            mol_mass = neutral_mass_calculator(ion_mz,
                                adduct_table_base_pos.loc[adduct_idx, "Adduct_mass"],
                                adduct_table_base_pos.loc[adduct_idx, "Mol_multiplier"],
                                adduct_table_base_pos.loc[adduct_idx, "Charge"])
            Species_rule = Validator_choice(ion_adduct_code, "POS")
            ion_rule_points = Species_rule(prec_mass_error, spec_id, pos_specturm_list)
            neutral_table.append((neutral_idx, mol_mass, ion_mz, ion_rt, i, ion_adduct, ion_tic, ion_rule_points, ion_samples, "POS"))
            
        for i in tmp_table_1['neg_ion'].unique():
            idx = tmp_table_1.index[tmp_table_1['neg_ion'] == i][0]
            spec_id = int(merged_node_table.loc[i, "spec_id"])
            ion_mz = merged_node_table.loc[i, mz_field]
            ion_rt = merged_node_table.loc[i, rt_field]
            ion_tic = merged_node_table.loc[i, "TIC"]
            ion_adduct = tmp_table_1.loc[idx, "neg_adduct"]
            ion_samples = merged_node_table.loc[i, "samples"]
            adduct_idx = adduct_table_base_neg.index[adduct_table_base_neg["Adduct"] == ion_adduct][0]
            ion_adduct_code = adduct_table_base_neg.loc[adduct_idx, "Adduct_code"]
            mol_mass = neutral_mass_calculator(ion_mz,
                                adduct_table_base_neg.loc[adduct_idx, "Adduct_mass"],
                                adduct_table_base_neg.loc[adduct_idx, "Mol_multiplier"],
                                adduct_table_base_neg.loc[adduct_idx, "Charge"])
            Species_rule = Validator_choice(ion_adduct_code, "NEG")
            ion_rule_points = Species_rule(prec_mass_error, spec_id, neg_spectrum_list)
            neutral_table.append((neutral_idx, mol_mass, ion_mz,  ion_rt, i, ion_adduct, ion_tic, ion_rule_points, ion_samples, "NEG"))            
        
        candidates_table.drop(idx_pool, inplace = True)
        neutral_idx += 1
    
    neutral_table = pd.DataFrame(neutral_table, columns = ['neutral_idx', 'neutral_mass', 'ion_mz', 'ion_rt', 'ion_idx', 'adduct', 'TIC', 'rule_points', 'samples', 'ion_mode'])         

    # Report the results: 
    print('Reporting results for POS singletons paired to NEG singletons...')
    for i in tqdm(neutral_table['neutral_idx'].unique()):
        tmp_table_1 = neutral_table[neutral_table['neutral_idx'] == i].copy()
        mol_mass = tmp_table_1['neutral_mass'].mean()
        mol_rt = tmp_table_1['ion_rt'].mean()
        mol_tic = tmp_table_1['TIC'].sum()
        mol_cluster = merged_node_table['cluster_id'].max() + 1
        adduct_count = len(tmp_table_1)
        tmp_table_1['mz_gap'] = abs(mol_mass - tmp_table_1['ion_mz'])
        tmp_table_1['rt_gap'] = abs(mol_rt - tmp_table_1['ion_rt'])
        mol_samples = '|'.join(tmp_table_1['samples'])
        mol_samples = list(set(mol_samples.split('|')))
        mol_samples.sort()
        mol_samples = '|'.join(mol_samples)
        new_idx = merged_node_table.index.max() + 1
        merged_node_table.loc[new_idx] = [None]*len(merged_node_table.columns)
        merged_node_table.loc[new_idx, [mz_field, rt_field, "TIC", "rule_points",
                                        'status', 'Adnotation', 'adduct_count',
                                        'ion_mode', 'samples', 'cluster_id']] = [mol_mass,
                                        mol_rt, mol_tic, 0, "neutral", None, adduct_count,
                                        "MIX", mol_samples, mol_cluster]
        
        for j in tmp_table_1.index:
            ion_idx = tmp_table_1.loc[j, 'ion_idx']
            adduct = tmp_table_1.loc[j, 'adduct']
            rule_points = tmp_table_1.loc[j, 'rule_points']
            mz_gap = tmp_table_1.loc[j, 'mz_gap']
            rt_gap = tmp_table_1.loc[j, 'rt_gap']
            tmp_mode = tmp_table_1.loc[j, "ion_mode"]
            merged_node_table.loc[ion_idx, 'rule_points'] = rule_points
            merged_node_table.loc[ion_idx, 'status'] = "adduct"
            merged_node_table.loc[ion_idx, 'Adnotation'] = adduct
            merged_node_table.loc[ion_idx, 'cluster_id'] = mol_cluster
            new_edge = merged_edge_table.index.max() + 1

            merged_edge_table.loc[new_idx] = [None]*len(merged_edge_table.columns)
            merged_edge_table.loc[new_idx, ['node_1', 'node_2', 'matched_peaks',
                                            'total_peaks', 'matching_score',
                   'rt_gap', 'mz_gap', 'status', 'Fragnotation', 'Adnotation',
                   'All_annotations', 'ion_mode']] = [new_idx, ion_idx, 0 ,0 ,0 ,
                    rt_gap, mz_gap, "add_edge", None, adduct, adduct, tmp_mode]
            del_edge = merged_edge_table[merged_edge_table['node_1'] == ion_idx]
            del_edge = del_edge[del_edge['node_2'] == ion_idx].index
            if len(del_edge) > 0 :
                merged_edge_table.drop(del_edge[0], inplace = True)
    
    # Update status node nodes (pos neutrals, neg neutrals, pos adducts_ neg_adducts)
    merged_node_table.insert(merged_node_table.columns.get_loc('status') + 1, 
                             "status_universal", merged_node_table['status'].copy())
    merged_node_table['status'] = merged_node_table['ion_mode'].str.lower() + "_" + merged_node_table['status']
    
    # Update status for edges (pos add, neg add, pos frag, neg frag, pos single, neg single)
    merged_edge_table.insert(merged_edge_table.columns.get_loc('status') + 1, 
                             "status_universal", merged_edge_table['status'].copy())
    merged_edge_table['status'] = merged_edge_table['ion_mode'].str.lower() + "_" + merged_edge_table['status']
    
    # Round values:
    merged_node_table[mz_field] = merged_node_table[mz_field].round(4)
    merged_node_table[rt_field] = merged_node_table[rt_field].round(3)
    merged_edge_table['mz_gap'] = merged_edge_table['mz_gap'].round(4)
    merged_node_table.drop("samples", axis = 1, inplace = True)

    # Adding csv columns to the node table (sample intensities)   
    neutrals_idx = list(merged_node_table.index[merged_node_table['status'] == "neg_neutral"])
    neutrals_idx += list(merged_node_table.index[merged_node_table['status'] == "pos_neutral"])
    neutrals_idx += list(merged_node_table.index[merged_node_table['status'] == "mix_neutral"])
    pos_ions_idx = list(set(merged_node_table.index[merged_node_table['ion_mode'] == "POS"]) - set(neutrals_idx))
    neg_ions_idx = list(set(merged_node_table.index[merged_node_table['ion_mode'] == "NEG"]) - set(neutrals_idx))
    samples = list(neg_table_global.columns) + list(pos_table_global.columns)
    samples = list(set(samples))
    samples.sort()
    
    for sample in samples:
        merged_node_table[sample] = [0.0]*len(merged_node_table)
        
    print("Adding POS sample intensities to the merged node table...")
    for i in tqdm(pos_ions_idx):
        ion_id = merged_node_table.loc[i, idx_column]
        for sample in samples:
            merged_node_table.loc[i, sample] = pos_table_global.loc[ion_id, sample]
    
    print("Adding NEG sample intensities to the merged node table...")
    for i in tqdm(neg_ions_idx):
        ion_id = merged_node_table.loc[i, idx_column]
        for sample in samples:
            merged_node_table.loc[i, sample] = neg_table_global.loc[ion_id, sample]
    
    print("Adding neutral sample intensities to the merged node table...")
    for i in tqdm(neutrals_idx):
        ion_ids = list(merged_edge_table['node_2'][merged_edge_table['node_1'] == i])
        for sample in samples:
            merged_node_table.loc[i, sample] = merged_node_table.loc[ion_ids, sample].sum()

    # Export the data
    merged_node_table.to_csv(out_full + 'node_table.csv', index_label = "Index")
    merged_edge_table.to_csv(out_full + 'edge_table.csv', index_label = "Index")
    
    if params['mm_export_samples'] : 
        mm_samplewise_export(neg_csv_file = in_path_full_neg + neg_csv,
                          pos_csv_file = in_path_full_pos + pos_csv,
                          out_path = out_samples,
                          edge_table = merged_edge_table,
                          node_table = merged_node_table,
                          params = params)
    return



