import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters
from matchms.similarity import ModifiedCosine
from MolNotator.others.global_functions import *

def dereplicator(params : dict, db_params : dict):

    # Load parameters
    mz_field = params['mz_field']
    rt_field = params['rt_field']
    mgf_file_neg= params['neg_mgf']
    mgf_file_pos= params['pos_mgf']
    neg_csv_file= params['neg_csv']
    pos_csv_file= params['pos_csv']
    input_mgf_neg_path= params['neg_out_0']
    input_mgf_pos_path= params['pos_out_0']
    node_table_path= params['mix_out_4_1']
    out_path_full= params['mix_out_5_1']
    out_path_samples= params['mix_out_5_2']
    
    db_prefix= db_params['db_prefix']
    db_type= db_params['db_type']
    db_folder= db_params['db_folder']
    db_file= db_params['db_file']
    db_cosine= db_params['db_cosine']
    db_matched_peaks= db_params['db_matched_peaks']
    db_prec_error= db_params['db_prec_error']
    db_mass_error= db_params['db_mass_error']
    db_rt_error= db_params['db_rt_error']
    db_hits= db_params['db_hits']
    db_adduct_filter= db_params['db_adduct_filter']
    db_adduct_field= db_params['db_adduct_field']
    db_rt_filter= db_params['db_rt_filter']
    db_mode_field= db_params['db_mode_field']
    db_unique_field= db_params['db_unique_field']
    db_export_fields= db_params['db_export_fields']
    db_prefix= db_params['db_prefix']

    # If the retention time is in minutes
    if params['rt_unit'] == "m":
        db_rt_error = db_rt_error/60

    #Create folders and load tables
    if not os.path.isdir(out_path_full) :
        os.mkdir(out_path_full)
        edge_table = pd.read_csv(node_table_path + 'edge_table.csv', index_col = "Index")
        node_table = pd.read_csv(node_table_path + 'node_table.csv', index_col = "Index",
                                 low_memory = False)
    else:
        edge_table = pd.read_csv(out_path_full + 'edge_table.csv', index_col = "Index")
        node_table = pd.read_csv(out_path_full + 'node_table.csv', index_col = "Index",
                                 low_memory = False)

    if not os.path.isdir(out_path_samples) :
        os.mkdir(out_path_samples) 

    # Process tables and remove samples to be stored in a samples DF
    node_table['Adnotation'] = node_table['Adnotation'].replace({np.nan: None})
    
    ###########################################################################
    if db_type == 'ion': # If database is MGF to dereplicate ion (MS/MS) data
    ###########################################################################
        modified_cosine = ModifiedCosine(tolerance=db_mass_error)
        Filter_fields_outer, filter_cols = Filter_choices_outer(db_params)
        
        # Load MGF files
        if params['process_mode'] == "NEG":
            print('Loading NEG MGF file...')
            mgf_neg = list(load_from_mgf(input_mgf_neg_path + mgf_file_neg))
            mgf_neg = [Spectrum_processing(s) for s in mgf_neg]
        elif params['process_mode'] == "POS":
            print('Loading POS MGF file...')
            mgf_pos = list(load_from_mgf(input_mgf_pos_path + mgf_file_pos))
            mgf_pos = [Spectrum_processing(s) for s in mgf_pos]
        else:
            print('Loading NEG MGF file...')
            mgf_neg = list(load_from_mgf(input_mgf_neg_path + mgf_file_neg))
            mgf_neg = [Spectrum_processing(s) for s in mgf_neg]
            print('Loading POS MGF file...')
            mgf_pos = list(load_from_mgf(input_mgf_pos_path + mgf_file_pos))
            mgf_pos = [Spectrum_processing(s) for s in mgf_pos]

        # Load the database file
        print('Loading database file and extracting data...')
        database_mgf = list(load_from_mgf(db_folder + db_file))
        database_mgf = [Float_prec_mz(s, db_params) for s in database_mgf]
        database_table = Database_table_mgf(database_mgf, db_params)
        
        # Start dereplication (cosine similarity)
        derep_table = list()
        print('Starting ion dereplication (cosine similarity)...')
        for i in tqdm(node_table.index):
            status = node_table.loc[i, "status_universal"]
            if status == "neutral" : 
                new_row = [i, None, None] + [None]*len(filter_cols) + [None]*len(db_export_fields) + [None, None, None]
                derep_table.append(new_row)                    
                continue
            ion_mz = node_table.loc[i, mz_field]
            ion_mgf_idx = int(node_table.loc[i, "spec_id"])
            ion_mode = node_table.loc[i, "ion_mode"]
            hits = database_table[database_table["mz"].between(ion_mz - db_prec_error, ion_mz + db_prec_error, inclusive = "both")].copy()
            
            # Ion mode filter
            if ion_mode == "NEG":
                hits = hits[hits['ion_mode'] == "negative"]
                exp_mgf = mgf_neg
            else:
                hits = hits[hits['ion_mode'] == "positive"]
                exp_mgf = mgf_pos
            
            # Optional filters : Adduct and RT
            if db_adduct_filter and (status == 'adduct'):
                adduct = node_table.loc[i, 'Adnotation']
                if adduct != None:
                    hits = hits[hits['adduct'] == adduct]
            if db_rt_filter:
                rt = node_table.loc[i, rt_field]
                hits = hits[hits["rt"].between(rt - db_rt_error, rt + db_rt_error, inclusive = "both")]
            
            # Calculate cosine similarity if hit table is not empty
            similarity_list = list()
            for j in hits.index:
                score, n_matches = modified_cosine.pair(exp_mgf[ion_mgf_idx], database_mgf[j])
                mass_error = abs(ion_mz - hits.loc[j, "mz"])*1000
                prod = score * n_matches
                similarity_list.append((j, score, n_matches, prod, mass_error))
            similarity_list = pd.DataFrame(similarity_list, columns = ["index", "cos", "matches", "prod", "error_mDa"])
            similarity_list.set_index('index', inplace = True)
            hits = hits.join(similarity_list, hits.index)
            
            # Filter using cosine and matched peaks thresholds
            hits = hits[hits['cos'] >= db_cosine]
            hits = hits[hits['matches'] >= db_matched_peaks]

            # If no hits at this point, export None and continue
            if len(hits) == 0 : 
                new_row = [i, None, None] + [None]*len(filter_cols) + [None]*len(db_export_fields) + [None, None, None]
                derep_table.append(new_row)    
                continue
            
            # Filter using unique field:
            unique_mols = hits['unique_field'].dropna().unique()
            mol_idx = hits.index[hits['unique_field'].isnull()].tolist()
            for mol in unique_mols:
                tmp_table_1 = hits[hits['unique_field'] == mol].copy()
                tmp_table_1.sort_values('prod', ascending = False, inplace = True)
                mol_idx.append(tmp_table_1.index[0])
            hits = hits.loc[mol_idx]
            
            # Keep the top db_hits hit in the hits table
            hits.sort_values('prod', ascending = False, inplace = True)
            hits = hits.iloc[:db_hits]
            hits = hits.fillna('')
            
            # Report the results in the derep_table list
            tmp_name = '|'.join(hits['name'])
            tmp_mz = '|'.join(hits['mz'].round(4).astype(str))
            tmp_filter_fields = ['|'.join(hits[col].astype(str)) for col in filter_cols]
            tmp_export_fields = ['|'.join(hits[f].astype(str)) for f in db_export_fields]
            tmp_cos = '|'.join(hits['cos'].round(2).astype(str))
            tmp_matches = '|'.join(hits['matches'].astype(str))
            tmp_error = '|'.join(hits['error_mDa'].round(1).astype(str))
            new_row = [i, tmp_name, tmp_mz] + tmp_filter_fields + tmp_export_fields + [tmp_cos, tmp_matches, tmp_error]
            derep_table.append(new_row)
        
        cols = ['index', 'name', 'mass'] + filter_cols + db_export_fields + ['cos', 'matches', 'error_mDa']
        cols = [db_prefix + '_' + col for col in cols]
        derep_table = pd.DataFrame(derep_table, columns = cols)
        derep_table.set_index(db_prefix + '_index', drop = True, inplace = True)
        
        # Report dereplication data to the node table
        cols.remove(db_prefix + '_index')
        if len(set(cols).intersection(node_table.columns)) > 0 :
            node_table.drop(cols, axis = 1, inplace = True)
        node_table = node_table.join(derep_table, node_table.index)
        
        # Propagate ion dereplications to neutrals
        neutrals_idx = node_table.index[node_table['status_universal'] == "neutral"]
        print('Propagating ion annotations to neutral nodes...')
        for i in tqdm(neutrals_idx):
            tmp_ions = edge_table['node_2'][edge_table['node_1'] == i].tolist()
            tmp_nodes = node_table.loc[tmp_ions]
            tmp_nodes = tmp_nodes[~tmp_nodes[db_prefix + '_mass'].isnull()]
            if len(tmp_nodes) == 0 : continue
            derep_table = list()
            for j in tmp_nodes.index:
                for k in range(len(tmp_nodes.loc[j, db_prefix + '_mass'].split('|'))):
                    new_row = list()
                    for col in cols:
                        new_row.append(tmp_nodes.loc[j, col].split('|')[k])
                    derep_table.append(new_row)
            derep_table = pd.DataFrame(derep_table, columns = cols)

            derep_table['prod'] = derep_table[db_prefix + '_cos'].astype(float) * derep_table[db_prefix + '_matches'].astype(int)
            unique_mols = derep_table[db_prefix + '_' + db_unique_field].dropna().unique()
            mol_idx = derep_table.index[derep_table[db_prefix + '_' + db_unique_field].isnull()].tolist()
            for mol in unique_mols:
                tmp_table_1 = derep_table[derep_table[db_prefix + '_' + db_unique_field] == mol].copy()
                tmp_table_1.sort_values('prod', ascending = False, inplace = True)
                mol_idx.append(tmp_table_1.index[0])
            derep_table = derep_table.loc[mol_idx]
            derep_table.sort_values("prod", ascending = False, inplace = True)
            derep_table = derep_table.iloc[:db_hits]
            derep_table = derep_table.fillna('')
            for col in cols:
                node_table.loc[i, col] = '|'.join(derep_table[col])
            if db_adduct_field != None:
                node_table.loc[i, db_adduct_field] = None
            
    #############################################################################
    elif db_type == "neutral": # If database is CSV to dereplicate molecule nodes
    #############################################################################
        # Load the database file
        print('Loading database file and extracting data...')
        if db_file[-3:] == 'tsv':
            sep = '\t'
        elif db_file[-3:] == 'csv':
            sep = ','
        if db_params['db_index_field'] != None :
            database_csv = pd.read_csv(db_folder + db_file, sep = sep, index_col = db_params['db_index_field'], low_memory = False)
        else:
            database_csv = pd.read_csv(db_folder + db_file, sep = sep, low_memory = False)
        database_table = Database_table_csv(database_csv, db_params)
        
        # Start dereplication (on neutrals only)
        derep_table = list()
        print('Starting neutral dereplication...')
        for i in tqdm(node_table.index):
            if node_table.loc[i, "status_universal"] == "neutral":
                mass = node_table.loc[i, mz_field]
                hits = database_table[database_table['mass'].between(mass - db_prec_error, mass + db_prec_error, inclusive = "both")].copy()
                if db_adduct_filter :
                    rt = node_table.loc[i, rt_field]
                    hits = hits[hits['rt'].between(rt - db_rt_error, rt + db_rt_error, inclusive = "both")]
                if len(hits) == 0 : 
                    new_row = [i, None, None, None] + [None]*len(db_params['db_export_fields'])
                    derep_table.append(new_row)    
                    continue
                hits['error_mDa'] = abs(mass - hits['mass'])*1000
                hits.sort_values('error_mDa', ascending = True, inplace = True)
                hits = hits.iloc[:db_hits]
                tmp_names = '|'.join(hits['name'])
                tmp_masses = '|'.join(hits['mass'].round(4).astype(str))
                tmp_error = '|'.join(hits['error_mDa'].round(1).astype(str))
                other_fields = ['|'.join(hits[field].astype(str)) for field in db_params['db_export_fields']]
                new_row = [i, tmp_names, tmp_masses, tmp_error]  + other_fields
                derep_table.append(new_row)
            else:
                new_row = [i, None, None, None] + [None]*len(db_params['db_export_fields'])
                derep_table.append(new_row)
        cols = ['index', 'name', 'mass', 'error_mDa'] + db_params['db_export_fields']
        cols = [db_prefix + '_' + col for col in cols]
        derep_table = pd.DataFrame(derep_table, columns = cols)
        derep_table.set_index(db_prefix + '_index', drop = True, inplace = True)
        
        # Report dereplication data to the node table
        node_table = node_table.join(derep_table, node_table.index)
        
    #############################################################################
    elif db_type == "formula": # CSV database dereplication based on chem formulas for neutrals and some ions
    #############################################################################

        # Load the database file
        print('Loading database file and extracting data...')
        if db_file[-3:] == 'tsv':
            sep = '\t'
        elif db_file[-3:] == 'csv':
            sep = ','
        if db_params['db_index_field'] != None :
            database_csv = pd.read_csv(db_folder + db_file, sep = sep, index_col = db_params['db_index_field'])
        else:
            database_csv = pd.read_csv(db_folder + db_file, sep = sep)
        database_table = Database_table_csv(database_csv, db_params)
        
        # Start dereplication (on neutrals only)
        derep_table = list()
        print('Starting formula dereplication...')
        for i in tqdm(node_table.index):
            formula = node_table.loc[i, "formula"]
            if type(formula) == float:
                new_row = [i, None, None, None] + [None]*len(db_params['db_export_fields'])
                derep_table.append(new_row)
                continue
            hits = database_table[database_table['formula'] == formula].copy()
            hits = hits.iloc[:db_hits]
            if len(hits) == 0 : 
                new_row = [i, None, None, None] + [None]*len(db_params['db_export_fields'])
                derep_table.append(new_row)
                continue                
            tmp_names = '|'.join(hits['name'])
            tmp_masses = str(hits['mass'].iloc[0].round(4))
            tmp_formula = hits['formula'].iloc[0]
            other_fields = ['|'.join(hits[field].astype(str)) for field in db_params['db_export_fields']]
            new_row = [i, tmp_names, tmp_masses, tmp_formula]  + other_fields
            derep_table.append(new_row)

        cols = ['index', 'name', 'mass', 'formula'] + db_params['db_export_fields']
        cols = [db_prefix + '_' + col for col in cols]
        derep_table = pd.DataFrame(derep_table, columns = cols)
        derep_table.set_index(db_prefix + '_index', drop = True, inplace = True)
        
        # Report dereplication data to the node table
        node_table = node_table.join(derep_table, node_table.index)
    
    ###########################################################################
    # Prepare data to be exported
    ###########################################################################
    
    # Round values where necessary and return the samples to the dataframe
    node_table[mz_field] = node_table[mz_field].round(4)
    node_table[rt_field] = node_table[rt_field].round(3)
    
    # Export data
    node_table.to_csv(out_path_full + 'node_table.csv', index_label = "Index")
    edge_table.to_csv(out_path_full + 'edge_table.csv', index_label = "Index")
    
    if db_params['db_export_samples'] : 
        if params['process_mode'] == "BOTH" :
            dr_samplewise_export(neg_csv_file = input_mgf_neg_path + neg_csv_file,
                              pos_csv_file = input_mgf_pos_path + pos_csv_file,
                              out_path = out_path_samples,
                              edge_table = edge_table,
                              node_table = node_table,
                              params = params)
        else :
            samplewise_export(merged_edge_table = edge_table,
                              merged_node_table = node_table,
                              step = "dereplicator",
                              ion_mode = params['process_mode'],
                              params = params)
    return
